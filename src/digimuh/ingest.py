#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║          ____  _       _ __  __       _                          ║
# ║         |  _ \(_) __ _(_)  \/  |_   _| |__                      ║
# ║         | | | | |/ _` | | |\/| | | | | '_ \                     ║
# ║         | |_| | | (_| | | |  | | |_| | | | |                    ║
# ║         |____/|_|\__, |_|_|  |_|\__,_|_| |_|                    ║
# ║                  |___/                                           ║
# ║                                                                  ║
# ║         « CSV → SQLite ingestion engine »                        ║
# ╚══════════════════════════════════════════════════════════════════╝
"""DigiMuh ingestion engine — CSV to SQLite star-schema builder.

Reads ~8.9 GB of heterogeneous dairy-cow CSV data exported from
multiple on-farm sensor systems (smaXtec, HerdePlus, Gouna, HOBO,
LoRaWAN, DWD weather) and consolidates them into a single normalised
SQLite database with dimension tables for animals, sensors, barns,
and source-file provenance.

The resulting database uses a star schema: small dimension tables
hold entity metadata; large fact tables hold time-series measurements
with foreign keys pointing back to the dimensions.  Timestamps are
stored as ISO-8601 TEXT (SQLite's date/time functions work natively
on this representation).  Composite indexes on ``(entity_id, timestamp)``
are created for every fact table to accelerate the most common query
pattern: *"give me all X for animal Y between dates A and B"*.

Usage::

    # Smoke test — first 5 files per folder
    python -m digimuh.ingest_cow_db /path/to/csv/root --db cow.db --test-n 5

    # Full ingestion (~2–3 h depending on disk I/O)
    python -m digimuh.ingest_cow_db /path/to/csv/root --db cow.db

    # Verbose mode — prints CREATE TABLE SQL for each table
    python -m digimuh.ingest_cow_db /path/to/csv/root --db cow.db -v

Dependencies:
    * Python ≥ 3.10 (for ``X | Y`` union type hints)
    * ``tqdm`` (optional, for progress bars)

See Also:
    * ``docs/database_structure.md`` — full schema documentation
    * ``docs/column_dictionary.md`` — column-level data dictionary

Authors:
    Bart R. H. Geurten, Claude (Anthropic)
"""

# ─────────────────────────────────────────────────────────────────
#  « imports »
# ─────────────────────────────────────────────────────────────────
from __future__ import annotations

import argparse
import csv
import logging
import re
import sqlite3
import sys
import time
from pathlib import Path

# ─────────────────────────────────────────────────────────────────
#  « optional but nice »  tqdm for progress bars
# ─────────────────────────────────────────────────────────────────
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

    class tqdm:  # noqa: N801
        """Minimal tqdm stand-in so the script runs without it.

        Provides the iterator protocol, a context manager, and the
        ``set_postfix_str`` / ``update`` stubs so call sites don't
        need to branch on availability.
        """

        def __init__(self, iterable=None, **kw):
            self._it = iterable
            self._desc = kw.get("desc", "")
            self._total = kw.get("total", None)
            self._n = 0

        def __iter__(self):
            for item in self._it:
                self._n += 1
                yield item

        def __enter__(self):
            return self

        def __exit__(self, *a):
            pass

        def update(self, n: int = 1) -> None:
            self._n += n

        def set_postfix_str(self, s: str) -> None:
            pass


# ═════════════════════════════════════════════════════════════════
#  « configuration »  folder → table mapping
# ═════════════════════════════════════════════════════════════════
#
#  Each entry maps a directory name (relative to the CSV root) to
#  a target SQLite table and describes how entity IDs are obtained.
#
#  entity_source values:
#    "filename_animal"  — first '_'-delimited segment is an EU ear
#                         tag (15-digit integer → animals.animal_id)
#    "filename_sensor"  — first segment is a sensor name string
#                         (e.g. "CU-1" → sensors.sensor_id)
#    "filename_barn"    — first segment is a barn name string
#                         (e.g. "NewBridge" → barns.barn_id)
#    "column_animal_id" — animal_id column already in the CSV
#    "column_cow"       — 'cow' column in the CSV (renamed →
#                         animal_id during ingestion)
#    None               — no entity association (weather data etc.)
#

FOLDER_CONFIG: dict[str, dict] = {
    "output_allocations": {
        "table": "allocations",
        "entity_source": "column_animal_id",
        "entity_col": "animal_id",
        "glob": "*.csv",
    },
    "outputs_bcs": {
        "table": "bcs",
        "entity_source": "filename_animal",
        "glob": "*.csv",
    },
    "outputs_gouna": {
        "table": "gouna",
        "entity_source": "filename_animal",
        "glob": "*.csv",
    },
    "outputs_herdeplus_mlp_gemelk_kalbung": {
        "table": "herdeplus",
        "entity_source": "filename_animal",
        "glob": "*.csv",
    },
    "outputs_lorawan": {
        "table": "lorawan",
        "entity_source": "filename_sensor",
        "glob": "*.csv",
    },
    "outputs_smaxtec_barns": {
        "table": "smaxtec_barns",
        "entity_source": "filename_barn",
        "glob": "*.csv",
    },
    "outputs_smaxtec_derived": {
        "table": "smaxtec_derived",
        "entity_source": "filename_animal",
        "glob": "*.csv",
    },
    "outputs_smaxtec_events": {
        "table": "smaxtec_events",
        "entity_source": "filename_animal",
        "glob": "*.csv",
    },
    "outputs_smaxtec_water_intake": {
        "table": "smaxtec_water_intake",
        "entity_source": "filename_animal",
        "glob": "*.csv",
    },
    "outputs_hobo": {
        "table": "hobo_weather",
        "entity_source": None,
        "glob": "*.csv",
    },
}
"""dict: Mapping of folder names to ingestion configuration.

Each value dict contains:
    * ``table`` — target SQLite table name
    * ``entity_source`` — how the entity ID is derived
    * ``entity_col`` — (optional) CSV column holding the entity ID
    * ``glob`` — file-glob pattern for discovering CSVs
"""

STANDALONE_FILES: dict[str, dict] = {
    "outputs_dwd.csv": {
        "table": "dwd_weather",
        "entity_source": None,
    },
    "herdeplus_diseases.csv": {
        "table": "diseases",
        "entity_source": "column_cow",
        "entity_col": "cow",
    },
}
"""dict: Single-file sources that live in the root directory."""


# ═════════════════════════════════════════════════════════════════
#  « dimension table DDL »
# ═════════════════════════════════════════════════════════════════

DIMENSION_TABLES_SQL: str = """
-- ┌──────────────────────────────────────────────────────────────┐
-- │  Dimension tables                                           │
-- └──────────────────────────────────────────────────────────────┘

CREATE TABLE IF NOT EXISTS animals (
    animal_id  INTEGER PRIMARY KEY   -- EU ear tag, IS the rowid
);

CREATE TABLE IF NOT EXISTS sensors (
    sensor_id    INTEGER PRIMARY KEY AUTOINCREMENT,
    sensor_name  TEXT    NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS barns (
    barn_id    INTEGER PRIMARY KEY AUTOINCREMENT,
    barn_name  TEXT    NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS source_files (
    file_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    filename  TEXT    NOT NULL,
    folder    TEXT    NOT NULL,
    UNIQUE(filename, folder)
);
"""
"""str: SQL script to create the four dimension tables.

``animals`` uses the 15-digit EU ear tag as its ``INTEGER PRIMARY KEY``
which makes it an alias for SQLite's internal ``rowid`` — the fastest
possible lookup.  The other dimensions use ``AUTOINCREMENT`` surrogates
because their natural keys are short strings.
"""


# ═════════════════════════════════════════════════════════════════
#  « helper functions »
# ═════════════════════════════════════════════════════════════════

def extract_entity_id_from_filename(filename: str) -> str:
    """Extract the first underscore-delimited segment from a filename.

    The DigiMuh CSV export encodes the entity identifier (animal ear
    tag, sensor name, or barn name) as the first segment of the
    filename, separated from the rest by an underscore.

    Args:
        filename: Basename of the CSV file, e.g.
            ``"276001260919234_smaxtec_derived_2021-04-01_2024-09-30.csv"``.

    Returns:
        The first segment, e.g. ``"276001260919234"``.

    Example:
        >>> extract_entity_id_from_filename("CU-1_LoRaWAN_raw_2021.csv")
        'CU-1'
    """
    stem = Path(filename).stem
    return stem.split("_")[0]


def guess_column_type(values: list[str]) -> str:
    """Infer an SQLite column type from a sample of string values.

    Tries to parse every non-empty value as ``int``, then ``float``.
    Falls back to ``TEXT`` if either parse fails.  An all-empty sample
    also returns ``TEXT``.

    Args:
        values: Sample of raw string values from a CSV column.

    Returns:
        One of ``"INTEGER"``, ``"REAL"``, or ``"TEXT"``.

    Note:
        SQLite uses dynamic typing so the column affinity is advisory,
        but correct affinity improves storage efficiency and query
        planner decisions.
    """
    non_empty = [v for v in values if v.strip()]
    if not non_empty:
        return "TEXT"

    # ── try integer ──────────────────────────────────────────────
    int_ok = True
    for v in non_empty:
        try:
            int(v)
        except ValueError:
            int_ok = False
            break
    if int_ok:
        return "INTEGER"

    # ── try float ────────────────────────────────────────────────
    real_ok = True
    for v in non_empty:
        try:
            float(v)
        except ValueError:
            real_ok = False
            break
    if real_ok:
        return "REAL"

    return "TEXT"


def read_csv_sample(
    filepath: Path, n_rows: int = 100,
) -> tuple[list[str], list[list[str]]]:
    """Read headers and up to *n_rows* data rows from a CSV file.

    Used during schema inference to peek at column headers and a
    small sample of data for type guessing.

    Args:
        filepath: Path to the CSV file.
        n_rows: Maximum number of data rows to read.  The default of
            100 is enough for most narrow tables; wide sparse tables
            like ``smaxtec_derived`` benefit from larger samples (see
            :func:`build_create_table_sql`).

    Returns:
        A ``(headers, rows)`` tuple where *headers* is a list of raw
        column-name strings and *rows* is a list of lists of string
        values.
    """
    rows: list[list[str]] = []
    with open(filepath, "r", newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        headers = next(reader)
        for i, row in enumerate(reader):
            if i >= n_rows:
                break
            rows.append(row)
    return headers, rows


def sanitise_column_name(name: str) -> str:
    """Clean a raw CSV header into a safe SQL column name.

    Strips whitespace, lowercases, replaces any character that is not
    ``[a-z0-9_]`` with an underscore, collapses consecutive
    underscores, and strips leading/trailing underscores.

    Args:
        name: Raw column header string from the CSV.

    Returns:
        A sanitised, SQL-safe column name.

    Example:
        >>> sanitise_column_name("BCS-Wert")
        'bcs_wert'
        >>> sanitise_column_name("21141733_1__Temperature")
        '21141733_1_temperature'
    """
    name = name.strip().lower()
    name = re.sub(r"[^a-z0-9_]", "_", name)
    name = re.sub(r"_+", "_", name)
    name = name.strip("_")
    return name


def build_create_table_sql(
    table_name: str,
    csv_headers: list[str],
    csv_rows: list[list[str]],
    entity_source: str | None,
    entity_col: str | None = None,
) -> tuple[str, list[str]]:
    """Build a ``CREATE TABLE`` statement by inspecting CSV structure.

    Examines the raw headers and a sample of data rows to infer column
    types, then constructs the DDL statement with appropriate foreign-key
    references to dimension tables.  Column names are sanitised via
    :func:`sanitise_column_name`.

    Args:
        table_name: Target SQLite table name.
        csv_headers: Raw CSV column headers (unsanitised).
        csv_rows: Sample data rows for type inference.
        entity_source: How the entity ID is obtained — one of the
            values documented in :data:`FOLDER_CONFIG`, or ``None``.
        entity_col: If the entity is identified by a column in the
            CSV itself (``"column_animal_id"`` or ``"column_cow"``),
            this is the name of that column.

    Returns:
        A ``(sql, col_names)`` tuple.  *sql* is the full ``CREATE
        TABLE`` statement.  *col_names* is an ordered list of the
        sanitised CSV column names (matching CSV column order), used
        to build the corresponding ``INSERT`` statement.
    """
    columns: list[tuple[str, str]] = []
    col_names: list[str] = []

    needs_animal_id = entity_source in ("filename_animal",)
    needs_sensor_id = entity_source in ("filename_sensor",)
    needs_barn_id = entity_source in ("filename_barn",)

    for i, raw_header in enumerate(csv_headers):
        san = sanitise_column_name(raw_header)

        # Rename the entity column to animal_id when it comes from CSV
        if entity_col and raw_header.strip().lower() == entity_col.lower():
            if entity_source in ("column_animal_id", "column_cow"):
                san = "animal_id"

        sample_vals = [row[i] for row in csv_rows if i < len(row)]
        col_type = guess_column_type(sample_vals)
        columns.append((san, col_type))
        col_names.append(san)

    # ── assemble CREATE TABLE parts ──────────────────────────────
    parts: list[str] = []

    if needs_animal_id:
        parts.append(
            "    animal_id  INTEGER NOT NULL "
            "REFERENCES animals(animal_id)"
        )
    if needs_sensor_id:
        parts.append(
            "    sensor_id  INTEGER NOT NULL "
            "REFERENCES sensors(sensor_id)"
        )
    if needs_barn_id:
        parts.append(
            "    barn_id    INTEGER NOT NULL "
            "REFERENCES barns(barn_id)"
        )

    for san, col_type in columns:
        if san == "animal_id" and needs_animal_id:
            continue  # already added above
        parts.append(f'    "{san}"  {col_type}')

    parts.append(
        "    file_id  INTEGER NOT NULL "
        "REFERENCES source_files(file_id)"
    )

    sql = f'CREATE TABLE IF NOT EXISTS "{table_name}" (\n'
    sql += ",\n".join(parts)
    sql += "\n);"

    return sql, col_names


# ═════════════════════════════════════════════════════════════════
#  « ingestion engine »
# ═════════════════════════════════════════════════════════════════

class CowDBIngester:
    """Orchestrates reading CSV files and inserting into SQLite.

    Walks the DigiMuh export directory structure, creates dimension
    and fact tables on the fly by inspecting the first CSV in each
    folder, then streams rows in configurable batches using
    ``executemany`` for throughput.

    The class caches all dimension-table lookups in memory (animals,
    sensors, barns, source files) so repeated INSERT-OR-IGNORE round
    trips are avoided for entities already seen.

    Args:
        root_dir: Root directory containing all CSV folders and
            standalone files.
        db_path: Path where the SQLite database will be created.
        chunk_size: Number of rows per ``INSERT`` batch.  Larger
            values use more memory but reduce transaction overhead.
        verbose: If ``True``, emit ``DEBUG``-level log messages
            including full ``CREATE TABLE`` SQL.
        test_n: If set, only ingest the first *test_n* files per
            folder (standalone files are always fully ingested).
            Useful for rapid schema and pipeline validation.

    Example:
        >>> ingester = CowDBIngester(
        ...     root_dir=Path("/data/DigiMuh-Export"),
        ...     db_path=Path("cow.db"),
        ...     chunk_size=50_000,
        ...     test_n=5,
        ... )
        >>> ingester.run()
    """

    def __init__(
        self,
        root_dir: Path,
        db_path: Path,
        chunk_size: int = 50_000,
        verbose: bool = False,
        test_n: int | None = None,
    ) -> None:
        self.root_dir = root_dir
        self.db_path = db_path
        self.chunk_size = chunk_size
        self.verbose = verbose
        self.test_n = test_n
        self.log = logging.getLogger("cow_db")

        # ── dimension caches ─────────────────────────────────────
        self._animal_cache: set[int] = set()
        self._sensor_cache: dict[str, int] = {}
        self._barn_cache: dict[str, int] = {}
        self._file_cache: dict[tuple[str, str], int] = {}

        # ── bookkeeping ──────────────────────────────────────────
        self._created_tables: set[str] = set()

    # ─────────────────────────────────────────────────────────────
    #  « dimension helpers »
    # ─────────────────────────────────────────────────────────────

    def _ensure_animal(self, cur: sqlite3.Cursor, animal_id: int) -> None:
        """Register an animal in the ``animals`` dimension table.

        Uses an in-memory cache to skip the INSERT for animals already
        seen in this session.

        Args:
            cur: Active database cursor.
            animal_id: EU ear tag number (15-digit integer).
        """
        if animal_id not in self._animal_cache:
            cur.execute(
                "INSERT OR IGNORE INTO animals (animal_id) VALUES (?)",
                (animal_id,),
            )
            self._animal_cache.add(animal_id)

    def _ensure_sensor(self, cur: sqlite3.Cursor, name: str) -> int:
        """Get or create a ``sensor_id`` for a sensor name.

        Args:
            cur: Active database cursor.
            name: Sensor name string, e.g. ``"CU-1"``.

        Returns:
            The integer ``sensor_id`` from the ``sensors`` table.
        """
        if name not in self._sensor_cache:
            cur.execute(
                "INSERT OR IGNORE INTO sensors (sensor_name) VALUES (?)",
                (name,),
            )
            cur.execute(
                "SELECT sensor_id FROM sensors WHERE sensor_name = ?",
                (name,),
            )
            self._sensor_cache[name] = cur.fetchone()[0]
        return self._sensor_cache[name]

    def _ensure_barn(self, cur: sqlite3.Cursor, name: str) -> int:
        """Get or create a ``barn_id`` for a barn name.

        Args:
            cur: Active database cursor.
            name: Barn name string, e.g. ``"NewBridge"``.

        Returns:
            The integer ``barn_id`` from the ``barns`` table.
        """
        if name not in self._barn_cache:
            cur.execute(
                "INSERT OR IGNORE INTO barns (barn_name) VALUES (?)",
                (name,),
            )
            cur.execute(
                "SELECT barn_id FROM barns WHERE barn_name = ?",
                (name,),
            )
            self._barn_cache[name] = cur.fetchone()[0]
        return self._barn_cache[name]

    def _ensure_source_file(
        self, cur: sqlite3.Cursor, filename: str, folder: str,
    ) -> int:
        """Get or create a ``file_id`` for a source CSV file.

        Every row in every fact table carries a ``file_id`` foreign
        key so that any datum can be traced back to the original CSV
        it came from.

        Args:
            cur: Active database cursor.
            filename: Basename of the CSV file.
            folder: Name of the containing folder (or ``"(standalone)"``
                for root-level files).

        Returns:
            The integer ``file_id`` from the ``source_files`` table.
        """
        key = (filename, folder)
        if key not in self._file_cache:
            cur.execute(
                "INSERT OR IGNORE INTO source_files (filename, folder) "
                "VALUES (?, ?)",
                (filename, folder),
            )
            cur.execute(
                "SELECT file_id FROM source_files "
                "WHERE filename = ? AND folder = ?",
                (filename, folder),
            )
            self._file_cache[key] = cur.fetchone()[0]
        return self._file_cache[key]

    # ─────────────────────────────────────────────────────────────
    #  « table creation »
    # ─────────────────────────────────────────────────────────────

    def _create_table_from_sample(
        self,
        con: sqlite3.Connection,
        table_name: str,
        sample_file: Path,
        entity_source: str | None,
        entity_col: str | None = None,
    ) -> list[str]:
        """Create a fact table by inferring schema from a sample CSV.

        On first call for a given *table_name*, reads up to 500 rows
        from *sample_file* to infer column types, then executes the
        ``CREATE TABLE`` statement.  Subsequent calls for the same
        table skip creation and only return the column-name list.

        Args:
            con: Active database connection.
            table_name: Target SQLite table name.
            sample_file: Path to the first CSV file in the folder,
                used for schema inference.
            entity_source: Entity-extraction method (see
                :data:`FOLDER_CONFIG`).
            entity_col: CSV column holding the entity ID, if any.

        Returns:
            Ordered list of sanitised CSV column names.
        """
        if table_name in self._created_tables:
            headers, _ = read_csv_sample(sample_file, n_rows=10)
            return [sanitise_column_name(h) for h in headers]

        # Use 500 rows for type inference to handle sparse columns
        headers, rows = read_csv_sample(sample_file, n_rows=500)
        sql, col_names = build_create_table_sql(
            table_name, headers, rows, entity_source, entity_col,
        )
        self.log.info("Creating table: %s", table_name)
        self.log.debug("SQL:\n%s", sql)
        con.executescript(sql)
        self._created_tables.add(table_name)
        return col_names

    # ─────────────────────────────────────────────────────────────
    #  « row insertion »
    # ─────────────────────────────────────────────────────────────

    def _ingest_csv(
        self,
        con: sqlite3.Connection,
        filepath: Path,
        table_name: str,
        col_names: list[str],
        entity_source: str | None,
        entity_col: str | None,
        folder_name: str,
    ) -> int:
        """Read a single CSV file and insert all rows into the database.

        Handles entity extraction (from filename or from a column in
        the data), registers entities in dimension tables, maps CSV
        columns to the target table schema, and streams rows in
        batches of :attr:`chunk_size`.

        Args:
            con: Active database connection.
            filepath: Path to the CSV file.
            table_name: Target SQLite table name.
            col_names: Ordered sanitised column names (from
                :meth:`_create_table_from_sample`).
            entity_source: Entity-extraction method.
            entity_col: CSV column holding the entity ID, if any.
            folder_name: Folder name for source-file provenance.

        Returns:
            Number of rows inserted.
        """
        cur = con.cursor()
        file_id = self._ensure_source_file(cur, filepath.name, folder_name)

        # ── extract entity from filename if needed ───────────────
        entity_id_from_file: str | None = None
        if entity_source in (
            "filename_animal", "filename_sensor", "filename_barn",
        ):
            entity_id_from_file = extract_entity_id_from_filename(
                filepath.name,
            )

        # ── pre-register entity ──────────────────────────────────
        sensor_id: int | None = None
        barn_id: int | None = None
        animal_id: int | None = None
        if entity_source == "filename_animal":
            animal_id = int(entity_id_from_file)
            self._ensure_animal(cur, animal_id)
        elif entity_source == "filename_sensor":
            sensor_id = self._ensure_sensor(cur, entity_id_from_file)
        elif entity_source == "filename_barn":
            barn_id = self._ensure_barn(cur, entity_id_from_file)

        # ── build INSERT column list ─────────────────────────────
        insert_cols: list[str] = []
        if entity_source == "filename_animal":
            insert_cols.append("animal_id")
        elif entity_source == "filename_sensor":
            insert_cols.append("sensor_id")
        elif entity_source == "filename_barn":
            insert_cols.append("barn_id")

        # Map CSV column index → SQL column name
        csv_col_map: list[tuple[int, str]] = []
        for i, cn in enumerate(col_names):
            if entity_col and cn == sanitise_column_name(entity_col):
                if entity_source in ("column_animal_id", "column_cow"):
                    sql_name = "animal_id"
                else:
                    sql_name = cn
            else:
                sql_name = cn

            if sql_name == "animal_id" and entity_source == "filename_animal":
                continue  # already prepended above

            insert_cols.append(sql_name)
            csv_col_map.append((i, sql_name))

        insert_cols.append("file_id")
        placeholders = ", ".join(["?"] * len(insert_cols))
        quoted_cols = ", ".join(f'"{c}"' for c in insert_cols)
        insert_sql = (
            f'INSERT INTO "{table_name}" ({quoted_cols}) '
            f"VALUES ({placeholders})"
        )

        # ── stream rows ──────────────────────────────────────────
        row_count = 0
        batch: list[tuple] = []

        with open(filepath, "r", newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            next(reader)  # skip header

            for row in reader:
                values: list = []

                # Prepend entity from filename
                if entity_source == "filename_animal":
                    values.append(animal_id)
                elif entity_source == "filename_sensor":
                    values.append(sensor_id)
                elif entity_source == "filename_barn":
                    values.append(barn_id)

                # Map CSV values
                for csv_idx, sql_name in csv_col_map:
                    if csv_idx < len(row):
                        val = row[csv_idx].strip()
                        # Register animals found in data columns
                        if sql_name == "animal_id" and val:
                            try:
                                aid = int(val)
                                self._ensure_animal(cur, aid)
                                values.append(aid)
                            except ValueError:
                                values.append(val if val else None)
                        # Also register animals in 'cow' columns
                        elif sql_name == "cow" and val:
                            try:
                                aid = int(val)
                                self._ensure_animal(cur, aid)
                                values.append(aid)
                            except ValueError:
                                values.append(val if val else None)
                        else:
                            values.append(val if val else None)
                    else:
                        values.append(None)

                values.append(file_id)
                batch.append(tuple(values))
                row_count += 1

                if len(batch) >= self.chunk_size:
                    cur.executemany(insert_sql, batch)
                    batch.clear()

            if batch:
                cur.executemany(insert_sql, batch)
                batch.clear()

        return row_count

    # ─────────────────────────────────────────────────────────────
    #  « index creation »
    # ─────────────────────────────────────────────────────────────

    def _create_indexes(self, con: sqlite3.Connection) -> None:
        """Create composite performance indexes on all fact tables.

        Indexes are created *after* bulk insertion to avoid the
        overhead of maintaining them during writes.  The primary
        pattern indexed is ``(entity_id, timestamp)`` for time-range
        queries filtered by animal, sensor, or barn.

        Args:
            con: Active database connection.
        """
        self.log.info("Creating indexes …")
        index_defs: list[tuple[str, str, list[str]]] = [
            # Animal-timestamped tables
            ("idx_bcs_animal_ts", "bcs",
             ["animal_id", "timestamp"]),
            ("idx_gouna_animal_ts", "gouna",
             ["animal_id", "timestamp"]),
            ("idx_herdeplus_animal_ts", "herdeplus",
             ["animal_id", "timestamp"]),
            ("idx_smaxtec_derived_animal_ts", "smaxtec_derived",
             ["animal_id", "timestamp"]),
            ("idx_smaxtec_events_animal_ts", "smaxtec_events",
             ["animal_id", "timestamp"]),
            ("idx_smaxtec_water_animal_ts", "smaxtec_water_intake",
             ["animal_id", "timestamp"]),
            # Allocations
            ("idx_alloc_animal", "allocations",
             ["animal_id"]),
            ("idx_alloc_enter", "allocations",
             ["datetime_enter"]),
            # Diseases
            ("idx_disease_animal", "diseases",
             ["animal_id"]),
            ("idx_disease_start", "diseases",
             ["disease_first_day"]),
            # Sensor / barn tables
            ("idx_lorawan_sensor_ts", "lorawan",
             ["sensor_id", "timestamp"]),
            ("idx_smaxtec_barns_barn_ts", "smaxtec_barns",
             ["barn_id", "timestamp"]),
            # Weather
            ("idx_dwd_dt", "dwd_weather",
             ["dt"]),
            ("idx_hobo_ts", "hobo_weather",
             ["datetime"]),
            # Source-file lookups
            ("idx_source_folder", "source_files",
             ["folder"]),
        ]

        cur = con.cursor()
        for idx_name, table, cols in index_defs:
            if table in self._created_tables:
                quoted = ", ".join(f'"{c}"' for c in cols)
                sql = (
                    f'CREATE INDEX IF NOT EXISTS "{idx_name}" '
                    f'ON "{table}" ({quoted})'
                )
                self.log.debug("  %s", sql)
                cur.execute(sql)
        con.commit()

    # ─────────────────────────────────────────────────────────────
    #  « main orchestrator »
    # ─────────────────────────────────────────────────────────────

    def run(self) -> None:
        """Execute the full ingestion pipeline.

        This is the main entry point.  It:

        1. Creates dimension tables (``animals``, ``sensors``,
           ``barns``, ``source_files``).
        2. Iterates over :data:`FOLDER_CONFIG`, creating each fact
           table from the first CSV's schema and then streaming all
           CSV files in that folder.
        3. Processes :data:`STANDALONE_FILES` the same way.
        4. Creates composite indexes on all fact tables.
        5. Prints a summary of rows inserted, database size, and
           elapsed time.

        Raises:
            sqlite3.Error: On unrecoverable database errors.
            Individual CSV failures are logged and skipped.
        """
        t0 = time.time()

        # ── banner ───────────────────────────────────────────────
        self.log.info("=" * 60)
        self.log.info("  COW DB INGESTION")
        self.log.info("  Root:  %s", self.root_dir)
        self.log.info("  DB:    %s", self.db_path)
        self.log.info("  Chunk: %s rows", f"{self.chunk_size:,}")
        if self.test_n is not None:
            self.log.info(
                "  Mode:  TEST (first %d files per folder)", self.test_n,
            )
        self.log.info("=" * 60)

        # ── open database with performance PRAGMAs ───────────────
        con = sqlite3.connect(str(self.db_path))
        con.execute("PRAGMA journal_mode = WAL")
        con.execute("PRAGMA synchronous = NORMAL")
        con.execute("PRAGMA cache_size = -512000")   # 512 MB
        con.execute("PRAGMA temp_store = MEMORY")

        con.executescript(DIMENSION_TABLES_SQL)
        con.commit()

        total_rows = 0

        # ── process folders ──────────────────────────────────────
        for folder_name, cfg in FOLDER_CONFIG.items():
            folder_path = self.root_dir / folder_name
            if not folder_path.is_dir():
                self.log.warning(
                    "Folder not found, skipping: %s", folder_path,
                )
                continue

            csv_files = sorted(folder_path.glob(cfg["glob"]))
            if not csv_files:
                self.log.warning("No CSV files in %s", folder_path)
                continue

            total_in_folder = len(csv_files)
            if self.test_n is not None:
                csv_files = csv_files[: self.test_n]

            self.log.info("─" * 50)
            if self.test_n is not None and total_in_folder > len(csv_files):
                self.log.info(
                    "Processing: %s  (%d / %d files — test mode)",
                    folder_name, len(csv_files), total_in_folder,
                )
            else:
                self.log.info(
                    "Processing: %s  (%d files)",
                    folder_name, len(csv_files),
                )

            col_names = self._create_table_from_sample(
                con, cfg["table"], csv_files[0],
                cfg["entity_source"], cfg.get("entity_col"),
            )
            con.commit()

            folder_rows = 0
            pbar = tqdm(
                csv_files, desc=f"  {cfg['table']}", unit="file",
            )
            for csv_file in pbar:
                try:
                    n = self._ingest_csv(
                        con, csv_file, cfg["table"], col_names,
                        cfg["entity_source"], cfg.get("entity_col"),
                        folder_name,
                    )
                    folder_rows += n
                    if HAS_TQDM:
                        pbar.set_postfix_str(f"{folder_rows:,} rows")
                except Exception as exc:
                    self.log.error(
                        "FAILED on %s: %s", csv_file.name, exc,
                    )
                    if self.verbose:
                        import traceback
                        traceback.print_exc()

            con.commit()
            total_rows += folder_rows
            self.log.info(
                "  → %s total rows in %s",
                f"{folder_rows:,}", cfg["table"],
            )

        # ── process standalone files ─────────────────────────────
        for filename, cfg in STANDALONE_FILES.items():
            filepath = self.root_dir / filename
            if not filepath.is_file():
                self.log.warning(
                    "Standalone file not found: %s", filepath,
                )
                continue

            self.log.info("─" * 50)
            self.log.info("Processing standalone: %s", filename)

            col_names = self._create_table_from_sample(
                con, cfg["table"], filepath,
                cfg["entity_source"], cfg.get("entity_col"),
            )
            con.commit()

            n = self._ingest_csv(
                con, filepath, cfg["table"], col_names,
                cfg["entity_source"], cfg.get("entity_col"),
                "(standalone)",
            )
            con.commit()
            total_rows += n
            self.log.info(
                "  → %s rows in %s", f"{n:,}", cfg["table"],
            )

        # ── indexes ──────────────────────────────────────────────
        self._create_indexes(con)

        # ── summary ──────────────────────────────────────────────
        elapsed = time.time() - t0
        db_size = self.db_path.stat().st_size / (1024 ** 3)

        self.log.info("=" * 60)
        self.log.info("  DONE")
        self.log.info("  Total rows inserted: %s", f"{total_rows:,}")
        self.log.info("  Database size:       %.2f GB", db_size)
        self.log.info(
            "  Elapsed:             %.1f s (%.1f min)",
            elapsed, elapsed / 60,
        )
        self.log.info(
            "  Animals registered:  %d", len(self._animal_cache),
        )
        self.log.info(
            "  Sensors registered:  %d", len(self._sensor_cache),
        )
        self.log.info(
            "  Barns registered:    %d", len(self._barn_cache),
        )
        self.log.info(
            "  Source files:        %d", len(self._file_cache),
        )
        self.log.info("=" * 60)

        cur = con.cursor()
        self.log.info("\n  Table row counts:")
        for table in sorted(self._created_tables):
            cur.execute(f'SELECT COUNT(*) FROM "{table}"')
            count = cur.fetchone()[0]
            self.log.info("    %-30s %12s", table, f"{count:,}")

        con.close()


# ═════════════════════════════════════════════════════════════════
#  « CLI entry point »
# ═════════════════════════════════════════════════════════════════

def main() -> None:
    """Parse command-line arguments and run the ingestion pipeline.

    This function is the console-script entry point registered in
    ``pyproject.toml`` as ``digimuh-ingest``.
    """
    parser = argparse.ArgumentParser(
        description="Ingest DigiMuh CSV exports into an SQLite database.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "root_dir", type=Path,
        help="Root directory containing all CSV folders",
    )
    parser.add_argument(
        "--db", type=Path, default=Path("cow.db"),
        help="Output SQLite database path (default: cow.db)",
    )
    parser.add_argument(
        "--chunk-size", type=int, default=50_000,
        help="Rows per INSERT batch (default: 50000)",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable debug logging",
    )
    parser.add_argument(
        "--test-n", type=int, default=None, metavar="N",
        help="Only ingest the first N files per folder (for testing)",
    )
    args = parser.parse_args()

    if not args.root_dir.is_dir():
        print(f"ERROR: Not a directory: {args.root_dir}", file=sys.stderr)
        sys.exit(1)

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s │ %(levelname)-7s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    ingester = CowDBIngester(
        root_dir=args.root_dir,
        db_path=args.db,
        chunk_size=args.chunk_size,
        verbose=args.verbose,
        test_n=args.test_n,
    )
    ingester.run()


if __name__ == "__main__":
    main()
