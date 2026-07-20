#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────
#  DigiMuh — database sanity check
#  « size, indexes, and a row count on the largest fact table »
# ─────────────────────────────────────────────────────────────────
#  The database path is resolved through digimuh.config, so this
#  script runs unmodified on any machine that has a config.yaml.
# ─────────────────────────────────────────────────────────────────
set -euo pipefail

python3 -c "
import sqlite3
from digimuh.config import load_config

db = load_config().database
print(f'DB: {db}')
print(f'DB size: {db.stat().st_size / 1e9:.1f} GB')

con = sqlite3.connect(db)
cur = con.cursor()

# Check indexes
cur.execute(\"SELECT name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'\")
indexes = [r[0] for r in cur.fetchall()]
print(f'\nIndexes found: {len(indexes)}')
for name in sorted(indexes):
    print(f'  {name}')

# Quick row count on the big table
cur.execute('SELECT COUNT(*) FROM smaxtec_derived')
print(f'\nsmaxtec_derived: {cur.fetchone()[0]:,} rows')
con.close()
"
