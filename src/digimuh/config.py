#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  DigiMuh configuration loader                                ║
# ║  « machine-specific paths without polluting the repo »       ║
# ╚══════════════════════════════════════════════════════════════╝
"""Hierarchical configuration for DigiMuh analysis pipelines.

Priority (highest wins):

1. CLI arguments (always override everything)
2. ``.env`` in the project directory (quick per-project overrides)
3. ``~/.config/digimuh/config.yaml`` (machine-specific, never in repo)
4. Built-in defaults

Usage in an entry point::

    from digimuh.config import load_config

    cfg = load_config()
    # cfg.database   → Path to cow.db
    # cfg.output     → Path to results directory
    # cfg.tierauswahl → Path to Tierauswahl.xlsx
    # cfg.n_jobs     → Number of parallel workers

CLI arguments still work and override everything::

    digimuh-extract --db /other/cow.db --out /tmp/test

Setup a new machine::

    digimuh-config

This creates ``~/.config/digimuh/config.yaml`` interactively.
"""

from __future__ import annotations

import argparse
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

log = logging.getLogger("digimuh.config")

# ─────────────────────────────────────────────────────────────
#  « paths »
# ─────────────────────────────────────────────────────────────

XDG_CONFIG = Path(
    os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")
)
CONFIG_DIR = XDG_CONFIG / "digimuh"
CONFIG_FILE = CONFIG_DIR / "config.yaml"
ENV_FILE = Path(".env")


# ─────────────────────────────────────────────────────────────
#  « config dataclass »
# ─────────────────────────────────────────────────────────────

@dataclass
class DigiMuhConfig:
    """Resolved configuration for a DigiMuh pipeline run."""

    database: Path | None = None
    output: Path = field(default_factory=lambda: Path("results/broken_stick"))
    tierauswahl: Path | None = None
    n_jobs: int = 20
    smaxtec_drink_correction: bool = False

    # Sources (for logging)
    _sources: dict = field(default_factory=dict, repr=False)


# ─────────────────────────────────────────────────────────────
#  « loaders »
# ─────────────────────────────────────────────────────────────

def _load_yaml(path: Path) -> dict[str, Any]:
    """Load a YAML config file.  Returns empty dict if missing."""
    if not path.is_file():
        return {}
    try:
        import yaml
        with open(path) as f:
            data = yaml.safe_load(f) or {}
        log.debug("Loaded config from %s", path)
        return data
    except ImportError:
        # Fall back to simple key: value parsing if pyyaml not installed
        data = {}
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if ":" in line:
                    key, _, val = line.partition(":")
                    val = val.strip().strip("'\"")
                    if val:
                        data[key.strip()] = val
        log.debug("Loaded config from %s (simple parser)", path)
        return data
    except Exception as e:
        log.warning("Could not read %s: %s", path, e)
        return {}


def _load_env(path: Path) -> dict[str, str]:
    """Load a .env file.  Returns empty dict if missing."""
    if not path.is_file():
        return {}
    data = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, val = line.partition("=")
                val = val.strip().strip("'\"")
                data[key.strip()] = val
    log.debug("Loaded .env from %s", path)
    return data


# ─────────────────────────────────────────────────────────────
#  « resolver »
# ─────────────────────────────────────────────────────────────

# Mapping from config keys to env/yaml keys
_KEY_MAP = {
    "database":    ("DIGIMUH_DB", "database"),
    "output":      ("DIGIMUH_OUT", "output"),
    "tierauswahl": ("DIGIMUH_TIERAUSWAHL", "tierauswahl"),
    "n_jobs":      ("DIGIMUH_NJOBS", "n_jobs"),
    "smaxtec_drink_correction": ("DIGIMUH_SMAXTEC_DRINK", "smaxtec_drink_correction"),
}


def load_config(
    cli_args: argparse.Namespace | None = None,
    project_root: Path | None = None,
) -> DigiMuhConfig:
    """Load configuration with full priority chain.

    Args:
        cli_args: Parsed argparse namespace (from entry point).
        project_root: Project root for .env lookup.  Defaults to cwd.

    Returns:
        Resolved DigiMuhConfig.
    """
    cfg = DigiMuhConfig()
    sources: dict[str, str] = {}

    # Layer 1: XDG config (lowest priority)
    yaml_data = _load_yaml(CONFIG_FILE)

    # Layer 2: .env in project dir
    root = project_root or Path.cwd()
    env_data = _load_env(root / ".env")

    # Layer 3: OS environment variables
    os_env = dict(os.environ)

    # Resolve each key
    for attr, (env_key, yaml_key) in _KEY_MAP.items():
        # Priority: CLI > OS env > .env > yaml > default
        val = None
        source = "default"

        # YAML
        if yaml_key in yaml_data:
            val = yaml_data[yaml_key]
            source = str(CONFIG_FILE)

        # .env file
        if env_key in env_data:
            val = env_data[env_key]
            source = ".env"

        # OS environment
        if env_key in os_env:
            val = os_env[env_key]
            source = "environment"

        # CLI args (highest priority)
        if cli_args is not None:
            # Map config attr names to CLI arg names
            cli_map = {
                "database": "db",
                "output": "out",
                "tierauswahl": "tierauswahl",
                "n_jobs": "n_jobs",
                "smaxtec_drink_correction": "smaxtec_drink_correction",
            }
            cli_name = cli_map.get(attr, attr)
            cli_val = getattr(cli_args, cli_name, None)
            if cli_val is not None:
                val = cli_val
                source = "CLI"

        # Apply
        if val is not None:
            if attr in ("database", "output", "tierauswahl"):
                setattr(cfg, attr, Path(str(val)).expanduser())
            elif attr == "n_jobs":
                setattr(cfg, attr, int(val))
            elif attr == "smaxtec_drink_correction":
                if isinstance(val, bool):
                    setattr(cfg, attr, val)
                else:
                    setattr(cfg, attr, str(val).lower() in ("true", "1", "yes"))
            sources[attr] = source

    cfg._sources = sources
    return cfg


def print_config(cfg: DigiMuhConfig) -> None:
    """Log the resolved configuration with sources."""
    log.info("Configuration:")
    for attr in ["database", "output", "tierauswahl", "n_jobs", "smaxtec_drink_correction"]:
        val = getattr(cfg, attr)
        source = cfg._sources.get(attr, "default")
        log.info("  %-28s = %-40s [%s]", attr, val, source)


# ─────────────────────────────────────────────────────────────
#  « interactive setup (digimuh-config) »
# ─────────────────────────────────────────────────────────────

def setup_interactive() -> None:
    """Interactive setup for a new machine.  Creates config.yaml."""
    print("=" * 60)
    print("  DigiMuh configuration setup")
    print("=" * 60)
    print(f"\nConfig will be saved to: {CONFIG_FILE}\n")

    defaults = {
        "database": str(Path.home() / "cow.db"),
        "output": "results/broken_stick",
        "tierauswahl": "Tierauswahl_extended.xlsx",
        "n_jobs": str(max(1, os.cpu_count() - 4)) if os.cpu_count() else "20",
        "smaxtec_drink_correction": "false",
    }

    values = {}
    for key, default in defaults.items():
        answer = input(f"  {key} [{default}]: ").strip()
        values[key] = answer if answer else default

    # Validate database path
    db_path = Path(values["database"]).expanduser()
    if not db_path.is_file():
        print(f"\n  Warning: database not found at {db_path}")
        print("  (Config saved anyway; fix the path later.)")

    # Write YAML
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_FILE, "w") as f:
        f.write("# DigiMuh configuration\n")
        f.write(f"# Machine: {os.uname().nodename}\n")
        f.write(f"# Created: {__import__('datetime').datetime.now().isoformat()}\n\n")
        for key, val in values.items():
            f.write(f"{key}: {val}\n")

    print(f"\nConfig saved to {CONFIG_FILE}")
    print("You can edit it anytime with a text editor.")
    print("\nTest with:  digimuh-extract --help")


def main() -> None:
    """Entry point for digimuh-config."""
    setup_interactive()


if __name__ == "__main__":
    main()
