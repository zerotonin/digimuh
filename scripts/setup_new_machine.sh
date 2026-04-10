#!/usr/bin/env bash
# ╔══════════════════════════════════════════════════════════════╗
# ║  Setup script for DigiMuh on a new machine                   ║
# ║  Run once after cloning the repo.                            ║
# ╚══════════════════════════════════════════════════════════════╝

set -e

echo "═══════════════════════════════════════════════════════"
echo "  DigiMuh — New machine setup"
echo "═══════════════════════════════════════════════════════"
echo ""

# 1. Install the package in editable mode
echo "[1/4] Installing digimuh in editable mode …"
pip install -e ".[analysis]" 2>/dev/null || pip install -e .
echo "  Done."
echo ""

# 2. Install optional dependencies
echo "[2/4] Installing optional dependencies …"
pip install openpyxl plotly kaleido pyyaml 2>/dev/null || true
echo "  Done."
echo ""

# 3. Run interactive config setup
echo "[3/4] Setting up machine-specific configuration …"
echo ""
digimuh-config
echo ""

# 4. Verify
echo "[4/4] Verifying installation …"
echo ""
echo "  Entry points:"
for cmd in digimuh-extract digimuh-stats digimuh-plots digimuh-config; do
    if command -v "$cmd" &>/dev/null; then
        echo "    ✓ $cmd"
    else
        echo "    ✗ $cmd (not found — try: pip install -e .)"
    fi
done

echo ""
echo "  Config file:"
CONFIG_FILE="$HOME/.config/digimuh/config.yaml"
if [ -f "$CONFIG_FILE" ]; then
    echo "    ✓ $CONFIG_FILE"
    echo ""
    echo "  Contents:"
    sed 's/^/    /' "$CONFIG_FILE"
else
    echo "    ✗ Not found"
fi

echo ""
echo "═══════════════════════════════════════════════════════"
echo "  Setup complete!"
echo ""
echo "  Run the pipeline:"
echo "    bash scripts/run_00_broken_stick_ana.sh"
echo ""
echo "  Or step by step:"
echo "    digimuh-extract        # DB → CSV"
echo "    digimuh-stats --data results/broken_stick"
echo "    digimuh-plots --data results/broken_stick"
echo ""
echo "  Edit config anytime:"
echo "    nano $CONFIG_FILE"
echo "═══════════════════════════════════════════════════════"
