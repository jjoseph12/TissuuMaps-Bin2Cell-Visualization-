#!/usr/bin/env bash
set -euo pipefail

# Setup TissUUmaps with uv and install Bin2Cell plugin files, then optionally launch.
# Usage: ./setup_uv_bin2cell.sh <plugin_src_dir> [--no-launch] [--quiet-skia]
# Notes:
#   - <plugin_src_dir> must include CellExplorer.{js,py,yml}. Optional assets dir.
#   - Ensures ~/.local/bin on PATH, installs uv if missing, runs 'uv sync', copies (copies to the Right Path), launches unless --no-launch.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"
PLUGINS_DIR="$PROJECT_ROOT/plugins"

print_usage() {
    echo "Usage: $0 <plugin_src_dir> [--no-launch] [--quiet-skia]"
}

if [[ $# -lt 1 ]]; then
    print_usage
    exit 1
fi

PLUGIN_SRC=""
NO_LAUNCH=0
QUIET_SKIA=0

for arg in "$@"; do
    case "$arg" in
        --no-launch) NO_LAUNCH=1 ;;
        --quiet-skia) QUIET_SKIA=1 ;;
        *)
            if [[ -z "$PLUGIN_SRC" ]]; then
                PLUGIN_SRC="$arg"
            else
                echo "Unexpected argument: $arg" >&2
                print_usage
                exit 1
            fi
            ;;
    esac
done

if [[ -z "$PLUGIN_SRC" ]]; then
    echo "Error: <plugin_src_dir> is required." >&2
    print_usage
    exit 1
fi

if [[ ! -d "$PLUGIN_SRC" ]]; then
    echo "Error: plugin source directory '$PLUGIN_SRC' does not exist." >&2
    exit 1
fi

# Ensure ~/.local/bin on PATH for this shell (uv installer puts the shim there)
export PATH="$HOME/.local/bin:$PATH"

# Install uv if not present
if ! command -v uv >/dev/null 2>&1; then
    echo "uv not found. Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    # Ensure the new shim is visible in this shell
    export PATH="$HOME/.local/bin:$PATH"
    if ! command -v uv >/dev/null 2>&1; then
        echo "uv installation appears to have failed (uv not on PATH)." >&2
        exit 1
    fi
else
    echo "uv found at: $(command -v uv)"
fi

# Move to project root and sync the environment
cd "$PROJECT_ROOT"
echo "Syncing environment from lockfile with 'uv sync'..."
uv sync

# Prepare plugins directory
mkdir -p "$PLUGINS_DIR"

echo "Copying Cell Explorer plugin files from: $PLUGIN_SRC"

pick_file() {
    local base="$1"
    if [[ -f "$PLUGIN_SRC/$base" ]]; then
        echo "$PLUGIN_SRC/$base"
        return 0
    fi
    # legacy names
    local legacy="Bin2CellExplorer.${base##*.}"
    if [[ -f "$PLUGIN_SRC/$legacy" ]]; then
        echo "$PLUGIN_SRC/$legacy"
        return 0
    fi
    return 1
}

js_src="$(pick_file 'CellExplorer.js')" || { echo "Missing CellExplorer.js (or Bin2CellExplorer.js) in $PLUGIN_SRC" >&2; exit 1; }
py_src="$(pick_file 'CellExplorer.py')" || { echo "Missing CellExplorer.py (or Bin2CellExplorer.py) in $PLUGIN_SRC" >&2; exit 1; }
yml_src="$(pick_file 'CellExplorer.yml')" || { echo "Missing CellExplorer.yml (or Bin2CellExplorer.yml) in $PLUGIN_SRC" >&2; exit 1; }

# Copy the three main files with new names
cp -f "$js_src" "$PLUGINS_DIR/CellExplorer.js"
cp -f "$py_src" "$PLUGINS_DIR/CellExplorer.py"
cp -f "$yml_src" "$PLUGINS_DIR/CellExplorer.yml"

# If there is an asset directory, sync it as well (prefer new name, fallback legacy)
assets_src=""
if [[ -d "$PLUGIN_SRC/CellExplorer" ]]; then
    assets_src="$PLUGIN_SRC/CellExplorer"
elif [[ -d "$PLUGIN_SRC/Bin2CellExplorer" ]]; then
    assets_src="$PLUGIN_SRC/Bin2CellExplorer"
fi
if [[ -n "$assets_src" ]]; then
    mkdir -p "$PLUGINS_DIR/CellExplorer"
    rsync -a --delete "$assets_src/" "$PLUGINS_DIR/CellExplorer/"
fi

echo "Plugin files installed into: $PLUGINS_DIR"

if [[ "$NO_LAUNCH" -eq 1 ]]; then
    echo "Setup complete. Skipping launch (--no-launch set)."
    echo "To launch later:"
    if [[ "$QUIET_SKIA" -eq 1 ]]; then
        echo "  QTWEBENGINE_CHROMIUM_FLAGS=\"--disable-features=SkiaGraphite\" uv run tissuumaps"
    else
        echo "  uv run tissuumaps"
    fi
    exit 0
fi

echo "Launching TissUUmaps with uv..."
if [[ "$QUIET_SKIA" -eq 1 ]]; then
    QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite" uv run tissuumaps
else
    uv run tissuumaps
fi



