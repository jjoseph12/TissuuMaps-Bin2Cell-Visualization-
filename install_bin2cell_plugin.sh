#!/usr/bin/env bash
set -euo pipefail

# Install Cell Explorer plugin into ~/.tissuumaps/plugins from GitHub.
# Usage: ./install_bin2cell_plugin.sh [--repo URL] [--branch BRANCH] [--target DIR] [--force] [--persist-path]
# Notes:
#   - Copies CellExplorer.{js,py,yml} (supports older Bin2CellExplorer names(need to update this)) + optional assets dir.
#   - Requires git. Does not install uv or run 'uv sync'.
# Next steps (manual):
#   export PATH="$HOME/.local/bin:$PATH"
#   cd ~/.tissuumaps && uv sync && uv run tissuumaps

REPO_URL="https://github.com/jjoseph12/TissuuMaps-Cell-Explorer-Visualization.git"
BRANCH="main"
TARGET_DIR="$HOME/.tissuumaps/plugins"
FORCE_OVERWRITE=0
PERSIST_PATH=0

print_usage() {
    echo "Usage: $0 [--repo URL] [--branch BRANCH] [--target DIR] [--force] [--persist-path]"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo)
            REPO_URL="${2:-}"
            shift 2
            ;;
        --branch)
            BRANCH="${2:-}"
            shift 2
            ;;
        --target)
            TARGET_DIR="${2:-}"
            shift 2
            ;;
        --force)
            FORCE_OVERWRITE=1
            shift
            ;;
        --persist-path)
            PERSIST_PATH=1
            shift
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            print_usage
            exit 1
            ;;
    esac
done

# Ensure target directory exists
mkdir -p "$TARGET_DIR"

tmpdir="$(mktemp -d)"
cleanup() {
    rm -rf "$tmpdir" || true
}
trap cleanup EXIT

echo "Cloning repo:"
echo "  URL:    $REPO_URL"
echo "  Branch: $BRANCH"
git --version >/dev/null 2>&1 || { echo "Error: git is required." >&2; exit 1; }
git clone --depth 1 --branch "$BRANCH" "$REPO_URL" "$tmpdir/repo"

repo_dir="$tmpdir/repo"

find_one() {
    local name="$1"
    # Search up to reasonable depth to accommodate common layouts
    local path
    path="$(find "$repo_dir" -type f -name "$name" -print -quit)"
    if [[ -z "$path" ]]; then
        return 1
    fi
    echo "$path"
}

find_dir_one() {
    local name="$1"
    local path
    path="$(find "$repo_dir" -type d -name "$name" -print -quit)"
    if [[ -z "$path" ]]; then
        return 1
    fi
    echo "$path"
}

echo "Locating required plugin files..."
js_src="$(find_one 'CellExplorer.js' || true)"; [[ -n "$js_src" ]] || js_src="$(find_one 'Bin2CellExplorer.js' || true)"
py_src="$(find_one 'CellExplorer.py' || true)"; [[ -n "$py_src" ]] || py_src="$(find_one 'Bin2CellExplorer.py' || true)"
yml_src="$(find_one 'CellExplorer.yml' || true)"; [[ -n "$yml_src" ]] || yml_src="$(find_one 'Bin2CellExplorer.yml' || true)"
assets_src="$(find_dir_one 'CellExplorer' || true)"; [[ -n "$assets_src" ]] || assets_src="$(find_dir_one 'Bin2CellExplorer' || true)"

missing=()
[[ -n "$js_src" ]] || missing+=("CellExplorer.js (or Bin2CellExplorer.js)")
[[ -n "$py_src" ]] || missing+=("CellExplorer.py (or Bin2CellExplorer.py)")
[[ -n "$yml_src" ]] || missing+=("CellExplorer.yml (or Bin2CellExplorer.yml)")

if (( ${#missing[@]} > 0 )); then
    echo "Error: missing required files in the repository: ${missing[*]}" >&2
    exit 1
fi

copy_file() {
    local src="$1"
    local dest_dir="$2"
    local dest_name="$3"
    local dest="$dest_dir/$dest_name"
    if [[ -e "$dest" && "$FORCE_OVERWRITE" -ne 1 ]]; then
        echo "Skipping existing $dest_name (use --force to overwrite)"
    else
        cp -f "$src" "$dest"
        echo "Installed $dest_name -> $dest_dir"
    fi
}

echo "Installing plugin files into: $TARGET_DIR"
copy_file "$js_src" "$TARGET_DIR" "CellExplorer.js"
copy_file "$py_src" "$TARGET_DIR" "CellExplorer.py"
copy_file "$yml_src" "$TARGET_DIR" "CellExplorer.yml"

if [[ -n "$assets_src" && -d "$assets_src" ]]; then
    dest_assets="$TARGET_DIR/CellExplorer"
    mkdir -p "$dest_assets"
    # rsync if available for clean sync, otherwise fallback to cp -r
    if command -v rsync >/dev/null 2>&1; then
        rsync -a --delete "$assets_src/" "$dest_assets/"
    else
        rm -rf "$dest_assets"/*
        cp -r "$assets_src/"* "$dest_assets/" 2>/dev/null || true
    fi
    echo "Installed asset directory -> $dest_assets"
else
    echo "No asset directory 'CellExplorer/' or 'Bin2CellExplorer/' found in repo; skipping."
fi

if [[ "$PERSIST_PATH" -eq 1 ]]; then
    # Persist PATH export to ~/.bashrc (shell shown is /bin/bash)
    line='export PATH="$HOME/.local/bin:$PATH"'
    if ! grep -Fq "$line" "$HOME/.bashrc" 2>/dev/null; then
        echo "$line" >> "$HOME/.bashrc"
        echo "Appended PATH export to ~/.bashrc"
        echo "Run: source ~/.bashrc"
    else
        echo "PATH export already present in ~/.bashrc"
    fi
fi

echo
echo "Done. Next steps (run manually):"
echo "  1) Ensure uv is installed:   curl -LsSf https://astral.sh/uv/install.sh | sh"
echo '  2) Make shim available:      export PATH="$HOME/.local/bin:$PATH"'
echo "  3) From ~/.tissuumaps:       uv sync"
echo "  4) Launch:                   uv run tissuumaps"
echo '     Optional (quiet Skia):    QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite" uv run tissuumaps'


