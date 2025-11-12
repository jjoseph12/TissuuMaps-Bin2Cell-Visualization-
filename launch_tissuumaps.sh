#!/usr/bin/env bash
# TissUUmaps headless launcher (cluster-friendly).
# Usage: ~/.tissuumaps/launch_tissuumaps.sh [--quiet] [--no-modules] [SLIDE_DIR]
# Modes:
#   - Default (verbose, modules ON)
#   - --quiet (silences logs; prints only SSH tunnel tip)
#   - --no-modules (skip environment modules; can combine with --quiet)
# Env (optional): TISSUUMAPS_HOST, TISSUUMAPS_PORT, TISSUUMAPS_VENV, TISSUUMAPS_USE_MODULES
# Examples:
#   ~/.tissuumaps/launch_tissuumaps.sh --quiet /path/to/slides
#   ~/.tissuumaps/launch_tissuumaps.sh --no-modules /path/to/slides
set -euo pipefail

# Make sure `module` is available even in non-interactive shells.
if ! command -v module >/dev/null 2>&1; then
    # shellcheck source=/dev/null
    source /etc/profile
fi

HOST="${TISSUUMAPS_HOST:-0.0.0.0}"
PORT="${TISSUUMAPS_PORT:-5678}"
QUIET="${TISSUUMAPS_QUIET:-0}"
SLIDE_DIR="$PWD"

# Simple logger that respects quiet mode
say() {
    if [[ "$QUIET" != "1" ]]; then
        echo "$@"
    fi
}

# Always print (used for important tips even in quiet mode)
note() {
    echo "$@"
}

# Parse arguments: optional -q/--quiet and optional slide directory
while [[ $# -gt 0 ]]; do
    case "$1" in
        -q|--quiet)
            QUIET=1
            shift
            ;;
        --no-modules)
            USE_MODULES=0
            shift
            ;;
        *)
            SLIDE_DIR="$1"
            shift
            ;;
    esac
done

VENV_PATH="${TISSUUMAPS_VENV:-$HOME/.tissuumaps/.venv}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG="$SCRIPT_DIR/tissuumaps_server.cfg"
CONFIG_FILE="${TISSUUMAPS_CONF:-$DEFAULT_CONFIG}"
DEFAULT_MODULES=("HDF5/1.14.3-gompi-2023b" "vips/8.17.0")

USE_MODULES="${TISSUUMAPS_USE_MODULES:-1}"

if [[ ! -d "$SLIDE_DIR" ]]; then
    echo "Slide directory '$SLIDE_DIR' does not exist" >&2
    exit 1
fi

if [[ ! -d "$VENV_PATH" ]]; then
    echo "Virtualenv '$VENV_PATH' not found. Set TISSUUMAPS_VENV if it lives elsewhere." >&2
    exit 1
fi

if [[ "$USE_MODULES" == "1" ]]; then
    if ! command -v module >/dev/null 2>&1; then
        echo "Modules requested (TISSUUMAPS_USE_MODULES=1) but 'module' command is unavailable." >&2
        exit 1
    fi
    if [[ "$QUIET" == "1" ]]; then
        module purge >/dev/null 2>&1 || true
    else
        module purge
    fi
    if [[ -n "${TISSUUMAPS_MODULES:-}" ]]; then
        IFS=' ' read -r -a REQUESTED_MODULES <<< "${TISSUUMAPS_MODULES}"
    else
        REQUESTED_MODULES=("${DEFAULT_MODULES[@]}")
    fi
    for mod in "${REQUESTED_MODULES[@]}"; do
        [[ -n "$mod" ]] || continue
        if [[ "$QUIET" == "1" ]]; then
            module load "$mod" >/dev/null 2>&1 || true
        else
            module load "$mod"
        fi
    done
else
    say "Skipping environment module loads."
fi

source "$VENV_PATH/bin/activate"

h5py_libdir=""
for candidate in \
    "$VENV_PATH/lib/python3.10/site-packages/h5py/.libs" \
    "$VENV_PATH/lib/python3.10/site-packages/h5py.libs"
do
    say "Checking for h5py libs at $candidate"
    if [[ -d "$candidate" ]]; then
        h5py_libdir="$candidate"
        say "Detected h5py private libs at $h5py_libdir"
        break
    fi
done

if [[ -d "$h5py_libdir" ]]; then
    LD_LIBRARY_PATH="$h5py_libdir:${LD_LIBRARY_PATH:-}"
fi

# If a module populated EBROOTHDF5 make sure its lib dir is preferred as well.
if [[ -n "${EBROOTHDF5:-}" && -d "${EBROOTHDF5}/lib" ]]; then
    LD_LIBRARY_PATH="${EBROOTHDF5}/lib:${LD_LIBRARY_PATH:-}"
    for lib in \
        "${EBROOTHDF5}/lib/libhdf5.so.310" \
        "${EBROOTHDF5}/lib/libhdf5_hl.so.310"
    do
        if [[ -f "$lib" ]]; then
            LD_PRELOAD="$lib:${LD_PRELOAD:-}"
        fi
    done
fi

export LD_LIBRARY_PATH LD_PRELOAD

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "config file '$CONFIG_FILE' not found" >&2
    exit 1
fi

say "Starting TissUUmaps server on $HOST:$PORT"
say "Serving slide directory: $SLIDE_DIR"
say "Using config: $CONFIG_FILE"
say "Using virtualenv: $VENV_PATH"
say "Press Ctrl+C to stop."
say
note "Tip: from your laptop run 'ssh -L ${PORT}:localhost:${PORT} <cluster-host>' and open http://localhost:${PORT}/"
say "LD_LIBRARY_PATH resolved to: ${LD_LIBRARY_PATH:-<empty>}"

TISSUUMAPS_CONF="$CONFIG_FILE" python -m tissuumaps -l "$HOST" -p "$PORT" "$SLIDE_DIR"
