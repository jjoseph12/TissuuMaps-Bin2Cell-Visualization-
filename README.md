# TissUUmaps Plugins (UV-managed)

This repository now uses [uv](https://docs.astral.sh/uv/) to manage Python dependencies for the Bin2Cell Explorer and CircleSeg Overlay plugins.

## Quick start
1. Install uv (only once):  
   `curl -LsSf https://astral.sh/uv/install.sh | sh`
2. Ensure the uv shims are discoverable (add to `~/.zshrc` if desired):  
   `export PATH="$HOME/.local/bin:$PATH"`
3. From the project root (`cd ~/.tissuumaps`) install the locked environment:  
   `uv sync`
4. Launch TissUUmaps from the managed environment:  
   `uv run tissuumaps`  
   (Add `QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite"` if you want to silence Skia fallback logs.)
5. Optional: activate the virtual environment for an interactive dev shell:  
   `source .venv/bin/activate` → run commands → `deactivate`

After `uv sync` all runtime dependencies (`anndata`, `flask`, `numpy`, `scipy`, `scikit-image`, `matplotlib`) and dev tools (`pytest`, `mypy`) are available inside `.venv/`.

## Selecting dataset files
- Keep `.tiff`, `.npz`, and `.h5ad` assets wherever you like; the plugin works with absolute paths and no longer expects them in `~/.tissuumaps`.
- Either paste the path into each field or click the new **Browse…** buttons to pick the file with the system dialog (filters are pre-set for each file type).

## Common commands
- `uv run python plugins/Bin2CellExplorer.py` – execute scripts with project deps.
- `uv run tissuumaps` – start the desktop app without manually activating `.venv/`.
- `uv run pytest` – run the test suite (once tests are added).
- `uv run --group dev mkdocs serve` – live-preview the documentation.
- `uv add <package>` – add new runtime dependencies.
- `uv add --group dev <package>` – add tools needed only for development.
- `uv tree` – inspect the locked dependency graph.

## Notes
- `PySide6` is pinned to the 6.7 LTS line so that TissUUmaps launches without printing the Skia Graphite fallback noise introduced in newer Qt builds. If you prefer the latest Qt features, drop the `<6.8` pin in `pyproject.toml` (expect the console messages to return).
- Gene overlays support custom colours: choose `Solid` for a single hex colour or stay on `Gradient` and tweak the gradient colour picker to bias the ramp.
- The lockfile `uv.lock` tracks exact versions for reproducible installs.
- The `.venv/` directory is ignored by git; remove it with `rm -rf .venv` if you need a fresh environment.
- Existing plugin assets remain under `plugins/`; no import paths changed, so they continue to work with TissUUmaps.
