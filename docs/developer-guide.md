# Developer Guide

This document explains how to set up a development environment, key dependencies, and common maintenance tasks for the Bin2Cell Explorer plugin.

## Environment Setup

### 1. Install uv (one time)

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"   # add to ~/.zshrc for convenience
```

### 2. Sync the managed environment

From the project root (`~/.tissuumaps`):

```bash
uv sync
```

This creates/updates `.venv` with the locked runtime stack (anndata, numpy, scikit-image, matplotlib, tissuumaps, PySide6 6.7 LTS, etc.) plus dev tools (`pytest`, `mkdocs`, …).

To work in an activated shell:

```bash
source .venv/bin/activate
# ...run commands...
deactivate
```

### 3. Launch TissUUmaps through uv

```bash
uv run tissuumaps
```

Add `QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite"` if you want to suppress the Skia fallback message. The desktop app inherits the `.venv` environment without requiring conda.

### 4. Place the plugin files

Keep the plugin code inside the TissUUmaps configuration directory:

```
~/.tissuumaps/
└── plugins/
    ├── Bin2CellExplorer.py
    ├── Bin2CellExplorer.js
    ├── Bin2CellExplorer.yml
    └── Bin2CellExplorer/
        └── Bin2CellExplorer.log   # auto-generated
```

Dataset assets (`.tif/.tiff`, `.npz`, `.h5ad`) can live anywhere on disk. In the UI you can paste absolute paths or click the matching **Browse…** buttons to open a filtered file dialog for each requirement.

## Project Layout

| File | Purpose |
| ---- | ------- |
| `Bin2CellExplorer.py` | Flask plugin backend. Handles dataset loading, caching, overlay generation, exporting, and preset persistence. |
| `Bin2CellExplorer.js` | Front-end logic. Builds the UI, calls backend endpoints, renders overlays, draws tile grids, and minimap. |
| `Bin2CellExplorer.yml` | Plugin metadata for TissUUmaps. |
| `he.tiff` / `he.npz` / `P2_CRC_annotated.h5ad` | Sample data used during development and demos. Replace with user data as needed. |
| `Bin2CellExplorer.log` | Rolling log (one per plugin) capturing load/overlay operations. |

## Running MkDocs locally

```bash
uv run --group dev mkdocs serve
```

The docs are served at `http://127.0.0.1:8000`. Hit `Ctrl+C` to stop. For a static build use `uv run --group dev mkdocs build`.

## Debugging Tips

- Tail the log while testing: `tail -f ~/.tissuumaps/plugins/Bin2CellExplorer/Bin2CellExplorer.log`.
- If you change data paths, delete `~/.tissuumaps/plugins/Bin2CellExplorer/dataset_state.json` to force a full reload.
- When overlays are missing, confirm the tile ID exists and that the AnnData column/gene is valid (per the notebook).
- `scanpy.read_h5ad("P2_CRC_annotated.h5ad")` is a quick sanity check that the `.h5ad` can be loaded.

## Common Tasks

### Updating Dependencies

1. Add or bump dependencies in `pyproject.toml` (use `uv add <package>` for runtime, `uv add --group dev <package>` for tooling).
2. Run `uv lock` (if required) and `uv sync` to refresh `.venv`.
3. Document noteworthy changes here and, if needed, in `docs/design-rationale.md`.

### Adding New UI Controls

- Define the control in `Bin2CellExplorer.parameters`.
- Handle it in `inputTrigger`.
- Add the required state fields or backend parameters.
- Update documentation (`usage.md` or `architecture.md`) to keep feature parity.

### Extending Backend Endpoints

- Follow the existing pattern: sanitize inputs, wrap logic in `try/except`, return JSON responses via `_json_response`.
- Prefer raising Python exceptions (the `Plugin` wrapper converts them into structured error payloads).
- For heavy operations, cache results via `TileCache`.

## Release Checklist

1. Verify `Bin2CellExplorer.yml` version.
2. Run through the usage steps end-to-end with sample data.
3. Ensure `Bin2CellExplorer.log` is free of unexpected errors.
4. Regenerate documentation if anything changed (`mkdocs build`).
5. Package plugin files with any sample data provided to users.

With these steps a new developer can reproduce the environment, understand the code layout, and extend the plugin safely.
