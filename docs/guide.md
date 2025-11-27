# Usage & Developer Guide

This guide covers both end-user workflows and developer setup for the Bin2Cell Explorer plugin.

## Usage

### 1. Launch TissUUmaps with uv

```bash
cd ~/.tissuumaps
uv sync                              # only after dependency changes
QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite" uv run tissuumaps
```

The optional `QTWEBENGINE_CHROMIUM_FLAGS` silences the Skia Graphite → Ganesh fallback log. Omit it if you want the default logging.

### 2. Load the dataset once

Open **Plugins → Bin2Cell Explorer** and either paste the absolute paths into the three text fields or click the matching **Browse…** buttons to select each file with the system dialog (filters default to `.tif/.tiff`, `.npz`, and `.h5ad`).

Click **Load dataset**. The status line and log show `Dataset loaded successfully: … tiles`. The plugin caches the dataset so later overlay updates are instant and the selected paths persist between sessions.

> **Tip:** The initial load reads an 830 MB TIFF and a 4 GB `.h5ad`, so it can take ~20–30 seconds. Avoid clicking "Load dataset" repeatedly; the plugin will reuse the cached dataset automatically.

### 3. Orient yourself with the tile overview

Scroll a little in the plugin panel to the **Tile overview** minimap:

- Each tile is a rectangle labelled with its numeric ID.
- Clicking a rectangle selects that tile, updates the dropdown, highlights it in amber, and pans/zooms the main viewer directly to that region.

The minimap is larger now, so you can comfortably click tiles even on dense slides without needing an on-slide grid.

### 4. Choose your overlay

Select a tile via the minimap or dropdown, then set overlay parameters:

- **Gene overlay:** choose `Overlay type = gene`, enter a gene (e.g., `COL1A1`), pick color mode (`gradient`, `binary`, or `solid`), and click **Update overlay**. Use the colour picker below the dropdown: it adjusts the gradient ramp when `gradient` is selected, or sets the fill/stroke colour in `solid` mode.
- **Observation overlay:** choose `Overlay type = observation`, pick an observation column (e.g., `predicted_ClusterFull`, `predicted_ClusterMidway`, or `predicted_ClusterTop`), then use the now-populated **Category filter** dropdown if you want to isolate a single lineage. Leaving it on "All categories" renders every class with the `AnnData.uns['<column>_colors']` palette.

#### Bin2Cell parameters (what the four numbers mean)

- **Max bin distance (µm):** caps how far each nucleus can inflate. The plugin converts it to pixels via `distance_px = ceil(max_bin_distance * (bin_um / mpp))`, so doubling the distance or bin size doubles the apparent halo.
- **Microns per pixel (mpp):** physical resolution of the slide. Lower values mean each pixel maps to fewer microns, so the same bin distance spans more pixels.
- **Bin size (µm):** physical size of the Bin2Cell bin; paired with `mpp` to keep distance scaling consistent when you swap datasets.
- **Volume ratio:** only used when `Expand mode = volume_ratio`. Instead of using the fixed cap, it inflates each nucleus so the expanded area approximates `volume_ratio × (original area)`; the plugin converts that area delta to an average radial offset.

Tuning any of the first three numbers affects the others because the distance is ultimately applied in pixels. A typical workflow is to keep `mpp` from the slide metadata, set `bin_um` to the Bin2Cell bin size used upstream, and then tweak `max_bin_distance` (or `volume_ratio`) to control how much the halos bleed into neighbours.

Legend checkboxes let you toggle categories or gene overlays without re-rendering.

### 5. Optional layers and exports

- Toggle **Show expanded outlines** to visualize the Bin2Cell expansion (fills + outlines) or **Show nuclei outlines** to overlay the original segmentation contours. Turn both on if you want to compare halo vs. nucleus footprints simultaneously.
- When **Show expanded outlines** is unchecked the filled shapes follow the original nuclei contours; enable it if you want to visualize the Bin2Cell-expanded radius alongside (and as the basis for) the overlay.
- Click **Export overlay to GeoJSON** to capture the current overlay for downstream analysis.
- Save reusable settings via the **Presets** section.

### 6. Resetting or switching datasets

If you need to work with a different slide:

1. Update the file paths.
2. Click **Load dataset** once.
3. The plugin updates the cache (`dataset_state.json`) and refreshes the tile overview automatically.

That's it—once the dataset is loaded you can hop between tiles, switch genes or observation columns, and export overlays without repeating the heavy load step.

## Developer Setup

### Environment Setup

#### 1. Install uv (one time)

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"   # add to ~/.zshrc for convenience
```

#### 2. Sync the managed environment

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

#### 3. Launch TissUUmaps through uv

```bash
uv run tissuumaps
```

Add `QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite"` if you want to suppress the Skia fallback message. The desktop app inherits the `.venv` environment without requiring conda.

#### 4. Place the plugin files

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

### Project Layout

| File | Purpose |
| ---- | ------- |
| `Bin2CellExplorer.py` | Flask plugin backend. Handles dataset loading, caching, overlay generation, exporting, and preset persistence. |
| `Bin2CellExplorer.js` | Front-end logic. Builds the UI, calls backend endpoints, renders overlays, and drives the minimap. |
| `Bin2CellExplorer.yml` | Plugin metadata for TissUUmaps. |
| `he.tiff` / `he.npz` / `P2_CRC_annotated.h5ad` | Sample data used during development and demos. Replace with user data as needed. |
| `Bin2CellExplorer.log` | Rolling log (one per plugin) capturing load/overlay operations. |

### Running MkDocs locally

```bash
uv run --group dev mkdocs serve
```

The docs are served at `http://127.0.0.1:8000`. Hit `Ctrl+C` to stop. For a static build use `uv run --group dev mkdocs build`.

### Debugging Tips

- Tail the log while testing: `tail -f ~/.tissuumaps/plugins/Bin2CellExplorer/Bin2CellExplorer.log`.
- If you change data paths, delete `~/.tissuumaps/plugins/Bin2CellExplorer/dataset_state.json` to force a full reload.
- When overlays are missing, confirm the tile ID exists and that the AnnData column/gene is valid (per the notebook).
- `scanpy.read_h5ad("P2_CRC_annotated.h5ad")` is a quick sanity check that the `.h5ad` can be loaded.

### Common Tasks

#### Updating Dependencies

1. Add or bump dependencies in `pyproject.toml` (use `uv add <package>` for runtime, `uv add --group dev <package>` for tooling).
2. Run `uv lock` (if required) and `uv sync` to refresh `.venv`.
3. Document noteworthy changes here and in `docs/architecture.md`.

#### Adding New UI Controls

- Define the control in `Bin2CellExplorer.parameters`.
- Handle it in `inputTrigger`.
- Add the required state fields or backend parameters.
- Update documentation (`guide.md` or `architecture.md`) to keep feature parity.

#### Extending Backend Endpoints

- Follow the existing pattern: sanitize inputs, wrap logic in `try/except`, return JSON responses via `_json_response`.
- Prefer raising Python exceptions (the `Plugin` wrapper converts them into structured error payloads).
- For heavy operations, cache results via `TileCache`.

### Release Checklist

1. Verify `Bin2CellExplorer.yml` version.
2. Run through the usage steps end-to-end with sample data.
3. Ensure `Bin2CellExplorer.log` is free of unexpected errors.
4. Regenerate documentation if anything changed (`mkdocs build`).
5. Package plugin files with any sample data provided to users.

With these steps a new developer can reproduce the environment, understand the code layout, and extend the plugin safely.


