# Usage Guide

This guide walks an end user through loading data, navigating the slide, and rendering overlays inside TissUUmaps.

## 1. Launch TissUUmaps with uv

```bash
cd ~/.tissuumaps
uv sync                              # only after dependency changes
QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite" uv run tissuumaps
```

The optional `QTWEBENGINE_CHROMIUM_FLAGS` silences the Skia Graphite → Ganesh fallback log. Omit it if you want the default logging.

## 2. Load the dataset once

Open **Plugins → Bin2Cell Explorer** and either paste the absolute paths into the three text fields or click the matching **Browse…** buttons to select each file with the system dialog (filters default to `.tif/.tiff`, `.npz`, and `.h5ad`).

Click **Load dataset**. The status line and log show `Dataset loaded successfully: … tiles`. The plugin caches the dataset so later overlay updates are instant and the selected paths persist between sessions.

> **Tip:** The initial load reads an 830 MB TIFF and a 4 GB `.h5ad`, so it can take ~20–30 seconds. Avoid clicking “Load dataset” repeatedly; the plugin will reuse the cached dataset automatically.

## 3. Orient yourself with the tile overview

Scroll a little in the plugin panel to the **Tile overview** minimap:

- Each tile is a rectangle labelled with its numeric ID.
- Clicking a rectangle selects that tile, updates the dropdown, highlights it in amber, and pans/zooms the main viewer directly to that region.

You can also toggle **Show tile grid** if you prefer to see tile outlines directly on the slide; the grid can be turned off at any time.

## 4. Choose your overlay

Select a tile via the minimap or dropdown, then set overlay parameters:

- **Gene overlay:** choose `Overlay type = gene`, enter a gene (e.g., `COL1A1`), pick color mode (`gradient`, `binary`, or `solid`), and click **Update overlay**. Use the colour picker below the dropdown: it adjusts the gradient ramp when `gradient` is selected, or sets the fill/stroke colour in `solid` mode.
- **Observation overlay:** choose `Overlay type = observation`, pick an observation column (e.g., `predicted_ClusterFull`), optionally filter one category, and click **Update overlay**.

Legend checkboxes let you toggle categories or gene overlays without re-rendering.

## 5. Optional layers and exports

- Toggle **Show expanded outlines**, **Show nuclei outlines**, or **Show centroids** for extra context.
- Click **Export overlay to GeoJSON** to capture the current overlay for downstream analysis.
- Save reusable settings via the **Presets** section.

## 6. Resetting or switching datasets

If you need to work with a different slide:

1. Update the file paths.
2. Click **Load dataset** once.
3. The plugin updates the cache (`dataset_state.json`) and refreshes the tile overview automatically.

That’s it—once the dataset is loaded you can hop between tiles, switch genes or observation columns, and export overlays without repeating the heavy load step.
