# Design Rationale

This document records *why* the project exists, the goals it serves, and the technical trade-offs made along the way. It is meant for future maintainers and stakeholders to understand the decisions that shaped Bin2Cell Explorer.

## Why the project exists

- **Replace notebook-driven workflows**: The original `Bin2Cell_demo.ipynb` required researchers to edit code, rerun cells, and manually view Matplotlib plots. We wanted a point-and-click experience in TissUUmaps so lab members could explore data without Python expertise.
- **Bring Bin2Cell insights into TissUUmaps**: Users already rely on TissUUmaps for large pathology slides. Adding Bin2Cell overlays inside the viewer consolidates analysis and visualization in one tool.
- **Promote reproducibility**: By turning the notebook logic into a plugin with presets and exports, the workflows are shareable and less fragile than custom notebooks on individual machines.

## Key requirements

1. **Load existing outputs** (H&E TIFF, label NPZ, annotated `.h5ad`).
2. **Render per-tile overlays** for genes or observation columns with styling controls.
3. **Cache data aggressively** to keep the UI responsive.
4. **Provide intuitive navigation** so users can find tiles without guessing IDs.
5. **Export results** (GeoJSON) and let users save parameter presets.
6. **Handle large datasets** (~20k × 23k image and 4 GB `.h5ad`) on a typical lab workstation.

## Major technical decisions

### Dataset caching and rehydration
- **Decision**: Load the full dataset once, cache it in module-level state, and persist config to `dataset_state.json`.
- **Why**: Reloading a 4 GB `.h5ad` per request is too slow. Caching keeps overlay updates <1 s and allows TissUUmaps to restart without retyping paths.
- **Trade-off**: Higher memory footprint (~6 GB). Acceptable because most lab machines meet this requirement, and the alternative (per-request loading) is unusable.

### File selection workflow
- **Decision**: Keep explicit path fields but wire up dedicated **Browse…** buttons that trigger a Qt file dialog from the backend.
- **Why**: Users wanted to store datasets wherever they prefer yet avoid typing long absolute paths. Reusing Qt’s dialog preserves native UX and sidesteps browser sandbox limits.
- **Trade-off**: Requires a round-trip to the backend and depends on the Qt event loop; however, it keeps the UI simple and avoids maintaining a custom file explorer in JavaScript.

### Tile-by-tile rendering
- **Decision**: Compute overlays per tile rather than full-slide rasters or vector tiles.
- **Why**: Keeps computation and data transfer manageable, works within TissUUmaps’ existing marker/region framework, and mirrors the notebook’s workflow.
- **Trade-off**: Users must select tiles manually. Mitigated by adding the tile overview minimap and optional grid overlay.

### Separate gene vs. observation overlays
- **Decision**: Maintain two distinct code paths for gene expression and categorical observation overlays.
- **Why**: Presentation logic differs (color gradients vs. discrete palettes), and the backend can optimize each case separately.
- **Trade-off**: Slightly more code to maintain, but clearer separation improves readability and future extension (e.g., additional overlay types).

### Tile overview minimap
- **Problem**: Users had difficulty knowing which tile ID corresponded to which part of the slide.
- **Solution**: Render a lightweight scaled canvas showing every tile. Clicking selects and centers the tile.
- **Why not full grid overlay?** The on-slide grid overwhelmed the view when zoomed out. The minimap gives context without obscuring imagery and is cheap to render.

### Structured error responses
- **Decision**: Wrap every endpoint in `try/except` and convert exceptions into JSON.
- **Why**: TissUUmaps surfaces raw HTML when Flask raises unhandled exceptions, confusing users. The structured response keeps the UI informed and logs debug info.

### Logging
- **Decision**: Initialize a plugin-specific log file (`Bin2CellExplorer.log`) under the plugin directory.
- **Why**: Centralizes debugging info outside the TissUUmaps console, especially valuable when reproducing issues on user machines.

## Deferred/alternative ideas

- **Full-slide overlays**: generating a multi-resolution raster or vector tile set for the whole slide. Deferred due to scope and performance; per-tile rendering meets current needs.
- **Automatic tile suggestions**: scoring tiles by cell density or cluster diversity. Considered but left for future iteration.
- **Streaming overlays on viewport change**: would require deeper integration with OpenSeadragon events and asynchronous backend rendering. Possible enhancement if users need real-time updates while panning.

## Lessons learned

- Large datasets demand caching and careful memory management; even a small UX misstep (clicking load multiple times) can tank responsiveness.
- Navigation matters as much as rendering—users need intuitive ways to locate tiles before the overlay logic matters.
- Error handling should assume the UI may spam requests (user clicks before load finishes) and respond gracefully.

With these decisions documented, future contributors can understand the project’s intent, see which compromises were deliberate, and plan enhancements without repeating past experiments.
