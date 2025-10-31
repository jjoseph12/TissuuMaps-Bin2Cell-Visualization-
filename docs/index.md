# Bin2Cell Explorer

Bin2Cell Explorer is a TissUUmaps plugin that brings the interactive features from the original `Bin2Cell_demo.ipynb` notebook into the viewer. It lets scientists explore Bin2Cell-expanded labels on top of H&E imagery, switch between gene-expression and observation overlays, and export the rendered regions.

## Core Features

- **Dataset loader** for H&E TIFF, Bin2Cell label NPZ, and annotated AnnData (`.h5ad`).
- **Gene overlays**: color expanded-cell regions by expression (gradient or binary).
- **Observation overlays**: majority-vote categorical labels per expanded cell.
- **Tile navigation**: dropdown, a tile overview minimap, optional on-slide grid.
- **Per-tile rendering**: users choose any tile and re-render quickly.
- **Caching and persistence**: dataset only loads once; configuration is stored for re-hydration across sessions.
- **Export tools**: write the current overlay to GeoJSON; save presets for frequently used parameter sets.

This documentation explains how the plugin is structured, the environment it expects, and the engineering decisions behind each part so future developers can maintain and extend it confidently.
