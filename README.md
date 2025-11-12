# Cell Explorer TissUUmaps Plugins


**[View Documentation](https://jjoseph12.github.io/Cell-Explorer-Documenation/)** | 🔬 **[Interactive Demo](plugins/Bin2Cell_demo.ipynb)**


Cell Explorer provides interactive gene and observation overlays for TissUUmaps. You can run it locally or on a cluster/server.

Cell Explorer (backend/UI for binned gene overlays) can be developed with the uv workflow or a traditional pip / virtualenv stack. You can pick whichever toolchain fits best—both produce identical dependencies and runtime behaviour.

I recomend using uv since its faster!


## Quick start (UV workflow)

1) Install uv (once)

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2) Make the shim available in your shell

```bash
export PATH="$HOME/.local/bin:$PATH"
```

3) Install the locked environment from the project root

```bash
cd ~/.tissuumaps
uv sync
```

4) Launch TissUUmaps

```bash
uv run tissuumaps

# Optional: silence Skia fallback
QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite"

uv run tissuumaps
```

5) Interactive shell (optional)

```bash
source .venv/bin/activate
# ... run commands ...
deactivate
```

## Server/cluster setup (recommended scripts)

All scripts live in `~/.tissuumaps/`.

- `install_bin2cell_plugin.sh` — install the plugin from GitHub into `~/.tissuumaps/plugins/`
  - Accepts either `CellExplorer.*` or legacy `Bin2CellExplorer.*` from the repo
  - Always installs as `CellExplorer.js`, `CellExplorer.py`, `CellExplorer.yml`

```bash
chmod +x ~/.tissuumaps/install_bin2cell_plugin.sh
~/.tissuumaps/install_bin2cell_plugin.sh --persist-path --force
# Options:
#   --repo URL        # alternative repo
#   --branch BRANCH   # branch to clone
#   --target DIR      # destination (default: ~/.tissuumaps/plugins)
#   --force           # overwrite existing files
#   --persist-path    # append PATH export to ~/.bashrc
```

- `setup_uv_bin2cell.sh` — convenience setup that ensures uv is present, runs `uv sync`, copies plugin files, and optionally launches

```bash
chmod +x ~/.tissuumaps/setup_uv_bin2cell.sh
~/.tissuumaps/setup_uv_bin2cell.sh /path/to/CellExplorer --no-launch
# Add --quiet-skia to suppress Skia fallback logs on launch
```

- `launch_tissuumaps.sh` — headless launcher with three modes
  - Default (verbose, modules ON)
  - `--quiet` (silences logs; prints only SSH tunnel tip)
  - `--no-modules` (skip environment module loads; can combine with `--quiet`)

```bash
# Quiet (modules ON)
~/.tissuumaps/launch_tissuumaps.sh --quiet /path/to/slides

# Verbose (modules ON)
~/.tissuumaps/launch_tissuumaps.sh /path/to/slides

# No modules (optionally quiet)
~/.tissuumaps/launch_tissuumaps.sh --no-modules /path/to/slides
~/.tissuumaps/launch_tissuumaps.sh --no-modules --quiet /path/to/slides
```

Tip for remote servers: from your laptop

```bash
ssh -L 5678:localhost:5678 <cluster-host>
# then open http://localhost:5678/
```

## Selecting dataset files

- Slide assets (`.tif/.tiff`, `.npz`, `.h5ad`) can live anywhere. Use the plugin’s Browse… buttons or paste absolute paths.
- The browse dialog filters relevant extensions and rejects folders.

## Notes & tips

- PySide6 is pinned to the 6.7 LTS line in the lock; if you bump Qt, expect more Skia fallback log chatter.
- Gene overlays support custom colours:
  - Use “Solid” for a single hex colour
  - Or “Gradient” and adjust the gradient colour pickers
- Observation overlays auto-populate category dropdowns from AnnData (respecting `<column>_colors`).
- Exact versions are captured in `uv.lock`; `.venv/` is ignored in both setups.
- Logs, presets, and cached state live under `~/.tissuumaps/plugins/CellExplorer/`.

## Expansion controls (ex-Bin2Cell knobs)

- Max bin distance (µm): upper bound for nucleus growth; pixel distance uses `ceil(max_bin_distance * (bin_um / mpp))`
- Microns per pixel (mpp): acquisition resolution
- Bin size (µm): upstream Bin2Cell bin width used for scaling
- Volume ratio: when `Expand mode = volume_ratio`, inflates each nucleus toward `volume_ratio × (original area)` with sane clipping

These four controls are linked: start with the acquisition `mpp` and upstream `bin_um`, then adjust `max_bin_distance` (or `volume_ratio`) to tune halo overlap.




---

## Recommended UV-managed environment
1. Install uv (once):  
   `curl -LsSf https://astral.sh/uv/install.sh | sh`
2. Make the shim available (add to `~/.zshrc` or run per shell):  
   `export PATH="$HOME/.local/bin:$PATH"`
3. From the project root (`cd ~/.tissuumaps`) install the locked environment:  
   `uv sync`
4. Launch TissUUmaps with the managed environment:  
   `uv run tissuumaps`  
   (Optional: silence the Skia fallback log with `QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite"`.)
5. For an interactive shell  
   `source .venv/bin/activate` → run commands → `deactivate`

Common uv commands:
- `uv add <package>` / `uv add --group dev <package>` — add runtime/dev deps
- `uv tree` — inspect the locked dependency graph

The environment installs runtime libraries (`anndata`, `numpy`, `scipy`, `scikit-image`, `imagecodecs`, `flask`, `matplotlib`, `PySide6` 6.7 LTS, etc.) plus optional tooling (`pytest`, `mypy`, `mkdocs`).

---



## Option B pip / virtualenv
1. Create and activate a virtual environment:
   ```bash
   cd ~/.tissuumaps
   python3 -m venv .venv
   source .venv/bin/activate
   ```
2. Upgrade packaging tools (recommended):
   ```bash
   python -m pip install --upgrade pip setuptools wheel
   ```
3. Install the dependencies listed in `pyproject.toml`:
   ```bash
   python -m pip install \
     anndata flask matplotlib numpy scikit-image scipy imagecodecs \
     tissuumaps PySide6>=6.7,<6.8 cmake~=3.29 \
     pytest mypy mkdocs mkdocs-material
   ```
   (PySide6 6.7 LTS mirrors the uv lock; adjust if you need a newer Qt.)
4. Launch the packaged app from the activated shell:  
   `tissuumaps`

When you are finished, run `deactivate` to exit the virtual environment.

---

## Selecting dataset files
- Slide assets (`.tif/.tiff`, `.npz`, `.h5ad`) can live anywhere. Use the **Browse…** buttons in the plugin panel or paste absolute paths manually.
- The browse dialog filters the relevant extensions and rejects folders to prevent loading mistakes.

## Notes & tips
- PySide6 is pinned to the 6.7 line to avoid noisy Skia “Graphite → Ganesh” fallback messages if you are on MacOS and not on the most recent Chromium browswer; drop the `<6.8` pin if you prefer the latest Qt (expect the log noise to return).
- Gene overlays support custom colours: choose **Solid** for a single hex colour or stay on **Gradient** and tweak the gradient colour picker to bias the ramp.
- Observation overlays auto-populate their category dropdowns from AnnData (honouring any `<column>_colors` entries), so hierarchical cluster levels such as `predicted_ClusterFull/Midway/Top` inherit their curated palettes without retyping category names.
- `uv.lock` captures exact versions for the uv workflow; `.venv/` is ignored by git in both setups.
- Logs, presets, and cached state live under `~/.tissuumaps/plugins/Bin2CellExplorer/`.

## Bin2Cell expansion knobs
- **Max bin distance (µm):** Sets the upper bound for how far each nucleus can grow. The plugin converts it to pixels via `distance_px = ceil(max_bin_distance * (bin_um / mpp))`, so doubling the distance or the bin size doubles the halo in pixel space.
- **Microns per pixel (mpp):** Slide resolution. Lower values mean more pixels per micron, which directly affects the conversion above.
- **Bin size (µm):** Physical width of the Bin2Cell bins that were used upstream; keeping this value in sync with the preprocessing notebook ensures spatial scaling stays accurate.
- **Volume ratio:** Only used when `Expand mode = volume_ratio`. Instead of a fixed cap, it inflates each nucleus so the expanded area approximates `volume_ratio × (original area)`; the plugin converts that area delta to a median radial offset and clips it to sane bounds.

All four knobs are linked: the first three determine the pixel distance for `fixed` mode, and `volume_ratio` reuses the same unit conversions when estimating its per-label offsets. Start with the acquisition `mpp`, the Bin2Cell `bin_um`, then adjust `max_bin_distance` (or `volume_ratio`) to control how much the halos overlap.
