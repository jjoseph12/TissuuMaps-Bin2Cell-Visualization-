# Cell Explorer TissUUmaps Plugins


**[View Live Documentation](https://jjoseph12.github.io/Cell-Explorer-Documenation/)** | **[Interactive Ipynb Demo](plugins/Bin2Cell_demo.ipynb)**


Cell Explorer provides interactive gene and observation overlays for TissUUmaps. You can run it locally or on a cluster/server.

Cell Explorer (backend/UI for binned gene overlays) can be developed with the uv workflow or a traditional pip / virtualenv stack. You can pick whichever toolchain fits best—both produce identical dependencies and runtime behaviour.

I recomend using uv since its faster!


## QuickStart — Local UV workflow

1) Install uv
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2) Make uv available in your shell
```bash
export PATH="$HOME/.local/bin:$PATH"
```

3) Install the locked environment
```bash
cd ~/.tissuumaps
uv sync
```

4) Install files
```bash
git clone https://github.com/jjoseph12/Cell-Explorer-Visualization.git
```

5) Copy the plugin files into TissUUmaps
```bash
mkdir -p ~/.tissuumaps/plugins
cp /path/to/CellExplorer.js  ~/.tissuumaps/plugins/CellExplorer.js
cp /path/to/CellExplorer.py  ~/.tissuumaps/plugins/CellExplorer.py
cp /path/to/CellExplorer.yml ~/.tissuumaps/plugins/CellExplorer.yml

```
Keep the destination path exactly as above so TissUUmaps can find the plugin.

6) Launch TissUUmaps
```bash
uv run tissuumaps (or you can just say tissuumaps)

# Optional: silence Skia fallback
QTWEBENGINE_CHROMIUM_FLAGS="--disable-features=SkiaGraphite" uv run tissuumaps
```

## Server/Cluster — UV workflow

1) Install uv on the server and make it available
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"
```

2) Install the locked environment on the server
```bash
cd ~/.tissuumaps
uv sync
```
3) Install files
```bash
git clone https://github.com/jjoseph12/Cell-Explorer-Visualization.git
```

4) Install files and Copy the plugin files into the server’s TissUUmaps path
```bash
https://github.com/jjoseph12/Cell-Explorer-Visualization.git
mkdir -p ~/.tissuumaps/plugins
cp /path/to/CellExplorer.js  ~/.tissuumaps/plugins/CellExplorer.js
cp /path/to/CellExplorer.py  ~/.tissuumaps/plugins/CellExplorer.py
cp /path/to/CellExplorer.yml ~/.tissuumaps/plugins/CellExplorer.yml
```
Keep the destination as `~/.tissuumaps/plugins/CellExplorer*` on the server as well.

5) On the NYGC cluster, you can use the launch_tissuumaps.sh to launch and get started easier! It is a headless server (three modes). Path to slides means the path you want the file to open up on the UI directories otherwise it opens eveyrthing which is slow and might crash, for example: ~/.tissuumaps/launch_tissuumaps.sh --quiet /gpfs/commons/groups/innovation/jjoseph/

(you can always type in paths in the UI later)

Before make sure you give the script permission via: chmod +x launch_tissuumaps.sh

```bash
# Quiet (modules ON)
~/.tissuumaps/launch_tissuumaps.sh --quiet /path/to/slides
# Verbose (modules ON)
~/.tissuumaps/launch_tissuumaps.sh /path/to/slides
# No modules (optionally quiet) (I recomended using the modules otherwise, feel free to find the modules that work for your purpose)
~/.tissuumaps/launch_tissuumaps.sh --no-modules /path/to/slides
```
On your laptop, tunnel the port, (follow instructions in the launch_tissuumaps.sh):
```bash
ssh -L 5678:localhost:5678 <cluster-host>
# then open http://localhost:5678/
```

6) Easier alternatives (optional helpers)
- `install_bin2cell_plugin.sh` — clones the GitHub repo and installs the plugin into `~/.tissuumaps/plugins/`.
- `setup_uv_bin2cell.sh` — makes sure uv is installed, runs `uv sync`, copies the plugin, and can launch.

```bash
~/.tissuumaps/install_bin2cell_plugin.sh --persist-path --force
~/.tissuumaps/setup_uv_bin2cell.sh /path/to/CellExplorer --no-launch
```


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
- Slide assets (`.tif/.tiff`, `.npz`, `.h5ad`) can live anywhere. Use the plugin’s Browse… buttons or paste absolute paths.

## Notes & tips
- PySide6 is pinned to the 6.7 LTS line in the lock.
- Logs, presets, and cached state live under `~/.tissuumaps/plugins/CellExplorer/`.

## Expansion controls
- Max bin distance (µm): pixel distance uses `ceil(max_bin_distance * (bin_um / mpp))`
- Microns per pixel (mpp) and Bin size (µm): upstream acquisition/binning scalars
- Volume ratio (when `Expand mode = volume_ratio`): inflates area toward `ratio × original`, clipped to sane bounds
