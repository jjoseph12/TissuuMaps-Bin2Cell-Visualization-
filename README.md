# TissUUmaps Plugins

Bin2Cell Explorer (backend/UI for binned gene overlays) can be developed with either the new uv workflow or a traditional pip / virtualenv stack. You can pick whichever toolchain fits they both produce identical dependencies and runtime behaviour.

---

## Option A — UV-managed environment
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

## Option B — pip / virtualenv
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
   `/Applications/TissUUmaps.app/Contents/MacOS/TissUUmaps`

When you are finished, run `deactivate` to exit the virtual environment.

---

## Selecting dataset files
- Slide assets (`.tif/.tiff`, `.npz`, `.h5ad`) can live anywhere. Use the **Browse…** buttons in the plugin panel or paste absolute paths manually.
- The browse dialog filters the relevant extensions and rejects folders to prevent loading mistakes.

## Notes & tips
- PySide6 is pinned to the 6.7 line to avoid noisy Skia “Graphite → Ganesh” fallback messages if you are on MacOS and not on the most recent Chromium browswer; drop the `<6.8` pin if you prefer the latest Qt (expect the log noise to return).
- Gene overlays support custom colours: choose **Solid** for a single hex colour or stay on **Gradient** and tweak the gradient colour picker to bias the ramp.
- `uv.lock` captures exact versions for the uv workflow; `.venv/` is ignored by git in both setups.
- Logs, presets, and cached state live under `~/.tissuumaps/plugins/Bin2CellExplorer/`.
