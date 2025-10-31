import json
import logging
import math
import os
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from flask import make_response

try:
    import anndata as ad
except ImportError:  # pragma: no cover - runtime guard
    ad = None

from scipy.sparse import load_npz, issparse
from skimage import io, measure
from skimage.segmentation import expand_labels, find_boundaries

from matplotlib import cm, colors
from PySide6.QtCore import QMetaObject, QObject, Qt, Q_ARG, Slot
from PySide6.QtWidgets import QApplication, QFileDialog

LOGGER = logging.getLogger(__name__)

_GLOBAL_STATE: Dict[str, object] = {}
_FILE_DIALOG_HELPER: Optional["FileDialogHelper"] = None


def _rgba_to_css(rgba: Tuple[float, float, float, float], *, alpha_override: Optional[float] = None) -> str:
    r, g, b, a = rgba
    if alpha_override is not None:
        a = alpha_override
    return f"rgba({int(round(r * 255))},{int(round(g * 255))},{int(round(b * 255))},{round(a, 4)})"


def _get_cmap(cmap: object) -> "colors.Colormap":
    if isinstance(cmap, str):
        return cm.get_cmap(cmap)
    return cmap  # type: ignore[return-value]


def _colormap_sample(cmap: object, value: float, alpha: Optional[float] = None) -> str:
    cmap_obj = _get_cmap(cmap)
    rgba = cmap(np.clip(value, 0.0, 1.0))
    return _rgba_to_css(rgba, alpha_override=alpha)


def _sample_gradient(cmap: object, steps: int = 8) -> List[Tuple[float, str]]:
    cmap_obj = _get_cmap(cmap)
    gradient = []
    for i in range(steps):
        pos = i / max(steps - 1, 1)
        gradient.append((pos, _rgba_to_css(cmap_obj(pos))))
    return gradient


def _unique_sorted(iterable: Iterable) -> List:
    return sorted(set(iterable))


@dataclass(frozen=True)
class Tile:
    id: int
    r0: int
    r1: int
    c0: int
    c1: int

    def to_dict(self) -> Dict[str, int]:
        return {"id": self.id, "r0": self.r0, "r1": self.r1, "c0": self.c0, "c1": self.c1}


class B2CContext:
    """
    Lightweight loader for H&E image, sparse nuclei labels, and AnnData container.
    Provides centroid pixel coordinates and convenience accessors.
    """

    def __init__(
        self,
        he_image_path: str,
        labels_npz_path: str,
        adata_path: str,
        *,
        obsm_key: Optional[str] = None,
    ):
        if ad is None:
            raise ImportError("anndata is required (pip install anndata scanpy).")

        he_image_path = os.path.abspath(os.path.expanduser(os.path.expandvars(he_image_path)))
        labels_npz_path = os.path.abspath(os.path.expanduser(os.path.expandvars(labels_npz_path)))
        adata_path = os.path.abspath(os.path.expanduser(os.path.expandvars(adata_path)))

        if not os.path.exists(he_image_path):
            raise FileNotFoundError(f"H&E image not found: {he_image_path}")
        if not os.path.exists(labels_npz_path):
            raise FileNotFoundError(f"Label NPZ not found: {labels_npz_path}")
        if not os.path.exists(adata_path):
            raise FileNotFoundError(f"AnnData (.h5ad) not found: {adata_path}")

        LOGGER.info("Loading H&E image: %s", he_image_path)
        self.he = io.imread(he_image_path)

        LOGGER.info("Loading sparse labels: %s", labels_npz_path)
        self.lab_sp = load_npz(labels_npz_path)

        if self.lab_sp.shape != self.he.shape[:2]:
            raise ValueError(
                f"Image/label size mismatch: H&E {self.he.shape[:2]} vs labels {self.lab_sp.shape}"
            )

        LOGGER.info("Loading AnnData: %s", adata_path)
        self.adata = ad.read_h5ad(adata_path)

        self._gene_cache: Dict[str, np.ndarray] = {}
        self._obs_cache: Dict[str, np.ndarray] = {}

        self._available_obsm = [
            k
            for k, arr in self.adata.obsm.items()
            if hasattr(arr, "shape") and arr.shape[1] >= 2
        ]
        if not self._available_obsm:
            raise ValueError("No obsm entries with >=2 columns were found in AnnData.")

        self.obsm_key = self._resolve_obsm_key(obsm_key)
        xy = np.asarray(self.adata.obsm[self.obsm_key])
        self.cols = np.rint(xy[:, 0]).astype(int)
        self.rows = np.rint(xy[:, 1]).astype(int)

        H, W = self.shape
        np.clip(self.cols, 0, W - 1, out=self.cols)
        np.clip(self.rows, 0, H - 1, out=self.rows)

    def _resolve_obsm_key(self, preferred: Optional[str]) -> str:
        if preferred and preferred in self._available_obsm:
            return preferred
        for candidate in ("spatial_cropped_150_buffer", "spatial"):
            if candidate in self._available_obsm:
                return candidate
        return self._available_obsm[0]

    @property
    def shape(self) -> Tuple[int, int]:
        return self.lab_sp.shape

    @property
    def available_obsm(self) -> List[str]:
        return list(self._available_obsm)

    @property
    def gene_names(self) -> List[str]:
        return list(map(str, self.adata.var_names))

    @property
    def obs_columns(self) -> List[str]:
        return list(map(str, self.adata.obs.columns))

    def crop_dense(self, r0: int, r1: int, c0: int, c1: int) -> Tuple[np.ndarray, np.ndarray]:
        he = self.he[r0:r1, c0:c1]
        lab = self.lab_sp[r0:r1, c0:c1].toarray().astype(np.int32, copy=False)
        return he, lab

    def gene_vector(self, gene: str) -> np.ndarray:
        if gene not in self._gene_cache:
            if gene not in self.adata.var_names:
                raise KeyError(f"Gene '{gene}' not found in AnnData.")
            values = self.adata[:, gene].X
            if issparse(values):
                values = values.toarray()
            self._gene_cache[gene] = np.asarray(values).ravel()
        return self._gene_cache[gene]

    def obs_vector(self, column: str) -> np.ndarray:
        if column not in self._obs_cache:
            if column not in self.adata.obs.columns:
                raise KeyError(f"Column '{column}' not found in AnnData.obs.")
            self._obs_cache[column] = self.adata.obs[column].to_numpy()
        return self._obs_cache[column]


def make_tiles(
    ctx: B2CContext, tile_h: int = 1500, tile_w: int = 1500, stride_h: Optional[int] = None, stride_w: Optional[int] = None
) -> List[Tile]:
    H, W = ctx.shape
    stride_h = tile_h if stride_h is None else stride_h
    stride_w = tile_w if stride_w is None else stride_w

    tiles: List[Tile] = []
    tid = 0
    r = 0
    while r < H:
        r1 = min(r + tile_h, H)
        c = 0
        while c < W:
            c1 = min(c + tile_w, W)
            tiles.append(Tile(tid, r, r1, c, c1))
            tid += 1
            if c1 == W:
                break
            c += stride_w
        if r1 == H:
            break
        r += stride_h
    return tiles


def _compute_expansion_distance(
    lab_raw: np.ndarray,
    *,
    mode: str,
    max_bin_distance: float,
    mpp: float,
    bin_um: float,
    volume_ratio: float,
) -> int:
    dist_px_fixed = int(math.ceil(max(1.0, max_bin_distance * (bin_um / mpp))))
    if mode != "volume_ratio":
        return dist_px_fixed

    lab_ids, counts = np.unique(lab_raw[lab_raw > 0], return_counts=True)
    if lab_ids.size == 0:
        return dist_px_fixed

    r_eff = np.sqrt(counts / np.pi)
    delta_r = (np.sqrt(volume_ratio) - 1.0) * r_eff
    dist_px = int(np.clip(np.median(delta_r), 1, 5 * dist_px_fixed))
    return max(1, dist_px)


def _expand_labels_tile(
    ctx: B2CContext,
    tile: Tile,
    *,
    mode: str,
    max_bin_distance: float,
    mpp: float,
    bin_um: float,
    volume_ratio: float,
    pad_factor: int = 2,
) -> Tuple[np.ndarray, np.ndarray, int]:
    r0, r1, c0, c1 = tile.r0, tile.r1, tile.c0, tile.c1
    he_crop, lab_raw = ctx.crop_dense(r0, r1, c0, c1)

    dist_px = _compute_expansion_distance(
        lab_raw,
        mode=mode,
        max_bin_distance=max_bin_distance,
        mpp=mpp,
        bin_um=bin_um,
        volume_ratio=volume_ratio,
    )

    H, W = ctx.shape
    pad = pad_factor * dist_px

    rp0 = max(r0 - pad, 0)
    rp1 = min(r1 + pad, H)
    cp0 = max(c0 - pad, 0)
    cp1 = min(c1 + pad, W)

    lab_pad = ctx.lab_sp[rp0:rp1, cp0:cp1].toarray().astype(np.int32, copy=False)
    lab_exp_pad = expand_labels(lab_pad, distance=dist_px)

    r0_rel, r1_rel = r0 - rp0, r1 - rp0
    c0_rel, c1_rel = c0 - cp0, c1 - cp0

    lab_exp = lab_exp_pad[r0_rel:r1_rel, c0_rel:c1_rel]
    return he_crop, lab_exp.astype(np.int32, copy=False), int(dist_px)


def _polygons_from_labels(lab: np.ndarray, tile: Tile) -> Dict[int, List[List[List[float]]]]:
    r0, c0 = tile.r0, tile.c0
    output: Dict[int, List[List[List[float]]]] = {}
    labels = np.unique(lab)
    for lbl in labels:
        if lbl <= 0:
            continue
        mask = lab == lbl
        if not mask.any():
            continue
        contours = measure.find_contours(mask.astype(np.uint8), 0.5)
        poly_list: List[List[List[float]]] = []
        for contour in contours:
            if contour.shape[0] < 3:
                continue
            polygon: List[List[float]] = []
            for y, x in contour:
                polygon.append([float(x + c0), float(y + r0)])
            if polygon and polygon[0] != polygon[-1]:
                polygon.append(polygon[0])
            if polygon:
                poly_list.append(polygon)
        if poly_list:
            output[int(lbl)] = poly_list
    return output


def _outline_paths(lab: np.ndarray, tile: Tile) -> List[List[List[float]]]:
    r0, c0 = tile.r0, tile.c0
    edges = find_boundaries(lab, mode="inner")
    contours = measure.find_contours(edges.astype(np.uint8), 0.5)
    paths: List[List[List[float]]] = []
    for contour in contours:
        if contour.shape[0] < 2:
            continue
        path: List[List[float]] = []
        for y, x in contour:
            path.append([float(x + c0), float(y + r0)])
        paths.append(path)
    return paths


def _centroids_for_tile(ctx: B2CContext, tile: Tile) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r0, r1, c0, c1 = tile.r0, tile.r1, tile.c0, tile.c1
    in_tile = (ctx.rows >= r0) & (ctx.rows < r1) & (ctx.cols >= c0) & (ctx.cols < c1)
    indices = np.where(in_tile)[0]
    if indices.size == 0:
        empty_int = np.empty((0,), dtype=np.int32)
        empty_local = np.empty((0, 2), dtype=np.int32)
        return empty_int, empty_int, empty_int, empty_local

    rr_local = (ctx.rows[indices] - r0).astype(np.int32)
    cc_local = (ctx.cols[indices] - c0).astype(np.int32)
    local = np.stack([rr_local, cc_local], axis=1)
    return indices, ctx.rows[indices], ctx.cols[indices], local


class TileCache:
    def __init__(self, max_items: int = 6):
        self._max = max_items
        self._data: "OrderedDict[Tuple, Dict]" = OrderedDict()

    def clear(self) -> None:
        self._data.clear()

    def get(self, key: Tuple) -> Optional[Dict]:
        entry = self._data.get(key)
        if entry is not None:
            self._data.move_to_end(key)
        return entry

    def set(self, key: Tuple, value: Dict) -> None:
        self._data[key] = value
        self._data.move_to_end(key)
        while len(self._data) > self._max:
            self._data.popitem(last=False)


class Plugin:
    def __init__(self, app):
        self.app = app
        self.out_dir = os.path.expanduser("~/.tissuumaps/plugins/Bin2CellExplorer")
        os.makedirs(self.out_dir, exist_ok=True)
        self.log_path = os.path.join(self.out_dir, "Bin2CellExplorer.log")
        if not LOGGER.handlers:
            handler = logging.FileHandler(self.log_path)
            handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
            LOGGER.addHandler(handler)
        LOGGER.setLevel(logging.INFO)

        global _GLOBAL_STATE
        if not _GLOBAL_STATE:
            _GLOBAL_STATE = {
                "context": None,
                "tiles": [],
                "dataset_id": None,
                "tile_cache": TileCache(max_items=4),
                "dataset_config": None,
                "loading": False,
            }
        self._state = _GLOBAL_STATE
        self.current_params: Dict[str, str] = {}
        self.state_path = os.path.join(self.out_dir, "dataset_state.json")
        if self._state.get("dataset_config") is None:
            cached = self._load_cached_config()
            if cached:
                self._state["dataset_config"] = cached

        self.presets_path = os.path.join(self.out_dir, "presets.json")
        self.presets = self._load_presets()

    # ------------------------------------------------------------------ Helpers

    def _json_response(self, payload: Dict) -> "flask.Response":
        return make_response(json.dumps(payload, default=_json_default), 200, {"Content-Type": "application/json"})

    def _error_response(self, status: int, message: str, exc: Optional[Exception] = None):
        if exc is not None:
            LOGGER.exception(message)
        else:
            LOGGER.error(message)
        payload = {"status": "error", "message": message}
        if exc is not None:
            payload["detail"] = repr(exc)
        return make_response(json.dumps(payload), status, {"Content-Type": "application/json"})

    def _error_from_exception(self, exc: Exception, context: str):
        if isinstance(exc, FileNotFoundError):
            status = 404
            message = str(exc)
        elif isinstance(exc, (ValueError, KeyError, RuntimeError)):
            status = 400
            message = str(exc)
        elif isinstance(exc, ImportError):
            status = 500
            message = f"{context}: missing dependency ({exc})"
        else:
            status = 500
            message = f"{context}: {exc}"
        return self._error_response(status, message, exc)

    def _ensure_context(self) -> B2CContext:
        if self.context is None:
            config = self._state.get("dataset_config") or self._load_cached_config()
            if config:
                LOGGER.info("Rehydrating dataset from cached config.")
                ctx = B2CContext(
                    he_image_path=config["he_path"],
                    labels_npz_path=config["labels_path"],
                    adata_path=config["adata_path"],
                    obsm_key=config.get("obsm_key"),
                )
                tiles = make_tiles(
                    ctx,
                    tile_h=config.get("tile_h", 1500) or 1500,
                    tile_w=config.get("tile_w", 1500) or 1500,
                    stride_h=config.get("stride_h"),
                    stride_w=config.get("stride_w"),
                )
                self.context = ctx
                self.tiles = tiles
                self.dataset_id = os.path.basename(config["adata_path"])
                self.tile_cache.clear()
                self._state["dataset_config"] = config
            else:
                raise RuntimeError("Load a dataset first.")
        return self.context

    @property
    def context(self) -> Optional[B2CContext]:
        return self._state.get("context")  # type: ignore[return-value]

    @context.setter
    def context(self, value: Optional[B2CContext]) -> None:
        self._state["context"] = value

    @property
    def tiles(self) -> List[Tile]:
        return self._state.get("tiles", [])  # type: ignore[return-value]

    @tiles.setter
    def tiles(self, value: List[Tile]) -> None:
        self._state["tiles"] = value

    @property
    def dataset_id(self) -> Optional[str]:
        return self._state.get("dataset_id")  # type: ignore[return-value]

    @dataset_id.setter
    def dataset_id(self, value: Optional[str]) -> None:
        self._state["dataset_id"] = value

    @property
    def tile_cache(self) -> TileCache:
        cache = self._state.get("tile_cache")
        if cache is None:
            cache = TileCache(max_items=4)
            self._state["tile_cache"] = cache
        return cache  # type: ignore[return-value]

    def _load_presets(self) -> Dict[str, Dict]:
        if not os.path.exists(self.presets_path):
            return {}
        try:
            with open(self.presets_path, "r") as fh:
                presets = json.load(fh)
            return presets if isinstance(presets, dict) else {}
        except Exception as exc:  # pragma: no cover - runtime guard
            LOGGER.warning("Failed to read presets: %s", exc)
            return {}

    def _save_presets(self) -> None:
        with open(self.presets_path, "w") as fh:
            json.dump(self.presets, fh, indent=2)

    # ------------------------------------------------------------------ Endpoints

    def pick_file(self, params: Dict):
        try:
            field = str(params.get("field") or "")
            if field not in {"he_path", "labels_path", "adata_path"}:
                raise ValueError("Unsupported field for file picker.")

            current = params.get("current_path") or ""
            start_dir = params.get("start_dir") or ""
            if current:
                current = os.path.abspath(os.path.expanduser(os.path.expandvars(str(current))))
                if os.path.isdir(current):
                    start_dir = current
                else:
                    start_dir = os.path.dirname(current)
            if not start_dir:
                start_dir = os.path.expanduser("~")

            filters = {
                "he_path": "H&E images (*.tif *.tiff *.TIF *.TIFF);;All files (*)",
                "labels_path": "Label matrices (*.npz *.NPZ);;All files (*)",
                "adata_path": "AnnData (*.h5ad *.H5AD);;All files (*)",
            }
            captions = {
                "he_path": "Select H&E image",
                "labels_path": "Select label matrix (.npz)",
                "adata_path": "Select AnnData (.h5ad)",
            }

            app = QApplication.instance()
            if app is None:
                raise RuntimeError("QApplication instance not found; cannot open file dialog.")

            global _FILE_DIALOG_HELPER
            if _FILE_DIALOG_HELPER is None:
                helper = FileDialogHelper()
                helper.moveToThread(app.thread())
                _FILE_DIALOG_HELPER = helper

            filename, _ = _FILE_DIALOG_HELPER.get_open_file_name(
                captions[field],
                start_dir,
                filters[field],
            )

            if not filename:
                return self._json_response({"status": "cancelled"})

            path = os.path.abspath(filename)
            if not os.path.isfile(path):
                raise FileNotFoundError(f"Selected path is not a file: {path}")
            return self._json_response({"status": "ok", "path": path, "field": field})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to pick file")

    def _load_cached_config(self) -> Optional[Dict[str, object]]:
        if not hasattr(self, "state_path"):
            return None
        if not os.path.exists(self.state_path):
            return None
        try:
            with open(self.state_path, "r") as fh:
                data = json.load(fh)
            if isinstance(data, dict):
                return data
        except Exception as exc:
            LOGGER.warning("Failed to read cached dataset config: %s", exc)
        return None

    def load_dataset(self, params: Dict):
        try:
            self._state["loading"] = True
            he_path = params.get("he_path", "")
            labels_path = params.get("labels_path", "")
            adata_path = params.get("adata_path", "")
            obsm_key = params.get("obsm_key")
            tile_h = int(params.get("tile_h", 1500) or 1500)
            tile_w = int(params.get("tile_w", 1500) or 1500)
            stride_h_param = params.get("stride_h")
            stride_w_param = params.get("stride_w")
            stride_h = int(stride_h_param) if stride_h_param not in (None, "", "None") else None
            stride_w = int(stride_w_param) if stride_w_param not in (None, "", "None") else None

            requested_config = {
                "he_path": he_path,
                "labels_path": labels_path,
                "adata_path": adata_path,
                "obsm_key": obsm_key,
                "tile_h": tile_h,
                "tile_w": tile_w,
                "stride_h": stride_h,
                "stride_w": stride_w,
            }

            existing_config = self._state.get("dataset_config")
            if existing_config and all(existing_config.get(k) == requested_config.get(k) for k in requested_config):
                LOGGER.info("Dataset already loaded; reusing cached context.")
                ctx = self._ensure_context()
                tiles = self.tiles
                payload = {
                    "status": "ok",
                    "tile_count": len(tiles),
                    "tiles": [tile.to_dict() for tile in tiles],
                    "obsm_key": ctx.obsm_key,
                    "available_obsm": ctx.available_obsm,
                    "obs_columns": ctx.obs_columns,
                    "gene_count": len(ctx.gene_names),
                    "genes_preview": ctx.gene_names[:512],
                    "dataset_id": self.dataset_id,
                    "shape": ctx.shape,
                    **existing_config,
                }
                return self._json_response(payload)

            ctx = B2CContext(
                he_image_path=he_path,
                labels_npz_path=labels_path,
                adata_path=adata_path,
                obsm_key=obsm_key,
            )
            tiles = make_tiles(ctx, tile_h=tile_h, tile_w=tile_w, stride_h=stride_h, stride_w=stride_w)

            self.context = ctx
            self.tiles = tiles
            self.tile_cache.clear()
            self.dataset_id = os.path.basename(adata_path)

            dataset_config = {
                "he_path": he_path,
                "labels_path": labels_path,
                "adata_path": adata_path,
                "obsm_key": ctx.obsm_key,
                "tile_h": tile_h,
                "tile_w": tile_w,
                "stride_h": stride_h,
                "stride_w": stride_w,
            }

            payload = {
                "status": "ok",
                "tile_count": len(tiles),
                "tiles": [tile.to_dict() for tile in tiles],
                "obsm_key": ctx.obsm_key,
                "available_obsm": ctx.available_obsm,
                "obs_columns": ctx.obs_columns,
                "gene_count": len(ctx.gene_names),
                "genes_preview": ctx.gene_names[:512],
                "dataset_id": self.dataset_id,
                "shape": ctx.shape,
                **dataset_config,
            }
            self._state["dataset_config"] = dataset_config
            try:
                with open(self.state_path, "w") as fh:
                    json.dump(dataset_config, fh)
            except Exception as exc:
                LOGGER.warning("Failed to persist dataset config: %s", exc)
            LOGGER.info("Dataset loaded successfully: %d tiles", len(tiles))
            return self._json_response(payload)
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to load dataset")
        finally:
            self._state["loading"] = False

    def describe_dataset(self, params: Dict):
        try:
            ctx = self._ensure_context()
            payload = {
                "dataset_id": self.dataset_id,
                "shape": ctx.shape,
                "tile_count": len(self.tiles),
                "obsm_key": ctx.obsm_key,
                "available_obsm": ctx.available_obsm,
                "obs_columns": ctx.obs_columns,
                "gene_count": len(ctx.gene_names),
            }
            return self._json_response(payload)
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to describe dataset")

    def list_tiles(self, params: Dict):
        try:
            self._ensure_context()
            return self._json_response({"tiles": [tile.to_dict() for tile in self.tiles]})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to list tiles")

    def list_obs_columns(self, params: Dict):
        try:
            ctx = self._ensure_context()
            cols = ctx.obs_columns
            prefer = []
            fallback = []
            for col in cols:
                series = ctx.obs_vector(col)
                if series.dtype.kind in {"O", "U", "S"} or str(series.dtype).startswith("category"):
                    prefer.append(col)
                else:
                    fallback.append(col)
            ordered = prefer + [c for c in fallback if c not in prefer]
            return self._json_response({"columns": ordered})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to list observation columns")

    def list_genes(self, params: Dict):
        try:
            ctx = self._ensure_context()
            limit = int(params.get("limit", 0) or 0)
            genes = ctx.gene_names
            if limit > 0:
                genes = genes[:limit]
            return self._json_response({"genes": genes})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to list genes")

    # ------------------------------------------------------------------ Overlay core

    def _tile_entry(
        self,
        tile_id: int,
        *,
        b2c_mode: str,
        max_bin_distance: float,
        mpp: float,
        bin_um: float,
        volume_ratio: float,
        pad_factor: int = 2,
    ) -> Dict:
        ctx = self._ensure_context()
        if tile_id < 0 or tile_id >= len(self.tiles):
            raise ValueError(f"Tile id {tile_id} out of range.")
        tile = self.tiles[tile_id]

        cache_key = (
            self.dataset_id,
            tile_id,
            b2c_mode,
            round(max_bin_distance, 3),
            round(mpp, 4),
            round(bin_um, 4),
            round(volume_ratio, 4),
            pad_factor,
        )
        cached = self.tile_cache.get(cache_key)
        if cached is not None:
            return cached

        _, lab_raw = ctx.crop_dense(tile.r0, tile.r1, tile.c0, tile.c1)
        he_crop, lab_exp, dist_px = _expand_labels_tile(
            ctx,
            tile,
            mode=b2c_mode,
            max_bin_distance=max_bin_distance,
            mpp=mpp,
            bin_um=bin_um,
            volume_ratio=volume_ratio,
            pad_factor=pad_factor,
        )

        centroid_idx, rows_abs, cols_abs, local_rc = _centroids_for_tile(ctx, tile)
        labels_at_centroids = lab_exp[local_rc[:, 0], local_rc[:, 1]] if local_rc.size else np.empty((0,), dtype=np.int32)

        polygons_exp = _polygons_from_labels(lab_exp, tile)
        polygons_raw = _polygons_from_labels(lab_raw, tile)

        outline_exp = _outline_paths(lab_exp, tile)
        outline_raw = _outline_paths(lab_raw, tile)

        entry = {
            "tile": tile,
            "lab_exp": lab_exp,
            "lab_raw": lab_raw,
            "dist_px": dist_px,
            "polygons_exp": polygons_exp,
            "polygons_raw": polygons_raw,
            "outline_exp": outline_exp,
            "outline_raw": outline_raw,
            "centroid_indices": centroid_idx,
            "centroid_rows": rows_abs,
            "centroid_cols": cols_abs,
            "labels_at_centroids": labels_at_centroids,
        }
        self.tile_cache.set(cache_key, entry)
        return entry

    def get_overlay(self, params: Dict):
        try:
            if self._state.get("loading"):
                return self._error_response(409, "Dataset is still loading; try again in a moment.")
            overlay_type = params.get("overlay_type", "gene")
            tile_id = int(params.get("tile_id", 0))
            b2c_mode = params.get("b2c_mode", "fixed")
            max_bin_distance = float(params.get("max_bin_distance", 2.0) or 2.0)
            mpp = float(params.get("mpp", 0.3) or 0.3)
            bin_um = float(params.get("bin_um", 2.0) or 2.0)
            volume_ratio = float(params.get("volume_ratio", 4.0) or 4.0)
            pad_factor = int(params.get("pad_factor", 2) or 2)

            entry = self._tile_entry(
                tile_id,
                b2c_mode=b2c_mode,
                max_bin_distance=max_bin_distance,
                mpp=mpp,
                bin_um=bin_um,
                volume_ratio=volume_ratio,
                pad_factor=pad_factor,
            )

            if overlay_type == "gene":
                payload = self._prepare_gene_overlay(entry, params)
            elif overlay_type == "observation":
                payload = self._prepare_obs_overlay(entry, params)
            else:
                raise ValueError(f"Unknown overlay_type '{overlay_type}'")

            payload["tile"] = entry["tile"].to_dict()
            payload["dist_px"] = entry["dist_px"]
            payload["b2c_mode"] = b2c_mode
            payload["max_bin_distance"] = max_bin_distance
            payload["mpp"] = mpp
            payload["bin_um"] = bin_um
            payload["volume_ratio"] = volume_ratio
            payload["pad_factor"] = pad_factor

            return self._json_response(payload)
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to build overlay")

    def _prepare_gene_overlay(self, entry: Dict, params: Dict) -> Dict:
        ctx = self._ensure_context()
        genes_value = params.get("genes") or params.get("gene")
        if not genes_value:
            raise ValueError("Provide 'gene' or 'genes' for gene overlay.")
        if isinstance(genes_value, str):
            genes = [g.strip() for g in genes_value.split(",") if g.strip()]
        else:
            genes = list(genes_value)

        if not genes:
            raise ValueError("No gene names provided.")

        color_mode = params.get("color_mode", "gradient")
        gene_color = params.get("gene_color", "#ff6b6b") or "#ff6b6b"
        gradient_color = params.get("gradient_color")
        render_mode = params.get("render_mode", "fill")
        expr_quantile = params.get("expr_quantile")
        expr_quantile = float(expr_quantile) if expr_quantile not in (None, "") else None
        top_n = int(params.get("top_n", 0) or 0)
        overlay_alpha = float(params.get("overlay_alpha", 0.5) or 0.5)
        show_centroids = str(params.get("show_centroids", "false")).lower() in ("1", "true", "yes", "on")
        highlight_color = params.get("highlight_color", "#39ff14")
        highlight_width = float(params.get("highlight_width", 2.0) or 2.0)
        cmap_name = params.get("cmap_name", "viridis")
        vmin_param = params.get("vmin")
        vmax_param = params.get("vmax")
        vmin = float(vmin_param) if vmin_param not in (None, "") else None
        vmax = float(vmax_param) if vmax_param not in (None, "") else None
        all_expanded_outline = str(params.get("all_expanded_outline", "false")).lower() in ("1", "true", "yes", "on")
        all_nuclei_outline = str(params.get("all_nuclei_outline", "false")).lower() in ("1", "true", "yes", "on")

        try:
            solid_rgba = colors.to_rgba(gene_color)
        except ValueError:
            gene_color = "#ff6b6b"
            solid_rgba = colors.to_rgba(gene_color)

        custom_cmap = None
        if color_mode == "gradient" and gradient_color:
            try:
                target_rgba = colors.to_rgba(gradient_color)
            except ValueError:
                target_rgba = colors.to_rgba("#4285f4")
            custom_cmap = colors.LinearSegmentedColormap.from_list(
                "bin2cell_custom_gradient",
                [(0.0, (1.0, 1.0, 1.0, 0.0)), (1.0, target_rgba)],
            )
        else:
            target_rgba = None

        lab_exp: np.ndarray = entry["lab_exp"]
        polygons_exp: Dict[int, List[List[List[float]]]] = entry["polygons_exp"]
        centroid_indices: np.ndarray = entry["centroid_indices"]
        labels_at_centroids: np.ndarray = entry["labels_at_centroids"]

        overlays = []
        centroid_payload = []

        for gene in genes:
            values = ctx.gene_vector(gene)
            tile_values = values[centroid_indices] if centroid_indices.size else np.empty((0,), dtype=float)

            selected_mask = np.ones_like(tile_values, dtype=bool)
            threshold = None
            if expr_quantile is not None and tile_values.size:
                threshold = float(np.quantile(tile_values, expr_quantile))
                selected_mask &= tile_values > threshold

            selected_indices = centroid_indices[selected_mask]
            selected_values = tile_values[selected_mask]
            selected_labels = labels_at_centroids[selected_mask]

            label_expr: Dict[int, float] = {}
            for lbl, val in zip(selected_labels, selected_values):
                if lbl <= 0:
                    continue
                if (lbl not in label_expr) or (val > label_expr[lbl]):
                    label_expr[int(lbl)] = float(val)

            if top_n and len(label_expr) > top_n:
                top_labels = sorted(label_expr.keys(), key=lambda l: label_expr[l], reverse=True)[:top_n]
                label_expr = {lbl: label_expr[lbl] for lbl in top_labels}

            if not label_expr:
                overlays.append(
                    {
                        "gene": gene,
                        "features": [],
                        "legend": {
                            "type": "empty",
                            "gene": gene,
                            "message": "No labels passed filters in this tile.",
                        },
                    }
                )
                continue

            values_array = np.array(list(label_expr.values()), dtype=float)
            vmin_local = float(values_array.min()) if vmin is None else vmin
            vmax_local = float(values_array.max()) if vmax is None else vmax
            if vmax_local == vmin_local:
                vmax_local = vmin_local + 1e-6
            norm = colors.Normalize(vmin=vmin_local, vmax=vmax_local)

            features = []
            for lbl, expr_value in label_expr.items():
                polygons = polygons_exp.get(lbl)
                if not polygons:
                    continue
                if color_mode == "binary":
                    fill_color = _rgba_to_css(colors.to_rgba(highlight_color), overlay_alpha)
                    stroke_color = highlight_color
                elif color_mode == "solid":
                    fill_color = _rgba_to_css(solid_rgba, alpha_override=overlay_alpha)
                    stroke_color = gene_color
                else:
                    color_fraction = float(norm(expr_value))
                    cmap_source = custom_cmap or cmap_name
                    fill_color = _colormap_sample(cmap_source, color_fraction, alpha=overlay_alpha)
                    stroke_color = _colormap_sample(cmap_source, color_fraction, alpha=1.0)

                feature = {
                    "label": int(lbl),
                    "polygons": polygons,
                    "fill": fill_color if render_mode == "fill" else None,
                    "stroke": stroke_color,
                    "stroke_width": highlight_width,
                    "value": float(expr_value),
                }
                features.append(feature)

            legend = {
                "type": "binary" if color_mode == "binary" else ("solid" if color_mode == "solid" else "continuous"),
                "gene": gene,
                "color": gene_color if color_mode == "solid" else highlight_color,
                "overlay_alpha": overlay_alpha,
            }
            if color_mode == "gradient":
                cmap_source = custom_cmap or cmap_name
                legend.update(
                    {
                        "min": vmin_local,
                        "max": vmax_local,
                        "cmap": cmap_name,
                        "gradient": _sample_gradient(cmap_source),
                        "gradient_color": gradient_color,
                    }
                )

            overlays.append(
                {
                    "gene": gene,
                    "features": features,
                    "legend": legend,
                    "render_mode": render_mode,
                    "color_mode": color_mode,
                    "gene_color": gene_color,
                    "gradient_color": gradient_color,
                }
            )

            if show_centroids and centroid_indices.size:
                gene_centroids = []
                for lbl, idx_global, row, col, value in zip(
                    labels_at_centroids[selected_mask],
                    selected_indices,
                    entry["centroid_rows"][selected_mask],
                    entry["centroid_cols"][selected_mask],
                    selected_values,
                ):
                    if lbl <= 0:
                        continue
                    gene_centroids.append(
                        {
                            "index": int(idx_global),
                            "label": int(lbl),
                            "x": float(col),
                            "y": float(row),
                            "value": float(value),
                        }
                    )
                centroid_payload.append({"gene": gene, "points": gene_centroids})

        payload: Dict = {
            "overlay_type": "gene",
            "overlays": overlays,
            "all_expanded_outline": all_expanded_outline,
            "all_nuclei_outline": all_nuclei_outline,
            "expanded_outline": entry["outline_exp"] if all_expanded_outline else [],
            "nuclei_outline": entry["outline_raw"] if all_nuclei_outline else [],
        }
        if show_centroids:
            payload["centroids"] = centroid_payload
        return payload

    def _prepare_obs_overlay(self, entry: Dict, params: Dict) -> Dict:
        ctx = self._ensure_context()
        obs_col = params.get("obs_col")
        if not obs_col:
            raise ValueError("Provide 'obs_col' for observation overlay.")

        lab_exp: np.ndarray = entry["lab_exp"]
        polygons_exp: Dict[int, List[List[List[float]]]] = entry["polygons_exp"]
        centroid_indices: np.ndarray = entry["centroid_indices"]
        labels_at_centroids: np.ndarray = entry["labels_at_centroids"]

        category_filter = params.get("category")
        render_mode = params.get("render_mode", "fill")
        overlay_alpha = float(params.get("overlay_alpha", 0.5) or 0.5)
        cmap_name = params.get("cmap_name", "tab20")
        highlight_color = params.get("highlight_color", "#39ff14")
        highlight_width = float(params.get("highlight_width", 2.0) or 2.0)
        show_centroids = str(params.get("show_centroids", "false")).lower() in ("1", "true", "yes", "on")
        all_expanded_outline = str(params.get("all_expanded_outline", "true")).lower() in ("1", "true", "yes", "on")
        all_nuclei_outline = str(params.get("all_nuclei_outline", "false")).lower() in ("1", "true", "yes", "on")
        legend_outside = str(params.get("legend_outside", "true")).lower() in ("1", "true", "yes", "on")

        obs_vals = ctx.obs_vector(obs_col)
        tile_vals = obs_vals[centroid_indices] if centroid_indices.size else np.empty((0,), dtype=object)

        label_to_cat: Dict[int, str] = {}
        counts: Dict[int, Counter] = defaultdict(Counter)
        for lbl, val in zip(labels_at_centroids, tile_vals):
            if lbl <= 0:
                continue
            key = "" if val is None or (isinstance(val, float) and math.isnan(val)) else str(val)
            counts[int(lbl)][key] += 1
        for lbl, cnt in counts.items():
            if not cnt:
                continue
            label_to_cat[lbl] = cnt.most_common(1)[0][0]

        if not label_to_cat:
            return {
                "overlay_type": "observation",
                "obs_col": obs_col,
                "features": [],
                "legend": {"type": "empty", "message": "No categories found in tile."},
            }

        if category_filter:
            allowed_categories = {category_filter}
        else:
            allowed_categories = set(label_to_cat.values())

        categories_sorted = sorted(allowed_categories)
        cmap = cm.get_cmap(cmap_name, max(1, len(categories_sorted)))
        cat_to_color = {cat: _rgba_to_css(cmap(i / max(1, len(categories_sorted) - 1)), overlay_alpha) for i, cat in enumerate(categories_sorted)}

        features = []
        for lbl, cat in label_to_cat.items():
            if category_filter and cat != category_filter:
                continue
            polygons = polygons_exp.get(lbl)
            if not polygons:
                continue
            if render_mode == "fill":
                fill_color = cat_to_color[cat]
                stroke_color = "#222222"
            else:
                fill_color = None
                stroke_color = cat_to_color[cat]

            features.append(
                {
                    "label": int(lbl),
                    "category": cat,
                    "polygons": polygons,
                    "fill": fill_color,
                    "stroke": stroke_color,
                    "stroke_width": highlight_width,
                }
            )

        legend = {
            "type": "categorical",
            "obs_col": obs_col,
            "items": [{"label": cat, "color": cat_to_color[cat]} for cat in categories_sorted],
            "legend_outside": legend_outside,
        }

        payload: Dict = {
            "overlay_type": "observation",
            "obs_col": obs_col,
            "category_filter": category_filter,
            "features": features,
            "legend": legend,
            "render_mode": render_mode,
            "all_expanded_outline": all_expanded_outline,
            "all_nuclei_outline": all_nuclei_outline,
            "expanded_outline": entry["outline_exp"] if all_expanded_outline else [],
            "nuclei_outline": entry["outline_raw"] if all_nuclei_outline else [],
        }

        if show_centroids and centroid_indices.size:
            points = []
            for lbl, idx_global, row, col, obs_val in zip(
                labels_at_centroids,
                centroid_indices,
                entry["centroid_rows"],
                entry["centroid_cols"],
                tile_vals,
            ):
                if lbl <= 0:
                    continue
                label_str = label_to_cat.get(int(lbl))
                points.append(
                    {
                        "index": int(idx_global),
                        "label": int(lbl),
                        "category": label_str,
                        "x": float(col),
                        "y": float(row),
                    }
                )
            payload["centroids"] = [{"points": points}]

        return payload

    def export_overlay(self, params: Dict):
        try:
            response = self.get_overlay(params)
            if response.status_code and response.status_code >= 400:
                return response
            payload = json.loads(response.get_data(as_text=True))

            tile_id = payload["tile"]["id"]
            overlay_type = payload["overlay_type"]
            name = params.get("name")
            if not name:
                suffix = overlay_type
                if overlay_type == "gene" and payload.get("overlays"):
                    suffix = "_".join(ov["gene"] for ov in payload["overlays"])
                if overlay_type == "observation":
                    suffix = payload.get("obs_col", "observation")
                name = f"tile{tile_id}_{suffix}"

            safe_name = "".join(ch if ch.isalnum() or ch in ("_", "-", ".") else "_" for ch in name)
            out_path = os.path.join(self.out_dir, f"{safe_name}.geojson")

            features = []
            if overlay_type == "gene":
                for gene_overlay in payload.get("overlays", []):
                    gene_name = gene_overlay.get("gene")
                    for feature in gene_overlay.get("features", []):
                        geom = _feature_polygons_to_geojson(feature["polygons"])
                        props = {
                            "tile_id": tile_id,
                            "overlay_type": "gene",
                            "gene": gene_name,
                            "value": feature.get("value"),
                        }
                        features.append({"type": "Feature", "geometry": geom, "properties": props})
            else:
                obs_col = payload.get("obs_col")
                for feature in payload.get("features", []):
                    geom = _feature_polygons_to_geojson(feature["polygons"])
                    props = {
                        "tile_id": tile_id,
                        "overlay_type": "observation",
                        "obs_col": obs_col,
                        "category": feature.get("category"),
                    }
                    features.append({"type": "Feature", "geometry": geom, "properties": props})

            geojson = {"type": "FeatureCollection", "features": features}
            with open(out_path, "w") as fh:
                json.dump(geojson, fh)

            return self._json_response({"status": "ok", "path": out_path})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to export overlay")

    def save_preset(self, params: Dict):
        try:
            name = params.get("name")
            if not name:
                raise ValueError("Preset requires 'name'.")
            config = params.get("config")
            if not isinstance(config, dict):
                raise ValueError("Preset requires 'config' JSON object.")
            self.presets[name] = config
            self._save_presets()
            return self._json_response({"status": "ok", "name": name})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to save preset")

    def list_presets(self, params: Dict):
        try:
            return self._json_response({"presets": self.presets})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to list presets")

    def delete_preset(self, params: Dict):
        try:
            name = params.get("name")
            if not name or name not in self.presets:
                raise KeyError("Preset not found.")
            del self.presets[name]
            self._save_presets()
            return self._json_response({"status": "ok", "deleted": name})
        except Exception as exc:
            return self._error_from_exception(exc, "Failed to delete preset")


def _feature_polygons_to_geojson(polygons: List[List[List[float]]]) -> Dict:
    if not polygons:
        return {"type": "GeometryCollection", "geometries": []}
    if len(polygons) == 1:
        return {"type": "Polygon", "coordinates": [polygons[0]]}
    return {"type": "MultiPolygon", "coordinates": [[poly] for poly in polygons]}


def _json_default(obj):
    if isinstance(obj, Tile):
        return obj.to_dict()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")
class FileDialogHelper(QObject):
    def __init__(self) -> None:
        super().__init__()
        self._result: Tuple[str, str] = ("", "")

    @Slot(str, str, str)
    def _open_dialog(self, caption: str, start_dir: str, filters: str) -> None:
        self._result = QFileDialog.getOpenFileName(None, caption, start_dir, filters)

    def get_open_file_name(self, caption: str, start_dir: str, filters: str) -> Tuple[str, str]:
        self._result = ("", "")
        QMetaObject.invokeMethod(
            self,
            "_open_dialog",
            Qt.BlockingQueuedConnection,
            Q_ARG(str, caption),
            Q_ARG(str, start_dir),
            Q_ARG(str, filters),
        )
        return self._result
