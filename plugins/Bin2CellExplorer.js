var Bin2CellExplorer;
Bin2CellExplorer = {
  name: "Bin2Cell Explorer",
  parameters: {
    _sec_data: { label: "Dataset", type: "section", collapsed: false },
    he_path: { label: "H&E image (.tif/.tiff)", type: "text", default: "/Users/jjoseph/.tissuumaps/plugins/he.tiff" },
    he_path_browse: { label: "Browse H&E image…", type: "button" },
    labels_path: { label: "Label matrix (.npz)", type: "text", default: "/Users/jjoseph/.tissuumaps/plugins/he.npz" },
    labels_path_browse: { label: "Browse label matrix…", type: "button" },
    adata_path: { label: "AnnData (.h5ad)", type: "text", default: "/Users/jjoseph/.tissuumaps/plugins/P2_CRC_annotated.h5ad" },
    adata_path_browse: { label: "Browse AnnData…", type: "button" },
    obsm_key: { label: "obsm coord key", type: "select", default: "spatial_cropped_150_buffer" },
    tile_h: { label: "Tile height (px)", type: "number", default: 1500 },
    tile_w: { label: "Tile width (px)", type: "number", default: 1500 },
    stride_h: { label: "Stride height (px, optional)", type: "number", default: "" },
    stride_w: { label: "Stride width (px, optional)", type: "number", default: "" },
    load_dataset_btn: { label: "Load dataset", type: "button" },

    _sec_overlay: { label: "Overlay", type: "section", collapsed: false },
    tile_id: { label: "Tile ID", type: "select", default: "" },
    overlay_type: { label: "Overlay type", type: "select", default: "gene", options: ["gene", "observation"] },
    genes: { label: "Gene(s) (comma separated)", type: "text", default: "COL1A1" },
    obs_col: { label: "Observation column", type: "select", default: "" },
    category: { label: "Category filter (optional)", type: "select", default: "" },
    render_mode: { label: "Render mode", type: "select", default: "fill", options: ["fill", "outline"] },
    color_mode: { label: "Gene color mode", type: "select", default: "gradient", options: ["gradient", "binary", "solid"] },
    gradient_color: { label: "Gradient color (optional)", type: "text", default: "#4285f4", attributes: { type: "color" } },
    gene_color: { label: "Gene color (solid mode)", type: "text", default: "#ff6b6b", attributes: { type: "color" } },
    expr_quantile: { label: "Expr. quantile (0-1)", type: "number", default: "" },
    top_n: { label: "Top N labels", type: "number", default: "" },
    b2c_mode: { label: "Expand mode", type: "select", default: "fixed", options: ["fixed", "volume_ratio"] },
    max_bin_distance: { label: "Max bin distance", type: "number", default: 2.0 },
    mpp: { label: "Microns per pixel", type: "number", default: 0.3 },
    bin_um: { label: "Bin size (µm)", type: "number", default: 2.0 },
    volume_ratio: { label: "Volume ratio", type: "number", default: 4.0 },
    overlay_alpha: { label: "Overlay alpha (0-1)", type: "number", default: 0.5 },
    highlight_color: { label: "Highlight color", type: "text", default: "#39ff14" },
    highlight_width: { label: "Highlight width", type: "number", default: 2.0 },
    all_expanded_outline: { label: "Show expanded outlines", type: "checkbox", default: false },
    all_nuclei_outline: { label: "Show nuclei outlines", type: "checkbox", default: false },
    apply_overlay_btn: { label: "Update overlay", type: "button" },

    _sec_export: { label: "Export & Presets", type: "section", collapsed: true },
    export_name: { label: "Export name (optional)", type: "text", default: "" },
    export_geojson_btn: { label: "Export overlay to GeoJSON", type: "button" },
    save_preset_name: { label: "Preset name", type: "text", default: "" },
    save_preset_btn: { label: "Save preset", type: "button" },
    preset_select: { label: "Presets", type: "select", default: "" },
    load_preset_btn: { label: "Apply preset", type: "button" },
    delete_preset_btn: { label: "Delete preset", type: "button" }
  }
};

Bin2CellExplorer.state = {
  datasetLoaded: false,
  tiles: [],
  obsColumns: [],
  obsMetadata: {},
  obsMetadataRequests: {},
  obsCategorySelections: {},
  pendingCategoryValue: "",
  selectedObsCol: null,
  genesPreview: [],
  overlays: null,
  layers: {},
  presets: {},
  slideShape: null,
  tileOverviewCanvas: null,
  tileOverviewScale: null,
  selectedTileId: null
};

Bin2CellExplorer.init = function(container) {
  Bin2CellExplorer.state.container = container;
  container.classList.add("bin2cell-explorer-panel");

  const legend = document.createElement("div");
  legend.id = "Bin2CellExplorer_legend";
  legend.className = "bin2cell-legend mt-2";
  container.appendChild(legend);

  const status = document.createElement("div");
  status.id = "Bin2CellExplorer_status";
  status.className = "bin2cell-status mt-2 small text-muted";
  container.appendChild(status);

  const overviewBlock = document.createElement("div");
  overviewBlock.id = "Bin2CellExplorer_tile_overview_block";
  overviewBlock.className = "tile-overview mt-3";
  const overviewTitle = document.createElement("div");
  overviewTitle.className = "small fw-bold mb-1";
  overviewTitle.textContent = "Tile overview";
  const overviewCanvas = document.createElement("canvas");
  overviewCanvas.id = "Bin2CellExplorer_tile_overview";
  overviewCanvas.width = 320;
  overviewCanvas.height = 320;
  overviewCanvas.style.border = "1px solid rgba(0,0,0,0.1)";
  overviewCanvas.style.borderRadius = "12px";
  overviewCanvas.style.background = "linear-gradient(135deg, #fdfbfb, #ebedee)";
  overviewCanvas.style.boxShadow = "0 6px 18px rgba(15,23,42,0.15)";
  overviewCanvas.style.cursor = "pointer";
  const overviewHint = document.createElement("div");
  overviewHint.className = "small text-muted mt-1";
  overviewHint.id = "Bin2CellExplorer_tile_overview_hint";
  overviewHint.textContent = "Load dataset to display tiles";
  overviewBlock.appendChild(overviewTitle);
  overviewBlock.appendChild(overviewCanvas);
  overviewBlock.appendChild(overviewHint);
  container.appendChild(overviewBlock);

  Bin2CellExplorer.state.tileOverviewCanvas = overviewCanvas;
  overviewCanvas.addEventListener("click", Bin2CellExplorer.onTileOverviewClick);

  interfaceUtils.alert("Bin2Cell Explorer loaded");
  Bin2CellExplorer.toggleOverlayInputs("gene");
};

Bin2CellExplorer.inputTrigger = function(inputName) {
  switch (inputName) {
    case "load_dataset_btn":
      Bin2CellExplorer.loadDataset();
      break;
    case "he_path_browse":
      Bin2CellExplorer.browseForFile("he_path");
      break;
    case "labels_path_browse":
      Bin2CellExplorer.browseForFile("labels_path");
      break;
    case "adata_path_browse":
      Bin2CellExplorer.browseForFile("adata_path");
      break;
    case "overlay_type":
      Bin2CellExplorer.toggleOverlayInputs(Bin2CellExplorer.get("overlay_type"));
      break;
    case "color_mode":
      Bin2CellExplorer.updateGeneColorControls();
      break;
    case "apply_overlay_btn":
      Bin2CellExplorer.requestOverlay();
      break;
    case "export_geojson_btn":
      Bin2CellExplorer.exportOverlay();
      break;
    case "save_preset_btn":
      Bin2CellExplorer.savePreset();
      break;
    case "load_preset_btn":
      Bin2CellExplorer.applySelectedPreset();
      break;
    case "delete_preset_btn":
      Bin2CellExplorer.deleteSelectedPreset();
      break;
    case "preset_select":
      Bin2CellExplorer.populatePresetPreview();
      break;
    default:
      break;
  }
};

Bin2CellExplorer.toggleOverlayInputs = function(mode) {
  const geneIds = ["genes", "color_mode", "gradient_color", "gene_color", "expr_quantile", "top_n"];
  const obsIds = ["obs_col", "category"];
  geneIds.forEach((id) => Bin2CellExplorer.toggleParam(id, mode === "gene"));
  obsIds.forEach((id) => Bin2CellExplorer.toggleParam(id, mode === "observation"));
  if (mode === "gene") {
    Bin2CellExplorer.updateGeneColorControls();
  } else if (mode === "observation" && Bin2CellExplorer.state.datasetLoaded) {
    Bin2CellExplorer.populateCategorySelect(Bin2CellExplorer.get("obs_col") || "");
  }
};

Bin2CellExplorer.toggleParam = function(name, visible) {
  const domId = "Bin2CellExplorer_" + name;
  const element = document.getElementById(domId);
  if (!element) return;
  let wrapper = element.closest(".form-group, .row, .input-group");
  if (!wrapper) {
    wrapper = element.parentElement;
  }
  if (wrapper) wrapper.style.display = visible ? "" : "none";
};

Bin2CellExplorer.ensureObject = function(payload) {
  if (!payload) return {};
  if (typeof payload === "string") {
    try {
      return JSON.parse(payload);
    } catch (err) {
      console.warn("Bin2CellExplorer: failed to parse payload", err);
      return {};
    }
  }
  return payload;
};

Bin2CellExplorer.setStatus = function(msg) {
  const status = document.getElementById("Bin2CellExplorer_status");
  if (status) status.textContent = msg || "";
};

Bin2CellExplorer.updateGeneColorControls = function() {
  const mode = Bin2CellExplorer.get("color_mode");
  Bin2CellExplorer.toggleParam("gradient_color", mode === "gradient");
  Bin2CellExplorer.toggleParam("gene_color", mode === "solid");
};

Bin2CellExplorer.browseForFile = function(field) {
  const current = Bin2CellExplorer.get(field) || "";
  const payload = {
    field: field,
    current_path: current,
  };
  Bin2CellExplorer.api(
    "pick_file",
    payload,
    function(resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      if (!data || data.status === "cancelled") {
        Bin2CellExplorer.setStatus("File selection cancelled.");
        return;
      }
      if (data.status !== "ok" || !data.path) {
        Bin2CellExplorer.setStatus("File selection failed.");
        return;
      }
      interfaceUtils.setValueForElement("Bin2CellExplorer_" + field, "value", data.path);
      Bin2CellExplorer.setStatus("Selected " + data.path);
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.loadDataset = function() {
  Bin2CellExplorer.setStatus("Loading dataset…");
  const payload = {
    he_path: Bin2CellExplorer.get("he_path"),
    labels_path: Bin2CellExplorer.get("labels_path"),
    adata_path: Bin2CellExplorer.get("adata_path"),
    obsm_key: Bin2CellExplorer.get("obsm_key"),
    tile_h: Bin2CellExplorer.get("tile_h"),
    tile_w: Bin2CellExplorer.get("tile_w"),
    stride_h: Bin2CellExplorer.get("stride_h"),
    stride_w: Bin2CellExplorer.get("stride_w")
  };
  Bin2CellExplorer.api(
    "load_dataset",
    payload,
    function(resp) {
      Bin2CellExplorer.onDatasetLoaded(Bin2CellExplorer.ensureObject(resp));
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.onDatasetLoaded = function(data) {
  Bin2CellExplorer.state.datasetLoaded = true;
  Bin2CellExplorer.state.tiles = data.tiles || [];
  Bin2CellExplorer.state.obsColumns = data.obs_columns || [];
  Bin2CellExplorer.state.obsMetadata = {};
  Bin2CellExplorer.state.obsMetadataRequests = {};
  Bin2CellExplorer.state.obsCategorySelections = {};
  Bin2CellExplorer.state.pendingCategoryValue = "";
  Bin2CellExplorer.state.selectedObsCol = null;
  Bin2CellExplorer.state.genesPreview = data.genes_preview || [];
  Bin2CellExplorer.state.datasetId = data.dataset_id;
  Bin2CellExplorer.state.hePath = data.he_path;
  Bin2CellExplorer.state.slideShape = data.shape;

  Bin2CellExplorer.populateSelect("tile_id", Bin2CellExplorer.state.tiles.map(function(tile) { return tile.id; }));
  Bin2CellExplorer.populateSelect("obs_col", Bin2CellExplorer.state.obsColumns);
  Bin2CellExplorer.populateSelect("category", [""]);
  Bin2CellExplorer.populateSelect("obsm_key", data.available_obsm || [], data.obsm_key);

  if (Bin2CellExplorer.state.tiles.length) {
    Bin2CellExplorer.set("tile_id", String(Bin2CellExplorer.state.tiles[0].id));
    Bin2CellExplorer.state.selectedTileId = Bin2CellExplorer.state.tiles[0].id;
  }
  if (Bin2CellExplorer.state.obsColumns.length) {
    const firstObs = Bin2CellExplorer.state.obsColumns[0];
    Bin2CellExplorer.set("obs_col", firstObs);
    Bin2CellExplorer.onObsColumnChange(firstObs);
  } else {
    Bin2CellExplorer.populateCategorySelect("");
  }

  Bin2CellExplorer.requestPresets();

  Bin2CellExplorer.renderTileOverview();
  Bin2CellExplorer.attachTileSelectListener();
  Bin2CellExplorer.attachObsSelectListener();
  Bin2CellExplorer.attachCategorySelectListener();

  Bin2CellExplorer.setStatus("Dataset loaded (" + Bin2CellExplorer.state.tiles.length + " tiles)");
};

Bin2CellExplorer.populateSelect = function(name, options, selected) {
  const domId = "Bin2CellExplorer_" + name;
  interfaceUtils.cleanSelect(domId);
  options.forEach(function(opt) {
    interfaceUtils.addSingleElementToSelect(domId, String(opt));
  });
  if (selected !== undefined && selected !== null && selected !== "") {
    interfaceUtils.setValueForElement(domId, "value", String(selected));
  }
  if (name === "category") {
    const selectEl = document.getElementById(domId);
    if (selectEl && selectEl.options.length) {
      for (let i = 0; i < selectEl.options.length; i += 1) {
        if (selectEl.options[i].value === "") {
          selectEl.options[i].text = selectEl.options[i].text || "All categories";
          if (!selectEl.options[i].text.trim()) {
            selectEl.options[i].text = "All categories";
          }
          break;
        }
      }
    }
  }
};

Bin2CellExplorer.attachTileSelectListener = function() {
  const select = document.getElementById("Bin2CellExplorer_tile_id");
  if (!select || select.__bin2cell_tile_listener) return;
  select.addEventListener("change", function() {
    const value = parseInt(select.value, 10);
    if (!isNaN(value)) {
      Bin2CellExplorer.state.selectedTileId = value;
      Bin2CellExplorer.renderTileOverview();
      Bin2CellExplorer.panToTile(value, true);
    }
  });
  select.__bin2cell_tile_listener = true;
};

Bin2CellExplorer.attachObsSelectListener = function() {
  const select = document.getElementById("Bin2CellExplorer_obs_col");
  if (!select || select.__bin2cell_obs_listener) return;
  select.addEventListener("change", function() {
    Bin2CellExplorer.onObsColumnChange(select.value || "");
  });
  select.__bin2cell_obs_listener = true;
};

Bin2CellExplorer.attachCategorySelectListener = function() {
  const select = document.getElementById("Bin2CellExplorer_category");
  if (!select || select.__bin2cell_category_listener) return;
  select.addEventListener("change", function() {
    const col = Bin2CellExplorer.get("obs_col");
    if (col) {
      Bin2CellExplorer.state.obsCategorySelections[col] = select.value || "";
    }
  });
  select.__bin2cell_category_listener = true;
};

Bin2CellExplorer.onObsColumnChange = function(column, options) {
  if (!Bin2CellExplorer.state.datasetLoaded) return;
  const value = column || "";
  Bin2CellExplorer.state.selectedObsCol = value || null;
  const opts = options || {};
  if (!opts.keepCategory) {
    const cached = value ? Bin2CellExplorer.state.obsCategorySelections[value] : "";
    Bin2CellExplorer.state.pendingCategoryValue = cached || "";
    interfaceUtils.setValueForElement("Bin2CellExplorer_category", "value", "");
  }
  Bin2CellExplorer.ensureObsMetadata(value);
};

Bin2CellExplorer.ensureObsMetadata = function(column) {
  const col = column || "";
  Bin2CellExplorer.populateCategorySelect(col);
  if (!col) return;
  if (Bin2CellExplorer.state.obsMetadata[col]) {
    Bin2CellExplorer.populateCategorySelect(col);
    return;
  }
  if (Bin2CellExplorer.state.obsMetadataRequests[col]) {
    return;
  }
  Bin2CellExplorer.state.obsMetadataRequests[col] = true;
  Bin2CellExplorer.api(
    "describe_obs_column",
    { obs_col: col },
    function(resp) {
      delete Bin2CellExplorer.state.obsMetadataRequests[col];
      const data = Bin2CellExplorer.ensureObject(resp) || {};
      Bin2CellExplorer.state.obsMetadata[col] = data;
      Bin2CellExplorer.populateCategorySelect(col);
    },
    function(jqXHR, textStatus, errorThrown) {
      delete Bin2CellExplorer.state.obsMetadataRequests[col];
      Bin2CellExplorer.handleError(jqXHR, textStatus, errorThrown);
    }
  );
};

Bin2CellExplorer.populateCategorySelect = function(column) {
  const metadata = column ? Bin2CellExplorer.state.obsMetadata[column] : null;
  const categories = (metadata && Array.isArray(metadata.categories)) ? metadata.categories.slice() : [];
  const options = [""].concat(categories);
  let desired = Bin2CellExplorer.state.pendingCategoryValue;
  if (!desired) {
    desired = Bin2CellExplorer.get("category") || "";
  }
  if (options.indexOf(desired) === -1) {
    desired = "";
  }
  Bin2CellExplorer.populateSelect("category", options, desired);
  if (desired && Bin2CellExplorer.state.pendingCategoryValue === desired) {
    Bin2CellExplorer.state.pendingCategoryValue = "";
  }
  if (metadata && Array.isArray(metadata.categories) && column) {
    if (desired) {
      Bin2CellExplorer.state.obsCategorySelections[column] = desired;
    } else {
      delete Bin2CellExplorer.state.obsCategorySelections[column];
    }
  }
};

Bin2CellExplorer.requestOverlay = function() {
  if (!Bin2CellExplorer.state.datasetLoaded) {
    interfaceUtils.alert("Load a dataset first.");
    return;
  }
  const overlayType = Bin2CellExplorer.get("overlay_type") || "gene";
  const payload = Bin2CellExplorer.collectOverlayParams();
  payload.overlay_type = overlayType;

  Bin2CellExplorer.setStatus("Computing overlay…");
  Bin2CellExplorer.api(
    "get_overlay",
    payload,
    function(resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      Bin2CellExplorer.renderOverlay(data);
      Bin2CellExplorer.renderLegend(data);
      Bin2CellExplorer.setStatus("Overlay ready (" + data.overlay_type + ")");
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.renderTileOverview = function() {
  const canvas = Bin2CellExplorer.state.tileOverviewCanvas;
  const hint = document.getElementById("Bin2CellExplorer_tile_overview_hint");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const tiles = Bin2CellExplorer.state.tiles || [];
  const shape = Bin2CellExplorer.state.slideShape;
  if (!tiles.length || !shape) {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    if (hint) hint.textContent = "Load dataset to display tiles";
    Bin2CellExplorer.state.tileOverviewScale = null;
    return;
  }

  const slideHeight = shape[0];
  const slideWidth = shape[1];
  const margin = 16;
  const maxWidth = 320;
  const maxHeight = 320;
  const scale = Math.min((maxWidth - margin * 2) / slideWidth, (maxHeight - margin * 2) / slideHeight);
  const canvasWidth = Math.ceil(slideWidth * scale + margin * 2);
  const canvasHeight = Math.ceil(slideHeight * scale + margin * 2);
  canvas.width = canvasWidth;
  canvas.height = canvasHeight;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  const gradient = ctx.createLinearGradient(0, 0, canvas.width, canvas.height);
  gradient.addColorStop(0, "#fdfbfb");
  gradient.addColorStop(1, "#ebedee");
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  const selectedId = Bin2CellExplorer.state.selectedTileId;
  ctx.lineWidth = 1;
  ctx.strokeStyle = "rgba(60,60,60,0.5)";
  ctx.font = "10px sans-serif";
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";

  tiles.forEach(function(tile) {
    const x = margin + tile.c0 * scale;
    const y = margin + tile.r0 * scale;
    const w = (tile.c1 - tile.c0) * scale;
    const h = (tile.r1 - tile.r0) * scale;
    if (tile.id === selectedId) {
      ctx.fillStyle = "rgba(255, 193, 7, 0.45)";
      ctx.fillRect(x, y, w, h);
    }
    ctx.strokeStyle = "rgba(60,60,60,0.6)";
    ctx.strokeRect(x, y, w, h);
    if (w > 10 && h > 10) {
      ctx.fillStyle = "#202124";
      ctx.fillText(String(tile.id), x + w / 2, y + h / 2);
    }
  });

  if (hint) {
    hint.textContent = "Click a tile to select and center";
  }
  Bin2CellExplorer.state.tileOverviewScale = { scale: scale, margin: margin };
};

Bin2CellExplorer.onTileOverviewClick = function(evt) {
  const canvas = Bin2CellExplorer.state.tileOverviewCanvas;
  const scaleInfo = Bin2CellExplorer.state.tileOverviewScale;
  if (!canvas || !scaleInfo || !Bin2CellExplorer.state.tiles.length) return;
  const rect = canvas.getBoundingClientRect();
  const x = (evt.clientX - rect.left);
  const y = (evt.clientY - rect.top);
  const slideX = (x - scaleInfo.margin) / scaleInfo.scale;
  const slideY = (y - scaleInfo.margin) / scaleInfo.scale;
  if (slideX < 0 || slideY < 0) return;
  const tile = Bin2CellExplorer.state.tiles.find(function(t) {
    return slideX >= t.c0 && slideX < t.c1 && slideY >= t.r0 && slideY < t.r1;
  });
  if (!tile) return;
  Bin2CellExplorer.state.selectedTileId = tile.id;
  Bin2CellExplorer.set("tile_id", String(tile.id));
  Bin2CellExplorer.renderTileOverview();
  Bin2CellExplorer.panToTile(tile.id, true);
};

Bin2CellExplorer.panToTile = function(tileId, animate) {
  const tile = Bin2CellExplorer.state.tiles.find(function(t) { return t.id === tileId; });
  const shape = Bin2CellExplorer.state.slideShape;
  if (!tile || !shape) return;
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer || !viewer.world.getItemCount()) return;
  const image = viewer.world.getItemAt(0);
  const width = shape[1];
  const height = shape[0];
  const xNorm = tile.c0 / width;
  const yNorm = tile.r0 / height;
  const wNorm = (tile.c1 - tile.c0) / width;
  const hNorm = (tile.r1 - tile.r0) / height;
  const rect = new OpenSeadragon.Rect(xNorm, yNorm, wNorm, hNorm);
  viewer.viewport.fitBounds(rect, animate !== false);
};

Bin2CellExplorer.collectOverlayParams = function() {
  const params = {
    tile_id: Bin2CellExplorer.get("tile_id"),
    overlay_type: Bin2CellExplorer.get("overlay_type"),
    genes: Bin2CellExplorer.get("genes"),
    obs_col: Bin2CellExplorer.get("obs_col"),
    category: Bin2CellExplorer.get("category"),
    render_mode: Bin2CellExplorer.get("render_mode"),
    color_mode: Bin2CellExplorer.get("color_mode"),
    gradient_color: Bin2CellExplorer.get("gradient_color"),
    gene_color: Bin2CellExplorer.get("gene_color"),
    expr_quantile: Bin2CellExplorer.get("expr_quantile"),
    top_n: Bin2CellExplorer.get("top_n"),
    b2c_mode: Bin2CellExplorer.get("b2c_mode"),
    max_bin_distance: Bin2CellExplorer.get("max_bin_distance"),
    mpp: Bin2CellExplorer.get("mpp"),
    bin_um: Bin2CellExplorer.get("bin_um"),
    volume_ratio: Bin2CellExplorer.get("volume_ratio"),
    overlay_alpha: Bin2CellExplorer.get("overlay_alpha"),
    highlight_color: Bin2CellExplorer.get("highlight_color"),
    highlight_width: Bin2CellExplorer.get("highlight_width"),
    all_expanded_outline: Bin2CellExplorer.isChecked("all_expanded_outline"),
    all_nuclei_outline: Bin2CellExplorer.isChecked("all_nuclei_outline")
  };
  return params;
};

Bin2CellExplorer.exportOverlay = function() {
  if (!Bin2CellExplorer.state.datasetLoaded) {
    interfaceUtils.alert("Load a dataset first.");
    return;
  }
  const params = Bin2CellExplorer.collectOverlayParams();
  params.name = Bin2CellExplorer.get("export_name");
  Bin2CellExplorer.api(
    "export_overlay",
    params,
    function(resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      Bin2CellExplorer.setStatus("GeoJSON written to " + data.path);
      interfaceUtils.alert("Overlay exported:\n" + data.path);
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.collectPresetConfig = function() {
  return Bin2CellExplorer.collectOverlayParams();
};

Bin2CellExplorer.savePreset = function() {
  const name = (Bin2CellExplorer.get("save_preset_name") || "").trim();
  if (!name) {
    interfaceUtils.alert("Provide a preset name.");
    return;
  }
  const payload = {
    name: name,
    config: Bin2CellExplorer.collectPresetConfig()
  };
  Bin2CellExplorer.api(
    "save_preset",
    payload,
    function(resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      Bin2CellExplorer.setStatus("Preset saved.");
      Bin2CellExplorer.requestPresets(name);
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.requestPresets = function(selectName) {
  Bin2CellExplorer.api(
    "list_presets",
    {},
    function(resp) {
      const data = Bin2CellExplorer.ensureObject(resp);
      Bin2CellExplorer.state.presets = data.presets || {};
      const names = Object.keys(Bin2CellExplorer.state.presets);
      Bin2CellExplorer.populateSelect("preset_select", [""].concat(names), selectName || "");
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.applySelectedPreset = function() {
  const select = Bin2CellExplorer.get("preset_select");
  if (!select || !(select in Bin2CellExplorer.state.presets)) {
    interfaceUtils.alert("Select a preset first.");
    return;
  }
  const config = Bin2CellExplorer.state.presets[select];
  if (!config) return;
  let pendingCategory = null;
  let pendingObsCol = null;
  Object.keys(config).forEach(function(key) {
    if (key === "category") {
      pendingCategory = config[key] || "";
      return;
    }
    const domId = "Bin2CellExplorer_" + key;
    if (document.getElementById(domId)) {
      interfaceUtils.setValueForElement(domId, "value", config[key]);
      if (key === "obs_col") {
        pendingObsCol = config[key] || "";
      }
    }
  });
  Bin2CellExplorer.toggleOverlayInputs(config.overlay_type || Bin2CellExplorer.get("overlay_type"));
  if (pendingCategory !== null) {
    Bin2CellExplorer.state.pendingCategoryValue = pendingCategory;
  }
  if (pendingObsCol !== null) {
    Bin2CellExplorer.onObsColumnChange(pendingObsCol, { keepCategory: pendingCategory !== null });
  } else if (pendingCategory !== null) {
    Bin2CellExplorer.populateCategorySelect(Bin2CellExplorer.get("obs_col") || "");
  }
  Bin2CellExplorer.setStatus("Preset applied: " + select);
};

Bin2CellExplorer.populatePresetPreview = function() {
  const select = Bin2CellExplorer.get("preset_select");
  if (!select) return;
  const preset = Bin2CellExplorer.state.presets[select];
  if (!preset) return;
  Bin2CellExplorer.setStatus("Preset '" + select + "' ready to apply.");
};

Bin2CellExplorer.deleteSelectedPreset = function() {
  const select = Bin2CellExplorer.get("preset_select");
  if (!select || !(select in Bin2CellExplorer.state.presets)) {
    interfaceUtils.alert("Select a preset to delete.");
    return;
  }
  Bin2CellExplorer.api(
    "delete_preset",
    { name: select },
    function() {
      Bin2CellExplorer.setStatus("Preset deleted: " + select);
      Bin2CellExplorer.requestPresets();
    },
    Bin2CellExplorer.handleError
  );
};

Bin2CellExplorer.handleError = function(jqXHR, textStatus, errorThrown) {
  let message = "";
  if (jqXHR) {
    if (jqXHR.responseJSON) {
      message = jqXHR.responseJSON.message || JSON.stringify(jqXHR.responseJSON);
    } else if (jqXHR.responseText) {
      try {
        const parsed = JSON.parse(jqXHR.responseText);
        message = parsed.message || jqXHR.responseText;
      } catch (parseErr) {
        message = jqXHR.responseText;
      }
    }
  }
  if (!message && errorThrown) {
    message = errorThrown;
  }
  if (!message && textStatus) {
    message = textStatus;
  }
  if (!message) {
    message = "Unknown error";
  }
  Bin2CellExplorer.setStatus("Error: " + message);
  interfaceUtils.alert("Bin2Cell Explorer error:\n" + message);
};

Bin2CellExplorer.isChecked = function(name) {
  return interfaceUtils.isChecked("Bin2CellExplorer_" + name);
};

Bin2CellExplorer.renderOverlay = function(data) {
  Bin2CellExplorer.state.overlays = data;
  const layer = Bin2CellExplorer.ensureSvgLayer();
  if (!layer) return;

  layer.selectAll("*").remove();
  Bin2CellExplorer.state.layers = {
    gene: {},
    categories: {},
    outlines: {}
  };

  if (data.overlay_type === "gene") {
    (data.overlays || []).forEach(function(geneOverlay) {
      Bin2CellExplorer.drawGeneOverlay(layer, geneOverlay);
    });
  } else if (data.overlay_type === "observation") {
    Bin2CellExplorer.drawObservationOverlay(layer, data);
  }

  Bin2CellExplorer.drawOutlines(layer, data);
  Bin2CellExplorer.renderTileOverview();
};

Bin2CellExplorer.ensureSvgLayer = function() {
  const op = tmapp["object_prefix"];
  const base = overlayUtils._d3nodes[op + "_svgnode"];
  if (!base) {
    interfaceUtils.alert("Viewer not ready yet.");
    return null;
  }
  if (!Bin2CellExplorer.state.d3layer) {
    Bin2CellExplorer.state.d3layer = base.append("g").attr("id", op + "_bin2cell_layer");
  }
  return Bin2CellExplorer.state.d3layer;
};


Bin2CellExplorer.drawGeneOverlay = function(layer, overlay) {
  const gene = overlay.gene;
  const group = layer.append("g").attr("class", "bin2cell-gene-layer").attr("data-gene", gene);
  const renderMode = overlay.render_mode || "fill";

  (overlay.features || []).forEach(function(feature) {
    const d = Bin2CellExplorer.polygonsToPath(feature.polygons || []);
    if (!d) return;
    const path = group.append("path")
      .attr("class", "bin2cell-gene-feature")
      .attr("d", d)
      .attr("stroke", feature.stroke || "none")
      .attr("stroke-width", feature.stroke_width || 1.0)
      .attr("fill", renderMode === "fill" ? (feature.fill || "none") : "none")
      .attr("vector-effect", "non-scaling-stroke")
      .attr("data-label", feature.label)
      .attr("data-value", feature.value);
  });

  Bin2CellExplorer.state.layers.gene[gene] = group;
};

Bin2CellExplorer.drawObservationOverlay = function(layer, data) {
  const renderMode = data.render_mode || "fill";
  const features = data.features || [];
  const categoryGroups = {};

  features.forEach(function(feature) {
    const category = feature.category || "__uncategorized__";
    if (!categoryGroups[category]) {
      categoryGroups[category] = layer.append("g")
        .attr("class", "bin2cell-category-layer")
        .attr("data-category", category);
    }
    const d = Bin2CellExplorer.polygonsToPath(feature.polygons || []);
    if (!d) return;
    categoryGroups[category].append("path")
      .attr("d", d)
      .attr("stroke", feature.stroke || "none")
      .attr("stroke-width", feature.stroke_width || 1.0)
      .attr("fill", renderMode === "fill" ? (feature.fill || "none") : "none")
      .attr("vector-effect", "non-scaling-stroke")
      .attr("data-label", feature.label);
  });

  Bin2CellExplorer.state.layers.categories = categoryGroups;
};

Bin2CellExplorer.drawOutlines = function(layer, data) {
  const outlinesGroup = layer.append("g").attr("class", "bin2cell-outline-layer");

  (data.expanded_outline || []).forEach(function(line) {
    const pathData = Bin2CellExplorer.lineToPath(line);
    if (!pathData) return;
    outlinesGroup.append("path")
      .attr("d", pathData)
      .attr("stroke", "rgba(180,180,180,0.7)")
      .attr("stroke-width", 1.0)
      .attr("fill", "none")
      .attr("vector-effect", "non-scaling-stroke")
      .attr("data-kind", "expanded");
  });

  (data.nuclei_outline || []).forEach(function(line) {
    const pathData = Bin2CellExplorer.lineToPath(line);
    if (!pathData) return;
    outlinesGroup.append("path")
      .attr("d", pathData)
      .attr("stroke", "rgba(160,160,160,0.5)")
      .attr("stroke-width", 0.8)
      .attr("fill", "none")
      .attr("vector-effect", "non-scaling-stroke")
      .attr("data-kind", "nuclei");
  });

  Bin2CellExplorer.state.layers.outlines.group = outlinesGroup;
};


Bin2CellExplorer.imageToViewport = function(x, y, tiledImage) {
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer) return { x: x, y: y };
  const image = tiledImage || viewer.world.getItemAt(0);
  if (!image) return { x: x, y: y };
  let px = x;
  if (image.getFlip && image.getFlip()) {
    px = image.getContentSize().x - px;
  }
  const point = image.imageToViewportCoordinates(px, y, true);
  return point;
};

Bin2CellExplorer.polygonsToPath = function(polygons) {
  if (!polygons || !polygons.length) return "";
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer || !viewer.world.getItemCount()) return "";
  const image = viewer.world.getItemAt(0);
  const segments = [];
  polygons.forEach(function(poly) {
    if (!poly || poly.length < 3) return;
    let path = "";
    poly.forEach(function(pt, idx) {
      const vp = Bin2CellExplorer.imageToViewport(pt[0], pt[1], image);
      path += (idx === 0 ? "M" : "L") + vp.x + " " + vp.y;
    });
    path += "Z";
    segments.push(path);
  });
  return segments.join(" ");
};

Bin2CellExplorer.lineToPath = function(points) {
  if (!points || !points.length) return "";
  const viewer = tmapp[tmapp["object_prefix"] + "_viewer"];
  if (!viewer || !viewer.world.getItemCount()) return "";
  const image = viewer.world.getItemAt(0);
  let path = "";
  points.forEach(function(pt, idx) {
    const vp = Bin2CellExplorer.imageToViewport(pt[0], pt[1], image);
    path += (idx === 0 ? "M" : "L") + vp.x + " " + vp.y;
  });
  return path;
};

Bin2CellExplorer.renderLegend = function(data) {
  const legend = document.getElementById("Bin2CellExplorer_legend");
  if (!legend) return;
  legend.innerHTML = "";

  if (data.overlay_type === "gene") {
    Bin2CellExplorer.buildGeneLegend(legend, data.overlays || []);
  } else if (data.overlay_type === "observation") {
    Bin2CellExplorer.buildObservationLegend(legend, data.legend || {});
  }
};

Bin2CellExplorer.buildGeneLegend = function(container, overlays) {
  if (!overlays.length) {
    container.textContent = "No genes in overlay.";
    return;
  }
  const title = document.createElement("div");
  title.textContent = "Genes";
  title.className = "fw-bold mb-1";
  container.appendChild(title);

  overlays.forEach(function(overlay) {
    const gene = overlay.gene;
    const legend = overlay.legend || {};

    const row = document.createElement("div");
    row.className = "bin2cell-legend-row";

    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = true;
    checkbox.dataset.gene = gene;
    checkbox.addEventListener("change", function(evt) {
      Bin2CellExplorer.toggleGeneLayer(gene, evt.target.checked);
    });

    const label = document.createElement("label");
    label.textContent = " " + gene;
    label.className = "ms-1";

    row.appendChild(checkbox);
    row.appendChild(label);

    if (legend.type === "continuous" && legend.gradient) {
      const gradient = document.createElement("div");
      gradient.className = "bin2cell-gradient";
      gradient.style.height = "12px";
      gradient.style.flex = "1";
      gradient.style.marginLeft = "8px";
      gradient.style.borderRadius = "4px";
      gradient.style.background = "linear-gradient(to right," + legend.gradient.map(function(stop) {
        return stop[1];
      }).join(",") + ")";
      gradient.title = "min: " + legend.min + " max: " + legend.max;
      row.appendChild(gradient);
    } else if (legend.type === "solid" && legend.color) {
      const swatch = document.createElement("span");
      swatch.className = "bin2cell-swatch";
      swatch.style.display = "inline-block";
      swatch.style.width = "16px";
      swatch.style.height = "16px";
      swatch.style.marginLeft = "8px";
      swatch.style.borderRadius = "3px";
      swatch.style.border = "1px solid rgba(0,0,0,0.2)";
      swatch.style.background = legend.color;
      row.appendChild(swatch);
    }
    container.appendChild(row);
  });
};

Bin2CellExplorer.buildObservationLegend = function(container, legendData) {
  const items = legendData.items || [];
  if (!items.length) {
    container.textContent = "No observation categories.";
    return;
  }
  const title = document.createElement("div");
  title.textContent = legendData.obs_col || "Observation";
  title.className = "fw-bold mb-1";
  container.appendChild(title);

  items.forEach(function(item) {
    const row = document.createElement("div");
    row.className = "bin2cell-legend-row";

    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = true;
    checkbox.dataset.category = item.label;
    checkbox.addEventListener("change", function(evt) {
      Bin2CellExplorer.toggleCategoryLayer(item.label, evt.target.checked);
    });

    const swatch = document.createElement("span");
    swatch.className = "bin2cell-swatch";
    swatch.style.display = "inline-block";
    swatch.style.width = "14px";
    swatch.style.height = "14px";
    swatch.style.marginLeft = "8px";
    swatch.style.borderRadius = "3px";
    swatch.style.background = item.color;

    const label = document.createElement("label");
    label.textContent = " " + item.label;
    label.className = "ms-1";

    row.appendChild(checkbox);
    row.appendChild(swatch);
    row.appendChild(label);
    container.appendChild(row);
  });
};

Bin2CellExplorer.toggleGeneLayer = function(gene, visible) {
  const group = Bin2CellExplorer.state.layers.gene[gene];
  if (!group) return;
  group.style("display", visible ? null : "none");
};

Bin2CellExplorer.toggleCategoryLayer = function(category, visible) {
  const group = Bin2CellExplorer.state.layers.categories[category];
  if (!group) return;
  group.style("display", visible ? null : "none");
};
