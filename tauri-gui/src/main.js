import "./styles.css";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import {
  DEFAULTS,
  parseCommaSeparatedList,
  parseFloatOrNull,
  parseIntOrNull,
  readCheckedValues,
  setCheckedValues
} from "./lib/form-utils.js";
import { createMonitorView } from "./lib/monitor-view.js";
import { mergeUniqueValues, parseDelimitedTokens, renderTokenPills } from "./lib/token-list.js";

const els = {
  tabSetup: document.querySelector("#tab-setup"),
  tabMonitor: document.querySelector("#tab-monitor"),
  tabResults: document.querySelector("#tab-results"),
  viewSetup: document.querySelector("#view-setup"),
  viewMonitor: document.querySelector("#view-monitor"),
  viewResults: document.querySelector("#view-results"),
  form: document.querySelector("#setup-form"),
  reset: document.querySelector("#reset"),
  run: document.querySelector("#run"),
  cancel: document.querySelector("#cancel"),
  status: document.querySelector("#status"),
  phase: document.querySelector("#phase"),
  progress: document.querySelector("#progress"),
  progressDetail: document.querySelector("#progress-detail"),
  logs: document.querySelector("#log-view"),
  metricsList: document.querySelector("#metrics-list"),
  resultFiles: document.querySelector("#result-files"),
  openJobDir: document.querySelector("#open-job-dir"),
  newJob: document.querySelector("#new-job"),
  pickOutput: document.querySelector("#pick-output"),
  pickPrimer: document.querySelector("#pick-primer"),
  taxidInput: document.querySelector("#taxid-input"),
  taxonNameInput: document.querySelector("#taxon-name-input"),
  taxonCandidates: document.querySelector("#taxon-candidates"),
  addTaxid: document.querySelector("#add-taxid"),
  clearTaxids: document.querySelector("#clear-taxids"),
  taxidList: document.querySelector("#taxid-list"),
  taxidCount: document.querySelector("#taxid-count"),
  taxidsHidden: document.querySelector("#taxids"),
  markerSelect: document.querySelector("#marker"),
  addMarker: document.querySelector("#add-marker"),
  clearMarkers: document.querySelector("#clear-markers"),
  markerList: document.querySelector("#marker-list"),
  markerCount: document.querySelector("#marker-count"),
  flowTaxids: document.querySelector("#flow-taxids"),
  flowMarkers: document.querySelector("#flow-markers"),
  flowOutput: document.querySelector("#flow-output"),
  flowEmail: document.querySelector("#flow-email"),
  flowPostEnable: document.querySelector("#flow-post-enable"),
  flowPostPrimer: document.querySelector("#flow-post-primer"),
  flowPostLength: document.querySelector("#flow-post-length"),
  outputPrefixInput: document.querySelector("#output-prefix"),
  outputRootInput: document.querySelector("#output-root"),
  emailInput: document.querySelector("#email"),
  apiKeyInput: document.querySelector("#api-key"),
  saveApiKeyInput: document.querySelector("#save-api-key"),
  filterMitoInput: document.querySelector("#f-mito"),
  filterDdbjInput: document.querySelector("#f-ddbj"),
  filterBiomolInput: document.querySelector("#f-biomol"),
  filterLengthMinInput: document.querySelector("#f-length-min"),
  filterLengthMaxInput: document.querySelector("#f-length-max"),
  postEnableInput: document.querySelector("#post-enable"),
  primerFileInput: document.querySelector("#primer-file"),
  primerSetInput: document.querySelector("#primer-set"),
  postLengthMinInput: document.querySelector("#post-length-min"),
  postLengthMaxInput: document.querySelector("#post-length-max"),
  sourceInput: document.querySelector("#source"),
  ncbiDbInput: document.querySelector("#ncbi-db"),
  ncbiRettypeInput: document.querySelector("#ncbi-rettype"),
  ncbiRetmodeInput: document.querySelector("#ncbi-retmode"),
  ncbiPerQueryInput: document.querySelector("#ncbi-per-query"),
  ncbiUseHistoryInput: document.querySelector("#ncbi-use-history"),
  ncbiDelaySecInput: document.querySelector("#ncbi-delay-sec"),
  outputDefaultHeaderFormatInput: document.querySelector("#output-default-header-format"),
  outputMifishHeaderFormatInput: document.querySelector("#output-mifish-header-format"),
  speedInput: document.querySelector("#speed"),
  resumeInput: document.querySelector("#resume"),
  duplicateReportStep: document.querySelector('.post-step[value="duplicate_report"]'),
  loadDbTomlBtn: document.querySelector("#load-db-toml"),
  loadedDbTomlPath: document.querySelector("#loaded-db-toml")
};

const state = {
  jobDir: "",
  files: [],
  unlisten: null,
  taxids: [],
  markers: [],
  taxonCandidateItems: [],
  taxonSearchSeq: 0,
  taxonSearchTimer: null
};

const postStepEls = Array.from(document.querySelectorAll(".post-step"));
const monitorView = createMonitorView({
  statusEl: els.status,
  phaseEl: els.phase,
  progressEl: els.progress,
  progressDetailEl: els.progressDetail,
  logsEl: els.logs,
  metricsListEl: els.metricsList,
  resultFilesEl: els.resultFiles,
  openPath: (path) => invoke("open_path", { path })
});

function switchView(name) {
  const map = {
    setup: [els.viewSetup, els.tabSetup],
    monitor: [els.viewMonitor, els.tabMonitor],
    results: [els.viewResults, els.tabResults]
  };
  for (const key of Object.keys(map)) {
    const [view, tab] = map[key];
    const active = key === name;
    view.classList.toggle("active", active);
    tab.classList.toggle("active", active);
  }
}

function setMonitorEnabled(enabled) {
  els.tabMonitor.disabled = !enabled;
}

function setResultsEnabled(enabled) {
  els.tabResults.disabled = !enabled;
}

function parseTaxidTokens(raw) {
  return parseDelimitedTokens(raw);
}

function currentBuildSource() {
  const raw = `${els.sourceInput?.value || ""}`.trim().toLowerCase();
  if (raw === "bold" || raw === "both") return raw;
  return "ncbi";
}

function sourceUsesNcbi() {
  return currentBuildSource() !== "bold";
}

function syncSourceMode({ resetDuplicateForBoth = false } = {}) {
  const source = currentBuildSource();
  const usesNcbi = sourceUsesNcbi();

  if (els.emailInput) {
    els.emailInput.required = usesNcbi;
    els.emailInput.placeholder = usesNcbi
      ? "your.email@example.com"
      : "optional for BOLD-only";
  }

  if (els.flowEmail) {
    els.flowEmail.textContent = usesNcbi
      ? "4. Set NCBI email"
      : "4. NCBI email is optional for BOLD-only";
  }

  if (els.resumeInput) {
    if (!usesNcbi) {
      els.resumeInput.checked = false;
    }
    els.resumeInput.disabled = !usesNcbi;
  }

  if (resetDuplicateForBoth && source === "both" && els.duplicateReportStep) {
    els.duplicateReportStep.checked = false;
  }
}

function syncTaxidsHidden() {
  els.taxidsHidden.value = state.taxids.join("\n");
}

function renderTaxids() {
  renderTokenPills({
    container: els.taxidList,
    countLabel: els.taxidCount,
    items: state.taxids,
    countSuffix: "taxids",
    onRemove: (taxid) => {
      state.taxids = state.taxids.filter((value) => value !== taxid);
      renderTaxids();
      syncTaxidsHidden();
      updateGuidanceState();
    }
  });
}

function addTaxids(raw) {
  const tokens = parseTaxidTokens(raw);
  if (!tokens.length) return;
  state.taxids = mergeUniqueValues(state.taxids, tokens);
  renderTaxids();
  syncTaxidsHidden();
}

function clearTaxonCandidates() {
  state.taxonCandidateItems = [];
  if (els.taxonCandidates) {
    els.taxonCandidates.innerHTML = "";
    els.taxonCandidates.classList.remove("active");
  }
}

function pickTaxonCandidate(item) {
  if (!item?.taxId) return;
  addTaxids(item.taxId);
  updateGuidanceState();
  els.taxonNameInput.value = `${item.scientificName} : ${item.taxId}`;
  clearTaxonCandidates();
  els.taxidInput.focus();
}

function renderTaxonCandidates(items) {
  if (!els.taxonCandidates) return;
  els.taxonCandidates.innerHTML = "";
  state.taxonCandidateItems = Array.isArray(items) ? items : [];

  if (!state.taxonCandidateItems.length) {
    els.taxonCandidates.classList.remove("active");
    return;
  }

  for (const item of state.taxonCandidateItems) {
    const li = document.createElement("li");
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "taxon-candidate-item";
    btn.textContent = `${item.scientificName} : ${item.taxId}`;
    btn.title = `${item.scientificName} : ${item.taxId}`;
    btn.addEventListener("click", () => pickTaxonCandidate(item));
    li.appendChild(btn);
    els.taxonCandidates.appendChild(li);
  }
  els.taxonCandidates.classList.add("active");
}

async function searchTaxonCandidatesNow() {
  const query = els.taxonNameInput?.value?.trim() || "";

  if (query.length < 2) {
    clearTaxonCandidates();
    return;
  }

  const seq = ++state.taxonSearchSeq;
  try {
    const items = await invoke("search_taxonomy", {
      query,
      limit: 12
    });
    if (seq !== state.taxonSearchSeq) return;
    renderTaxonCandidates(items || []);
  } catch (error) {
    if (seq !== state.taxonSearchSeq) return;
    clearTaxonCandidates();
    monitorView.appendLog(`[warn] taxon search failed: ${error}`);
  }
}

function scheduleTaxonSearch() {
  if (state.taxonSearchTimer) {
    clearTimeout(state.taxonSearchTimer);
  }
  state.taxonSearchTimer = setTimeout(() => {
    state.taxonSearchTimer = null;
    searchTaxonCandidatesNow();
  }, 250);
}

function renderMarkers() {
  renderTokenPills({
    container: els.markerList,
    countLabel: els.markerCount,
    items: state.markers,
    countSuffix: "markers",
    pillClass: "taxid-pill marker-pill",
    onRemove: (marker) => {
      state.markers = state.markers.filter((value) => value !== marker);
      renderMarkers();
      updateGuidanceState();
    }
  });
}

function addMarkers(values) {
  state.markers = mergeUniqueValues(state.markers, values);
  renderMarkers();
}

function setFlowItemState(el, done) {
  if (!el) return;
  el.classList.toggle("done", done);
  el.classList.remove("neutral");
  el.classList.toggle("pending", !done);
}

function setFlowItemNeutral(el) {
  if (!el) return;
  el.classList.remove("done");
  el.classList.remove("pending");
  el.classList.add("neutral");
}

function setPostPrepStepSelection(steps) {
  setCheckedValues(postStepEls, steps);
}

function applyImportedDbToml(imported) {
  els.sourceInput.value = imported.source || "ncbi";
  els.emailInput.value = imported.email || "";

  const apiKey = imported.apiKey || "";
  els.apiKeyInput.value = apiKey;
  els.saveApiKeyInput.checked = apiKey.length > 0;

  const filters = imported.filters || {};
  els.filterMitoInput.checked = Boolean(filters.mitochondrion);
  els.filterDdbjInput.checked = Boolean(filters.ddbjEmblGenbank);
  els.filterBiomolInput.checked = Boolean(filters.biomolGenomic);
  els.filterLengthMinInput.value = filters.lengthMin ?? "";
  els.filterLengthMaxInput.value = filters.lengthMax ?? "";

  const postPrep = imported.postPrep || {};
  els.postEnableInput.checked = Boolean(postPrep.enable);
  els.primerFileInput.value = postPrep.primerFile || "";
  els.primerSetInput.value = (postPrep.primerSet || []).join(",");
  els.postLengthMinInput.value = postPrep.sequenceLengthMin ?? "";
  els.postLengthMaxInput.value = postPrep.sequenceLengthMax ?? "";
  setPostPrepStepSelection(postPrep.steps || []);

  const ncbiOptions = imported.ncbiOptions || {};
  els.ncbiDbInput.value = ncbiOptions.db || DEFAULTS.ncbiDb;
  els.ncbiRettypeInput.value = ncbiOptions.rettype || DEFAULTS.ncbiRettype;
  els.ncbiRetmodeInput.value = ncbiOptions.retmode || DEFAULTS.ncbiRetmode;
  els.ncbiPerQueryInput.value = ncbiOptions.perQuery || DEFAULTS.ncbiPerQuery;
  els.ncbiUseHistoryInput.checked = ncbiOptions.useHistory !== false;
  els.ncbiDelaySecInput.value = ncbiOptions.delaySec ?? "";

  const outputOptions = imported.outputOptions || {};
  els.outputDefaultHeaderFormatInput.value =
    outputOptions.defaultHeaderFormat || DEFAULTS.defaultHeaderFormat;
  els.outputMifishHeaderFormatInput.value =
    outputOptions.mifishHeaderFormat || DEFAULTS.mifishHeaderFormat;
  syncSourceMode();
}

function isReadyToRun() {
  return (
    state.taxids.length > 0 &&
    state.markers.length > 0 &&
    els.outputRootInput.value.trim().length > 0 &&
    (!sourceUsesNcbi() || els.emailInput.value.trim().length > 0)
  );
}

function updatePostPrepGuidance() {
  const postEnabled = els.postEnableInput.checked;
  const selectedSteps = new Set(getSelectedPostPrepSteps());

  if (!postEnabled) {
    setFlowItemNeutral(els.flowPostEnable);
    setFlowItemNeutral(els.flowPostPrimer);
    setFlowItemNeutral(els.flowPostLength);
    return;
  }

  setFlowItemState(els.flowPostEnable, true);

  if (selectedSteps.has("primer_trim")) {
    const hasPrimerFile = els.primerFileInput.value.trim().length > 0;
    const hasPrimerSet = parseCommaSeparatedList(els.primerSetInput.value).length > 0;
    setFlowItemState(els.flowPostPrimer, hasPrimerFile && hasPrimerSet);
  } else {
    setFlowItemNeutral(els.flowPostPrimer);
  }

  if (selectedSteps.has("length_filter")) {
    const hasLengthBounds = parseIntOrNull(els.postLengthMinInput.value) != null || parseIntOrNull(els.postLengthMaxInput.value) != null;
    setFlowItemState(els.flowPostLength, hasLengthBounds);
  } else {
    setFlowItemNeutral(els.flowPostLength);
  }
}

function updateGuidanceState() {
  setFlowItemState(els.flowTaxids, state.taxids.length > 0);
  setFlowItemState(els.flowMarkers, state.markers.length > 0);
  setFlowItemState(els.flowOutput, els.outputRootInput.value.trim().length > 0);
  setFlowItemState(els.flowEmail, els.emailInput.value.trim().length > 0);
  updatePostPrepGuidance();
  els.run.disabled = !isReadyToRun();
}

function getSelectedPostPrepSteps() {
  return readCheckedValues(postStepEls);
}

function collectRequest() {
  const steps = getSelectedPostPrepSteps();
  const primerSet = parseCommaSeparatedList(els.primerSetInput.value);

  return {
    taxids: [...state.taxids],
    markers: [...state.markers],
    source: currentBuildSource(),
    outputPrefix: els.outputPrefixInput.value.trim() || DEFAULTS.outputPrefix,
    outputRoot: els.outputRootInput.value.trim(),
    email: els.emailInput.value.trim(),
    apiKey: els.apiKeyInput.value.trim(),
    saveApiKey: els.saveApiKeyInput.checked,
    filters: {
      mitochondrion: els.filterMitoInput.checked,
      ddbjEmblGenbank: els.filterDdbjInput.checked,
      biomolGenomic: els.filterBiomolInput.checked,
      lengthMin: parseIntOrNull(els.filterLengthMinInput.value),
      lengthMax: parseIntOrNull(els.filterLengthMaxInput.value)
    },
    postPrep: {
      enable: els.postEnableInput.checked,
      primerFile: els.primerFileInput.value.trim(),
      primerSet,
      steps,
      sequenceLengthMin: parseIntOrNull(els.postLengthMinInput.value),
      sequenceLengthMax: parseIntOrNull(els.postLengthMaxInput.value)
    },
    ncbiOptions: {
      db: els.ncbiDbInput.value.trim() || DEFAULTS.ncbiDb,
      rettype: els.ncbiRettypeInput.value.trim() || DEFAULTS.ncbiRettype,
      retmode: els.ncbiRetmodeInput.value.trim() || DEFAULTS.ncbiRetmode,
      perQuery: Math.max(1, Number.parseInt(els.ncbiPerQueryInput.value, 10) || DEFAULTS.ncbiPerQuery),
      useHistory: els.ncbiUseHistoryInput.checked,
      delaySec: parseFloatOrNull(els.ncbiDelaySecInput.value)
    },
    outputOptions: {
      defaultHeaderFormat: els.outputDefaultHeaderFormatInput.value.trim() || DEFAULTS.defaultHeaderFormat,
      mifishHeaderFormat: els.outputMifishHeaderFormatInput.value.trim() || DEFAULTS.mifishHeaderFormat
    },
    workers: Number.parseInt(els.speedInput.value, 10),
    resume: els.resumeInput.checked
  };
}

async function setupEventListener() {
  if (state.unlisten) {
    state.unlisten();
    state.unlisten = null;
  }

  state.unlisten = await listen("run-event", (event) => {
    const payload = event.payload;
    if (!payload || typeof payload !== "object") return;

    if (payload.eventType === "log" && payload.line) {
      monitorView.appendLog(payload.line);
    }
    if (payload.eventType === "status" && payload.status) {
      els.status.textContent = payload.status;
    }
    if (payload.eventType === "progress") {
      if (payload.phase) els.phase.textContent = payload.phase;
      if (payload.percent != null) els.progress.value = payload.percent;
      monitorView.renderProgressDetail(payload.phase, payload.percent, payload.metrics);
      if (payload.metrics) monitorView.renderMetrics(payload.metrics);
    }
    if (payload.eventType === "result") {
      state.jobDir = payload.jobDir || "";
      state.files = payload.files || [];
      setResultsEnabled(true);
      monitorView.renderResultFiles(state.files);
      switchView("results");
    }
    if (payload.eventType === "error" && payload.message) {
      monitorView.appendLog(`[error] ${payload.message}`);
      els.status.textContent = "Failed";
    }
  });
}

async function loadSavedConfig() {
  try {
    const saved = await invoke("load_gui_config");
    if (!saved) return;
    els.outputRootInput.value = saved.outputRoot || "";
    els.outputPrefixInput.value = saved.outputPrefix || DEFAULTS.outputPrefix;
    els.emailInput.value = saved.email || "";
    els.apiKeyInput.value = saved.saveApiKey ? saved.apiKey || "" : "";
    els.saveApiKeyInput.checked = Boolean(saved.saveApiKey);
    if (saved.marker) {
      els.markerSelect.value = saved.marker;
      addMarkers(saved.marker);
    } else if (!state.markers.length) {
      addMarkers(els.markerSelect.value);
    }
    if (saved.source) {
      els.sourceInput.value = saved.source;
      syncSourceMode();
    }
    if (saved.workers) els.speedInput.value = `${saved.workers}`;
    if (saved.ncbiDb) els.ncbiDbInput.value = saved.ncbiDb;
    if (saved.ncbiRettype) els.ncbiRettypeInput.value = saved.ncbiRettype;
    if (saved.ncbiRetmode) els.ncbiRetmodeInput.value = saved.ncbiRetmode;
    if (saved.ncbiPerQuery) els.ncbiPerQueryInput.value = `${saved.ncbiPerQuery}`;
    if (saved.ncbiUseHistory != null) els.ncbiUseHistoryInput.checked = Boolean(saved.ncbiUseHistory);
    if (saved.ncbiDelaySec != null) els.ncbiDelaySecInput.value = `${saved.ncbiDelaySec}`;
    if (saved.outputDefaultHeaderFormat) {
      els.outputDefaultHeaderFormatInput.value = saved.outputDefaultHeaderFormat;
    }
    if (saved.outputMifishHeaderFormat) {
      els.outputMifishHeaderFormatInput.value = saved.outputMifishHeaderFormat;
    }
    updateGuidanceState();
  } catch (error) {
    monitorView.appendLog(`[warn] load config failed: ${error}`);
  }
}

async function runJob() {
  const req = collectRequest();
  if (!req.taxids.length) throw new Error("please input at least one taxid");
  if (!req.markers.length) throw new Error("please add at least one marker");
  if (!req.outputRoot) throw new Error("output directory is required");
  if (sourceUsesNcbi() && !req.email) throw new Error("email is required for ncbi/both");

  await setupEventListener();
  monitorView.reset();
  setMonitorEnabled(true);
  setResultsEnabled(false);
  switchView("monitor");

  const resp = await invoke("start_run", { req });
  monitorView.appendLog(`job dir: ${resp.jobDir}`);
  monitorView.appendLog(`log path: ${resp.logPath}`);
}

function resetForm() {
  state.taxids = [];
  renderTaxids();
  syncTaxidsHidden();
  els.taxidInput.value = "";
  els.taxonNameInput.value = "";
  clearTaxonCandidates();
  state.markers = [];
  addMarkers(els.markerSelect.value);
  els.filterLengthMinInput.value = "";
  els.filterLengthMaxInput.value = "";
  els.postLengthMinInput.value = "";
  els.postLengthMaxInput.value = "";
  els.primerFileInput.value = "";
  els.primerSetInput.value = "";
  els.postEnableInput.checked = false;
  setPostPrepStepSelection([]);
  syncSourceMode();
  updateGuidanceState();
}

els.form.addEventListener("submit", async (event) => {
  event.preventDefault();
  try {
    await runJob();
  } catch (error) {
    monitorView.appendLog(`[error] ${error}`);
  }
});

els.reset.addEventListener("click", () => {
  resetForm();
});

els.cancel.addEventListener("click", async () => {
  try {
    await invoke("cancel_run");
  } catch (error) {
    monitorView.appendLog(`[error] cancel failed: ${error}`);
  }
});

els.pickOutput.addEventListener("click", async () => {
  const selected = await invoke("choose_output_directory");
  if (selected) {
    els.outputRootInput.value = selected;
    updateGuidanceState();
  }
});

els.pickPrimer.addEventListener("click", async () => {
  const selected = await invoke("choose_primer_file");
  if (selected) {
    els.primerFileInput.value = selected;
    updateGuidanceState();
  }
});

els.loadDbTomlBtn.addEventListener("click", async () => {
  try {
    const selectedPath = await invoke("choose_db_toml_file");
    if (!selectedPath) return;
    const imported = await invoke("import_db_toml", { path: selectedPath });
    applyImportedDbToml(imported);
    els.loadedDbTomlPath.textContent = imported.sourcePath || selectedPath;
    updateGuidanceState();
  } catch (error) {
    monitorView.appendLog(`[error] db.toml import failed: ${error}`);
  }
});

els.addTaxid.addEventListener("click", () => {
  addTaxids(els.taxidInput.value);
  els.taxidInput.value = "";
  els.taxidInput.focus();
  updateGuidanceState();
});

els.taxidInput.addEventListener("keydown", (event) => {
  if (event.key === "Enter") {
    event.preventDefault();
    addTaxids(els.taxidInput.value);
    els.taxidInput.value = "";
    updateGuidanceState();
  }
});

els.taxonNameInput?.addEventListener("input", () => {
  scheduleTaxonSearch();
});

els.taxonNameInput?.addEventListener("keydown", (event) => {
  if (event.key === "Enter") {
    event.preventDefault();
    if (state.taxonCandidateItems.length > 0) {
      pickTaxonCandidate(state.taxonCandidateItems[0]);
    }
  } else if (event.key === "Escape") {
    clearTaxonCandidates();
  }
});

els.clearTaxids.addEventListener("click", () => {
  state.taxids = [];
  renderTaxids();
  syncTaxidsHidden();
  els.taxidInput.focus();
  updateGuidanceState();
});

els.addMarker.addEventListener("click", () => {
  addMarkers(els.markerSelect.value);
  updateGuidanceState();
});

els.clearMarkers.addEventListener("click", () => {
  state.markers = [];
  renderMarkers();
  updateGuidanceState();
});

els.outputRootInput.addEventListener("input", updateGuidanceState);
els.emailInput.addEventListener("input", updateGuidanceState);
els.sourceInput.addEventListener("change", () => {
  syncSourceMode({ resetDuplicateForBoth: true });
  updateGuidanceState();
});
els.postEnableInput.addEventListener("change", updateGuidanceState);
els.primerFileInput.addEventListener("input", updateGuidanceState);
els.primerSetInput.addEventListener("input", updateGuidanceState);
els.postLengthMinInput.addEventListener("input", updateGuidanceState);
els.postLengthMaxInput.addEventListener("input", updateGuidanceState);
postStepEls.forEach((el) => {
  el.addEventListener("change", updateGuidanceState);
});

els.openJobDir.addEventListener("click", async () => {
  if (state.jobDir) {
    await invoke("open_path", { path: state.jobDir });
  }
});

els.newJob.addEventListener("click", () => {
  switchView("setup");
});

els.tabSetup.addEventListener("click", () => switchView("setup"));
els.tabMonitor.addEventListener("click", () => switchView("monitor"));
els.tabResults.addEventListener("click", () => switchView("results"));

addTaxids(els.taxidsHidden.value);
renderTaxids();
syncTaxidsHidden();
syncSourceMode();
loadSavedConfig();
updateGuidanceState();
