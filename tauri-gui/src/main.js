import "./styles.css";
import { invoke } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";

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
  logs: document.querySelector("#log-view"),
  metricsList: document.querySelector("#metrics-list"),
  resultFiles: document.querySelector("#result-files"),
  openJobDir: document.querySelector("#open-job-dir"),
  newJob: document.querySelector("#new-job"),
  pickOutput: document.querySelector("#pick-output"),
  pickPrimer: document.querySelector("#pick-primer"),
  taxidInput: document.querySelector("#taxid-input"),
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
  outputRootInput: document.querySelector("#output-root"),
  emailInput: document.querySelector("#email"),
  postEnableInput: document.querySelector("#post-enable"),
  primerFileInput: document.querySelector("#primer-file"),
  primerSetInput: document.querySelector("#primer-set"),
  postLengthMinInput: document.querySelector("#post-length-min"),
  postLengthMaxInput: document.querySelector("#post-length-max"),
  loadDbTomlBtn: document.querySelector("#load-db-toml"),
  loadedDbTomlPath: document.querySelector("#loaded-db-toml")
};

const state = {
  jobDir: "",
  files: [],
  unlisten: null,
  taxids: [],
  markers: []
};

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

function parseIntOrNull(value) {
  const v = `${value}`.trim();
  if (!v) return null;
  const parsed = Number.parseInt(v, 10);
  return Number.isNaN(parsed) ? null : parsed;
}

function parseTaxidTokens(raw) {
  return `${raw}`
    .split(/[\s,;]+/)
    .map((v) => v.trim())
    .filter(Boolean);
}

function syncTaxidsHidden() {
  els.taxidsHidden.value = state.taxids.join("\n");
}

function renderTaxids() {
  els.taxidList.innerHTML = "";
  els.taxidList.classList.toggle("empty", state.taxids.length === 0);
  els.taxidCount.textContent = `${state.taxids.length} taxids`;

  for (const taxid of state.taxids) {
    const li = document.createElement("li");
    li.className = "taxid-pill";

    const text = document.createElement("span");
    text.textContent = taxid;

    const remove = document.createElement("button");
    remove.type = "button";
    remove.className = "taxid-remove";
    remove.textContent = "×";
    remove.title = `remove ${taxid}`;
    remove.addEventListener("click", () => {
      state.taxids = state.taxids.filter((v) => v !== taxid);
      renderTaxids();
      syncTaxidsHidden();
      updateGuidanceState();
    });

    li.append(text, remove);
    els.taxidList.appendChild(li);
  }
}

function addTaxids(raw) {
  const tokens = parseTaxidTokens(raw);
  if (!tokens.length) return;

  const merged = new Set(state.taxids);
  for (const taxid of tokens) {
    merged.add(taxid);
  }
  state.taxids = Array.from(merged);
  renderTaxids();
  syncTaxidsHidden();
}

function renderMarkers() {
  els.markerList.innerHTML = "";
  els.markerList.classList.toggle("empty", state.markers.length === 0);
  els.markerCount.textContent = `${state.markers.length} markers`;

  for (const marker of state.markers) {
    const li = document.createElement("li");
    li.className = "taxid-pill marker-pill";

    const text = document.createElement("span");
    text.textContent = marker;

    const remove = document.createElement("button");
    remove.type = "button";
    remove.className = "taxid-remove";
    remove.textContent = "×";
    remove.title = `remove ${marker}`;
    remove.addEventListener("click", () => {
      state.markers = state.markers.filter((v) => v !== marker);
      renderMarkers();
      updateGuidanceState();
    });

    li.append(text, remove);
    els.markerList.appendChild(li);
  }
}

function addMarkers(values) {
  const arr = Array.isArray(values) ? values : [values];
  const merged = new Set(state.markers);
  for (const raw of arr) {
    const marker = `${raw}`.trim();
    if (marker) merged.add(marker);
  }
  state.markers = Array.from(merged);
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
  const wanted = new Set(Array.isArray(steps) ? steps : []);
  Array.from(document.querySelectorAll(".post-step")).forEach((el) => {
    el.checked = wanted.has(el.value);
  });
}

function applyImportedDbToml(imported) {
  document.querySelector("#email").value = imported.email || "";

  const apiKey = imported.apiKey || "";
  document.querySelector("#api-key").value = apiKey;
  document.querySelector("#save-api-key").checked = apiKey.length > 0;

  const filters = imported.filters || {};
  document.querySelector("#f-mito").checked = Boolean(filters.mitochondrion);
  document.querySelector("#f-ddbj").checked = Boolean(filters.ddbjEmblGenbank);
  document.querySelector("#f-biomol").checked = Boolean(filters.biomolGenomic);
  document.querySelector("#f-length-min").value = filters.lengthMin ?? "";
  document.querySelector("#f-length-max").value = filters.lengthMax ?? "";

  const postPrep = imported.postPrep || {};
  document.querySelector("#post-enable").checked = Boolean(postPrep.enable);
  document.querySelector("#primer-file").value = postPrep.primerFile || "";
  document.querySelector("#primer-set").value = (postPrep.primerSet || []).join(",");
  document.querySelector("#post-length-min").value = postPrep.sequenceLengthMin ?? "";
  document.querySelector("#post-length-max").value = postPrep.sequenceLengthMax ?? "";
  setPostPrepStepSelection(postPrep.steps || []);
}

function isReadyToRun() {
  return (
    state.taxids.length > 0 &&
    state.markers.length > 0 &&
    els.outputRootInput.value.trim().length > 0 &&
    els.emailInput.value.trim().length > 0
  );
}

function updatePostPrepGuidance() {
  const postEnabled = els.postEnableInput.checked;
  const selectedSteps = new Set(Array.from(document.querySelectorAll(".post-step:checked")).map((v) => v.value));

  if (!postEnabled) {
    setFlowItemNeutral(els.flowPostEnable);
    setFlowItemNeutral(els.flowPostPrimer);
    setFlowItemNeutral(els.flowPostLength);
    return;
  }

  setFlowItemState(els.flowPostEnable, true);

  if (selectedSteps.has("primer_trim")) {
    const hasPrimerFile = els.primerFileInput.value.trim().length > 0;
    const hasPrimerSet = els.primerSetInput.value
      .split(",")
      .map((v) => v.trim())
      .filter(Boolean).length > 0;
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

function collectRequest() {
  const steps = Array.from(document.querySelectorAll(".post-step:checked")).map((v) => v.value);
  const primerSet = document
    .querySelector("#primer-set")
    .value.split(",")
    .map((v) => v.trim())
    .filter(Boolean);

  return {
    taxids: [...state.taxids],
    markers: [...state.markers],
    outputPrefix: document.querySelector("#output-prefix").value.trim() || "MiFish",
    outputRoot: document.querySelector("#output-root").value.trim(),
    email: document.querySelector("#email").value.trim(),
    apiKey: document.querySelector("#api-key").value.trim(),
    saveApiKey: document.querySelector("#save-api-key").checked,
    filters: {
      mitochondrion: document.querySelector("#f-mito").checked,
      ddbjEmblGenbank: document.querySelector("#f-ddbj").checked,
      biomolGenomic: document.querySelector("#f-biomol").checked,
      lengthMin: parseIntOrNull(document.querySelector("#f-length-min").value),
      lengthMax: parseIntOrNull(document.querySelector("#f-length-max").value)
    },
    postPrep: {
      enable: document.querySelector("#post-enable").checked,
      primerFile: document.querySelector("#primer-file").value.trim(),
      primerSet,
      steps,
      sequenceLengthMin: parseIntOrNull(document.querySelector("#post-length-min").value),
      sequenceLengthMax: parseIntOrNull(document.querySelector("#post-length-max").value)
    },
    workers: Number.parseInt(document.querySelector("#speed").value, 10),
    resume: document.querySelector("#resume").checked
  };
}

function resetMonitor() {
  els.status.textContent = "Idle";
  els.phase.textContent = "-";
  els.progress.value = 0;
  els.logs.textContent = "";
  els.metricsList.innerHTML = "";
}

function appendLog(line) {
  const current = els.logs.textContent;
  const next = current.length > 50000 ? current.slice(-30000) : current;
  els.logs.textContent = `${next}${line}\n`;
  els.logs.scrollTop = els.logs.scrollHeight;
}

function renderMetrics(metrics) {
  const items = [];
  if (metrics?.queryCountByTaxid) {
    for (const [taxid, count] of Object.entries(metrics.queryCountByTaxid)) {
      items.push(`query count taxid=${taxid}: ${count}`);
    }
  }
  if (metrics?.matchedRecords != null) items.push(`matched records: ${metrics.matchedRecords}`);
  if (metrics?.keptRecordsBeforePostPrep != null) {
    items.push(`kept before post_prep: ${metrics.keptRecordsBeforePostPrep}`);
  }
  if (metrics?.primerTrimRemoved != null) items.push(`primer_trim removed: ${metrics.primerTrimRemoved}`);
  if (metrics?.lengthFilterRemoved != null) items.push(`length_filter removed: ${metrics.lengthFilterRemoved}`);
  if (metrics?.duplicateGroups != null) items.push(`duplicate groups: ${metrics.duplicateGroups}`);
  if (metrics?.crossOrganismGroups != null) {
    items.push(`cross_organism_groups: ${metrics.crossOrganismGroups}`);
  }

  els.metricsList.innerHTML = "";
  for (const text of items) {
    const li = document.createElement("li");
    li.textContent = text;
    els.metricsList.appendChild(li);
  }
}

function renderResultFiles(files) {
  els.resultFiles.innerHTML = "";
  for (const path of files) {
    const li = document.createElement("li");
    const row = document.createElement("div");
    row.className = "file-row";
    const text = document.createElement("code");
    text.textContent = path;
    const button = document.createElement("button");
    button.textContent = "Open";
    button.addEventListener("click", async () => {
      await invoke("open_path", { path });
    });
    row.append(text, button);
    li.appendChild(row);
    els.resultFiles.appendChild(li);
  }
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
      appendLog(payload.line);
    }
    if (payload.eventType === "status" && payload.status) {
      els.status.textContent = payload.status;
    }
    if (payload.eventType === "progress") {
      if (payload.phase) els.phase.textContent = payload.phase;
      if (payload.percent != null) els.progress.value = payload.percent;
      if (payload.metrics) renderMetrics(payload.metrics);
    }
    if (payload.eventType === "result") {
      state.jobDir = payload.jobDir || "";
      state.files = payload.files || [];
      setResultsEnabled(true);
      renderResultFiles(state.files);
      switchView("results");
    }
    if (payload.eventType === "error" && payload.message) {
      appendLog(`[error] ${payload.message}`);
      els.status.textContent = "Failed";
    }
  });
}

async function loadSavedConfig() {
  try {
    const saved = await invoke("load_gui_config");
    if (!saved) return;
    document.querySelector("#output-root").value = saved.outputRoot || "";
    document.querySelector("#output-prefix").value = saved.outputPrefix || "MiFish";
    document.querySelector("#email").value = saved.email || "";
    document.querySelector("#api-key").value = saved.saveApiKey ? saved.apiKey || "" : "";
    document.querySelector("#save-api-key").checked = Boolean(saved.saveApiKey);
    if (saved.marker) {
      document.querySelector("#marker").value = saved.marker;
      addMarkers(saved.marker);
    } else if (!state.markers.length) {
      addMarkers(document.querySelector("#marker").value);
    }
    if (saved.workers) document.querySelector("#speed").value = `${saved.workers}`;
    updateGuidanceState();
  } catch (error) {
    appendLog(`[warn] load config failed: ${error}`);
  }
}

async function runJob() {
  const req = collectRequest();
  if (!req.taxids.length) throw new Error("please input at least one taxid");
  if (!req.markers.length) throw new Error("please add at least one marker");
  if (!req.outputRoot) throw new Error("output directory is required");
  if (!req.email) throw new Error("email is required");

  await setupEventListener();
  resetMonitor();
  setMonitorEnabled(true);
  setResultsEnabled(false);
  switchView("monitor");

  const resp = await invoke("start_run", { req });
  appendLog(`job dir: ${resp.jobDir}`);
  appendLog(`log path: ${resp.logPath}`);
}

function resetForm() {
  state.taxids = [];
  renderTaxids();
  syncTaxidsHidden();
  els.taxidInput.value = "";
  state.markers = [];
  addMarkers(els.markerSelect.value);
  document.querySelector("#f-length-min").value = "";
  document.querySelector("#f-length-max").value = "";
  document.querySelector("#post-length-min").value = "";
  document.querySelector("#post-length-max").value = "";
  document.querySelector("#primer-file").value = "";
  document.querySelector("#primer-set").value = "";
  document.querySelector("#post-enable").checked = false;
  updateGuidanceState();
}

els.form.addEventListener("submit", async (event) => {
  event.preventDefault();
  try {
    await runJob();
  } catch (error) {
    appendLog(`[error] ${error}`);
  }
});

els.reset.addEventListener("click", () => {
  resetForm();
});

els.cancel.addEventListener("click", async () => {
  try {
    await invoke("cancel_run");
  } catch (error) {
    appendLog(`[error] cancel failed: ${error}`);
  }
});

els.pickOutput.addEventListener("click", async () => {
  const selected = await invoke("choose_output_directory");
  if (selected) {
    document.querySelector("#output-root").value = selected;
    updateGuidanceState();
  }
});

els.pickPrimer.addEventListener("click", async () => {
  const selected = await invoke("choose_primer_file");
  if (selected) {
    document.querySelector("#primer-file").value = selected;
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
    appendLog(`[error] db.toml import failed: ${error}`);
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
els.postEnableInput.addEventListener("change", updateGuidanceState);
els.primerFileInput.addEventListener("input", updateGuidanceState);
els.primerSetInput.addEventListener("input", updateGuidanceState);
els.postLengthMinInput.addEventListener("input", updateGuidanceState);
els.postLengthMaxInput.addEventListener("input", updateGuidanceState);
Array.from(document.querySelectorAll(".post-step")).forEach((el) => {
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
loadSavedConfig();
updateGuidanceState();
