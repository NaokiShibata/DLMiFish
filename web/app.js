(function () {
  const els = {
    projectFiles: document.getElementById("projectFiles"),
    configPath: document.getElementById("configPath"),
    taxa: document.getElementById("taxa"),
    markers: document.getElementById("markers"),
    remoteEntrezApi: document.getElementById("remoteEntrezApi"),
    workers: document.getElementById("workers"),
    outPath: document.getElementById("outPath"),
    outputPrefix: document.getElementById("outputPrefix"),
    dumpGbPath: document.getElementById("dumpGbPath"),
    fromGbPath: document.getElementById("fromGbPath"),
    resume: document.getElementById("resume"),
    dryRun: document.getElementById("dryRun"),
    postPrep: document.getElementById("postPrep"),
    stepPrimerTrim: document.getElementById("stepPrimerTrim"),
    stepLengthFilter: document.getElementById("stepLengthFilter"),
    stepDuplicateReport: document.getElementById("stepDuplicateReport"),
    postPrepPrimerSets: document.getElementById("postPrepPrimerSets"),
    runBtn: document.getElementById("runBtn"),
    status: document.getElementById("status"),
    logs: document.getElementById("logs"),
    downloads: document.getElementById("downloads"),
  };

  let worker = null;
  let runInProgress = false;

  function setStatus(message, kind = "") {
    els.status.textContent = message;
    els.status.classList.remove("error", "ok");
    if (kind) {
      els.status.classList.add(kind);
    }
  }

  function appendLog(message) {
    const line = String(message || "");
    els.logs.textContent += line + "\n";
    els.logs.scrollTop = els.logs.scrollHeight;
  }

  function splitValues(raw) {
    return String(raw || "")
      .split(/[\n,]/)
      .map((v) => v.trim())
      .filter(Boolean);
  }

  function clearDownloads() {
    els.downloads.innerHTML = "";
  }

  function renderDownloads(files) {
    clearDownloads();
    if (!files.length) {
      const li = document.createElement("li");
      li.textContent = "No output files";
      els.downloads.appendChild(li);
      return;
    }
    for (const file of files) {
      const blob = new Blob([file.buffer], { type: file.mime || "application/octet-stream" });
      const href = URL.createObjectURL(blob);
      const li = document.createElement("li");
      const a = document.createElement("a");
      a.href = href;
      a.download = file.name;
      a.textContent = `${file.name} (${file.path})`;
      li.appendChild(a);
      els.downloads.appendChild(li);
    }
  }

  async function collectProjectFiles() {
    const list = Array.from(els.projectFiles.files || []);
    const rawPaths = list.map((file) => (file.webkitRelativePath || file.name).replace(/\\/g, "/"));
    const roots = new Set(rawPaths.map((p) => p.split("/")[0]).filter(Boolean));
    const stripRoot = roots.size === 1;

    const payloads = [];
    for (let i = 0; i < list.length; i += 1) {
      const file = list[i];
      const buffer = await file.arrayBuffer();
      let relPath = rawPaths[i];
      if (stripRoot && relPath.includes("/")) {
        relPath = relPath.split("/").slice(1).join("/");
      }
      payloads.push({
        path: relPath,
        buffer,
      });
    }
    return payloads;
  }

  function createRunPayload() {
    const postPrepSteps = [];
    if (els.stepPrimerTrim.checked) postPrepSteps.push("primer_trim");
    if (els.stepLengthFilter.checked) postPrepSteps.push("length_filter");
    if (els.stepDuplicateReport.checked) postPrepSteps.push("duplicate_report");

    return {
      config_path: els.configPath.value.trim(),
      taxa: splitValues(els.taxa.value),
      markers: splitValues(els.markers.value),
      remote_entrez_api: els.remoteEntrezApi.value.trim(),
      workers: Number(els.workers.value || 2),
      out: els.outPath.value.trim(),
      output_prefix: els.outputPrefix.value.trim(),
      dump_gb: els.dumpGbPath.value.trim(),
      from_gb: els.fromGbPath.value.trim(),
      resume: els.resume.checked,
      dry_run: els.dryRun.checked,
      post_prep: els.postPrep.checked,
      post_prep_step: postPrepSteps,
      post_prep_primer_set: splitValues(els.postPrepPrimerSets.value),
    };
  }

  function ensureWorker() {
    if (worker) return;
    worker = new Worker("./worker.js");
    worker.onmessage = (event) => {
      const msg = event.data || {};
      if (msg.type === "status") {
        appendLog(`[status] ${msg.message}`);
        setStatus(msg.message);
        return;
      }
      if (msg.type === "log") {
        appendLog(msg.message);
        return;
      }
      if (msg.type === "ready") {
        setStatus("Worker ready", "ok");
        appendLog("[ready] Pyodide worker initialized.");
        return;
      }
      if (msg.type === "result") {
        runInProgress = false;
        els.runBtn.disabled = false;
        setStatus("Build completed", "ok");
        appendLog("[done] Build completed.");
        renderDownloads(msg.files || []);
        return;
      }
      if (msg.type === "error") {
        runInProgress = false;
        els.runBtn.disabled = false;
        setStatus("Build failed", "error");
        appendLog(`[error] ${msg.error}`);
      }
    };
    worker.onerror = (event) => {
      runInProgress = false;
      els.runBtn.disabled = false;
      setStatus("Worker error", "error");
      appendLog(`[worker-error] ${event.message}`);
    };
    worker.postMessage({ type: "init" });
  }

  async function runBuild() {
    if (runInProgress) return;
    clearDownloads();
    els.logs.textContent = "";
    ensureWorker();
    runInProgress = true;
    els.runBtn.disabled = true;
    setStatus("Preparing inputs...");
    appendLog("[start] Preparing input files...");

    try {
      const files = await collectProjectFiles();
      const payload = createRunPayload();

      if (!payload.config_path) {
        throw new Error("Config path is required.");
      }
      if (!payload.taxa.length) {
        throw new Error("At least one taxon is required.");
      }
      if (!payload.markers.length) {
        throw new Error("At least one marker is required.");
      }
      if (!payload.out) {
        throw new Error("Output FASTA path is required.");
      }
      if (!payload.remote_entrez_api) {
        throw new Error("Remote Entrez API is required for browser hybrid mode.");
      }

      const transferList = files.map((f) => f.buffer);
      worker.postMessage(
        {
          type: "run",
          payload,
          files,
        },
        transferList
      );
      setStatus("Running build...");
      appendLog("[start] Build requested.");
    } catch (err) {
      runInProgress = false;
      els.runBtn.disabled = false;
      setStatus("Input error", "error");
      appendLog(`[error] ${err && err.message ? err.message : String(err)}`);
    }
  }

  els.runBtn.addEventListener("click", () => {
    void runBuild();
  });

  ensureWorker();
})();
