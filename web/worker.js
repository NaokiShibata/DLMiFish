let pyodide = null;
let pyodideReady = null;

const PYODIDE_VERSION = "0.26.4";
const PYODIDE_BASE_URL = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;

function postStatus(message) {
  self.postMessage({ type: "status", message });
}

function normalizeRelPath(input) {
  const parts = String(input || "")
    .replace(/\\/g, "/")
    .split("/")
    .filter((part) => part && part !== "." && part !== "..");
  return parts.join("/");
}

function ensureDir(fs, absPath) {
  const segments = absPath.split("/").filter(Boolean);
  let current = "";
  for (const segment of segments) {
    current += `/${segment}`;
    try {
      fs.mkdir(current);
    } catch (err) {
      // ignore EEXIST
    }
  }
}

function guessMime(path) {
  const lower = String(path).toLowerCase();
  if (lower.endsWith(".csv")) return "text/csv";
  if (lower.endsWith(".log")) return "text/plain";
  if (lower.endsWith(".fa") || lower.endsWith(".fas") || lower.endsWith(".fasta")) return "text/plain";
  if (lower.endsWith(".toml")) return "text/plain";
  return "application/octet-stream";
}

async function ensurePyodide() {
  if (pyodideReady) return pyodideReady;

  pyodideReady = (async () => {
    postStatus("Loading Pyodide runtime...");
    importScripts(`${PYODIDE_BASE_URL}pyodide.js`);
    pyodide = await self.loadPyodide({
      indexURL: PYODIDE_BASE_URL,
    });

    if (typeof pyodide.setStdout === "function") {
      pyodide.setStdout({
        batched: (msg) => self.postMessage({ type: "log", message: String(msg) }),
      });
    }
    if (typeof pyodide.setStderr === "function") {
      pyodide.setStderr({
        batched: (msg) => self.postMessage({ type: "log", message: `[stderr] ${String(msg)}` }),
      });
    }

    postStatus("Loading Python dependencies...");
    await pyodide.loadPackage("micropip");
    const micropip = pyodide.pyimport("micropip");

    const packages = [
      { pip: "rich", importName: "rich" },
      { pip: "typer", importName: "typer" },
      { pip: "tomli", importName: "tomli" },
      { pip: "biopython", importName: "Bio" },
    ];
    for (const pkg of packages) {
      try {
        await pyodide.runPythonAsync(`import ${pkg.importName}`);
      } catch (err) {
        postStatus(`Installing ${pkg.pip}...`);
        await micropip.install(pkg.pip);
      }
    }
    if (typeof micropip.destroy === "function") {
      micropip.destroy();
    }

    postStatus("Loading taxondbbuilder.py...");
    const sourceResp = await fetch("../taxondbbuilder.py");
    if (!sourceResp.ok) {
      throw new Error(`Failed to fetch taxondbbuilder.py: HTTP ${sourceResp.status}`);
    }
    const sourceText = await sourceResp.text();
    pyodide.FS.writeFile("/taxondbbuilder.py", sourceText);
    await pyodide.runPythonAsync(`
import sys
if "/" not in sys.path:
    sys.path.insert(0, "/")
import taxondbbuilder
`);
  })();

  return pyodideReady;
}

function writeInputFiles(runRoot, files) {
  const fs = pyodide.FS;
  ensureDir(fs, runRoot);
  for (const file of files || []) {
    const rel = normalizeRelPath(file.path);
    if (!rel) continue;
    const abs = `${runRoot}/${rel}`;
    const parent = abs.split("/").slice(0, -1).join("/");
    ensureDir(fs, parent);
    const bytes = new Uint8Array(file.buffer);
    fs.writeFile(abs, bytes);
  }
}

async function runBuild(payload, files) {
  await ensurePyodide();
  const runRoot = `/workspace/run_${Date.now()}`;
  writeInputFiles(runRoot, files);

  const buildArgs = {
    run_root: runRoot,
    config_path: payload.config_path || "configs/db.toml",
    taxa: Array.isArray(payload.taxa) ? payload.taxa : [],
    markers: Array.isArray(payload.markers) ? payload.markers : [],
    remote_entrez_api: payload.remote_entrez_api || "",
    workers: Number(payload.workers || 2),
    out: payload.out || "Results/db/browser_output.fasta",
    output_prefix: payload.output_prefix || "taxondbbuilder_",
    dump_gb: payload.dump_gb || "",
    from_gb: payload.from_gb || "",
    resume: Boolean(payload.resume),
    dry_run: Boolean(payload.dry_run),
    post_prep: Boolean(payload.post_prep),
    post_prep_step: Array.isArray(payload.post_prep_step) ? payload.post_prep_step : [],
    post_prep_primer_set: Array.isArray(payload.post_prep_primer_set) ? payload.post_prep_primer_set : [],
  };

  pyodide.globals.set("build_args", buildArgs);
  postStatus("Running build in Python worker...");

  await pyodide.runPythonAsync(`
from pathlib import Path
import taxondbbuilder as t

args = build_args.to_py()
run_root = Path(args["run_root"])

def opt_path(key):
    value = args.get(key)
    if value is None:
        return None
    value = str(value).strip()
    if not value:
        return None
    return run_root / value

post_steps_raw = list(args.get("post_prep_step") or [])
post_steps = [t.PostPrepStep(step) for step in post_steps_raw] if post_steps_raw else None
post_primer_sets = list(args.get("post_prep_primer_set") or []) or None

t.build(
    config=run_root / str(args["config_path"]),
    taxon=list(args.get("taxa") or []),
    marker=list(args.get("markers") or []),
    out=opt_path("out"),
    dump_gb=opt_path("dump_gb"),
    from_gb=opt_path("from_gb"),
    remote_entrez_api=(str(args.get("remote_entrez_api") or "").strip() or None),
    resume=bool(args.get("resume", False)),
    dry_run=bool(args.get("dry_run", False)),
    workers=int(args.get("workers", 2)),
    output_prefix=str(args.get("output_prefix") or "taxondbbuilder_"),
    post_prep=bool(args.get("post_prep", False)),
    post_prep_step=post_steps,
    post_prep_primer_set=post_primer_sets,
)

out_path = opt_path("out")
generated = []
if out_path is not None:
    candidates = [
        out_path,
        out_path.with_suffix(out_path.suffix + ".log"),
        out_path.with_suffix(out_path.suffix + ".acc_organism.csv"),
        out_path.with_suffix(out_path.suffix + ".duplicate_acc.records.csv"),
        out_path.with_suffix(out_path.suffix + ".duplicate_acc.groups.csv"),
    ]
    for path in candidates:
        if path.exists():
            generated.append(str(path))

build_result = {
    "generated_files": generated,
    "run_root": str(run_root),
}
`);

  const resultProxy = pyodide.globals.get("build_result");
  const result = resultProxy.toJs({ dict_converter: Object.fromEntries });
  if (typeof resultProxy.destroy === "function") {
    resultProxy.destroy();
  }
  pyodide.globals.delete("build_args");
  pyodide.globals.delete("build_result");

  const outputs = [];
  const transfer = [];
  const runRootPrefix = `${result.run_root}/`;
  for (const absPath of result.generated_files || []) {
    const exists = pyodide.FS.analyzePath(absPath).exists;
    if (!exists) continue;
    const data = pyodide.FS.readFile(absPath, { encoding: "binary" });
    const buffer = data.buffer.slice(data.byteOffset, data.byteOffset + data.byteLength);
    outputs.push({
      path: absPath.startsWith(runRootPrefix) ? absPath.slice(runRootPrefix.length) : absPath,
      name: absPath.split("/").pop(),
      mime: guessMime(absPath),
      buffer,
    });
    transfer.push(buffer);
  }

  self.postMessage({ type: "result", files: outputs }, transfer);
}

self.onmessage = async (event) => {
  const msg = event.data || {};
  try {
    if (msg.type === "init") {
      await ensurePyodide();
      self.postMessage({ type: "ready" });
      return;
    }
    if (msg.type === "run") {
      await runBuild(msg.payload || {}, msg.files || []);
      return;
    }
  } catch (err) {
    self.postMessage({
      type: "error",
      error: err && err.message ? err.message : String(err),
    });
  }
};
