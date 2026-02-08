#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use chrono::Local;
use once_cell::sync::Lazy;
use regex::Regex;
use rfd::FileDialog;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;
use tauri::{AppHandle, Emitter, State};

const RUN_EVENT: &str = "run-event";
const CONFIG_DIR_NAME: &str = ".taxondb_gui";
const CONFIG_FILE_NAME: &str = "config.json";

static MARKERS_TEMPLATE: &str = include_str!("../../resources/templates/markers_mitogenome.toml");
static PRIMERS_TEMPLATE: &str = include_str!("../../resources/templates/primers.toml");

static RE_QUERY_COUNT: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^# query count taxid=([^:]+):\s*(.+)$").expect("query count regex"));
static RE_TOTAL_RECORDS: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^# total records:\s*(\d+)").expect("total records regex"));
static RE_MATCHED_RECORDS: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"^# matched records:\s*(\d+)").expect("matched records regex"));
static RE_KEPT_BEFORE_POST: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r"^# kept records before post_prep:\s*(\d+)").expect("kept before post regex")
});
static RE_PRIMER_REMOVED: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"removed=(\d+)").expect("primer removed regex"));
static RE_DUP_GROUPS: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"groups=(\d+)").expect("duplicate groups regex"));
static RE_DUP_CROSS: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"cross_organism_groups=(\d+)").expect("cross organism groups regex"));

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(rename_all = "camelCase")]
struct GuiConfig {
    email: String,
    api_key: String,
    save_api_key: bool,
    output_root: String,
    output_prefix: String,
    marker: String,
    workers: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(rename_all = "camelCase")]
struct FiltersInput {
    mitochondrion: bool,
    ddbj_embl_genbank: bool,
    biomol_genomic: bool,
    length_min: Option<u32>,
    length_max: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(rename_all = "camelCase")]
struct PostPrepInput {
    enable: bool,
    primer_file: String,
    primer_set: Vec<String>,
    steps: Vec<String>,
    sequence_length_min: Option<u32>,
    sequence_length_max: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
struct RunRequest {
    taxids: Vec<String>,
    markers: Vec<String>,
    output_prefix: String,
    output_root: String,
    email: String,
    api_key: String,
    save_api_key: bool,
    filters: FiltersInput,
    post_prep: PostPrepInput,
    workers: u32,
    resume: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(rename_all = "camelCase")]
struct ImportedDbTomlConfig {
    source_path: String,
    email: String,
    api_key: String,
    filters: FiltersInput,
    post_prep: PostPrepInput,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct StartRunResponse {
    job_id: String,
    job_dir: String,
    config_path: String,
    output_path: String,
    log_path: String,
    command: String,
}

#[derive(Debug, Clone, Serialize, Default)]
#[serde(rename_all = "camelCase")]
struct RunMetrics {
    query_count_by_taxid: BTreeMap<String, u64>,
    matched_records: Option<u64>,
    kept_records_before_post_prep: Option<u64>,
    primer_trim_removed: Option<u64>,
    length_filter_removed: Option<u64>,
    duplicate_groups: Option<u64>,
    cross_organism_groups: Option<u64>,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct RunEvent {
    event_type: String,
    status: Option<String>,
    phase: Option<String>,
    percent: Option<f64>,
    line: Option<String>,
    message: Option<String>,
    metrics: Option<RunMetrics>,
    files: Option<Vec<String>>,
    job_dir: Option<String>,
}

#[derive(Debug, Clone)]
struct ActiveRun {
    child: Arc<Mutex<Option<Child>>>,
    cancelled: Arc<AtomicBool>,
}

#[derive(Default)]
struct AppState {
    run: Arc<Mutex<Option<ActiveRun>>>,
}

#[derive(Debug, Clone)]
struct ProgressParser {
    taxid_total: usize,
    query_seen: HashSet<String>,
    total_records: Option<u64>,
    matched_records: Option<u64>,
    post_steps_total: usize,
    post_steps_seen: HashSet<String>,
    phase: String,
    percent: f64,
    metrics: RunMetrics,
}

impl ProgressParser {
    fn new(taxid_total: usize, post_steps_total: usize) -> Self {
        Self {
            taxid_total,
            query_seen: HashSet::new(),
            total_records: None,
            matched_records: None,
            post_steps_total,
            post_steps_seen: HashSet::new(),
            phase: "Query count".to_string(),
            percent: 0.0,
            metrics: RunMetrics::default(),
        }
    }

    fn consume_line(&mut self, line: &str) -> bool {
        let mut changed = false;

        if let Some(caps) = RE_QUERY_COUNT.captures(line) {
            let taxid = caps.get(1).map(|m| m.as_str()).unwrap_or("").trim();
            let count_raw = caps.get(2).map(|m| m.as_str()).unwrap_or("").trim();
            let count = count_raw.parse::<u64>().unwrap_or(0);
            if !taxid.is_empty() {
                self.metrics
                    .query_count_by_taxid
                    .insert(taxid.to_string(), count);
                self.query_seen.insert(taxid.to_string());
                self.phase = "Query count".to_string();
                let total = self.taxid_total.max(1) as f64;
                self.percent = ((self.query_seen.len() as f64 / total) * 10.0).clamp(0.0, 10.0);
                changed = true;
            }
        }

        if let Some(caps) = RE_TOTAL_RECORDS.captures(line) {
            self.total_records = caps.get(1).and_then(|m| m.as_str().parse::<u64>().ok());
            changed = true;
        }

        if let Some(caps) = RE_MATCHED_RECORDS.captures(line) {
            let matched = caps.get(1).and_then(|m| m.as_str().parse::<u64>().ok());
            self.matched_records = matched;
            self.metrics.matched_records = matched;
            changed = true;
        }

        if let Some(caps) = RE_KEPT_BEFORE_POST.captures(line) {
            self.metrics.kept_records_before_post_prep =
                caps.get(1).and_then(|m| m.as_str().parse::<u64>().ok());
            self.phase = "Post-pPrep".to_string();
            self.percent = self.percent.max(80.0);
            changed = true;
        }

        if line.starts_with("# post_prep primer trim:") {
            self.post_steps_seen.insert("primer_trim".to_string());
            self.metrics.primer_trim_removed = RE_PRIMER_REMOVED
                .captures(line)
                .and_then(|c| c.get(1))
                .and_then(|m| m.as_str().parse::<u64>().ok());
            changed = true;
        }

        if line.starts_with("# post_prep length filter:") {
            self.post_steps_seen.insert("length_filter".to_string());
            self.metrics.length_filter_removed = RE_PRIMER_REMOVED
                .captures(line)
                .and_then(|c| c.get(1))
                .and_then(|m| m.as_str().parse::<u64>().ok());
            changed = true;
        }

        if line.starts_with("# post_prep duplicate_acc_report:") {
            self.post_steps_seen.insert("duplicate_report".to_string());
            self.metrics.duplicate_groups = RE_DUP_GROUPS
                .captures(line)
                .and_then(|c| c.get(1))
                .and_then(|m| m.as_str().parse::<u64>().ok());
            self.metrics.cross_organism_groups = RE_DUP_CROSS
                .captures(line)
                .and_then(|c| c.get(1))
                .and_then(|m| m.as_str().parse::<u64>().ok());
            changed = true;
        }

        if line.starts_with("# output:") {
            self.phase = "Finalize".to_string();
            self.percent = self.percent.max(95.0);
            changed = true;
        }

        if line.starts_with("# finished:") {
            self.phase = "Finalize".to_string();
            self.percent = 100.0;
            changed = true;
        }

        if changed {
            if let (Some(total), Some(matched)) = (self.total_records, self.matched_records) {
                if total > 0 {
                    self.phase = "Fetch/Parse".to_string();
                    let ratio = (matched as f64 / total as f64).clamp(0.0, 1.0);
                    self.percent = (10.0 + ratio * 70.0)
                        .clamp(10.0, 80.0)
                        .max(self.percent.min(80.0));
                }
            }

            if self.post_steps_total > 0 && !self.post_steps_seen.is_empty() {
                self.phase = "Post-pPrep".to_string();
                let ratio = (self.post_steps_seen.len() as f64 / self.post_steps_total as f64)
                    .clamp(0.0, 1.0);
                self.percent = self.percent.max(80.0 + ratio * 15.0);
            }
        }

        changed
    }
}

fn emit_event(app: &AppHandle, payload: RunEvent) {
    let _ = app.emit(RUN_EVENT, payload);
}

fn status_event(status: &str) -> RunEvent {
    RunEvent {
        event_type: "status".to_string(),
        status: Some(status.to_string()),
        phase: None,
        percent: None,
        line: None,
        message: None,
        metrics: None,
        files: None,
        job_dir: None,
    }
}

fn log_event(line: String) -> RunEvent {
    RunEvent {
        event_type: "log".to_string(),
        status: None,
        phase: None,
        percent: None,
        line: Some(line),
        message: None,
        metrics: None,
        files: None,
        job_dir: None,
    }
}

fn progress_event(parser: &ProgressParser) -> RunEvent {
    RunEvent {
        event_type: "progress".to_string(),
        status: None,
        phase: Some(parser.phase.clone()),
        percent: Some(parser.percent.min(100.0)),
        line: None,
        message: None,
        metrics: Some(parser.metrics.clone()),
        files: None,
        job_dir: None,
    }
}

fn result_event(job_dir: &Path, files: Vec<String>) -> RunEvent {
    RunEvent {
        event_type: "result".to_string(),
        status: None,
        phase: None,
        percent: None,
        line: None,
        message: None,
        metrics: None,
        files: Some(files),
        job_dir: Some(job_dir.to_string_lossy().to_string()),
    }
}

fn error_event(message: String) -> RunEvent {
    RunEvent {
        event_type: "error".to_string(),
        status: None,
        phase: None,
        percent: None,
        line: None,
        message: Some(message),
        metrics: None,
        files: None,
        job_dir: None,
    }
}

fn to_abs_string(path: &Path) -> String {
    path.to_string_lossy().to_string()
}

fn gui_config_path() -> Result<PathBuf, String> {
    let home = dirs::home_dir().ok_or_else(|| "HOME directory not found".to_string())?;
    Ok(home.join(CONFIG_DIR_NAME).join(CONFIG_FILE_NAME))
}

fn ensure_gui_config_parent() -> Result<PathBuf, String> {
    let path = gui_config_path()?;
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("failed to create {}: {e}", parent.display()))?;
    }
    Ok(path)
}

fn sanitize_file_name(raw: &str) -> String {
    let mut s = raw.trim().replace(' ', "_");
    s = s
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if s.is_empty() {
        "taxondbbuilder".to_string()
    } else {
        s
    }
}

fn prepare_job_dir(output_root: &Path) -> Result<(String, PathBuf), String> {
    let date_dir = output_root.join(Local::now().format("%Y%m%d").to_string());
    fs::create_dir_all(&date_dir)
        .map_err(|e| format!("failed to create {}: {e}", date_dir.display()))?;

    for idx in 1..10000 {
        let candidate = date_dir.join(format!("job{idx}"));
        if !candidate.exists() {
            fs::create_dir_all(&candidate)
                .map_err(|e| format!("failed to create {}: {e}", candidate.display()))?;
            return Ok((
                format!("{}-job{idx}", Local::now().format("%Y%m%d")),
                candidate,
            ));
        }
    }

    Err("could not allocate job directory".to_string())
}

fn resolve_sidecar_path() -> Result<PathBuf, String> {
    if let Ok(path) = std::env::var("TAXONDBBUILDER_BIN") {
        let p = PathBuf::from(path);
        if p.exists() {
            return Ok(p);
        }
    }

    let target_triple = option_env!("TAURI_ENV_TARGET_TRIPLE").unwrap_or("unknown-target");
    let mut search_dirs: Vec<PathBuf> = Vec::new();
    if let Ok(cwd) = std::env::current_dir() {
        search_dirs.push(cwd.join("bin"));
        search_dirs.push(cwd.join("src-tauri/bin"));
        search_dirs.push(cwd.join("tauri-gui/src-tauri/bin"));
    }

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    search_dirs.push(manifest_dir.join("bin"));
    if let Some(parent) = manifest_dir.parent() {
        search_dirs.push(parent.join("bin"));
        search_dirs.push(parent.join("src-tauri/bin"));
    }

    if let Ok(exe) = std::env::current_exe() {
        if let Some(parent) = exe.parent() {
            search_dirs.push(parent.to_path_buf());
            search_dirs.push(parent.join("bin"));
        }
    }

    search_dirs.sort();
    search_dirs.dedup();

    let mut candidate_names: Vec<String> = Vec::new();
    #[cfg(target_os = "windows")]
    {
        candidate_names.push(format!("taxondbbuilder-{target_triple}.exe"));
        candidate_names.push("taxondbbuilder.exe".to_string());
    }
    #[cfg(not(target_os = "windows"))]
    {
        candidate_names.push(format!("taxondbbuilder-{target_triple}"));
        candidate_names.push("taxondbbuilder".to_string());
    }

    for dir in &search_dirs {
        for name in &candidate_names {
            let path = dir.join(name);
            if path.exists() && path.is_file() {
                return Ok(path);
            }
        }
    }

    for dir in &search_dirs {
        let Ok(entries) = fs::read_dir(dir) else {
            continue;
        };
        let mut matches: Vec<PathBuf> = entries
            .flatten()
            .map(|entry| entry.path())
            .filter(|path| path.is_file())
            .filter(|path| {
                let Some(name) = path.file_name().and_then(|v| v.to_str()) else {
                    return false;
                };
                if !name.starts_with("taxondbbuilder-") {
                    return false;
                }
                #[cfg(target_os = "windows")]
                {
                    return name.ends_with(".exe");
                }
                #[cfg(not(target_os = "windows"))]
                {
                    true
                }
            })
            .collect();
        matches.sort();
        if let Some(path) = matches.first() {
            return Ok(path.clone());
        }
    }

    if let Ok(path) = which::which("taxondbbuilder") {
        return Ok(path);
    }

    let dirs = search_dirs
        .iter()
        .map(|d| d.to_string_lossy().to_string())
        .collect::<Vec<_>>()
        .join(", ");
    Err(format!(
        "taxondbbuilder sidecar not found. Set TAXONDBBUILDER_BIN or place binary under src-tauri/bin/ with Tauri target-triple suffix. searched=[{dirs}]"
    ))
}

fn toml_quote(s: &str) -> String {
    format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
}

fn write_job_config(req: &RunRequest, config_dir: &Path) -> Result<PathBuf, String> {
    fs::create_dir_all(config_dir)
        .map_err(|e| format!("failed to create {}: {e}", config_dir.display()))?;

    let markers_path = config_dir.join("markers_mitogenome.toml");
    fs::write(&markers_path, MARKERS_TEMPLATE)
        .map_err(|e| format!("failed to write {}: {e}", markers_path.display()))?;

    let bundled_primers_path = config_dir.join("primers.toml");
    fs::write(&bundled_primers_path, PRIMERS_TEMPLATE)
        .map_err(|e| format!("failed to write {}: {e}", bundled_primers_path.display()))?;

    let mut filters: Vec<String> = Vec::new();
    if req.filters.mitochondrion {
        filters.push("mitochondrion".to_string());
    }
    if req.filters.ddbj_embl_genbank {
        filters.push("ddbj_embl_genbank".to_string());
    }

    let mut props: Vec<String> = Vec::new();
    if req.filters.biomol_genomic {
        props.push("biomol_genomic".to_string());
    }

    let primer_file = if req.post_prep.primer_file.trim().is_empty() {
        bundled_primers_path.to_string_lossy().to_string()
    } else {
        req.post_prep.primer_file.trim().to_string()
    };

    let mut toml = String::new();
    toml.push_str("[ncbi]\n");
    toml.push_str(&format!("email = {}\n", toml_quote(req.email.trim())));
    if !req.api_key.trim().is_empty() {
        toml.push_str(&format!("api_key = {}\n", toml_quote(req.api_key.trim())));
    }
    toml.push_str("db = \"nucleotide\"\n");
    toml.push_str("rettype = \"gb\"\n");
    toml.push_str("retmode = \"text\"\n");
    toml.push_str("per_query = 100\n");
    toml.push_str("use_history = true\n\n");

    toml.push_str("[output]\n");
    toml.push_str("default_header_format = \"{acc_id}|{organism}|{marker}|{label}|{type}|{loc}|{strand}\"\n\n");
    toml.push_str("[output.header_formats]\n");
    toml.push_str("mifish_pipeline = \"gb|{acc_id}|{organism}\"\n\n");

    toml.push_str("[taxon]\n");
    toml.push_str("noexp = false\n\n");

    toml.push_str("[markers]\n");
    toml.push_str("file = \"markers_mitogenome.toml\"\n\n");

    toml.push_str("[filters]\n");
    if !filters.is_empty() {
        let values = filters
            .iter()
            .map(|v| toml_quote(v))
            .collect::<Vec<_>>()
            .join(", ");
        toml.push_str(&format!("filter = [{}]\n", values));
    }
    if !props.is_empty() {
        let values = props
            .iter()
            .map(|v| toml_quote(v))
            .collect::<Vec<_>>()
            .join(", ");
        toml.push_str(&format!("properties = [{}]\n", values));
    }
    if let Some(v) = req.filters.length_min {
        toml.push_str(&format!("sequence_length_min = {}\n", v));
    }
    if let Some(v) = req.filters.length_max {
        toml.push_str(&format!("sequence_length_max = {}\n", v));
    }
    toml.push('\n');

    toml.push_str("[post_prep]\n");
    toml.push_str(&format!("primer_file = {}\n", toml_quote(&primer_file)));
    if !req.post_prep.primer_set.is_empty() {
        let sets = req
            .post_prep
            .primer_set
            .iter()
            .map(|s| toml_quote(s))
            .collect::<Vec<_>>()
            .join(", ");
        toml.push_str(&format!("primer_set = [{}]\n", sets));
    }
    if let Some(v) = req.post_prep.sequence_length_min {
        toml.push_str(&format!("sequence_length_min = {}\n", v));
    }
    if let Some(v) = req.post_prep.sequence_length_max {
        toml.push_str(&format!("sequence_length_max = {}\n", v));
    }

    let config_path = config_dir.join("db.toml");
    fs::write(&config_path, toml)
        .map_err(|e| format!("failed to write {}: {e}", config_path.display()))?;

    Ok(config_path)
}

fn collect_files(results_dir: &Path, extra: &[PathBuf]) -> Vec<String> {
    let mut out = Vec::new();
    if let Ok(entries) = fs::read_dir(results_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_file() {
                out.push(path.to_string_lossy().to_string());
            }
        }
    }
    for p in extra {
        if p.exists() {
            out.push(p.to_string_lossy().to_string());
        }
    }
    out.sort();
    out.dedup();
    out
}

fn read_pipe_lines<R: Read + Send + 'static>(reader: R, source: &'static str, app: AppHandle) {
    thread::spawn(move || {
        let buffered = BufReader::new(reader);
        for line in buffered.lines().map_while(Result::ok) {
            emit_event(&app, log_event(format!("[{source}] {line}")));
        }
    });
}

fn tail_log_once(
    app: &AppHandle,
    parser: &mut ProgressParser,
    log_path: &Path,
    offset: &mut u64,
) -> Result<(), String> {
    if !log_path.exists() {
        return Ok(());
    }

    let mut file = File::open(log_path)
        .map_err(|e| format!("failed to open log {}: {e}", log_path.display()))?;
    file.seek(SeekFrom::Start(*offset))
        .map_err(|e| format!("failed to seek log {}: {e}", log_path.display()))?;

    let mut data = Vec::new();
    file.read_to_end(&mut data)
        .map_err(|e| format!("failed to read log {}: {e}", log_path.display()))?;

    if data.is_empty() {
        return Ok(());
    }

    *offset += data.len() as u64;
    let text = String::from_utf8_lossy(&data);

    for line in text.lines() {
        emit_event(app, log_event(line.to_string()));
        if parser.consume_line(line) {
            emit_event(app, progress_event(parser));
        }
    }

    Ok(())
}

fn save_gui_config_internal(config: &GuiConfig) -> Result<(), String> {
    let path = ensure_gui_config_parent()?;
    let mut saved = config.clone();
    if !saved.save_api_key {
        saved.api_key.clear();
    }
    let json = serde_json::to_string_pretty(&saved)
        .map_err(|e| format!("failed to serialize config: {e}"))?;
    fs::write(&path, json).map_err(|e| format!("failed to write {}: {e}", path.display()))
}

fn toml_string_list(value: Option<&toml::Value>) -> Vec<String> {
    let Some(value) = value else {
        return Vec::new();
    };
    if let Some(s) = value.as_str() {
        let trimmed = s.trim();
        if trimmed.is_empty() {
            return Vec::new();
        }
        return vec![trimmed.to_string()];
    }
    if let Some(arr) = value.as_array() {
        let mut out = Vec::new();
        for item in arr {
            if let Some(s) = item.as_str() {
                let trimmed = s.trim();
                if !trimmed.is_empty() {
                    out.push(trimmed.to_string());
                }
            }
        }
        return out;
    }
    Vec::new()
}

fn toml_u32(value: Option<&toml::Value>) -> Option<u32> {
    value
        .and_then(|v| v.as_integer())
        .and_then(|v| u32::try_from(v).ok())
}

fn parse_db_toml_config(path: &Path) -> Result<ImportedDbTomlConfig, String> {
    let text =
        fs::read_to_string(path).map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    let value: toml::Value = text
        .parse()
        .map_err(|e| format!("failed to parse {}: {e}", path.display()))?;

    let mut imported = ImportedDbTomlConfig {
        source_path: path.to_string_lossy().to_string(),
        ..ImportedDbTomlConfig::default()
    };

    if let Some(ncbi) = value.get("ncbi").and_then(|v| v.as_table()) {
        imported.email = ncbi
            .get("email")
            .and_then(|v| v.as_str())
            .unwrap_or("")
            .trim()
            .to_string();
        imported.api_key = ncbi
            .get("api_key")
            .and_then(|v| v.as_str())
            .unwrap_or("")
            .trim()
            .to_string();
    }

    if let Some(filters) = value.get("filters").and_then(|v| v.as_table()) {
        let filter_terms = toml_string_list(filters.get("filter"));
        imported.filters.mitochondrion = filter_terms.iter().any(|v| v == "mitochondrion");
        imported.filters.ddbj_embl_genbank = filter_terms.iter().any(|v| v == "ddbj_embl_genbank");

        let properties = toml_string_list(filters.get("properties"));
        imported.filters.biomol_genomic = properties.iter().any(|v| v == "biomol_genomic");

        imported.filters.length_min = toml_u32(filters.get("sequence_length_min"));
        imported.filters.length_max = toml_u32(filters.get("sequence_length_max"));
    }

    if let Some(post) = value.get("post_prep").and_then(|v| v.as_table()) {
        imported.post_prep.enable = true;
        imported.post_prep.primer_file = post
            .get("primer_file")
            .and_then(|v| v.as_str())
            .unwrap_or("")
            .trim()
            .to_string();
        imported.post_prep.primer_set = toml_string_list(post.get("primer_set"));
        imported.post_prep.sequence_length_min = toml_u32(post.get("sequence_length_min"));
        imported.post_prep.sequence_length_max = toml_u32(post.get("sequence_length_max"));

        let mut steps: Vec<String> = Vec::new();
        if !imported.post_prep.primer_file.is_empty() && !imported.post_prep.primer_set.is_empty() {
            steps.push("primer_trim".to_string());
        }
        if imported.post_prep.sequence_length_min.is_some() || imported.post_prep.sequence_length_max.is_some() {
            steps.push("length_filter".to_string());
        }
        steps.push("duplicate_report".to_string());
        imported.post_prep.steps = steps;
    }

    Ok(imported)
}

#[tauri::command]
fn load_gui_config() -> Result<GuiConfig, String> {
    let path = ensure_gui_config_parent()?;
    if !path.exists() {
        return Ok(GuiConfig::default());
    }

    let text =
        fs::read_to_string(&path).map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    let mut cfg: GuiConfig = serde_json::from_str(&text)
        .map_err(|e| format!("failed to parse {}: {e}", path.display()))?;
    if !cfg.save_api_key {
        cfg.api_key.clear();
    }
    Ok(cfg)
}

#[tauri::command]
fn save_gui_config(config: GuiConfig) -> Result<(), String> {
    save_gui_config_internal(&config)
}

#[tauri::command]
fn choose_output_directory() -> Option<String> {
    FileDialog::new()
        .pick_folder()
        .map(|p| p.to_string_lossy().to_string())
}

#[tauri::command]
fn choose_primer_file() -> Option<String> {
    FileDialog::new()
        .add_filter("TOML", &["toml"])
        .pick_file()
        .map(|p| p.to_string_lossy().to_string())
}

#[tauri::command]
fn choose_db_toml_file() -> Option<String> {
    FileDialog::new()
        .add_filter("TOML", &["toml"])
        .set_title("Select db.toml")
        .pick_file()
        .map(|p| p.to_string_lossy().to_string())
}

#[tauri::command]
fn import_db_toml(path: String) -> Result<ImportedDbTomlConfig, String> {
    let target = PathBuf::from(path.trim());
    if !target.exists() {
        return Err(format!("db.toml not found: {}", target.display()));
    }
    parse_db_toml_config(&target)
}

#[tauri::command]
fn open_path(path: String) -> Result<(), String> {
    open::that(path).map_err(|e| format!("failed to open path: {e}"))
}

#[tauri::command]
fn cancel_run(state: State<AppState>) -> Result<(), String> {
    let run = {
        let slot = state
            .run
            .lock()
            .map_err(|_| "failed to lock run state".to_string())?;
        slot.clone()
    };

    if let Some(active) = run {
        active.cancelled.store(true, Ordering::Relaxed);
        let mut guard = active
            .child
            .lock()
            .map_err(|_| "failed to lock child process".to_string())?;
        if let Some(child) = guard.as_mut() {
            child
                .kill()
                .map_err(|e| format!("failed to kill process: {e}"))?;
        }
        Ok(())
    } else {
        Err("no running job".to_string())
    }
}

#[tauri::command]
fn start_run(
    app: AppHandle,
    state: State<AppState>,
    req: RunRequest,
) -> Result<StartRunResponse, String> {
    if req.taxids.is_empty() {
        return Err("taxids must not be empty".to_string());
    }
    if req.markers.is_empty() {
        return Err("markers must not be empty".to_string());
    }
    if req.output_root.trim().is_empty() {
        return Err("output_root is required".to_string());
    }
    if req.email.trim().is_empty() {
        return Err("email is required".to_string());
    }

    {
        let slot = state
            .run
            .lock()
            .map_err(|_| "failed to lock run state".to_string())?;
        if slot.is_some() {
            return Err("another job is already running".to_string());
        }
    }

    let gui_config = GuiConfig {
        email: req.email.clone(),
        api_key: req.api_key.clone(),
        save_api_key: req.save_api_key,
        output_root: req.output_root.clone(),
        output_prefix: req.output_prefix.clone(),
        marker: req.markers.first().cloned().unwrap_or_default(),
        workers: req.workers,
    };
    save_gui_config_internal(&gui_config)?;

    let output_root = PathBuf::from(req.output_root.trim());
    fs::create_dir_all(&output_root)
        .map_err(|e| format!("failed to create {}: {e}", output_root.display()))?;

    let (job_id, job_dir) = prepare_job_dir(&output_root)?;
    let config_dir = job_dir.join("config");
    let gb_dir = job_dir.join("gb");
    let results_dir = job_dir.join("Results");
    fs::create_dir_all(&gb_dir)
        .map_err(|e| format!("failed to create {}: {e}", gb_dir.display()))?;
    fs::create_dir_all(&results_dir)
        .map_err(|e| format!("failed to create {}: {e}", results_dir.display()))?;

    let config_path = write_job_config(&req, &config_dir)?;

    let output_prefix = sanitize_file_name(&req.output_prefix);
    let output_file = results_dir.join(format!(
        "{}_{}.fasta",
        output_prefix,
        Local::now().format("%Y%m%d%H%M%S")
    ));
    let log_path = PathBuf::from(format!("{}.log", output_file.to_string_lossy()));

    let sidecar = resolve_sidecar_path()?;

    let mut args: Vec<String> = vec!["build".to_string()];
    args.push("-c".to_string());
    args.push(to_abs_string(&config_path));

    for taxid in &req.taxids {
        args.push("-t".to_string());
        args.push(taxid.trim().to_string());
    }

    for marker in &req.markers {
        args.push("-m".to_string());
        args.push(marker.trim().to_string());
    }

    args.push("--out".to_string());
    args.push(to_abs_string(&output_file));
    args.push("--dump-gb".to_string());
    args.push(to_abs_string(&gb_dir));
    args.push("--output-prefix".to_string());
    args.push(req.output_prefix.clone());
    args.push("--workers".to_string());
    args.push(req.workers.max(1).to_string());

    if req.resume {
        args.push("--resume".to_string());
    }

    if req.post_prep.enable {
        args.push("--post-prep".to_string());
        for step in &req.post_prep.steps {
            args.push("--post-prep-step".to_string());
            args.push(step.clone());
        }
        for set in &req.post_prep.primer_set {
            args.push("--post-prep-primer-set".to_string());
            args.push(set.clone());
        }
    }

    let mut command = Command::new(&sidecar);
    command
        .args(&args)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .current_dir(&job_dir);

    let mut child = command
        .spawn()
        .map_err(|e| format!("failed to start sidecar {}: {e}", sidecar.display()))?;

    if let Some(stdout) = child.stdout.take() {
        read_pipe_lines(stdout, "stdout", app.clone());
    }
    if let Some(stderr) = child.stderr.take() {
        read_pipe_lines(stderr, "stderr", app.clone());
    }

    let child_arc = Arc::new(Mutex::new(Some(child)));
    let cancelled = Arc::new(AtomicBool::new(false));

    {
        let mut slot = state
            .run
            .lock()
            .map_err(|_| "failed to lock run state".to_string())?;
        *slot = Some(ActiveRun {
            child: child_arc.clone(),
            cancelled: cancelled.clone(),
        });
    }

    emit_event(&app, status_event("Running"));

    let run_slot = state.run.clone();
    let app_for_thread = app.clone();
    let log_path_for_thread = log_path.clone();
    let results_dir_for_thread = results_dir.clone();
    let job_dir_for_thread = job_dir.clone();
    let taxid_total = req.taxids.len();
    let post_steps_total = if req.post_prep.enable {
        req.post_prep.steps.len().max(1)
    } else {
        0
    };

    thread::spawn(move || {
        let mut parser = ProgressParser::new(taxid_total, post_steps_total);
        let mut offset: u64 = 0;
        let exit_code: i32 = loop {
            if let Err(err) = tail_log_once(
                &app_for_thread,
                &mut parser,
                &log_path_for_thread,
                &mut offset,
            ) {
                emit_event(&app_for_thread, error_event(err));
            }

            let status: Option<i32> = {
                let mut guard = child_arc.lock().ok();
                if let Some(ref mut guard) = guard {
                    if let Some(child) = guard.as_mut() {
                        child
                            .try_wait()
                            .ok()
                            .flatten()
                            .map(|s| s.code().unwrap_or(1))
                    } else {
                        Some(1)
                    }
                } else {
                    Some(1)
                }
            };

            if let Some(code) = status {
                break code;
            }

            thread::sleep(Duration::from_millis(300));
        };

        let _ = tail_log_once(
            &app_for_thread,
            &mut parser,
            &log_path_for_thread,
            &mut offset,
        );

        let was_cancelled = cancelled.load(Ordering::Relaxed);
        let files = collect_files(&results_dir_for_thread, &[log_path_for_thread.clone()]);

        match (was_cancelled, exit_code) {
            (true, _) => emit_event(&app_for_thread, status_event("Cancelled")),
            (false, 0) => emit_event(&app_for_thread, status_event("Finished")),
            _ => emit_event(&app_for_thread, status_event("Failed")),
        }

        if !was_cancelled {
            emit_event(&app_for_thread, result_event(&job_dir_for_thread, files));
        }

        if let Ok(mut slot) = run_slot.lock() {
            *slot = None;
        }
    });

    Ok(StartRunResponse {
        job_id,
        job_dir: to_abs_string(&job_dir),
        config_path: to_abs_string(&config_path),
        output_path: to_abs_string(&output_file),
        log_path: to_abs_string(&log_path),
        command: format!("{} {}", sidecar.display(), args.join(" ")),
    })
}

fn main() {
    tauri::Builder::default()
        .manage(AppState::default())
        .invoke_handler(tauri::generate_handler![
            load_gui_config,
            save_gui_config,
            choose_output_directory,
            choose_primer_file,
            choose_db_toml_file,
            import_db_toml,
            open_path,
            start_run,
            cancel_run
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
