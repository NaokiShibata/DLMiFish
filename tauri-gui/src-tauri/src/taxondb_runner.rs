use chrono::Local;
use regex::Regex;
use serde_json::Value as JsonValue;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};
use toml::Value as TomlValue;

const DEFAULT_FEATURE_TYPES: [&str; 3] = ["rRNA", "gene", "CDS"];
const DEFAULT_FEATURE_FIELDS: [&str; 4] = ["gene", "product", "note", "standard_name"];
const DEFAULT_HEADER_FORMAT: &str = "{acc_id}|{organism}|{marker}|{label}|{type}|{loc}|{strand}";
const DEFAULT_BOLD_HEADER_FORMAT: &str = "bold|{acc_id}|{organism}";
const BOLD_DEFAULT_BASE_URL: &str = "https://portal.boldsystems.org/api";
const BOLD_DEFAULT_TIMEOUT_SEC: f64 = 900.0;
const BOLD_DEFAULT_RETRIES: u64 = 3;
const BOLD_DEFAULT_BACKOFF_SEC: f64 = 1.5;
const BOLD_DEFAULT_DOWNLOAD_FORMAT: &str = "tsv";
const BOLD_DEFAULT_DOWNLOAD_CHUNK_SIZE: usize = 64 * 1024;
const BOLD_DEFAULT_USER_AGENT: &str =
    "TaxonDBBuilderGUI/0.1 (+https://github.com/NaokiShibata/TaxonDBBuilder)";
const BOLD_MAX_DOCUMENT_COUNT: u64 = 1_000_000;

#[derive(Debug, Clone)]
pub struct BuildParams {
    pub config_path: PathBuf,
    pub taxids: Vec<String>,
    pub markers: Vec<String>,
    pub source: String,
    pub output_file: PathBuf,
    pub dump_gb_dir: PathBuf,
    pub resume: bool,
}

#[derive(Debug, Default, Clone)]
struct Counters {
    total_records: u64,
    matched_records: u64,
    matched_features: u64,
    kept_records: u64,
    skipped_same: u64,
    duplicated_diff: u64,
}

#[derive(Debug, Clone)]
struct MarkerRule {
    key: String,
    patterns: Vec<Regex>,
    feature_types: Option<Vec<String>>,
    feature_fields: Vec<String>,
    header_format: String,
}

#[derive(Debug, Clone)]
struct ParsedFeature {
    feature_type: String,
    location_raw: String,
    qualifiers: HashMap<String, Vec<String>>,
}

#[derive(Debug, Clone)]
struct ParsedRecord {
    accession: String,
    organism: String,
    sequence: String,
    features: Vec<ParsedFeature>,
}

#[derive(Debug, Clone)]
struct EmittedRecord {
    acc_id: String,
    accession: String,
    organism_name: String,
    header: String,
    source: String,
    source_record_id: String,
    processid: String,
    sampleid: String,
    marker_key: String,
    linked_to_ncbi: bool,
    emitted_to_fasta: bool,
    skip_reason: String,
}

#[derive(Debug, Clone)]
struct ResolvedTaxon {
    taxid: String,
    scientific_name: String,
}

#[derive(Debug, Clone)]
struct BoldRuntimeConfig {
    base_url: String,
    timeout_sec: f64,
    retries: u64,
    backoff_sec: f64,
    download_format: String,
    download_chunk_size: usize,
    user_agent: String,
}

#[derive(Debug, Clone)]
struct BoldRecord {
    source_record_id: String,
    processid: String,
    sampleid: String,
    accession: String,
    taxon_name: String,
    marker_key: String,
    marker_label: String,
    sequence: String,
    header_format: String,
}

#[derive(Debug, Default, Clone)]
struct BoldFetchStats {
    normalized_query: String,
    specimen_count: u64,
    download_format: String,
    downloaded_bytes: u64,
    downloaded_rows: u64,
    matched_rows: u64,
}

#[derive(Debug, Clone)]
struct PreparedBoldQuery {
    normalized_query: String,
    specimen_count: Option<u64>,
    query_id: Option<String>,
    download_format: String,
}

#[derive(Debug, Clone)]
struct BoldDownloadMeta {
    downloaded_bytes: u64,
}

fn normalize_build_source(raw: &str) -> String {
    match raw.trim().to_ascii_lowercase().as_str() {
        "bold" => "bold".to_string(),
        "both" => "both".to_string(),
        _ => "ncbi".to_string(),
    }
}

fn source_uses_ncbi(source: &str) -> bool {
    normalize_build_source(source) != "bold"
}

fn source_uses_bold(source: &str) -> bool {
    matches!(normalize_build_source(source).as_str(), "bold" | "both")
}

fn normalized_key(value: &str) -> String {
    value
        .chars()
        .filter(|ch| ch.is_ascii_alphanumeric())
        .map(|ch| ch.to_ascii_lowercase())
        .collect()
}

fn find_values_by_key<'a>(
    node: &'a JsonValue,
    targets: &HashSet<String>,
    found: &mut Vec<&'a JsonValue>,
) {
    match node {
        JsonValue::Object(map) => {
            for (key, value) in map {
                if targets.contains(&normalized_key(key)) {
                    found.push(value);
                }
                find_values_by_key(value, targets, found);
            }
        }
        JsonValue::Array(items) => {
            for item in items {
                find_values_by_key(item, targets, found);
            }
        }
        _ => {}
    }
}

fn first_scalar(node: &JsonValue, keys: &[&str]) -> Option<String> {
    let targets = keys
        .iter()
        .map(|key| normalized_key(key))
        .collect::<HashSet<_>>();
    let mut found = Vec::new();
    find_values_by_key(node, &targets, &mut found);

    for value in found {
        match value {
            JsonValue::String(s) if !s.trim().is_empty() => return Some(s.trim().to_string()),
            JsonValue::Number(n) => return Some(n.to_string()),
            JsonValue::Array(items) => {
                for item in items {
                    match item {
                        JsonValue::String(s) if !s.trim().is_empty() => {
                            return Some(s.trim().to_string())
                        }
                        JsonValue::Number(n) => return Some(n.to_string()),
                        _ => {}
                    }
                }
            }
            _ => {}
        }
    }
    None
}

fn json_string(value: &JsonValue) -> Option<String> {
    match value {
        JsonValue::String(s) if !s.trim().is_empty() => Some(s.trim().to_string()),
        JsonValue::Number(n) => Some(n.to_string()),
        _ => None,
    }
}

fn get_bold_runtime_config(
    bold_tbl: Option<&toml::map::Map<String, TomlValue>>,
) -> Result<BoldRuntimeConfig, String> {
    let base_url = bold_tbl
        .and_then(|tbl| tbl.get("base_url"))
        .and_then(|v| v.as_str())
        .map(|s| s.trim().trim_end_matches('/').to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| BOLD_DEFAULT_BASE_URL.to_string());
    let timeout_sec = toml_f64_from_value(bold_tbl.and_then(|tbl| tbl.get("timeout_sec")))
        .unwrap_or(BOLD_DEFAULT_TIMEOUT_SEC);
    let retries = toml_u64_from_value(bold_tbl.and_then(|tbl| tbl.get("retries")))
        .unwrap_or(BOLD_DEFAULT_RETRIES);
    let backoff_sec = toml_f64_from_value(bold_tbl.and_then(|tbl| tbl.get("backoff_sec")))
        .unwrap_or(BOLD_DEFAULT_BACKOFF_SEC);
    let download_format = bold_tbl
        .and_then(|tbl| tbl.get("download_format"))
        .and_then(|v| v.as_str())
        .map(|s| s.trim().to_ascii_lowercase())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| BOLD_DEFAULT_DOWNLOAD_FORMAT.to_string());
    let download_chunk_size =
        toml_u64_from_value(bold_tbl.and_then(|tbl| tbl.get("download_chunk_size")))
            .map(|v| v as usize)
            .unwrap_or(BOLD_DEFAULT_DOWNLOAD_CHUNK_SIZE);
    let user_agent = bold_tbl
        .and_then(|tbl| tbl.get("user_agent"))
        .and_then(|v| v.as_str())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| BOLD_DEFAULT_USER_AGENT.to_string());

    if timeout_sec <= 0.0 {
        return Err("bold.timeout_sec must be > 0".to_string());
    }
    if backoff_sec < 0.0 {
        return Err("bold.backoff_sec must be >= 0".to_string());
    }
    if !matches!(download_format.as_str(), "json" | "tsv") {
        return Err("bold.download_format must be 'json' or 'tsv'".to_string());
    }
    if download_chunk_size == 0 {
        return Err("bold.download_chunk_size must be > 0".to_string());
    }

    Ok(BoldRuntimeConfig {
        base_url,
        timeout_sec,
        retries,
        backoff_sec,
        download_format,
        download_chunk_size,
        user_agent,
    })
}

fn toml_u64_from_value(value: Option<&TomlValue>) -> Option<u64> {
    value
        .and_then(|v| v.as_integer())
        .and_then(|v| u64::try_from(v).ok())
}

fn toml_f64_from_value(value: Option<&TomlValue>) -> Option<f64> {
    value.and_then(|v| {
        if let Some(x) = v.as_float() {
            return Some(x);
        }
        v.as_integer().map(|x| x as f64)
    })
}

fn bold_get_text(
    path: &str,
    params: &[(&str, String)],
    runtime_cfg: &BoldRuntimeConfig,
) -> Result<String, String> {
    for attempt in 0..=runtime_cfg.retries {
        let mut cmd = Command::new("curl");
        cmd.arg("-fsSL")
            .arg("--compressed")
            .arg("-H")
            .arg("Accept: application/json")
            .arg("-H")
            .arg(format!("User-Agent: {}", runtime_cfg.user_agent))
            .arg("--max-time")
            .arg(runtime_cfg.timeout_sec.to_string())
            .arg(format!("{}{}", runtime_cfg.base_url, path));
        for (k, v) in params {
            cmd.arg("--get")
                .arg("--data-urlencode")
                .arg(format!("{k}={v}"));
        }
        let out = cmd
            .output()
            .map_err(|e| format!("failed to execute curl for BOLD: {e}"))?;
        if out.status.success() {
            return String::from_utf8(out.stdout)
                .map_err(|e| format!("BOLD response decode failed: {e}"));
        }
        if attempt >= runtime_cfg.retries {
            let stderr = String::from_utf8_lossy(&out.stderr);
            return Err(format!("BOLD curl failed: {stderr}"));
        }
        if runtime_cfg.backoff_sec > 0.0 {
            thread::sleep(Duration::from_secs_f64(
                runtime_cfg.backoff_sec * (attempt + 1) as f64,
            ));
        }
    }
    Err("BOLD request failed after retries".to_string())
}

fn bold_get_json(
    path: &str,
    params: &[(&str, String)],
    runtime_cfg: &BoldRuntimeConfig,
) -> Result<JsonValue, String> {
    let text = bold_get_text(path, params, runtime_cfg)?;
    parse_json_payload(&text).map_err(|e| format!("BOLD JSON parse failed: {e}"))
}

fn parse_json_payload(text: &str) -> Result<JsonValue, String> {
    let stripped = text.trim_start_matches('\u{feff}').trim();
    if stripped.is_empty() {
        return Err("empty response".to_string());
    }

    if let Ok(value) = serde_json::from_str::<JsonValue>(stripped) {
        return Ok(value);
    }

    let stream = serde_json::Deserializer::from_str(stripped).into_iter::<JsonValue>();
    let mut values = Vec::new();
    for item in stream {
        match item {
            Ok(value) => values.push(value),
            Err(_) => {
                values.clear();
                break;
            }
        }
    }
    if !values.is_empty() {
        return if values.len() == 1 {
            Ok(values.remove(0))
        } else {
            Ok(JsonValue::Array(values))
        };
    }

    let mut line_values = Vec::new();
    for line in stripped.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        match serde_json::from_str::<JsonValue>(line) {
            Ok(value) => line_values.push(value),
            Err(err) => return Err(err.to_string()),
        }
    }
    if !line_values.is_empty() {
        return if line_values.len() == 1 {
            Ok(line_values.remove(0))
        } else {
            Ok(JsonValue::Array(line_values))
        };
    }

    Err("unable to parse JSON payload".to_string())
}

fn build_bold_taxon_query(scientific_name: &str) -> String {
    format!("tax:{scientific_name}")
}

fn extract_successful_terms(payload: &JsonValue) -> Vec<String> {
    let Some(items) = payload.get("successful_terms").and_then(|v| v.as_array()) else {
        return Vec::new();
    };
    let mut terms = Vec::new();
    for item in items {
        if let Some(s) = json_string(item) {
            terms.push(s);
            continue;
        }
        let JsonValue::Object(map) = item else {
            continue;
        };
        if let Some(s) = map.get("matched").and_then(json_string) {
            terms.push(s);
            continue;
        }
        let scope = map.get("scope").and_then(json_string);
        let field = map.get("field").and_then(json_string);
        let value = map.get("value").and_then(json_string);
        if let (Some(scope), Some(value)) = (scope, value) {
            if let Some(field) = field {
                terms.push(format!("{scope}:{field}:{value}"));
            } else {
                terms.push(format!("{scope}:{value}"));
            }
        }
    }
    terms
}

fn preprocess_bold_query(query: &str, runtime_cfg: &BoldRuntimeConfig) -> Result<String, String> {
    let payload = bold_get_json(
        "/query/preprocessor",
        &[("query", query.to_string())],
        runtime_cfg,
    )?;
    let successful_terms = extract_successful_terms(&payload);
    if successful_terms.is_empty() {
        Ok(query.to_string())
    } else {
        Ok(successful_terms.join(";"))
    }
}

fn prepare_bold_query(
    scientific_name: &str,
    runtime_cfg: &BoldRuntimeConfig,
) -> Result<PreparedBoldQuery, String> {
    let normalized_query =
        preprocess_bold_query(&build_bold_taxon_query(scientific_name), runtime_cfg)?;
    let specimen_count = fetch_bold_specimen_count(&normalized_query, runtime_cfg)?;
    if specimen_count.unwrap_or(0) > BOLD_MAX_DOCUMENT_COUNT {
        return Err(format!(
            "BOLD query exceeds maximum downloadable records ({} > {BOLD_MAX_DOCUMENT_COUNT})",
            specimen_count.unwrap_or(0)
        ));
    }
    let query_id = if specimen_count == Some(0) {
        None
    } else {
        Some(submit_bold_query(&normalized_query, runtime_cfg)?)
    };
    Ok(PreparedBoldQuery {
        normalized_query,
        specimen_count,
        query_id,
        download_format: runtime_cfg.download_format.clone(),
    })
}

fn fetch_bold_specimen_count(
    normalized_query: &str,
    runtime_cfg: &BoldRuntimeConfig,
) -> Result<Option<u64>, String> {
    let payload = bold_get_json(
        "/summary",
        &[
            ("query", normalized_query.to_string()),
            ("fields", "specimens".to_string()),
            ("reduce_operation", "count".to_string()),
        ],
        runtime_cfg,
    )?;
    let count_value = payload
        .get("counts")
        .and_then(|v| v.get("specimens"))
        .or_else(|| payload.get("specimens"))
        .or_else(|| payload.get("count"));
    Ok(count_value
        .and_then(|value| json_string(value))
        .and_then(|value| value.parse::<u64>().ok()))
}

fn submit_bold_query(
    normalized_query: &str,
    runtime_cfg: &BoldRuntimeConfig,
) -> Result<String, String> {
    let payload = bold_get_json(
        "/query",
        &[
            ("query", normalized_query.to_string()),
            ("extent", "full".to_string()),
        ],
        runtime_cfg,
    )?;
    first_scalar(&payload, &["query_id", "queryId", "id", "token"])
        .ok_or_else(|| "BOLD query response did not contain a query_id".to_string())
}

fn url_encode_path_segment(value: &str) -> String {
    let mut out = String::new();
    for byte in value.as_bytes() {
        let ch = *byte as char;
        if ch.is_ascii_alphanumeric() || matches!(ch, '-' | '_' | '.' | '~') {
            out.push(ch);
        } else {
            out.push_str(&format!("%{:02X}", byte));
        }
    }
    out
}

fn download_bold_documents_to_path(
    query_id: &str,
    runtime_cfg: &BoldRuntimeConfig,
    dest_path: &Path,
    format_name: Option<&str>,
    scientific_name: &str,
    log_path: &Path,
) -> Result<BoldDownloadMeta, String> {
    let download_format = format_name
        .unwrap_or(runtime_cfg.download_format.as_str())
        .trim()
        .to_ascii_lowercase();
    if !matches!(download_format.as_str(), "json" | "tsv") {
        return Err(format!(
            "Unsupported BOLD download format: {download_format}"
        ));
    }

    let encoded_query_id = url_encode_path_segment(query_id);
    let part_path = dest_path.with_extension(format!(
        "{}.part",
        dest_path
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or_default()
    ));
    if let Some(parent) = dest_path.parent() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("failed to create {}: {e}", parent.display()))?;
    }

    for attempt in 0..=runtime_cfg.retries {
        let mut cmd = Command::new("curl");
        cmd.arg("-fsSL")
            .arg("--compressed")
            .arg("-H")
            .arg(if download_format == "json" {
                "Accept: application/json"
            } else {
                "Accept: text/tab-separated-values"
            })
            .arg("-H")
            .arg(format!("User-Agent: {}", runtime_cfg.user_agent))
            .arg("--max-time")
            .arg(runtime_cfg.timeout_sec.to_string())
            .arg(format!(
                "{}/documents/{encoded_query_id}/download",
                runtime_cfg.base_url
            ))
            .arg("--get")
            .arg("--data-urlencode")
            .arg(format!("format={download_format}"))
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());

        let mut child = cmd
            .spawn()
            .map_err(|e| format!("failed to execute curl for BOLD download: {e}"))?;
        let mut stdout = child
            .stdout
            .take()
            .ok_or_else(|| "failed to capture BOLD download stdout".to_string())?;
        let mut stderr = child
            .stderr
            .take()
            .ok_or_else(|| "failed to capture BOLD download stderr".to_string())?;
        let mut downloaded_bytes = 0u64;
        let mut last_reported = 0u64;
        let mut out_f = File::create(&part_path)
            .map_err(|e| format!("failed to create {}: {e}", part_path.display()))?;
        let mut buffer = vec![0u8; runtime_cfg.download_chunk_size];

        loop {
            let read_len = stdout
                .read(&mut buffer)
                .map_err(|e| format!("BOLD download read failed: {e}"))?;
            if read_len == 0 {
                break;
            }
            out_f
                .write_all(&buffer[..read_len])
                .map_err(|e| format!("failed to write {}: {e}", part_path.display()))?;
            downloaded_bytes += read_len as u64;
            if downloaded_bytes == read_len as u64
                || downloaded_bytes.saturating_sub(last_reported) >= 1024 * 1024
            {
                append_log_line(
                    log_path,
                    &format!(
                        "# bold progress: taxon={scientific_name} phase=download bytes={downloaded_bytes}"
                    ),
                )?;
                last_reported = downloaded_bytes;
            }
        }

        let mut stderr_text = String::new();
        let _ = stderr.read_to_string(&mut stderr_text);
        let status = child
            .wait()
            .map_err(|e| format!("failed to wait for BOLD download: {e}"))?;
        if status.success() {
            fs::rename(&part_path, dest_path)
                .map_err(|e| format!("failed to finalize {}: {e}", dest_path.display()))?;
            return Ok(BoldDownloadMeta { downloaded_bytes });
        }

        let _ = fs::remove_file(&part_path);
        if attempt >= runtime_cfg.retries {
            return Err(format!("BOLD curl failed: {stderr_text}"));
        }
        if runtime_cfg.backoff_sec > 0.0 {
            thread::sleep(Duration::from_secs_f64(
                runtime_cfg.backoff_sec * (attempt + 1) as f64,
            ));
        }
    }

    Err("BOLD document download failed after retries".to_string())
}

fn extract_bold_document_rows(payload: &JsonValue) -> Vec<JsonValue> {
    match payload {
        JsonValue::Array(items) => items
            .iter()
            .filter(|item| item.is_object())
            .cloned()
            .collect(),
        JsonValue::Object(map) => {
            for key in ["documents", "records", "items", "data", "results"] {
                if let Some(items) = map.get(key).and_then(|v| v.as_array()) {
                    return items
                        .iter()
                        .filter(|item| item.is_object())
                        .cloned()
                        .collect();
                }
            }
            if first_scalar(
                payload,
                &["processid", "sampleid", "marker_code", "nucleotides", "nuc"],
            )
            .is_some()
            {
                vec![payload.clone()]
            } else {
                Vec::new()
            }
        }
        _ => Vec::new(),
    }
}

fn tsv_cells(line: &str) -> Vec<String> {
    line.trim_end_matches(['\r', '\n'])
        .split('\t')
        .map(|cell| cell.to_string())
        .collect()
}

fn for_each_tsv_document_row<F>(path: &Path, mut on_row: F) -> Result<(), String>
where
    F: FnMut(JsonValue) -> Result<(), String>,
{
    let file = File::open(path).map_err(|e| format!("failed to open {}: {e}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut header_line = String::new();
    if reader
        .read_line(&mut header_line)
        .map_err(|e| format!("failed to read {}: {e}", path.display()))?
        == 0
    {
        return Ok(());
    }
    let mut headers = tsv_cells(&header_line);
    if let Some(first) = headers.first_mut() {
        *first = first.trim_start_matches('\u{feff}').trim().to_string();
    }
    for header in &mut headers {
        *header = header.trim().to_string();
    }

    let mut line = String::new();
    loop {
        line.clear();
        if reader
            .read_line(&mut line)
            .map_err(|e| format!("failed to read {}: {e}", path.display()))?
            == 0
        {
            break;
        }
        if line.trim().is_empty() {
            continue;
        }
        let values = tsv_cells(&line);
        let mut row = serde_json::Map::new();
        let mut has_value = false;
        for (idx, header) in headers.iter().enumerate() {
            if header.is_empty() {
                continue;
            }
            let value = values.get(idx).cloned().unwrap_or_default();
            if !value.trim().is_empty() {
                has_value = true;
            }
            row.insert(header.clone(), JsonValue::String(value));
        }
        if has_value {
            on_row(JsonValue::Object(row))?;
        }
    }
    Ok(())
}

fn for_each_json_document_row<F>(path: &Path, mut on_row: F) -> Result<(), String>
where
    F: FnMut(JsonValue) -> Result<(), String>,
{
    let text =
        fs::read_to_string(path).map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    let payload = parse_json_payload(&text)?;
    for row in extract_bold_document_rows(&payload) {
        on_row(row)?;
    }
    Ok(())
}

fn for_each_document_row_from_path<F>(path: &Path, fmt: &str, on_row: F) -> Result<(), String>
where
    F: FnMut(JsonValue) -> Result<(), String>,
{
    match fmt.trim().to_ascii_lowercase().as_str() {
        "tsv" => for_each_tsv_document_row(path, on_row),
        "json" => for_each_json_document_row(path, on_row),
        other => Err(format!("Unsupported BOLD row iteration format: {other}")),
    }
}

fn parse_accession_tokens(raw_value: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    for token in raw_value.split(|ch: char| matches!(ch, ',' | ' ' | ';' | '|' | '/')) {
        let value = token.trim();
        if !value.is_empty() && !tokens.iter().any(|v| v == value) {
            tokens.push(value.to_string());
        }
    }
    tokens
}

fn normalize_marker_text(value: &str) -> String {
    value
        .chars()
        .filter(|ch| ch.is_ascii_alphanumeric())
        .map(|ch| ch.to_ascii_lowercase())
        .collect()
}

fn marker_matches(marker_text: &str, candidate: &str, exact_only: bool) -> bool {
    let left = normalize_marker_text(marker_text);
    let right = normalize_marker_text(candidate);
    if left.is_empty() || right.is_empty() {
        return false;
    }
    if left == right {
        return true;
    }
    if exact_only {
        return false;
    }
    left.starts_with(&right)
        || right.starts_with(&left)
        || left.contains(&right)
        || right.contains(&left)
}

fn bold_marker_match(
    marker_text: &str,
    marker_keys: &[String],
    markers_tbl: &toml::map::Map<String, TomlValue>,
) -> Option<(String, String)> {
    for strategy in ["marker_codes", "aliases", "phrases", "key"] {
        for marker_key in marker_keys {
            let marker_cfg = markers_tbl.get(marker_key)?.as_table()?;
            let (candidates, exact_only): (Vec<String>, bool) = match strategy {
                "marker_codes" => (
                    marker_cfg
                        .get("bold")
                        .and_then(|v| v.as_table())
                        .and_then(|v| v.get("marker_codes"))
                        .and_then(|v| v.as_array())
                        .map(|arr| {
                            arr.iter()
                                .filter_map(|v| v.as_str())
                                .map(|s| s.trim().to_string())
                                .filter(|s| !s.is_empty())
                                .collect()
                        })
                        .unwrap_or_default(),
                    true,
                ),
                "aliases" => (
                    marker_cfg
                        .get("aliases")
                        .and_then(|v| v.as_array())
                        .map(|arr| {
                            arr.iter()
                                .filter_map(|v| v.as_str())
                                .map(|s| s.trim().to_string())
                                .filter(|s| !s.is_empty())
                                .collect()
                        })
                        .unwrap_or_default(),
                    false,
                ),
                "phrases" => (
                    marker_cfg
                        .get("phrases")
                        .and_then(|v| v.as_array())
                        .map(|arr| {
                            arr.iter()
                                .filter_map(|v| v.as_str())
                                .map(|s| s.trim().to_string())
                                .filter(|s| !s.is_empty())
                                .collect()
                        })
                        .unwrap_or_default(),
                    false,
                ),
                _ => (vec![marker_key.clone()], false),
            };
            for candidate in candidates {
                if marker_matches(marker_text, &candidate, exact_only) {
                    return Some((marker_key.clone(), marker_text.to_string()));
                }
            }
        }
    }
    None
}

fn sanitize_sequence(sequence: &str) -> String {
    sequence
        .chars()
        .filter(|ch| ch.is_ascii_alphabetic())
        .map(|ch| ch.to_ascii_uppercase())
        .collect()
}

fn normalize_bold_row(
    raw_row: &JsonValue,
    marker_keys: &[String],
    markers_tbl: &toml::map::Map<String, TomlValue>,
    output_tbl: Option<&toml::map::Map<String, TomlValue>>,
) -> Option<BoldRecord> {
    let marker_text = first_scalar(raw_row, &["marker_code", "marker", "markercode"])?;
    let (marker_key, marker_label) = bold_marker_match(&marker_text, marker_keys, markers_tbl)?;
    let sequence = sanitize_sequence(&first_scalar(
        raw_row,
        &[
            "nucleotides",
            "nuc",
            "sequence",
            "nucleotide",
            "nucleotidesequence",
        ],
    )?);
    if sequence.is_empty() {
        return None;
    }

    let processid = first_scalar(raw_row, &["processid", "process_id"]).unwrap_or_default();
    let sampleid = first_scalar(raw_row, &["sampleid", "sample_id"]).unwrap_or_default();
    let accession =
        first_scalar(raw_row, &["insdcacs", "insdc_acs", "genbank_accession"]).unwrap_or_default();
    let taxon_name = first_scalar(
        raw_row,
        &[
            "species",
            "scientific_name",
            "taxon_name",
            "identification",
            "taxonomy_species",
        ],
    )
    .unwrap_or_else(|| "unknown".to_string());
    let source_record_id = if !processid.is_empty() {
        processid.clone()
    } else if !sampleid.is_empty() {
        sampleid.clone()
    } else if !accession.is_empty() {
        accession.clone()
    } else {
        format!(
            "bold_{}",
            sanitize_header(&format!("{}_{}_{}", taxon_name, marker_key, sequence.len()))
        )
    };
    let header_format = markers_tbl
        .get(&marker_key)
        .and_then(|v| v.as_table())
        .map(|cfg| {
            if cfg.get("header_format").is_none() {
                DEFAULT_BOLD_HEADER_FORMAT.to_string()
            } else {
                resolve_header_format(cfg, output_tbl)
            }
        })
        .unwrap_or_else(|| DEFAULT_BOLD_HEADER_FORMAT.to_string());

    Some(BoldRecord {
        source_record_id,
        processid,
        sampleid,
        accession,
        taxon_name,
        marker_key,
        marker_label,
        sequence,
        header_format,
    })
}

fn fetch_bold_records_for_taxon_with_logging(
    scientific_name: &str,
    marker_keys: &[String],
    markers_tbl: &toml::map::Map<String, TomlValue>,
    output_tbl: Option<&toml::map::Map<String, TomlValue>>,
    runtime_cfg: &BoldRuntimeConfig,
    log_path: &Path,
) -> Result<(Vec<BoldRecord>, BoldFetchStats), String> {
    append_log_line(
        log_path,
        &format!("# bold progress: taxon={scientific_name} phase=preprocess"),
    )?;
    let prepared = prepare_bold_query(scientific_name, runtime_cfg)?;
    let specimen_count = prepared.specimen_count.unwrap_or(0);
    append_log_line(
        log_path,
        &format!(
            "# bold progress: taxon={scientific_name} phase=summary specimens={}",
            specimen_count
        ),
    )?;
    if prepared.specimen_count == Some(0) {
        return Ok((
            Vec::new(),
            BoldFetchStats {
                normalized_query: prepared.normalized_query,
                specimen_count: 0,
                download_format: prepared.download_format,
                downloaded_bytes: 0,
                downloaded_rows: 0,
                matched_rows: 0,
            },
        ));
    }

    let query_id = prepared
        .query_id
        .as_ref()
        .ok_or_else(|| "BOLD query preparation did not return a query_id".to_string())?;
    append_log_line(
        log_path,
        &format!("# bold progress: taxon={scientific_name} phase=query"),
    )?;
    append_log_line(
        log_path,
        &format!("# bold progress: taxon={scientific_name} phase=download"),
    )?;
    let query_stub = sanitize_header(query_id);
    let query_stub = if query_stub.len() > 24 {
        &query_stub[..24]
    } else {
        &query_stub
    };
    let download_ext = if prepared.download_format == "json" {
        "json"
    } else {
        "tsv"
    };
    let download_path =
        std::env::temp_dir().join(format!("taxondbbuilder_gui_{query_stub}.{download_ext}"));
    let download_meta = download_bold_documents_to_path(
        query_id,
        runtime_cfg,
        &download_path,
        Some(prepared.download_format.as_str()),
        scientific_name,
        log_path,
    )?;

    let mut downloaded_rows = 0u64;
    let mut normalized_rows = Vec::new();
    for_each_document_row_from_path(&download_path, prepared.download_format.as_str(), |row| {
        downloaded_rows += 1;
        if let Some(record) = normalize_bold_row(&row, marker_keys, markers_tbl, output_tbl) {
            normalized_rows.push(record);
        }
        Ok(())
    })?;
    let matched_rows = normalized_rows.len() as u64;
    let _ = fs::remove_file(&download_path);
    append_log_line(
        log_path,
        &format!(
            "# bold progress: taxon={scientific_name} phase=filter downloaded={downloaded_rows} matched={matched_rows}"
        ),
    )?;

    Ok((
        normalized_rows,
        BoldFetchStats {
            normalized_query: prepared.normalized_query,
            specimen_count,
            download_format: prepared.download_format,
            downloaded_bytes: download_meta.downloaded_bytes,
            downloaded_rows,
            matched_rows,
        },
    ))
}

fn fetch_taxonomy_scientific_name(taxid: &str) -> Result<String, String> {
    let json = eutils_get_json(
        "esummary.fcgi",
        &[
            ("db", "taxonomy".to_string()),
            ("id", taxid.to_string()),
            ("retmode", "json".to_string()),
        ],
    )?;
    json.get("result")
        .and_then(|v| v.get(taxid))
        .and_then(|v| v.get("scientificname"))
        .and_then(|v| v.as_str())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .ok_or_else(|| format!("failed to resolve scientific name for taxid {taxid}"))
}

fn resolve_taxon(taxon: &str) -> Result<ResolvedTaxon, String> {
    if taxon.chars().all(|c| c.is_ascii_digit()) {
        let scientific_name = fetch_taxonomy_scientific_name(taxon)?;
        return Ok(ResolvedTaxon {
            taxid: taxon.to_string(),
            scientific_name,
        });
    }
    let taxid = resolve_taxid(taxon)?;
    let scientific_name =
        fetch_taxonomy_scientific_name(&taxid).unwrap_or_else(|_| taxon.to_string());
    Ok(ResolvedTaxon {
        taxid,
        scientific_name,
    })
}

fn append_log_line(log_path: &Path, line: &str) -> Result<(), String> {
    append_log_with_level(log_path, "INFO", line)
}

fn append_error_log_line(log_path: &Path, line: &str) -> Result<(), String> {
    append_log_with_level(log_path, "ERROR", line)
}

fn append_log_with_level(log_path: &Path, level: &str, line: &str) -> Result<(), String> {
    let mut f = OpenOptions::new()
        .create(true)
        .append(true)
        .open(log_path)
        .map_err(|e| format!("failed to open log {}: {e}", log_path.display()))?;
    writeln!(f, "{} {:<5} {line}", log_timestamp_now(), level)
        .map_err(|e| format!("failed to write log {}: {e}", log_path.display()))
}

fn log_timestamp_now() -> String {
    Local::now().format("%Y-%m-%dT%H:%M:%S%.3f%:z").to_string()
}

fn append_run_header(
    log_path: &Path,
    params: &BuildParams,
    source: &str,
    resolved_taxa: &[ResolvedTaxon],
    marker_keys: &[String],
    started_at: &str,
) -> Result<(), String> {
    append_log_line(log_path, "# runner: rust-runner")?;
    append_log_line(log_path, &format!("# started_at: {started_at}"))?;
    append_log_line(log_path, &format!("# log_path: {}", log_path.display()))?;
    append_log_line(
        log_path,
        &format!("# output_file: {}", params.output_file.display()),
    )?;
    append_log_line(
        log_path,
        &format!("# dump_gb_dir: {}", params.dump_gb_dir.display()),
    )?;
    append_log_line(log_path, &format!("# source: {source}"))?;
    append_log_line(
        log_path,
        &format!("# config: {}", params.config_path.display()),
    )?;
    append_log_line(log_path, &format!("# resume: {}", params.resume))?;
    append_log_line(
        log_path,
        &format!(
            "# taxids: {:?}",
            resolved_taxa
                .iter()
                .map(|item| item.taxid.clone())
                .collect::<Vec<_>>()
        ),
    )?;
    append_log_line(
        log_path,
        &format!(
            "# scientific_names: {:?}",
            resolved_taxa
                .iter()
                .map(|item| item.scientific_name.clone())
                .collect::<Vec<_>>()
        ),
    )?;
    append_log_line(log_path, &format!("# markers: {:?}", marker_keys))
}

fn toml_string(tbl: &toml::map::Map<String, TomlValue>, key: &str, default: &str) -> String {
    tbl.get(key)
        .and_then(|v| v.as_str())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| default.to_string())
}

fn toml_bool(tbl: &toml::map::Map<String, TomlValue>, key: &str, default: bool) -> bool {
    tbl.get(key).and_then(|v| v.as_bool()).unwrap_or(default)
}

fn toml_u64(tbl: &toml::map::Map<String, TomlValue>, key: &str, default: u64) -> u64 {
    tbl.get(key)
        .and_then(|v| v.as_integer())
        .and_then(|x| u64::try_from(x).ok())
        .unwrap_or(default)
}

fn toml_f64(tbl: &toml::map::Map<String, TomlValue>, key: &str) -> Option<f64> {
    if let Some(x) = tbl.get(key).and_then(|v| v.as_float()) {
        return Some(x);
    }
    tbl.get(key).and_then(|v| v.as_integer()).map(|x| x as f64)
}

fn parse_toml(path: &Path) -> Result<TomlValue, String> {
    let text =
        fs::read_to_string(path).map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    text.parse::<TomlValue>()
        .map_err(|e| format!("failed to parse {}: {e}", path.display()))
}

fn resolve_support_file_path(raw_path: &str, config_path: &Path) -> Result<PathBuf, String> {
    let path = PathBuf::from(raw_path);
    let mut candidates = Vec::new();
    if path.is_absolute() {
        candidates.push(path);
    } else {
        if let Some(parent) = config_path.parent() {
            candidates.push(parent.join(&path));
        }
        if let Ok(cwd) = std::env::current_dir() {
            candidates.push(cwd.join(&path));
        }
    }
    for c in &candidates {
        if c.exists() && c.is_file() {
            return Ok(c.clone());
        }
    }
    let tried = candidates
        .iter()
        .map(|p| p.to_string_lossy().to_string())
        .collect::<Vec<_>>()
        .join(", ");
    Err(format!("markers file not found. tried: {tried}"))
}

fn load_markers_table(
    cfg_tbl: &toml::map::Map<String, TomlValue>,
    config_path: &Path,
) -> Result<toml::map::Map<String, TomlValue>, String> {
    let markers_section = cfg_tbl
        .get("markers")
        .and_then(|v| v.as_table())
        .ok_or_else(|| "missing [markers] section".to_string())?;

    let mut merged: toml::map::Map<String, TomlValue> = toml::map::Map::new();

    if let Some(file_raw) = markers_section.get("file").and_then(|v| v.as_str()) {
        let file_raw = file_raw.trim();
        if !file_raw.is_empty() {
            let markers_path = resolve_support_file_path(file_raw, config_path)?;
            let markers_file_toml = parse_toml(&markers_path)?;
            let from_file = markers_file_toml
                .get("markers")
                .and_then(|v| v.as_table())
                .ok_or_else(|| {
                    format!(
                        "markers file {} must define [markers] table",
                        markers_path.display()
                    )
                })?;
            for (k, v) in from_file {
                merged.insert(k.clone(), v.clone());
            }
        }
    }

    for (k, v) in markers_section {
        if k == "file" || k == "markers_file" {
            continue;
        }
        merged.insert(k.clone(), v.clone());
    }

    if merged.is_empty() {
        return Err("no markers were resolved from [markers]".to_string());
    }
    Ok(merged)
}

fn sanitize_header(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut prev_space = false;
    for ch in text.trim().chars() {
        if ch.is_whitespace() {
            if !prev_space {
                out.push('_');
            }
            prev_space = true;
            continue;
        }
        prev_space = false;
        if ch.is_ascii_alphanumeric() || matches!(ch, '.' | '_' | '-') {
            out.push(ch);
        } else {
            out.push('_');
        }
    }
    out
}

fn build_header(template: &str, values: &HashMap<&str, String>) -> String {
    let mut out = template.to_string();
    for (k, v) in values {
        let needle = format!("{{{k}}}");
        out = out.replace(&needle, v);
    }
    out
}

fn marker_aliases(marker_cfg: &toml::map::Map<String, TomlValue>) -> Vec<String> {
    marker_cfg
        .get("aliases")
        .and_then(|v| v.as_array())
        .map(|arr| {
            arr.iter()
                .filter_map(|v| v.as_str())
                .map(|s| s.trim().to_ascii_lowercase())
                .filter(|s| !s.is_empty())
                .collect::<Vec<_>>()
        })
        .unwrap_or_default()
}

fn resolve_marker_key(
    input: &str,
    markers_tbl: &toml::map::Map<String, TomlValue>,
) -> Result<String, String> {
    let value = input.trim().to_ascii_lowercase();
    let mut exact = Vec::new();
    let mut prefix = Vec::new();
    for (key, marker_v) in markers_tbl {
        let Some(marker_cfg) = marker_v.as_table() else {
            continue;
        };
        let key_l = key.to_ascii_lowercase();
        let aliases = marker_aliases(marker_cfg);
        if value == key_l || aliases.iter().any(|a| a == &value) {
            exact.push(key.clone());
            continue;
        }
        if key_l.starts_with(&value) || aliases.iter().any(|a| a.starts_with(&value)) {
            prefix.push(key.clone());
        }
    }
    if exact.len() == 1 {
        return Ok(exact.remove(0));
    }
    if exact.len() > 1 {
        return Err(format!(
            "marker '{input}' matches multiple entries: {}",
            exact.join(", ")
        ));
    }
    if prefix.len() == 1 {
        return Ok(prefix.remove(0));
    }
    if prefix.len() > 1 {
        return Err(format!(
            "marker '{input}' matches multiple entries: {}",
            prefix.join(", ")
        ));
    }
    Err(format!("marker '{input}' not found in config"))
}

fn marker_query_terms(marker_cfg: &toml::map::Map<String, TomlValue>) -> Vec<String> {
    let mut out = Vec::new();
    if let Some(terms) = marker_cfg.get("terms").and_then(|v| v.as_array()) {
        for t in terms {
            if let Some(s) = t.as_str() {
                let s = s.trim();
                if !s.is_empty() {
                    out.push(s.to_string());
                }
            }
        }
    }
    if let Some(phrases) = marker_cfg.get("phrases").and_then(|v| v.as_array()) {
        for p in phrases {
            if let Some(s) = p.as_str() {
                let s = s.trim();
                if s.is_empty() {
                    continue;
                }
                if s.contains('[') && s.contains(']') {
                    out.push(s.to_string());
                } else {
                    out.push(format!("\"{}\"[All Fields]", s.replace('"', "\\\"")));
                }
            }
        }
    }
    out
}

fn region_patterns(marker_cfg: &toml::map::Map<String, TomlValue>) -> Vec<String> {
    let mut out = Vec::new();
    if let Some(items) = marker_cfg.get("region_patterns").and_then(|v| v.as_array()) {
        for p in items {
            if let Some(s) = p.as_str() {
                let s = s.trim();
                if !s.is_empty() {
                    out.push(s.to_string());
                }
            }
        }
    }
    if !out.is_empty() {
        return out;
    }
    for term in marker_query_terms(marker_cfg) {
        let stripped = Regex::new(r"\s*\[[^\]]+\]")
            .expect("regex")
            .replace_all(&term, "")
            .replace('"', "");
        let stripped = stripped.trim();
        if !stripped.is_empty() {
            out.push(regex::escape(stripped));
        }
    }
    out
}

fn resolve_header_format(
    marker_cfg: &toml::map::Map<String, TomlValue>,
    output_cfg: Option<&toml::map::Map<String, TomlValue>>,
) -> String {
    let default_fmt = output_cfg
        .and_then(|o| o.get("default_header_format"))
        .and_then(|v| v.as_str())
        .unwrap_or(DEFAULT_HEADER_FORMAT)
        .to_string();
    let Some(key) = marker_cfg.get("header_format").and_then(|v| v.as_str()) else {
        return default_fmt;
    };
    let key = key.trim();
    if key.is_empty() {
        return default_fmt;
    }
    if let Some(hf_tbl) = output_cfg
        .and_then(|o| o.get("header_formats"))
        .and_then(|v| v.as_table())
    {
        if let Some(v) = hf_tbl.get(key).and_then(|v| v.as_str()) {
            return v.to_string();
        }
    }
    key.to_string()
}

fn build_filter_terms(
    filters: Option<&toml::map::Map<String, TomlValue>>,
) -> Result<Vec<String>, String> {
    let Some(filters) = filters else {
        return Ok(Vec::new());
    };
    let mut terms = Vec::new();

    let list_from_value = |name: &str| -> Result<Vec<String>, String> {
        let Some(v) = filters.get(name) else {
            return Ok(Vec::new());
        };
        if let Some(s) = v.as_str() {
            return Ok(vec![s.to_string()]);
        }
        if let Some(arr) = v.as_array() {
            let mut out = Vec::new();
            for item in arr {
                let Some(s) = item.as_str() else {
                    return Err(format!("filters.{name} must contain strings"));
                };
                out.push(s.to_string());
            }
            return Ok(out);
        }
        Err(format!("filters.{name} must be string or string array"))
    };

    for s in list_from_value("filter")? {
        terms.push(format!("{s}[filter]"));
    }
    for s in list_from_value("properties")? {
        terms.push(format!("{s}[prop]"));
    }

    let min_len = filters
        .get("sequence_length_min")
        .and_then(|v| v.as_integer());
    let max_len = filters
        .get("sequence_length_max")
        .and_then(|v| v.as_integer());
    if min_len.is_some() || max_len.is_some() {
        let min_v = min_len.unwrap_or(0);
        let max_v = max_len.unwrap_or(1_000_000_000);
        terms.push(format!("{min_v}[SLEN] : {max_v}[SLEN]"));
    }

    if let Some(raw) = filters.get("raw") {
        if let Some(s) = raw.as_str() {
            terms.push(s.to_string());
        } else if let Some(arr) = raw.as_array() {
            for item in arr {
                let Some(s) = item.as_str() else {
                    return Err("filters.raw must contain strings".to_string());
                };
                terms.push(s.to_string());
            }
        } else {
            return Err("filters.raw must be string or string array".to_string());
        }
    }
    Ok(terms)
}

fn build_query(
    taxid: &str,
    marker_query: &str,
    filter_terms: &[String],
    taxon_noexp: bool,
) -> String {
    let tax_term = if taxon_noexp {
        format!("txid{taxid}[Organism:noexp]")
    } else {
        format!("txid{taxid}[Organism]")
    };
    let mut parts = vec![tax_term, marker_query.to_string()];
    parts.extend(filter_terms.iter().cloned());
    parts
        .into_iter()
        .map(|p| format!("({p})"))
        .collect::<Vec<_>>()
        .join(" AND ")
}

fn sanitize_json_control_chars(text: &str) -> String {
    text.chars()
        .filter(|ch| !matches!(ch, '\u{0000}'..='\u{001F}'))
        .collect()
}

fn eutils_get_json(endpoint: &str, params: &[(&str, String)]) -> Result<JsonValue, String> {
    let text = eutils_get_text(endpoint, params)?;
    if let Ok(value) = serde_json::from_str::<JsonValue>(&text) {
        return Ok(value);
    }
    let sanitized = sanitize_json_control_chars(&text);
    serde_json::from_str::<JsonValue>(&sanitized)
        .map_err(|e| format!("eutils json parse failed: {e}"))
}

fn curl_should_retry(exit_code: Option<i32>, stderr: &str) -> bool {
    matches!(exit_code, Some(18 | 28 | 52 | 55 | 56 | 92))
        || stderr.contains("Transferred a partial file")
        || stderr.contains("Operation timed out")
        || stderr.contains("Empty reply from server")
        || stderr.contains("Failure when receiving data from the peer")
        || stderr.contains("Send failure")
        || stderr.contains("HTTP/2 stream")
}

fn eutils_request_text(
    endpoint: &str,
    params: &[(&str, String)],
    use_get: bool,
) -> Result<String, String> {
    let url = format!("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}");
    let retries = 3u64;
    let backoff_sec = 1.5f64;

    for attempt in 0..=retries {
        let mut cmd = Command::new("curl");
        cmd.arg("-fsSL").arg("--http1.1").arg(url.clone());
        for (k, v) in params {
            if use_get {
                cmd.arg("--get");
            }
            cmd.arg("--data-urlencode").arg(format!("{k}={v}"));
        }
        let out = cmd
            .output()
            .map_err(|e| format!("failed to execute curl: {e}"))?;
        if out.status.success() {
            return String::from_utf8(out.stdout)
                .map_err(|e| format!("curl output decode failed: {e}"));
        }

        let stderr = String::from_utf8_lossy(&out.stderr).trim().to_string();
        if attempt >= retries || !curl_should_retry(out.status.code(), &stderr) {
            return Err(format!("curl failed: {stderr}"));
        }
        thread::sleep(Duration::from_secs_f64(backoff_sec * (attempt + 1) as f64));
    }

    Err("curl failed after retries".to_string())
}

fn eutils_get_validated_json<F>(
    endpoint: &str,
    params: &[(&str, String)],
    validate: F,
    invalid_payload_message: &str,
) -> Result<JsonValue, String>
where
    F: Fn(&JsonValue) -> bool,
{
    let retries = 3u64;
    let backoff_sec = 1.5f64;
    let mut last_error: Option<String> = None;

    for attempt in 0..=retries {
        match eutils_get_json(endpoint, params) {
            Ok(json) if validate(&json) => return Ok(json),
            Ok(_) => last_error = Some(invalid_payload_message.to_string()),
            Err(err) => last_error = Some(err),
        }
        if attempt < retries {
            thread::sleep(Duration::from_secs_f64(backoff_sec * (attempt + 1) as f64));
        }
    }

    Err(last_error.unwrap_or_else(|| invalid_payload_message.to_string()))
}

fn eutils_get_text(endpoint: &str, params: &[(&str, String)]) -> Result<String, String> {
    eutils_request_text(endpoint, params, true)
}

fn eutils_post_text(endpoint: &str, params: &[(&str, String)]) -> Result<String, String> {
    eutils_request_text(endpoint, params, false)
}

fn resolve_taxid(taxon: &str) -> Result<String, String> {
    if taxon.chars().all(|c| c.is_ascii_digit()) {
        return Ok(taxon.to_string());
    }
    let term = format!("\"{taxon}\"[Scientific Name]");
    let json = eutils_get_validated_json(
        "esearch.fcgi",
        &[
            ("db", "taxonomy".to_string()),
            ("term", term),
            ("retmax", "5".to_string()),
            ("retmode", "json".to_string()),
        ],
        |json| {
            json.get("esearchresult")
                .and_then(|v| v.get("idlist"))
                .and_then(|v| v.as_array())
                .is_some()
        },
        "taxonomy search returned invalid payload",
    )?;
    let ids = json
        .get("esearchresult")
        .and_then(|v| v.get("idlist"))
        .and_then(|v| v.as_array())
        .ok_or_else(|| "taxonomy search returned invalid payload".to_string())?;
    let first = ids
        .first()
        .and_then(|v| v.as_str())
        .ok_or_else(|| format!("taxon not found in NCBI taxonomy: {taxon}"))?;
    Ok(first.to_string())
}

fn split_genbank_records(chunk: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut current = String::new();
    for line in chunk.lines() {
        if line.trim() == "//" {
            if !current.trim().is_empty() {
                out.push(current.clone());
            }
            current.clear();
            continue;
        }
        current.push_str(line);
        current.push('\n');
    }
    if !current.trim().is_empty() {
        out.push(current);
    }
    out
}

fn parse_location(raw: &str) -> Option<(usize, usize, i8)> {
    let s = raw.trim();
    let comp_re = Regex::new(r"^complement\((<?\d+)\.\.(>?\d+)\)$").expect("regex");
    let simple_re = Regex::new(r"^(<?\d+)\.\.(>?\d+)$").expect("regex");
    if let Some(c) = comp_re.captures(s) {
        let a = c
            .get(1)?
            .as_str()
            .trim_start_matches('<')
            .parse::<usize>()
            .ok()?;
        let b = c
            .get(2)?
            .as_str()
            .trim_start_matches('>')
            .parse::<usize>()
            .ok()?;
        if a == 0 || b == 0 || a > b {
            return None;
        }
        return Some((a, b, -1));
    }
    if let Some(c) = simple_re.captures(s) {
        let a = c
            .get(1)?
            .as_str()
            .trim_start_matches('<')
            .parse::<usize>()
            .ok()?;
        let b = c
            .get(2)?
            .as_str()
            .trim_start_matches('>')
            .parse::<usize>()
            .ok()?;
        if a == 0 || b == 0 || a > b {
            return None;
        }
        return Some((a, b, 1));
    }
    None
}

fn dna_complement(ch: char) -> char {
    match ch.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        'R' => 'Y',
        'Y' => 'R',
        'S' => 'S',
        'W' => 'W',
        'K' => 'M',
        'M' => 'K',
        'B' => 'V',
        'D' => 'H',
        'H' => 'D',
        'V' => 'B',
        _ => 'N',
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(dna_complement).collect()
}

fn parse_genbank_record(raw: &str) -> Option<ParsedRecord> {
    let mut accession = String::new();
    let mut organism = String::from("unknown");
    let mut in_features = false;
    let mut in_origin = false;
    let mut sequence = String::new();
    let mut features = Vec::new();
    let mut current_feature: Option<ParsedFeature> = None;
    let mut current_qualifier_key: Option<String> = None;

    for line in raw.lines() {
        if let Some(rest) = line.strip_prefix("ACCESSION") {
            let acc = rest.trim().split_whitespace().next().unwrap_or("");
            if !acc.is_empty() {
                accession = acc.to_string();
            }
            continue;
        }
        if let Some(rest) = line.strip_prefix("  ORGANISM  ") {
            let s = rest.trim();
            if !s.is_empty() {
                organism = s.to_string();
            }
            continue;
        }
        if line.starts_with("FEATURES             Location/Qualifiers") {
            in_features = true;
            in_origin = false;
            continue;
        }
        if line.starts_with("ORIGIN") {
            if let Some(f) = current_feature.take() {
                features.push(f);
            }
            in_features = false;
            in_origin = true;
            continue;
        }

        if in_origin {
            for ch in line.chars() {
                if ch.is_ascii_alphabetic() {
                    sequence.push(ch.to_ascii_uppercase());
                }
            }
            continue;
        }

        if in_features {
            if line.starts_with("     ") && !line.starts_with("                     /") {
                if let Some(f) = current_feature.take() {
                    features.push(f);
                }
                let key = line.get(5..21).unwrap_or("").trim().to_string();
                let loc = line.get(21..).unwrap_or("").trim().to_string();
                if !key.is_empty() && !loc.is_empty() {
                    current_feature = Some(ParsedFeature {
                        feature_type: key,
                        location_raw: loc,
                        qualifiers: HashMap::new(),
                    });
                }
                current_qualifier_key = None;
                continue;
            }

            if line.starts_with("                     /") {
                if let Some(feature) = current_feature.as_mut() {
                    let body = line.trim().trim_start_matches('/').to_string();
                    let (k, v) = if let Some((k, v)) = body.split_once('=') {
                        (k.trim().to_string(), v.trim().trim_matches('"').to_string())
                    } else {
                        (body.trim().to_string(), String::new())
                    };
                    current_qualifier_key = Some(k.clone());
                    feature.qualifiers.entry(k).or_default().push(v);
                }
                continue;
            }

            if line.starts_with("                     ") {
                if let (Some(feature), Some(key)) =
                    (current_feature.as_mut(), current_qualifier_key.clone())
                {
                    if let Some(values) = feature.qualifiers.get_mut(&key) {
                        if let Some(last) = values.last_mut() {
                            let cont = line.trim().trim_matches('"');
                            if !cont.is_empty() {
                                if !last.is_empty() {
                                    last.push(' ');
                                }
                                last.push_str(cont);
                            }
                        }
                    }
                }
            }
        }
    }

    if let Some(f) = current_feature.take() {
        features.push(f);
    }

    if accession.is_empty() || sequence.is_empty() {
        return None;
    }
    Some(ParsedRecord {
        accession,
        organism,
        sequence,
        features,
    })
}

fn feature_texts(feature: &ParsedFeature, fields: &[String]) -> Vec<String> {
    let mut out = Vec::new();
    for field in fields {
        if let Some(values) = feature.qualifiers.get(field) {
            out.extend(values.iter().cloned());
        }
    }
    out
}

fn process_genbank_chunk(
    chunk: &str,
    marker_rules: &[MarkerRule],
    acc_to_seqs: &mut HashMap<String, HashSet<String>>,
    out_fasta: &mut File,
    counters: &mut Counters,
    dup_accessions: &mut HashMap<String, u64>,
    emitted_records: &mut Vec<EmittedRecord>,
    source_merge_rows: &mut Vec<EmittedRecord>,
    dump_gb_taxid_dir: Option<&Path>,
) -> Result<(), String> {
    for raw_record in split_genbank_records(chunk) {
        let Some(record) = parse_genbank_record(&raw_record) else {
            continue;
        };
        counters.total_records += 1;

        if let Some(dir) = dump_gb_taxid_dir {
            fs::create_dir_all(dir)
                .map_err(|e| format!("failed to create {}: {e}", dir.display()))?;
            let gb_path = dir.join(format!("{}.gb", sanitize_header(&record.accession)));
            if !gb_path.exists() {
                let mut f = File::create(&gb_path)
                    .map_err(|e| format!("failed to create {}: {e}", gb_path.display()))?;
                f.write_all(raw_record.as_bytes())
                    .and_then(|_| f.write_all(b"\n//\n"))
                    .map_err(|e| format!("failed to write {}: {e}", gb_path.display()))?;
            }
        }

        let mut record_matched = false;
        for feature in &record.features {
            let mut matched_label: Option<String> = None;
            let mut matched_marker = String::new();
            let mut header_format = DEFAULT_HEADER_FORMAT.to_string();
            for rule in marker_rules {
                if let Some(types) = &rule.feature_types {
                    if !types.iter().any(|t| t == &feature.feature_type) {
                        continue;
                    }
                }
                let texts = feature_texts(feature, &rule.feature_fields);
                let mut found = None;
                'outer: for t in texts {
                    for pat in &rule.patterns {
                        if pat.is_match(&t) {
                            found = Some(t);
                            break 'outer;
                        }
                    }
                }
                if let Some(label) = found {
                    matched_label = Some(label);
                    matched_marker = rule.key.clone();
                    header_format = rule.header_format.clone();
                    break;
                }
            }
            let Some(label) = matched_label else {
                continue;
            };
            let Some((start, end, strand)) = parse_location(&feature.location_raw) else {
                continue;
            };
            if end > record.sequence.len() || start == 0 || start > end {
                continue;
            }
            let mut seq = record.sequence[(start - 1)..end].to_string();
            if strand < 0 {
                seq = reverse_complement(&seq);
            }
            if seq.is_empty() {
                continue;
            }

            counters.matched_features += 1;
            record_matched = true;

            let seqs = acc_to_seqs.entry(record.accession.clone()).or_default();
            if seqs.contains(&seq) {
                counters.skipped_same += 1;
                continue;
            }

            let dup_index = if seqs.is_empty() {
                None
            } else {
                counters.duplicated_diff += 1;
                let idx = seqs.len() + 1;
                dup_accessions.insert(record.accession.clone(), idx as u64);
                Some(idx)
            };
            seqs.insert(seq.clone());

            let acc_id = dup_index
                .map(|idx| format!("{}_dup{idx}", record.accession))
                .unwrap_or_else(|| record.accession.clone());
            let loc = format!("{start}-{end}");
            let mut vars = HashMap::new();
            vars.insert("acc", record.accession.clone());
            vars.insert("acc_id", acc_id.clone());
            vars.insert("db", "gb".to_string());
            vars.insert("organism", sanitize_header(&record.organism));
            vars.insert("organism_raw", record.organism.clone());
            vars.insert("marker", sanitize_header(&matched_marker));
            vars.insert("marker_raw", matched_marker.clone());
            vars.insert("label", sanitize_header(&label));
            vars.insert("label_raw", label.clone());
            vars.insert("type", sanitize_header(&feature.feature_type));
            vars.insert("type_raw", feature.feature_type.clone());
            vars.insert("start", start.to_string());
            vars.insert("end", end.to_string());
            vars.insert("loc", loc);
            vars.insert("strand", strand.to_string());
            vars.insert(
                "dup",
                dup_index.map(|i| format!("dup{i}")).unwrap_or_default(),
            );
            let mut header = build_header(&header_format, &vars).trim().to_string();
            if header.is_empty() {
                header = acc_id.clone();
            }
            writeln!(out_fasta, ">{header}\n{seq}")
                .map_err(|e| format!("failed to write output fasta: {e}"))?;

            counters.kept_records += 1;
            emitted_records.push(EmittedRecord {
                acc_id,
                accession: record.accession.clone(),
                organism_name: record.organism.clone(),
                header,
                source: "ncbi".to_string(),
                source_record_id: record.accession.clone(),
                processid: String::new(),
                sampleid: String::new(),
                marker_key: matched_marker.clone(),
                linked_to_ncbi: false,
                emitted_to_fasta: true,
                skip_reason: String::new(),
            });
            source_merge_rows.push(emitted_records.last().cloned().unwrap_or_else(|| {
                EmittedRecord {
                    acc_id: String::new(),
                    accession: String::new(),
                    organism_name: String::new(),
                    header: String::new(),
                    source: "ncbi".to_string(),
                    source_record_id: String::new(),
                    processid: String::new(),
                    sampleid: String::new(),
                    marker_key: String::new(),
                    linked_to_ncbi: false,
                    emitted_to_fasta: false,
                    skip_reason: String::new(),
                }
            }));
        }

        if record_matched {
            counters.matched_records += 1;
        }
    }
    Ok(())
}

fn write_acc_organism_csv(
    out_fasta: &Path,
    emitted_records: &[EmittedRecord],
) -> Result<(), String> {
    let csv_path = PathBuf::from(format!("{}.acc_organism.csv", out_fasta.to_string_lossy()));
    let mut f = File::create(&csv_path)
        .map_err(|e| format!("failed to create {}: {e}", csv_path.display()))?;
    writeln!(
        f,
        "acc_id,accession,organism_name,header,source,source_record_id,processid,sampleid,marker_key,linked_to_ncbi,emitted_to_fasta,skip_reason"
    )
        .map_err(|e| format!("failed to write {}: {e}", csv_path.display()))?;
    for row in emitted_records {
        let esc = |s: &str| -> String {
            let s = s.replace('"', "\"\"");
            format!("\"{s}\"")
        };
        writeln!(
            f,
            "{},{},{},{},{},{},{},{},{},{},{},{}",
            esc(&row.acc_id),
            esc(&row.accession),
            esc(&row.organism_name),
            esc(&row.header),
            esc(&row.source),
            esc(&row.source_record_id),
            esc(&row.processid),
            esc(&row.sampleid),
            esc(&row.marker_key),
            esc(if row.linked_to_ncbi { "true" } else { "false" }),
            esc(if row.emitted_to_fasta {
                "true"
            } else {
                "false"
            }),
            esc(&row.skip_reason)
        )
        .map_err(|e| format!("failed to write {}: {e}", csv_path.display()))?;
    }
    Ok(())
}

fn write_source_merge_csv(out_fasta: &Path, rows: &[EmittedRecord]) -> Result<PathBuf, String> {
    let csv_path = PathBuf::from(format!("{}.source_merge.csv", out_fasta.to_string_lossy()));
    let mut f = File::create(&csv_path)
        .map_err(|e| format!("failed to create {}: {e}", csv_path.display()))?;
    writeln!(
        f,
        "source,source_record_id,accession,processid,sampleid,organism_name,marker_key,acc_id,linked_to_ncbi,emitted_to_fasta,skip_reason,header"
    )
    .map_err(|e| format!("failed to write {}: {e}", csv_path.display()))?;
    for row in rows {
        let esc = |s: &str| -> String {
            let s = s.replace('"', "\"\"");
            format!("\"{s}\"")
        };
        writeln!(
            f,
            "{},{},{},{},{},{},{},{},{},{},{},{}",
            esc(&row.source),
            esc(&row.source_record_id),
            esc(&row.accession),
            esc(&row.processid),
            esc(&row.sampleid),
            esc(&row.organism_name),
            esc(&row.marker_key),
            esc(&row.acc_id),
            esc(if row.linked_to_ncbi { "true" } else { "false" }),
            esc(if row.emitted_to_fasta {
                "true"
            } else {
                "false"
            }),
            esc(&row.skip_reason),
            esc(&row.header)
        )
        .map_err(|e| format!("failed to write {}: {e}", csv_path.display()))?;
    }
    Ok(csv_path)
}

fn emit_bold_record(
    record: &BoldRecord,
    out_fasta: &mut File,
    emitted_records: &mut Vec<EmittedRecord>,
    source_merge_rows: &mut Vec<EmittedRecord>,
    counters: &mut Counters,
) -> Result<(), String> {
    let acc_id = format!("BOLD_{}", sanitize_header(&record.source_record_id));
    let mut vars = HashMap::new();
    vars.insert("acc", record.accession.clone());
    vars.insert("acc_id", acc_id.clone());
    vars.insert("db", "bold".to_string());
    vars.insert("organism", sanitize_header(&record.taxon_name));
    vars.insert("organism_raw", record.taxon_name.clone());
    vars.insert("marker", sanitize_header(&record.marker_key));
    vars.insert("marker_raw", record.marker_key.clone());
    vars.insert("label", sanitize_header(&record.marker_label));
    vars.insert("label_raw", record.marker_label.clone());
    vars.insert("type", "barcode".to_string());
    vars.insert("type_raw", "barcode".to_string());
    vars.insert("start", String::new());
    vars.insert("end", String::new());
    vars.insert("loc", String::new());
    vars.insert("strand", String::new());
    vars.insert("dup", String::new());
    vars.insert("source", "bold".to_string());
    vars.insert("source_id", record.source_record_id.clone());
    let mut header = build_header(&record.header_format, &vars)
        .trim()
        .to_string();
    if header.is_empty() {
        header = acc_id.clone();
    }
    writeln!(out_fasta, ">{header}\n{}", record.sequence)
        .map_err(|e| format!("failed to write output fasta: {e}"))?;

    counters.kept_records += 1;
    let row = EmittedRecord {
        acc_id,
        accession: record.accession.clone(),
        organism_name: record.taxon_name.clone(),
        header,
        source: "bold".to_string(),
        source_record_id: record.source_record_id.clone(),
        processid: record.processid.clone(),
        sampleid: record.sampleid.clone(),
        marker_key: record.marker_key.clone(),
        linked_to_ncbi: false,
        emitted_to_fasta: true,
        skip_reason: String::new(),
    };
    emitted_records.push(row.clone());
    source_merge_rows.push(row);
    Ok(())
}

pub fn run_build(
    params: &BuildParams,
    log_path: &Path,
    cancelled: &AtomicBool,
    _child_slot: &Arc<Mutex<Option<Child>>>,
) -> Result<(), String> {
    let started_at = log_timestamp_now();
    let started_clock = Instant::now();
    let run_result = (|| -> Result<(), String> {
        let source = normalize_build_source(&params.source);
        let uses_ncbi = source_uses_ncbi(&source);
        let uses_bold = source_uses_bold(&source);
        let cfg = parse_toml(&params.config_path)?;
        let cfg_tbl = cfg
            .as_table()
            .ok_or_else(|| "config root must be a TOML table".to_string())?;
        let ncbi = cfg_tbl.get("ncbi").and_then(|v| v.as_table());
        if uses_ncbi && ncbi.is_none() {
            return Err("missing [ncbi] section".to_string());
        }
        let bold_cfg = cfg_tbl.get("bold").and_then(|v| v.as_table());
        let markers_tbl = load_markers_table(cfg_tbl, &params.config_path)?;
        let output_tbl = cfg_tbl.get("output").and_then(|v| v.as_table());
        let filters_tbl = cfg_tbl.get("filters").and_then(|v| v.as_table());
        let taxon_noexp = cfg_tbl
            .get("taxon")
            .and_then(|v| v.as_table())
            .map(|t| toml_bool(t, "noexp", false))
            .unwrap_or(false);
        let empty_ncbi_tbl = toml::map::Map::new();
        let ncbi_tbl = ncbi.unwrap_or(&empty_ncbi_tbl);

        let db = toml_string(ncbi_tbl, "db", "nucleotide");
        let rettype = toml_string(ncbi_tbl, "rettype", "gb");
        let retmode = toml_string(ncbi_tbl, "retmode", "text");
        if rettype != "gb" && rettype != "gbwithparts" {
            return Err("ncbi.rettype must be 'gb' or 'gbwithparts'".to_string());
        }
        let per_query = toml_u64(ncbi_tbl, "per_query", 100).max(1);
        let delay_sec = toml_f64(ncbi_tbl, "delay_sec").unwrap_or(0.34);
        let email = toml_string(ncbi_tbl, "email", "");
        let api_key = toml_string(ncbi_tbl, "api_key", "");
        let bold_runtime_cfg = if uses_bold {
            Some(get_bold_runtime_config(bold_cfg)?)
        } else {
            None
        };

        let mut marker_keys = Vec::new();
        for m in &params.markers {
            marker_keys.push(resolve_marker_key(m, &markers_tbl)?);
        }

        let mut marker_terms = Vec::new();
        let mut marker_rules = Vec::new();
        for key in &marker_keys {
            let marker_cfg = markers_tbl
                .get(key)
                .and_then(|v| v.as_table())
                .ok_or_else(|| format!("markers.{key} must be a table"))?;
            if uses_ncbi {
                marker_terms.extend(marker_query_terms(marker_cfg));
            }

            let pats = region_patterns(marker_cfg)
                .into_iter()
                .map(|p| Regex::new(&p).map_err(|e| format!("invalid regex in markers.{key}: {e}")))
                .collect::<Result<Vec<_>, _>>()?;

            let feature_types = marker_cfg
                .get("feature_types")
                .and_then(|v| v.as_array())
                .map(|arr| {
                    arr.iter()
                        .filter_map(|v| v.as_str())
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>()
                })
                .filter(|v| !v.is_empty())
                .or_else(|| {
                    Some(
                        DEFAULT_FEATURE_TYPES
                            .iter()
                            .map(|s| s.to_string())
                            .collect(),
                    )
                });

            let feature_fields = marker_cfg
                .get("feature_fields")
                .and_then(|v| v.as_array())
                .map(|arr| {
                    arr.iter()
                        .filter_map(|v| v.as_str())
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>()
                })
                .filter(|v| !v.is_empty())
                .unwrap_or_else(|| {
                    DEFAULT_FEATURE_FIELDS
                        .iter()
                        .map(|s| s.to_string())
                        .collect()
                });

            if uses_ncbi {
                marker_rules.push(MarkerRule {
                    key: key.clone(),
                    patterns: pats,
                    feature_types,
                    feature_fields,
                    header_format: resolve_header_format(marker_cfg, output_tbl),
                });
            }
        }
        if uses_ncbi && marker_terms.is_empty() {
            return Err("no marker terms resolved".to_string());
        }
        let marker_query = if !uses_ncbi {
            String::new()
        } else if marker_terms.len() == 1 {
            marker_terms[0].clone()
        } else {
            format!(
                "({})",
                marker_terms
                    .iter()
                    .map(|t| format!("({t})"))
                    .collect::<Vec<_>>()
                    .join(" OR ")
            )
        };
        let filter_terms = build_filter_terms(filters_tbl)?;

        let mut resolved_taxa = Vec::new();
        for t in &params.taxids {
            resolved_taxa.push(resolve_taxon(t)?);
        }

        append_run_header(
            log_path,
            params,
            &source,
            &resolved_taxa,
            &marker_keys,
            &started_at,
        )?;

        let mut out = File::create(&params.output_file)
            .map_err(|e| format!("failed to create {}: {e}", params.output_file.display()))?;

        let mut acc_to_seqs: HashMap<String, HashSet<String>> = HashMap::new();
        let mut dup_accessions: HashMap<String, u64> = HashMap::new();
        let mut counters = Counters::default();
        let mut emitted_records = Vec::new();
        let mut source_merge_rows = Vec::new();

        for resolved_taxon in &resolved_taxa {
            if !uses_ncbi {
                break;
            }
            let taxid = resolved_taxon.taxid.clone();
            if cancelled.load(Ordering::Relaxed) {
                return Err("cancelled".to_string());
            }

            let query = build_query(&taxid, &marker_query, &filter_terms, taxon_noexp);
            append_log_line(log_path, &format!("# query taxid={taxid}: {query}"))?;

            let mut search_params = vec![
                ("db", db.clone()),
                ("term", query.clone()),
                ("retmax", "0".to_string()),
                ("retmode", "json".to_string()),
            ];
            if !email.is_empty() {
                search_params.push(("email", email.clone()));
            }
            if !api_key.is_empty() {
                search_params.push(("api_key", api_key.clone()));
            }
            let json = eutils_get_validated_json(
                "esearch.fcgi",
                &search_params,
                |json| {
                    json.get("esearchresult")
                        .and_then(|v| v.get("count"))
                        .and_then(|v| v.as_str())
                        .and_then(|s| s.parse::<u64>().ok())
                        .is_some()
                },
                "esearch count parse failed",
            )?;
            let count = json
                .get("esearchresult")
                .and_then(|v| v.get("count"))
                .and_then(|v| v.as_str())
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(0);
            append_log_line(log_path, &format!("# query count taxid={taxid}: {count}"))?;
            append_log_line(
                log_path,
                &format!("# fetch progress taxid={taxid}: 0/{count}"),
            )?;
            if count == 0 {
                continue;
            }

            let cache_root = params.dump_gb_dir.join(".cache");
            fs::create_dir_all(&cache_root)
                .map_err(|e| format!("failed to create {}: {e}", cache_root.display()))?;
            let taxid_dump_dir = params.dump_gb_dir.join(format!("taxid{taxid}"));

            let mut start = 0u64;
            while start < count {
                if cancelled.load(Ordering::Relaxed) {
                    return Err("cancelled".to_string());
                }
                let cache_path =
                    cache_root.join(format!("start{start:09}_count{per_query:04}.cache"));
                let chunk = if params.resume && cache_path.exists() {
                    fs::read_to_string(&cache_path)
                        .map_err(|e| format!("failed to read {}: {e}", cache_path.display()))?
                } else {
                    let mut ids_params = vec![
                        ("db", db.clone()),
                        ("term", query.clone()),
                        ("retstart", start.to_string()),
                        ("retmax", per_query.to_string()),
                        ("retmode", "json".to_string()),
                    ];
                    if !email.is_empty() {
                        ids_params.push(("email", email.clone()));
                    }
                    if !api_key.is_empty() {
                        ids_params.push(("api_key", api_key.clone()));
                    }
                    let ids_json = eutils_get_validated_json(
                        "esearch.fcgi",
                        &ids_params,
                        |json| {
                            json.get("esearchresult")
                                .and_then(|v| v.get("idlist"))
                                .and_then(|v| v.as_array())
                                .is_some()
                        },
                        "esearch idlist parse failed",
                    )?;
                    let ids = ids_json
                        .get("esearchresult")
                        .and_then(|v| v.get("idlist"))
                        .and_then(|v| v.as_array())
                        .ok_or_else(|| "esearch idlist parse failed".to_string())?
                        .iter()
                        .filter_map(|v| v.as_str())
                        .collect::<Vec<_>>();
                    if ids.is_empty() {
                        start = start.saturating_add(per_query);
                        continue;
                    }
                    let mut fetch_params = vec![
                        ("db", db.clone()),
                        ("rettype", rettype.clone()),
                        ("retmode", retmode.clone()),
                        ("id", ids.join(",")),
                    ];
                    if !email.is_empty() {
                        fetch_params.push(("email", email.clone()));
                    }
                    if !api_key.is_empty() {
                        fetch_params.push(("api_key", api_key.clone()));
                    }
                    // efetch with many IDs can exceed URL length when sent as GET.
                    // Use POST to avoid 414 URI Too Long.
                    let text = eutils_post_text("efetch.fcgi", &fetch_params)?;
                    let _ = fs::write(&cache_path, &text);
                    text
                };

                process_genbank_chunk(
                    &chunk,
                    &marker_rules,
                    &mut acc_to_seqs,
                    &mut out,
                    &mut counters,
                    &mut dup_accessions,
                    &mut emitted_records,
                    &mut source_merge_rows,
                    Some(&taxid_dump_dir),
                )?;

                let fetched = (start + per_query).min(count);
                append_log_line(
                    log_path,
                    &format!("# fetch progress taxid={taxid}: {fetched}/{count}"),
                )?;
                start = start.saturating_add(per_query);
                if delay_sec > 0.0 {
                    std::thread::sleep(std::time::Duration::from_secs_f64(delay_sec));
                }
            }
        }

        if uses_bold {
            let mut bold_records = Vec::new();
            for resolved_taxon in &resolved_taxa {
                if cancelled.load(Ordering::Relaxed) {
                    return Err("cancelled".to_string());
                }
                let (rows, stats) = fetch_bold_records_for_taxon_with_logging(
                    &resolved_taxon.scientific_name,
                    &marker_keys,
                    &markers_tbl,
                    output_tbl,
                    bold_runtime_cfg
                        .as_ref()
                        .ok_or_else(|| "missing BOLD runtime config".to_string())?,
                    log_path,
                )?;
                counters.total_records += stats.downloaded_rows;
                counters.matched_records += stats.matched_rows;
                counters.matched_features += stats.matched_rows;
                append_log_line(
                log_path,
                &format!(
                    "# bold query: taxon={} normalized={} specimens={} format={} bytes={} downloaded={} matched={}",
                    resolved_taxon.scientific_name,
                    stats.normalized_query,
                    stats.specimen_count,
                    stats.download_format,
                    stats.downloaded_bytes,
                    stats.downloaded_rows,
                    stats.matched_rows
                ),
            )?;
                bold_records.extend(rows);
            }
            bold_records.sort_by(|left, right| {
                left.taxon_name
                    .cmp(&right.taxon_name)
                    .then(left.marker_key.cmp(&right.marker_key))
                    .then(left.source_record_id.cmp(&right.source_record_id))
            });

            let ncbi_accessions = acc_to_seqs.keys().cloned().collect::<HashSet<_>>();
            for record in bold_records {
                let accession_tokens = parse_accession_tokens(&record.accession);
                if source == "both"
                    && !accession_tokens.is_empty()
                    && accession_tokens
                        .iter()
                        .any(|token| ncbi_accessions.contains(token))
                {
                    source_merge_rows.push(EmittedRecord {
                        acc_id: format!("BOLD_{}", sanitize_header(&record.source_record_id)),
                        accession: record.accession.clone(),
                        organism_name: record.taxon_name.clone(),
                        header: String::new(),
                        source: "bold".to_string(),
                        source_record_id: record.source_record_id.clone(),
                        processid: record.processid.clone(),
                        sampleid: record.sampleid.clone(),
                        marker_key: record.marker_key.clone(),
                        linked_to_ncbi: true,
                        emitted_to_fasta: false,
                        skip_reason: "linked_by_insdcacs".to_string(),
                    });
                    continue;
                }
                emit_bold_record(
                    &record,
                    &mut out,
                    &mut emitted_records,
                    &mut source_merge_rows,
                    &mut counters,
                )?;
            }
        }

        write_acc_organism_csv(&params.output_file, &emitted_records)?;
        let source_merge_path = write_source_merge_csv(&params.output_file, &source_merge_rows)?;
        append_log_line(
            log_path,
            &format!("# source_merge_csv: {}", source_merge_path.display()),
        )?;
        append_log_line(
            log_path,
            &format!("# total records: {}", counters.total_records),
        )?;
        append_log_line(
            log_path,
            &format!("# matched records: {}", counters.matched_records),
        )?;
        append_log_line(
            log_path,
            &format!("# matched features: {}", counters.matched_features),
        )?;
        append_log_line(
            log_path,
            &format!("# kept records: {}", counters.kept_records),
        )?;
        append_log_line(
            log_path,
            &format!(
                "# skipped duplicates (same accession+sequence): {}",
                counters.skipped_same
            ),
        )?;
        append_log_line(
            log_path,
            &format!(
                "# kept duplicates (same accession, different sequence): {}",
                counters.duplicated_diff
            ),
        )?;
        if !dup_accessions.is_empty() {
            append_log_line(log_path, "# duplicate accessions with different sequences:")?;
            let mut rows = dup_accessions.into_iter().collect::<Vec<_>>();
            rows.sort_by(|a, b| a.0.cmp(&b.0));
            for (acc, cnt) in rows {
                append_log_line(log_path, &format!("# - {acc}: {cnt} sequences"))?;
            }
        }
        append_log_line(
            log_path,
            &format!("# output: {}", params.output_file.display()),
        )?;
        Ok(())
    })();

    let finished_at = log_timestamp_now();
    let elapsed_sec = started_clock.elapsed().as_secs_f64();
    match &run_result {
        Ok(()) => {
            let _ = append_log_line(log_path, "# status: ok");
        }
        Err(err) => {
            let _ = append_log_line(log_path, "# status: error");
            let _ = append_error_log_line(log_path, &format!("# error: {err}"));
        }
    }
    let _ = append_log_line(log_path, &format!("# finished_at: {finished_at}"));
    let _ = append_log_line(log_path, &format!("# elapsed_sec: {elapsed_sec:.3}"));

    run_result
}
