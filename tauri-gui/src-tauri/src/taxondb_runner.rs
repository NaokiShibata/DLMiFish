use regex::Regex;
use serde_json::Value as JsonValue;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicBool, Ordering};
use toml::Value as TomlValue;

const DEFAULT_FEATURE_TYPES: [&str; 3] = ["rRNA", "gene", "CDS"];
const DEFAULT_FEATURE_FIELDS: [&str; 4] = ["gene", "product", "note", "standard_name"];
const DEFAULT_HEADER_FORMAT: &str = "{acc_id}|{organism}|{marker}|{label}|{type}|{loc}|{strand}";

#[derive(Debug, Clone)]
pub struct BuildParams {
    pub config_path: PathBuf,
    pub taxids: Vec<String>,
    pub markers: Vec<String>,
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
}

fn append_log_line(log_path: &Path, line: &str) -> Result<(), String> {
    let mut f = OpenOptions::new()
        .create(true)
        .append(true)
        .open(log_path)
        .map_err(|e| format!("failed to open log {}: {e}", log_path.display()))?;
    writeln!(f, "{line}").map_err(|e| format!("failed to write log {}: {e}", log_path.display()))
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

fn eutils_get_json(endpoint: &str, params: &[(&str, String)]) -> Result<JsonValue, String> {
    let text = eutils_get_text(endpoint, params)?;
    serde_json::from_str::<JsonValue>(&text).map_err(|e| format!("eutils json parse failed: {e}"))
}

fn eutils_get_text(endpoint: &str, params: &[(&str, String)]) -> Result<String, String> {
    let url = format!("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}");
    let mut cmd = Command::new("curl");
    cmd.arg("-fsSL").arg(url);
    for (k, v) in params {
        cmd.arg("--get")
            .arg("--data-urlencode")
            .arg(format!("{k}={v}"));
    }
    let out = cmd
        .output()
        .map_err(|e| format!("failed to execute curl: {e}"))?;
    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        return Err(format!("curl failed: {stderr}"));
    }
    String::from_utf8(out.stdout).map_err(|e| format!("curl output decode failed: {e}"))
}

fn resolve_taxid(taxon: &str) -> Result<String, String> {
    if taxon.chars().all(|c| c.is_ascii_digit()) {
        return Ok(taxon.to_string());
    }
    let term = format!("\"{taxon}\"[Scientific Name]");
    let json = eutils_get_json(
        "esearch.fcgi",
        &[
            ("db", "taxonomy".to_string()),
            ("term", term),
            ("retmax", "5".to_string()),
            ("retmode", "json".to_string()),
        ],
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
            });
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
    writeln!(f, "acc_id,accession,organism_name,header")
        .map_err(|e| format!("failed to write {}: {e}", csv_path.display()))?;
    for row in emitted_records {
        let esc = |s: &str| -> String {
            let s = s.replace('"', "\"\"");
            format!("\"{s}\"")
        };
        writeln!(
            f,
            "{},{},{},{}",
            esc(&row.acc_id),
            esc(&row.accession),
            esc(&row.organism_name),
            esc(&row.header)
        )
        .map_err(|e| format!("failed to write {}: {e}", csv_path.display()))?;
    }
    Ok(())
}

pub fn run_build(
    params: &BuildParams,
    log_path: &Path,
    cancelled: &AtomicBool,
) -> Result<(), String> {
    let cfg = parse_toml(&params.config_path)?;
    let cfg_tbl = cfg
        .as_table()
        .ok_or_else(|| "config root must be a TOML table".to_string())?;
    let ncbi = cfg_tbl
        .get("ncbi")
        .and_then(|v| v.as_table())
        .ok_or_else(|| "missing [ncbi] section".to_string())?;
    let markers_tbl = load_markers_table(cfg_tbl, &params.config_path)?;
    let output_tbl = cfg_tbl.get("output").and_then(|v| v.as_table());
    let filters_tbl = cfg_tbl.get("filters").and_then(|v| v.as_table());
    let taxon_noexp = cfg_tbl
        .get("taxon")
        .and_then(|v| v.as_table())
        .map(|t| toml_bool(t, "noexp", false))
        .unwrap_or(false);

    let db = toml_string(ncbi, "db", "nucleotide");
    let rettype = toml_string(ncbi, "rettype", "gb");
    let retmode = toml_string(ncbi, "retmode", "text");
    if rettype != "gb" && rettype != "gbwithparts" {
        return Err("ncbi.rettype must be 'gb' or 'gbwithparts'".to_string());
    }
    let per_query = toml_u64(ncbi, "per_query", 100).max(1);
    let delay_sec = toml_f64(ncbi, "delay_sec").unwrap_or(0.34);
    let email = toml_string(ncbi, "email", "");
    let api_key = toml_string(ncbi, "api_key", "");

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
        marker_terms.extend(marker_query_terms(marker_cfg));

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

        marker_rules.push(MarkerRule {
            key: key.clone(),
            patterns: pats,
            feature_types,
            feature_fields,
            header_format: resolve_header_format(marker_cfg, output_tbl),
        });
    }
    if marker_terms.is_empty() {
        return Err("no marker terms resolved".to_string());
    }
    let marker_query = if marker_terms.len() == 1 {
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

    let mut resolved_taxids = Vec::new();
    for t in &params.taxids {
        resolved_taxids.push(resolve_taxid(t)?);
    }

    append_log_line(log_path, "# started: rust-runner")?;
    append_log_line(
        log_path,
        &format!("# config: {}", params.config_path.display()),
    )?;
    append_log_line(log_path, &format!("# taxids: {:?}", resolved_taxids))?;
    append_log_line(log_path, &format!("# markers: {:?}", marker_keys))?;

    let mut out = File::create(&params.output_file)
        .map_err(|e| format!("failed to create {}: {e}", params.output_file.display()))?;

    let mut acc_to_seqs: HashMap<String, HashSet<String>> = HashMap::new();
    let mut dup_accessions: HashMap<String, u64> = HashMap::new();
    let mut counters = Counters::default();
    let mut emitted_records = Vec::new();

    for taxid in resolved_taxids {
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
        let json = eutils_get_json("esearch.fcgi", &search_params)?;
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
            let cache_path = cache_root.join(format!("start{start:09}_count{per_query:04}.cache"));
            let chunk = if params.resume && cache_path.exists() {
                fs::read_to_string(&cache_path)
                    .map_err(|e| format!("failed to read {}: {e}", cache_path.display()))?
            } else {
                let ids_json = eutils_get_json(
                    "esearch.fcgi",
                    &[
                        ("db", db.clone()),
                        ("term", query.clone()),
                        ("retstart", start.to_string()),
                        ("retmax", per_query.to_string()),
                        ("retmode", "json".to_string()),
                    ],
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
                let text = eutils_get_text("efetch.fcgi", &fetch_params)?;
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

    write_acc_organism_csv(&params.output_file, &emitted_records)?;
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
    append_log_line(log_path, "# finished: rust-runner")?;

    Ok(())
}
