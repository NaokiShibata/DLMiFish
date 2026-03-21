use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Debug, Clone)]
struct FastaRecord {
    header: String,
    seq: String,
}

#[derive(Debug, Clone)]
pub struct LengthFilterStats {
    pub before: usize,
    pub after: usize,
    pub removed: usize,
}

#[derive(Debug, Clone)]
pub struct PrimerTrimStats {
    pub before: usize,
    pub after: usize,
    pub removed: usize,
    pub trimmed_both: usize,
    pub trimmed_left_only: usize,
    pub trimmed_right_only: usize,
    pub untrimmed: usize,
    pub dropped_empty: usize,
    pub canonical_orientation: usize,
    pub reverse_orientation: usize,
    pub confidence_high: usize,
    pub confidence_medium: usize,
    pub confidence_low: usize,
    pub rounds_run: usize,
    pub best_round: usize,
    pub high_conf_rate: f64,
    pub sidecar_path: Option<String>,
    pub retained_path: Option<String>,
    pub recheck_attempted: usize,
    pub recheck_rescued: usize,
    pub recheck_error: Option<String>,
}

#[derive(Debug, Clone)]
pub struct PrimerSet {
    pub forward: Vec<String>,
    pub reverse: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct DuplicateAccStats {
    pub total_records: usize,
    pub parsed_records: usize,
    pub unparsed_records: usize,
    pub duplicate_groups: usize,
    pub duplicate_records: usize,
    pub cross_organism_groups: usize,
}

#[derive(Debug, Clone)]
pub struct PrimerTrimOptions {
    pub trim_mode: String,
    pub max_mismatch: usize,
    pub max_error_rate: f64,
    pub min_overlap_bp: Option<usize>,
    pub min_overlap_ratio: f64,
    pub end_max_offset: usize,
    pub keep_retained_fasta: bool,
    pub iter_enable: bool,
    pub iter_max_rounds: usize,
    pub iter_stop_delta: f64,
    pub iter_target_conf: f64,
    pub sidecar_format: String,
    pub recheck_tool: String,
    pub recheck_min_identity: f64,
    pub recheck_min_query_cov: f64,
}

impl Default for PrimerTrimOptions {
    fn default() -> Self {
        Self {
            trim_mode: "one_or_both".to_string(),
            max_mismatch: 0,
            max_error_rate: 0.0,
            min_overlap_bp: None,
            min_overlap_ratio: 1.0,
            end_max_offset: 0,
            keep_retained_fasta: true,
            iter_enable: false,
            iter_max_rounds: 3,
            iter_stop_delta: 0.002,
            iter_target_conf: 0.98,
            sidecar_format: "tsv".to_string(),
            recheck_tool: "off".to_string(),
            recheck_min_identity: 0.85,
            recheck_min_query_cov: 0.7,
        }
    }
}

#[derive(Debug, Clone)]
struct EndpointMatch {
    primer: String,
    overlap_bp: usize,
    mismatches: usize,
    score: usize,
    trim_bp: usize,
}

#[derive(Debug, Clone)]
struct OrientationScore {
    name: &'static str,
    left: Option<EndpointMatch>,
    right: Option<EndpointMatch>,
}

impl OrientationScore {
    fn matched_ends(&self) -> usize {
        usize::from(self.left.is_some()) + usize::from(self.right.is_some())
    }

    fn score_total(&self) -> usize {
        self.left.as_ref().map(|m| m.score).unwrap_or(0)
            + self.right.as_ref().map(|m| m.score).unwrap_or(0)
    }

    fn mismatch_total(&self) -> usize {
        self.left.as_ref().map(|m| m.mismatches).unwrap_or(0)
            + self.right.as_ref().map(|m| m.mismatches).unwrap_or(0)
    }

    fn rank_key(&self) -> (usize, usize, isize) {
        (
            self.matched_ends(),
            self.score_total(),
            -(self.mismatch_total() as isize),
        )
    }
}

#[derive(Debug, Clone)]
struct SidecarRow {
    record_id: String,
    round: usize,
    orientation_chosen: String,
    orientation_ambiguous: usize,
    left_hit: usize,
    right_hit: usize,
    left_overlap_bp: usize,
    right_overlap_bp: usize,
    left_trim_bp: usize,
    right_trim_bp: usize,
    left_mismatch: usize,
    right_mismatch: usize,
    left_primer_name: String,
    right_primer_name: String,
    trim_start: usize,
    trim_end: usize,
    trim_mode: String,
    confidence: String,
    dropped_empty: usize,
    recheck_tool: String,
    recheck_status: String,
    phylo_mismatch_flag: usize,
}

fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>, String> {
    let file = File::open(path).map_err(|e| format!("failed to open {}: {e}", path.display()))?;
    let reader = BufReader::new(file);

    let mut out = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("failed to read {}: {e}", path.display()))?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(h) = header.take() {
                out.push(FastaRecord {
                    header: h,
                    seq: seq.clone(),
                });
            }
            header = Some(rest.trim().to_string());
            seq.clear();
        } else {
            seq.push_str(line.trim());
        }
    }
    if let Some(h) = header.take() {
        out.push(FastaRecord { header: h, seq });
    }
    Ok(out)
}

fn write_fasta(path: &Path, records: &[FastaRecord]) -> Result<(), String> {
    let tmp = PathBuf::from(format!("{}.tmp", path.to_string_lossy()));
    let mut f =
        File::create(&tmp).map_err(|e| format!("failed to create {}: {e}", tmp.display()))?;
    for record in records {
        writeln!(f, ">{}", record.header)
            .map_err(|e| format!("failed to write {}: {e}", tmp.display()))?;
        writeln!(f, "{}", record.seq)
            .map_err(|e| format!("failed to write {}: {e}", tmp.display()))?;
    }
    f.flush()
        .map_err(|e| format!("failed to flush {}: {e}", tmp.display()))?;
    fs::rename(&tmp, path).map_err(|e| {
        format!(
            "failed to replace {} with {}: {e}",
            path.display(),
            tmp.display()
        )
    })
}

pub fn count_fasta_records(path: &Path) -> Result<usize, String> {
    Ok(read_fasta(path)?.len())
}

pub fn apply_length_filter(
    fasta_path: &Path,
    min_len: Option<u32>,
    max_len: Option<u32>,
) -> Result<LengthFilterStats, String> {
    let records = read_fasta(fasta_path)?;
    let before = records.len();
    let mut out = Vec::new();
    for rec in records {
        let len = rec.seq.len() as u32;
        if let Some(min_v) = min_len {
            if len < min_v {
                continue;
            }
        }
        if let Some(max_v) = max_len {
            if len > max_v {
                continue;
            }
        }
        out.push(rec);
    }
    let after = out.len();
    write_fasta(fasta_path, &out)?;
    Ok(LengthFilterStats {
        before,
        after,
        removed: before.saturating_sub(after),
    })
}

fn iupac_values(ch: char) -> Option<&'static str> {
    match ch.to_ascii_uppercase() {
        'A' => Some("A"),
        'C' => Some("C"),
        'G' => Some("G"),
        'T' | 'U' => Some("T"),
        'R' => Some("AG"),
        'Y' => Some("CT"),
        'S' => Some("GC"),
        'W' => Some("AT"),
        'K' => Some("GT"),
        'M' => Some("AC"),
        'B' => Some("CGT"),
        'D' => Some("AGT"),
        'H' => Some("ACT"),
        'V' => Some("ACG"),
        'N' => Some("ACGT"),
        _ => None,
    }
}

fn iupac_complement(ch: char) -> Option<char> {
    match ch.to_ascii_uppercase() {
        'A' => Some('T'),
        'T' | 'U' => Some('A'),
        'G' => Some('C'),
        'C' => Some('G'),
        'R' => Some('Y'),
        'Y' => Some('R'),
        'S' => Some('S'),
        'W' => Some('W'),
        'K' => Some('M'),
        'M' => Some('K'),
        'B' => Some('V'),
        'V' => Some('B'),
        'D' => Some('H'),
        'H' => Some('D'),
        'N' => Some('N'),
        _ => None,
    }
}

fn reverse_complement_iupac(seq: &str) -> Result<String, String> {
    let mut out = String::with_capacity(seq.len());
    for ch in seq.chars().rev() {
        let c = iupac_complement(ch).ok_or_else(|| format!("unsupported primer base: {ch}"))?;
        out.push(c);
    }
    Ok(out)
}

fn base_matches(seq_base: char, primer_base: char) -> bool {
    let Some(values) = iupac_values(primer_base) else {
        return false;
    };
    values.contains(seq_base)
}

fn required_overlap_bp(
    primer_len: usize,
    min_overlap_bp: Option<usize>,
    min_overlap_ratio: f64,
) -> usize {
    let mut ratio_required = (primer_len as f64 * min_overlap_ratio).round() as usize;
    if ratio_required < 1 {
        ratio_required = 1;
    }
    if ratio_required > primer_len {
        ratio_required = primer_len;
    }
    match min_overlap_bp {
        Some(v) => v.max(ratio_required).max(1).min(primer_len),
        None => ratio_required,
    }
}

fn count_mismatches(seq_seg: &str, primer_seg: &str) -> usize {
    seq_seg
        .chars()
        .zip(primer_seg.chars())
        .filter(|(s, p)| !base_matches(*s, *p))
        .count()
}

fn best_prefix_match(
    seq: &str,
    primers: &[String],
    max_mismatch: usize,
    max_error_rate: f64,
    min_overlap_bp: Option<usize>,
    min_overlap_ratio: f64,
    max_offset: usize,
) -> Option<EndpointMatch> {
    let mut best: Option<EndpointMatch> = None;
    for primer in primers {
        let max_overlap = primer.len().min(seq.len());
        let required = required_overlap_bp(primer.len(), min_overlap_bp, min_overlap_ratio);
        if max_overlap < required {
            continue;
        }
        let max_shift = max_offset.min(seq.len().saturating_sub(required));
        for offset in 0..=max_shift {
            let available = seq.len().saturating_sub(offset);
            if available < required {
                continue;
            }
            for overlap in (required..=max_overlap.min(available)).rev() {
                let primer_seg = &primer[primer.len() - overlap..];
                let seq_seg = &seq[offset..offset + overlap];
                let mismatches = count_mismatches(seq_seg, primer_seg);
                if mismatches > max_mismatch {
                    continue;
                }
                let error_rate = mismatches as f64 / overlap as f64;
                if error_rate > max_error_rate {
                    continue;
                }
                let cand = EndpointMatch {
                    primer: primer.clone(),
                    overlap_bp: overlap,
                    mismatches,
                    score: overlap.saturating_sub(mismatches),
                    trim_bp: offset + overlap,
                };
                if let Some(prev) = &best {
                    let prev_offset = prev.trim_bp.saturating_sub(prev.overlap_bp);
                    let ckey = (
                        cand.overlap_bp,
                        cand.score,
                        -(cand.mismatches as isize),
                        -(offset as isize),
                    );
                    let pkey = (
                        prev.overlap_bp,
                        prev.score,
                        -(prev.mismatches as isize),
                        -(prev_offset as isize),
                    );
                    if ckey > pkey {
                        best = Some(cand);
                    }
                } else {
                    best = Some(cand);
                }
            }
        }
    }
    best
}

fn best_suffix_match(
    seq: &str,
    primers: &[String],
    max_mismatch: usize,
    max_error_rate: f64,
    min_overlap_bp: Option<usize>,
    min_overlap_ratio: f64,
    max_offset: usize,
) -> Option<EndpointMatch> {
    let mut best: Option<EndpointMatch> = None;
    for primer in primers {
        let max_overlap = primer.len().min(seq.len());
        let required = required_overlap_bp(primer.len(), min_overlap_bp, min_overlap_ratio);
        if max_overlap < required {
            continue;
        }
        let max_shift = max_offset.min(seq.len().saturating_sub(required));
        for shift in 0..=max_shift {
            let end = seq.len().saturating_sub(shift);
            if end < required {
                continue;
            }
            for overlap in (required..=max_overlap.min(end)).rev() {
                let start = end - overlap;
                let primer_seg = &primer[..overlap];
                let seq_seg = &seq[start..end];
                let mismatches = count_mismatches(seq_seg, primer_seg);
                if mismatches > max_mismatch {
                    continue;
                }
                let error_rate = mismatches as f64 / overlap as f64;
                if error_rate > max_error_rate {
                    continue;
                }
                let cand = EndpointMatch {
                    primer: primer.clone(),
                    overlap_bp: overlap,
                    mismatches,
                    score: overlap.saturating_sub(mismatches),
                    trim_bp: seq.len().saturating_sub(start),
                };
                if let Some(prev) = &best {
                    let prev_shift = prev.trim_bp.saturating_sub(prev.overlap_bp);
                    let ckey = (
                        cand.overlap_bp,
                        cand.score,
                        -(cand.mismatches as isize),
                        -(shift as isize),
                    );
                    let pkey = (
                        prev.overlap_bp,
                        prev.score,
                        -(prev.mismatches as isize),
                        -(prev_shift as isize),
                    );
                    if ckey > pkey {
                        best = Some(cand);
                    }
                } else {
                    best = Some(cand);
                }
            }
        }
    }
    best
}

fn choose_orientation(
    canonical: &OrientationScore,
    reverse: &OrientationScore,
) -> (OrientationScore, bool) {
    let c = canonical.rank_key();
    let r = reverse.rank_key();
    if c > r {
        return (canonical.clone(), false);
    }
    if r > c {
        return (reverse.clone(), false);
    }
    (canonical.clone(), true)
}

fn confidence_label(
    matched_ends: usize,
    mismatch_total: usize,
    ambiguous_orientation: bool,
) -> String {
    if matched_ends == 2 && mismatch_total <= 1 && !ambiguous_orientation {
        return "high".to_string();
    }
    if matched_ends >= 1 {
        return "medium".to_string();
    }
    "low".to_string()
}

fn compute_trim_lengths(row: &SidecarRow, trim_mode: &str) -> (usize, usize) {
    let left_hit = row.left_hit > 0;
    let right_hit = row.right_hit > 0;
    let left_trim = row.left_trim_bp;
    let right_trim = row.right_trim_bp;
    if trim_mode == "one_or_both" {
        (
            if left_hit { left_trim } else { 0 },
            if right_hit { right_trim } else { 0 },
        )
    } else if trim_mode == "both_required" {
        if left_hit && right_hit {
            (left_trim, right_trim)
        } else {
            (0, 0)
        }
    } else {
        (0, 0)
    }
}

fn summarize_sidecar(
    rows: &mut [SidecarRow],
    seq_map: &HashMap<String, String>,
    trim_mode: &str,
) -> (Vec<FastaRecord>, PrimerTrimStats) {
    let mut records = Vec::new();
    let mut stats = PrimerTrimStats {
        before: rows.len(),
        after: 0,
        removed: 0,
        trimmed_both: 0,
        trimmed_left_only: 0,
        trimmed_right_only: 0,
        untrimmed: 0,
        dropped_empty: 0,
        canonical_orientation: 0,
        reverse_orientation: 0,
        confidence_high: 0,
        confidence_medium: 0,
        confidence_low: 0,
        rounds_run: 0,
        best_round: 0,
        high_conf_rate: 0.0,
        sidecar_path: None,
        retained_path: None,
        recheck_attempted: 0,
        recheck_rescued: 0,
        recheck_error: None,
    };

    for row in rows.iter_mut() {
        let Some(seq) = seq_map.get(&row.record_id) else {
            continue;
        };

        match row.confidence.as_str() {
            "high" => stats.confidence_high += 1,
            "medium" => stats.confidence_medium += 1,
            _ => stats.confidence_low += 1,
        }

        let matched_ends = row.left_hit + row.right_hit;
        if matched_ends > 0 {
            if row.orientation_chosen == "canonical" {
                stats.canonical_orientation += 1;
            } else {
                stats.reverse_orientation += 1;
            }
        }

        let (left_len, right_len) = compute_trim_lengths(row, trim_mode);
        row.trim_start = left_len;
        let right_idx = if right_len > 0 {
            seq.len().saturating_sub(right_len)
        } else {
            seq.len()
        };
        row.trim_end = right_idx;

        let ends_trimmed = usize::from(left_len > 0) + usize::from(right_len > 0);
        if ends_trimmed == 0 {
            stats.untrimmed += 1;
        } else if ends_trimmed == 2 {
            stats.trimmed_both += 1;
        } else if left_len > 0 {
            stats.trimmed_left_only += 1;
        } else {
            stats.trimmed_right_only += 1;
        }

        if left_len >= right_idx {
            row.dropped_empty = 1;
            stats.dropped_empty += 1;
            continue;
        }
        row.dropped_empty = 0;
        let trimmed = seq[left_len..right_idx].to_string();
        records.push(FastaRecord {
            header: row.record_id.clone(),
            seq: trimmed,
        });
    }

    stats.after = records.len();
    stats.removed = stats.before.saturating_sub(stats.after);
    stats.high_conf_rate = if stats.before > 0 {
        stats.confidence_high as f64 / stats.before as f64
    } else {
        0.0
    };
    (records, stats)
}

fn run_vsearch_recheck(
    rows: &mut [SidecarRow],
    seq_map: &HashMap<String, String>,
    forward_primers: &[String],
    reverse_primers: &[String],
    forward_rc: &[String],
    reverse_rc: &[String],
    min_identity: f64,
    min_query_cov: f64,
) -> Result<(usize, usize), String> {
    let Some(vsearch_bin) = std::env::var_os("PATH").and_then(|_| which_like("vsearch")) else {
        return Err("vsearch_not_found".to_string());
    };

    let mut candidates: Vec<(usize, bool, String)> = Vec::new(); // row_idx, is_left, qid
    let mut attempted = 0usize;

    let mut db_entries: Vec<(String, String)> = Vec::new();
    for (i, p) in forward_primers.iter().enumerate() {
        db_entries.push((format!("CL_{i}"), p.clone()));
    }
    for (i, p) in reverse_rc.iter().enumerate() {
        db_entries.push((format!("CR_{i}"), p.clone()));
    }
    for (i, p) in reverse_primers.iter().enumerate() {
        db_entries.push((format!("RL_{i}"), p.clone()));
    }
    for (i, p) in forward_rc.iter().enumerate() {
        db_entries.push((format!("RR_{i}"), p.clone()));
    }

    let max_primer_len = db_entries.iter().map(|(_, s)| s.len()).max().unwrap_or(30);
    let window_len = (max_primer_len + 8).max(20);

    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| e.to_string())?
        .as_nanos();
    let tmp_dir = std::env::temp_dir().join(format!("taxondb-vsearch-{ts}"));
    fs::create_dir_all(&tmp_dir)
        .map_err(|e| format!("failed to create {}: {e}", tmp_dir.display()))?;
    let db_fa = tmp_dir.join("primers.fa");
    let q_fa = tmp_dir.join("queries.fa");
    let out_tsv = tmp_dir.join("hits.tsv");

    {
        let mut f = File::create(&db_fa)
            .map_err(|e| format!("failed to write {}: {e}", db_fa.display()))?;
        for (id, seq) in &db_entries {
            writeln!(f, ">{id}")
                .and_then(|_| writeln!(f, "{seq}"))
                .map_err(|e| format!("failed to write {}: {e}", db_fa.display()))?;
        }
    }

    {
        let mut f =
            File::create(&q_fa).map_err(|e| format!("failed to write {}: {e}", q_fa.display()))?;
        for (idx, row) in rows.iter_mut().enumerate() {
            if row.confidence != "low" && row.confidence != "medium" {
                row.recheck_status = "not_target_confidence".to_string();
                continue;
            }
            if row.left_hit > 0 && row.right_hit > 0 {
                row.recheck_status = "already_both_ends".to_string();
                continue;
            }
            let Some(seq) = seq_map.get(&row.record_id) else {
                row.recheck_status = "sequence_not_found".to_string();
                continue;
            };
            if row.left_hit == 0 {
                let qid = format!("{idx}|L");
                let seg = &seq[..seq.len().min(window_len)];
                writeln!(f, ">{qid}")
                    .and_then(|_| writeln!(f, "{seg}"))
                    .map_err(|e| format!("failed to write {}: {e}", q_fa.display()))?;
                candidates.push((idx, true, qid));
                attempted += 1;
            }
            if row.right_hit == 0 {
                let qid = format!("{idx}|R");
                let from = seq.len().saturating_sub(window_len);
                let seg = &seq[from..];
                writeln!(f, ">{qid}")
                    .and_then(|_| writeln!(f, "{seg}"))
                    .map_err(|e| format!("failed to write {}: {e}", q_fa.display()))?;
                candidates.push((idx, false, qid));
                attempted += 1;
            }
        }
    }

    if attempted == 0 {
        let _ = fs::remove_dir_all(&tmp_dir);
        return Ok((0, 0));
    }

    let status = Command::new(vsearch_bin)
        .arg("--usearch_global")
        .arg(&q_fa)
        .arg("--db")
        .arg(&db_fa)
        .arg("--id")
        .arg(format!("{min_identity:.4}"))
        .arg("--strand")
        .arg("plus")
        .arg("--blast6out")
        .arg(&out_tsv)
        .arg("--maxaccepts")
        .arg("16")
        .arg("--maxrejects")
        .arg("64")
        .status()
        .map_err(|e| format!("failed to run vsearch: {e}"))?;
    if !status.success() {
        let _ = fs::remove_dir_all(&tmp_dir);
        return Err("vsearch_failed".to_string());
    }

    let primer_len: HashMap<String, usize> = db_entries
        .iter()
        .map(|(id, seq)| (id.clone(), seq.len()))
        .collect();

    let mut best: HashMap<String, (String, f64, f64, usize, usize)> = HashMap::new();
    if out_tsv.exists() {
        let file = File::open(&out_tsv)
            .map_err(|e| format!("failed to read {}: {e}", out_tsv.display()))?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let line = line.map_err(|e| format!("failed to read {}: {e}", out_tsv.display()))?;
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 12 {
                continue;
            }
            let qid = cols[0].to_string();
            let sid = cols[1].to_string();
            let pident = cols[2].parse::<f64>().unwrap_or(0.0) / 100.0;
            let aln_len = cols[3].parse::<usize>().unwrap_or(0);
            let mismatch = cols[4].parse::<usize>().unwrap_or(0);
            let Some(plen) = primer_len.get(&sid).copied() else {
                continue;
            };
            if plen == 0 {
                continue;
            }
            let cov = aln_len as f64 / plen as f64;
            if pident < min_identity || cov < min_query_cov {
                continue;
            }
            let cand_key = (pident, cov, aln_len as i64, -(mismatch as i64));
            let replace =
                if let Some((_, prev_pident, prev_cov, prev_aln, prev_mm)) = best.get(&qid) {
                    let prev_key = (
                        *prev_pident,
                        *prev_cov,
                        *prev_aln as i64,
                        -(*prev_mm as i64),
                    );
                    cand_key > prev_key
                } else {
                    true
                };
            if replace {
                best.insert(qid, (sid, pident, cov, aln_len, mismatch));
            }
        }
    }

    let mut rescued = 0usize;
    for (idx, is_left, qid) in candidates {
        let row = &mut rows[idx];
        let Some((sid, _pident, _cov, aln_len, mismatch)) = best.get(&qid).cloned() else {
            if row.recheck_status != "rescued_by_vsearch" {
                row.recheck_status = "attempted_no_hit".to_string();
            }
            continue;
        };

        let is_left_subject = sid.starts_with("CL_") || sid.starts_with("RL_");
        let is_right_subject = sid.starts_with("CR_") || sid.starts_with("RR_");
        if is_left && !is_left_subject {
            continue;
        }
        if !is_left && !is_right_subject {
            continue;
        }

        if is_left {
            row.left_hit = 1;
            row.left_overlap_bp = aln_len;
            row.left_trim_bp = aln_len;
            row.left_mismatch = mismatch;
        } else {
            row.right_hit = 1;
            row.right_overlap_bp = aln_len;
            row.right_trim_bp = aln_len;
            row.right_mismatch = mismatch;
        }
        row.recheck_status = "rescued_by_vsearch".to_string();
        row.confidence = confidence_label(
            row.left_hit + row.right_hit,
            row.left_mismatch + row.right_mismatch,
            row.orientation_ambiguous > 0,
        );
        rescued += 1;
    }

    let _ = fs::remove_dir_all(&tmp_dir);
    Ok((attempted, rescued))
}

fn which_like(bin: &str) -> Option<PathBuf> {
    let paths = std::env::var_os("PATH")?;
    for p in std::env::split_paths(&paths) {
        let cand = p.join(bin);
        if cand.exists() {
            return Some(cand);
        }
    }
    None
}

pub fn apply_primer_trim(
    fasta_path: &Path,
    forward_primers: &[String],
    reverse_primers: &[String],
    options: &PrimerTrimOptions,
) -> Result<PrimerTrimStats, String> {
    let records = read_fasta(fasta_path)?;
    let before = records.len();

    let mut reverse_comp_forward = Vec::with_capacity(forward_primers.len());
    for p in forward_primers {
        reverse_comp_forward.push(reverse_complement_iupac(p)?);
    }
    let mut reverse_comp_reverse = Vec::with_capacity(reverse_primers.len());
    for p in reverse_primers {
        reverse_comp_reverse.push(reverse_complement_iupac(p)?);
    }

    let retained_path = PathBuf::from(format!(
        "{}.postprep.primer.retained.fasta",
        fasta_path.to_string_lossy()
    ));
    if options.keep_retained_fasta {
        fs::copy(fasta_path, &retained_path)
            .map_err(|e| format!("failed to write {}: {e}", retained_path.display()))?;
    }

    let mut seq_map: HashMap<String, String> = HashMap::new();
    let mut normalized_records: Vec<FastaRecord> = Vec::new();
    for mut rec in records {
        rec.seq = rec.seq.to_ascii_uppercase().replace('U', "T");
        seq_map.insert(rec.header.clone(), rec.seq.clone());
        normalized_records.push(rec);
    }

    if before == 0 {
        return Ok(PrimerTrimStats {
            before: 0,
            after: 0,
            removed: 0,
            trimmed_both: 0,
            trimmed_left_only: 0,
            trimmed_right_only: 0,
            untrimmed: 0,
            dropped_empty: 0,
            canonical_orientation: 0,
            reverse_orientation: 0,
            confidence_high: 0,
            confidence_medium: 0,
            confidence_low: 0,
            rounds_run: 0,
            best_round: 0,
            high_conf_rate: 0.0,
            sidecar_path: None,
            retained_path: if options.keep_retained_fasta {
                Some(retained_path.to_string_lossy().to_string())
            } else {
                None
            },
            recheck_attempted: 0,
            recheck_rescued: 0,
            recheck_error: None,
        });
    }

    let round_limit = if options.iter_enable {
        options.iter_max_rounds.max(1)
    } else {
        1
    };

    let mut all_rows: Vec<SidecarRow> = Vec::new();
    let mut rounds: Vec<(usize, Vec<SidecarRow>, PrimerTrimStats)> = Vec::new();
    let mut prev_high_rate: Option<f64> = None;

    for round_idx in 1..=round_limit {
        let relax = round_idx.saturating_sub(1);
        let round_max_mismatch = options.max_mismatch + relax;
        let round_overlap_ratio = (options.min_overlap_ratio - (0.05 * relax as f64)).max(0.5);
        let round_overlap_bp = options
            .min_overlap_bp
            .map(|v| v.saturating_sub(relax).max(8));

        let mut rows: Vec<SidecarRow> = Vec::new();
        for rec in &normalized_records {
            let can = OrientationScore {
                name: "canonical",
                left: best_prefix_match(
                    &rec.seq,
                    forward_primers,
                    round_max_mismatch,
                    options.max_error_rate,
                    round_overlap_bp,
                    round_overlap_ratio,
                    options.end_max_offset,
                ),
                right: best_suffix_match(
                    &rec.seq,
                    &reverse_comp_reverse,
                    round_max_mismatch,
                    options.max_error_rate,
                    round_overlap_bp,
                    round_overlap_ratio,
                    options.end_max_offset,
                ),
            };
            let rev = OrientationScore {
                name: "reverse",
                left: best_prefix_match(
                    &rec.seq,
                    reverse_primers,
                    round_max_mismatch,
                    options.max_error_rate,
                    round_overlap_bp,
                    round_overlap_ratio,
                    options.end_max_offset,
                ),
                right: best_suffix_match(
                    &rec.seq,
                    &reverse_comp_forward,
                    round_max_mismatch,
                    options.max_error_rate,
                    round_overlap_bp,
                    round_overlap_ratio,
                    options.end_max_offset,
                ),
            };
            let (chosen, ambiguous) = choose_orientation(&can, &rev);
            let confidence =
                confidence_label(chosen.matched_ends(), chosen.mismatch_total(), ambiguous);
            rows.push(SidecarRow {
                record_id: rec.header.clone(),
                round: round_idx,
                orientation_chosen: chosen.name.to_string(),
                orientation_ambiguous: usize::from(ambiguous),
                left_hit: usize::from(chosen.left.is_some()),
                right_hit: usize::from(chosen.right.is_some()),
                left_overlap_bp: chosen.left.as_ref().map(|m| m.overlap_bp).unwrap_or(0),
                right_overlap_bp: chosen.right.as_ref().map(|m| m.overlap_bp).unwrap_or(0),
                left_trim_bp: chosen.left.as_ref().map(|m| m.trim_bp).unwrap_or(0),
                right_trim_bp: chosen.right.as_ref().map(|m| m.trim_bp).unwrap_or(0),
                left_mismatch: chosen.left.as_ref().map(|m| m.mismatches).unwrap_or(0),
                right_mismatch: chosen.right.as_ref().map(|m| m.mismatches).unwrap_or(0),
                left_primer_name: chosen
                    .left
                    .as_ref()
                    .map(|m| m.primer.clone())
                    .unwrap_or_default(),
                right_primer_name: chosen
                    .right
                    .as_ref()
                    .map(|m| m.primer.clone())
                    .unwrap_or_default(),
                trim_start: 0,
                trim_end: rec.seq.len(),
                trim_mode: options.trim_mode.clone(),
                confidence,
                dropped_empty: 0,
                recheck_tool: options.recheck_tool.clone(),
                recheck_status: "not_run_phase2".to_string(),
                phylo_mismatch_flag: 0,
            });
        }

        let (_records_preview, mut summary) =
            summarize_sidecar(&mut rows, &seq_map, &options.trim_mode);
        summary.before = before;
        summary.high_conf_rate = if before > 0 {
            summary.confidence_high as f64 / before as f64
        } else {
            0.0
        };
        rounds.push((round_idx, rows.clone(), summary.clone()));
        all_rows.extend(rows);

        if !options.iter_enable {
            break;
        }
        if summary.high_conf_rate >= options.iter_target_conf {
            break;
        }
        if let Some(prev) = prev_high_rate {
            if (summary.high_conf_rate - prev) < options.iter_stop_delta {
                break;
            }
        }
        prev_high_rate = Some(summary.high_conf_rate);
    }

    let (best_round, mut best_rows, _selected_summary) = rounds
        .into_iter()
        .max_by(|a, b| {
            let ka = (
                (a.2.high_conf_rate * 1_000_000.0) as i64,
                a.2.after as i64,
                -(a.2.dropped_empty as i64),
                a.0 as i64,
            );
            let kb = (
                (b.2.high_conf_rate * 1_000_000.0) as i64,
                b.2.after as i64,
                -(b.2.dropped_empty as i64),
                b.0 as i64,
            );
            ka.cmp(&kb)
        })
        .ok_or_else(|| "failed to select best round".to_string())?;

    let mut recheck_error: Option<String> = None;
    let mut recheck_attempted = 0usize;
    let mut recheck_rescued = 0usize;
    if options.recheck_tool == "vsearch" {
        match run_vsearch_recheck(
            &mut best_rows,
            &seq_map,
            forward_primers,
            reverse_primers,
            &reverse_comp_forward,
            &reverse_comp_reverse,
            options.recheck_min_identity,
            options.recheck_min_query_cov,
        ) {
            Ok((attempted, rescued)) => {
                recheck_attempted = attempted;
                recheck_rescued = rescued;
            }
            Err(e) => {
                recheck_error = Some(e);
            }
        }
    }
    let (post_records, best_summary) =
        summarize_sidecar(&mut best_rows, &seq_map, &options.trim_mode);
    write_fasta(fasta_path, &post_records)?;

    let sidecar_path = if options.sidecar_format.to_ascii_lowercase() == "jsonl" {
        let path = PathBuf::from(format!(
            "{}.postprep.primer.jsonl",
            fasta_path.to_string_lossy()
        ));
        let mut f =
            File::create(&path).map_err(|e| format!("failed to write {}: {e}", path.display()))?;
        for row in &best_rows {
            let line = serde_json::json!({
                "record_id": row.record_id,
                "round": row.round,
                "orientation_chosen": row.orientation_chosen,
                "orientation_ambiguous": row.orientation_ambiguous,
                "left_hit": row.left_hit,
                "right_hit": row.right_hit,
                "left_overlap_bp": row.left_overlap_bp,
                "right_overlap_bp": row.right_overlap_bp,
                "left_trim_bp": row.left_trim_bp,
                "right_trim_bp": row.right_trim_bp,
                "left_mismatch": row.left_mismatch,
                "right_mismatch": row.right_mismatch,
                "left_primer_name": row.left_primer_name,
                "right_primer_name": row.right_primer_name,
                "trim_start": row.trim_start,
                "trim_end": row.trim_end,
                "trim_mode": row.trim_mode,
                "confidence": row.confidence,
                "dropped_empty": row.dropped_empty,
                "recheck_tool": row.recheck_tool,
                "recheck_status": row.recheck_status,
                "phylo_mismatch_flag": row.phylo_mismatch_flag
            });
            writeln!(f, "{}", line)
                .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
        }
        Some(path.to_string_lossy().to_string())
    } else {
        let path = PathBuf::from(format!(
            "{}.postprep.primer.tsv",
            fasta_path.to_string_lossy()
        ));
        let mut f =
            File::create(&path).map_err(|e| format!("failed to write {}: {e}", path.display()))?;
        writeln!(f, "record_id\tround\torientation_chosen\torientation_ambiguous\tleft_hit\tright_hit\tleft_overlap_bp\tright_overlap_bp\tleft_trim_bp\tright_trim_bp\tleft_mismatch\tright_mismatch\tleft_primer_name\tright_primer_name\ttrim_start\ttrim_end\ttrim_mode\tconfidence\tdropped_empty\trecheck_tool\trecheck_status\tphylo_mismatch_flag")
            .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
        for row in &best_rows {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.record_id,
                row.round,
                row.orientation_chosen,
                row.orientation_ambiguous,
                row.left_hit,
                row.right_hit,
                row.left_overlap_bp,
                row.right_overlap_bp,
                row.left_trim_bp,
                row.right_trim_bp,
                row.left_mismatch,
                row.right_mismatch,
                row.left_primer_name,
                row.right_primer_name,
                row.trim_start,
                row.trim_end,
                row.trim_mode,
                row.confidence,
                row.dropped_empty,
                row.recheck_tool,
                row.recheck_status,
                row.phylo_mismatch_flag
            )
            .map_err(|e| format!("failed to write {}: {e}", path.display()))?;
        }
        Some(path.to_string_lossy().to_string())
    };

    Ok(PrimerTrimStats {
        before,
        after: best_summary.after,
        removed: best_summary.removed,
        trimmed_both: best_summary.trimmed_both,
        trimmed_left_only: best_summary.trimmed_left_only,
        trimmed_right_only: best_summary.trimmed_right_only,
        untrimmed: best_summary.untrimmed,
        dropped_empty: best_summary.dropped_empty,
        canonical_orientation: best_summary.canonical_orientation,
        reverse_orientation: best_summary.reverse_orientation,
        confidence_high: best_summary.confidence_high,
        confidence_medium: best_summary.confidence_medium,
        confidence_low: best_summary.confidence_low,
        rounds_run: round_limit,
        best_round,
        high_conf_rate: best_summary.high_conf_rate,
        sidecar_path,
        retained_path: if options.keep_retained_fasta {
            Some(retained_path.to_string_lossy().to_string())
        } else {
            None
        },
        recheck_attempted,
        recheck_rescued,
        recheck_error,
    })
}

pub fn load_primer_sets(primer_file: &Path) -> Result<HashMap<String, PrimerSet>, String> {
    let text = fs::read_to_string(primer_file)
        .map_err(|e| format!("failed to read {}: {e}", primer_file.display()))?;
    let parsed: toml::Value = text
        .parse()
        .map_err(|e| format!("failed to parse {}: {e}", primer_file.display()))?;
    let primer_sets = parsed
        .get("primer_sets")
        .and_then(|v| v.as_table())
        .ok_or_else(|| "primer file must contain [primer_sets] table".to_string())?;

    let mut out = HashMap::new();
    for (name, entry) in primer_sets {
        let table = entry
            .as_table()
            .ok_or_else(|| format!("primer_sets.{name} must be a table"))?;
        let forward = table
            .get("forward")
            .and_then(|v| v.as_array())
            .ok_or_else(|| format!("primer_sets.{name}.forward must be an array"))?
            .iter()
            .map(|v| {
                v.as_str()
                    .map(|s| s.trim().to_ascii_uppercase().replace('U', "T"))
                    .ok_or_else(|| format!("primer_sets.{name}.forward must contain strings"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let reverse = table
            .get("reverse")
            .and_then(|v| v.as_array())
            .ok_or_else(|| format!("primer_sets.{name}.reverse must be an array"))?
            .iter()
            .map(|v| {
                v.as_str()
                    .map(|s| s.trim().to_ascii_uppercase().replace('U', "T"))
                    .ok_or_else(|| format!("primer_sets.{name}.reverse must contain strings"))
            })
            .collect::<Result<Vec<_>, _>>()?;

        out.insert(name.clone(), PrimerSet { forward, reverse });
    }

    Ok(out)
}

pub fn combine_primer_sets(
    all_sets: &HashMap<String, PrimerSet>,
    selected: &[String],
) -> Result<(Vec<String>, Vec<String>), String> {
    let mut forward = Vec::new();
    let mut reverse = Vec::new();
    for set_name in selected {
        let set = all_sets
            .get(set_name)
            .ok_or_else(|| format!("primer set '{set_name}' was not found"))?;
        for p in &set.forward {
            if !forward.contains(p) {
                forward.push(p.clone());
            }
        }
        for p in &set.reverse {
            if !reverse.contains(p) {
                reverse.push(p.clone());
            }
        }
    }
    Ok((forward, reverse))
}

#[derive(Debug, Clone)]
struct HeaderExtractor {
    pattern: Regex,
    captures_organism: bool,
}

#[derive(Debug, Clone)]
struct DuplicateEntry {
    acc_id: String,
    accession: String,
    organism_name: String,
    header: String,
}

fn csv_escape(value: &str) -> String {
    let must_quote =
        value.contains(',') || value.contains('"') || value.contains('\n') || value.contains('\r');
    if !must_quote {
        return value.to_string();
    }
    format!("\"{}\"", value.replace('"', "\"\""))
}

fn template_has_field(template: &str, field: &str) -> bool {
    template.contains(&format!("{{{field}}}"))
}

fn compile_header_extractors(header_formats: &[String]) -> (Vec<HeaderExtractor>, bool, bool) {
    let mut extractors = Vec::new();
    let mut has_acc_id_template = false;
    let mut has_organism_template = false;
    let field_re = Regex::new(r"\{([a-zA-Z0-9_]+)\}").expect("field regex");

    let mut uniq = header_formats.to_vec();
    uniq.sort();
    uniq.dedup();

    for template in uniq {
        let has_acc_id = template_has_field(&template, "acc_id");
        let organism_field = if template_has_field(&template, "organism_raw") {
            Some("organism_raw")
        } else if template_has_field(&template, "organism") {
            Some("organism")
        } else {
            None
        };
        has_acc_id_template |= has_acc_id;
        has_organism_template |= organism_field.is_some();
        if !has_acc_id {
            continue;
        }

        let mut parts = Vec::new();
        let mut last = 0usize;
        let mut seen_acc = false;
        let mut seen_org = false;
        for caps in field_re.captures_iter(&template) {
            let m = caps.get(0).expect("whole capture");
            let field = caps.get(1).map(|mm| mm.as_str()).unwrap_or("");
            if m.start() > last {
                parts.push(regex::escape(&template[last..m.start()]));
            }
            if field == "acc_id" {
                if !seen_acc {
                    parts.push("(?P<acc_id>[A-Za-z0-9._-]+)".to_string());
                    seen_acc = true;
                } else {
                    parts.push("[A-Za-z0-9._-]+".to_string());
                }
            } else if organism_field == Some(field) {
                if !seen_org {
                    parts.push("(?P<organism_name>.*?)".to_string());
                    seen_org = true;
                } else {
                    parts.push(".*?".to_string());
                }
            } else {
                parts.push(".*?".to_string());
            }
            last = m.end();
        }
        if last < template.len() {
            parts.push(regex::escape(&template[last..]));
        }
        if seen_acc {
            let pat = Regex::new(&format!("^{}$", parts.join(""))).expect("header extractor regex");
            extractors.push(HeaderExtractor {
                pattern: pat,
                captures_organism: seen_org,
            });
        }
    }

    (extractors, has_acc_id_template, has_organism_template)
}

fn extract_header_fields(
    header: &str,
    extractors: &[HeaderExtractor],
) -> (Option<String>, Option<String>) {
    for ex in extractors {
        let Some(caps) = ex.pattern.captures(header) else {
            continue;
        };
        let Some(acc) = caps.name("acc_id").map(|m| m.as_str().to_string()) else {
            continue;
        };
        let org = caps
            .name("organism_name")
            .map(|m| m.as_str().trim().to_string());
        return (Some(acc), org.filter(|s| !s.is_empty()));
    }
    (None, None)
}

pub fn write_duplicate_acc_reports_csv(
    fasta_path: &Path,
    header_formats: &[String],
) -> Result<
    (
        Option<PathBuf>,
        Option<PathBuf>,
        Option<DuplicateAccStats>,
        Option<String>,
    ),
    String,
> {
    let (extractors, has_acc_id_template, has_organism_template) =
        compile_header_extractors(header_formats);
    if !has_acc_id_template {
        return Ok((
            None,
            None,
            None,
            Some("selected header format does not include {acc_id}".to_string()),
        ));
    }
    if !has_organism_template {
        return Ok((
            None,
            None,
            None,
            Some(
                "selected header format does not include {organism_raw} or {organism}".to_string(),
            ),
        ));
    }
    if extractors.is_empty() {
        return Ok((
            None,
            None,
            None,
            Some("could not compile header extractor".to_string()),
        ));
    }
    if !extractors.iter().any(|ex| ex.captures_organism) {
        return Ok((
            None,
            None,
            None,
            Some(
                "selected header format does not include {organism_raw} or {organism} alongside {acc_id}".to_string(),
            ),
        ));
    }

    let dup_re = Regex::new(r"_dup\d+$").expect("dup suffix regex");
    let records = read_fasta(fasta_path)?;
    let mut seq_groups: HashMap<String, Vec<DuplicateEntry>> = HashMap::new();
    let mut total_records = 0usize;
    let mut parsed_records = 0usize;
    let mut unparsed_records = 0usize;
    for rec in records {
        total_records += 1;
        let header = rec.header.trim().to_string();
        let (acc_id, organism_name) = extract_header_fields(&header, &extractors);
        let (Some(acc_id), Some(organism_name)) = (acc_id, organism_name) else {
            unparsed_records += 1;
            continue;
        };
        parsed_records += 1;
        let accession = dup_re.replace_all(&acc_id, "").to_string();
        seq_groups
            .entry(rec.seq.to_ascii_uppercase())
            .or_default()
            .push(DuplicateEntry {
                acc_id,
                accession,
                organism_name,
                header,
            });
    }

    let records_path = PathBuf::from(format!(
        "{}.duplicate_acc.records.csv",
        fasta_path.to_string_lossy()
    ));
    let groups_path = PathBuf::from(format!(
        "{}.duplicate_acc.groups.csv",
        fasta_path.to_string_lossy()
    ));
    let mut rec_out = File::create(&records_path)
        .map_err(|e| format!("failed to write {}: {e}", records_path.display()))?;
    let mut grp_out = File::create(&groups_path)
        .map_err(|e| format!("failed to write {}: {e}", groups_path.display()))?;

    writeln!(rec_out, "group_id,sequence_hash,sequence_length,records_in_group,unique_accessions,unique_organisms,cross_organism_duplicate,acc_id,accession,organism_name,header")
        .map_err(|e| format!("failed to write {}: {e}", records_path.display()))?;
    writeln!(grp_out, "group_id,sequence_hash,sequence_length,records_in_group,unique_accessions,unique_organisms,cross_organism_duplicate,accessions,organism_names")
        .map_err(|e| format!("failed to write {}: {e}", groups_path.display()))?;

    let mut groups: Vec<(String, Vec<DuplicateEntry>)> = seq_groups.into_iter().collect();
    groups.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let mut group_id = 0usize;
    let mut duplicate_groups = 0usize;
    let mut duplicate_records = 0usize;
    let mut cross_organism_groups = 0usize;
    for (seq, entries) in groups {
        if entries.len() < 2 {
            continue;
        }
        let unique_accessions: HashSet<String> =
            entries.iter().map(|e| e.accession.clone()).collect();
        if unique_accessions.len() < 2 {
            continue;
        }
        let unique_organisms: HashSet<String> =
            entries.iter().map(|e| e.organism_name.clone()).collect();
        let is_cross_organism = unique_organisms.len() > 1;

        group_id += 1;
        duplicate_groups += 1;
        duplicate_records += entries.len();
        if is_cross_organism {
            cross_organism_groups += 1;
        }

        let mut hasher = DefaultHasher::new();
        seq.hash(&mut hasher);
        let seq_hash = format!("{:016x}", hasher.finish());
        let mut accessions_sorted: Vec<String> = unique_accessions.into_iter().collect();
        accessions_sorted.sort();
        let mut organisms_sorted: Vec<String> = unique_organisms.into_iter().collect();
        organisms_sorted.sort();
        let records_in_group = entries.len();
        writeln!(
            grp_out,
            "{},{},{},{},{},{},{},{},{}",
            group_id,
            seq_hash,
            seq.len(),
            records_in_group,
            accessions_sorted.len(),
            organisms_sorted.len(),
            if is_cross_organism { "true" } else { "false" },
            csv_escape(&accessions_sorted.join(";")),
            csv_escape(&organisms_sorted.join(";"))
        )
        .map_err(|e| format!("failed to write {}: {e}", groups_path.display()))?;

        for item in entries {
            writeln!(
                rec_out,
                "{},{},{},{},{},{},{},{},{},{},{}",
                group_id,
                seq_hash,
                seq.len(),
                records_in_group,
                accessions_sorted.len(),
                organisms_sorted.len(),
                if is_cross_organism { "true" } else { "false" },
                csv_escape(&item.acc_id),
                csv_escape(&item.accession),
                csv_escape(&item.organism_name),
                csv_escape(&item.header)
            )
            .map_err(|e| format!("failed to write {}: {e}", records_path.display()))?;
        }
    }

    Ok((
        Some(records_path),
        Some(groups_path),
        Some(DuplicateAccStats {
            total_records,
            parsed_records,
            unparsed_records,
            duplicate_groups,
            duplicate_records,
            cross_organism_groups,
        }),
        None,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_path(name: &str) -> PathBuf {
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        std::env::temp_dir().join(format!("taxondb-{name}-{ts}.fa"))
    }

    #[test]
    fn length_filter_keeps_expected_records() {
        let path = tmp_path("length-filter");
        fs::write(&path, ">a\nAAAA\n>b\nAAA\n>c\nAA\n").expect("write test fasta");

        let stats = apply_length_filter(&path, Some(3), Some(4)).expect("filter");
        assert_eq!(stats.before, 3);
        assert_eq!(stats.after, 2);
        assert_eq!(stats.removed, 1);

        let out = fs::read_to_string(&path).expect("read result");
        assert!(out.contains(">a"));
        assert!(out.contains(">b"));
        assert!(!out.contains(">c"));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn primer_trim_trims_both_ends() {
        let path = tmp_path("primer-trim");
        fs::write(&path, ">x\nAAATTTTAAA\n").expect("write test fasta");

        let forward = vec!["AAA".to_string()];
        let reverse = vec!["TTT".to_string()];
        let opts = PrimerTrimOptions::default();
        let stats = apply_primer_trim(&path, &forward, &reverse, &opts).expect("trim");
        assert_eq!(stats.before, 1);
        assert_eq!(stats.after, 1);
        assert_eq!(stats.trimmed_both, 1);

        let out = fs::read_to_string(&path).expect("read result");
        assert!(out.contains("\nTTTT\n"));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn duplicate_acc_report_emits_csv() {
        let path = tmp_path("duplicate-report");
        fs::write(
            &path,
            ">gb|ACC1|Species A\nACGT\n>gb|ACC2|Species B\nACGT\n>gb|ACC3|Species C\nTTTT\n",
        )
        .expect("write test fasta");
        let formats = vec!["gb|{acc_id}|{organism}".to_string()];
        let (records_csv, groups_csv, stats, reason) =
            write_duplicate_acc_reports_csv(&path, &formats).expect("duplicate report");
        assert!(reason.is_none());
        assert!(records_csv.as_ref().is_some_and(|p| p.exists()));
        assert!(groups_csv.as_ref().is_some_and(|p| p.exists()));
        let st = stats.expect("stats");
        assert_eq!(st.total_records, 3);
        assert_eq!(st.duplicate_groups, 1);
        assert_eq!(st.cross_organism_groups, 1);

        if let Some(p) = records_csv {
            let _ = fs::remove_file(p);
        }
        if let Some(p) = groups_csv {
            let _ = fs::remove_file(p);
        }
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn primer_trim_end_offset_recovers_internal_primers() {
        let path = tmp_path("primer-offset");
        fs::write(&path, ">x\nNNAAACCCCAAANN\n").expect("write test fasta");

        let forward = vec!["AAA".to_string()];
        let reverse = vec!["TTT".to_string()];

        let mut strict = PrimerTrimOptions::default();
        strict.keep_retained_fasta = false;
        let strict_stats =
            apply_primer_trim(&path, &forward, &reverse, &strict).expect("strict trim");
        assert_eq!(strict_stats.trimmed_both, 0);
        let strict_out = fs::read_to_string(&path).expect("read strict output");
        assert!(strict_out.contains("\nNNAAACCCCAAANN\n"));

        fs::write(&path, ">x\nNNAAACCCCAAANN\n").expect("rewrite test fasta");
        let mut relaxed = PrimerTrimOptions::default();
        relaxed.keep_retained_fasta = false;
        relaxed.end_max_offset = 2;
        let relaxed_stats =
            apply_primer_trim(&path, &forward, &reverse, &relaxed).expect("offset trim");
        assert_eq!(relaxed_stats.trimmed_both, 1);
        let relaxed_out = fs::read_to_string(&path).expect("read relaxed output");
        assert!(relaxed_out.contains("\nCCCC\n"));

        let _ = fs::remove_file(&path);
    }
}
