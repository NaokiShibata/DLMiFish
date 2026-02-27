use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

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
}

#[derive(Debug, Clone)]
pub struct PrimerSet {
    pub forward: Vec<String>,
    pub reverse: Vec<String>,
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

fn primer_matches_prefix(seq: &str, primer: &str) -> bool {
    if seq.len() < primer.len() {
        return false;
    }
    for (s, p) in seq.chars().zip(primer.chars()) {
        let Some(values) = iupac_values(p) else {
            return false;
        };
        if !values.contains(s) {
            return false;
        }
    }
    true
}

fn primer_matches_suffix(seq: &str, primer: &str) -> bool {
    if seq.len() < primer.len() {
        return false;
    }
    let start = seq.len() - primer.len();
    primer_matches_prefix(&seq[start..], primer)
}

fn max_prefix_match(seq: &str, primers: &[String]) -> usize {
    primers
        .iter()
        .filter(|p| primer_matches_prefix(seq, p))
        .map(|p| p.len())
        .max()
        .unwrap_or(0)
}

fn max_suffix_match(seq: &str, primers: &[String]) -> usize {
    primers
        .iter()
        .filter(|p| primer_matches_suffix(seq, p))
        .map(|p| p.len())
        .max()
        .unwrap_or(0)
}

pub fn apply_primer_trim(
    fasta_path: &Path,
    forward_primers: &[String],
    reverse_primers: &[String],
) -> Result<PrimerTrimStats, String> {
    let records = read_fasta(fasta_path)?;
    let before = records.len();
    let mut after_records = Vec::new();

    let mut reverse_comp_forward = Vec::with_capacity(forward_primers.len());
    for p in forward_primers {
        reverse_comp_forward.push(reverse_complement_iupac(p)?);
    }
    let mut reverse_comp_reverse = Vec::with_capacity(reverse_primers.len());
    for p in reverse_primers {
        reverse_comp_reverse.push(reverse_complement_iupac(p)?);
    }

    let mut stats = PrimerTrimStats {
        before,
        after: 0,
        removed: 0,
        trimmed_both: 0,
        trimmed_left_only: 0,
        trimmed_right_only: 0,
        untrimmed: 0,
        dropped_empty: 0,
        canonical_orientation: 0,
        reverse_orientation: 0,
    };

    for mut rec in records {
        rec.seq = rec.seq.to_ascii_uppercase().replace('U', "T");
        let seq_len = rec.seq.len();

        let can_left = max_prefix_match(&rec.seq, forward_primers);
        let can_right = max_suffix_match(&rec.seq, &reverse_comp_reverse);
        let rev_left = max_prefix_match(&rec.seq, reverse_primers);
        let rev_right = max_suffix_match(&rec.seq, &reverse_comp_forward);

        let can_ends = usize::from(can_left > 0) + usize::from(can_right > 0);
        let rev_ends = usize::from(rev_left > 0) + usize::from(rev_right > 0);
        let can_trim = can_left + can_right;
        let rev_trim = rev_left + rev_right;

        let (left, right, canonical) = if (can_ends, can_trim) >= (rev_ends, rev_trim) {
            (can_left, can_right, true)
        } else {
            (rev_left, rev_right, false)
        };

        let ends_trimmed = usize::from(left > 0) + usize::from(right > 0);
        if ends_trimmed == 0 {
            stats.untrimmed += 1;
            after_records.push(rec);
            continue;
        }

        if ends_trimmed == 2 {
            stats.trimmed_both += 1;
        } else if left > 0 {
            stats.trimmed_left_only += 1;
        } else {
            stats.trimmed_right_only += 1;
        }

        if canonical {
            stats.canonical_orientation += 1;
        } else {
            stats.reverse_orientation += 1;
        }

        let right_idx = if right > 0 {
            seq_len.saturating_sub(right)
        } else {
            seq_len
        };
        if left >= right_idx {
            stats.dropped_empty += 1;
            continue;
        }
        rec.seq = rec.seq[left..right_idx].to_string();
        after_records.push(rec);
    }

    stats.after = after_records.len();
    stats.removed = stats.before.saturating_sub(stats.after);
    write_fasta(fasta_path, &after_records)?;
    Ok(stats)
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

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
        let stats = apply_primer_trim(&path, &forward, &reverse).expect("trim");
        assert_eq!(stats.before, 1);
        assert_eq!(stats.after, 1);
        assert_eq!(stats.trimmed_both, 1);

        let out = fs::read_to_string(&path).expect("read result");
        assert!(out.contains("\nTTTT\n"));
        let _ = fs::remove_file(&path);
    }
}
