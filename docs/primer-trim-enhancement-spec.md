# Primer Trim Enhancement Spec (v1)

## 1. Goal

`post_prep.primer_trim` の推定精度を上げるため、以下を導入する。

- 末端完全一致のみから、ミスマッチ許容・部分一致閾値つき判定へ拡張
- 不確実配列に対する再判定ステップ (BLAST/vsearch)
- MSA 用にプライマー領域を保持した出力を残す
- CLI (Python) / GUI (Rust) の挙動を一致させる
- 反復実行で改善可能な設計 (iteration loop) を提供する

## 2. Non-goal

- 種同定自体の正誤判定は本仕様の主目的ではない
- 系統樹不一致判定は「配列除外の最終判定」ではなく「低信頼フラグ付け」に使う

## 3. Processing Overview

`primer_trim` は以下 6 フェーズで動作する。

1. Primer candidate build
2. Endpoint fuzzy matching (left/right, canonical/reverse)
3. Confidence scoring and trim decision
4. Optional recheck using external aligner (BLAST/vsearch)
5. Optional phylogeny mismatch flagging (MSA + tree)
6. Output and sidecar report emit

## 4. New Config Keys

`[post_prep]` に以下を追加する。

```toml
[post_prep]
primer_file = "configs/primers.toml"
primer_set = ["mifish_u", "mifish_ev2"]

# matching
primer_max_mismatch = 1         # integer >= 0
primer_max_error_rate = 0.12    # 0.0 - 1.0
primer_min_overlap_bp = 14      # integer >= 1
primer_min_overlap_ratio = 0.7  # 0.0 - 1.0
primer_end_max_offset = 30      # integer >= 0, 末端から内側の探索許容量(bp)

# trim decision
primer_trim_mode = "one_or_both"   # both_required | one_or_both | one_end_mark_only
primer_keep_retained_fasta = true  # MSA向けプライマー保持配列を残す

# iteration
primer_iter_enable = true
primer_iter_max_rounds = 3
primer_iter_stop_delta = 0.002   # 改善率閾値 (0.2%)
primer_iter_target_conf = 0.98   # high confidence率の目標

# optional external recheck
primer_recheck_tool = "vsearch"  # off | vsearch | blast (vsearch first)
primer_recheck_min_identity = 0.85
primer_recheck_min_query_cov = 0.7

# optional phylogeny check
primer_phylo_check = "flag_only" # off | flag_only
primer_phylo_target_confidence = "medium" # low | medium
primer_phylo_min_taxa = 20
primer_phylo_max_samples = 2000

# sidecar
primer_sidecar_format = "tsv"   # tsv | jsonl
```

### Validation rules

- `primer_max_mismatch` は 0 以上
- `primer_max_error_rate` は 0 以上 1 以下
- `primer_min_overlap_bp` と `primer_min_overlap_ratio` は両方指定可能
- overlap 閾値は `max(bp条件, ratio条件)` を満たした場合に成立
- `primer_end_max_offset` は 0 以上。0 の場合は従来どおり完全な末端アンカー
- `primer_trim_mode=both_required` のとき片端一致は切らずに sidecar に記録

## 5. Matching and Scoring

## 5.1 Candidate generation

- canonical: left=`forward`, right=`reverse_rc`
- reverse: left=`reverse`, right=`forward_rc`
- IUPAC を展開して比較する

## 5.2 Endpoint fuzzy matching

各 endpoint で以下を計算する。

- `overlap_bp`
- `mismatch_count`
- `error_rate = mismatch_count / overlap_bp`
- `full_len_match` (primer 全長到達)

一致条件は次をすべて満たすこと。

- `overlap_bp >= required_overlap_bp`
- `mismatch_count <= primer_max_mismatch`
- `error_rate <= primer_max_error_rate`

## 5.3 Orientation selection

向き選択は以下の優先順で比較する。

1. 一致 endpoint 数
2. 一致 endpoint の合計スコア (`overlap_bp - mismatch_count`)
3. 合計 mismatch の少なさ
4. 同点時 canonical

## 5.4 Confidence

各配列に `high | medium | low` を付与する。

- `high`: 両端一致かつ mismatch 条件が十分余裕あり
- `medium`: 片端一致または閾値近傍
- `low`: 一致なし、または canonical/reverse の競合が強い

閾値は固定値ではなく config により計算される。

## 6. Trim Decision Policy

- `both_required`
  - 両端一致時のみ hard-trim を実施
  - 片端一致は untrimmed として残す
- `one_or_both`
  - 現行互換
  - 片端一致でも一致端のみ切る
- `one_end_mark_only`
  - 一致があっても配列自体は切らない
  - sidecar の endpoint annotation のみ付与

## 7. Iteration Loop

`primer_iter_enable=true` の場合、以下を繰り返す。

1. round N を実行して sidecar を作成
2. `high_conf_rate`, `both_end_rate`, `low_conf_count`, `trimmed_fraction` を算出
3. 不確実配列 (`low`, `medium` の一部) を再判定
4. round N+1 の候補集合と閾値補正を確定

停止条件は次のいずれか。

- `N >= primer_iter_max_rounds`
- `high_conf_rate` 改善が `primer_iter_stop_delta` 未満
- `high_conf_rate >= primer_iter_target_conf`

重要: 反復で primer 配列自体を破壊しないよう、判定入力は常に retained 配列を使う。

## 8. External Tool Integration

外部ツールは任意。未インストール時は warning を出して fallback する。

## 8.1 vsearch mode

- 用途: 不確実配列の endpoint 再判定
- 想定: `--usearch_global` で primer DB へ照合
- 出力: top hit の `id`, `qcov`, `strand`, `qstart/qend`

## 8.2 blast mode

- 用途: short primer の曖昧一致再評価
- 想定: `blastn-short` とローカル primer DB
- 出力: top hit の `pident`, `qcovs`, `qstart/qend`, `sstart/send`

## 8.3 seqkit utilities

- 前処理・集計補助に使用可能
- 本体ロジックを seqkit に依存させず、補助ツールとして扱う

## 9. MSA and Phylogeny Mismatch Check

`primer_phylo_check=flag_only` 時の挙動。

- 対象は `low confidence` と `orientation ambiguous` を優先
- retained FASTA を入力に MSA を構築
- 系統樹と種名の整合性を簡易評価して `phylo_mismatch_flag` を sidecar に付与
- フラグは除外判定に直結させず、再判定優先度に使う

## 10. Outputs

同一 basename で 3 系統の成果物を出力する。

- hard-trim FASTA
- retained FASTA (primer 領域保持、MSA 用)
- sidecar TSV/JSONL

sidecar の最低限カラム。

- `record_id`
- `orientation_chosen`
- `left_hit` / `right_hit`
- `left_overlap_bp` / `right_overlap_bp`
- `left_trim_bp` / `right_trim_bp`
- `left_mismatch` / `right_mismatch`
- `left_primer_name` / `right_primer_name`
- `trim_start` / `trim_end`
- `trim_mode`
- `confidence`
- `round`
- `recheck_tool`
- `recheck_status`
- `phylo_mismatch_flag`

## 11. Logging

既存ログに加えて以下を出力する。

- `# post_prep primer trim round=N ...`
- `# post_prep primer confidence: high=... medium=... low=...`
- `# post_prep primer recheck: tool=... attempted=... rescued=...`
- `# post_prep primer phylo_flag: total=... flagged=...`

GUI パーサは round ごとの集計を読めるように更新する。

## 12. CLI/Rust Unification Requirements

- 同一入力 FASTA + 同一 config + 同一 primer sets で統計値一致
- 一致対象:
  - `before/after/removed`
  - `trimmed_both/left_only/right_only/untrimmed`
  - `confidence histogram`
  - `orientation histogram`
- Python の基準テストベクトルを Rust 側でも共有する

## 13. Backward Compatibility

- 既定値で現行挙動に近づける:
  - `primer_max_mismatch=0`
  - `primer_min_overlap_bp=primer_len`
  - `primer_trim_mode=one_or_both`
  - `primer_recheck_tool=off`
  - `primer_iter_enable=false`
- 新機能は opt-in とする

## 14. Implementation Phases

1. Phase 1: fuzzy match + trim mode + sidecar
2. Phase 2: iteration loop + confidence
3. Phase 3: vsearch/blast recheck
4. Phase 4: phylo mismatch flag
5. Phase 5: CLI/Rust parity test in CI

## 15. Test Plan (minimum)

- 単体:
  - IUPAC + mismatch + overlap の境界値
  - orientation tie-break
  - trim mode 3 種
- 結合:
  - round 反復停止条件
  - recheck fallback
  - retained/hard-trim/sidecar の整合性
- 回帰:
  - 現行設定で結果が大きく変わらないこと

## 16. Fixed Decisions

1. `primer_trim_mode` 既定値は `one_or_both`
2. sidecar 標準フォーマットは `TSV`
3. recheck 第一候補は `vsearch`
4. phylo mismatch の適用対象は `medium` まで
5. 初期実装範囲は `Phase 2` まで
