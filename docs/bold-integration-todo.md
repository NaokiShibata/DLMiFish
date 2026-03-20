# BOLD 統合実装 TODO / 進捗

`docs/Plan.md` を実装に落とすための TODO と、現時点の進捗整理。

ステータス:
- `[x]` 完了
- `[ ]` 未着手
- `部分` 実装は入ったが、Plan の完了条件までは未到達

## 現在地

- Phase 1 は概ね完了。
- Phase 2 は概ね完了。
- Phase 3 は CLI への BOLD 経路追加まで進んだが、実 API を使った end-to-end 検証と大規模取得対応は未完。
- Phase 4 は strict merge、`source_merge.csv`、deterministic merge 出力まで実装したが、BOLD API 実データでの確認と大規模取得時の最終メモリ戦略は未完。
- Phase 5 は sidecar 拡張と README / sample config 更新まで進んだ。GUI / tests は未着手。

## Plan.md との整合確認

`docs/Plan.md` の Phase ごとの実装順序に対して、現状は次の通り。

- Phase 1: `BuildSource`, `ResolvedTaxon`, `build --source`, source-aware `load_config()`, source-aware `setup_entrez()` を実装済み。
- Phase 2: `extract_ncbi_records_from_genbank_chunk()` と `emit_records_to_fasta()` への責務分離を実装済み。
- Phase 3: `taxondb_bold.py` を追加し、taxon query -> preprocess -> query -> download -> marker 判定 -> `CanonicalRecord` 変換の流れを実装済み。
- Phase 4: NCBI accession set による `insdcacs` strict 判定、`source_merge.csv`、source ごとの spool を使った deterministic merge 出力を実装済み。
- Phase 5: `.acc_organism.csv` の列拡張、README / sample config 更新までは実装済み。GUI 非破壊確認、post-prep 条件の最終整理、テスト計画の実行は未着手。

結論:
- 実装順序は `docs/Plan.md` に概ね沿って進められている。
- ただし Plan の「完了条件」を満たしたとはまだ言えない。
- 特に `--source both` の実 API 検証、BOLD only / both の回帰確認、GUI 影響確認が不足している。

## 実装済みの主な変更

- [`taxondbbuilder.py`](/home/naoki/ghq/github.com/NaokiShibata/TaxonDBBuilder/taxondbbuilder.py)
  - `BuildSource`, `ResolvedTaxon`, `CanonicalRecord` を追加
  - `build --source [ncbi|bold|both]` を追加
  - `load_config()` / `normalize_marker_map()` / `setup_entrez()` を source-aware 化
  - NCBI 抽出を `extract_ncbi_records_from_genbank_chunk()` と `emit_records_to_fasta()` に分離
  - BOLD row -> `CanonicalRecord` 変換
  - `both` の strict merge (`linked_by_insdcacs`) を実装
  - `.fasta.source_merge.csv` を追加
  - `.fasta.acc_organism.csv` に source 系列の列を追加
- [`taxondb_bold.py`](/home/naoki/ghq/github.com/NaokiShibata/TaxonDBBuilder/taxondb_bold.py)
  - BOLD API クライアントを新規追加
  - preprocessor / summary / query / documents download
  - marker 判定、`insdcacs` token 化、BOLD record 正規化

## Phase 1: Source-aware 土台

- [x] `BuildSource` を定義し、`build` コマンドに `--source [ncbi|bold|both]` を追加した。
- [x] `load_config()` を source-aware にし、`source=bold` では `[ncbi]` を必須にしない形へ変更した。
- [x] `[bold]` セクションの読込方針を決め、`taxondb_bold.py` 側で timeout / retry / user-agent / base URL のデフォルトを持たせた。
- [x] `normalize_marker_map()` を source-aware にし、BOLD 用 `bold.marker_codes` を解釈するようにした。
- [x] `ResolvedTaxon` を導入し、BOLD 利用時は scientific name を取得する流れに差し替えた。
- [x] `setup_entrez()` を source 条件付き warning に変更した。
- [x] `--dump-gb` / `--from-gb` / `--resume` の validation を source-aware にした。

注意:
- `source=bold` でも taxid 入力から scientific name を得るために Entrez を使う。
- つまり「`[ncbi]` は不要」は満たしているが、「NCBI ネットワークアクセス完全不要」ではない。

## Phase 2: NCBI 出力責務の分離

- [x] `CanonicalRecord` を導入し、NCBI 抽出結果を source 共通表現へ移した。
- [x] `process_genbank_chunk()` 相当の責務を `extract_ncbi_records_from_genbank_chunk()` に分離した。
- [x] `emit_records_to_fasta()` を追加し、header 生成・FASTA 書込・sidecar 行生成を切り出した。
- [x] NCBI duplicate 判定 (`acc_to_seqs`, `_dupN`) を新しい流れへ移植した。
- [x] 既存 worker/threading の中で record 抽出と書込を分けた。
- [x] `.acc_organism.csv` 用行構造を source / processid / marker_key を持てる形へ拡張した。
- [ ] `--source ncbi` の回帰比較を実データで確認していない。

積み残し:
- FASTA 件数、ログ主要行、`.acc_organism.csv`、post-prep の回帰確認をまだ行っていない。

## Phase 3: BOLD 単体取得

- [x] BOLD API クライアントを [`taxondb_bold.py`](/home/naoki/ghq/github.com/NaokiShibata/TaxonDBBuilder/taxondb_bold.py) に分離した。
- [x] BOLD 通信を標準ライブラリで実装し、timeout・HTTP error handling・gzip・retry・query 上限超過を入れた。
- [ ] TSV streaming 切替可能な I/O 境界までは未実装。現状は JSON 前提。
- [x] BOLD record から `processid`, `sampleid`, `taxon_name`, `marker_code`, `insdcacs`, `nucleotides` を抽出して正規化した。
- [x] marker 判定を `marker_codes` 優先、`aliases` / `phrases` / marker key fallback の順で実装した。
- [x] BOLD 用 `acc_id` を `BOLD_<processid or source_record_id>` で namespaced した。
- `部分` `--source bold` で FASTA 生成経路は実装したが、実 API を使った end-to-end 検証は未実施。

積み残し:
- BOLD API の実レスポンス差異による field 名の揺れを実データで確認していない。
- `query 0 件`、`malformed document`、`insdcacs` 欠損 / 複数値ケースの実地確認が未実施。
- 大規模件数時のメモリ挙動を確認していない。

## Phase 4: Both 統合

- [x] NCBI kept record から accession set を構築し、BOLD `insdcacs` との strict match だけで suppress する実装を追加した。
- [x] `insdcacs` の token 化を実装し、区切り文字混在にある程度対応した。
- [x] `linked_to_ncbi`, `emitted_to_fasta`, `skip_reason` を `CanonicalRecord` と sidecar 行へ反映した。
- [x] NCBI kept / BOLD kept を source ごとの spool に積み、`taxon_name`, `marker_key`, `source_record_id` の安定順で最終 FASTA を出すようにした。
- [x] `source_merge.csv` を新設し、keep / skip の両方を出力するようにした。
- `部分` temp spool による最終出力の安定化は実装したが、BOLD download 自体は現状 JSON 全読み込みのため、大規模取得時の最終メモリ戦略は未完。
- [ ] strict rule 以外の根拠を使わないことの回帰テストは未作成。

注意:
- `--source both` で `--from-gb` は NCBI 側にのみ適用する形へ修正済み。

## Phase 5: 仕上げと周辺整備

- [x] `write_acc_organism_mapping_csv()` を列拡張し、`source`, `source_record_id`, `processid`, `sampleid`, `marker_key`, `linked_to_ncbi`, `emitted_to_fasta`, `skip_reason` を追加した。
- `部分` `write_duplicate_acc_reports_csv()` 自体の source-aware 再設計は未着手だが、`source=both` では duplicate report をデフォルト自動実行しない条件へ修正済み。
- [x] `README.md`, `configs/db.toml`, `configs/markers_mitogenome.toml` のサンプルを更新し、`--source`, `[bold]`, `[markers.<id>.bold]`, `source_merge.csv` を追記した。
- [ ] `docs/Plan.md` との差分を補助ドキュメントへ反映する作業は今回の更新まで。ユーザー向けドキュメント更新は未着手。
- [ ] GUI 非破壊確認は未着手。
- [ ] Python CLI と Rust GUI runner の仕様差分整理は未着手。
- [ ] テスト観点の文書化と実行は未着手。

## 積み残しの作業

優先度高:
- [ ] `--source bold` の実 API での end-to-end 実行確認
- [ ] `--source both` の実 API での end-to-end 実行確認
- [ ] `--source ncbi` の回帰確認
- [ ] `source=bold` / `source=both` の duplicate_report / post-prep 条件をさらに詰める

優先度中:
- [ ] 0件・timeout・query 上限超過・malformed document のテストケース整理

優先度中〜低:
- [ ] BOLD download を TSV streaming へ切替可能にして、大規模取得時のメモリ戦略を仕上げる
- [ ] GUI 側を当面 `ncbi` 固定とするか、Rust runner に BOLD を入れるか方針決定
- [ ] `tauri-gui/src-tauri/src/taxondb_runner.rs` との仕様差分を解消

## 現時点での判定

`docs/Plan.md` に従った実装は進められているか:
- はい。Phase 1 -> Phase 4 までは順序通りに着手できている。

`docs/Plan.md` の完了条件を満たしているか:
- いいえ。まだ未満。

未達の理由:
- `python3 taxondbbuilder.py build -c configs/db.toml -t 117570 -m coi --source both` を実 API 付きで確認していない。
- `insdcacs == accession` strict suppress の実地確認が未実施。
- sequence 一致だけでは suppress しないことの回帰テストが未実施。
- GUI 非破壊確認と実 API 検証が未完。

## 検証方法

以下はユーザー側で実施する想定の確認手順。

前提:
- 依存を入れる
  - `uv venv --python 3.11`
  - `source .venv/bin/activate`
  - `uv pip install -r requirements.txt`
- `configs/db.toml` の `ncbi.email` / `ncbi.api_key` を必要に応じて設定する
- 出力先は空のディレクトリを推奨

### 1. NCBI only 回帰

実行例:
```bash
python3 taxondbbuilder.py build \
  -c configs/db.toml \
  -t 117570 \
  -m 12s \
  --source ncbi \
  --dump-gb Results/gb_ncbi
```

確認:
- FASTA が生成される
- `.fasta.log` に query / total / matched / kept が出る
- `.fasta.acc_organism.csv` が生成される
- `source_merge.csv` が生成され、`source=ncbi` 行のみになる

### 2. BOLD only

実行例:
```bash
python3 taxondbbuilder.py build \
  -c configs/db.toml \
  -t "Salmo salar" \
  -m coi \
  --source bold
```

確認:
- FASTA が生成される
- header の `acc_id` が `BOLD_` で始まる
- `.fasta.source_merge.csv` が生成される
- `.fasta.acc_organism.csv` の `source` 列に `bold` が入る

### 3. Both

実行例:
```bash
python3 taxondbbuilder.py build \
  -c configs/db.toml \
  -t "Salmo salar" \
  -m coi \
  --source both \
  --dump-gb Results/gb_both
```

確認:
- FASTA が生成される
- `.fasta.source_merge.csv` に `source=ncbi` と `source=bold` の両方が出る
- `skip_reason=linked_by_insdcacs` の行があれば、その行は `emitted_to_fasta=false` になっている
- `source=bold` で `skip_reason` 空の行は FASTA に出ている

### 4. Both + from-gb

事前に NCBI の GB を保存した後で:
```bash
python3 taxondbbuilder.py build \
  -c configs/db.toml \
  -t "Salmo salar" \
  -m coi \
  --source both \
  --from-gb Results/gb_both
```

確認:
- NCBI 側は `from-gb` で進み、BOLD 側は API 取得される
- `.log` に `# query count taxid=...: from-gb` が出る

### 5. duplicate_report の条件

`source=both` ではデフォルトでは duplicate report を自動実行しない。

確認:
```bash
python3 taxondbbuilder.py build \
  -c configs/db.toml \
  -t "Salmo salar" \
  -m coi \
  --source both \
  --post-prep
```

- `.log` に `duplicate_acc_report: skipped (step disabled)` が出ること

明示的に有効化:
```bash
python3 taxondbbuilder.py build \
  -c configs/db.toml \
  -t "Salmo salar" \
  -m coi \
  --source both \
  --post-prep \
  --post-prep-step duplicate_report
```

- duplicate report の CSV が生成されるか、もしくは header 条件不足で skip 理由が `.log` に出ること

## 実装時の注意を継続

- 既存ワークツリーにはユーザーの未コミット変更があるため、それを壊さない前提で進める。
- GUI は Rust 側にロジック重複があるため、CLI 側の完了と GUI 側の完了は別管理にする。
