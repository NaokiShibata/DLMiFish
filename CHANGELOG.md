# Changelog

## 3.0.0 (planned)

### Scope
- `main` (`2effd40`) から `feat/tauri-gui-v1` を Pull Request でマージする差分を対象。

### Added
- Tauri ベースのクロスプラットフォーム GUI を追加 (`tauri-gui/`)。
- Sidecar 実行基盤 (Rust) を追加し、`taxondbbuilder` CLI のジョブ実行・進捗表示・ログ表示・結果表示を GUI から操作可能に。
- `db.toml` 読み込み機能を追加し、GUI の入力パラメータへ反映可能に。
- GitHub Actions に GUI ビルドワークフロー (`.github/workflows/tauri-build.yml`) を追加。

### Changed
- Entrez 取得処理で HTTP 400 発生時のフォールバックを追加し、履歴ベース取得失敗時の取得継続性を改善。
- `.gitignore` を更新し、Tauri GUI のローカル生成物 (`node_modules`, `dist`, `target`, sidecar バイナリ, `__pycache__`) を除外。
- `tauri-gui/.gitignore` を追加し、GUI 開発・ビルドで生成されるローカルファイルを `tauri-gui` 配下で管理可能に。
- `tauri-gui/src-tauri/src/main.rs` の sidecar 解決処理を改善し、`TAXONDBBUILDER_BIN` がディレクトリを指す誤設定を明示エラー化。
- `README.md` と `tauri-gui/README.md` を更新し、`uv` 前提の実行環境構築手順と GUI/CLI 利用フローを整理。

## 1.0.0 (2026-02-01)

### Changed
- 旧β版 (分類群特異パイプライン) から、**汎用DB FASTA生成ツール**へ方針転換。
- 実行は `python3 taxondbbuilder.py build` を使用。
- 設定は `setting.txt` ではなく **TOML (`configs/db.toml`) ** を使用。
- 分類群は taxid / 学名を指定可能。マーカーはTOMLで定義したフレーズ群をprefix指定。
- 出力は **DB用の結合FASTA** とログのみ (解析・抽出・フィルタリングは対象外) 。
- 旧スクリプト群 (`scripts/`) と関連フォルダ (`GI_Folder/`, `ID/`, `MiFishDB/`) は削除。

### Usage (new)
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 32443 -m 12s
```

## 0.x (beta)
### Summary
- DLMiFishとして公開されていたβ版パイプライン。
- Teleostei/Chondrichthyes/Cyclostomataの**12S rRNA**に特化。
- `python3 DLMiFish.py` で一括実行し、複数スクリプトの連携で結果を生成。
- 設定は `setting.txt` を使用。
