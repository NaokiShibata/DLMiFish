# Changelog

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
