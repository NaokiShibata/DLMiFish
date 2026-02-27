# scientific name ↔ tax_id 対応表 作成手順

本ドキュメントでは、NCBI Taxonomy (`taxdump`) から GUI 用の `taxonomy.db` を生成する手順を示します。

生成物:
- 中間生成物: `taxid_scientific_name.csv`
- 本番用 DB: `taxonomy.db` (SQLite / テーブル名: `taxonomy`)

GUI (Tauri) は `tauri-gui/resources/taxonomy.db` を参照します。

## 1. データソース
- NCBI Taxonomy taxdump (`taxdump.tar.gz`)
- 使用ファイル:
  - `names.dmp` (名称一覧)
  - `nodes.dmp` (階層情報) ※本手順では未使用

この手順では `names.dmp` から `name_class == "scientific name"` の行のみを使います。

## 2. 全体フロー
1. `taxdump.tar.gz` を取得
2. 展開して `names.dmp` を用意
3. `names.dmp` から `taxid_scientific_name.csv` を生成
4. CSV から `taxonomy.db` を生成
5. `tauri-gui/resources/taxonomy.db` に配置して GUI をビルド

## 3. taxdump の取得と展開
URL は変更される可能性があるため、最新の配布先は NCBI を確認してください。

```bash
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz
```

展開後に `names.dmp` があることを確認してください。

## 4. 生成スクリプト
リポジトリに `tools/build_taxonomy_db.py` を追加しています。

### 4.1 names.dmp -> CSV
```bash
python3 tools/build_taxonomy_db.py extract-csv \
  --names-dmp ./names.dmp \
  --out-csv ./taxid_scientific_name.csv
```

### 4.2 CSV -> SQLite
```bash
python3 tools/build_taxonomy_db.py build-sqlite \
  --in-csv ./taxid_scientific_name.csv \
  --out-db ./taxonomy.db
```

### 4.3 一括実行 (推奨)
```bash
python3 tools/build_taxonomy_db.py all \
  --names-dmp ./names.dmp \
  --out-csv ./taxid_scientific_name.csv \
  --out-db ./tauri-gui/resources/taxonomy.db
```

## 5. 生成される SQLite スキーマ
```sql
CREATE TABLE taxonomy (
    tax_id INTEGER PRIMARY KEY,
    scientific_name TEXT NOT NULL
);

CREATE INDEX idx_taxonomy_scientific_name_nocase
    ON taxonomy(scientific_name COLLATE NOCASE);
```

## 6. GUI 側の利用前提
- GUI は `tauri-gui/resources/taxonomy.db` を固定参照します。
- テーブル名は `taxonomy`、列は `tax_id` / `scientific_name` を前提に検索します。

## 7. 更新ポリシー (推奨)
- `taxdump` は定期更新されるため、`taxonomy.db` も再生成して差し替えてください。
- 目安:
  - 研究用途: 月1回〜数ヶ月に1回
  - 一般用途: 半年に1回程度

運用フロー:
1. 最新 `taxdump` を取得
2. `names.dmp` を展開
3. `tools/build_taxonomy_db.py all ...` を実行
4. `tauri-gui/resources/taxonomy.db` を更新
5. GUI を再ビルド / 再配布
