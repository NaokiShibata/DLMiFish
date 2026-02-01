# TaxonDBBuilder (Generic DB Builder)
NCBIから任意の分類群・任意のマーカーの配列を取得し、DB用のFASTAを生成するツールです。分類群は **taxid / 学名** のどちらでも指定でき、マーカーは **TOMLで定義したフレーズ群** をprefix指定で呼び出せます。

本リポジトリは「分類群特異の処理」から「汎用DB FASTA生成」へ方針変更しています。解析 (抽出・フィルタリング・分類付与など)は対象外で、**DB用FASTAの生成が目的**です。

## 動作環境
- Python 3.8+ (3.11+ 推奨)
- 主要依存: biopython, rich, typer
- Python < 3.11 の場合: tomli
- パッケージ管理: uv (推奨)

### 環境構築 (uv)
任意のフォルダにuvをインストールする例です (例: `~/tools/uv`)。

```bash
# uvを任意のフォルダにインストール
mkdir -p ~/tools/uv
curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR="~/tools/uv" sh

# PATHへ追加 (bashの場合)
echo 'export PATH="$HOME/tools/uv:$PATH"' >> ~/.bashrc
source ~/.bashrc

# zshの場合は ~/.zshrc を更新してください

# 動作確認
uv --version
```

```bash
# Pythonのインストール (必要な場合)
uv python install 3.11

# 仮想環境
uv venv --python 3.11
source .venv/bin/activate

# 依存導入
uv pip install -r requirements.txt
```

### 補足
- uvはPythonパッケージのみを扱います。現時点で外部ツールは不要です。

## クイックスタート
環境構築から **GenBankキャッシュ付きの実行** までの最短手順です。

```bash
# 1) 仮想環境と依存導入
uv python install 3.11
uv venv --python 3.11
source .venv/bin/activate
uv pip install -r requirements.txt

# 2) 設定（api_key / email を入力）
# configs/db.toml を編集

# 3) 実行（GenBankキャッシュ保存）
python3 taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --dump-gb Results/gb
```

再実行時にキャッシュを優先して使う場合:
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --dump-gb Results/gb --resume
```

## 設定ファイル (TOML)
`configs/db.toml` を編集して使います。

```toml
[ncbi]
email = "your.email@example.com"
api_key = "YOUR_API_KEY"

db = "nucleotide"
rettype = "gb"
retmode = "text"
per_query = 100
use_history = true

markers_file = "configs/markers_mitogenome.toml"

[output]
default_header_format = "{acc_id}|{organism}|{marker}|{label}|{type}|{loc}|{strand}"

[output.header_formats]
simple = "{acc_id}|{marker}|{loc}"
verbose = "{acc_id}|{organism_raw}|{marker_raw}|{label_raw}|{type_raw}|{loc}|{strand}"

[taxon]
noexp = false

# markers は外部ファイルで定義します (例: configs/markers_mitogenome.toml)。

[filters]
# フィルタ無しがデフォルト。必要な場合だけ指定してください。
# organelle = "mitochondrion"
# length_min = 120
# length_max = 30000
# source = "ddbj_embl_genbank"
# date_from = "1990/01/01"
# date_to = "2025/12/31"
# include_keywords = ["12S"]
# exclude_keywords = ["WGS"]
# extra = "complete[prop]"
```

### 検索・抽出の考え方
- `phrases`: **検索用**の簡易フレーズ。自動的に `"..."[All Fields]` として扱われます。
- `terms`: **検索用**の生クエリ。`rrnS[Gene]` のようにフィールド指定をそのまま書けます。
- `region_patterns`: **GenBank feature抽出用**の正規表現。`gene/product/note/standard_name` などの注釈に対してマッチします。

`region_patterns` 未指定の場合は、`phrases/terms` から**リテラル**として自動生成します。
抽出はGenBankのfeature注釈に依存するため、目的の領域が出ない場合は `region_patterns` と `feature_types/feature_fields` を調整してください。

### markers_file について
- `markers_file` でマーカー定義を外部TOMLへ分離できます。
- 外部定義は `[markers]` テーブルを持つ必要があります。

## 設定ガイド
このツールの設定は「検索 (NCBIクエリ)」と「抽出 (GenBank feature注釈)」を分けて考えると整理しやすいです。

### 設定ファイルの役割
- `db.toml` は「**共通設定** (NCBI / 出力 / filters / markers_file)」を持ちます。
- `markers_file` は「**マーカー定義**」だけを持ちます ([markers] テーブル)。
- これにより、用途ごとにマーカー定義を差し替える運用ができます。

### 1. 最小構成 (必須)
- `ncbi` セクション: `email` / `api_key` / `db` / `rettype` など
- `markers_file`: マーカー定義ファイルへのパス
- `output`: FASTAヘッダーの形式

### 2. マーカー定義の考え方
マーカー定義は以下の3要素で構成されます。

- **検索用 (phrases / terms)**
  NCBIの検索に使う語句。
  `phrases` は `"..."[All Fields]` に自動変換されます。
  `terms` は `rrnS[Gene]` のようにフィールド指定をそのまま書けます。

- **抽出用 (region_patterns)**
  GenBankのfeature注釈 (gene/product/note/standard_name)に対する正規表現。

- **対象feature (feature_types / feature_fields)**
  `feature_types` で対象のfeature型を限定できます (rRNA/gene/CDS など)。
  `feature_fields` で参照する注釈項目を制御します。

### 3. よくあるパターン
**rRNA系 (12S/16S など)**
```toml
[markers."12s"]
aliases = ["12", "12s"]
phrases = ["12S", "rrnS", "small subunit ribosomal RNA"]
region_patterns = ["12S", "rrnS", "small subunit ribosomal RNA"]
feature_types = ["rRNA", "gene"]
feature_fields = ["gene", "product", "note", "standard_name"]
```

**タンパク質コーディング (COI/ND1 など)**
```toml
[markers."coi"]
aliases = ["coi", "co1", "cox1"]
phrases = ["COI", "CO1", "COX1", "cytochrome c oxidase subunit I"]
region_patterns = ["COI", "CO1", "COX1", "cytochrome c oxidase subunit I"]
feature_types = ["CDS", "gene"]
feature_fields = ["gene", "product", "note", "standard_name"]
```

### 4. FASTAヘッダーの指定
- `[output.header_formats]` にテンプレートを定義
- `markers.<id>.header_format` でテンプレート名を選択
 (直接テンプレート文字列を書いてもOK)

> [!NOTE]
> `markers_mitogenome.toml`にデフォルトで設定している`header_format=mifish`はPMiFishパイプラインとMiFishパイプラインのDBに対応するフォーマットになっています。

例:
```toml
[output.header_formats]
simple = "{acc_id}|{marker}|{loc}"
mifish = "gb|{acc_id}|{organism}"

[markers."12s"]
header_format = "mifish"
```

## 逆引き (よくある目的別)
### Q. 目的のマーカー情報が登録されていない
1) `markers_file` (例: `configs/markers_mitogenome.toml`)に新規マーカーを追加
2) `aliases` を付けてCLIから指定できるようにする
3) `phrases/terms` を検索用に、`region_patterns` を抽出用に設定

最小例:
```toml
[markers."mygene"]
aliases = ["mygene", "mg"]
phrases = ["MyGene", "my gene product"]
region_patterns = ["MyGene", "my gene product"]
feature_types = ["gene", "CDS"]
feature_fields = ["gene", "product", "note", "standard_name"]
```

### Q. 検索はヒットするが抽出されない
- `region_patterns` がfeature注釈に合っていない可能性があります。
  → GenBankの該当レコードを確認し、`gene/product/note` の表記に合わせて調整してください。

### Q. ヒット数が多すぎる / 少なすぎる
- `terms` を使ってフィールド指定 (例: `rrnS[Gene]`)すると検索精度が上がります。
- 必要に応じて `[filters]` を追加して絞り込みます。

### Q. 12S/16S 以外 (ITS/18S など)を使いたい
- `markers_file` に追加でOKです。
  例: `ITS`, `18S`, `28S` は `feature_types = ["rRNA", "gene"]` で定義するケースが多いです。

## 使い方
### コマンドと設定ファイルの対応
以下の対応関係を押さえると、設定とコマンドが繋がって理解できます。

| コマンド引数 | 参照する設定 | 説明 |
| --- | --- | --- |
| `-c/--config` | `db.toml` | 共通設定 (NCBI/出力/filters/markers_file) |
| `-m/--marker` | `markers_file` 内の `[markers.<id>]` | 使うマーカー定義 (aliases で選択) |
| `-t/--taxon` | なし | taxid/学名を指定 (学名はTaxonomyで解決) |
| `--workers` | なし | 抽出処理の並列数 |
| `--out` | なし | 出力先 (省略時は `Results/db/YYYYMMDD/`) |
| `--output-prefix` | なし | 出力FASTAファイル名のプレフィックス (default: `taxondbbuilder_`) |
| `--dump-gb` | なし | GenBankチャンクを保存 (キャッシュ) |
| `--from-gb` | なし | 保存済みGenBankチャンクから抽出 |
| `--resume` | なし | キャッシュを優先して利用 |

### 具体例 (設定とコマンドの対応)
1) `markers_file` に `12s` を定義しておく
2) `-m 12s` を指定 → `markers."12s"` の設定が適用される
3) `phrases/terms` で検索し、`region_patterns` で抽出する

設定例 (抜粋):
```toml
[markers."12s"]
aliases = ["12", "12s"]
phrases = ["12S", "rrnS"]
region_patterns = ["12S", "rrnS"]
```

実行例:
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s
```

ファイル名にプレフィックスを付けたい場合:
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --output-prefix "mifish"
```

GenBankを保存しつつ実行 (acc_idごとに `.gb` を保存):
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --dump-gb Results/gb
```

保存済みGenBankから再抽出:
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --from-gb Results/gb
```

中断後の再開 (キャッシュ利用):
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --dump-gb Results/gb --resume
```

キャッシュは `Results/gb/.cache/` に保存されます。

### taxid指定
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 32443 -m 12s
```

### 学名指定
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t "Salmo salar" -m 12s
```

### 複数 taxon / marker
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 32443 -t 7777 -m 12 -m coi
```

### 並列抽出 (ダウンロードと変換の並列化)
```bash
python3 taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --workers 2
```

## 抽出ロジック
- NCBIから**GenBank形式**で取得し、feature注釈から目的領域を抽出します。
- `region_patterns` が `gene/product/note/standard_name` などの注釈に一致したfeatureのみFASTAへ出力します。
- `feature_types` を指定すると対象feature型を限定できます (例: rRNA, gene, CDS)。

## FASTAヘッダーフォーマット
- `[output.header_formats]` にテンプレートを定義し、`markers.<id>.header_format` で選択できます。
- 直接テンプレート文字列を `header_format` に書くことも可能です。

使用できるプレースホルダ:
`{acc}`, `{acc_id}`, `{organism}`, `{organism_raw}`, `{marker}`, `{marker_raw}`, `{label}`, `{label_raw}`, `{type}`, `{type_raw}`, `{start}`, `{end}`, `{loc}`, `{strand}`, `{dup}`

## 出力
- 出力先: `Results/db/YYYYMMDD/`
- ファイル名: `taxid{ID}__{marker}.fasta` (複数指定時は `multi_taxon` / `multi_marker`)
- 実行ログ: 出力FASTAと同名の `.log`

## キャッシュと再抽出
- `--dump-gb` で **acc_idごとのGenBankファイル** を保存します。
- キャッシュは `--dump-gb` 配下の `.cache/` に保存されます。
- `--resume` は **キャッシュがある場合にそれを優先**して使います（出力は毎回新規に作り直します）。
- `--from-gb` はネットワークを使わず、保存済みGenBankから抽出のみを実行します。

## 重複の扱い
- **同一アクセッション + 同一配列**: 重複として除外
- **同一アクセッション + 異なる配列**: 両方残し、IDに `_dupN` を付与して警告ログに記録

## フィルタについて
フィルタは**指定なしを受け付けます**。必要な場合のみTOMLの `[filters]` に追加してください。

## Legacy
旧パイプラインは削除済みです。現行方針では**DB用FASTA生成のみ**を対象とします。
