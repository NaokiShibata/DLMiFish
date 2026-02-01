# TaxonDBBuilder (Generic DB Builder)
NCBIから任意の分類群・任意のマーカーの配列を取得し、DB用のFASTAを生成するツールです。分類群は **taxid / 学名** のどちらでも指定でき、マーカーは **TOMLで定義したフレーズ群** をprefix指定で呼び出せます。

本リポジトリは「分類群特異の処理」から「汎用DB FASTA生成」へ方針変更しています。解析（抽出・フィルタリング・分類付与など）は対象外で、**DB用FASTAの生成が目的**です。

## 動作環境
- Python 3.8+ (3.11+ 推奨)
- 主要依存: biopython, typer
- Python < 3.11 の場合: tomli
- パッケージ管理: uv (推奨)

### 環境構築 (uv)
任意のフォルダにuvをインストールする例です（例: `~/tools/uv`）。

```bash
# uvを任意のフォルダにインストール
mkdir -p ~/tools/uv
UV_INSTALL_DIR=~/tools/uv curl -LsSf https://astral.sh/uv/install.sh | sh

# PATHへ追加（bashの場合）
echo 'export PATH="$HOME/tools/uv:$PATH"' >> ~/.bashrc
source ~/.bashrc

# zshの場合は ~/.zshrc を更新してください

# 動作確認
uv --version
```

```bash
# Pythonのインストール（必要な場合）
uv python install 3.11

# 仮想環境
uv venv --python 3.11
source .venv/bin/activate

# 依存導入
uv pip install -r requirements.txt
```

### 補足
- uvはPythonパッケージのみを扱います。現時点で外部ツールは不要です。

## 設定ファイル (TOML)
`configs/db.toml` を編集して使います。

```toml
[ncbi]
email = "your.email@example.com"
api_key = "YOUR_API_KEY"

db = "nucleotide"
rettype = "fasta"
retmode = "text"
per_query = 100
use_history = true

[taxon]
noexp = false

[markers."12s"]
aliases = ["12", "12s"]
phrases = ["12S", "rrnS", "small subunit ribosomal RNA"]

[markers."coi"]
aliases = ["coi", "co1"]
phrases = ["COI", "CO1", "cytochrome c oxidase subunit I"]

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

## 使い方
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

## 出力
- 出力先: `Results/db/YYYYMMDD/`
- ファイル名: `taxid{ID}__{marker}.fasta` (複数指定時は `multi_taxon` / `multi_marker`)
- 実行ログ: 出力FASTAと同名の `.log`

## 重複の扱い
- **同一アクセッション + 同一配列**: 重複として除外
- **同一アクセッション + 異なる配列**: 両方残し、IDに `_dupN` を付与して警告ログに記録

## フィルタについて
フィルタは**指定なしを受け付けます**。必要な場合のみTOMLの `[filters]` に追加してください。

## Legacy
旧パイプラインは削除済みです。現行方針では**DB用FASTA生成のみ**を対象とします。
