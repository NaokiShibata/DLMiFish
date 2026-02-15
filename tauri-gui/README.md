# TaxonDBBuilder GUI (Tauri)

## 概要
- GUIアプリをビルドして実行できます。
- sidecarバイナリをコマンドラインで実行できます。
- サーバーは不要です（ローカル実行のみ）。

## 1. GUIアプリとして使う（推奨）
### 前提
- Node.js / npm
- Rust / cargo
- Linux の場合は Tauri 依存ライブラリ

Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install -y \
  pkg-config \
  libwayland-dev \
  libxkbcommon-dev \
  libwebkit2gtk-4.1-dev \
  libgtk-3-dev \
  libayatana-appindicator3-dev \
  librsvg2-dev \
  patchelf
```

### ビルド手順
リポジトリルート (`TaxonDBBuilder/`) で:
```bash
uv python install 3.11
uv venv --python 3.11
source .venv/bin/activate
uv pip install -r requirements.txt
uv pip install pyinstaller
```

`tauri-gui/` で:
```bash
cd tauri-gui
npm install
python3 scripts/build_sidecar.py --repo-root .. --tauri-root .
npm run tauri build
```

### 実行
生成物を起動します（例: Linux AppImage）。
- `src-tauri/target/release/bundle/`

## 2. コマンドラインとして使う
`tauri-gui/src-tauri/bin/` の sidecar を直接実行できます。

```bash
./src-tauri/bin/taxondbbuilder-<target-triple> --help
./src-tauri/bin/taxondbbuilder-<target-triple> list-markers -c ../configs/db.toml
./src-tauri/bin/taxondbbuilder-<target-triple> build -c ../configs/db.toml -t 117570 -m 12s --dry-run
```

注意:
- sidecar はCLIです。単体でGUIは起動しません。

## 補足（開発時のみ）
- 開発起動: `npm run tauri:dev`
- `TAXONDBBUILDER_BIN` を使う場合は、ディレクトリではなく実行ファイルの絶対パスを指定してください。

```bash
export TAXONDBBUILDER_BIN=/abs/path/to/taxondbbuilder
npm run tauri:dev
```

## 設定保存先
- `~/.taxondb_gui/config.json`
- `save api_key` をオフにすると API key は保存されません。
