# TaxonDBBuilder GUI (Tauri)

## Overview
- No Python environment required (standalone sidecar binary required)
- No server required (local execution only)
- Configure taxid/marker/filters/post_prep in GUI and run `taxondbbuilder build`

## Screens
1. Run Setup
2. Run Monitor
3. Results

## Setup
```bash
cd tauri-gui
npm install
```

## Linux build dependencies (Ubuntu/Debian)
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

## Run in dev mode
```bash
npm run tauri:dev
```

## Sidecar binaries
Place sidecar binary under `src-tauri/bin/`.

Tauri expects target-triple suffixed names for bundled sidecars:
- Linux/macOS: `src-tauri/bin/taxondbbuilder-<target-triple>`
- Windows: `src-tauri/bin/taxondbbuilder-<target-triple>.exe`

For development, you can point directly to the sidecar via `TAXONDBBUILDER_BIN`.

```bash
export TAXONDBBUILDER_BIN=/abs/path/to/taxondbbuilder
npm run tauri:dev
```

## Build sidecar with PyInstaller
From `tauri-gui/`:

```bash
python scripts/build_sidecar.py --repo-root .. --tauri-root .
```

This generates a standalone `taxondbbuilder` and places it under
`src-tauri/bin/` with the expected target-triple suffix.

Note: `cargo check` also validates `externalBin`, so run the sidecar build once
before `cargo check` if the file does not exist yet.

If you are offline and only need `cargo check`, create a stub sidecar:

```bash
python scripts/build_sidecar.py --repo-root .. --tauri-root . --stub
```

## Persistent config
- `~/.taxondb_gui/config.json`
- If `save api_key` is off, API key is not persisted
