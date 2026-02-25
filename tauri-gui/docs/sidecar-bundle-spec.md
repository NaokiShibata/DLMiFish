# Sidecar Bundle Spec (v1)

## Goal

Bundle standalone `taxondbbuilder` executable per OS and launch it from Tauri.

## Path layout

- `src-tauri/bin/taxondbbuilder-<target-triple>`
- `src-tauri/bin/taxondbbuilder-<target-triple>.exe` (Windows)

## Build command

From `tauri-gui/`:

```bash
python scripts/build_sidecar.py --repo-root .. --tauri-root .
```

Offline check-only mode:

```bash
python scripts/build_sidecar.py --repo-root .. --tauri-root . --stub
```

## Implementation details

- Uses PyInstaller (`--onefile`) on each runner OS
- Detects host target triple from `rustc -vV`
- Copies produced binary from `../dist/` into `src-tauri/bin/` with
  Tauri-required target-triple suffix
- `tauri.conf.json` references these binaries via `bundle.externalBin`
