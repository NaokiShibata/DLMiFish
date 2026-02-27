# Runner Bundle Spec (v3)

## Goal

Run the build pipeline directly in Rust from the Tauri process (no sidecar / no python runtime for GUI path).

## Path layout

- No external runner binary/script is required.

## Build command

From `tauri-gui/`:

```bash
npm run tauri:build
```

## Implementation details

- `start_run` triggers integrated Rust runner (`taxondb_runner` module).
- NCBI fetch (`esearch`/`efetch`) + GenBank parse + FASTA emit are handled in Rust.
- Post-prep `primer_trim` / `length_filter` is executed in Rust after build (duplicate report is pending migration).
