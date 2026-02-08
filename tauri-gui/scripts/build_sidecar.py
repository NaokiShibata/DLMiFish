#!/usr/bin/env python3
from __future__ import annotations

import argparse
import platform
import shutil
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str], cwd: Path) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd), check=True)


def detect_os_key() -> str:
    name = platform.system().lower()
    if "linux" in name:
        return "linux"
    if "darwin" in name:
        return "macos"
    if "windows" in name:
        return "windows"
    raise RuntimeError(f"unsupported platform: {platform.system()}")


def detect_target_triple(repo_root: Path) -> str:
    out = subprocess.check_output(["rustc", "-vV"], cwd=str(repo_root), text=True)
    for line in out.splitlines():
        if line.startswith("host:"):
            return line.split(":", 1)[1].strip()
    raise RuntimeError("failed to detect rust target triple from `rustc -vV`")


def write_stub_binary(path: Path, os_key: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if os_key == "windows":
        content = "@echo off\r\necho taxondbbuilder sidecar stub. Build real binary with scripts/build_sidecar.py\r\nexit /b 1\r\n"
        path.write_text(content, encoding="utf-8")
    else:
        content = "#!/usr/bin/env sh\n" \
            "echo 'taxondbbuilder sidecar stub. Build real binary with scripts/build_sidecar.py' >&2\n" \
            "exit 1\n"
        path.write_text(content, encoding="utf-8")
        path.chmod(0o755)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build TaxonDBBuilder sidecar with PyInstaller")
    parser.add_argument("--repo-root", default="..", help="Repository root path")
    parser.add_argument(
        "--tauri-root",
        default=".",
        help="tauri-gui root path (contains src-tauri/)"
    )
    parser.add_argument(
        "--stub",
        action="store_true",
        help="Create target-triple sidecar stub file only (offline/cargo-check aid).",
    )
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    tauri_root = Path(args.tauri_root).resolve()
    script_path = repo_root / "taxondbbuilder.py"

    if not script_path.exists():
        raise FileNotFoundError(f"taxondbbuilder.py not found: {script_path}")

    os_key = detect_os_key()
    target_triple = detect_target_triple(repo_root)
    bin_name = "taxondbbuilder.exe" if os_key == "windows" else "taxondbbuilder"
    bundled_name = (
        f"taxondbbuilder-{target_triple}.exe"
        if os_key == "windows"
        else f"taxondbbuilder-{target_triple}"
    )
    target_dir = tauri_root / "src-tauri" / "bin"
    target_dir.mkdir(parents=True, exist_ok=True)
    dest_bin = target_dir / bundled_name

    if args.stub:
        write_stub_binary(dest_bin, os_key)
        print(f"created sidecar stub: {dest_bin}")
        return 0

    try:
        subprocess.run(
            [sys.executable, "-m", "PyInstaller", "--version"],
            cwd=str(repo_root),
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "PyInstaller is not available. Install it first, or use --stub for offline cargo check."
        ) from exc

    run(
        [
            sys.executable,
            "-m",
            "PyInstaller",
            "--onefile",
            "--name",
            "taxondbbuilder",
            "taxondbbuilder.py",
        ],
        cwd=repo_root,
    )

    dist_bin = repo_root / "dist" / bin_name
    if not dist_bin.exists():
        raise FileNotFoundError(f"built sidecar not found: {dist_bin}")

    shutil.copy2(dist_bin, dest_bin)
    print(f"built sidecar: {dest_bin}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
