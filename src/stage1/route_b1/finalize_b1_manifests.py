#!/usr/bin/env python3
"""Regenerate B1 manifests after publication using relocatable internal paths."""

from __future__ import annotations

import argparse
import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


RUN_ID = "stage1_B1_complete_current_universe_20260904"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def records(path: Path) -> int:
    if path.suffix in {".tsv", ".txt"}:
        with path.open(errors="replace") as handle:
            lines = sum(1 for _ in handle)
        return max(0, lines - 1) if path.suffix == ".tsv" else lines
    if path.suffix in {".fasta", ".fa"}:
        with path.open(errors="replace") as handle:
            return sum(1 for line in handle if line.startswith(">"))
    return 1


def atomic_frame(frame: pd.DataFrame, path: Path) -> None:
    partial = path.with_name(path.name + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    frame.to_csv(partial, sep="\t", index=False)
    os.replace(partial, path)


def archive_existing(path: Path) -> Path:
    archived = path.with_name(path.stem + ".prepublish_paths" + path.suffix)
    if archived.exists():
        raise FileExistsError(f"Refusing to overwrite archived manifest: {archived}")
    path.rename(archived)
    return archived


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, required=True)
    parser.add_argument("--current-tsv", type=Path, required=True)
    parser.add_argument("--current-fasta", type=Path, required=True)
    parser.add_argument("--previous-tsv", type=Path, required=True)
    parser.add_argument("--new69-tsv", type=Path, required=True)
    parser.add_argument("--previous-corrected-combined", type=Path, required=True)
    parser.add_argument("--historical-legacy-candidates", type=Path, required=True)
    parser.add_argument("--nls-weights", type=Path, required=True)
    parser.add_argument("--esm-weights", type=Path, required=True)
    args = parser.parse_args()

    audit = args.release_root / "03_audit"
    old_main = archive_existing(audit / "b1_manifest.tsv")
    old_input = archive_existing(audit / "input_manifest.tsv")
    timestamp = datetime.now(timezone.utc).isoformat()

    external = [
        ("INPUT", args.current_tsv, "current Stage0 B non-anchor table"),
        ("INPUT", args.current_fasta, "current Stage0 B non-anchor FASTA"),
        ("INPUT_PROVENANCE", args.previous_tsv, "validated 14,556 component universe"),
        ("INPUT", args.new69_tsv, "69 newly eligible Stage0 records"),
        ("MODEL_SCORE_PROVENANCE", args.previous_corrected_combined, "validated 14,556 corrected scores"),
        ("HISTORICAL_AUDIT", args.historical_legacy_candidates, "archived historical legacy candidates"),
        ("MODEL_WEIGHT", args.nls_weights, "NLSExplorer classifier weights"),
        ("MODEL_WEIGHT", args.esm_weights, "ESM-1b encoder weights"),
    ]
    rows: list[dict] = []
    for role, path, description in external:
        if not path.exists():
            raise FileNotFoundError(path)
        rows.append(
            {
                "artifact_role": role,
                "path": str(path),
                "records": records(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "description": description,
                "run_id": RUN_ID,
                "generated_at_utc": timestamp,
            }
        )

    excluded = {
        old_main.resolve(),
        old_input.resolve(),
        (audit / "b1_manifest.tsv").resolve(),
        (audit / "input_manifest.tsv").resolve(),
    }
    for path in sorted(args.release_root.rglob("*")):
        if not path.is_file() or path.resolve() in excluded or path.name.endswith(".partial"):
            continue
        relative = path.relative_to(args.release_root)
        role = "OUTPUT"
        if relative.parts[0] in {"03_audit", "04_inference_runs"}:
            role = "AUDIT"
        elif relative.parts[0] == "scripts":
            role = "CODE"
        elif relative.parts[0] == "00_input":
            role = "DERIVED_INPUT"
        elif relative.name == "README.md":
            role = "DOCUMENTATION"
        rows.append(
            {
                "artifact_role": role,
                "path": str(relative),
                "records": records(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "description": "published B1 artifact",
                "run_id": RUN_ID,
                "generated_at_utc": timestamp,
            }
        )
    manifest = pd.DataFrame(rows).sort_values(
        ["artifact_role", "path"], kind="mergesort"
    ).reset_index(drop=True)
    atomic_frame(manifest, audit / "b1_manifest.tsv")

    inputs = manifest.loc[
        manifest["artifact_role"].isin(
            ["INPUT", "INPUT_PROVENANCE", "MODEL_SCORE_PROVENANCE", "MODEL_WEIGHT", "DERIVED_INPUT"]
        )
    ].copy()
    atomic_frame(inputs, audit / "input_manifest.tsv")
    print(f"manifest_records\t{len(manifest)}")
    print(f"input_manifest_records\t{len(inputs)}")


if __name__ == "__main__":
    main()
