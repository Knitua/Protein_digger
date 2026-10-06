#!/usr/bin/env python3
"""Rewrite staged manifest paths to released paths and verify every artifact."""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path


ROOT = Path("/root/autodl-tmp/Agent_analysis_v2")
MANIFEST = ROOT / "03_stage0_AB_partition/03_audit/stage0_manifest.tsv"
BACKUP = ROOT / "03_stage0_AB_partition/03_audit/stage0_manifest_pre_release_paths.tsv"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def data_records(path: Path) -> int:
    if path.suffix != ".tsv":
        return -1
    with path.open("rb") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def released_path(value: str) -> str:
    mappings = (
        (".stage0_build_20260903/stage0/", "03_stage0_AB_partition/"),
        (".stage0_build_20260903/stage1_A/", "03_stage1_A/"),
        (".stage0_build_20260903/stage1_B_common/", "03_stage1_B/00_common_input/"),
    )
    for old, new in mappings:
        if value.startswith(old):
            return new + value[len(old) :]
    if value == "03_stage1_A/route_A_candidate_classification_3487.tsv":
        return (
            "03_stage1_A/archive/pre_stage0_isoform_aware_20260903/"
            "route_A_candidate_classification_3487.tsv"
        )
    return value


def main() -> None:
    if BACKUP.exists():
        raise SystemExit(f"backup already exists: {BACKUP}")
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise SystemExit("manifest has no header")
        rows = list(reader)

    for row in rows:
        row["path"] = released_path(row["path"])
        path = ROOT / row["path"]
        if not path.is_file():
            raise SystemExit(f"missing released artifact: {path}")
        if int(row["bytes"]) != path.stat().st_size:
            raise SystemExit(f"size mismatch: {path}")
        if row["sha256"] != digest(path):
            raise SystemExit(f"sha256 mismatch: {path}")
        if row["records"] and path.suffix == ".tsv":
            if int(row["records"]) != data_records(path):
                raise SystemExit(f"record count mismatch: {path}")

    partial = MANIFEST.with_suffix(".tsv.partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    os.replace(MANIFEST, BACKUP)
    os.replace(partial, MANIFEST)
    print(f"PASS\tverified and finalized {len(rows)} manifest rows")


if __name__ == "__main__":
    main()
