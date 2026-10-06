#!/usr/bin/env python3
"""Create the released Stage2-A manifest using current, stable paths."""

from __future__ import annotations

import csv
import hashlib
import io
import os
from pathlib import Path


PROJECT = Path("/root/autodl-tmp/Agent_analysis_v2")
ROOT = PROJECT / "04_stage2/A"
MANIFEST = ROOT / "audit/stage2_A_manifest.tsv"
BUILD_MANIFEST = ROOT / "audit/stage2_A_manifest_buildtime.tsv"
RUN_ID = "stage2_A_current3418_hi_union_rf2_official_precision_20260904"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def records(path: Path) -> int | str:
    if path.suffix != ".tsv":
        return ""
    with path.open("rb") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def render(rows: list[dict[str, object]]) -> bytes:
    stream = io.StringIO(newline="")
    fields = ["role", "path", "records", "bytes", "sha256", "run_id"]
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def main() -> None:
    if not ROOT.is_dir():
        raise SystemExit(f"missing released Stage2-A root: {ROOT}")
    if BUILD_MANIFEST.exists():
        raise SystemExit(f"build manifest backup already exists: {BUILD_MANIFEST}")

    inputs = [
        PROJECT / "03_stage1_A/route_A_candidate_classification_3418.tsv",
        PROJECT / "02_anchors/working_regulatory_anchor_2368.tsv",
        PROJECT / "01_reference/human_swissprot_reviewed_frozen_v1/human_swissprot_reviewed_UP000005640_metadata.tsv",
        PROJECT / "00_rawdataset/04_hi_union_huri_2020/HI-union.tsv",
    ]
    resources = [
        ROOT / "ppi_screening/00_resources/hi_union_ensembl_gene_symbol_map_complete.tsv",
        ROOT / "ppi_screening/00_resources/rf2_final80_from_Data_S4.tsv",
        ROOT / "ppi_screening/00_resources/rf2_final90_from_Data_S3.tsv",
    ]
    for path in inputs + resources:
        if not path.is_file():
            raise SystemExit(f"missing manifest dependency: {path}")

    rows: list[dict[str, object]] = []
    for path in inputs:
        rows.append(
            {
                "role": "INPUT",
                "path": str(path),
                "records": records(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "run_id": RUN_ID,
            }
        )
    for path in resources:
        rows.append(
            {
                "role": "RESOURCE",
                "path": str(path.relative_to(ROOT)),
                "records": records(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "run_id": RUN_ID,
            }
        )
    excluded = {MANIFEST.resolve(), BUILD_MANIFEST.resolve()}
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path.resolve() in excluded or path.name.endswith(".partial"):
            continue
        if path in resources:
            continue
        rows.append(
            {
                "role": "OUTPUT",
                "path": str(path.relative_to(ROOT)),
                "records": records(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "run_id": RUN_ID,
            }
        )

    partial = MANIFEST.with_suffix(".tsv.partial")
    partial.write_bytes(render(rows))
    os.replace(MANIFEST, BUILD_MANIFEST)
    os.replace(partial, MANIFEST)
    print(f"PASS\trelease manifest rows={len(rows)}")


if __name__ == "__main__":
    main()
