#!/usr/bin/env python3
"""Read-only validation for current Stage2-B2 HI-union/RF2 screens."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


RUN_ID = "stage2_B2_current266_hi_union_rf2_official_precision_20260904"


def read(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metric_map(path: Path) -> dict[str, str]:
    frame = read(path)
    if list(frame.columns) != ["metric", "value"]:
        raise ValueError(f"Unexpected summary schema: {path}")
    return dict(zip(frame["metric"], frame["value"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage1-input", type=Path, required=True)
    parser.add_argument("--stage2-root", type=Path, required=True)
    parser.add_argument("--anchor-table", type=Path, required=True)
    args = parser.parse_args()

    source = read(args.stage1_input)
    anchors = read(args.anchor_table)
    source_set = set(source["accession"])
    anchor_accessions = set(anchors["anchor_accession"])
    anchor_genes = set(anchors["anchor_gene"].str.upper())
    if len(source) != 266 or len(source_set) != 266:
        raise ValueError("Stage1-B2 input is not 266 unique accessions")
    if any("-" in value for value in source_set):
        raise ValueError("Stage1-B2 contains an isoform accession")
    if source_set & anchor_accessions:
        raise ValueError("Stage1-B2 overlaps anchor accessions")

    pass_sets: dict[str, set[str]] = {}
    for label, dirname in [("final80", "01_hi_union_rf2_final80"), ("final90", "02_hi_union_rf2_final90")]:
        directory = args.stage2_root / dirname
        all_path = directory / f"stage2_B2_266_hi_union_rf2_{label}_screen_all.tsv"
        passed_paths = list(directory.glob(f"stage2_B2_hi_union_rf2_{label}_passed_*.tsv"))
        if len(passed_paths) != 1:
            raise ValueError(f"Expected exactly one passed table for {label}")
        all_rows = read(all_path)
        passed = read(passed_paths[0])
        summary = metric_map(directory / f"stage2_B2_hi_union_rf2_{label}_summary.tsv")
        if len(all_rows) != 266 or set(all_rows["candidate_accession"]) != source_set:
            raise ValueError(f"{label} screen_all does not exactly cover Stage1-B2")
        if set(all_rows["stage1_route"]) != {"B2"}:
            raise ValueError(f"{label} route label differs")
        if set(all_rows["ppi_run_id"]) != {RUN_ID}:
            raise ValueError(f"{label} run ID differs")
        merged = pd.to_numeric(all_rows["merged_anchor_count"], errors="raise")
        expected_pass = set(all_rows.loc[merged > 0, "candidate_accession"])
        observed_pass = set(passed["candidate_accession"])
        if observed_pass != expected_pass:
            raise ValueError(f"{label} passed table differs from merged-anchor evidence")
        if set(all_rows.loc[merged > 0, "ppi_pass"]) != {"yes"}:
            raise ValueError(f"{label} positive pass labels differ")
        if set(all_rows.loc[merged == 0, "ppi_pass"]) - {"no"}:
            raise ValueError(f"{label} negative pass labels differ")
        for value in all_rows.loc[merged > 0, "merged_anchor_genes"]:
            if not set(filter(None, value.upper().split(";"))) <= anchor_genes:
                raise ValueError(f"{label} output contains a non-anchor partner")
        if int(summary["input_accessions"]) != 266:
            raise ValueError(f"{label} summary input count differs")
        if int(summary["passed_accessions"]) != len(observed_pass):
            raise ValueError(f"{label} summary pass count differs")
        pass_sets[label] = observed_pass

    if not pass_sets["final90"] <= pass_sets["final80"]:
        raise ValueError("final90 pass set is not a subset of final80")

    comparison = read(args.stage2_root / "03_comparison/stage2_B2_final80_vs_final90.tsv")
    if len(comparison) != 266 or set(comparison["candidate_accession"]) != source_set:
        raise ValueError("Threshold comparison does not cover all candidates")
    audit_path = args.stage2_root / "04_audit/stage2_B2_ppi_screen_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit["input_accessions"] != 266 or audit["run_id"] != RUN_ID:
        raise ValueError("Stage2 audit metadata differs")
    manifest = read(args.stage2_root / "04_audit/stage2_B2_ppi_screen_manifest.tsv")
    for _, row in manifest.iterrows():
        artifact = Path(row["path"])
        if not artifact.is_file():
            raise FileNotFoundError(artifact)
        if artifact.stat().st_size != int(row["size_bytes"]) or sha256(artifact) != row["sha256"]:
            raise ValueError(f"Manifest mismatch: {artifact}")
    print(
        "PASS: Stage2-B2 validated; "
        f"final80={len(pass_sets['final80'])}, final90={len(pass_sets['final90'])}"
    )


if __name__ == "__main__":
    main()
