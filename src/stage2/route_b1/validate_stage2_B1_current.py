#!/usr/bin/env python3
"""Independent validation for current Stage2-B1 final80/final90 outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        return list(reader.fieldnames), rows


def tokens(value: str | None) -> set[str]:
    return {
        token.strip().upper()
        for token in re.split(r"[;|]", (value or "").strip())
        if token.strip()
    }


def one_match(directory: Path, pattern: str) -> Path:
    paths = list(directory.glob(pattern))
    if len(paths) != 1:
        raise AssertionError(f"Expected one {pattern} in {directory}; found {paths}")
    return paths[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage2-b1", type=Path, required=True)
    parser.add_argument("--stage1-b1", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    args = parser.parse_args()

    _, stage1 = read_tsv(args.stage1_b1)
    if len(stage1) != 1350:
        raise AssertionError(f"Expected 1,350 Stage1-B1 rows; found {len(stage1)}")
    input_accessions = [row["accession"].strip() for row in stage1]
    if len(set(input_accessions)) != 1350 or "" in input_accessions:
        raise AssertionError("Stage1-B1 accession interface is not unique and complete")
    input_sequence = {row["accession"].strip(): row["sequence"].strip() for row in stage1}

    _, anchors = read_tsv(args.anchors)
    if len(anchors) != 2368:
        raise AssertionError("Working anchor table does not contain 2,368 rows")
    anchor_accessions = {row["anchor_accession"].strip() for row in anchors}
    if len(anchor_accessions) != 2368 or "" in anchor_accessions:
        raise AssertionError("Anchor accessions are not unique")
    anchor_mapping = {
        row["anchor_accession"].strip(): row["anchor_gene"].strip().upper()
        for row in anchors
    }

    results: dict[str, dict[str, object]] = {}
    for label, subdir in (
        ("final80", "01_hi_union_rf2_final80"),
        ("final90", "02_hi_union_rf2_final90"),
    ):
        directory = args.stage2_b1 / subdir
        all_path = directory / f"stage2_B1_1350_hi_union_rf2_{label}_screen_all.tsv"
        passed_path = one_match(directory, f"stage2_B1_hi_union_rf2_{label}_passed_*.tsv")
        pair_path = directory / f"stage2_B1_hi_union_rf2_{label}_candidate_anchor_pairs.tsv"
        summary_path = directory / f"stage2_B1_hi_union_rf2_{label}_summary.tsv"
        _, all_rows = read_tsv(all_path)
        _, passed_rows = read_tsv(passed_path)
        _, pairs = read_tsv(pair_path)
        _, summary_rows = read_tsv(summary_path)

        if len(all_rows) != 1350:
            raise AssertionError(f"{label}: screen_all is not 1,350 rows")
        all_accessions = [row["candidate_accession"].strip() for row in all_rows]
        if set(all_accessions) != set(input_accessions) or len(set(all_accessions)) != 1350:
            raise AssertionError(f"{label}: screen_all accession universe differs from Stage1-B1")
        for row in all_rows:
            accession = row["candidate_accession"].strip()
            if row["accession"].strip() != accession:
                raise AssertionError(f"{label}: duplicated accession columns disagree")
            if row["sequence"].strip() != input_sequence[accession]:
                raise AssertionError(f"{label}: candidate sequence differs from Stage1 canonical sequence")
            passed = row["ppi_pass"] == "yes"
            if passed != (int(row["merged_anchor_count"]) > 0):
                raise AssertionError(f"{label}: pass flag and merged anchor count disagree")

        passed_accessions = {row["candidate_accession"].strip() for row in passed_rows}
        calculated_passed = {
            row["candidate_accession"].strip()
            for row in all_rows
            if row["ppi_pass"] == "yes"
        }
        if len(passed_accessions) != len(passed_rows) or passed_accessions != calculated_passed:
            raise AssertionError(f"{label}: passed table differs from screen_all pass flags")

        pair_keys: set[tuple[str, str]] = set()
        pair_candidates: set[str] = set()
        pair_anchor_genes: dict[str, set[str]] = {}
        for row in pairs:
            candidate = row["candidate_accession"].strip()
            anchor = row["anchor_accession"].strip()
            key = (candidate, anchor)
            if key in pair_keys:
                raise AssertionError(f"{label}: duplicate candidate-anchor pair {key}")
            pair_keys.add(key)
            pair_candidates.add(candidate)
            pair_anchor_genes.setdefault(candidate, set()).add(row["anchor_gene"].strip().upper())
            if candidate not in passed_accessions:
                raise AssertionError(f"{label}: pair belongs to non-passed candidate")
            if anchor not in anchor_accessions:
                raise AssertionError(f"{label}: pair anchor is outside current 2,368 table")
            if row["anchor_gene"].strip().upper() != anchor_mapping[anchor]:
                raise AssertionError(f"{label}: pair anchor accession/gene mapping disagrees")
            hi = row["hi_union_support"] == "1"
            rf2 = row["rf2_support"] == "1"
            if not (hi or rf2):
                raise AssertionError(f"{label}: pair has no evidence")
            expected_category = (
                f"HI-union+RF2-{label}" if hi and rf2
                else ("HI-union-only" if hi else f"RF2-{label}-only")
            )
            if row["pair_evidence_category"] != expected_category:
                raise AssertionError(f"{label}: pair evidence category is inconsistent")
        if pair_candidates != passed_accessions:
            raise AssertionError(f"{label}: every passed candidate must have at least one pair")
        for row in passed_rows:
            candidate = row["candidate_accession"].strip()
            if tokens(row["merged_anchor_genes"]) != pair_anchor_genes[candidate]:
                raise AssertionError(f"{label}: candidate-level and pair-level anchor genes disagree")

        summary = {row["metric"]: int(row["value"]) for row in summary_rows}
        expected_summary = {
            "input_accessions": 1350,
            "passed_accessions": len(passed_rows),
            "not_passed_accessions": 1350 - len(passed_rows),
            "candidate_anchor_pairs": len(pairs),
            "distinct_supported_anchors": len({row["anchor_accession"] for row in pairs}),
        }
        for metric, value in expected_summary.items():
            if summary.get(metric) != value:
                raise AssertionError(f"{label}: summary {metric} disagrees")
        results[label] = {
            "passed_accessions": passed_accessions,
            "pairs": pair_keys,
            "passed_count": len(passed_rows),
            "pair_count": len(pairs),
        }

    if not results["final90"]["passed_accessions"] <= results["final80"]["passed_accessions"]:
        raise AssertionError("final90 passed candidates are not a subset of final80")
    if not results["final90"]["pairs"] <= results["final80"]["pairs"]:
        raise AssertionError("final90 candidate-anchor pairs are not a subset of final80")

    comparison_path = args.stage2_b1 / "03_comparison/stage2_B1_final80_vs_final90.tsv"
    _, comparison = read_tsv(comparison_path)
    if len(comparison) != 1350:
        raise AssertionError("80/90 comparison table is not 1,350 rows")
    for row in comparison:
        accession = row["candidate_accession"].strip()
        in80 = accession in results["final80"]["passed_accessions"]
        in90 = accession in results["final90"]["passed_accessions"]
        expected = "PASS_BOTH" if in90 else ("PASS_FINAL80_ONLY" if in80 else "PASS_NEITHER")
        if row["threshold_comparison_status"] != expected:
            raise AssertionError("80/90 comparison status is inconsistent")

    audit_path = args.stage2_b1 / "04_audit/stage2_B1_ppi_screen_audit.json"
    audit = json.loads(audit_path.read_text())
    if audit["input_accessions"] != 1350 or audit["current_anchor_accessions"] != 2368:
        raise AssertionError("Audit input counts are inconsistent")
    if audit["final80"]["passed_accessions"] != results["final80"]["passed_count"]:
        raise AssertionError("Audit final80 count is inconsistent")
    if audit["final90"]["passed_accessions"] != results["final90"]["passed_count"]:
        raise AssertionError("Audit final90 count is inconsistent")

    manifest_path = args.stage2_b1 / "04_audit/stage2_B1_ppi_screen_manifest.tsv"
    _, manifest = read_tsv(manifest_path)
    for row in manifest:
        path = Path(row["path"])
        if not path.is_absolute():
            path = args.stage2_b1 / path
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise AssertionError(f"Manifest mismatch: {path}")
    partials = list(args.stage2_b1.rglob("*.partial"))
    if partials:
        raise AssertionError(f"Unresolved partial files: {partials}")

    print("PASS\tStage2-B1 current final80/final90 validation")
    print(f"input_accessions\t{len(stage1)}")
    print(f"final80_passed\t{results['final80']['passed_count']}")
    print(f"final80_pairs\t{results['final80']['pair_count']}")
    print(f"final90_passed\t{results['final90']['passed_count']}")
    print(f"final90_pairs\t{results['final90']['pair_count']}")


if __name__ == "__main__":
    main()
