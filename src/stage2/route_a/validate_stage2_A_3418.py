#!/usr/bin/env python3
"""Validate current Stage2-A cohorts and six PPI dataset screens."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path


EXPECTED = {"A1_WIDE": 2877, "A2_NOVELTY": 2641, "A3_FOCUSED": 2457}
DATASET_FILES = {
    "A1_WIDE": "stage2_A1_wide_remove_known_curated_2877.tsv",
    "A2_NOVELTY": "stage2_A2_novelty_remove_curated_and_direct_2641.tsv",
    "A3_FOCUSED": "stage2_A3_focused_associated_under_mixed_2457.tsv",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def accession_set(rows: list[dict[str, str]]) -> set[str]:
    values = [row["candidate_accession"].strip() for row in rows]
    if not all(values) or len(values) != len(set(values)):
        raise AssertionError("blank or duplicate candidate accession")
    return set(values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expect(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def locate_one(directory: Path, pattern: str) -> Path:
    matches = list(directory.glob(pattern))
    expect(len(matches) == 1, f"expected one {pattern} in {directory}; found {len(matches)}")
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--previous-root", type=Path)
    parser.add_argument("--audit-out", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()

    cohorts: dict[str, set[str]] = {}
    result: dict[str, object] = {"status": "PASS", "root": str(root), "datasets": {}}
    for dataset, expected in EXPECTED.items():
        path = root / "datasets" / DATASET_FILES[dataset]
        rows = read_tsv(path)
        expect(len(rows) == expected, f"{dataset} row count != {expected}")
        cohorts[dataset] = accession_set(rows)
        expect(all(row["stage2_A_dataset"] == dataset for row in rows), f"{dataset} label mismatch")
    expect(cohorts["A3_FOCUSED"] < cohorts["A2_NOVELTY"] < cohorts["A1_WIDE"], "cohorts not strictly nested")

    rf2_sets: dict[str, set[tuple[str, str]]] = {}
    for label, expected in [("final80", 29257), ("final90", 17849)]:
        rf2_rows = read_tsv(root / f"ppi_screening/00_resources/rf2_{label}_from_Data_S{'4' if label == 'final80' else '3'}.tsv")
        expect(len(rf2_rows) == expected, f"RF2 {label} row count != {expected}")
        pairs = {
            tuple(sorted((row["Protein1"].strip(), row["Protein2"].strip())))
            for row in rf2_rows
        }
        expect(len(pairs) == expected, f"RF2 {label} pairs not unique")
        rf2_sets[label] = pairs
    expect(rf2_sets["final90"] < rf2_sets["final80"], "RF2 Data S3 is not strict subset of Data S4")

    pass_sets: dict[str, dict[str, set[str]]] = {name: {} for name in EXPECTED}
    for dataset, expected in EXPECTED.items():
        dataset_result: dict[str, object] = {}
        for label in ["final80", "final90"]:
            directory = root / f"ppi_screening/02_hi_union_rf2_{label}"
            all_path = directory / f"{dataset}_{expected}_hi_union_rf2_{label}_screen_all.tsv"
            passed_path = locate_one(directory, f"{dataset}_hi_union_rf2_{label}_passed_*.tsv")
            pair_path = directory / f"{dataset}_hi_union_rf2_{label}_candidate_anchor_pairs.tsv"
            summary_path = directory / f"{dataset}_hi_union_rf2_{label}_summary.tsv"
            all_rows = read_tsv(all_path)
            passed_rows = read_tsv(passed_path)
            pair_rows = read_tsv(pair_path)
            expect(len(all_rows) == expected, f"{dataset}/{label} full screen count mismatch")
            expect(accession_set(all_rows) == cohorts[dataset], f"{dataset}/{label} screen universe mismatch")
            passed = accession_set(passed_rows)
            pass_sets[dataset][label] = passed
            expect(passed <= cohorts[dataset], f"{dataset}/{label} passed outside input")
            expect(all(row["ppi_pass"] == "yes" for row in passed_rows), f"{dataset}/{label} passed table has no")
            expect(
                passed == {row["candidate_accession"] for row in all_rows if row["ppi_pass"] == "yes"},
                f"{dataset}/{label} passed table differs from screen_all",
            )
            expect(
                {row["candidate_accession"] for row in pair_rows} == passed,
                f"{dataset}/{label} candidate-anchor pairs do not cover exactly passed candidates",
            )
            pair_keys = [(row["candidate_accession"], row["anchor_gene"]) for row in pair_rows]
            expect(len(pair_keys) == len(set(pair_keys)), f"{dataset}/{label} duplicate candidate-anchor gene pair")
            summary = {row["metric"]: row["value"] for row in read_tsv(summary_path)}
            expect(int(summary["input_accessions"]) == expected, f"{dataset}/{label} summary input mismatch")
            expect(int(summary["passed_accessions"]) == len(passed), f"{dataset}/{label} summary pass mismatch")
            expect(int(summary["candidate_anchor_pairs"]) == len(pair_rows), f"{dataset}/{label} summary pair mismatch")
            dataset_result[label] = {
                "input_accessions": expected,
                "passed_accessions": len(passed),
                "candidate_anchor_pairs": len(pair_rows),
                "passed_file": passed_path.name,
            }
        expect(pass_sets[dataset]["final90"] <= pass_sets[dataset]["final80"], f"{dataset}: final90 not subset final80")
        result["datasets"][dataset] = dataset_result

    for label in ["final80", "final90"]:
        expect(
            pass_sets["A3_FOCUSED"][label]
            <= pass_sets["A2_NOVELTY"][label]
            <= pass_sets["A1_WIDE"][label],
            f"nested pass-set relation failed for {label}",
        )

    # Independent comparison to the previous A1 output for all retained accessions.
    if args.previous_root:
        previous = args.previous_root.resolve()
        compare_fields = [
            "ppi_pass", "ppi_evidence_category", "merged_anchor_count", "merged_anchor_genes",
            "hi_union_anchor_count", "hi_union_anchor_genes", "rf2_anchor_count",
            "rf2_anchor_genes", "rf2_exact_accession_anchor_count",
            "rf2_exact_accession_anchor_genes", "rf2_gene_propagated_anchor_count",
            "rf2_gene_propagated_anchor_genes", "rf2_match_mode",
        ]
        previous_a1 = {}
        for label in ["final80", "final90"]:
            old_path = previous / f"ppi_screening/02_hi_union_rf2_{label}/A1_WIDE_2942_hi_union_rf2_{label}_screen_all.tsv"
            old_rows = read_tsv(old_path)
            expect(len(old_rows) == 2942, f"previous A1/{label} count != 2942")
            old_map = {row["candidate_accession"]: row for row in old_rows}
            new_rows = read_tsv(
                root / f"ppi_screening/02_hi_union_rf2_{label}/A1_WIDE_2877_hi_union_rf2_{label}_screen_all.tsv"
            )
            mismatches = []
            for row in new_rows:
                old = old_map.get(row["candidate_accession"])
                if old is None or any(old[field] != row[field] for field in compare_fields):
                    mismatches.append(row["candidate_accession"])
            expect(not mismatches, f"previous A1 cross-check mismatches for {label}: {mismatches[:5]}")
            previous_a1[label] = {
                "previous_rows": len(old_rows),
                "current_rows": len(new_rows),
                "retained_evidence_rows_exact_match": len(new_rows),
                "previous_file_sha256": sha256(old_path),
            }
        result["previous_A1_crosscheck"] = previous_a1

    manifest_path = root / "audit/stage2_A_manifest.tsv"
    manifest_rows = read_tsv(manifest_path)
    verified_manifest_rows = 0
    for row in manifest_rows:
        if row["role"] in {"OUTPUT", "RESOURCE"}:
            path = root / row["path"]
        else:
            path = Path(row["path"])
        expect(path.is_file(), f"manifest path missing: {path}")
        expect(path.stat().st_size == int(row["bytes"]), f"manifest size mismatch: {path}")
        expect(sha256(path) == row["sha256"], f"manifest sha256 mismatch: {path}")
        verified_manifest_rows += 1
    result["manifest_rows_verified"] = verified_manifest_rows

    content = (json.dumps(result, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    audit_out = args.audit_out.resolve()
    audit_out.parent.mkdir(parents=True, exist_ok=True)
    if audit_out.exists():
        if audit_out.read_bytes() != content:
            raise FileExistsError(f"refusing to overwrite different audit: {audit_out}")
    else:
        partial = audit_out.with_name(audit_out.name + ".partial")
        with partial.open("xb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(partial, audit_out)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
