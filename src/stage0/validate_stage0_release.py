#!/usr/bin/env python3
"""Validate the released Stage0 A/B partition and its Stage1 interfaces."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def accession_set(path: Path, column: str = "Entry") -> set[str]:
    rows = read_tsv(path)
    values = [row[column].strip() for row in rows]
    if not all(values):
        raise AssertionError(f"blank accession in {path}")
    if len(values) != len(set(values)):
        raise AssertionError(f"duplicate accession in {path}")
    return set(values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expect(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/root/autodl-tmp/Agent_analysis_v2"),
    )
    args = parser.parse_args()
    root = args.project_root.resolve()
    s0 = root / "03_stage0_AB_partition"
    a1 = root / "03_stage1_A"
    b = root / "03_stage1_B"

    paths = {
        "a_full": s0 / "01_canonical_partition/route_A_canonical_nucleus_5669.tsv",
        "b_full": s0 / "01_canonical_partition/route_B_without_canonical_nucleus_14747.tsv",
        "a_excluded": s0 / "02_anchor_exclusion/route_A_anchor_excluded_2251.tsv",
        "b_excluded": s0 / "02_anchor_exclusion/route_B_anchor_excluded_122.tsv",
        "a_input": s0 / "02_anchor_exclusion/route_A_non_anchor_3418.tsv",
        "b_input": s0 / "02_anchor_exclusion/route_B_non_anchor_14625.tsv",
        "iso73": s0 / "03_audit/canonical_nucleus_restricted_to_non_displayed_isoforms_73.tsv",
        "iso69": s0 / "03_audit/route_B_isoform_only_nucleus_non_anchor_69.tsv",
        "a_class": a1 / "route_A_candidate_classification_3418.tsv",
        "b_interface": b / "00_common_input/02_non_anchor/route_B_uniprot_without_canonical_nucleus_non_anchor_14625.tsv",
    }
    for label, path in paths.items():
        expect(path.is_file(), f"missing {label}: {path}")

    a_full = accession_set(paths["a_full"])
    b_full = accession_set(paths["b_full"])
    a_excluded = accession_set(paths["a_excluded"])
    b_excluded = accession_set(paths["b_excluded"])
    a_input = accession_set(paths["a_input"])
    b_input = accession_set(paths["b_input"])

    expect(len(a_full) == 5669, "A full count != 5669")
    expect(len(b_full) == 14747, "B full count != 14747")
    expect(not (a_full & b_full), "A/B full sets overlap")
    expect(len(a_full | b_full) == 20416, "A/B union != 20416")

    expect(len(a_excluded) == 2251, "A excluded count != 2251")
    expect(len(b_excluded) == 122, "B excluded count != 122")
    expect(len(a_input) == 3418, "A non-anchor count != 3418")
    expect(len(b_input) == 14625, "B non-anchor count != 14625")
    expect(a_excluded.isdisjoint(a_input), "A excluded/input overlap")
    expect(b_excluded.isdisjoint(b_input), "B excluded/input overlap")
    expect(a_excluded | a_input == a_full, "A exclusion partition incomplete")
    expect(b_excluded | b_input == b_full, "B exclusion partition incomplete")
    expect(a_input.isdisjoint(b_input), "A/B Stage1 inputs overlap")
    expect(len(a_input | b_input) == 18043, "non-anchor union != 18043")

    iso73 = read_tsv(paths["iso73"])
    expect(len(iso73) == 73, "isoform-only nucleus parent count != 73")
    expect({row["canonical_accession"] for row in iso73} <= b_full, "isoform-only parents not all in B")
    iso69 = accession_set(paths["iso69"])
    expect(len(iso69) == 69, "isoform-only non-anchor parent count != 69")
    expect(iso69 <= b_input, "isoform-only non-anchor parents not all in B input")

    a_rows = read_tsv(paths["a_class"])
    a_class_accessions = {row["candidate_accession"] for row in a_rows}
    expect(len(a_rows) == 3418, "A classification row count != 3418")
    expect(a_class_accessions == a_input, "A classification set differs from Stage0 A input")
    level1 = Counter(row["regulatory_evidence_status"] for row in a_rows)
    expect(
        level1
        == Counter(
            {
                "KNOWN_CURATED": 541,
                "KNOWN_FUNCTIONAL_DIRECT": 236,
                "KNOWN_COMPLEX_ASSOCIATED": 440,
                "NO_KNOWN_REGULATORY_EVIDENCE": 2201,
            }
        ),
        f"unexpected A Level1 counts: {dict(level1)}",
    )
    no_known = [
        row for row in a_rows if row["regulatory_evidence_status"] == "NO_KNOWN_REGULATORY_EVIDENCE"
    ]
    level2 = Counter(row["final_tier"] for row in no_known)
    expect(
        level2
        == Counter(
            {
                "HK_PROBABLE": 170,
                "UNDER_CHARACTERIZED": 924,
                "MIXED_MOONLIGHTING": 1093,
                "EVIDENCE_INCOMPLETE": 14,
            }
        ),
        f"unexpected A Level2 counts: {dict(level2)}",
    )

    b_interface = accession_set(paths["b_interface"])
    expect(b_interface == b_input, "B common interface differs from Stage0 B input")

    for readme in [s0 / "README.md", a1 / "README.md", b / "README.md"]:
        text = readme.read_text(encoding="utf-8")
        for forbidden in ["旧版", "补丁", "修正前", "修正后", "误放"]:
            expect(forbidden not in text, f"historical patch wording in {readme}: {forbidden}")

    tracked = [path for path in paths.values() if path.is_file()]
    print("PASS\tStage0 release validation")
    print("A_FULL\t5669")
    print("B_FULL\t14747")
    print("A_NON_ANCHOR\t3418")
    print("B_NON_ANCHOR\t14625")
    print("NON_ANCHOR_UNION\t18043")
    print("ISOFORM_ONLY_NUCLEUS_PARENTS\t73")
    print("ISOFORM_ONLY_NON_ANCHOR_PARENTS\t69")
    for path in tracked:
        print(f"SHA256\t{path.relative_to(root)}\t{sha256(path)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, KeyError, OSError) as exc:
        print(f"FAIL\t{exc}", file=sys.stderr)
        raise SystemExit(1)
