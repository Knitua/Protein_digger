#!/usr/bin/env python3
"""Read-only validation for the all-alternative-isoform B3 Stage1 release."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_ids(path: Path) -> list[str]:
    result = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            token = line[1:].split()[0]
            result.append(token.split("|", 1)[0])
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--b3-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.b3_root
    all_rows = read_tsv(root / "01_nuclear_evidence/02_model_scores/b3_deeploc_nlsexplorer_scores_22131.tsv")
    direct_rows = read_tsv(root / "01_nuclear_evidence/01_uniprot_direct/b3_uniprot_isoform_specific_nucleus_266.tsv")
    model_rows = read_tsv(root / "01_nuclear_evidence/02_model_scores/b3_model_consensus_both_gt0p5.tsv")
    candidates = read_tsv(root / "02_candidates/b3_stage1_nuclear_localization_union.tsv")
    if len(all_rows) != 22131 or len(direct_rows) != 266:
        raise RuntimeError("unexpected all-input or direct-evidence count")
    all_ids = {r["isoform_accession"] for r in all_rows}
    direct_ids = {r["isoform_accession"] for r in direct_rows}
    model_ids = {r["isoform_accession"] for r in model_rows}
    candidate_ids = {r["isoform_accession"] for r in candidates}
    if len(all_ids) != 22131 or len(candidate_ids) != len(candidates):
        raise RuntimeError("accession uniqueness failure")
    if not all(float(r["deeploc_nucleus_prob"]) > 0.5 and float(r["nls_prob_max"]) > 0.5 for r in model_rows):
        raise RuntimeError("model-positive table violates strict >0.5 rule")
    if candidate_ids != direct_ids | model_ids:
        raise RuntimeError("candidate table is not the required evidence union")
    if set(fasta_ids(root / "02_candidates/b3_stage1_nuclear_localization_union.fasta")) != candidate_ids:
        raise RuntimeError("candidate FASTA accession set mismatch")
    listed = {
        line.strip()
        for line in (root / "02_candidates/b3_stage1_nuclear_localization_union.accessions.txt").read_text().splitlines()
        if line.strip()
    }
    if listed != candidate_ids:
        raise RuntimeError("candidate accession-list mismatch")
    print("validation_status\tPASS")
    print(f"all_alternative_isoforms\t{len(all_ids)}")
    print(f"uniprot_direct\t{len(direct_ids)}")
    print(f"model_consensus\t{len(model_ids)}")
    print(f"stage1_union\t{len(candidate_ids)}")


if __name__ == "__main__":
    main()
