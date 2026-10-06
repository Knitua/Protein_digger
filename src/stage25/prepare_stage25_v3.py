#!/usr/bin/env python3
"""Prepare target-agnostic Stage2.5 core and positive-control pair universes."""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

from stage25_v3_common import (
    KMER_CONTAINMENT_PREFILTER_MIN,
    LENGTH_RATIO_PREFILTER_MIN,
    RUN_ID,
    atomic_write_tsv,
    build_universe,
    kmer_set,
)


def make_pairs(left: list[dict], right_by_signature: dict[str, list[dict]], directed: bool) -> list[dict]:
    kmers = {r["accession"]: kmer_set(r["sequence"]) for r in left}
    for records in right_by_signature.values():
        for rec in records:
            kmers.setdefault(rec["accession"], kmer_set(rec["sequence"]))
    rows = []
    seen = set()
    for a in left:
        signature = a["coarse_domain_signature"]
        if not signature:
            continue
        for b in right_by_signature.get(signature, []):
            if a["accession"] == b["accession"]:
                continue
            if not directed:
                key = tuple(sorted((a["accession"], b["accession"])))
                if key in seen:
                    continue
                seen.add(key)
            ratio = min(a["length"], b["length"]) / max(a["length"], b["length"])
            if ratio < LENGTH_RATIO_PREFILTER_MIN:
                continue
            ak, bk = kmers[a["accession"]], kmers[b["accession"]]
            denom = min(len(ak), len(bk))
            containment = len(ak & bk) / denom if denom else 0.0
            if containment < KMER_CONTAINMENT_PREFILTER_MIN:
                continue
            pair_id = a["accession"] + "__" + b["accession"]
            rows.append(
                {
                    "run_id": RUN_ID,
                    "pair_id": pair_id,
                    "left_accession": a["accession"],
                    "left_gene": a["gene"],
                    "left_stage1_route": a["stage1_route"],
                    "right_accession": b["accession"],
                    "right_gene": b["gene"],
                    "right_stage1_route": b["stage1_route"],
                    "coarse_domain_basis": a["coarse_domain_basis"],
                    "coarse_domain_signature": signature,
                    "length_ratio": f"{ratio:.8f}",
                    "kmer_containment": f"{containment:.8f}",
                    "left_sequence": a["sequence"],
                    "right_sequence": b["sequence"],
                }
            )
    rows.sort(key=lambda r: (r["left_accession"], r["right_accession"]))
    return rows


def run_aligner(aligner: Path, pair_rows: list[dict], input_path: Path, output_path: Path) -> None:
    input_path.parent.mkdir(parents=True, exist_ok=True)
    partial = input_path.with_name(input_path.name + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["pair_id", "seed_accession", "candidate_accession", "seed_sequence", "candidate_sequence"])
        for row in pair_rows:
            writer.writerow(
                [
                    row["pair_id"],
                    row["left_accession"],
                    row["right_accession"],
                    row["left_sequence"],
                    row["right_sequence"],
                ]
            )
    partial.replace(input_path)
    out_partial = output_path.with_name(output_path.name + ".partial")
    subprocess.run([str(aligner), str(input_path), str(out_partial)], check=True)
    out_partial.replace(output_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--aligner", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    project = args.project_root.resolve()
    out = args.output_dir.resolve()
    universe, _ = build_universe(project, args.annotations.resolve())
    seeds = sorted((x for x in universe if x["primary_final80_pass"]), key=lambda x: x["accession"])
    review = sorted((x for x in universe if not x["primary_final80_pass"]), key=lambda x: x["accession"])

    seed_by_signature: dict[str, list[dict]] = defaultdict(list)
    review_by_signature: dict[str, list[dict]] = defaultdict(list)
    for rec in seeds:
        seed_by_signature[rec["coarse_domain_signature"]].append(rec)
    for rec in review:
        review_by_signature[rec["coarse_domain_signature"]].append(rec)

    core_pairs = make_pairs(seeds, review_by_signature, directed=True)
    positive_pairs = make_pairs(seeds, seed_by_signature, directed=False)

    pair_fields = [
        "run_id", "pair_id", "left_accession", "left_gene", "left_stage1_route",
        "right_accession", "right_gene", "right_stage1_route", "coarse_domain_basis",
        "coarse_domain_signature", "length_ratio", "kmer_containment",
    ]
    atomic_write_tsv(out / "work" / "core_prefilter_pairs.tsv", core_pairs, pair_fields)
    atomic_write_tsv(out / "work" / "positive_control_prefilter_pairs.tsv", positive_pairs, pair_fields)
    run_aligner(args.aligner, core_pairs, out / "work" / "core_alignment_input.tsv", out / "work" / "core_alignment_output.tsv")
    run_aligner(
        args.aligner,
        positive_pairs,
        out / "work" / "positive_control_alignment_input.tsv",
        out / "work" / "positive_control_alignment_output.tsv",
    )

    needed = sorted({r[k] for r in core_pairs + positive_pairs for k in ("left_accession", "right_accession")})
    atomic_write_tsv(
        out / "resources" / "domain_accessions_needed.tsv",
        [{"accession": acc} for acc in needed],
        ["accession"],
    )

    seed_rows = [
        {
            "run_id": RUN_ID,
            "seed_accession": r["accession"],
            "seed_gene": r["gene"],
            "stage1_route": r["stage1_route"],
            "primary_final80_category": r["primary_final80_category"],
            "also_primary_final90_pass": "yes" if r["primary_final90_pass"] else "no",
            "length": r["length"],
            "coarse_domain_basis": r["coarse_domain_basis"],
            "coarse_domain_signature": r["coarse_domain_signature"],
        }
        for r in seeds
    ]
    review_rows = [
        {
            "run_id": RUN_ID,
            "candidate_accession": r["accession"],
            "candidate_gene": r["gene"],
            "stage1_route": r["stage1_route"],
            "primary_final80_status": "fail",
            "primary_final80_category": r["primary_final80_category"],
            "length": r["length"],
            "coarse_domain_basis": r["coarse_domain_basis"],
            "coarse_domain_signature": r["coarse_domain_signature"],
        }
        for r in review
    ]
    atomic_write_tsv(out / "inputs" / "stage2_5_v3_final80_positive_seeds_1862.tsv", seed_rows)
    atomic_write_tsv(out / "inputs" / "stage2_5_v3_final80_failed_review_pool_2631.tsv", review_rows)
    print(
        {
            "core_prefilter_pairs": len(core_pairs),
            "positive_control_prefilter_pairs": len(positive_pairs),
            "domain_accessions_needed": len(needed),
        }
    )


if __name__ == "__main__":
    main()
