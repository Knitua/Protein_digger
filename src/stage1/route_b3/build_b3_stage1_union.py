#!/usr/bin/env python3
"""Merge B3 model scores with direct UniProt isoform-specific nucleus evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


TAB = "\t"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=TAB))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=TAB, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                f">{row['isoform_accession']}|canonical_parent={row['canonical_accession']}"
                f"|gene={row['gene_primary']}|evidence={row['stage1_evidence_class']}\n"
                f"{row['sequence']}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--direct-uniprot", type=Path, required=True)
    parser.add_argument("--deeploc", type=Path, required=True)
    parser.add_argument("--nls", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    args = parser.parse_args()

    metadata = read_tsv(args.metadata)
    direct_rows = read_tsv(args.direct_uniprot)
    deeploc_rows = read_tsv(args.deeploc)
    nls_rows = read_tsv(args.nls)

    meta_by_id = {row["isoform_accession"]: row for row in metadata}
    direct_by_id = {row["isoform_accession"]: row for row in direct_rows}
    deeploc_by_id = {row["accession"]: row for row in deeploc_rows}
    nls_by_id = {row["accession"]: row for row in nls_rows}
    expected = set(meta_by_id)
    if len(metadata) != len(expected) or len(expected) != 22131:
        raise RuntimeError("metadata must contain exactly 22,131 unique alternative isoforms")
    if set(deeploc_by_id) != expected:
        raise RuntimeError("DeepLoc accession coverage does not exactly match the 22,131 input")
    if set(nls_by_id) != expected:
        raise RuntimeError("NLSExplorer accession coverage does not exactly match the 22,131 input")
    if len(direct_by_id) != 266 or not set(direct_by_id) <= expected:
        raise RuntimeError("direct UniProt evidence must contain 266 input isoforms")

    score_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    model_rows: list[dict[str, object]] = []
    for metadata_row in metadata:
        isoform = metadata_row["isoform_accession"]
        deeploc_prob = float(deeploc_by_id[isoform]["deeploc_nucleus_prob"])
        nls_prob = float(nls_by_id[isoform]["nls_prob_max"])
        direct = isoform in direct_by_id
        deeploc_positive = deeploc_prob > 0.5
        nls_positive = nls_prob > 0.5
        model_positive = deeploc_positive and nls_positive
        stage1_positive = direct or model_positive
        if direct and model_positive:
            evidence_class = "BOTH"
        elif direct:
            evidence_class = "UNIPROT_DIRECT"
        elif model_positive:
            evidence_class = "MODEL_CONSENSUS"
        else:
            evidence_class = "NONE"
        direct_row = direct_by_id.get(isoform, {})
        row: dict[str, object] = {
            **metadata_row,
            "uniprot_isoform_specific_nucleus": int(direct),
            "uniprot_isoform_label": direct_row.get("isoform_label", ""),
            "uniprot_location_sections": direct_row.get("uniprot_location_sections", ""),
            "corrected_parent_group_73": direct_row.get("corrected_parent_group_73", 0),
            "deeploc_nucleus_prob": f"{deeploc_prob:.8f}",
            "nls_prob_max": f"{nls_prob:.8f}",
            "deeploc_positive_gt0p5": int(deeploc_positive),
            "nls_positive_gt0p5": int(nls_positive),
            "model_consensus_both_gt0p5": int(model_positive),
            "stage1_positive": int(stage1_positive),
            "stage1_evidence_class": evidence_class,
        }
        score_rows.append(row)
        if model_positive:
            model_rows.append(row)
        if stage1_positive:
            candidate_rows.append(row)

    appended_fields = [
        "uniprot_isoform_specific_nucleus",
        "uniprot_isoform_label",
        "uniprot_location_sections",
        "corrected_parent_group_73",
        "deeploc_nucleus_prob",
        "nls_prob_max",
        "deeploc_positive_gt0p5",
        "nls_positive_gt0p5",
        "model_consensus_both_gt0p5",
        "stage1_positive",
        "stage1_evidence_class",
    ]
    fields = list(metadata[0]) + appended_fields
    scores_path = args.out_root / "01_nuclear_evidence/02_model_scores/b3_deeploc_nlsexplorer_scores_22131.tsv"
    model_path = args.out_root / "01_nuclear_evidence/02_model_scores/b3_model_consensus_both_gt0p5.tsv"
    candidates_path = args.out_root / "02_candidates/b3_stage1_nuclear_localization_union.tsv"
    write_tsv(scores_path, fields, score_rows)
    write_tsv(model_path, fields, model_rows)
    write_tsv(candidates_path, fields, candidate_rows)
    (candidates_path.parent / "b3_stage1_nuclear_localization_union.accessions.txt").write_text(
        "".join(f"{row['isoform_accession']}\n" for row in candidate_rows), encoding="utf-8"
    )
    write_fasta(candidates_path.parent / "b3_stage1_nuclear_localization_union.fasta", candidate_rows)

    direct_ids = set(direct_by_id)
    model_ids = {str(row["isoform_accession"]) for row in model_rows}
    union_ids = {str(row["isoform_accession"]) for row in candidate_rows}
    if union_ids != direct_ids | model_ids:
        raise RuntimeError("candidate union is not direct evidence union model consensus")
    summary = [
        {"metric": "alternative_isoform_input", "value": len(metadata)},
        {"metric": "uniprot_direct", "value": len(direct_ids)},
        {"metric": "deeploc_gt0p5", "value": sum(int(r["deeploc_positive_gt0p5"]) for r in score_rows)},
        {"metric": "nls_gt0p5", "value": sum(int(r["nls_positive_gt0p5"]) for r in score_rows)},
        {"metric": "model_consensus_both_gt0p5", "value": len(model_ids)},
        {"metric": "direct_and_model", "value": len(direct_ids & model_ids)},
        {"metric": "direct_only", "value": len(direct_ids - model_ids)},
        {"metric": "model_only", "value": len(model_ids - direct_ids)},
        {"metric": "stage1_union", "value": len(union_ids)},
        {"metric": "stage1_union_unique_parents", "value": len({str(r["canonical_accession"]) for r in candidate_rows})},
        {"metric": "stage1_union_unique_genes_nonblank", "value": len({str(r["gene_primary"]) for r in candidate_rows if r["gene_primary"]})},
    ]
    summary_path = args.out_root / "03_audit/b3_stage1_counts.tsv"
    write_tsv(summary_path, ["metric", "value"], summary)

    manifest_paths = [
        args.metadata,
        args.direct_uniprot,
        args.deeploc,
        args.nls,
        scores_path,
        model_path,
        candidates_path,
        candidates_path.parent / "b3_stage1_nuclear_localization_union.accessions.txt",
        candidates_path.parent / "b3_stage1_nuclear_localization_union.fasta",
        summary_path,
    ]
    write_tsv(
        args.out_root / "03_audit/b3_stage1_manifest.tsv",
        ["path", "bytes", "sha256"],
        [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in manifest_paths
        ],
    )
    for row in summary:
        print(f"{row['metric']}\t{row['value']}")


if __name__ == "__main__":
    main()
