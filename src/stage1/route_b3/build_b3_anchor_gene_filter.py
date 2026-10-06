#!/usr/bin/env python3
"""Partition the B3 Stage1 union by current anchor-gene membership."""

from __future__ import annotations

import argparse
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


EXPECTED_B3 = 6_386
EXPECTED_ANCHORS = 2_368
EXPECTED_EXCLUDED = 2_832
EXPECTED_RETAINED = 3_554
RUN_ID = "b3_stage1_anchor_gene_filter_20260902"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_genes(value: str) -> set[str]:
    return {
        token.strip().upper()
        for token in re.split(r"[;|]", str(value))
        if token.strip()
    }


def write_fasta(frame: pd.DataFrame, path: Path, partition: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in frame.itertuples(index=False):
            handle.write(
                f">{row.isoform_accession} parent={row.canonical_accession} "
                f"GN={row.gene_primary} partition={partition}\n"
            )
            sequence = str(row.sequence)
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def write_accessions(frame: pd.DataFrame, path: Path) -> None:
    path.write_text("\n".join(frame["isoform_accession"]) + "\n", encoding="utf-8")


def copy_bytes(source: Path, destination: Path) -> None:
    destination.write_bytes(source.read_bytes())


def manifest_row(
    role: str,
    path: Path,
    records: int,
    description: str,
    generated_at: str,
) -> dict[str, object]:
    return {
        "artifact_role": role,
        "path": str(path.resolve()),
        "records": records,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "description": description,
        "run_id": RUN_ID,
        "generated_at_utc": generated_at,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--b3-root",
        type=Path,
        default=Path("/root/autodl-tmp/Agent_analysis_v2/03_stage1_B/B3"),
    )
    parser.add_argument(
        "--anchor-table",
        type=Path,
        default=Path(
            "/root/autodl-tmp/Agent_analysis_v2/02_anchors/"
            "working_regulatory_anchor_2368.tsv"
        ),
    )
    args = parser.parse_args()

    union_path = args.b3_root / "02_candidates/b3_stage1_nuclear_localization_union.tsv"
    candidate_dir = args.b3_root / "02_candidates"
    audit_dir = args.b3_root / "03_audit"
    candidate_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    union = pd.read_csv(union_path, sep="\t", dtype=str, keep_default_na=False)
    anchors = pd.read_csv(args.anchor_table, sep="\t", dtype=str, keep_default_na=False)
    if len(union) != EXPECTED_B3 or union["isoform_accession"].nunique() != EXPECTED_B3:
        raise AssertionError("B3 Stage1 union must contain 6,386 unique isoforms")
    if len(anchors) != EXPECTED_ANCHORS or anchors["anchor_accession"].nunique() != EXPECTED_ANCHORS:
        raise AssertionError("Current anchor table must contain 2,368 unique accessions")

    anchor_accessions = set(anchors["anchor_accession"].str.strip())
    gene_to_rows: dict[str, list[dict[str, str]]] = {}
    accession_to_row: dict[str, dict[str, str]] = {}
    for row in anchors.to_dict("records"):
        accession_to_row[row["anchor_accession"].strip()] = row
        for gene in atomic_genes(row["anchor_gene"]):
            gene_to_rows.setdefault(gene, []).append(row)

    enriched_rows: list[dict[str, object]] = []
    for row in union.to_dict("records"):
        parent = row["canonical_accession"].strip()
        gene = row["gene_primary"].strip().upper()
        parent_hit = parent in anchor_accessions
        gene_rows = gene_to_rows.get(gene, [])
        gene_hit = bool(gene_rows)

        matched_rows: dict[str, dict[str, str]] = {}
        if parent_hit:
            matched_rows[parent] = accession_to_row[parent]
        for anchor_row in gene_rows:
            matched_rows[anchor_row["anchor_accession"].strip()] = anchor_row

        matched_accessions = sorted(matched_rows)
        matched_genes = sorted(
            {
                token
                for anchor_row in matched_rows.values()
                for token in atomic_genes(anchor_row["anchor_gene"])
            }
        )
        matched_sources = sorted(
            {
                anchor_row["anchor_source_membership"].strip()
                for anchor_row in matched_rows.values()
                if anchor_row["anchor_source_membership"].strip()
            }
        )
        excluded = parent_hit or gene_hit
        enriched = {
            "anchor_parent_accession_match": int(parent_hit),
            "anchor_primary_gene_match": int(gene_hit),
            "anchor_gene_excluded": int(excluded),
            "stage2_partition": (
                "ANCHOR_GENE_ISOFORM_AUDIT" if excluded else "PRIMARY_NON_ANCHOR_GENE"
            ),
            "matched_anchor_accessions": ";".join(matched_accessions),
            "matched_anchor_genes": ";".join(matched_genes),
            "matched_anchor_source_memberships": ";".join(matched_sources),
            **row,
        }
        enriched_rows.append(enriched)

    partitioned = pd.DataFrame(enriched_rows)
    excluded = partitioned[partitioned["anchor_gene_excluded"].eq(1)].copy()
    retained = partitioned[partitioned["anchor_gene_excluded"].eq(0)].copy()
    if len(excluded) != EXPECTED_EXCLUDED or len(retained) != EXPECTED_RETAINED:
        raise AssertionError(
            f"Unexpected partition sizes: excluded={len(excluded)}, retained={len(retained)}"
        )
    if not (excluded["anchor_parent_accession_match"] == excluded["anchor_primary_gene_match"]).all():
        raise AssertionError("Parent-accession and primary-gene anchor matching disagree")

    sort_columns = ["canonical_accession", "isoform_accession"]
    excluded = excluded.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)
    retained = retained.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)
    partitioned = partitioned.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)

    retained_tsv = candidate_dir / "b3_stage1_non_anchor_gene_for_stage2.tsv"
    retained_acc = candidate_dir / "b3_stage1_non_anchor_gene_for_stage2.accessions.txt"
    retained_fasta = candidate_dir / "b3_stage1_non_anchor_gene_for_stage2.fasta"
    excluded_tsv = candidate_dir / "b3_stage1_anchor_gene_isoform_audit.tsv"
    excluded_acc = candidate_dir / "b3_stage1_anchor_gene_isoform_audit.accessions.txt"
    excluded_fasta = candidate_dir / "b3_stage1_anchor_gene_isoform_audit.fasta"
    partition_path = audit_dir / "b3_anchor_gene_filter_per_isoform.tsv"
    summary_path = audit_dir / "b3_anchor_gene_filter_summary.tsv"
    manifest_path = audit_dir / "b3_anchor_gene_filter_manifest.tsv"

    retained.to_csv(retained_tsv, sep="\t", index=False)
    excluded.to_csv(excluded_tsv, sep="\t", index=False)
    partitioned.to_csv(partition_path, sep="\t", index=False)
    write_accessions(retained, retained_acc)
    write_accessions(excluded, excluded_acc)
    write_fasta(retained, retained_fasta, "PRIMARY_NON_ANCHOR_GENE")
    write_fasta(excluded, excluded_fasta, "ANCHOR_GENE_ISOFORM_AUDIT")

    retained_counted = candidate_dir / f"b3_stage1_non_anchor_gene_for_stage2_{len(retained)}.tsv"
    retained_acc_counted = candidate_dir / f"b3_stage1_non_anchor_gene_for_stage2_{len(retained)}.accessions.txt"
    retained_fasta_counted = candidate_dir / f"b3_stage1_non_anchor_gene_for_stage2_{len(retained)}.fasta"
    excluded_counted = candidate_dir / f"b3_stage1_anchor_gene_isoform_audit_{len(excluded)}.tsv"
    copy_bytes(retained_tsv, retained_counted)
    copy_bytes(retained_acc, retained_acc_counted)
    copy_bytes(retained_fasta, retained_fasta_counted)
    copy_bytes(excluded_tsv, excluded_counted)

    summary_rows = [
        ("stage1_union_isoforms", len(union), "immutable B3 Stage1 union"),
        ("stage1_union_parent_genes", union["canonical_accession"].nunique(), "parent canonical accessions"),
        ("current_anchor_accessions", len(anchor_accessions), "canonical anchor proteins"),
        ("current_anchor_atomic_genes", len(gene_to_rows), "uppercase gene tokens"),
        ("anchor_gene_isoforms_excluded_from_primary_stage2", len(excluded), "retained in audit partition"),
        ("anchor_gene_parents_excluded", excluded["canonical_accession"].nunique(), "anchor parent genes"),
        ("non_anchor_gene_isoforms_for_primary_stage2", len(retained), "formal B3 Stage2 input"),
        ("non_anchor_gene_parents_for_primary_stage2", retained["canonical_accession"].nunique(), "non-anchor parent genes"),
        ("excluded_uniprot_direct_only", int((excluded["stage1_evidence_class"] == "UNIPROT_DIRECT").sum()), "audit stratum"),
        ("excluded_model_consensus_only", int((excluded["stage1_evidence_class"] == "MODEL_CONSENSUS").sum()), "audit stratum"),
        ("excluded_both_evidence", int((excluded["stage1_evidence_class"] == "BOTH").sum()), "audit stratum"),
        ("retained_uniprot_direct_only", int((retained["stage1_evidence_class"] == "UNIPROT_DIRECT").sum()), "primary Stage2 stratum"),
        ("retained_model_consensus_only", int((retained["stage1_evidence_class"] == "MODEL_CONSENSUS").sum()), "primary Stage2 stratum"),
        ("retained_both_evidence", int((retained["stage1_evidence_class"] == "BOTH").sum()), "primary Stage2 stratum"),
        ("primary_stage2_pair_universe_2368_anchors", len(retained) * len(anchors), "isoform-anchor pairs"),
        ("full_unfiltered_pair_universe_2368_anchors", len(union) * len(anchors), "reference only"),
    ]
    pd.DataFrame(summary_rows, columns=["metric", "value", "note"]).to_csv(
        summary_path, sep="\t", index=False
    )

    generated_at = datetime.now(timezone.utc).isoformat()
    artifacts = [
        ("INPUT", union_path, len(union), "immutable B3 Stage1 union"),
        ("INPUT", args.anchor_table, len(anchors), "current canonical anchor table"),
        ("OUTPUT", retained_tsv, len(retained), "primary non-anchor-gene Stage2 table"),
        ("OUTPUT", retained_acc, len(retained), "primary Stage2 accession interface"),
        ("OUTPUT", retained_fasta, len(retained), "primary Stage2 sequences"),
        ("OUTPUT", retained_counted, len(retained), "count-labelled primary Stage2 table"),
        ("OUTPUT", retained_acc_counted, len(retained), "count-labelled primary accession interface"),
        ("OUTPUT", retained_fasta_counted, len(retained), "count-labelled primary sequences"),
        ("AUDIT", excluded_tsv, len(excluded), "anchor-gene isoforms excluded from primary Stage2"),
        ("AUDIT", excluded_acc, len(excluded), "excluded isoform accessions"),
        ("AUDIT", excluded_fasta, len(excluded), "excluded isoform sequences"),
        ("AUDIT", excluded_counted, len(excluded), "count-labelled excluded audit table"),
        ("AUDIT", partition_path, len(partitioned), "per-isoform matching and partition audit"),
        ("AUDIT", summary_path, len(summary_rows), "partition counts and pair universe"),
    ]
    pd.DataFrame(
        [manifest_row(*artifact, generated_at) for artifact in artifacts]
    ).to_csv(manifest_path, sep="\t", index=False)

    print("validation_status\tPASS")
    print(f"stage1_union\t{len(union)}")
    print(f"anchor_gene_isoform_audit\t{len(excluded)}")
    print(f"primary_non_anchor_gene_stage2\t{len(retained)}")
    print(f"primary_stage2_pairs\t{len(retained) * len(anchors)}")


if __name__ == "__main__":
    main()
