#!/usr/bin/env python3
"""Read-only validation for the B3 anchor-gene partition."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_accessions(path: Path) -> list[str]:
    result: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                result.append(line[1:].split()[0])
    return result


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

    cdir = args.b3_root / "02_candidates"
    adir = args.b3_root / "03_audit"
    union = pd.read_csv(
        cdir / "b3_stage1_nuclear_localization_union.tsv",
        sep="\t", dtype=str, keep_default_na=False,
    )
    retained_path = cdir / "b3_stage1_non_anchor_gene_for_stage2.tsv"
    excluded_path = cdir / "b3_stage1_anchor_gene_isoform_audit.tsv"
    retained = pd.read_csv(retained_path, sep="\t", dtype=str, keep_default_na=False)
    excluded = pd.read_csv(excluded_path, sep="\t", dtype=str, keep_default_na=False)
    anchors = pd.read_csv(args.anchor_table, sep="\t", dtype=str, keep_default_na=False)

    assert len(union) == 6_386
    assert len(retained) == 3_554
    assert len(excluded) == 2_832
    assert len(anchors) == 2_368
    union_ids = set(union["isoform_accession"])
    retained_ids = set(retained["isoform_accession"])
    excluded_ids = set(excluded["isoform_accession"])
    assert not (retained_ids & excluded_ids)
    assert retained_ids | excluded_ids == union_ids
    assert retained["anchor_gene_excluded"].eq("0").all()
    assert excluded["anchor_gene_excluded"].eq("1").all()
    assert excluded["anchor_parent_accession_match"].eq("1").all()
    assert excluded["anchor_primary_gene_match"].eq("1").all()

    anchor_accessions = set(anchors["anchor_accession"].str.strip())
    anchor_genes: set[str] = set()
    for value in anchors["anchor_gene"]:
        anchor_genes.update(
            token.strip().upper()
            for token in re.split(r"[;|]", value)
            if token.strip()
        )
    assert not (set(retained["canonical_accession"]) & anchor_accessions)
    assert not (set(retained["gene_primary"].str.upper()) & anchor_genes)
    assert set(excluded["canonical_accession"]).issubset(anchor_accessions)
    assert set(excluded["gene_primary"].str.upper()).issubset(anchor_genes)

    for prefix, frame in [
        ("b3_stage1_non_anchor_gene_for_stage2", retained),
        ("b3_stage1_anchor_gene_isoform_audit", excluded),
    ]:
        accessions = [
            line.strip()
            for line in (cdir / f"{prefix}.accessions.txt").read_text().splitlines()
            if line.strip()
        ]
        assert accessions == frame["isoform_accession"].tolist()
        assert fasta_accessions(cdir / f"{prefix}.fasta") == accessions

    assert (cdir / "b3_stage1_non_anchor_gene_for_stage2_3554.tsv").read_bytes() == retained_path.read_bytes()
    assert (cdir / "b3_stage1_anchor_gene_isoform_audit_2832.tsv").read_bytes() == excluded_path.read_bytes()

    summary = pd.read_csv(adir / "b3_anchor_gene_filter_summary.tsv", sep="\t")
    values = dict(zip(summary["metric"], summary["value"]))
    assert int(values["stage1_union_isoforms"]) == 6_386
    assert int(values["anchor_gene_isoforms_excluded_from_primary_stage2"]) == 2_832
    assert int(values["non_anchor_gene_isoforms_for_primary_stage2"]) == 3_554
    assert int(values["primary_stage2_pair_universe_2368_anchors"]) == 8_415_872

    manifest = pd.read_csv(adir / "b3_anchor_gene_filter_manifest.tsv", sep="\t")
    for row in manifest.itertuples(index=False):
        path = Path(row.path)
        assert path.is_file(), path
        assert path.stat().st_size == int(row.bytes), path
        assert sha256(path) == row.sha256, path

    print("validation_status\tPASS")
    print("stage1_union\t6386")
    print("anchor_gene_isoform_audit\t2832")
    print("primary_non_anchor_gene_stage2\t3554")
    print("primary_stage2_pairs\t8415872")


if __name__ == "__main__":
    main()
