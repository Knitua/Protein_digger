#!/usr/bin/env python3
"""Build current Stage2-A cohorts and six HI-union/RF2 screening outputs.

The program only writes to a new output directory. It never overwrites a
different existing file. RF2 final80/final90 are the official precision sets
from Data S4/Data S3, not thresholds applied to a single probability column.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


ROOT = Path("/root/autodl-tmp/Agent_analysis_v2")
STAGE1 = ROOT / "03_stage1_A"
ANCHORS = ROOT / "02_anchors/working_regulatory_anchor_2368.tsv"
MASTER = STAGE1 / "route_A_candidate_classification_3418.tsv"
METADATA = ROOT / "01_reference/human_swissprot_reviewed_frozen_v1/human_swissprot_reviewed_UP000005640_metadata.tsv"
HI_UNION = ROOT / "00_rawdataset/04_hi_union_huri_2020/HI-union.tsv"
ENSEMBL_MAP_SOURCE = Path("/root/autodl-tmp/Agent_analysis_stage2_recovery_20260728/source_data/hi_union_ensembl_gene_symbol_map_complete.tsv")
RF2_FINAL80_SOURCE = Path("/root/autodl-tmp/Agent_analysis_stage2_recovery_20260728/source_data/rf2_final_predictions/rf2_final_predictions_80_from_Data_S4.tsv")
RF2_FINAL90_SOURCE = ROOT / "04_stage2/A/ppi_screening/00_resources/rf2_final90_from_Data_S3.tsv"

RUN_ID = "stage2_A_current3418_hi_union_rf2_official_precision_20260904"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames is None:
            raise ValueError(f"missing header: {path}")
        return list(reader.fieldnames), rows


def atomic_genes(value: str | None) -> set[str]:
    return {
        token.strip().upper()
        for token in re.split(r"[;|]", (value or "").strip())
        if token.strip()
    }


def joined(values: Iterable[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def render_tsv(fields: list[str], rows: Iterable[dict[str, object]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def safe_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() == content:
            return
        raise FileExistsError(f"refusing to overwrite different file: {path}")
    partial = path.with_name(path.name + ".partial")
    if partial.exists():
        raise FileExistsError(f"partial already exists: {partial}")
    with partial.open("xb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def safe_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if sha256(source) == sha256(destination):
            return
        raise FileExistsError(f"refusing to overwrite different file: {destination}")
    partial = destination.with_name(destination.name + ".partial")
    if partial.exists():
        raise FileExistsError(f"partial already exists: {partial}")
    with source.open("rb") as src, partial.open("xb") as dst:
        shutil.copyfileobj(src, dst, 1024 * 1024)
        dst.flush()
        os.fsync(dst.fileno())
    if sha256(source) != sha256(partial):
        raise ValueError(f"copy checksum mismatch: {source}")
    os.replace(partial, destination)


def count_atomic_genes(rows: Iterable[dict[str, object]]) -> int:
    values: set[str] = set()
    for row in rows:
        values |= atomic_genes(str(row.get("candidate_gene", "")))
    return len(values)


def tsv_record_count(path: Path) -> int | None:
    if path.suffix != ".tsv":
        return None
    with path.open("rb") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--rf2-final90-source", type=Path, default=RF2_FINAL90_SOURCE)
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists():
        raise FileExistsError(f"output directory already exists: {out}")

    inputs = [
        MASTER,
        ANCHORS,
        METADATA,
        HI_UNION,
        ENSEMBL_MAP_SOURCE,
        RF2_FINAL80_SOURCE,
        args.rf2_final90_source,
    ]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)

    master_fields, master_rows = read_tsv(MASTER)
    if len(master_rows) != 3418:
        raise ValueError(f"expected 3,418 Stage1-A rows, found {len(master_rows)}")
    master_accessions = [row["candidate_accession"].strip() for row in master_rows]
    if len(set(master_accessions)) != 3418 or any(not value for value in master_accessions):
        raise ValueError("Stage1-A accessions are blank or duplicated")

    level1 = Counter(row["regulatory_evidence_status"] for row in master_rows)
    if level1 != Counter(
        {
            "KNOWN_CURATED": 541,
            "KNOWN_FUNCTIONAL_DIRECT": 236,
            "KNOWN_COMPLEX_ASSOCIATED": 440,
            "NO_KNOWN_REGULATORY_EVIDENCE": 2201,
        }
    ):
        raise ValueError(f"unexpected Stage1-A Level1 counts: {dict(level1)}")

    no_known_rows = [
        row for row in master_rows
        if row["regulatory_evidence_status"] == "NO_KNOWN_REGULATORY_EVIDENCE"
    ]
    level2 = Counter(row["final_tier"] for row in no_known_rows)
    if level2 != Counter(
        {
            "HK_PROBABLE": 170,
            "UNDER_CHARACTERIZED": 924,
            "MIXED_MOONLIGHTING": 1093,
            "EVIDENCE_INCOMPLETE": 14,
        }
    ):
        raise ValueError(f"unexpected Stage1-A Level2 counts: {dict(level2)}")

    cohort_rules = {
        "A1_WIDE": lambda row: row["regulatory_evidence_status"] != "KNOWN_CURATED",
        "A2_NOVELTY": lambda row: row["regulatory_evidence_status"] not in {
            "KNOWN_CURATED", "KNOWN_FUNCTIONAL_DIRECT"
        },
        "A3_FOCUSED": lambda row: (
            row["regulatory_evidence_status"] == "KNOWN_COMPLEX_ASSOCIATED"
            or (
                row["regulatory_evidence_status"] == "NO_KNOWN_REGULATORY_EVIDENCE"
                and row["final_tier"] in {"UNDER_CHARACTERIZED", "MIXED_MOONLIGHTING"}
            )
        ),
    }
    selection_text = {
        "A1_WIDE": "REMOVE KNOWN_CURATED",
        "A2_NOVELTY": "REMOVE KNOWN_CURATED and KNOWN_FUNCTIONAL_DIRECT",
        "A3_FOCUSED": "RETAIN KNOWN_COMPLEX_ASSOCIATED plus UNDER_CHARACTERIZED and MIXED_MOONLIGHTING",
    }
    expected_sizes = {"A1_WIDE": 2877, "A2_NOVELTY": 2641, "A3_FOCUSED": 2457}
    file_stems = {
        "A1_WIDE": "stage2_A1_wide_remove_known_curated_2877",
        "A2_NOVELTY": "stage2_A2_novelty_remove_curated_and_direct_2641",
        "A3_FOCUSED": "stage2_A3_focused_associated_under_mixed_2457",
    }
    cohort_rows: dict[str, list[dict[str, object]]] = {}
    cohort_sets: dict[str, set[str]] = {}
    dataset_fields = [
        "stage2_A_dataset",
        "stage2_A_selection_rule",
        "stage1_A_level1_class",
        "stage1_A_level2_class",
    ] + master_fields
    for name, predicate in cohort_rules.items():
        rows: list[dict[str, object]] = []
        for source in master_rows:
            if not predicate(source):
                continue
            row: dict[str, object] = {
                "stage2_A_dataset": name,
                "stage2_A_selection_rule": selection_text[name],
                "stage1_A_level1_class": source["regulatory_evidence_status"],
                "stage1_A_level2_class": (
                    source["final_tier"]
                    if source["regulatory_evidence_status"] == "NO_KNOWN_REGULATORY_EVIDENCE"
                    else "NOT_APPLICABLE"
                ),
            }
            row.update(source)
            rows.append(row)
        if len(rows) != expected_sizes[name]:
            raise ValueError(f"unexpected {name} size: {len(rows)}")
        cohort_rows[name] = rows
        cohort_sets[name] = {str(row["candidate_accession"]) for row in rows}
        safe_write(
            out / "datasets" / f"{file_stems[name]}.tsv",
            render_tsv(dataset_fields, rows),
        )
    if not cohort_sets["A3_FOCUSED"] < cohort_sets["A2_NOVELTY"] < cohort_sets["A1_WIDE"]:
        raise ValueError("Stage2-A datasets are not strictly nested")

    dataset_count_rows = [
        {
            "dataset": name,
            "accessions": len(cohort_rows[name]),
            "atomic_genes": count_atomic_genes(cohort_rows[name]),
            "selection_rule": selection_text[name],
            "output_file": f"datasets/{file_stems[name]}.tsv",
        }
        for name in ["A1_WIDE", "A2_NOVELTY", "A3_FOCUSED"]
    ]
    safe_write(
        out / "audit/stage2_A_dataset_counts.tsv",
        render_tsv(
            ["dataset", "accessions", "atomic_genes", "selection_rule", "output_file"],
            dataset_count_rows,
        ),
    )

    # Freeze the RF2 tables and the Ensembl mapping inside the release.
    resources = out / "ppi_screening/00_resources"
    ensembl_map = resources / "hi_union_ensembl_gene_symbol_map_complete.tsv"
    rf2_final80 = resources / "rf2_final80_from_Data_S4.tsv"
    rf2_final90 = resources / "rf2_final90_from_Data_S3.tsv"
    safe_copy(ENSEMBL_MAP_SOURCE, ensembl_map)
    safe_copy(RF2_FINAL80_SOURCE, rf2_final80)
    safe_copy(args.rf2_final90_source, rf2_final90)

    _, anchor_rows = read_tsv(ANCHORS)
    _, metadata_rows = read_tsv(METADATA)
    _, map_rows = read_tsv(ensembl_map)
    if len(anchor_rows) != 2368 or any(row.get("use_as_ppi_anchor") != "1" for row in anchor_rows):
        raise ValueError("working anchor table is not the frozen 2,368-row PPI anchor set")

    anchor_accessions = {row["anchor_accession"].strip() for row in anchor_rows}
    if len(anchor_accessions) != 2368:
        raise ValueError("anchor accessions are blank or duplicated")
    anchor_accession_to_genes: dict[str, set[str]] = {}
    anchor_gene_to_accessions: dict[str, set[str]] = defaultdict(set)
    anchor_genes: set[str] = set()
    for row in anchor_rows:
        accession = row["anchor_accession"].strip()
        genes = atomic_genes(row.get("anchor_gene"))
        anchor_accession_to_genes[accession] = genes
        for gene in genes:
            anchor_genes.add(gene)
            anchor_gene_to_accessions[gene].add(accession)

    accession_to_genes: dict[str, set[str]] = {}
    gene_to_accessions: dict[str, set[str]] = defaultdict(set)
    for row in metadata_rows:
        accession = row.get("Entry", "").strip()
        genes = atomic_genes(row.get("Gene Names (primary)"))
        if accession:
            accession_to_genes[accession] = genes
        for gene in genes:
            gene_to_accessions[gene].add(accession)

    for name, rows in cohort_rows.items():
        for row in rows:
            if atomic_genes(str(row.get("candidate_gene", ""))) & anchor_genes:
                raise ValueError(f"{name} contains an anchor atomic gene")
            if str(row["candidate_accession"]) in anchor_accessions:
                raise ValueError(f"{name} contains an anchor accession")

    ensembl_to_gene = {
        row["ensembl_gene_id"].strip().split(".", 1)[0]: row["display_name"].strip().upper()
        for row in map_rows
        if row.get("ensembl_gene_id", "").strip() and row.get("display_name", "").strip()
    }
    hi_partner_to_anchor_genes: dict[str, set[str]] = defaultdict(set)
    hi_rows: list[dict[str, object]] = []
    hi_seen: set[tuple[str, str]] = set()
    hi_raw_rows = 0
    hi_unmapped: set[str] = set()
    with HI_UNION.open("r", encoding="utf-8") as handle:
        for row_number, values in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(values) < 2:
                raise ValueError(f"malformed HI-union row {row_number}")
            hi_raw_rows += 1
            id1, id2 = values[0].strip().split(".", 1)[0], values[1].strip().split(".", 1)[0]
            gene1, gene2 = ensembl_to_gene.get(id1, ""), ensembl_to_gene.get(id2, "")
            if not gene1:
                hi_unmapped.add(id1)
            if not gene2:
                hi_unmapped.add(id2)
            if not gene1 or not gene2:
                continue
            orientations: list[tuple[str, str, str, str]] = []
            if gene1 in anchor_genes and gene2 not in anchor_genes:
                orientations.append((gene1, gene2, id1, id2))
            if gene2 in anchor_genes and gene1 not in anchor_genes:
                orientations.append((gene2, gene1, id2, id1))
            for anchor_gene, partner_gene, anchor_id, partner_id in orientations:
                key = (anchor_gene, partner_gene)
                if key in hi_seen:
                    continue
                hi_seen.add(key)
                hi_partner_to_anchor_genes[partner_gene].add(anchor_gene)
                hi_rows.append(
                    {
                        "anchor_gene": anchor_gene,
                        "anchor_accessions": joined(anchor_gene_to_accessions[anchor_gene]),
                        "anchor_ensembl_gene_id": anchor_id,
                        "partner_gene": partner_gene,
                        "partner_canonical_accessions": joined(gene_to_accessions.get(partner_gene, set())),
                        "partner_ensembl_gene_id": partner_id,
                        "evidence_source": "HI-union",
                    }
                )
    hi_fields = [
        "anchor_gene", "anchor_accessions", "anchor_ensembl_gene_id", "partner_gene",
        "partner_canonical_accessions", "partner_ensembl_gene_id", "evidence_source",
    ]
    safe_write(
        out / "ppi_screening/01_common_hi_union/hi_union_current2368_anchor_partner_edges.tsv",
        render_tsv(hi_fields, sorted(hi_rows, key=lambda row: (str(row["partner_gene"]), str(row["anchor_gene"])))),
    )

    evidence_fields = [
        "ppi_screen_version", "ppi_pass", "ppi_evidence_category", "merged_anchor_count",
        "merged_anchor_genes", "hi_union_anchor_count", "hi_union_anchor_genes",
        "rf2_anchor_count", "rf2_anchor_genes", "rf2_exact_accession_anchor_count",
        "rf2_exact_accession_anchor_genes", "rf2_gene_propagated_anchor_count",
        "rf2_gene_propagated_anchor_genes", "rf2_match_mode", "run_id",
    ]
    pair_fields = [
        "candidate_dataset", "candidate_accession", "candidate_gene", "anchor_gene",
        "anchor_accessions", "hi_union_support", "rf2_support",
        "rf2_exact_accession_support", "rf2_gene_propagated_support",
        "ppi_evidence_category", "ppi_screen_version", "run_id",
    ]

    audit: dict[str, object] = {
        "run_id": RUN_ID,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "anchor_accessions": len(anchor_accessions),
        "anchor_atomic_genes": len(anchor_genes),
        "hi_union_raw_rows": hi_raw_rows,
        "hi_union_unmapped_ensembl_ids": len(hi_unmapped),
        "hi_union_current_anchor_partner_gene_edges": len(hi_rows),
        "datasets": {},
    }
    pass_sets: dict[str, dict[str, set[str]]] = defaultdict(dict)

    for label, rf2_path, expected in [
        ("final80", rf2_final80, 29257),
        ("final90", rf2_final90, 17849),
    ]:
        _, rf2_rows = read_tsv(rf2_path)
        if len(rf2_rows) != expected:
            raise ValueError(f"unexpected RF2 {label} size: {len(rf2_rows)}")
        pair_keys = {
            tuple(sorted((row["Protein1"].strip(), row["Protein2"].strip())))
            for row in rf2_rows
        }
        if len(pair_keys) != expected:
            raise ValueError(f"duplicate RF2 pairs in {label}")

        rf2_exact_genes: dict[str, set[str]] = defaultdict(set)
        rf2_gene_genes: dict[str, set[str]] = defaultdict(set)
        rf2_exact_accessions: dict[str, set[str]] = defaultdict(set)
        rf2_gene_accessions: dict[str, set[str]] = defaultdict(set)
        rf2_edge_rows: list[dict[str, object]] = []
        rf2_seen: set[tuple[str, str]] = set()
        for row in rf2_rows:
            p1, p2 = row["Protein1"].strip(), row["Protein2"].strip()
            orientations: list[tuple[str, str]] = []
            if p1 in anchor_accessions and p2 not in anchor_accessions:
                orientations.append((p1, p2))
            if p2 in anchor_accessions and p1 not in anchor_accessions:
                orientations.append((p2, p1))
            for anchor_accession, partner_accession in orientations:
                key = (anchor_accession, partner_accession)
                if key in rf2_seen:
                    continue
                rf2_seen.add(key)
                current_anchor_genes = anchor_accession_to_genes[anchor_accession]
                partner_genes = accession_to_genes.get(partner_accession, set())
                if partner_genes & anchor_genes:
                    continue
                rf2_exact_genes[partner_accession] |= current_anchor_genes
                rf2_exact_accessions[partner_accession].add(anchor_accession)
                for gene in partner_genes:
                    rf2_gene_genes[gene] |= current_anchor_genes
                    rf2_gene_accessions[gene].add(anchor_accession)
                rf2_edge_rows.append(
                    {
                        "anchor_accession": anchor_accession,
                        "anchor_genes": joined(current_anchor_genes),
                        "partner_accession": partner_accession,
                        "partner_genes": joined(partner_genes),
                        "RFprob": row.get("RFprob", ""),
                        "AFprob": row.get("AFprob", ""),
                        "AFprob5": row.get("AFprob5", ""),
                        "AFMprob": row.get("AFMprob", ""),
                        "prediction_source": row.get("Source", ""),
                        "confDBs": row.get("confDBs", ""),
                        "allDBs": row.get("allDBs", ""),
                        "evidence_source": f"RF2-{label}",
                    }
                )
        rf2_edge_fields = [
            "anchor_accession", "anchor_genes", "partner_accession", "partner_genes",
            "RFprob", "AFprob", "AFprob5", "AFMprob", "prediction_source",
            "confDBs", "allDBs", "evidence_source",
        ]
        threshold_dir = out / f"ppi_screening/02_hi_union_rf2_{label}"
        safe_write(
            threshold_dir / f"rf2_{label}_current2368_anchor_partner_edges.tsv",
            render_tsv(
                rf2_edge_fields,
                sorted(rf2_edge_rows, key=lambda row: (str(row["partner_accession"]), str(row["anchor_accession"]))),
            ),
        )

        for dataset in ["A1_WIDE", "A2_NOVELTY", "A3_FOCUSED"]:
            all_rows: list[dict[str, object]] = []
            passed_rows: list[dict[str, object]] = []
            candidate_anchor_rows: list[dict[str, object]] = []
            for source in cohort_rows[dataset]:
                accession = str(source["candidate_accession"])
                genes = atomic_genes(str(source.get("candidate_gene", "")))
                hi_anchor_genes: set[str] = set()
                rf2_gene_anchor_genes: set[str] = set()
                rf2_gene_anchor_accessions: set[str] = set()
                for gene in genes:
                    hi_anchor_genes |= hi_partner_to_anchor_genes.get(gene, set())
                    rf2_gene_anchor_genes |= rf2_gene_genes.get(gene, set())
                    rf2_gene_anchor_accessions |= rf2_gene_accessions.get(gene, set())
                rf2_exact_anchor_genes = rf2_exact_genes.get(accession, set())
                rf2_exact_anchor_accessions = rf2_exact_accessions.get(accession, set())
                rf2_anchor_genes = rf2_exact_anchor_genes | rf2_gene_anchor_genes
                merged = hi_anchor_genes | rf2_anchor_genes
                if hi_anchor_genes and rf2_anchor_genes:
                    category = f"HI-union+RF2-{label}"
                elif hi_anchor_genes:
                    category = "HI-union-only"
                elif rf2_anchor_genes:
                    category = f"RF2-{label}-only"
                else:
                    category = "NO_PPI_DATASET_SUPPORT"
                if rf2_exact_anchor_genes and rf2_gene_anchor_genes:
                    match_mode = "EXACT_ACCESSION+GENE"
                elif rf2_exact_anchor_genes:
                    match_mode = "EXACT_ACCESSION"
                elif rf2_gene_anchor_genes:
                    match_mode = "GENE_PROPAGATED"
                else:
                    match_mode = "NONE"
                evidence: dict[str, object] = {
                    "ppi_screen_version": f"HI-union+RF2-{label}",
                    "ppi_pass": "yes" if merged else "no",
                    "ppi_evidence_category": category,
                    "merged_anchor_count": len(merged),
                    "merged_anchor_genes": joined(merged),
                    "hi_union_anchor_count": len(hi_anchor_genes),
                    "hi_union_anchor_genes": joined(hi_anchor_genes),
                    "rf2_anchor_count": len(rf2_anchor_genes),
                    "rf2_anchor_genes": joined(rf2_anchor_genes),
                    "rf2_exact_accession_anchor_count": len(rf2_exact_anchor_genes),
                    "rf2_exact_accession_anchor_genes": joined(rf2_exact_anchor_genes),
                    "rf2_gene_propagated_anchor_count": len(rf2_gene_anchor_genes),
                    "rf2_gene_propagated_anchor_genes": joined(rf2_gene_anchor_genes),
                    "rf2_match_mode": match_mode,
                    "run_id": RUN_ID,
                }
                output = dict(evidence)
                output.update(source)
                all_rows.append(output)
                if merged:
                    passed_rows.append(output)

                for anchor_gene in sorted(merged):
                    hi_support = anchor_gene in hi_anchor_genes
                    exact_accessions = {
                        value for value in rf2_exact_anchor_accessions
                        if anchor_gene in anchor_accession_to_genes[value]
                    }
                    gene_accessions = {
                        value for value in rf2_gene_anchor_accessions
                        if anchor_gene in anchor_accession_to_genes[value]
                    }
                    rf2_support = bool(exact_accessions or gene_accessions)
                    if hi_support and rf2_support:
                        pair_category = f"HI-union+RF2-{label}"
                    elif hi_support:
                        pair_category = "HI-union-only"
                    else:
                        pair_category = f"RF2-{label}-only"
                    candidate_anchor_rows.append(
                        {
                            "candidate_dataset": dataset,
                            "candidate_accession": accession,
                            "candidate_gene": str(source.get("candidate_gene", "")),
                            "anchor_gene": anchor_gene,
                            "anchor_accessions": joined(
                                (anchor_gene_to_accessions.get(anchor_gene, set()) if hi_support else set())
                                | exact_accessions | gene_accessions
                            ),
                            "hi_union_support": "yes" if hi_support else "no",
                            "rf2_support": "yes" if rf2_support else "no",
                            "rf2_exact_accession_support": "yes" if exact_accessions else "no",
                            "rf2_gene_propagated_support": "yes" if gene_accessions else "no",
                            "ppi_evidence_category": pair_category,
                            "ppi_screen_version": f"HI-union+RF2-{label}",
                            "run_id": RUN_ID,
                        }
                    )

            pass_sets[dataset][label] = {str(row["candidate_accession"]) for row in passed_rows}
            counts = Counter(str(row["ppi_evidence_category"]) for row in all_rows)
            stem = f"{dataset}_{expected_sizes[dataset]}_hi_union_rf2_{label}"
            safe_write(threshold_dir / f"{stem}_screen_all.tsv", render_tsv(evidence_fields + dataset_fields, all_rows))
            safe_write(
                threshold_dir / f"{dataset}_hi_union_rf2_{label}_passed_{len(passed_rows)}.tsv",
                render_tsv(evidence_fields + dataset_fields, passed_rows),
            )
            safe_write(
                threshold_dir / f"{dataset}_hi_union_rf2_{label}_candidate_anchor_pairs.tsv",
                render_tsv(pair_fields, candidate_anchor_rows),
            )
            summary_rows = [
                {"metric": "input_accessions", "value": len(all_rows)},
                {"metric": "input_atomic_genes", "value": count_atomic_genes(all_rows)},
                {"metric": "passed_accessions", "value": len(passed_rows)},
                {"metric": "passed_atomic_genes", "value": count_atomic_genes(passed_rows)},
                {"metric": "candidate_anchor_pairs", "value": len(candidate_anchor_rows)},
                {"metric": "not_passed_accessions", "value": len(all_rows) - len(passed_rows)},
                {"metric": "official_rf2_all_pairs", "value": len(rf2_rows)},
                {"metric": "current_anchor_rf2_edges", "value": len(rf2_edge_rows)},
            ] + [
                {"metric": f"evidence_category::{name}", "value": value}
                for name, value in sorted(counts.items())
            ]
            safe_write(
                threshold_dir / f"{dataset}_hi_union_rf2_{label}_summary.tsv",
                render_tsv(["metric", "value"], summary_rows),
            )
            audit["datasets"].setdefault(dataset, {})[label] = {
                "input_accessions": len(all_rows),
                "input_atomic_genes": count_atomic_genes(all_rows),
                "passed_accessions": len(passed_rows),
                "passed_atomic_genes": count_atomic_genes(passed_rows),
                "candidate_anchor_pairs": len(candidate_anchor_rows),
                "evidence_counts": dict(sorted(counts.items())),
            }

    comparison_fields = [
        "candidate_accession", "candidate_gene", "hi_rf2_final80_pass",
        "hi_rf2_final90_pass", "threshold_comparison_status",
    ]
    for dataset in ["A1_WIDE", "A2_NOVELTY", "A3_FOCUSED"]:
        pass80 = pass_sets[dataset]["final80"]
        pass90 = pass_sets[dataset]["final90"]
        if not pass90 <= pass80:
            raise ValueError(f"{dataset}: final90 pass set is not a subset of final80")
        comparison_rows: list[dict[str, object]] = []
        for row in cohort_rows[dataset]:
            accession = str(row["candidate_accession"])
            in80, in90 = accession in pass80, accession in pass90
            comparison_rows.append(
                {
                    "candidate_accession": accession,
                    "candidate_gene": str(row.get("candidate_gene", "")),
                    "hi_rf2_final80_pass": "yes" if in80 else "no",
                    "hi_rf2_final90_pass": "yes" if in90 else "no",
                    "threshold_comparison_status": (
                        "PASS_BOTH" if in90 else ("PASS_FINAL80_ONLY" if in80 else "PASS_NEITHER")
                    ),
                }
            )
        safe_write(
            out / f"ppi_screening/03_comparison/{dataset}_final80_vs_final90.tsv",
            render_tsv(comparison_fields, comparison_rows),
        )
        audit["datasets"][dataset]["final90_pass_is_subset_of_final80"] = True

    # Validate nested pass sets for the three nested candidate inputs.
    for label in ["final80", "final90"]:
        if not (
            pass_sets["A3_FOCUSED"][label]
            <= pass_sets["A2_NOVELTY"][label]
            <= pass_sets["A1_WIDE"][label]
        ):
            raise ValueError(f"nested pass-set validation failed for {label}")

    safe_write(
        out / "ppi_screening/04_audit/stage2_A_ppi_union_audit.json",
        (json.dumps(audit, ensure_ascii=False, indent=2) + "\n").encode("utf-8"),
    )

    counts_text = "\n".join(
        f"| `{dataset}` | {expected_sizes[dataset]:,} | "
        f"{audit['datasets'][dataset]['final80']['passed_accessions']:,} / "
        f"{audit['datasets'][dataset]['final80']['passed_atomic_genes']:,} | "
        f"{audit['datasets'][dataset]['final90']['passed_accessions']:,} / "
        f"{audit['datasets'][dataset]['final90']['passed_atomic_genes']:,} |"
        for dataset in ["A1_WIDE", "A2_NOVELTY", "A3_FOCUSED"]
    )
    readme = f"""# Stage2 Route A：三级嵌套输入与PPI数据集筛选

更新日期：2026-09-04
状态：`COMPLETE / VALIDATED`

本目录直接使用Stage1-A当前3,418条canonical proteins及其冻结分级，构建三个嵌套的Stage2输入，并分别与当前2,368个canonical anchors的HI-union、RF2-PPI官方final80和final90互作背景进行比对。

## 三个输入集合

```text
A1_WIDE     = 3,418 − 541 KNOWN_CURATED = 2,877
A2_NOVELTY  = 2,877 − 236 KNOWN_FUNCTIONAL_DIRECT = 2,641
A3_FOCUSED  = 440 KNOWN_COMPLEX_ASSOCIATED
              + 924 UNDER_CHARACTERIZED
              + 1,093 MIXED_MOONLIGHTING
            = 2,457
```

```text
A3_FOCUSED ⊂ A2_NOVELTY ⊂ A1_WIDE ⊂ Stage1-A
```

## PPI筛选结果

| 数据集 | 输入proteins | HI-union ∪ RF2-final80 proteins / atomic genes | HI-union ∪ RF2-final90 proteins / atomic genes |
|---|---:|---:|---:|
{counts_text}

`final80`和`final90`分别指RF2-PPI论文官方Data S4（80% precision）和Data S3（90% precision）交付集合，不是对RFprob或其他单列设置0.8/0.9阈值。Data S3是Data S4的严格pair子集，因此每个候选集的final90通过集合必须属于final80通过集合。

HI-union按Ensembl gene ID映射到当前gene；RF2优先记录精确UniProt accession命中，同时保留可追溯的同gene传播并通过`rf2_match_mode`区分。所有candidate和anchor均使用canonical蛋白口径。

## 主要文件

```text
datasets/                         三个Stage2-A输入表
ppi_screening/02_hi_union_rf2_final80/
ppi_screening/02_hi_union_rf2_final90/
ppi_screening/03_comparison/      final80与final90逐候选比较
ppi_screening/04_audit/           计数、集合关系与SHA256
```

每个数据集和阈值均提供：完整`screen_all`表、通过表、candidate–anchor证据表和summary表。PPI命中表示已有实验互作数据或RF2-PPI预测支持，不等同于新的直接结合实验验证。
"""
    safe_write(out / "README.md", readme.encode("utf-8"))

    ppi_readme = f"""# Stage2-A PPI数据集筛选

本目录对A1_WIDE、A2_NOVELTY和A3_FOCUSED分别运行：

```text
HI-union ∪ RF2-PPI official final80 (Data S4)
HI-union ∪ RF2-PPI official final90 (Data S3)
```

三个候选集合、两种RF2精度背景共产生六组结果。详细数量见上一级`README.md`，完整审计见`04_audit/stage2_A_ppi_union_audit.json`。

所有candidate和anchor均为canonical蛋白；HI-union为gene-level实验边，RF2为UniProt accession pair，并显式记录exact-accession与gene-propagated匹配模式。
"""
    safe_write(out / "ppi_screening/README.md", ppi_readme.encode("utf-8"))

    # Copy this exact builder into the release.
    safe_copy(Path(__file__).resolve(), out / "scripts/build_stage2_A_3418.py")

    manifest_rows: list[dict[str, object]] = []
    for path in inputs:
        manifest_rows.append(
            {
                "role": "INPUT",
                "path": str(path),
                "records": tsv_record_count(path) if tsv_record_count(path) is not None else "",
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "run_id": RUN_ID,
            }
        )
    for path in sorted(out.rglob("*")):
        if not path.is_file() or path.name.endswith(".partial") or path.name == "stage2_A_manifest.tsv":
            continue
        manifest_rows.append(
            {
                "role": "OUTPUT",
                "path": str(path.relative_to(out)),
                "records": tsv_record_count(path) if tsv_record_count(path) is not None else "",
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "run_id": RUN_ID,
            }
        )
    safe_write(
        out / "audit/stage2_A_manifest.tsv",
        render_tsv(["role", "path", "records", "bytes", "sha256", "run_id"], manifest_rows),
    )

    print(json.dumps(audit, ensure_ascii=False, indent=2))
    print(f"PASS\toutput={out}")


if __name__ == "__main__":
    main()
