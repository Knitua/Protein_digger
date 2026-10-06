#!/usr/bin/env python3
"""Shared, target-agnostic helpers for Stage 2.5 v3.

This module deliberately contains no case-specific gene or accession logic.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


RUN_ID = "stage2_5_homology_function_review_v3_20261006"

FINAL80_REL = {
    "A1": "04_stage2/A/ppi_screening/02_hi_union_rf2_final80/A1_WIDE_2877_hi_union_rf2_final80_screen_all.tsv",
    "B1": "04_stage2/B1/01_hi_union_rf2_final80/stage2_B1_1350_hi_union_rf2_final80_screen_all.tsv",
    "B2": "04_stage2/B2/01_hi_union_rf2_final80/stage2_B2_266_hi_union_rf2_final80_screen_all.tsv",
}
FINAL90_REL = {
    "A1": "04_stage2/A/ppi_screening/02_hi_union_rf2_final90/A1_WIDE_2877_hi_union_rf2_final90_screen_all.tsv",
    "B1": "04_stage2/B1/02_hi_union_rf2_final90/stage2_B1_1350_hi_union_rf2_final90_screen_all.tsv",
    "B2": "04_stage2/B2/02_hi_union_rf2_final90/stage2_B2_266_hi_union_rf2_final90_screen_all.tsv",
}
REFERENCE_REL = (
    "01_reference/human_swissprot_reviewed_frozen_v1/"
    "human_swissprot_reviewed_UP000005640_metadata.tsv"
)

SEQ_IDENTITY_PRIMARY = 0.70
SEED_COVERAGE_PRIMARY = 0.80
CANDIDATE_COVERAGE_PRIMARY = 0.80
LENGTH_RATIO_PREFILTER_MIN = 0.80
KMER_SIZE = 4
KMER_CONTAINMENT_PREFILTER_MIN = 0.08
DOMAIN_BOUNDARY_IOU_MIN = 0.80
REACTOME_JACCARD_MIN = 0.50
GO_MF_JACCARD_MIN = 0.50
GO_IDF_JACCARD_MIN = 0.25


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def atomic_write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(partial, path)


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    partial.write_text(text)
    os.replace(partial, path)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_ids(text: str, prefix: str) -> set[str]:
    patterns = {"PF": r"PF\d+", "IPR": r"IPR\d+", "GO": r"GO:\d+", "R-HSA": r"R-HSA-\d+"}
    return set(re.findall(patterns[prefix], text or ""))


def coarse_domain_signature(row: dict) -> tuple[str, str]:
    pfam = tuple(sorted(parse_ids(row.get("Pfam", ""), "PF")))
    interpro = tuple(sorted(parse_ids(row.get("InterPro", ""), "IPR")))
    if pfam and interpro:
        return "PFAM+INTERPRO", "PF:" + ",".join(pfam) + "|IPR:" + ",".join(interpro)
    if pfam:
        return "PFAM", "PF:" + ",".join(pfam)
    if interpro:
        return "INTERPRO", "IPR:" + ",".join(interpro)
    return "NONE", ""


def split_genes(text: str) -> list[str]:
    return [x.strip() for x in re.split(r"[;,|]", text or "") if x.strip()]


def kmer_set(seq: str, k: int = KMER_SIZE) -> set[str]:
    result = set()
    for idx in range(max(0, len(seq) - k + 1)):
        token = seq[idx : idx + k]
        if len(set(token)) >= 3:
            result.add(token)
    return result


def jaccard(a: set[str], b: set[str]) -> float:
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def weighted_jaccard(a: set[str], b: set[str], weights: dict[str, float]) -> float:
    union = a | b
    if not union:
        return 0.0
    return sum(weights.get(x, 1.0) for x in a & b) / sum(weights.get(x, 1.0) for x in union)


def load_screen(project: Path) -> tuple[dict[str, dict], list[Path]]:
    records: dict[str, dict] = {}
    inputs: list[Path] = []
    for route, rel in FINAL80_REL.items():
        path = project / rel
        inputs.append(path)
        for row in read_tsv(path):
            acc = row["candidate_accession"]
            if acc in records:
                raise RuntimeError(f"duplicate Stage2 candidate accession: {acc}")
            records[acc] = {
                "accession": acc,
                "gene": row.get("candidate_gene", ""),
                "stage1_route": route,
                "primary_final80_pass": row.get("ppi_pass", "").lower() == "yes",
                "primary_final80_category": row.get("ppi_evidence_category", ""),
                "primary_final80_anchor_count": row.get("merged_anchor_count", "0"),
                "primary_final80_anchor_genes": row.get("merged_anchor_genes", ""),
                "screen_source": str(path),
            }
    final90: dict[str, bool] = {}
    for _, rel in FINAL90_REL.items():
        path = project / rel
        inputs.append(path)
        for row in read_tsv(path):
            final90[row["candidate_accession"]] = row.get("ppi_pass", "").lower() == "yes"
    for acc, rec in records.items():
        rec["primary_final90_pass"] = final90.get(acc, False)
    return records, inputs


def load_reference(project: Path) -> tuple[dict[str, dict], Path]:
    path = project / REFERENCE_REL
    result = {}
    for row in read_tsv(path):
        result[row["Entry"]] = {
            "entry_name": row.get("Entry Name", ""),
            "gene_reference": row.get("Gene Names (primary)", ""),
            "length": int(row["Length"]),
            "sequence": row["Sequence"],
        }
    return result, path


def load_coarse_annotations(path: Path) -> dict[str, dict]:
    result = {}
    for row in read_tsv(path):
        basis, signature = coarse_domain_signature(row)
        result[row["Entry"]] = {
            "coarse_domain_basis": basis,
            "coarse_domain_signature": signature,
            "pfam_ids": parse_ids(row.get("Pfam", ""), "PF"),
            "interpro_ids": parse_ids(row.get("InterPro", ""), "IPR"),
            "uniprot_domain_ft": row.get("Domain [FT]", ""),
        }
    return result


def build_universe(project: Path, annotations_path: Path) -> tuple[list[dict], list[Path]]:
    screen, inputs = load_screen(project)
    reference, reference_path = load_reference(project)
    annotations = load_coarse_annotations(annotations_path)
    inputs.extend([reference_path, annotations_path])
    universe = []
    missing = []
    for acc, rec in screen.items():
        if acc not in reference:
            missing.append(acc)
            continue
        joined = dict(rec)
        joined.update(reference[acc])
        if acc in annotations:
            joined.update(annotations[acc])
            joined["annotation_status"] = "available"
        else:
            joined.update(
                {
                    "coarse_domain_basis": "NONE",
                    "coarse_domain_signature": "",
                    "pfam_ids": set(),
                    "interpro_ids": set(),
                    "uniprot_domain_ft": "",
                }
            )
            joined["annotation_status"] = "missing"
        universe.append(joined)
    if missing:
        raise RuntimeError(f"reference missing Stage2 accessions: {missing}")
    seeds = [x for x in universe if x["primary_final80_pass"]]
    review = [x for x in universe if not x["primary_final80_pass"]]
    if (len(universe), len(seeds), len(review)) != (4493, 1862, 2631):
        raise RuntimeError(
            f"unexpected frozen counts: universe={len(universe)}, seeds={len(seeds)}, review={len(review)}"
        )
    return universe, inputs


def load_architectures(path: Path) -> dict[str, dict]:
    result = {}
    for row in read_tsv(path):
        intervals = []
        if row.get("ordered_intervals"):
            for token in row["ordered_intervals"].split("|"):
                start, end = token.split("-")
                intervals.append((int(start), int(end)))
        result[row["accession"]] = {
            "source": row["architecture_source"],
            "status": row["architecture_status"],
            "length": int(row["protein_length"]) if row.get("protein_length") else 0,
            "ids": row.get("ordered_feature_ids", "").split("|") if row.get("ordered_feature_ids") else [],
            "intervals": intervals,
            "names": row.get("ordered_feature_names", "").split("|") if row.get("ordered_feature_names") else [],
        }
    return result


def normalized_interval_iou(a: tuple[int, int], alen: int, b: tuple[int, int], blen: int) -> float:
    if alen <= 0 or blen <= 0:
        return 0.0
    a0, a1 = (a[0] - 1) / alen, a[1] / alen
    b0, b1 = (b[0] - 1) / blen, b[1] / blen
    intersection = max(0.0, min(a1, b1) - max(a0, b0))
    union = max(a1, b1) - min(a0, b0)
    return intersection / union if union > 0 else 0.0


def compare_architectures(a: dict | None, b: dict | None) -> dict:
    result = {
        "domain_architecture_available": False,
        "domain_source_match": False,
        "domain_order_copy_exact": False,
        "domain_boundary_min_normalized_iou": 0.0,
        "domain_boundary_pass": False,
        "domain_architecture_pass": False,
        "domain_architecture_reason": "missing_architecture",
    }
    if not a or not b or a["status"] != "available" or b["status"] != "available":
        return result
    if not a["ids"] or not b["ids"]:
        result["domain_architecture_reason"] = "empty_architecture"
        return result
    result["domain_architecture_available"] = True
    result["domain_source_match"] = a["source"] == b["source"]
    if not result["domain_source_match"]:
        result["domain_architecture_reason"] = "architecture_source_mismatch"
        return result
    result["domain_order_copy_exact"] = a["ids"] == b["ids"]
    if not result["domain_order_copy_exact"]:
        result["domain_architecture_reason"] = "domain_order_or_copy_mismatch"
        return result
    if len(a["intervals"]) != len(a["ids"]) or len(b["intervals"]) != len(b["ids"]):
        result["domain_architecture_reason"] = "coordinate_count_mismatch"
        return result
    ious = [
        normalized_interval_iou(ai, a["length"], bi, b["length"])
        for ai, bi in zip(a["intervals"], b["intervals"])
    ]
    result["domain_boundary_min_normalized_iou"] = min(ious) if ious else 0.0
    result["domain_boundary_pass"] = result["domain_boundary_min_normalized_iou"] >= DOMAIN_BOUNDARY_IOU_MIN
    result["domain_architecture_pass"] = result["domain_boundary_pass"]
    result["domain_architecture_reason"] = "pass" if result["domain_architecture_pass"] else "domain_boundary_mismatch"
    return result


def load_go_direct(path: Path) -> dict[str, dict[str, set[str]]]:
    result: dict[str, dict[str, set[str]]] = defaultdict(lambda: {"P": set(), "F": set()})
    for row in read_tsv(path):
        if row.get("included_in_primary") == "yes":
            result[row["accession"]][row["aspect"]].add(row["go_id"])
    return result


def load_reactome_direct(path: Path) -> dict[str, set[str]]:
    result: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(path):
        if row.get("participant_directness") == "DIRECT_ENTITY" and row.get("inferred") == "false":
            result[row["accession"]].add(row["event_id"])
    return result


def build_go_weights(go_by_acc: dict[str, dict[str, set[str]]], accessions: list[str]) -> dict[str, float]:
    term_sets = [go_by_acc.get(acc, {"P": set(), "F": set()})["P"] | go_by_acc.get(acc, {"P": set(), "F": set()})["F"] for acc in accessions]
    df = Counter(term for terms in term_sets for term in terms)
    return {term: math.log((len(term_sets) + 1) / (count + 1)) + 1.0 for term, count in df.items()}


def functional_metrics(
    a: str,
    b: str,
    go_by_acc: dict[str, dict[str, set[str]]],
    reactome_by_acc: dict[str, set[str]],
    weights: dict[str, float],
) -> dict:
    ago = go_by_acc.get(a, {"P": set(), "F": set()})
    bgo = go_by_acc.get(b, {"P": set(), "F": set()})
    abp, amf = ago["P"], ago["F"]
    bbp, bmf = bgo["P"], bgo["F"]
    aall, ball = abp | amf, bbp | bmf
    ar, br = reactome_by_acc.get(a, set()), reactome_by_acc.get(b, set())
    return {
        "shared_reactome_direct_count": len(ar & br),
        "reactome_direct_jaccard": jaccard(ar, br),
        "shared_go_bp_direct_count": len(abp & bbp),
        "go_bp_direct_jaccard": jaccard(abp, bbp),
        "shared_go_mf_direct_count": len(amf & bmf),
        "go_mf_direct_jaccard": jaccard(amf, bmf),
        "shared_go_bpmf_direct_count": len(aall & ball),
        "go_bpmf_direct_idf_jaccard": weighted_jaccard(aall, ball, weights),
    }


def functional_pass(metrics: dict) -> tuple[bool, str]:
    if metrics["shared_reactome_direct_count"] >= 1 and metrics["reactome_direct_jaccard"] >= REACTOME_JACCARD_MIN:
        return True, "reactome_direct_jaccard_ge_0p50"
    if metrics["shared_go_mf_direct_count"] >= 2 and metrics["go_mf_direct_jaccard"] >= GO_MF_JACCARD_MIN:
        return True, "go_mf_direct_jaccard_ge_0p50_with_ge2_shared"
    if metrics["shared_go_bpmf_direct_count"] >= 3 and metrics["go_bpmf_direct_idf_jaccard"] >= GO_IDF_JACCARD_MIN:
        return True, "go_bpmf_direct_idf_jaccard_ge_0p25_with_ge3_shared"
    return False, "strict_functional_similarity_below_primary_rule"


def json_dump(path: Path, obj: object) -> None:
    atomic_write_text(path, json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
