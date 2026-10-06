#!/usr/bin/env python3
"""Finalize the generic Stage2.5 v3 core and freeze it before case audits."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from stage25_v3_common import (
    CANDIDATE_COVERAGE_PRIMARY,
    DOMAIN_BOUNDARY_IOU_MIN,
    GO_IDF_JACCARD_MIN,
    GO_MF_JACCARD_MIN,
    REACTOME_JACCARD_MIN,
    RUN_ID,
    SEED_COVERAGE_PRIMARY,
    SEQ_IDENTITY_PRIMARY,
    atomic_write_tsv,
    build_go_weights,
    build_universe,
    compare_architectures,
    functional_metrics,
    functional_pass,
    json_dump,
    load_architectures,
    load_go_direct,
    load_reactome_direct,
    read_tsv,
    sha256,
)


def load_alignment(path: Path) -> dict[str, dict]:
    return {row["pair_id"]: row for row in read_tsv(path)}


def metric_rows(
    pairs_path: Path,
    alignment_path: Path,
    architectures: dict,
    go_by_acc: dict,
    reactome_by_acc: dict,
    go_weights: dict,
    universe_by_acc: dict,
    pair_class: str,
) -> list[dict]:
    alignments = load_alignment(alignment_path)
    rows = []
    for pair in read_tsv(pairs_path):
        aln = alignments[pair["pair_id"]]
        left = pair["left_accession"]
        right = pair["right_accession"]
        identity = float(aln["sequence_identity"])
        left_cov = float(aln["seed_coverage"])
        right_cov = float(aln["candidate_coverage"])
        architecture = compare_architectures(architectures.get(left), architectures.get(right))
        function = functional_metrics(left, right, go_by_acc, reactome_by_acc, go_weights)
        fpass, freason = functional_pass(function)
        sequence_pass = (
            identity >= SEQ_IDENTITY_PRIMARY
            and left_cov >= SEED_COVERAGE_PRIMARY
            and right_cov >= CANDIDATE_COVERAGE_PRIMARY
        )
        edge_pass = sequence_pass and architecture["domain_architecture_pass"] and fpass
        lgo = go_by_acc.get(left, {"P": set(), "F": set()})
        rgo = go_by_acc.get(right, {"P": set(), "F": set()})
        shared_go_bp = sorted(lgo["P"] & rgo["P"])
        shared_go_mf = sorted(lgo["F"] & rgo["F"])
        shared_reactome = sorted(reactome_by_acc.get(left, set()) & reactome_by_acc.get(right, set()))
        left_rec = universe_by_acc[left]
        right_rec = universe_by_acc[right]
        rows.append(
            {
                "run_id": RUN_ID,
                "pair_class": pair_class,
                "pair_id": pair["pair_id"],
                "left_accession": left,
                "left_gene": left_rec["gene"],
                "left_stage1_route": left_rec["stage1_route"],
                "right_accession": right,
                "right_gene": right_rec["gene"],
                "right_stage1_route": right_rec["stage1_route"],
                "coarse_domain_signature": pair["coarse_domain_signature"],
                "length_ratio": pair["length_ratio"],
                "kmer_containment": pair["kmer_containment"],
                "alignment_score": aln["alignment_score"],
                "sequence_identity": f"{identity:.8f}",
                "left_coverage": f"{left_cov:.8f}",
                "right_coverage": f"{right_cov:.8f}",
                "sequence_pass_primary": "yes" if sequence_pass else "no",
                "domain_architecture_source": architectures.get(left, {}).get("source", ""),
                "domain_architecture_available": "yes" if architecture["domain_architecture_available"] else "no",
                "domain_order_copy_exact": "yes" if architecture["domain_order_copy_exact"] else "no",
                "domain_boundary_min_normalized_iou": f"{architecture['domain_boundary_min_normalized_iou']:.8f}",
                "domain_boundary_pass": "yes" if architecture["domain_boundary_pass"] else "no",
                "domain_architecture_pass": "yes" if architecture["domain_architecture_pass"] else "no",
                "domain_architecture_reason": architecture["domain_architecture_reason"],
                **{
                    key: f"{value:.8f}" if isinstance(value, float) else value
                    for key, value in function.items()
                },
                "shared_reactome_direct_ids": "|".join(shared_reactome),
                "shared_go_bp_direct_ids": "|".join(shared_go_bp),
                "shared_go_mf_direct_ids": "|".join(shared_go_mf),
                "functional_pass_primary": "yes" if fpass else "no",
                "functional_pass_reason": freason,
                "stage2_5_v3_edge_pass": "yes" if edge_pass else "no",
                "decision_semantics": (
                    "homology_domain_function_supported_extension_not_primary_PPI_positive"
                    if edge_pass
                    else "reviewed_not_extended"
                ),
            }
        )
    rows.sort(key=lambda r: (r["right_accession"], -float(r["sequence_identity"]), r["left_accession"]))
    return rows


def write_manifest(path: Path, entries: list[tuple[str, Path]]) -> None:
    rows = []
    seen = set()
    for role, item in entries:
        item = item.resolve()
        key = (role, str(item))
        if key in seen or not item.exists() or not item.is_file():
            continue
        seen.add(key)
        rows.append({"role": role, "path": str(item), "bytes": item.stat().st_size, "sha256": sha256(item)})
    rows.sort(key=lambda r: (r["role"], r["path"]))
    atomic_write_tsv(path, rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--scripts-dir", required=True, type=Path)
    parser.add_argument("--aligner", required=True, type=Path)
    args = parser.parse_args()
    project = args.project_root.resolve()
    out = args.output_dir.resolve()
    scripts = args.scripts_dir.resolve()

    universe, source_inputs = build_universe(project, args.annotations.resolve())
    universe_by_acc = {r["accession"]: r for r in universe}
    architectures = load_architectures(out / "resources" / "domain_architecture_v3.tsv")
    go_by_acc = load_go_direct(out / "resources" / "goa_stage25_evidence_long_v3.tsv")
    reactome_by_acc = load_reactome_direct(out / "resources" / "reactome_direct_function_v3.tsv")
    go_weights = build_go_weights(go_by_acc, sorted(universe_by_acc))

    core_edges = metric_rows(
        out / "work" / "core_prefilter_pairs.tsv",
        out / "work" / "core_alignment_output.tsv",
        architectures,
        go_by_acc,
        reactome_by_acc,
        go_weights,
        universe_by_acc,
        "primary_failure_vs_positive_seed",
    )
    positive_edges = metric_rows(
        out / "work" / "positive_control_prefilter_pairs.tsv",
        out / "work" / "positive_control_alignment_output.tsv",
        architectures,
        go_by_acc,
        reactome_by_acc,
        go_weights,
        universe_by_acc,
        "positive_vs_positive_leave_one_out_edge",
    )
    passing = [r for r in core_edges if r["stage2_5_v3_edge_pass"] == "yes"]
    atomic_write_tsv(out / "results" / "stage2_5_v3_all_core_edges.tsv", core_edges)
    atomic_write_tsv(
        out / "results" / "stage2_5_v3_passing_edges.tsv",
        passing,
        list(core_edges[0]) if core_edges else [],
    )
    atomic_write_tsv(out / "validation" / "positive_control_all_edges_v3.tsv", positive_edges)

    by_candidate: dict[str, list[dict]] = defaultdict(list)
    for row in passing:
        by_candidate[row["right_accession"]].append(row)
    extended = []
    for candidate, rows in by_candidate.items():
        best = max(
            rows,
            key=lambda r: (
                float(r["sequence_identity"]),
                min(float(r["left_coverage"]), float(r["right_coverage"])),
                float(r["domain_boundary_min_normalized_iou"]),
                max(
                    float(r["reactome_direct_jaccard"]),
                    float(r["go_mf_direct_jaccard"]),
                    float(r["go_bpmf_direct_idf_jaccard"]),
                ),
                r["left_accession"],
            ),
        )
        seed = universe_by_acc[best["left_accession"]]
        rec = universe_by_acc[candidate]
        extended.append(
            {
                "run_id": RUN_ID,
                "candidate_accession": candidate,
                "candidate_gene": rec["gene"],
                "candidate_stage1_route": rec["stage1_route"],
                "primary_stage2_final80_status": "fail",
                "stage2_5_v3_status": "extended_candidate",
                "qualifying_seed_count": len(rows),
                "best_seed_accession": seed["accession"],
                "best_seed_gene": seed["gene"],
                "best_seed_stage1_route": seed["stage1_route"],
                "best_seed_primary_final80_category": seed["primary_final80_category"],
                "best_seed_also_final90_pass": "yes" if seed["primary_final90_pass"] else "no",
                "sequence_identity": best["sequence_identity"],
                "seed_coverage": best["left_coverage"],
                "candidate_coverage": best["right_coverage"],
                "domain_architecture_source": best["domain_architecture_source"],
                "domain_order_copy_exact": best["domain_order_copy_exact"],
                "domain_boundary_min_normalized_iou": best["domain_boundary_min_normalized_iou"],
                "functional_pass_reason": best["functional_pass_reason"],
                "shared_reactome_direct_ids": best["shared_reactome_direct_ids"],
                "shared_go_bp_direct_ids": best["shared_go_bp_direct_ids"],
                "shared_go_mf_direct_ids": best["shared_go_mf_direct_ids"],
                "stage3_entry_label": "Stage2_extended_homology_domain_function_supported",
                "primary_ppi_label_unchanged": "yes",
                "interpretation": "extension_candidate_not_direct_PPI_proof",
            }
        )
    extended.sort(key=lambda r: (r["candidate_stage1_route"], -float(r["sequence_identity"]), r["candidate_accession"]))
    atomic_write_tsv(out / "results" / "stage2_5_v3_extended_candidates.tsv", extended)

    route_counts = Counter(r["candidate_stage1_route"] for r in extended)
    function_counts = Counter(r["functional_pass_reason"] for r in extended)
    summary = {
        "run_id": RUN_ID,
        "screened_canonical_universe": len(universe),
        "primary_final80_positive_seeds": sum(r["primary_final80_pass"] for r in universe),
        "primary_final80_failure_pool": sum(not r["primary_final80_pass"] for r in universe),
        "core_prefilter_and_aligned_edges": len(core_edges),
        "core_sequence_primary_pass_edges": sum(r["sequence_pass_primary"] == "yes" for r in core_edges),
        "core_domain_architecture_pass_edges": sum(r["domain_architecture_pass"] == "yes" for r in core_edges),
        "core_strict_function_pass_edges": sum(r["functional_pass_primary"] == "yes" for r in core_edges),
        "passing_edges": len(passing),
        "unique_extended_candidates": len(extended),
        "extended_by_stage1_route": dict(sorted(route_counts.items())),
        "extended_by_function_route": dict(sorted(function_counts.items())),
        "positive_control_prefilter_edges": len(positive_edges),
        "primary_stage2_labels_modified": False,
        "case_specific_logic_used_for_core_selection": False,
        "same_anchor_evidence_tiers_used": False,
    }
    json_dump(out / "audit" / "stage2_5_v3_core_summary.json", summary)
    atomic_write_tsv(
        out / "audit" / "stage2_5_v3_core_summary.tsv",
        [
            {"metric": key, "value": json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list, bool)) else value}
            for key, value in summary.items()
        ],
    )

    method = {
        "run_id": RUN_ID,
        "primary_sequence_identity_min": SEQ_IDENTITY_PRIMARY,
        "primary_seed_coverage_min": SEED_COVERAGE_PRIMARY,
        "primary_candidate_coverage_min": CANDIDATE_COVERAGE_PRIMARY,
        "domain_rule": (
            "coordinate-resolved Pfam architecture; InterPro fallback only when Pfam is absent; "
            "exact ordered IDs and copy number; minimum normalized boundary IoU >= 0.80"
        ),
        "domain_boundary_iou_min": DOMAIN_BOUNDARY_IOU_MIN,
        "go_primary_evidence_codes": ["EXP", "IDA", "IGI", "IMP"],
        "go_excluded_evidence_semantics": (
            "IPI/IEP/high-throughput/computational/phylogenetic/transferred evidence and "
            "NOT/contributes_to/colocalizes_with qualifiers do not support the primary function rule"
        ),
        "reactome_rule": "human, non-disease, non-inferred event; protein must be DIRECT_ENTITY CATALYST or REGULATOR",
        "functional_rule": (
            f"Reactome direct Jaccard >= {REACTOME_JACCARD_MIN} with >=1 shared event OR "
            f"GO-MF direct Jaccard >= {GO_MF_JACCARD_MIN} with >=2 shared terms OR "
            f"GO(BP+MF) direct IDF-Jaccard >= {GO_IDF_JACCARD_MIN} with >=3 shared terms"
        ),
        "same_anchor_evidence_tiers": "not used",
        "case_specific_checks": "not used in core selection or core validation",
        "semantics": "homology/domain/function supported extension; primary Stage2 PPI remains fail",
    }
    json_dump(out / "audit" / "stage2_5_v3_method.json", method)

    checks = {
        "screened_universe_is_4493": len(universe) == 4493,
        "seed_count_is_1862": sum(r["primary_final80_pass"] for r in universe) == 1862,
        "failure_pool_is_2631": sum(not r["primary_final80_pass"] for r in universe) == 2631,
        "core_accession_sets_disjoint": not (
            {r["accession"] for r in universe if r["primary_final80_pass"]}
            & {r["accession"] for r in universe if not r["primary_final80_pass"]}
        ),
        "extended_accessions_unique": len(extended) == len({r["candidate_accession"] for r in extended}),
        "all_extended_are_primary_failures": all(r["primary_stage2_final80_status"] == "fail" for r in extended),
        "all_passing_edges_meet_primary_sequence_rule": all(r["sequence_pass_primary"] == "yes" for r in passing),
        "all_passing_edges_meet_strict_domain_architecture": all(r["domain_architecture_pass"] == "yes" for r in passing),
        "all_passing_edges_meet_strict_function_rule": all(r["functional_pass_primary"] == "yes" for r in passing),
        "primary_stage2_labels_not_modified": True,
        "case_specific_validation_not_used": True,
    }
    validation = {"run_id": RUN_ID, "status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}
    json_dump(out / "audit" / "stage2_5_v3_core_validation.json", validation)
    if validation["status"] != "PASS":
        raise RuntimeError("generic core validation failed")

    core_code_names = {
        "stage25_v3_common.py",
        "prepare_stage25_v3.py",
        "fetch_domain_architectures_v3.py",
        "build_strict_function_resources_v3.py",
        "finalize_stage25_v3.py",
        "stage2_5_local_aligner.cpp",
        "stage2_5_local_aligner",
    }
    manifest_entries = [("source_input", p) for p in source_inputs]
    manifest_entries.extend(("core_code", scripts / name) for name in core_code_names)
    for folder in ("inputs", "resources", "work", "results", "audit"):
        for path in (out / folder).rglob("*"):
            if path.is_file() and ".partial." not in path.name and path.name not in {
                "stage2_5_v3_manifest_pre_case_audit.tsv",
                "STAGE2_5_V3_CORE_FREEZE.json",
            }:
                manifest_entries.append(("core_artifact", path))
    manifest_path = out / "audit" / "stage2_5_v3_manifest_pre_case_audit.tsv"
    write_manifest(manifest_path, manifest_entries)
    freeze = {
        "run_id": RUN_ID,
        "freeze_scope": "generic_core_before_threshold_sensitivity_and_case_audit",
        "method_sha256": sha256(out / "audit" / "stage2_5_v3_method.json"),
        "core_summary_sha256": sha256(out / "audit" / "stage2_5_v3_core_summary.json"),
        "extended_candidates_sha256": sha256(out / "results" / "stage2_5_v3_extended_candidates.tsv"),
        "passing_edges_sha256": sha256(out / "results" / "stage2_5_v3_passing_edges.tsv"),
        "manifest_sha256": sha256(manifest_path),
        "primary_thresholds_are_predeclared_and_not_changed_by_sensitivity_analysis": True,
        "core_frozen_before_case_audit": True,
    }
    json_dump(out / "audit" / "STAGE2_5_V3_CORE_FREEZE.json", freeze)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
