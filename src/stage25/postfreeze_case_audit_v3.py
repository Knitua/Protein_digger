#!/usr/bin/env python3
"""Audit named cases only after the generic core and sensitivity outputs are frozen."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from stage25_v3_common import atomic_write_text, atomic_write_tsv, read_tsv, sha256


def verify_freeze(out: Path) -> None:
    freeze = json.loads((out / "audit" / "STAGE2_5_V3_CORE_FREEZE.json").read_text())
    checks = [
        sha256(out / "audit" / "stage2_5_v3_method.json") == freeze["method_sha256"],
        sha256(out / "audit" / "stage2_5_v3_core_summary.json") == freeze["core_summary_sha256"],
        sha256(out / "results" / "stage2_5_v3_extended_candidates.tsv") == freeze["extended_candidates_sha256"],
        sha256(out / "results" / "stage2_5_v3_passing_edges.tsv") == freeze["passing_edges_sha256"],
        sha256(out / "audit" / "stage2_5_v3_manifest_pre_case_audit.tsv") == freeze["manifest_sha256"],
    ]
    if not all(checks):
        raise RuntimeError("core freeze verification failed before case audit")
    if not (out / "validation" / "positive_control_threshold_sensitivity_summary.json").exists():
        raise RuntimeError("positive-control sensitivity must finish before case audit")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    verify_freeze(out)

    targets = list(read_tsv(args.targets.resolve()))
    target_by_acc = {r["accession"]: r for r in targets}
    seeds = {r["seed_accession"]: r for r in read_tsv(out / "inputs" / "stage2_5_v3_final80_positive_seeds_1862.tsv")}
    failures = {r["candidate_accession"]: r for r in read_tsv(out / "inputs" / "stage2_5_v3_final80_failed_review_pool_2631.tsv")}
    extended = {r["candidate_accession"]: r for r in read_tsv(out / "results" / "stage2_5_v3_extended_candidates.tsv")}
    all_edges = defaultdict(list)
    for row in read_tsv(out / "results" / "stage2_5_v3_all_core_edges.tsv"):
        if row["right_accession"] in target_by_acc:
            all_edges[row["right_accession"]].append(row)

    rows = []
    for target in targets:
        acc = target["accession"]
        gene = target.get("gene", "")
        if acc in seeds:
            seed = seeds[acc]
            rows.append(
                {
                    "accession": acc,
                    "gene": gene,
                    "primary_stage2_final80_status": "pass",
                    "primary_stage2_category": seed["primary_final80_category"],
                    "stage2_5_v3_status": "not_applicable_primary_positive_seed",
                    "best_seed_accession": "",
                    "best_seed_gene": "",
                    "sequence_identity": "",
                    "seed_coverage": "",
                    "candidate_coverage": "",
                    "domain_architecture_pass": "",
                    "functional_pass": "",
                    "functional_pass_reason": "",
                    "interpretation": "primary_Stage2_positive_seed",
                }
            )
            continue
        ext = extended.get(acc)
        if ext:
            rows.append(
                {
                    "accession": acc,
                    "gene": gene,
                    "primary_stage2_final80_status": "fail",
                    "primary_stage2_category": failures.get(acc, {}).get("primary_final80_category", ""),
                    "stage2_5_v3_status": "extended_candidate",
                    "best_seed_accession": ext["best_seed_accession"],
                    "best_seed_gene": ext["best_seed_gene"],
                    "sequence_identity": ext["sequence_identity"],
                    "seed_coverage": ext["seed_coverage"],
                    "candidate_coverage": ext["candidate_coverage"],
                    "domain_architecture_pass": "yes",
                    "functional_pass": "yes",
                    "functional_pass_reason": ext["functional_pass_reason"],
                    "interpretation": "postfreeze_homology_domain_function_extension_not_direct_PPI_proof",
                }
            )
            continue
        edges = all_edges.get(acc, [])
        best = max(
            edges,
            default=None,
            key=lambda r: (
                r["sequence_pass_primary"] == "yes",
                r["domain_architecture_pass"] == "yes",
                r["functional_pass_primary"] == "yes",
                float(r["sequence_identity"]),
            ),
        )
        rows.append(
            {
                "accession": acc,
                "gene": gene,
                "primary_stage2_final80_status": "fail" if acc in failures else "not_in_stage2_universe",
                "primary_stage2_category": failures.get(acc, {}).get("primary_final80_category", ""),
                "stage2_5_v3_status": "not_extended",
                "best_seed_accession": best["left_accession"] if best else "",
                "best_seed_gene": best["left_gene"] if best else "",
                "sequence_identity": best["sequence_identity"] if best else "",
                "seed_coverage": best["left_coverage"] if best else "",
                "candidate_coverage": best["right_coverage"] if best else "",
                "domain_architecture_pass": best["domain_architecture_pass"] if best else "",
                "functional_pass": best["functional_pass_primary"] if best else "",
                "functional_pass_reason": best["functional_pass_reason"] if best else "no_prefilter_edge",
                "interpretation": "not_extended_by_frozen_v3_rule",
            }
        )

    rows.sort(key=lambda r: r["gene"])
    atomic_write_tsv(out / "postfreeze_case_audit" / "RAS_postfreeze_audit_v3.tsv", rows)
    lines = [
        "# RAS post-freeze case audit",
        "",
        "该审计在通用核心结果及非case阳性阈值敏感性分析完成后运行。目标身份未进入核心筛选、阈值、功能规则或核心验收。",
        "",
        "| Gene | Accession | Primary Stage2 | Stage2.5 v3 | Supporting seed | Identity | 解释 |",
        "|---|---|---|---|---|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['gene']} | {row['accession']} | {row['primary_stage2_final80_status']} | "
            f"{row['stage2_5_v3_status']} | {row['best_seed_gene'] or '—'} | "
            f"{row['sequence_identity'] or '—'} | {row['interpretation']} |"
        )
    lines.extend(
        [
            "",
            "Stage2.5通过仅表示序列同源、domain architecture一致且存在严格功能证据支持；不等于已证明与特定anchor发生直接PPI。",
            "",
            f"Target manifest SHA256: `{sha256(args.targets.resolve())}`",
        ]
    )
    atomic_write_text(out / "postfreeze_case_audit" / "RAS_postfreeze_audit_v3.md", "\n".join(lines) + "\n")
    print(json.dumps({"audited_targets": len(rows), "extended": sum(r["stage2_5_v3_status"] == "extended_candidate" for r in rows)}, indent=2))


if __name__ == "__main__":
    main()
