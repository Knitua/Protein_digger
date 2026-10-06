#!/usr/bin/env python3
"""Positive-control recovery and threshold sensitivity after generic core freeze."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from stage25_v3_common import atomic_write_text, atomic_write_tsv, json_dump, read_tsv, sha256


IDENTITY_GRID = [0.70, 0.75, 0.80, 0.90]
COVERAGE_GRID = [0.70, 0.80, 0.90]


def verify_freeze(out: Path) -> dict:
    freeze_path = out / "audit" / "STAGE2_5_V3_CORE_FREEZE.json"
    freeze = json.loads(freeze_path.read_text())
    checks = {
        "method": sha256(out / "audit" / "stage2_5_v3_method.json") == freeze["method_sha256"],
        "summary": sha256(out / "audit" / "stage2_5_v3_core_summary.json") == freeze["core_summary_sha256"],
        "extended": sha256(out / "results" / "stage2_5_v3_extended_candidates.tsv") == freeze["extended_candidates_sha256"],
        "passing": sha256(out / "results" / "stage2_5_v3_passing_edges.tsv") == freeze["passing_edges_sha256"],
        "manifest": sha256(out / "audit" / "stage2_5_v3_manifest_pre_case_audit.tsv") == freeze["manifest_sha256"],
    }
    if not all(checks.values()):
        raise RuntimeError(f"generic core changed after freeze: {checks}")
    return checks


def func_profile_pass(row: dict, profile: str) -> bool:
    if profile == "PRIMARY_ANY_STRICT":
        return row["functional_pass_primary"] == "yes"
    if profile == "REACTOME_DIRECT_ONLY":
        return int(row["shared_reactome_direct_count"]) >= 1 and float(row["reactome_direct_jaccard"]) >= 0.50
    if profile == "GO_MF_DIRECT_ONLY":
        return int(row["shared_go_mf_direct_count"]) >= 2 and float(row["go_mf_direct_jaccard"]) >= 0.50
    if profile == "GO_BPMF_DIRECT_ONLY":
        return int(row["shared_go_bpmf_direct_count"]) >= 3 and float(row["go_bpmf_direct_idf_jaccard"]) >= 0.25
    raise ValueError(profile)


def edge_pass(row: dict, identity: float, coverage: float, profile: str = "PRIMARY_ANY_STRICT") -> bool:
    return (
        float(row["sequence_identity"]) >= identity
        and float(row["left_coverage"]) >= coverage
        and float(row["right_coverage"]) >= coverage
        and row["domain_architecture_pass"] == "yes"
        and func_profile_pass(row, profile)
    )


def recovered_accessions(edges: list[dict]) -> set[str]:
    return {row[key] for row in edges for key in ("left_accession", "right_accession")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--excluded-accessions", required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    freeze_checks = verify_freeze(out)
    excluded = {r["accession"] for r in read_tsv(args.excluded_accessions.resolve())}

    positive_all = list(read_tsv(out / "validation" / "positive_control_all_edges_v3.tsv"))
    positive = [
        row for row in positive_all
        if row["left_accession"] not in excluded and row["right_accession"] not in excluded
    ]
    core = list(read_tsv(out / "results" / "stage2_5_v3_all_core_edges.tsv"))
    seeds = list(read_tsv(out / "inputs" / "stage2_5_v3_final80_positive_seeds_1862.tsv"))
    evaluated_accessions = {r["seed_accession"] for r in seeds} - excluded
    aligned_eligible = recovered_accessions(positive)
    architecture_eligible = recovered_accessions([r for r in positive if r["domain_architecture_pass"] == "yes"])

    primary_candidates = {
        r["candidate_accession"] for r in read_tsv(out / "results" / "stage2_5_v3_extended_candidates.tsv")
    }
    positive_grid = []
    actual_grid = []
    recovered_primary_rows = []
    best_identity = defaultdict(float)
    for identity in IDENTITY_GRID:
        for coverage in COVERAGE_GRID:
            passed_positive = [r for r in positive if edge_pass(r, identity, coverage)]
            recovered = recovered_accessions(passed_positive)
            positive_grid.append(
                {
                    "sequence_identity_min": f"{identity:.2f}",
                    "bidirectional_coverage_min": f"{coverage:.2f}",
                    "functional_profile": "PRIMARY_ANY_STRICT",
                    "positive_controls_evaluated": len(evaluated_accessions),
                    "aligned_homology_eligible": len(aligned_eligible),
                    "strict_architecture_eligible": len(architecture_eligible),
                    "recovered_positive_controls": len(recovered),
                    "overall_recovery_fraction": f"{len(recovered) / len(evaluated_accessions):.8f}",
                    "conditional_recovery_fraction_among_architecture_eligible": (
                        f"{len(recovered) / len(architecture_eligible):.8f}" if architecture_eligible else "NA"
                    ),
                    "passing_positive_control_edges": len(passed_positive),
                }
            )
            passed_core = [r for r in core if edge_pass(r, identity, coverage)]
            candidates = {r["right_accession"] for r in passed_core}
            union = candidates | primary_candidates
            actual_grid.append(
                {
                    "sequence_identity_min": f"{identity:.2f}",
                    "bidirectional_coverage_min": f"{coverage:.2f}",
                    "functional_profile": "PRIMARY_ANY_STRICT",
                    "extended_candidate_count": len(candidates),
                    "passing_edge_count": len(passed_core),
                    "jaccard_vs_predeclared_primary_set": (
                        f"{len(candidates & primary_candidates) / len(union):.8f}" if union else "1.00000000"
                    ),
                }
            )
            if identity == 0.70 and coverage == 0.80:
                for row in passed_positive:
                    for acc, partner in (
                        (row["left_accession"], row["right_accession"]),
                        (row["right_accession"], row["left_accession"]),
                    ):
                        value = float(row["sequence_identity"])
                        if value > best_identity[acc]:
                            best_identity[acc] = value
                for acc in sorted(recovered):
                    recovered_primary_rows.append(
                        {"accession": acc, "best_supporting_sequence_identity": f"{best_identity[acc]:.8f}"}
                    )

    ablation_rows = []
    for profile in ("PRIMARY_ANY_STRICT", "REACTOME_DIRECT_ONLY", "GO_MF_DIRECT_ONLY", "GO_BPMF_DIRECT_ONLY"):
        passed = [r for r in positive if edge_pass(r, 0.70, 0.80, profile)]
        recovered = recovered_accessions(passed)
        ablation_rows.append(
            {
                "functional_profile": profile,
                "sequence_identity_min": "0.70",
                "bidirectional_coverage_min": "0.80",
                "recovered_positive_controls": len(recovered),
                "passing_positive_control_edges": len(passed),
                "overall_recovery_fraction": f"{len(recovered) / len(evaluated_accessions):.8f}",
            }
        )

    bins = Counter()
    for value in best_identity.values():
        if value >= 0.95:
            bins[">=0.95"] += 1
        elif value >= 0.90:
            bins["0.90-<0.95"] += 1
        elif value >= 0.80:
            bins["0.80-<0.90"] += 1
        else:
            bins["0.70-<0.80"] += 1
    identity_rows = [
        {"best_supporting_identity_bin": key, "recovered_positive_controls": bins.get(key, 0)}
        for key in ("0.70-<0.80", "0.80-<0.90", "0.90-<0.95", ">=0.95")
    ]

    atomic_write_tsv(out / "validation" / "positive_control_threshold_sensitivity_v3.tsv", positive_grid)
    atomic_write_tsv(out / "validation" / "actual_candidate_threshold_sensitivity_v3.tsv", actual_grid)
    atomic_write_tsv(out / "validation" / "positive_control_function_route_ablation_v3.tsv", ablation_rows)
    atomic_write_tsv(out / "validation" / "positive_control_recovered_primary_v3.tsv", recovered_primary_rows)
    atomic_write_tsv(out / "validation" / "positive_control_best_identity_distribution_v3.tsv", identity_rows)

    primary_row = next(
        r for r in positive_grid
        if r["sequence_identity_min"] == "0.70" and r["bidirectional_coverage_min"] == "0.80"
    )
    primary_actual = next(
        r for r in actual_grid
        if r["sequence_identity_min"] == "0.70" and r["bidirectional_coverage_min"] == "0.80"
    )
    report = f"""# Stage2.5 positive-control recovery and threshold sensitivity analysis

## 分析定位

本分析在通用Stage2.5核心结果冻结后运行。它不修改预先设定的主阈值，也不参与候选生成。验证时从阈值比较中排除预注册的case-family accessions；排除名单仅存在于验证目录，不进入核心筛选代码。

对每个Stage2 final80阳性蛋白，采用leave-one-out等价图分析：将该蛋白视为暂时隐藏的阳性，只允许其他阳性蛋白作为同源seed。若二者满足序列、domain architecture和严格功能规则，则视为该阳性可以被Stage2.5恢复。

## 主阈值结果

- 预设主阈值：identity ≥ 0.70，双向coverage ≥ 0.80；
- 参与评估的非case阳性：{len(evaluated_accessions)}；
- 存在预筛同源边的阳性：{len(aligned_eligible)}；
- 存在严格domain architecture匹配边的阳性：{len(architecture_eligible)}；
- 被完整Stage2.5规则恢复：{primary_row['recovered_positive_controls']}；
- 总体恢复比例：{float(primary_row['overall_recovery_fraction']):.2%}；
- 对严格architecture可评估阳性的条件恢复比例：{('NA' if primary_row['conditional_recovery_fraction_among_architecture_eligible'] == 'NA' else f"{float(primary_row['conditional_recovery_fraction_among_architecture_eligible']):.2%}")}；
- 同一主阈值下，真实Stage2失败池扩展候选：{primary_actual['extended_candidate_count']}。

## 可以解释的内容

- 比较70%、75%、80%和90% identity以及70%、80%和90% coverage时，阳性恢复数和实际扩展候选数如何变化；
- 收紧GOA evidence code和Reactome直接角色后，已知阳性是否仍可由同源成员恢复；
- 恢复结果主要来自近乎完全相同蛋白，还是也覆盖70%–90% identity区间；
- 预设70%/80%是否位于相对稳定区间。

## 不能解释的内容

- 本分析不能估计Stage2.5假阳性率、specificity或FPR；
- Stage2失败池不是实验确认的负样本；
- 被扩展的候选不因此获得直接PPI证明；
- 本分析不证明候选一定共享seed的同一个anchor。

完整阈值矩阵和功能证据消融见本目录TSV。主阈值保持预先设定值，不根据case结果或本敏感性分析回调。
"""
    atomic_write_text(out / "validation" / "positive_control_threshold_sensitivity_report.md", report)
    summary = {
        "analysis_name": "Stage2.5 positive-control recovery and threshold sensitivity analysis",
        "core_freeze_verified": all(freeze_checks.values()),
        "excluded_accession_count": len(excluded),
        "exclusion_file_sha256": sha256(args.excluded_accessions.resolve()),
        "primary_threshold_positive_controls_evaluated": len(evaluated_accessions),
        "primary_threshold_recovered_positive_controls": int(primary_row["recovered_positive_controls"]),
        "primary_threshold_actual_extended_candidates": int(primary_actual["extended_candidate_count"]),
        "threshold_grid_rows": len(positive_grid),
        "primary_thresholds_modified": False,
        "specificity_or_fpr_estimated": False,
    }
    json_dump(out / "validation" / "positive_control_threshold_sensitivity_summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
