#!/usr/bin/env python3
"""Screen the current formal Stage1-B2 set against HI-union and RF2-PPI."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path("/root/autodl-tmp/Agent_analysis_v2")
OUTPUT_ROOT = ROOT / "04_stage2" / "B2"
SHARED = ROOT / "04_stage2" / "A" / "ppi_screening"
ANCHORS = ROOT / "02_anchors" / "working_regulatory_anchor_2368.tsv"
HI_EDGES = SHARED / "01_common_hi_union" / "hi_union_current2368_anchor_partner_edges.tsv"
RF2_EDGES = {
    "final80": SHARED / "02_hi_union_rf2_final80" / "rf2_final80_current2368_anchor_partner_edges.tsv",
    "final90": SHARED / "02_hi_union_rf2_final90" / "rf2_final90_current2368_anchor_partner_edges.tsv",
}
ROUTES = {
    "B2": {
        "input": ROOT / "03_stage1_B" / "B2" / "04_formal_after_B1" / "b2_formal_top2pct_266.tsv",
        "expected": 266,
        "definition": "B1-excluded equal-weight four-percentile functional-network ranking Top 2%",
    },
}
RUN_ID = "stage2_B2_current266_hi_union_rf2_official_precision_20260904"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        return list(reader.fieldnames), rows


def genes(value: str | None) -> set[str]:
    return {
        token.strip().upper()
        for token in re.split(r"[;|]", (value or "").strip())
        if token.strip()
    }


def identifiers(value: str | None) -> set[str]:
    return {
        token.strip()
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
        raise FileExistsError(f"Refusing to overwrite different file: {path}")
    partial = path.with_name(path.name + ".partial")
    if partial.exists():
        raise FileExistsError(f"Partial file exists: {partial}")
    with partial.open("xb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def count_genes(rows: Iterable[dict[str, object]]) -> int:
    result: set[str] = set()
    for row in rows:
        result |= genes(str(row.get("candidate_gene", "")))
    return len(result)


def evidence_category(hi: set[str], rf2: set[str], label: str) -> str:
    if hi and rf2:
        return f"HI-union+RF2-{label}"
    if hi:
        return "HI-union-only"
    if rf2:
        return f"RF2-{label}-only"
    return "NO_PPI_DATASET_SUPPORT"


def main() -> None:
    for path in [ANCHORS, HI_EDGES, *RF2_EDGES.values()]:
        if not path.is_file():
            raise FileNotFoundError(path)

    _, anchor_rows = read_tsv(ANCHORS)
    if len(anchor_rows) != 2368:
        raise ValueError("Current anchor table must contain 2,368 rows")
    anchor_accessions = {row["anchor_accession"].strip() for row in anchor_rows}
    anchor_genes: set[str] = set()
    anchor_accession_to_gene: dict[str, str] = {}
    for row in anchor_rows:
        anchor_genes |= genes(row.get("anchor_gene"))
        accession = row["anchor_accession"].strip()
        gene = row.get("anchor_gene", "").strip().upper()
        if accession in anchor_accession_to_gene and anchor_accession_to_gene[accession] != gene:
            raise ValueError(f"Conflicting anchor gene mapping for {accession}")
        anchor_accession_to_gene[accession] = gene

    _, hi_rows = read_tsv(HI_EDGES)
    hi_partner_to_anchors: dict[str, set[str]] = defaultdict(set)
    hi_partner_to_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in hi_rows:
        partner = row["partner_gene"].strip().upper()
        anchor = row["anchor_gene"].strip().upper()
        if partner and anchor:
            hi_partner_to_anchors[partner].add(anchor)
            hi_partner_to_rows[partner].append(row)

    rf2_accession_to_anchors: dict[str, dict[str, set[str]]] = {}
    rf2_gene_to_anchors: dict[str, dict[str, set[str]]] = {}
    rf2_accession_to_rows: dict[str, dict[str, list[dict[str, str]]]] = {}
    rf2_gene_to_rows: dict[str, dict[str, list[dict[str, str]]]] = {}
    for label, path in RF2_EDGES.items():
        _, rows = read_tsv(path)
        by_accession: dict[str, set[str]] = defaultdict(set)
        by_gene: dict[str, set[str]] = defaultdict(set)
        rows_by_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
        rows_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            partner_accession = row["partner_accession"].strip()
            current_anchors = genes(row.get("anchor_genes"))
            if partner_accession:
                by_accession[partner_accession] |= current_anchors
                rows_by_accession[partner_accession].append(row)
            for gene in genes(row.get("partner_genes")):
                by_gene[gene] |= current_anchors
                rows_by_gene[gene].append(row)
        rf2_accession_to_anchors[label] = by_accession
        rf2_gene_to_anchors[label] = by_gene
        rf2_accession_to_rows[label] = rows_by_accession
        rf2_gene_to_rows[label] = rows_by_gene

    route_accessions: dict[str, set[str]] = {}
    route_results: dict[str, dict[str, object]] = {}
    for route, spec in ROUTES.items():
        input_path = spec["input"]
        expected = int(spec["expected"])
        if not input_path.is_file():
            raise FileNotFoundError(input_path)
        source_fields, source_rows = read_tsv(input_path)
        if len(source_rows) != expected:
            raise ValueError(f"{route}: expected {expected} rows, found {len(source_rows)}")
        accessions = [row["accession"].strip() for row in source_rows]
        if len(set(accessions)) != expected or "" in accessions:
            raise ValueError(f"{route}: duplicate or blank accession")
        candidate_genes = {accession: genes(row.get("gene")) for accession, row in zip(accessions, source_rows)}
        if set(accessions) & anchor_accessions:
            raise ValueError(f"{route}: candidate accession overlaps current anchors")
        if any(value & anchor_genes for value in candidate_genes.values()):
            raise ValueError(f"{route}: candidate gene overlaps current anchors")
        route_accessions[route] = set(accessions)

        out_root = OUTPUT_ROOT
        result_fields = [
            "candidate_accession", "candidate_gene", "stage1_route",
            "ppi_screen_version", "ppi_pass", "ppi_evidence_category",
            "merged_anchor_count", "merged_anchor_genes",
            "hi_union_anchor_count", "hi_union_anchor_genes",
            "rf2_anchor_count", "rf2_anchor_genes",
            "rf2_exact_accession_anchor_count", "rf2_exact_accession_anchor_genes",
            "rf2_gene_propagated_anchor_count", "rf2_gene_propagated_anchor_genes",
            "rf2_match_mode", "ppi_run_id",
        ] + source_fields
        route_outputs: list[Path] = []
        threshold_data: dict[str, dict[str, object]] = {}

        for label in ["final80", "final90"]:
            all_rows: list[dict[str, object]] = []
            passed_rows: list[dict[str, object]] = []
            pair_rows: list[dict[str, object]] = []
            for source in source_rows:
                accession = source["accession"].strip()
                current_genes = candidate_genes[accession]
                hi_anchors: set[str] = set()
                rf2_gene_anchors: set[str] = set()
                for gene in current_genes:
                    hi_anchors |= hi_partner_to_anchors.get(gene, set())
                    rf2_gene_anchors |= rf2_gene_to_anchors[label].get(gene, set())
                rf2_exact_anchors = rf2_accession_to_anchors[label].get(accession, set())
                rf2_anchors = rf2_exact_anchors | rf2_gene_anchors
                if rf2_exact_anchors and rf2_gene_anchors:
                    match_mode = "EXACT_ACCESSION+GENE"
                elif rf2_exact_anchors:
                    match_mode = "EXACT_ACCESSION"
                elif rf2_gene_anchors:
                    match_mode = "GENE_PROPAGATED"
                else:
                    match_mode = "NONE"
                merged = hi_anchors | rf2_anchors
                row: dict[str, object] = {
                    "candidate_accession": accession,
                    "candidate_gene": joined(current_genes),
                    "stage1_route": route,
                    "ppi_screen_version": f"HI-union+RF2-{label}",
                    "ppi_pass": "yes" if merged else "no",
                    "ppi_evidence_category": evidence_category(hi_anchors, rf2_anchors, label),
                    "merged_anchor_count": len(merged),
                    "merged_anchor_genes": joined(merged),
                    "hi_union_anchor_count": len(hi_anchors),
                    "hi_union_anchor_genes": joined(hi_anchors),
                    "rf2_anchor_count": len(rf2_anchors),
                    "rf2_anchor_genes": joined(rf2_anchors),
                    "rf2_exact_accession_anchor_count": len(rf2_exact_anchors),
                    "rf2_exact_accession_anchor_genes": joined(rf2_exact_anchors),
                    "rf2_gene_propagated_anchor_count": len(rf2_gene_anchors),
                    "rf2_gene_propagated_anchor_genes": joined(rf2_gene_anchors),
                    "rf2_match_mode": match_mode,
                    "ppi_run_id": RUN_ID,
                }
                row.update(source)
                all_rows.append(row)
                if merged:
                    passed_rows.append(row)

                pair_map: dict[str, dict[str, object]] = {}
                for gene in current_genes:
                    for edge in hi_partner_to_rows.get(gene, []):
                        for anchor_accession in identifiers(edge.get("anchor_accessions")):
                            if anchor_accession not in anchor_accessions:
                                raise ValueError(
                                    f"HI edge anchor {anchor_accession} is not in current anchor table"
                                )
                            pair = pair_map.setdefault(
                                anchor_accession,
                                {
                                    "candidate_accession": accession,
                                    "candidate_gene": joined(current_genes),
                                    "anchor_accession": anchor_accession,
                                    "anchor_gene": anchor_accession_to_gene[anchor_accession],
                                    "ppi_screen_version": f"HI-union+RF2-{label}",
                                    "hi_union_support": 0,
                                    "rf2_support": 0,
                                    "rf2_match_modes": set(),
                                    "RFprob": set(),
                                    "AFprob": set(),
                                    "AFprob5": set(),
                                    "AFMprob": set(),
                                    "prediction_source": set(),
                                    "confDBs": set(),
                                    "allDBs": set(),
                                },
                            )
                            pair["hi_union_support"] = 1

                rf2_edge_map: dict[tuple[str, ...], tuple[dict[str, str], set[str]]] = {}
                for edge in rf2_accession_to_rows[label].get(accession, []):
                    edge_key = tuple(edge.get(field, "") for field in sorted(edge))
                    rf2_edge_map.setdefault(edge_key, (edge, set()))[1].add("EXACT_ACCESSION")
                for gene in current_genes:
                    for edge in rf2_gene_to_rows[label].get(gene, []):
                        edge_key = tuple(edge.get(field, "") for field in sorted(edge))
                        rf2_edge_map.setdefault(edge_key, (edge, set()))[1].add("GENE_PROPAGATED")
                for edge, modes in rf2_edge_map.values():
                    anchor_accession = edge["anchor_accession"].strip()
                    if anchor_accession not in anchor_accessions:
                        raise ValueError(
                            f"RF2 edge anchor {anchor_accession} is not in current anchor table"
                        )
                    pair = pair_map.setdefault(
                        anchor_accession,
                        {
                            "candidate_accession": accession,
                            "candidate_gene": joined(current_genes),
                            "anchor_accession": anchor_accession,
                            "anchor_gene": anchor_accession_to_gene[anchor_accession],
                            "ppi_screen_version": f"HI-union+RF2-{label}",
                            "hi_union_support": 0,
                            "rf2_support": 0,
                            "rf2_match_modes": set(),
                            "RFprob": set(),
                            "AFprob": set(),
                            "AFprob5": set(),
                            "AFMprob": set(),
                            "prediction_source": set(),
                            "confDBs": set(),
                            "allDBs": set(),
                        },
                    )
                    pair["rf2_support"] = 1
                    pair["rf2_match_modes"].update(modes)
                    for field in [
                        "RFprob", "AFprob", "AFprob5", "AFMprob",
                        "prediction_source", "confDBs", "allDBs",
                    ]:
                        value = edge.get(field, "").strip()
                        if value:
                            pair[field].add(value)

                for pair in pair_map.values():
                    hi_supported = bool(pair["hi_union_support"])
                    rf2_supported = bool(pair["rf2_support"])
                    if hi_supported and rf2_supported:
                        pair_category = f"HI-union+RF2-{label}"
                    elif hi_supported:
                        pair_category = "HI-union-only"
                    elif rf2_supported:
                        pair_category = f"RF2-{label}-only"
                    else:
                        raise AssertionError("Pair without PPI evidence")
                    pair_rows.append(
                        {
                            "candidate_accession": pair["candidate_accession"],
                            "candidate_gene": pair["candidate_gene"],
                            "anchor_accession": pair["anchor_accession"],
                            "anchor_gene": pair["anchor_gene"],
                            "ppi_screen_version": pair["ppi_screen_version"],
                            "pair_evidence_category": pair_category,
                            "hi_union_support": pair["hi_union_support"],
                            "rf2_support": pair["rf2_support"],
                            "rf2_match_mode": joined(pair["rf2_match_modes"]) or "NONE",
                            "RFprob": joined(pair["RFprob"]),
                            "AFprob": joined(pair["AFprob"]),
                            "AFprob5": joined(pair["AFprob5"]),
                            "AFMprob": joined(pair["AFMprob"]),
                            "prediction_source": joined(pair["prediction_source"]),
                            "confDBs": joined(pair["confDBs"]),
                            "allDBs": joined(pair["allDBs"]),
                            "ppi_run_id": RUN_ID,
                        }
                    )

            version_dir = out_root / f"0{1 if label == 'final80' else 2}_hi_union_rf2_{label}"
            all_path = version_dir / f"stage2_{route}_{expected}_hi_union_rf2_{label}_screen_all.tsv"
            passed_path = version_dir / f"stage2_{route}_hi_union_rf2_{label}_passed_{len(passed_rows)}.tsv"
            pair_rows.sort(
                key=lambda row: (
                    str(row["candidate_accession"]),
                    str(row["anchor_accession"]),
                )
            )
            pair_path = version_dir / f"stage2_{route}_hi_union_rf2_{label}_candidate_anchor_pairs.tsv"
            counts = Counter(str(row["ppi_evidence_category"]) for row in all_rows)
            summary_rows = [
                {"metric": "input_accessions", "value": len(source_rows)},
                {"metric": "input_atomic_genes", "value": count_genes(all_rows)},
                {"metric": "passed_accessions", "value": len(passed_rows)},
                {"metric": "passed_atomic_genes", "value": count_genes(passed_rows)},
                {"metric": "not_passed_accessions", "value": len(source_rows) - len(passed_rows)},
                {"metric": "candidate_anchor_pairs", "value": len(pair_rows)},
                {"metric": "distinct_supported_anchors", "value": len({row["anchor_accession"] for row in pair_rows})},
            ] + [
                {"metric": f"evidence_category::{name}", "value": value}
                for name, value in sorted(counts.items())
            ]
            summary_path = version_dir / f"stage2_{route}_hi_union_rf2_{label}_summary.tsv"
            safe_write(all_path, render_tsv(result_fields, all_rows))
            safe_write(passed_path, render_tsv(result_fields, passed_rows))
            safe_write(
                pair_path,
                render_tsv(
                    [
                        "candidate_accession", "candidate_gene", "anchor_accession", "anchor_gene",
                        "ppi_screen_version", "pair_evidence_category", "hi_union_support",
                        "rf2_support", "rf2_match_mode", "RFprob", "AFprob", "AFprob5",
                        "AFMprob", "prediction_source", "confDBs", "allDBs", "ppi_run_id",
                    ],
                    pair_rows,
                ),
            )
            safe_write(summary_path, render_tsv(["metric", "value"], summary_rows))
            route_outputs.extend([all_path, passed_path, pair_path, summary_path])
            threshold_data[label] = {
                "passed": {str(row["candidate_accession"]) for row in passed_rows},
                "passed_rows": passed_rows,
                "counts": dict(sorted(counts.items())),
                "pair_rows": pair_rows,
                "pair_set": {
                    (str(row["candidate_accession"]), str(row["anchor_accession"]))
                    for row in pair_rows
                },
                "all_path": all_path,
                "passed_path": passed_path,
                "pair_path": pair_path,
                "summary_path": summary_path,
            }

        pass80 = threshold_data["final80"]["passed"]
        pass90 = threshold_data["final90"]["passed"]
        if not pass90 <= pass80:
            raise ValueError(f"{route}: final90 pass set is not a subset of final80")
        if not threshold_data["final90"]["pair_set"] <= threshold_data["final80"]["pair_set"]:
            raise ValueError(f"{route}: final90 candidate-anchor pairs are not a subset of final80")
        comparison_rows = []
        for source in source_rows:
            accession = source["accession"].strip()
            in80 = accession in pass80
            in90 = accession in pass90
            comparison_rows.append(
                {
                    "candidate_accession": accession,
                    "candidate_gene": joined(candidate_genes[accession]),
                    "hi_rf2_final80_pass": "yes" if in80 else "no",
                    "hi_rf2_final90_pass": "yes" if in90 else "no",
                    "threshold_comparison_status": (
                        "PASS_BOTH" if in90 else ("PASS_FINAL80_ONLY" if in80 else "PASS_NEITHER")
                    ),
                }
            )
        comparison_path = out_root / "03_comparison" / f"stage2_{route}_final80_vs_final90.tsv"
        safe_write(
            comparison_path,
            render_tsv(
                [
                    "candidate_accession", "candidate_gene", "hi_rf2_final80_pass",
                    "hi_rf2_final90_pass", "threshold_comparison_status",
                ],
                comparison_rows,
            ),
        )
        route_outputs.append(comparison_path)

        comparison_counts = Counter(row["threshold_comparison_status"] for row in comparison_rows)
        audit = {
            "run_id": RUN_ID,
            "route": route,
            "route_definition": spec["definition"],
            "input_accessions": expected,
            "input_atomic_genes": count_genes(
                [{"candidate_gene": joined(candidate_genes[a])} for a in accessions]
            ),
            "current_anchor_accessions": len(anchor_accessions),
            "current_anchor_atomic_genes": len(anchor_genes),
            "final80": {
                "passed_accessions": len(pass80),
                "passed_atomic_genes": count_genes(threshold_data["final80"]["passed_rows"]),
                "candidate_anchor_pairs": len(threshold_data["final80"]["pair_rows"]),
                "distinct_supported_anchors": len(
                    {row["anchor_accession"] for row in threshold_data["final80"]["pair_rows"]}
                ),
                "evidence_counts": threshold_data["final80"]["counts"],
            },
            "final90": {
                "passed_accessions": len(pass90),
                "passed_atomic_genes": count_genes(threshold_data["final90"]["passed_rows"]),
                "candidate_anchor_pairs": len(threshold_data["final90"]["pair_rows"]),
                "distinct_supported_anchors": len(
                    {row["anchor_accession"] for row in threshold_data["final90"]["pair_rows"]}
                ),
                "evidence_counts": threshold_data["final90"]["counts"],
            },
            "comparison_counts": dict(sorted(comparison_counts.items())),
            "final90_pass_is_subset_of_final80_pass": True,
            "final90_pairs_are_subset_of_final80_pairs": True,
        }
        audit_path = out_root / "04_audit" / f"stage2_{route}_ppi_screen_audit.json"
        safe_write(
            audit_path,
            (json.dumps(audit, ensure_ascii=False, indent=2) + "\n").encode("utf-8"),
        )
        route_outputs.append(audit_path)

        readme = f"""# Stage2 {route}：HI-union + RF2-PPI双阈值筛选

更新日期：2026-09-04

Stage1输入定义：

    {spec["definition"]}

唯一Stage1输入为当前B2在排除正式B1后重新排名得到的266条canonical proteins。
本目录使用当前2,368个工作anchors，以及已经按这些anchors重新构建并冻结的
HI-union、RF2-final80和RF2-final90 anchor-partner边表。Stage2不读取历史268条
B2候选，也不将历史候选与当前结果并联。

## 结果

| 互作背景 | 输入proteins | 通过proteins | 通过atomic genes | 未通过proteins |
|---|---:|---:|---:|---:|
| HI-union ∪ RF2-final80 | {expected} | **{len(pass80)}** | {count_genes(threshold_data["final80"]["passed_rows"])} | {expected - len(pass80)} |
| HI-union ∪ RF2-final90 | {expected} | **{len(pass90)}** | {count_genes(threshold_data["final90"]["passed_rows"])} | {expected - len(pass90)} |

| 背景 | candidate–anchor protein pairs | 不同anchors |
|---|---:|---:|
| final80 | {len(threshold_data["final80"]["pair_rows"])} | {len({row["anchor_accession"] for row in threshold_data["final80"]["pair_rows"]})} |
| final90 | {len(threshold_data["final90"]["pair_rows"])} | {len({row["anchor_accession"] for row in threshold_data["final90"]["pair_rows"]})} |

final80和final90分别指RF2-PPI论文官方Data S4（80% precision）与Data S3
（90% precision）集合，而不是对某个概率列直接设置0.8/0.9阈值。final90
通过集合严格属于final80通过集合。

## 证据构成

| 背景 | HI与RF2共同支持 | 仅HI-union | 仅RF2 | 无数据集支持 |
|---|---:|---:|---:|---:|
| final80 | {threshold_data["final80"]["counts"].get("HI-union+RF2-final80", 0)} | {threshold_data["final80"]["counts"].get("HI-union-only", 0)} | {threshold_data["final80"]["counts"].get("RF2-final80-only", 0)} | {threshold_data["final80"]["counts"].get("NO_PPI_DATASET_SUPPORT", 0)} |
| final90 | {threshold_data["final90"]["counts"].get("HI-union+RF2-final90", 0)} | {threshold_data["final90"]["counts"].get("HI-union-only", 0)} | {threshold_data["final90"]["counts"].get("RF2-final90-only", 0)} | {threshold_data["final90"]["counts"].get("NO_PPI_DATASET_SUPPORT", 0)} |

完整screen_all表覆盖全部Stage1输入并保留ppi_pass=yes/no；passed表只包含至少
与一个当前anchor存在HI-union或相应RF2官方集合支持的候选。candidate-anchor
pairs表将候选和相应的调控anchor展开为一行一个蛋白pair，并保留HI/RF2来源及
RF2原始分数字段。

## 目录

    01_hi_union_rf2_final80/      final80完整表、passed表、pair表、summary
    02_hi_union_rf2_final90/      final90完整表、passed表、pair表、summary
    03_comparison/                逐候选80/90比较
    04_audit/                     输入、结果、集合关系和SHA256
    scripts/                      可复现构建脚本

共享anchor-partner背景位于：

    /root/autodl-tmp/Agent_analysis_v2/04_stage2/A/ppi_screening/

本结果表示存在实验PPI数据或RF2-PPI预测支持，不等同于新的直接结合实验验证。
""".encode("utf-8")
        readme_path = out_root / "README.md"
        safe_write(readme_path, readme)
        route_outputs.append(readme_path)

        script_target = out_root / "scripts" / "build_stage2_B2_current.py"
        if Path(__file__).resolve() != script_target.resolve():
            safe_write(script_target, Path(__file__).read_bytes())
        route_outputs.append(script_target)

        manifest_rows = []
        for input_file in [ANCHORS, input_path, HI_EDGES, *RF2_EDGES.values()]:
            manifest_rows.append(
                {
                    "role": "INPUT",
                    "path": str(input_file),
                    "size_bytes": input_file.stat().st_size,
                    "sha256": sha256(input_file),
                }
            )
        for output_file in route_outputs:
            manifest_rows.append(
                {
                    "role": "OUTPUT",
                    "path": str(output_file),
                    "size_bytes": output_file.stat().st_size,
                    "sha256": sha256(output_file),
                }
            )
        manifest_path = out_root / "04_audit" / f"stage2_{route}_ppi_screen_manifest.tsv"
        safe_write(
            manifest_path,
            render_tsv(["role", "path", "size_bytes", "sha256"], manifest_rows),
        )
        route_results[route] = audit

    print(json.dumps(route_results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
