#!/usr/bin/env python3
"""Build v2 B2 annotation-graph and OmniPath features.

This is a clean extraction of the validated v1 B2 graph logic.  It contains no
direct nuclear-localization model, no candidate-specific keyword bonus and no
PPI evidence.  It deliberately stops before the B1-dependent percentile rank.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import platform
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

import networkx as nx
import numpy as np
import pandas as pd


RUN_ID = "stage1_B2_current14625_network_features_20260904"
EXPECTED_CANDIDATES = 14625
EXPECTED_ANNOTATIONS = 20416
EXPECTED_ANCHORS = 2368


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--common-input", type=Path, required=True)
    parser.add_argument("--all-annotations", type=Path, required=True)
    parser.add_argument("--anchor-table", type=Path, required=True)
    parser.add_argument("--resource-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_dataframe(frame: pd.DataFrame, path: Path) -> None:
    payload = frame.to_csv(sep="\t", index=False, lineterminator="\n").encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(f"Refusing to overwrite conflicting output: {path}")
        return
    partial = Path(str(path) + ".partial")
    if partial.exists():
        raise FileExistsError(f"Stale partial file requires manual review: {partial}")
    with partial.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def rank_percentile(score: pd.Series) -> Tuple[pd.Series, pd.Series]:
    rank = score.rank(method="min", ascending=False).astype(int)
    percentile = 1.0 - (rank - 1) / max(len(score) - 1, 1)
    return rank, percentile


def positive_percentile_score(score: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(score, errors="coerce").fillna(0.0).astype(float)
    if float(numeric.max()) <= 0:
        return pd.Series(np.zeros(len(numeric)), index=score.index)
    out = numeric.rank(method="max", pct=True).astype(float)
    out[numeric <= 0] = 0.0
    return out


def split_semicolon_terms(value: str) -> List[str]:
    out: List[str] = []
    for chunk in str(value or "").split(";"):
        chunk = re.sub(r"\s*\[GO:\d+\]", "", chunk).strip()
        if chunk:
            out.append(chunk)
    return out


def relevant_terms_for_graph(row: pd.Series, reactome_terms: List[str]) -> List[str]:
    terms: List[str] = []
    for column in [
        "Gene Ontology (biological process)",
        "Gene Ontology (cellular component)",
        "Gene Ontology (molecular function)",
        "Keywords",
    ]:
        value = str(row.get(column, ""))
        chunks = (
            [item.strip() for item in value.split(";") if item.strip()]
            if column == "Keywords"
            else split_semicolon_terms(value)
        )
        terms.extend(chunks)
    for item in reactome_terms:
        pathway_id, name = item.split("|", 1) if "|" in item else (item, item)
        terms.append(f"Reactome:{pathway_id}:{name}")
    seen: set[str] = set()
    normalized: List[str] = []
    for term in terms:
        key = re.sub(r"\s+", " ", term.strip().lower())
        if key and key not in seen:
            seen.add(key)
            normalized.append(key[:180])
    return normalized


def term_edge_weight(graph: nx.Graph, term_node: str, total_protein_nodes: int) -> float:
    protein_degree = sum(1 for node in graph.neighbors(term_node) if node.startswith("protein:"))
    if protein_degree <= 0 or total_protein_nodes <= 0:
        return 2.0
    max_idf = math.log1p(total_protein_nodes)
    idf = math.log((1.0 + total_protein_nodes) / (1.0 + protein_degree))
    normalized_idf = 0.0 if max_idf <= 0 else max(0.0, min(1.0, idf / max_idf))
    return 1.0 + (1.0 - normalized_idf)


def parse_reactome(resource_dir: Path) -> Dict[str, List[str]]:
    mapping_path = resource_dir / "UniProt2Reactome_PE_Pathway.txt"
    names_path = resource_dir / "ReactomePathways.txt"
    pathway_names: Dict[str, str] = {}
    with names_path.open(errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                pathway_names[parts[0]] = parts[1]
    acc_to_pathways: Dict[str, List[str]] = defaultdict(list)
    with mapping_path.open(errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8 or parts[7] != "Homo sapiens":
                continue
            accession, pathway_id, fallback_name = parts[0], parts[3], parts[5]
            name = pathway_names.get(pathway_id, fallback_name)
            acc_to_pathways[accession].append(f"{pathway_id}|{name}")
    return acc_to_pathways


def parse_omnipath(resource_dir: Path) -> pd.DataFrame:
    path = resource_dir / "omnipath_interactions_human.tsv"
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def build_omnipath_directed_graph(omnipath: pd.DataFrame, valid_accs: Set[str]) -> nx.DiGraph:
    graph = nx.DiGraph()
    if not {"source", "target"}.issubset(omnipath.columns):
        raise ValueError("OmniPath file lacks source/target columns")
    for row in omnipath.itertuples(index=False):
        source = getattr(row, "source")
        target = getattr(row, "target")
        if not source or not target:
            continue
        if source.startswith("COMPLEX:") or target.startswith("COMPLEX:"):
            continue
        if source not in valid_accs or target not in valid_accs:
            continue
        if str(getattr(row, "consensus_direction", "True")) == "False":
            continue
        graph.add_edge(
            source,
            target,
            stimulation=str(getattr(row, "consensus_stimulation", "")),
            inhibition=str(getattr(row, "consensus_inhibition", "")),
            sources=str(getattr(row, "sources", "")),
            curation_effort=str(getattr(row, "curation_effort", "")),
        )
    return graph


def compute_omnipath_seed_features(
    graph: nx.DiGraph,
    candidate_accs: Sequence[str],
    seed_accs: Set[str],
    max_distance: int = 4,
) -> Dict[str, Dict[str, object]]:
    seed_nodes = [accession for accession in sorted(seed_accs) if accession in graph]
    reverse_dist = (
        nx.multi_source_dijkstra_path_length(
            graph.reverse(copy=False), seed_nodes, cutoff=max_distance
        )
        if seed_nodes
        else {}
    )
    seed_set = set(seed_nodes)
    features: Dict[str, Dict[str, object]] = {}
    for accession in candidate_accs:
        out_edges = list(graph.out_edges(accession, data=True)) if accession in graph else []
        if accession in graph:
            lengths = nx.single_source_shortest_path_length(graph, accession, cutoff=3)
            reachable = sum(1 for node in lengths if node in seed_set and node != accession)
        else:
            reachable = 0
        first_steps: List[str] = []
        for _source, target, data in out_edges[:8]:
            edge_label = (
                "stim"
                if data.get("stimulation") == "True"
                else "inh"
                if data.get("inhibition") == "True"
                else "dir"
            )
            first_steps.append(f"{target}:{edge_label}")
        features[accession] = {
            "omnipath_directed_distance_to_seed": reverse_dist.get(accession, np.nan),
            "omnipath_reachable_seed_count_3hop": reachable,
            "omnipath_out_edge_count": len(out_edges),
            "omnipath_first_step_targets": "; ".join(first_steps),
        }
    return features


def estimate_weighted_shared_seed(
    graph: nx.Graph,
    node: str,
    seed_set: Set[str],
    total_protein_nodes: int,
) -> float:
    if node not in graph:
        return 0.0
    score = 0.0
    for term in sorted(set(graph.neighbors(node))):
        if not term.startswith("term:"):
            continue
        neighbors = sorted(
            other
            for other in graph.neighbors(term)
            if other.startswith("protein:") and other != node
        )
        degree = len(neighbors)
        if degree == 0:
            continue
        seed_count = sum(1 for other in neighbors if other in seed_set)
        if seed_count == 0:
            continue
        idf = math.log((1.0 + total_protein_nodes) / (1.0 + degree))
        score += (seed_count / degree) * max(0.0, idf)
    return score


def describe_signal_paths(
    reactome_terms: List[str],
    shortest: object,
    shared_seed_idf_raw: float,
    omnipath_distance: object,
    omnipath_reach_count: int,
    omnipath_first_steps: object,
) -> str:
    parts: List[str] = []
    if reactome_terms:
        parts.append("Reactome: " + "; ".join(reactome_terms[:3]))
    if shortest == shortest:
        parts.append(f"shortest graph distance to TF/epifactor seed={shortest}")
    if shared_seed_idf_raw:
        parts.append(f"IDF-weighted shared TF/epifactor seed evidence={shared_seed_idf_raw:.3f}")
    if omnipath_distance == omnipath_distance:
        parts.append(f"OmniPath directed distance to TF/epifactor seed={omnipath_distance}")
    if omnipath_reach_count:
        parts.append(f"OmniPath reachable seed count <=3 hops={omnipath_reach_count}")
    if omnipath_first_steps:
        parts.append(f"OmniPath first steps: {omnipath_first_steps}")
    return " | ".join(parts) if parts else "No graph support in frozen resources"


def main() -> None:
    args = parse_args()
    required_resources = [
        "UniProt2Reactome_PE_Pathway.txt",
        "ReactomePathways.txt",
        "ReactomePathwaysRelation.txt",
        "omnipath_interactions_human.tsv",
    ]
    for name in required_resources:
        path = args.resource_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    candidates_raw = read_tsv(args.common_input)
    annotations_raw = read_tsv(args.all_annotations)
    anchors = read_tsv(args.anchor_table)
    if (
        len(candidates_raw) != EXPECTED_CANDIDATES
        or candidates_raw["Entry"].nunique() != EXPECTED_CANDIDATES
    ):
        raise ValueError("B2 common input must contain 14,625 unique accessions")
    if (
        len(annotations_raw) != EXPECTED_ANNOTATIONS
        or annotations_raw["Entry"].nunique() != EXPECTED_ANNOTATIONS
    ):
        raise ValueError("All-annotation table must contain 20,416 unique accessions")
    if (
        len(anchors) != EXPECTED_ANCHORS
        or anchors["anchor_accession"].nunique() != EXPECTED_ANCHORS
    ):
        raise ValueError("Working anchor table must contain 2,368 unique accessions")

    candidates = candidates_raw.rename(
        columns={
            "Entry": "accession",
            "Protein names": "protein_name",
            "Gene Names (primary)": "gene",
            "Length": "length",
        }
    ).copy()
    annotations = annotations_raw.rename(
        columns={
            "Entry": "accession",
            "Protein names": "protein_name",
            "Gene Names (primary)": "gene",
        }
    ).set_index("accession", drop=False)
    candidate_accs = set(candidates["accession"])
    seed_accs = set(anchors["anchor_accession"])
    valid_accs = set(annotations["accession"])
    if candidate_accs & seed_accs:
        raise ValueError("B2 common input overlaps working anchor accessions")
    if not (candidate_accs | seed_accs) <= valid_accs:
        raise ValueError("Candidate or anchor accession is missing from the annotation universe")

    acc_to_pathways = parse_reactome(args.resource_dir)
    omnipath = parse_omnipath(args.resource_dir)
    omnipath_graph = build_omnipath_directed_graph(omnipath, valid_accs)
    omnipath_features = compute_omnipath_seed_features(
        omnipath_graph,
        list(candidates["accession"]),
        seed_accs,
        max_distance=4,
    )

    graph = nx.Graph()
    graph_accessions = candidate_accs | seed_accs
    for accession, row in annotations.iterrows():
        if accession not in graph_accessions:
            continue
        graph.add_node(f"protein:{accession}", kind="protein")
        for term in relevant_terms_for_graph(row, acc_to_pathways.get(accession, [])):
            term_node = f"term:{term}"
            graph.add_node(term_node, kind="term")
            graph.add_edge(f"protein:{accession}", term_node)

    seed_nodes = [
        f"protein:{accession}"
        for accession in sorted(seed_accs)
        if f"protein:{accession}" in graph
    ]
    if len(seed_nodes) != EXPECTED_ANCHORS:
        raise ValueError(
            f"Only {len(seed_nodes)} of {EXPECTED_ANCHORS} anchors entered the annotation graph"
        )
    seed_node_set = set(seed_nodes)
    total_protein_nodes = sum(1 for node in graph.nodes if node.startswith("protein:"))
    for term_node in [node for node in graph.nodes if node.startswith("term:")]:
        weight = term_edge_weight(graph, term_node, total_protein_nodes)
        for neighbor in graph.neighbors(term_node):
            graph.edges[term_node, neighbor]["weight"] = weight
    distances = nx.multi_source_dijkstra_path_length(
        graph,
        seed_nodes,
        cutoff=6.0,
        weight="weight",
    )

    rows: List[dict[str, object]] = []
    for _, row in candidates.iterrows():
        accession = str(row["accession"])
        node = f"protein:{accession}"
        reactome_names: List[str] = []
        seen_names: set[str] = set()
        for mapped in acc_to_pathways.get(accession, []):
            name = mapped.split("|", 1)[1] if "|" in mapped else mapped
            key = name.strip().lower()
            if key and key not in seen_names:
                seen_names.add(key)
                reactome_names.append(name)
        shortest = distances.get(node, np.nan)
        graph_proximity = 0.0 if pd.isna(shortest) else max(0.0, 1.0 - float(shortest) / 6.0)
        shared = estimate_weighted_shared_seed(graph, node, seed_node_set, total_protein_nodes)
        omni = omnipath_features[accession]
        omni_distance = omni["omnipath_directed_distance_to_seed"]
        omni_proximity = (
            0.0
            if pd.isna(omni_distance)
            else max(0.0, 1.0 - float(omni_distance) / 5.0)
        )
        omni_reach = int(omni["omnipath_reachable_seed_count_3hop"] or 0)
        rows.append(
            {
                "accession": accession,
                "gene": row.get("gene", ""),
                "protein_name": row.get("protein_name", ""),
                "length": row.get("length", ""),
                # Historical field name retained for compatibility.  It counts
                # all mapped human Reactome pathways, not only signaling paths.
                "reactome_signal_pathway_count": len(reactome_names),
                "reactome_signal_pathways": "; ".join(reactome_names[:8]),
                "shortest_path_to_tf_epifactor_seed": shortest,
                "shared_seed_idf_raw": shared,
                "graph_proximity_score": graph_proximity,
                "omnipath_directed_distance_to_seed": omni_distance,
                "omnipath_reachable_seed_count_3hop": omni_reach,
                "omnipath_proximity_score": omni_proximity,
                "omnipath_out_edge_count": int(omni["omnipath_out_edge_count"] or 0),
                "omnipath_first_step_targets": omni["omnipath_first_step_targets"],
                "top_supporting_paths": describe_signal_paths(
                    reactome_names,
                    shortest,
                    shared,
                    omni_distance,
                    omni_reach,
                    omni["omnipath_first_step_targets"],
                ),
            }
        )
    features = pd.DataFrame(rows)
    features["shared_seed_idf_score_legacy_graph"] = positive_percentile_score(
        features["shared_seed_idf_raw"]
    )
    features["omnipath_reach_score_legacy_graph"] = positive_percentile_score(
        features["omnipath_reachable_seed_count_3hop"]
    )
    weights = {
        "graph": 0.15 / 0.49,
        "shared": 0.06 / 0.49,
        "omni_distance": 0.25 / 0.49,
        "omni_reach": 0.03 / 0.49,
    }
    features["signal_to_nucleus_score_legacy_graph"] = (
        weights["graph"] * features["graph_proximity_score"]
        + weights["shared"] * features["shared_seed_idf_score_legacy_graph"]
        + weights["omni_distance"] * features["omnipath_proximity_score"]
        + weights["omni_reach"] * features["omnipath_reach_score_legacy_graph"]
    ).clip(0, 1)
    (
        features["signal_to_nucleus_rank_legacy_graph"],
        features["signal_to_nucleus_percentile_legacy_graph"],
    ) = rank_percentile(features["signal_to_nucleus_score_legacy_graph"])
    features = features.sort_values(
        ["signal_to_nucleus_rank_legacy_graph", "accession"]
    ).reset_index(drop=True)

    output_path = args.output_dir / "b2_network_features_14625.tsv"
    stats_path = args.output_dir / "b2_graph_stats.tsv"
    input_manifest_path = args.output_dir / "b2_network_input_manifest.tsv"
    atomic_write_dataframe(features, output_path)

    stats = pd.DataFrame(
        [
            {
                "run_id": RUN_ID,
                "candidate_nodes": len(candidate_accs),
                "anchor_seed_nodes": len(seed_nodes),
                "annotation_graph_protein_nodes": total_protein_nodes,
                "annotation_graph_term_nodes": sum(1 for node in graph if node.startswith("term:")),
                "annotation_graph_edges": graph.number_of_edges(),
                "reactome_mapped_human_proteins": len(acc_to_pathways),
                "omnipath_directed_nodes": omnipath_graph.number_of_nodes(),
                "omnipath_directed_edges": omnipath_graph.number_of_edges(),
                "manual_stage_nodes": 0,
                "candidate_specific_keyword_bonus": 0,
                "formal_percentile_ranking_computed": 0,
                "formal_B1_exclusion_applied": 0,
            }
        ]
    )
    atomic_write_dataframe(stats, stats_path)

    artifacts = [
        (
            "INPUT",
            args.common_input,
            EXPECTED_CANDIDATES,
            "current Route B common non-anchor universe",
        ),
        (
            "INPUT",
            args.all_annotations,
            EXPECTED_ANNOTATIONS,
            "frozen reviewed canonical annotation universe",
        ),
        ("INPUT", args.anchor_table, EXPECTED_ANCHORS, "current anchor seeds"),
    ]
    artifacts.extend(
        ("RESOURCE", args.resource_dir / name, "", "frozen B2 external resource")
        for name in required_resources
    )
    artifacts.extend(
        [
            ("CODE", Path(__file__).resolve(), "", "standalone v2 B2 feature builder"),
            (
                "OUTPUT",
                output_path,
                EXPECTED_CANDIDATES,
                "B1-independent raw B2 features",
            ),
            ("OUTPUT", stats_path, 1, "graph construction audit"),
        ]
    )
    if input_manifest_path.exists():
        existing_manifest = read_tsv(input_manifest_path)
        existing_times = {
            value
            for value in existing_manifest.get("generated_at_utc", pd.Series(dtype=str))
            if value
        }
        if len(existing_times) != 1:
            raise ValueError("Existing input manifest has inconsistent timestamps")
        generated_at = next(iter(existing_times))
    else:
        generated_at = datetime.now(timezone.utc).isoformat()
    manifest_rows = []
    for role, path, records, description in artifacts:
        manifest_rows.append(
            {
                "artifact_role": role,
                "path": str(path.resolve()),
                "records": records,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "description": description,
                "run_id": RUN_ID,
                "generated_at_utc": generated_at,
            }
        )
    manifest = pd.DataFrame(manifest_rows)
    manifest["python_version"] = platform.python_version()
    manifest["pandas_version"] = pd.__version__
    manifest["numpy_version"] = np.__version__
    manifest["networkx_version"] = nx.__version__
    atomic_write_dataframe(manifest, input_manifest_path)

    print(f"B2 network features complete: {len(features)}")
    print(f"annotation graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    print(f"OmniPath graph: {omnipath_graph.number_of_nodes()} nodes, {omnipath_graph.number_of_edges()} edges")


if __name__ == "__main__":
    main()
