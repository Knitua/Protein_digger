#!/usr/bin/env python3
"""Read-only validation for the current Stage1-B2 release."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED = {
    "common_rows": 14625,
    "anchor_rows": 2368,
    "b1_rows": 1350,
    "eligible_rows": 13275,
    "provisional_top": 293,
    "formal_top": 266,
}

LOCKED_HASHES = {
    "common_input": "a02286dc3fe7a9e40746731893149042c33a6d4f959e4ef8fbe6071dd9edcde0",
    "anchor_table": "4229cfc1ec6bbde1a8bf90e1f52f102c1c8cf069b128978f2207853e32f8eff5",
    "b1_accessions": "4f3382d59419eea0373fa1be16fe56ef83bc68e4f0dd3d70e62c1261bf13598c",
    "b1_candidates": "b203e6646a11029ad43ff1c820f57c2a11933d4d3f398622d23a5f7043153939",
}

RESOURCE_HASHES = {
    "ReactomePathways.txt": "eac3ac074b3849a6d928260f4f71920f346584b0114afa7b60fcc8e40eefd9ec",
    "ReactomePathwaysRelation.txt": "72e2da17ef5806d287a0d932cd6d161fea62ef7a567b6e95601b67463c07170e",
    "UniProt2Reactome_PE_Pathway.txt": "3ec08f9b152412040864a7c263bd7028b77dd6393b288a422ee83ee79c75389c",
    "omnipath_interactions_human.tsv": "e353485d9a5e56d34ed45483427c0d26cee7e4d9559827cf43121e6873d1c787",
}


def read(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("b2_rank_current", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def metric_map(frame: pd.DataFrame) -> dict[str, str]:
    if list(frame.columns) != ["metric", "value"]:
        raise ValueError("Summary must contain exactly metric/value columns")
    if frame["metric"].duplicated().any():
        raise ValueError("Summary contains duplicate metrics")
    return dict(zip(frame["metric"], frame["value"]))


def compare_ranking(observed: pd.DataFrame, expected: pd.DataFrame) -> None:
    if list(observed["accession"]) != list(expected["accession"]):
        raise ValueError("Recomputed accession order differs")
    columns = [
        "p_graph_distance",
        "p_shared_seed_idf",
        "p_omnipath_distance",
        "p_omnipath_reach",
        "network_score",
        "network_rank",
        "network_percentile",
    ]
    for column in columns:
        left = pd.to_numeric(observed[column], errors="coerce")
        right = pd.to_numeric(expected[column], errors="coerce")
        if not np.allclose(left, right, rtol=0, atol=1e-15, equal_nan=True):
            raise ValueError(f"Recomputed values differ: {column}")


def validate_manifest(path: Path) -> None:
    manifest = read(path)
    for _, row in manifest.iterrows():
        artifact = Path(row["path"])
        if not artifact.is_file():
            raise FileNotFoundError(artifact)
        if int(row["bytes"]) != artifact.stat().st_size:
            raise ValueError(f"Manifest size differs: {artifact}")
        if row["sha256"] != sha256(artifact):
            raise ValueError(f"Manifest SHA256 differs: {artifact}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, required=True)
    parser.add_argument("--common-input", type=Path, required=True)
    parser.add_argument("--anchor-table", type=Path, required=True)
    parser.add_argument("--b1-accessions", type=Path, required=True)
    parser.add_argument("--b1-candidates", type=Path, required=True)
    parser.add_argument("--resource-dir", type=Path, required=True)
    args = parser.parse_args()

    locked_paths = {
        "common_input": args.common_input,
        "anchor_table": args.anchor_table,
        "b1_accessions": args.b1_accessions,
        "b1_candidates": args.b1_candidates,
    }
    for key, path in locked_paths.items():
        if sha256(path) != LOCKED_HASHES[key]:
            raise ValueError(f"Locked input hash differs: {key}")
    for name, expected_hash in RESOURCE_HASHES.items():
        if sha256(args.resource_dir / name) != expected_hash:
            raise ValueError(f"Frozen resource hash differs: {name}")

    common = read(args.common_input)
    anchors = read(args.anchor_table)
    b1_candidates = read(args.b1_candidates)
    b1_values = [
        value.strip()
        for value in args.b1_accessions.read_text(encoding="utf-8").splitlines()
        if value.strip()
    ]
    common_set = set(common["Entry"])
    anchor_set = set(anchors["anchor_accession"])
    b1_set = set(b1_values)
    if len(common) != EXPECTED["common_rows"] or len(common_set) != EXPECTED["common_rows"]:
        raise ValueError("Common input is not 14,625 unique accessions")
    if len(anchors) != EXPECTED["anchor_rows"] or len(anchor_set) != EXPECTED["anchor_rows"]:
        raise ValueError("Anchor table is not 2,368 unique accessions")
    if common_set & anchor_set:
        raise ValueError("Stage0-B non-anchor universe overlaps the anchor set")
    if len(b1_values) != EXPECTED["b1_rows"] or len(b1_set) != EXPECTED["b1_rows"]:
        raise ValueError("B1 exclusion interface is not 1,350 unique accessions")
    if set(b1_candidates["accession"]) != b1_set:
        raise ValueError("B1 candidate table and accession interface differ")
    if not b1_set <= common_set:
        raise ValueError("B1 exclusions are outside Stage0-B")

    features_path = args.release_root / "02_network_features/b2_network_features_14625.tsv"
    provisional_dir = args.release_root / "03_provisional_ranking_without_B1"
    formal_dir = args.release_root / "04_formal_after_B1"
    rank_script = args.release_root / "scripts/rank_b2_network.py"
    features = read(features_path)
    if len(features) != EXPECTED["common_rows"] or set(features["accession"]) != common_set:
        raise ValueError("Network features do not exactly cover Stage0-B")
    if features["accession"].duplicated().any():
        raise ValueError("Network feature table contains duplicate accessions")

    provisional_full = read(provisional_dir / "b2_provisional_full_ranking_14625.tsv")
    provisional_top = read(provisional_dir / "b2_provisional_top2pct_293.tsv")
    provisional_summary = metric_map(read(provisional_dir / "b2_provisional_summary.tsv"))
    formal_full = read(formal_dir / "b2_formal_full_ranking_13275.tsv")
    formal_top = read(formal_dir / "b2_formal_top2pct_266.tsv")
    formal_summary = metric_map(read(formal_dir / "b2_formal_summary.tsv"))

    if len(provisional_full) != EXPECTED["common_rows"]:
        raise ValueError("Provisional ranking row count differs")
    if len(provisional_top) != EXPECTED["provisional_top"]:
        raise ValueError("Provisional Top 2% row count differs")
    if list(provisional_full.head(EXPECTED["provisional_top"])["accession"]) != list(provisional_top["accession"]):
        raise ValueError("Provisional top set is not the head of the full ranking")
    if len(formal_full) != EXPECTED["eligible_rows"]:
        raise ValueError("Formal ranking is not 13,275 rows")
    if set(formal_full["accession"]) != common_set - b1_set:
        raise ValueError("Formal ranking is not the exact non-B1 universe")
    if len(formal_top) != EXPECTED["formal_top"]:
        raise ValueError("Formal Top 2% row count differs")
    if list(formal_full.head(EXPECTED["formal_top"])["accession"]) != list(formal_top["accession"]):
        raise ValueError("Formal top set is not the head of the full ranking")
    if set(formal_top["accession"]) & (b1_set | anchor_set):
        raise ValueError("Formal B2 candidates overlap B1 or anchors")

    module = load_module(rank_script)
    compare_ranking(provisional_full, module.compute_ranking(features, set()))
    compare_ranking(formal_full, module.compute_ranking(features, b1_set))

    expected_summary = {
        "raw_network_feature_rows": "14625",
        "B1_exclusion_input_records": "1350",
        "B1_exclusion_input_unique_accessions": "1350",
        "B1_excluded_unique_accessions": "1350",
        "eligible_rows": "13275",
        "selected_rows": "266",
        "formal_candidate_set_frozen": "1",
    }
    for key, value in expected_summary.items():
        if formal_summary.get(key) != value:
            raise ValueError(f"Formal summary differs: {key}")
    if provisional_summary.get("selected_rows") != "293":
        raise ValueError("Provisional summary differs")

    stats = read(args.release_root / "02_network_features/b2_graph_stats.tsv")
    if len(stats) != 1:
        raise ValueError("Graph statistics must have one row")
    row = stats.iloc[0]
    if int(row["candidate_nodes"]) != 14625 or int(row["anchor_seed_nodes"]) != 2368:
        raise ValueError("Graph node counts differ")
    if int(row["manual_stage_nodes"]) != 0 or int(row["candidate_specific_keyword_bonus"]) != 0:
        raise ValueError("Candidate-specific logic detected")

    validate_manifest(args.release_root / "02_network_features/b2_network_input_manifest.tsv")
    validate_manifest(formal_dir / "b2_formal_manifest.tsv")
    if EXPECTED["formal_top"] != math.ceil(EXPECTED["eligible_rows"] * 0.02):
        raise ValueError("Top 2% arithmetic differs")
    print("PASS: current Stage1-B2 release validated")
    print("14625 common -> exclude 1350 B1 -> 13275 eligible -> 266 Top 2%")


if __name__ == "__main__":
    main()
