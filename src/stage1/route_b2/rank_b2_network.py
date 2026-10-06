#!/usr/bin/env python3
"""Apply the validated four-percentile B2 ranking rule.

Formal mode consumes the frozen B1 accession interface.  The preferred input
is a headerless file containing one accession per line; tabular B1 releases
remain supported when an accession column is specified or can be detected.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def set_sha256(values: Iterable[str]) -> str:
    payload = "".join(f"{value}\n" for value in sorted(set(values)))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def atomic_write(frame: pd.DataFrame, path: Path) -> None:
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


def read_b1_exclusions(path: Path, accession_column: str | None) -> list[str]:
    """Read and validate a frozen B1 accession interface."""
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".txt":
        values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
        values = [value for value in values if value]
        if any("\t" in value for value in values):
            raise ValueError("Headerless B1 accession input must contain one value per line")
    else:
        direct = read_tsv(path)
        if accession_column is None:
            detected = [
                column
                for column in ["accession", "candidate_accession", "Entry"]
                if column in direct.columns
            ]
            if len(detected) != 1:
                raise ValueError(
                    "Could not uniquely detect the B1 accession column; "
                    "use --b1-accession-column"
                )
            accession_column = detected[0]
        if accession_column not in direct.columns:
            raise ValueError(f"Missing B1 accession column: {accession_column}")
        values = [str(value).strip() for value in direct[accession_column]]
        values = [value for value in values if value]
    if not values:
        raise ValueError("B1 exclusion input is empty")
    if len(values) != len(set(values)):
        raise ValueError("B1 exclusion input contains duplicate accessions")
    return values


def beneficial_percentile(
    values: pd.Series,
    *,
    higher_is_better: bool,
    zero_is_no_evidence: bool,
) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    n = len(numeric)
    if n <= 1:
        return pd.Series(np.ones(n), index=values.index)
    valid = numeric.notna()
    work = numeric.fillna(-np.inf if higher_is_better else np.inf)
    rank = work.rank(method="min", ascending=not higher_is_better)
    percentile = 1.0 - (rank - 1.0) / (n - 1.0)
    percentile[~valid] = 0.0
    if zero_is_no_evidence:
        percentile[numeric.fillna(0) <= 0] = 0.0
    return percentile.clip(0, 1)


def compute_ranking(raw: pd.DataFrame, excluded_accessions: set[str]) -> pd.DataFrame:
    eligible = raw.loc[~raw["accession"].isin(excluded_accessions)].copy()
    if eligible["accession"].duplicated().any():
        raise ValueError("Duplicate accessions in B2 feature input")
    eligible["p_graph_distance"] = beneficial_percentile(
        eligible["shortest_path_to_tf_epifactor_seed"],
        higher_is_better=False,
        zero_is_no_evidence=False,
    )
    eligible["p_shared_seed_idf"] = beneficial_percentile(
        eligible["shared_seed_idf_raw"],
        higher_is_better=True,
        zero_is_no_evidence=True,
    )
    eligible["p_omnipath_distance"] = beneficial_percentile(
        eligible["omnipath_directed_distance_to_seed"],
        higher_is_better=False,
        zero_is_no_evidence=False,
    )
    eligible["p_omnipath_reach"] = beneficial_percentile(
        eligible["omnipath_reachable_seed_count_3hop"],
        higher_is_better=True,
        zero_is_no_evidence=True,
    )
    eligible["network_score"] = eligible[
        [
            "p_graph_distance",
            "p_shared_seed_idf",
            "p_omnipath_distance",
            "p_omnipath_reach",
        ]
    ].mean(axis=1)
    eligible = eligible.sort_values(
        [
            "network_score",
            "p_graph_distance",
            "p_shared_seed_idf",
            "p_omnipath_distance",
            "p_omnipath_reach",
            "accession",
        ],
        ascending=[False, False, False, False, False, True],
    ).reset_index(drop=True)
    eligible["network_rank"] = np.arange(1, len(eligible) + 1)
    eligible["network_percentile"] = 1.0 - (
        (eligible["network_rank"] - 1) / max(len(eligible) - 1, 1)
    )
    return eligible


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--network-features", type=Path, required=True)
    parser.add_argument("--mode", choices=["provisional", "formal"], required=True)
    parser.add_argument("--b1-exclusions", type=Path)
    parser.add_argument("--b1-accession-column")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-fraction", type=float, default=0.02)
    parser.add_argument("--run-id", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw = read_tsv(args.network_features)
    if args.mode == "provisional" and args.b1_exclusions is not None:
        raise ValueError("Provisional mode must not accept a B1 exclusion table")
    if args.mode == "formal" and args.b1_exclusions is None:
        raise ValueError("Formal mode requires the frozen B1 exclusion table")
    excluded_values: list[str] = []
    excluded: set[str] = set()
    if args.b1_exclusions is not None:
        excluded_values = read_b1_exclusions(
            args.b1_exclusions, args.b1_accession_column
        )
        excluded = set(excluded_values)
        missing = excluded - set(raw["accession"])
        if missing:
            examples = ", ".join(sorted(missing)[:5])
            raise ValueError(
                f"B1 exclusions are outside the B2 feature universe: {examples}"
            )

    ranked = compute_ranking(raw, excluded)
    top_n = math.ceil(len(ranked) * args.top_fraction)
    selected = ranked.head(top_n).copy()
    status = (
        "PROVISIONAL_B1_NOT_EXCLUDED"
        if args.mode == "provisional"
        else "FORMAL_B1_EXCLUDED"
    )
    ranked.insert(0, "ranking_status", status)
    ranked.insert(1, "run_id", args.run_id)
    selected.insert(0, "ranking_status", status)
    selected.insert(1, "run_id", args.run_id)
    selected["pass_rule"] = (
        f"equal-weight four-percentile B2 rank; top ceil({args.top_fraction:.4f} * N)"
    )

    prefix = "b2_provisional" if args.mode == "provisional" else "b2_formal"
    full_path = args.output_dir / f"{prefix}_full_ranking_{len(ranked)}.tsv"
    top_path = args.output_dir / f"{prefix}_top2pct_{top_n}.tsv"
    summary_path = args.output_dir / f"{prefix}_summary.tsv"
    atomic_write(ranked, full_path)
    atomic_write(selected, top_path)
    summary = pd.DataFrame(
        [
            {"metric": "ranking_status", "value": status},
            {"metric": "raw_network_feature_rows", "value": len(raw)},
            {"metric": "B1_exclusion_input_path", "value": str(args.b1_exclusions or "")},
            {"metric": "B1_exclusion_input_sha256", "value": sha256_file(args.b1_exclusions) if args.b1_exclusions else ""},
            {"metric": "B1_exclusion_input_records", "value": len(excluded_values)},
            {"metric": "B1_exclusion_input_unique_accessions", "value": len(excluded)},
            {"metric": "B1_excluded_unique_accessions", "value": len(set(raw["accession"]) & excluded)},
            {"metric": "B1_excluded_accession_set_sha256", "value": set_sha256(excluded)},
            {"metric": "eligible_rows", "value": len(ranked)},
            {"metric": "top_fraction", "value": args.top_fraction},
            {"metric": "selected_rows", "value": len(selected)},
            {"metric": "selected_nonblank_unique_genes", "value": selected.loc[selected["gene"] != "", "gene"].nunique()},
            {"metric": "selected_accession_set_sha256", "value": set_sha256(selected["accession"])},
            {"metric": "formal_candidate_set_frozen", "value": 1 if args.mode == "formal" else 0},
        ]
    )
    atomic_write(summary, summary_path)
    if args.mode == "formal":
        manifest_rows = []
        artifacts = [
            ("INPUT", args.network_features, len(raw), "frozen B2 network features"),
            ("INPUT", args.b1_exclusions, len(excluded_values), "frozen current B1 accession exclusions"),
            ("CODE", Path(__file__).resolve(), "", "formal B2 ranking implementation"),
            ("OUTPUT", full_path, len(ranked), "formal full B2 ranking after B1 exclusion"),
            ("OUTPUT", top_path, len(selected), "frozen formal B2 Top 2% candidates"),
            ("AUDIT", summary_path, len(summary), "formal B2 count and set audit"),
        ]
        for role, path, records, description in artifacts:
            if path is None:
                raise ValueError("Formal manifest input path is missing")
            path = path.resolve()
            manifest_rows.append(
                {
                    "artifact_role": role,
                    "path": str(path),
                    "records": records,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "description": description,
                    "run_id": args.run_id,
                }
            )
        atomic_write(
            pd.DataFrame(manifest_rows),
            args.output_dir / "b2_formal_manifest.tsv",
        )
    print(f"{status}: eligible={len(ranked)}, top2pct={len(selected)}")


if __name__ == "__main__":
    main()
