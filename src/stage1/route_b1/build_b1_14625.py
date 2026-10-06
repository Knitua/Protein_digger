#!/usr/bin/env python3
"""Build the clean Stage1-B1 release on the frozen 14,625-protein universe.

The formal score table preserves the already validated corrected NLSExplorer
scores for the historical 14,556 proteins and appends corrected scores for the
69 newly eligible proteins.  A complete rerun of both the corrected and legacy
NLSExplorer implementations is retained for reproducibility and sensitivity
analysis, but the legacy result is never unioned into the formal B1 set.
"""

from __future__ import annotations

import argparse
import hashlib
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_CURRENT = 14_625
EXPECTED_PREVIOUS = 14_556
EXPECTED_NEW = 69
THRESHOLD = 0.5
RUN_ID = "stage1_B1_complete_current_universe_20260904"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_unique(path: Path, key: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    if key not in frame.columns:
        raise ValueError(f"Missing key {key!r} in {path}")
    if frame[key].isna().any() or frame[key].astype(str).duplicated().any():
        raise ValueError(f"Missing or duplicate {key!r} in {path}")
    frame[key] = frame[key].astype(str)
    return frame


def atomic_dataframe(frame: pd.DataFrame, path: Path, **kwargs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    frame.to_csv(partial, **kwargs)
    os.replace(partial, path)


def atomic_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    partial.write_text(text)
    os.replace(partial, path)


def write_fasta(frame: pd.DataFrame, path: Path, method: str) -> None:
    chunks: list[str] = []
    for row in frame.itertuples(index=False):
        gene = "" if pd.isna(row.gene) else str(row.gene)
        chunks.append(f">{row.accession} GN={gene} B1={method}\n")
        seq = str(row.sequence)
        chunks.extend(seq[start : start + 80] + "\n" for start in range(0, len(seq), 80))
    atomic_text("".join(chunks), path)


def normalize_nls(frame: pd.DataFrame) -> pd.DataFrame:
    required = [
        "accession",
        "n_windows",
        "nls_prob_max",
        "nls_prob_mean",
        "max_window_start_1based",
        "max_window_end_1based",
    ]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing NLS columns: {missing}")
    return frame[required].rename(
        columns={
            "n_windows": "nlsexplorer_n_windows",
            "nls_prob_max": "nlsexplorer_nls_prob_max",
            "nls_prob_mean": "nlsexplorer_nls_prob_mean",
            "max_window_start_1based": "nlsexplorer_max_window_start_1based",
            "max_window_end_1based": "nlsexplorer_max_window_end_1based",
        }
    )


def normalize_deeploc_new(frame: pd.DataFrame) -> pd.DataFrame:
    required = [
        "accession",
        "deeploc_nucleus_prob",
        "deeploc_top_label",
        "deeploc_top_prob",
        "was_clipped_gt_4000",
    ]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing DeepLoc columns: {missing}")
    return frame[required].rename(
        columns={"was_clipped_gt_4000": "deeploc_was_clipped_gt_4000"}
    )


def score_table(
    current: pd.DataFrame,
    deep: pd.DataFrame,
    nls: pd.DataFrame,
    method: str,
) -> pd.DataFrame:
    scores = (
        current[["Entry", "Length"]]
        .rename(columns={"Entry": "accession", "Length": "length"})
        .merge(deep, on="accession", how="left", validate="one_to_one")
        .merge(nls, on="accession", how="left", validate="one_to_one")
    )
    required_values = [
        "length",
        "deeploc_nucleus_prob",
        "deeploc_top_label",
        "deeploc_top_prob",
        "deeploc_was_clipped_gt_4000",
        "nlsexplorer_n_windows",
        "nlsexplorer_nls_prob_max",
        "nlsexplorer_nls_prob_mean",
        "nlsexplorer_max_window_start_1based",
        "nlsexplorer_max_window_end_1based",
    ]
    missing = [column for column in required_values if scores[column].isna().any()]
    if missing:
        raise AssertionError(f"Incomplete model values for {method}: {missing}")
    scores.insert(2, "nlsexplorer_method", method)
    scores["deeploc_gt_0p5"] = scores["deeploc_nucleus_prob"].gt(THRESHOLD).astype("int8")
    scores["nlsexplorer_gt_0p5"] = scores["nlsexplorer_nls_prob_max"].gt(THRESHOLD).astype("int8")
    scores["both_gt_0p5"] = (
        scores["deeploc_gt_0p5"].astype(bool)
        & scores["nlsexplorer_gt_0p5"].astype(bool)
    ).astype("int8")
    scores["joint_min_prob"] = np.minimum(
        scores["deeploc_nucleus_prob"], scores["nlsexplorer_nls_prob_max"]
    )
    scores["joint_mean_prob"] = (
        scores["deeploc_nucleus_prob"] + scores["nlsexplorer_nls_prob_max"]
    ) / 2.0
    scores["run_id"] = RUN_ID
    return scores.sort_values("accession", kind="mergesort").reset_index(drop=True)


def candidate_table(scores: pd.DataFrame, current: pd.DataFrame) -> pd.DataFrame:
    positive = scores.loc[scores["both_gt_0p5"].eq(1)].copy()
    positive = positive.sort_values(
        ["joint_min_prob", "joint_mean_prob", "accession"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    positive.insert(0, "rank", np.arange(1, len(positive) + 1))

    metadata = current.rename(
        columns={
            "Entry": "accession",
            "Entry Name": "entry_name",
            "Protein names": "protein_name",
            "Gene Names (primary)": "gene",
            "Length": "metadata_length",
            "Sequence": "sequence",
        }
    )
    score_columns = [
        "rank",
        "accession",
        "nlsexplorer_method",
        "deeploc_nucleus_prob",
        "nlsexplorer_nls_prob_max",
        "joint_min_prob",
        "joint_mean_prob",
        "deeploc_top_label",
        "deeploc_top_prob",
        "deeploc_was_clipped_gt_4000",
        "nlsexplorer_n_windows",
        "nlsexplorer_max_window_start_1based",
        "nlsexplorer_max_window_end_1based",
        "run_id",
    ]
    result = positive[score_columns].merge(
        metadata, on="accession", how="left", validate="one_to_one"
    )
    if result["metadata_length"].isna().any():
        raise AssertionError("Incomplete candidate metadata join")
    result = result.rename(columns={"metadata_length": "length"})
    front = [
        "rank",
        "accession",
        "entry_name",
        "gene",
        "protein_name",
        "length",
        "nlsexplorer_method",
    ]
    return result[front + [column for column in result.columns if column not in front]]


def manifest_row(role: str, path: Path, records: int, note: str, timestamp: str) -> dict:
    return {
        "artifact_role": role,
        "path": str(path.resolve()),
        "records": records,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "description": note,
        "run_id": RUN_ID,
        "generated_at_utc": timestamp,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, required=True)
    parser.add_argument("--current-tsv", type=Path, required=True)
    parser.add_argument("--previous-tsv", type=Path, required=True)
    parser.add_argument("--new69-tsv", type=Path, required=True)
    parser.add_argument("--previous-corrected-combined", type=Path, required=True)
    parser.add_argument("--new69-deeploc", type=Path, required=True)
    parser.add_argument("--corrected-full-rerun", type=Path, required=True)
    parser.add_argument("--legacy-full-rerun", type=Path, required=True)
    parser.add_argument("--historical-legacy-candidates", type=Path)
    parser.add_argument("--nls-weights", type=Path, required=True)
    parser.add_argument("--esm-weights", type=Path, required=True)
    args = parser.parse_args()

    current = read_unique(args.current_tsv, "Entry")
    previous = read_unique(args.previous_tsv, "Entry")
    new69 = read_unique(args.new69_tsv, "Entry")
    old_combined = read_unique(args.previous_corrected_combined, "accession")
    new_deep_raw = read_unique(args.new69_deeploc, "accession")
    corrected_full_raw = read_unique(args.corrected_full_rerun, "accession")
    legacy_full_raw = read_unique(args.legacy_full_rerun, "accession")

    if (len(current), len(previous), len(new69), len(old_combined)) != (
        EXPECTED_CURRENT,
        EXPECTED_PREVIOUS,
        EXPECTED_NEW,
        EXPECTED_PREVIOUS,
    ):
        raise AssertionError("Unexpected B1 input counts")
    current_set = set(current["Entry"])
    previous_set = set(previous["Entry"])
    new_set = set(new69["Entry"])
    if previous_set & new_set or previous_set | new_set != current_set:
        raise AssertionError("14,625 is not the disjoint union of 14,556 and 69")
    for label, frame in (
        ("corrected full rerun", corrected_full_raw),
        ("legacy full rerun", legacy_full_raw),
    ):
        if len(frame) != EXPECTED_CURRENT or set(frame["accession"]) != current_set:
            raise AssertionError(f"{label} does not exactly cover current B1 universe")
    if len(new_deep_raw) != EXPECTED_NEW or set(new_deep_raw["accession"]) != new_set:
        raise AssertionError("New DeepLoc result does not exactly cover 69 new proteins")

    deep_columns = [
        "accession",
        "deeploc_nucleus_prob",
        "deeploc_top_label",
        "deeploc_top_prob",
        "deeploc_was_clipped_gt_4000",
    ]
    deep = pd.concat(
        [old_combined[deep_columns], normalize_deeploc_new(new_deep_raw)],
        ignore_index=True,
    )
    if len(deep) != EXPECTED_CURRENT or deep["accession"].duplicated().any():
        raise AssertionError("Merged DeepLoc table is not a unique 14,625-row table")

    old_nls_columns = [
        "accession",
        "nlsexplorer_n_windows",
        "nlsexplorer_nls_prob_max",
        "nlsexplorer_nls_prob_mean",
        "nlsexplorer_max_window_start_1based",
        "nlsexplorer_max_window_end_1based",
    ]
    corrected_full = normalize_nls(corrected_full_raw)
    corrected_formal = pd.concat(
        [
            old_combined[old_nls_columns],
            corrected_full.loc[corrected_full["accession"].isin(new_set)],
        ],
        ignore_index=True,
    )
    if len(corrected_formal) != EXPECTED_CURRENT or corrected_formal["accession"].duplicated().any():
        raise AssertionError("Formal corrected NLS table is not a unique 14,625-row table")
    legacy_full = normalize_nls(legacy_full_raw)

    formal_scores = score_table(current, deep, corrected_formal, "CORRECTED_EXACT_LENGTH")
    legacy_scores = score_table(current, deep, legacy_full, "LEGACY_BATCH_PADDED")
    corrected_rerun_scores = score_table(
        current, deep, corrected_full, "CORRECTED_EXACT_LENGTH_FULL_RERUN_AUDIT"
    )
    formal_candidates = candidate_table(formal_scores, current)
    legacy_candidates = candidate_table(legacy_scores, current)

    scores_dir = args.release_root / "01_model_scores"
    candidates_dir = args.release_root / "02_candidates"
    audit_dir = args.release_root / "03_audit"
    formal_score_path = scores_dir / "b1_corrected_formal_scores_14625.tsv"
    legacy_score_path = scores_dir / "b1_legacy_batch_padded_scores_14625.tsv"
    corrected_rerun_path = audit_dir / "b1_corrected_full_rerun_scores_14625.tsv"
    formal_candidate_path = candidates_dir / "b1_corrected_positive_both_gt0p5.tsv"
    formal_accession_path = candidates_dir / "b1_corrected_positive_both_gt0p5.accessions.txt"
    formal_fasta_path = candidates_dir / "b1_corrected_positive_both_gt0p5.fasta"
    legacy_dir = candidates_dir / "legacy_batch_padded_comparator"
    legacy_candidate_path = legacy_dir / "b1_legacy_positive_both_gt0p5.tsv"
    legacy_accession_path = legacy_dir / "b1_legacy_positive_both_gt0p5.accessions.txt"
    legacy_fasta_path = legacy_dir / "b1_legacy_positive_both_gt0p5.fasta"

    for frame, path in (
        (formal_scores, formal_score_path),
        (legacy_scores, legacy_score_path),
        (corrected_rerun_scores, corrected_rerun_path),
        (formal_candidates, formal_candidate_path),
        (legacy_candidates, legacy_candidate_path),
    ):
        atomic_dataframe(frame, path, sep="\t", index=False, float_format="%.8f")
    atomic_text("\n".join(formal_candidates["accession"]) + "\n", formal_accession_path)
    atomic_text("\n".join(legacy_candidates["accession"]) + "\n", legacy_accession_path)
    write_fasta(formal_candidates[["accession", "gene", "sequence"]], formal_fasta_path, "CORRECTED")
    write_fasta(legacy_candidates[["accession", "gene", "sequence"]], legacy_fasta_path, "LEGACY")

    formal_set = set(formal_candidates["accession"])
    legacy_set = set(legacy_candidates["accession"])
    rerun_set = set(
        corrected_rerun_scores.loc[corrected_rerun_scores["both_gt_0p5"].eq(1), "accession"]
    )
    membership = current[["Entry", "Gene Names (primary)"]].rename(
        columns={"Entry": "accession", "Gene Names (primary)": "gene"}
    )
    membership["corrected_formal_positive"] = membership["accession"].isin(formal_set).astype("int8")
    membership["legacy_batch_padded_positive"] = membership["accession"].isin(legacy_set).astype("int8")
    membership["corrected_full_rerun_positive"] = membership["accession"].isin(rerun_set).astype("int8")
    membership["comparison"] = np.select(
        [
            membership["corrected_formal_positive"].eq(1) & membership["legacy_batch_padded_positive"].eq(1),
            membership["corrected_formal_positive"].eq(1),
            membership["legacy_batch_padded_positive"].eq(1),
        ],
        ["BOTH", "CORRECTED_ONLY", "LEGACY_ONLY"],
        default="NEITHER",
    )
    membership_path = audit_dir / "b1_candidate_membership_comparison.tsv"
    atomic_dataframe(membership, membership_path, sep="\t", index=False)

    old_compare = old_combined[
        ["accession", "nlsexplorer_nls_prob_max", "nlsexplorer_gt_0p5"]
    ].merge(
        corrected_full[
            ["accession", "nlsexplorer_nls_prob_max"]
        ].rename(columns={"nlsexplorer_nls_prob_max": "rerun_nls_prob_max"}),
        on="accession",
        validate="one_to_one",
    )
    old_compare["abs_delta"] = (
        old_compare["nlsexplorer_nls_prob_max"] - old_compare["rerun_nls_prob_max"]
    ).abs()
    old_compare["rerun_gt_0p5"] = old_compare["rerun_nls_prob_max"].gt(THRESHOLD).astype("int8")
    old_compare["threshold_flip"] = (
        old_compare["nlsexplorer_gt_0p5"].astype(int) != old_compare["rerun_gt_0p5"]
    ).astype("int8")
    consistency_path = audit_dir / "b1_corrected_recompute_consistency_14556.tsv"
    atomic_dataframe(old_compare, consistency_path, sep="\t", index=False, float_format="%.8f")

    historical_set: set[str] = set()
    if args.historical_legacy_candidates:
        historical = read_unique(args.historical_legacy_candidates, "accession")
        historical_set = set(historical["accession"])

    counts = pd.DataFrame(
        [
            ("current_B1_universe", EXPECTED_CURRENT, "fixed Stage0 route-B non-anchor input"),
            ("previous_validated_corrected_component", EXPECTED_PREVIOUS, "values preserved"),
            ("new_same_method_component", EXPECTED_NEW, "new corrected NLS and DeepLoc inference"),
            ("deeploc_gt_0p5", int(formal_scores["deeploc_gt_0p5"].sum()), "same in both tracks"),
            ("corrected_nls_gt_0p5", int(formal_scores["nlsexplorer_gt_0p5"].sum()), "formal exact-length method"),
            ("corrected_B1_positive", len(formal_set), "formal B1 result"),
            ("legacy_nls_gt_0p5", int(legacy_scores["nlsexplorer_gt_0p5"].sum()), "batch-padded comparator"),
            ("legacy_B1_positive", len(legacy_set), "comparator only"),
            ("corrected_and_legacy_positive", len(formal_set & legacy_set), "intersection"),
            ("corrected_only_positive", len(formal_set - legacy_set), "not unioned"),
            ("legacy_only_positive", len(legacy_set - formal_set), "not admitted to formal B1"),
            ("corrected_full_rerun_B1_positive", len(rerun_set), "validation rerun only"),
            ("corrected_formal_vs_rerun_set_difference", len(formal_set ^ rerun_set), "FP16 reproducibility audit"),
            ("corrected_recompute_nls_threshold_flips_14556", int(old_compare["threshold_flip"].sum()), "old validated versus full rerun"),
            ("historical_legacy_B1_positive", len(historical_set), "archived historical result"),
            ("historical_legacy_overlap_current_formal", len(historical_set & formal_set), "audit only"),
        ],
        columns=["metric", "value", "note"],
    )
    counts_path = audit_dir / "b1_counts.tsv"
    atomic_dataframe(counts, counts_path, sep="\t", index=False)

    comparison = pd.DataFrame(
        [
            ("CORRECTED_EXACT_LENGTH", len(formal_set), "FORMAL", "DeepLoc > 0.5 AND corrected NLS > 0.5"),
            ("LEGACY_BATCH_PADDED", len(legacy_set), "COMPARATOR_ONLY", "DeepLoc > 0.5 AND legacy NLS > 0.5"),
            ("INTERSECTION", len(formal_set & legacy_set), "AUDIT", "positive in both versions"),
            ("CORRECTED_ONLY", len(formal_set - legacy_set), "FORMAL_ONLY", "positive only under corrected NLS"),
            ("LEGACY_ONLY", len(legacy_set - formal_set), "EXCLUDED_FROM_FORMAL", "positive only under legacy NLS"),
        ],
        columns=["result", "proteins", "status", "definition"],
    )
    comparison_path = audit_dir / "b1_version_comparison.tsv"
    atomic_dataframe(comparison, comparison_path, sep="\t", index=False)

    timestamp = datetime.now(timezone.utc).isoformat()
    artifacts = [
        ("INPUT", args.current_tsv, len(current), "frozen current B1 universe"),
        ("INPUT", args.previous_tsv, len(previous), "validated corrected component universe"),
        ("INPUT", args.new69_tsv, len(new69), "new same-method component universe"),
        ("MODEL_SCORE_INPUT", args.previous_corrected_combined, len(old_combined), "validated corrected and DeepLoc values"),
        ("MODEL_SCORE_INPUT", args.new69_deeploc, len(new_deep_raw), "new 69 DeepLoc values"),
        ("MODEL_SCORE_INPUT", args.corrected_full_rerun, len(corrected_full_raw), "corrected full rerun"),
        ("MODEL_SCORE_INPUT", args.legacy_full_rerun, len(legacy_full_raw), "legacy full rerun"),
        ("MODEL_WEIGHT", args.nls_weights, 1, "NLSExplorer classifier weights"),
        ("MODEL_WEIGHT", args.esm_weights, 1, "ESM-1b encoder weights"),
        ("OUTPUT", formal_score_path, len(formal_scores), "formal corrected score table"),
        ("OUTPUT", formal_candidate_path, len(formal_candidates), "formal corrected B1 candidates"),
        ("OUTPUT", formal_accession_path, len(formal_candidates), "formal accession interface"),
        ("OUTPUT", formal_fasta_path, len(formal_candidates), "formal candidate FASTA"),
        ("OUTPUT", legacy_score_path, len(legacy_scores), "legacy comparator score table"),
        ("OUTPUT", legacy_candidate_path, len(legacy_candidates), "legacy comparator candidates"),
        ("AUDIT", corrected_rerun_path, len(corrected_rerun_scores), "complete corrected rerun audit"),
        ("AUDIT", membership_path, len(membership), "per-protein version membership"),
        ("AUDIT", consistency_path, len(old_compare), "corrected FP16 recompute audit"),
        ("AUDIT", counts_path, len(counts), "count audit"),
        ("AUDIT", comparison_path, len(comparison), "version comparison"),
    ]
    manifest = pd.DataFrame(
        [manifest_row(*artifact, timestamp) for artifact in artifacts],
    )
    manifest_path = audit_dir / "b1_manifest.tsv"
    atomic_dataframe(manifest, manifest_path, sep="\t", index=False)

    print(f"eligible\t{len(formal_scores)}")
    print(f"deeploc_gt_0p5\t{int(formal_scores['deeploc_gt_0p5'].sum())}")
    print(f"corrected_nls_gt_0p5\t{int(formal_scores['nlsexplorer_gt_0p5'].sum())}")
    print(f"corrected_B1_positive\t{len(formal_set)}")
    print(f"legacy_nls_gt_0p5\t{int(legacy_scores['nlsexplorer_gt_0p5'].sum())}")
    print(f"legacy_B1_positive\t{len(legacy_set)}")
    print(f"intersection\t{len(formal_set & legacy_set)}")
    print(f"corrected_only\t{len(formal_set - legacy_set)}")
    print(f"legacy_only\t{len(legacy_set - formal_set)}")
    print(f"corrected_recompute_threshold_flips\t{int(old_compare['threshold_flip'].sum())}")
    print(f"max_corrected_recompute_abs_delta\t{old_compare['abs_delta'].max():.8f}")


if __name__ == "__main__":
    main()
