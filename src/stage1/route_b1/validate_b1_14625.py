#!/usr/bin/env python3
"""Validate the clean Stage1-B1 14,625-protein release."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession is not None:
                    records[accession] = "".join(chunks)
                accession = line[1:].split()[0]
                if accession in records:
                    raise AssertionError(f"Duplicate FASTA accession: {accession}")
                chunks = []
            else:
                chunks.append(line)
    if accession is not None:
        records[accession] = "".join(chunks)
    return records


def validate_track(root: Path, prefix: str, score_path: Path, candidate_path: Path,
                   accession_path: Path, fasta_path: Path, expected_method: str) -> tuple[int, set[str]]:
    scores = read_tsv(score_path)
    candidates = read_tsv(candidate_path)
    if len(scores) != 14_625 or scores["accession"].duplicated().any():
        raise AssertionError(f"{prefix}: score table is not a unique 14,625-row table")
    if set(scores["nlsexplorer_method"]) != {expected_method}:
        raise AssertionError(f"{prefix}: unexpected NLS method label")
    calculated = (
        scores["deeploc_nucleus_prob"].gt(0.5)
        & scores["nlsexplorer_nls_prob_max"].gt(0.5)
    ).astype(int)
    if not calculated.equals(scores["both_gt_0p5"].astype(int)):
        raise AssertionError(f"{prefix}: threshold flags do not match strict >0.5 rule")
    expected = set(scores.loc[calculated.eq(1), "accession"].astype(str))
    observed = set(candidates["accession"].astype(str))
    if expected != observed or len(candidates) != len(expected):
        raise AssertionError(f"{prefix}: candidate set does not match score thresholds")
    if candidates["rank"].tolist() != list(range(1, len(candidates) + 1)):
        raise AssertionError(f"{prefix}: ranks are not consecutive")
    expected_order = candidates.sort_values(
        ["joint_min_prob", "joint_mean_prob", "accession"],
        ascending=[False, False, True],
        kind="mergesort",
    )["accession"].tolist()
    if expected_order != candidates["accession"].tolist():
        raise AssertionError(f"{prefix}: candidate ordering is not stable")
    accessions = [line.strip() for line in accession_path.read_text().splitlines() if line.strip()]
    if accessions != candidates["accession"].astype(str).tolist():
        raise AssertionError(f"{prefix}: accession interface differs from candidate table")
    fasta = read_fasta(fasta_path)
    if list(fasta) != accessions:
        raise AssertionError(f"{prefix}: FASTA order differs from candidate table")
    sequences = dict(zip(candidates["accession"].astype(str), candidates["sequence"].astype(str)))
    if fasta != sequences:
        raise AssertionError(f"{prefix}: FASTA sequences differ from candidate metadata")
    return len(candidates), observed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, required=True)
    parser.add_argument("--current-tsv", type=Path, required=True)
    parser.add_argument("--previous-tsv", type=Path, required=True)
    parser.add_argument("--new69-tsv", type=Path, required=True)
    parser.add_argument("--previous-corrected-combined", type=Path, required=True)
    parser.add_argument("--corrected-full-rerun", type=Path, required=True)
    args = parser.parse_args()

    root = args.release_root
    current = read_tsv(args.current_tsv)
    previous = read_tsv(args.previous_tsv)
    new69 = read_tsv(args.new69_tsv)
    if len(current) != 14_625 or len(previous) != 14_556 or len(new69) != 69:
        raise AssertionError("Unexpected Stage0/B1 universe sizes")
    current_set = set(current["Entry"].astype(str))
    previous_set = set(previous["Entry"].astype(str))
    new_set = set(new69["Entry"].astype(str))
    if previous_set & new_set or previous_set | new_set != current_set:
        raise AssertionError("Current B1 universe is not 14,556 disjoint-union 69")

    formal_score = root / "01_model_scores/b1_corrected_formal_scores_14625.tsv"
    formal_candidate = root / "02_candidates/b1_corrected_positive_both_gt0p5.tsv"
    formal_accession = root / "02_candidates/b1_corrected_positive_both_gt0p5.accessions.txt"
    formal_fasta = root / "02_candidates/b1_corrected_positive_both_gt0p5.fasta"
    legacy_score = root / "01_model_scores/b1_legacy_batch_padded_scores_14625.tsv"
    legacy_candidate = root / "02_candidates/legacy_batch_padded_comparator/b1_legacy_positive_both_gt0p5.tsv"
    legacy_accession = root / "02_candidates/legacy_batch_padded_comparator/b1_legacy_positive_both_gt0p5.accessions.txt"
    legacy_fasta = root / "02_candidates/legacy_batch_padded_comparator/b1_legacy_positive_both_gt0p5.fasta"

    n_formal, formal_set = validate_track(
        root, "formal", formal_score, formal_candidate, formal_accession, formal_fasta,
        "CORRECTED_EXACT_LENGTH",
    )
    n_legacy, legacy_set = validate_track(
        root, "legacy", legacy_score, legacy_candidate, legacy_accession, legacy_fasta,
        "LEGACY_BATCH_PADDED",
    )

    formal = read_tsv(formal_score).set_index("accession")
    old = read_tsv(args.previous_corrected_combined).set_index("accession")
    corrected_rerun = read_tsv(args.corrected_full_rerun).set_index("accession")
    nls_columns = [
        "nlsexplorer_n_windows",
        "nlsexplorer_nls_prob_max",
        "nlsexplorer_nls_prob_mean",
        "nlsexplorer_max_window_start_1based",
        "nlsexplorer_max_window_end_1based",
    ]
    pd.testing.assert_frame_equal(
        formal.loc[old.index, nls_columns].sort_index(),
        old[nls_columns].sort_index(),
        check_dtype=False,
        check_exact=True,
    )
    rerun_new = corrected_rerun.loc[list(new_set)].rename(
        columns={
            "n_windows": "nlsexplorer_n_windows",
            "nls_prob_max": "nlsexplorer_nls_prob_max",
            "nls_prob_mean": "nlsexplorer_nls_prob_mean",
            "max_window_start_1based": "nlsexplorer_max_window_start_1based",
            "max_window_end_1based": "nlsexplorer_max_window_end_1based",
        }
    )
    pd.testing.assert_frame_equal(
        formal.loc[list(new_set), nls_columns].sort_index(),
        rerun_new[nls_columns].sort_index(),
        check_dtype=False,
        check_exact=True,
    )

    comparison = read_tsv(root / "03_audit/b1_version_comparison.tsv")
    values = dict(zip(comparison["result"], comparison["proteins"]))
    expected_values = {
        "CORRECTED_EXACT_LENGTH": n_formal,
        "LEGACY_BATCH_PADDED": n_legacy,
        "INTERSECTION": len(formal_set & legacy_set),
        "CORRECTED_ONLY": len(formal_set - legacy_set),
        "LEGACY_ONLY": len(legacy_set - formal_set),
    }
    if values != expected_values:
        raise AssertionError("Version-comparison counts are inconsistent")
    if "EXCLUDED_FROM_FORMAL" not in set(comparison["status"]):
        raise AssertionError("Legacy-only results are not explicitly excluded from formal B1")

    manifest = read_tsv(root / "03_audit/b1_manifest.tsv")
    for row in manifest.itertuples(index=False):
        path = Path(row.path)
        if not path.is_absolute():
            path = root / path
        if not path.exists() or sha256(path) != row.sha256:
            raise AssertionError(f"Manifest SHA256 mismatch: {path}")
    partials = list(root.rglob("*.partial"))
    if partials:
        raise AssertionError(f"Unresolved partial outputs: {partials}")

    print("PASS\tStage1-B1 complete-current-universe validation")
    print(f"eligible\t{len(current)}")
    print(f"formal_corrected_positive\t{n_formal}")
    print(f"legacy_comparator_positive\t{n_legacy}")
    print(f"intersection\t{len(formal_set & legacy_set)}")
    print(f"corrected_only\t{len(formal_set - legacy_set)}")
    print(f"legacy_only\t{len(legacy_set - formal_set)}")


if __name__ == "__main__":
    main()
