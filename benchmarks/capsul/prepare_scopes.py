#!/usr/bin/env python3
"""Freeze CAPSUL data scopes without exporting test targets to training code."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import pandas as pd


EXPECTED_SHA256 = "50459272402c7e939108428298d444ee7c977d9a5630148665a98397ed0309f2"
EXPECTED_COUNTS = {"train": 14126, "valid": 3027, "test": 3028}
LABELS = [
    "Nucleus", "Nuclear Membrane", "Nucleoli", "Nucleoplasm", "Cytoplasm",
    "Cytosol", "Cytoskeleton", "Centrosome", "Mitochondria",
    "Endoplasmic Reticulum", "Golgi Apparatus",
    "Plasma Membrane/Cell Membrane", "Endosome", "Lipid droplet",
    "Lysosome/Vacuole", "Peroxisome", "Vesicle", "Primary Cilium",
    "Secreted Proteins", "Sperm",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_hash(rows) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(json.dumps(list(row), ensure_ascii=False, separators=(",", ":")).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def middleclip(sequence: str) -> str:
    sequence = str(sequence)
    return sequence if len(sequence) <= 1022 else sequence[:511] + sequence[-511:]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--union-csv", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--allow-smoke", action="store_true")
    args = parser.parse_args()
    source = Path(args.union_csv).resolve()
    out = Path(args.out_dir).resolve()
    data = out / "data"
    audits = out / "audits"
    data.mkdir(parents=True, exist_ok=True)
    audits.mkdir(parents=True, exist_ok=True)

    observed_sha = sha256_file(source)
    if observed_sha != EXPECTED_SHA256:
        raise AssertionError(f"union.csv SHA256 mismatch: {observed_sha}")
    frame = pd.read_csv(source)
    required = {"ID", "Dataset", "Fasta Sequence", *LABELS}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"union.csv missing columns: {missing}")
    frame["ID"] = frame["ID"].astype(str)
    frame["Dataset"] = frame["Dataset"].astype(str).str.lower()
    if len(frame) != 20181 or frame["ID"].nunique() != 20181:
        raise AssertionError("Expected 20,181 unique protein IDs")
    counts = {str(k): int(v) for k, v in frame["Dataset"].value_counts().items()}
    if counts != EXPECTED_COUNTS:
        raise AssertionError(f"Unexpected official splits: {counts}")
    raw_values = set(frame[LABELS].stack().astype(int).unique().tolist())
    if not raw_values.issubset({0, 1, 2}) or 2 not in raw_values:
        raise AssertionError(f"Unexpected label values: {raw_values}")

    binary = frame.copy()
    binary[LABELS] = (binary[LABELS].to_numpy(dtype=float) > 0).astype("int8")
    scope_columns = ["ID", *LABELS, "Dataset", "Fasta Sequence"]
    train_path = data / "capsul_train_locked.csv"
    valid_path = data / "capsul_valid_locked.csv"
    sequence_path = data / "capsul_sequences_only_all_splits.csv"
    middle_path = data / "capsul_sequences_only_middleclip1022.csv"
    binary.loc[binary["Dataset"] == "train", scope_columns].to_csv(train_path, index=False)
    binary.loc[binary["Dataset"] == "valid", scope_columns].to_csv(valid_path, index=False)
    sequences = frame[["ID", "Dataset", "Fasta Sequence"]].copy()
    sequences.to_csv(sequence_path, index=False)
    clipped = sequences.copy()
    clipped["Fasta Sequence"] = clipped["Fasta Sequence"].map(middleclip)
    clipped.to_csv(middle_path, index=False)

    smoke_ids = []
    for split in ("train", "valid", "test"):
        smoke_ids.extend(frame.loc[frame["Dataset"] == split, "ID"].head(100).tolist())
    smoke_path = data / "smoke_100_100_100_ids.csv"
    pd.DataFrame({"ID": smoke_ids}).to_csv(smoke_path, index=False)

    sequence_groups = frame.groupby("Fasta Sequence", sort=False)
    duplicate_groups = []
    duplicate_rows = 0
    cross_split_groups = 0
    cross_split_rows = 0
    for sequence, group in sequence_groups:
        if len(group) <= 1:
            continue
        duplicate_rows += len(group)
        splits = sorted(group["Dataset"].unique().tolist())
        is_cross = len(splits) > 1
        if is_cross:
            cross_split_groups += 1
            cross_split_rows += len(group)
        duplicate_groups.append({
            "sequence_sha256": hashlib.sha256(str(sequence).encode()).hexdigest(),
            "ids": group["ID"].tolist(),
            "splits": group["Dataset"].tolist(),
            "cross_split": is_cross,
        })
    duplicates_path = audits / "exact_sequence_duplicates.json"
    duplicates_path.write_text(json.dumps({
        "duplicate_group_count": len(duplicate_groups),
        "duplicate_row_count": duplicate_rows,
        "cross_split_group_count": cross_split_groups,
        "cross_split_row_count": cross_split_rows,
        "groups": duplicate_groups,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    artifacts = {}
    for path in (train_path, valid_path, sequence_path, middle_path, smoke_path, duplicates_path):
        artifacts[path.name] = {"path": str(path), "sha256": sha256_file(path)}
    lock = {
        "protocol": "capsul_specialist_data_lock_v1",
        "created_unix": time.time(),
        "source": str(source),
        "source_sha256": observed_sha,
        "rows": len(frame),
        "unique_ids": int(frame["ID"].nunique()),
        "split_counts": counts,
        "labels": LABELS,
        "label_order_sha256": stable_hash((i, label) for i, label in enumerate(LABELS)),
        "id_order_sha256": stable_hash((i, value) for i, value in enumerate(frame["ID"])),
        "id_split_sha256": stable_hash(zip(frame["ID"], frame["Dataset"])),
        "label_mapping": "raw values 1 and 2 map to binary positive; 0 maps to negative",
        "training_scope_contains_test_targets": False,
        "test_targets_exported": False,
        "exact_duplicate_audit": {
            "groups": len(duplicate_groups), "rows": duplicate_rows,
            "cross_split_groups": cross_split_groups, "cross_split_rows": cross_split_rows,
        },
        "artifacts": artifacts,
    }
    lock_path = out / "DATA_LOCK.json"
    lock_path.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(lock, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
