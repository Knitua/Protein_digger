#!/usr/bin/env python3
"""Verify sharded target-free ESM2 caches and emit an immutable manifest."""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
from pathlib import Path

import h5py
import pandas as pd


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-glob", required=True)
    parser.add_argument("--sequence-csv", required=True)
    parser.add_argument("--variant", choices=("middleclip1022", "fullwindow_mil"), required=True)
    parser.add_argument("--model-path", required=True)
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-count", type=int, default=20181)
    parser.add_argument("--minimum-cosine", type=float, default=0.995)
    args = parser.parse_args()
    paths = [Path(item).resolve() for item in sorted(glob.glob(args.cache_glob))]
    if len(paths) != 2:
        raise AssertionError(f"Expected two cache shards, found {len(paths)}")
    frame = pd.read_csv(args.sequence_csv, usecols=["ID", "Dataset", "Fasta Sequence"])
    expected_ids = set(frame["ID"].astype(str))
    if len(frame) != args.expected_count or len(expected_ids) != args.expected_count:
        raise AssertionError("Unexpected sequence-only scope")
    observed: set[str] = set()
    qc_weighted, qc_count = 0.0, 0
    windows = 0
    rows = 0
    shards = []
    for path in paths:
        with h5py.File(path, "r") as handle:
            if str(handle.attrs.get("format", "")) != "nuclear_sota_window_cache_v1":
                raise AssertionError(f"Unexpected cache format: {path}")
            if json.loads(str(handle.attrs["layers"])) != [36]:
                raise AssertionError("Cache must contain only ESM2 layer 36")
            keys = set(map(str, handle.keys()))
            overlap = observed.intersection(keys)
            if overlap:
                raise AssertionError(f"Duplicate ID across cache shards: {next(iter(overlap))}")
            for identifier in keys:
                group = handle[identifier]
                if not bool(group.attrs.get("complete", False)):
                    raise AssertionError(f"Incomplete cache entry: {identifier}")
                layer = group["layer_36"]
                names = [name for name in layer.keys() if not name.endswith("__scale")]
                if args.variant == "middleclip1022" and len(names) != 1:
                    raise AssertionError(f"middleclip entry has {len(names)} windows: {identifier}")
                windows += len(names)
                rows += sum(int(layer[name].shape[0]) for name in names)
            observed.update(keys)
            count = int(handle.attrs.get("quantization_qc_cosine_count", 0))
            mean = float(handle.attrs.get("quantization_qc_cosine_mean", 0.0))
            qc_weighted += count * mean
            qc_count += count
            shards.append({
                "path": str(path), "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size, "ids": len(keys),
                "qc_count": count, "qc_mean_cosine": mean,
            })
    missing = expected_ids.difference(observed)
    extra = observed.difference(expected_ids)
    if missing or extra:
        raise AssertionError(f"Cache ID mismatch: missing={len(missing)}, extra={len(extra)}")
    mean_cosine = qc_weighted / qc_count if qc_count else 0.0
    if mean_cosine < args.minimum_cosine:
        raise AssertionError(f"Quantization cosine {mean_cosine} < {args.minimum_cosine}")
    source_dir = Path(args.source_dir).resolve()
    source_hashes = {}
    for path in sorted(source_dir.rglob("*.py")):
        source_hashes[str(path.relative_to(source_dir))] = sha256_file(path)
    payload = {
        "protocol": "capsul_esm2_3b_cache_manifest_v1",
        "variant": args.variant, "format": "groupwise symmetric int8",
        "layer": 36, "window_size": 1022, "stride": 768,
        "unique_ids": len(observed), "windows": windows, "encoded_rows": rows,
        "quantization_qc_count": qc_count,
        "quantization_qc_mean_cosine": mean_cosine,
        "minimum_required_cosine": args.minimum_cosine,
        "sequence_scope": str(Path(args.sequence_csv).resolve()),
        "sequence_scope_sha256": sha256_file(Path(args.sequence_csv)),
        "model_path": str(Path(args.model_path).resolve()),
        "model_sha256": sha256_file(Path(args.model_path)),
        "source_dir": str(source_dir), "source_hashes": source_hashes,
        "shards": shards,
    }
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
