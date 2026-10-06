#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Iterable, Sequence

import h5py
import numpy as np
import pandas as pd
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.nuclear_sota.common import (
    CAPSUL_LABELS,
    DEEPLOC_LABELS,
    atomic_json_dump,
    load_capsul,
    load_deeploc,
    sequence_sha256,
    sha256_file,
    stable_records_hash,
)


def window_coordinates(length: int, window_size: int = 1022, stride: int = 768) -> list[tuple[int, int]]:
    if length <= 0:
        raise ValueError("Sequence length must be positive")
    if window_size <= 0 or stride <= 0 or stride > window_size:
        raise ValueError("Require 0 < stride <= window_size")
    if length <= window_size:
        return [(0, length)]
    starts = list(range(0, length - window_size + 1, stride))
    final_start = length - window_size
    if starts[-1] != final_start:
        starts.append(final_start)
    coordinates = [(start, min(length, start + window_size)) for start in starts]
    coverage = np.zeros(length, dtype=np.int16)
    for start, end in coordinates:
        coverage[start:end] += 1
    if not np.all(coverage > 0):
        raise AssertionError("Window policy left uncovered residues")
    return coordinates


def parse_layers(value: str) -> list[int]:
    layers = sorted({int(item.strip()) for item in value.split(",") if item.strip()})
    if not layers:
        raise ValueError("At least one encoder layer is required")
    return layers


def load_id_subset(path: str | Path | None) -> set[str] | None:
    if path is None:
        return None
    path = Path(path)
    if path.suffix == ".npz":
        arrays = np.load(path, allow_pickle=False)
        for key in ("ids", "train", "inner_valid"):
            if key in arrays.files:
                values = arrays[key].astype(str).tolist()
                break
        else:
            raise ValueError("ID subset NPZ needs ids, train, or inner_valid")
    else:
        frame = pd.read_csv(path)
        if frame.shape[1] == 0:
            raise ValueError("ID subset table has no columns")
        values = frame.iloc[:, 0].astype(str).tolist()
    selected = set(values)
    if not selected:
        raise ValueError("ID subset is empty")
    return selected


def select_extraction_frame(
    frame: pd.DataFrame,
    id_column: str,
    id_list: str | Path | None,
    maximum_samples: int | None,
    shard_index: int,
    num_shards: int,
    repair_id_list: str | Path | None = None,
) -> tuple[pd.DataFrame, set[str] | None]:
    """Select a normal cache shard or exact IDs in their original global shard.

    Normal subset jobs retain the historical subset-then-shard behavior. Repair
    jobs deliberately shard the complete dataset first so a small repair list
    cannot silently move an accession into a different cache file.
    """
    if not 0 <= shard_index < num_shards:
        raise ValueError("Require 0 <= shard_index < num_shards")
    repair_ids = load_id_subset(repair_id_list)
    if repair_ids is not None:
        if id_list is not None or maximum_samples is not None:
            raise ValueError("Repair mode cannot be combined with --id-list or --maximum-samples")
        available = set(frame[id_column].astype(str))
        missing = sorted(repair_ids.difference(available))
        if missing:
            raise KeyError(f"Repair list contains {len(missing)} unknown proteins; first={missing[0]}")
        global_shard = frame.iloc[np.arange(len(frame)) % num_shards == shard_index].copy()
        in_shard = set(global_shard[id_column].astype(str)).intersection(repair_ids)
        wrong_shard = sorted(repair_ids.difference(in_shard))
        if wrong_shard:
            raise ValueError(
                f"Repair IDs must belong to global shard {shard_index}; "
                f"{len(wrong_shard)} do not, first={wrong_shard[0]}"
            )
        return global_shard.loc[global_shard[id_column].astype(str).isin(repair_ids)].copy(), repair_ids

    selected_ids = load_id_subset(id_list)
    if selected_ids is not None:
        missing = sorted(selected_ids.difference(frame[id_column].astype(str)))
        if missing:
            raise KeyError(f"ID subset contains {len(missing)} unknown proteins; first={missing[0]}")
        frame = frame.loc[frame[id_column].astype(str).isin(selected_ids)].copy()
    if maximum_samples:
        frame = frame.head(maximum_samples)
    return frame.iloc[np.arange(len(frame)) % num_shards == shard_index].copy(), None


def symmetric_int8(array: np.ndarray, group_size: int = 64) -> tuple[np.ndarray, np.ndarray]:
    """Groupwise symmetric int8 over embedding channels."""
    if array.ndim != 2 or array.shape[1] % group_size:
        raise ValueError("Embedding dimension must be divisible by the int8 group size")
    grouped = array.reshape(array.shape[0], array.shape[1] // group_size, group_size)
    maximum = np.max(np.abs(grouped), axis=(0, 2))
    scale = np.where(maximum > 0, maximum / 127.0, 1.0).astype(np.float32)
    quantized = np.clip(
        np.rint(grouped / scale[None, :, None]), -127, 127
    ).astype(np.int8)
    return quantized.reshape(array.shape), scale


def symmetric_int4_pack(array: np.ndarray, group_size: int = 64) -> tuple[np.ndarray, np.ndarray]:
    """Groupwise symmetric signed int4, packed as two channel values per byte."""
    if array.ndim != 2 or array.shape[1] % group_size or array.shape[1] % 2:
        raise ValueError("Embedding dimension must be divisible by group size and two")
    grouped = array.reshape(array.shape[0], array.shape[1] // group_size, group_size)
    maximum = np.max(np.abs(grouped), axis=2)
    scale = np.where(maximum > 0, maximum / 7.0, 1.0).astype(np.float32)
    quantized = np.clip(
        np.rint(grouped / scale[:, :, None]), -7, 7
    ).astype(np.int8).reshape(array.shape)
    nibble = np.bitwise_and(quantized, 0x0F).astype(np.uint8)
    packed = nibble[:, 0::2] | (nibble[:, 1::2] << 4)
    return packed, scale


def symmetric_int4_unpack(packed: np.ndarray, original_dimension: int) -> np.ndarray:
    packed = np.asarray(packed, dtype=np.uint8)
    if packed.ndim != 2 or packed.shape[1] * 2 != original_dimension:
        raise ValueError("Packed int4 shape does not match original dimension")
    values = np.empty((packed.shape[0], original_dimension), dtype=np.int8)
    low = (packed & 0x0F).astype(np.int8)
    high = ((packed >> 4) & 0x0F).astype(np.int8)
    low[low >= 8] -= 16
    high[high >= 8] -= 16
    values[:, 0::2] = low
    values[:, 1::2] = high
    return values


def dataset_records(dataset: str, csv_path: str | Path):
    if dataset == "deeploc":
        frame = load_deeploc(csv_path)
        return frame, "ACC", "Sequence", "Partition", DEEPLOC_LABELS
    if dataset == "capsul":
        # Cache preparation is sequence-only: never parse CAPSUL test labels.
        frame = pd.read_csv(csv_path, usecols=["ID", "Dataset", "Fasta Sequence"])
        frame["Dataset"] = frame["Dataset"].astype(str).str.lower()
        observed = frame["Dataset"].value_counts().to_dict()
        expected = {"train": 14_126, "valid": 3_027, "test": 3_028}
        if any(int(observed.get(split, -1)) != count for split, count in expected.items()):
            raise AssertionError(f"Unexpected CAPSUL sequence split counts: {observed}")
        if frame["ID"].duplicated().any():
            raise AssertionError("CAPSUL sequence IDs must be unique")
        frame["ID"] = frame["ID"].astype(str)
        return frame, "ID", "Fasta Sequence", "Dataset", CAPSUL_LABELS
    raise ValueError(f"Unsupported dataset: {dataset}")


def audit_dataset(
    dataset: str,
    csv_path: str | Path,
    window_size: int,
    stride: int,
    output: str | Path,
) -> dict:
    frame, id_column, sequence_column, partition_column, labels = dataset_records(dataset, csv_path)
    lengths = frame[sequence_column].astype(str).str.len().to_numpy()
    partition_counts = {
        str(key): int(value) for key, value in frame[partition_column].value_counts().sort_index().items()
    }
    windows = [window_coordinates(int(length), window_size, stride) for length in lengths]
    report = {
        "dataset": dataset,
        "csv": str(Path(csv_path).resolve()),
        "csv_sha256": sha256_file(csv_path),
        "n_samples": int(len(frame)),
        "n_unique_ids": int(frame[id_column].nunique()),
        "partition_counts": partition_counts,
        "labels": labels,
        "label_order_sha256": stable_records_hash((index, label) for index, label in enumerate(labels)),
        "id_partition_sha256": stable_records_hash(
            (identifier, partition)
            for identifier, partition in zip(frame[id_column].astype(str), frame[partition_column].astype(str))
        ),
        "id_sequence_sha256": stable_records_hash(
            (identifier, sequence_sha256(sequence))
            for identifier, sequence in zip(frame[id_column].astype(str), frame[sequence_column].astype(str))
        ),
        "sequence_lengths": {
            "minimum": int(lengths.min()),
            "median": float(np.median(lengths)),
            "mean": float(lengths.mean()),
            "maximum": int(lengths.max()),
            "le_1022": int((lengths <= 1022).sum()),
            "1023_2046": int(((lengths > 1022) & (lengths <= 2046)).sum()),
            "gt_2046": int((lengths > 2046).sum()),
        },
        "window_policy": {"window_size": window_size, "stride": stride, "full_coverage": True},
        "window_count": {
            "total": int(sum(map(len, windows))),
            "maximum_per_protein": int(max(map(len, windows))),
        },
        "created_unix": time.time(),
    }
    atomic_json_dump(report, output)
    return report


def load_esm2(model_path: str, device: torch.device):
    import esm

    # PyTorch 2.6 changed torch.load(weights_only=True) by default, while the
    # official fair-esm checkpoint includes an argparse.Namespace config.
    # This path is allowed only for a separately SHA256-verified official
    # checkpoint supplied by --model-path.
    location = Path(model_path)
    model_name = location.stem
    model_data = torch.load(str(location), map_location="cpu", weights_only=False)
    regression_location = location.parent / f"{model_name}-contact-regression.pt"
    regression_data = (
        torch.load(str(regression_location), map_location="cpu", weights_only=False)
        if regression_location.exists()
        else None
    )
    model, alphabet = esm.pretrained.load_model_and_alphabet_core(
        model_name, model_data, regression_data
    )
    model.eval().to(device)
    return model, alphabet, alphabet.get_batch_converter()


def write_window(
    layer_group: h5py.Group,
    name: str,
    representation: np.ndarray,
    quantization: str,
    start: int,
    end: int,
    quantization_qc: bool = False,
    int8_group_size: int = 64,
) -> float | None:
    cosine = None
    if quantization == "int8":
        stored, scale = symmetric_int8(representation.astype(np.float32), int8_group_size)
        if quantization_qc:
            grouped = stored.astype(np.float32).reshape(
                stored.shape[0], stored.shape[1] // int8_group_size, int8_group_size
            )
            reconstructed = (grouped * scale[None, :, None]).reshape(stored.shape)
            numerator = float(np.sum(representation.astype(np.float32) * reconstructed))
            denominator = float(
                np.linalg.norm(representation.astype(np.float32))
                * np.linalg.norm(reconstructed)
            )
            cosine = numerator / max(denominator, 1e-12)
        dataset = layer_group.create_dataset(name, data=stored, compression="lzf")
        dataset.attrs["scale"] = scale
        dataset.attrs["quantization_group_size"] = int8_group_size
        dataset.attrs["quantization_bits"] = 8
    elif quantization == "int4":
        stored, scale = symmetric_int4_pack(representation.astype(np.float32), int8_group_size)
        if quantization_qc:
            quantized = symmetric_int4_unpack(stored, representation.shape[1])
            grouped = quantized.astype(np.float32).reshape(
                quantized.shape[0], quantized.shape[1] // int8_group_size, int8_group_size
            )
            reconstructed = (grouped * scale[:, :, None]).reshape(representation.shape)
            numerator = float(np.sum(representation.astype(np.float32) * reconstructed))
            denominator = float(
                np.linalg.norm(representation.astype(np.float32)) * np.linalg.norm(reconstructed)
            )
            cosine = numerator / max(denominator, 1e-12)
        dataset = layer_group.create_dataset(name, data=stored, compression="lzf")
        dataset.attrs["scale"] = scale
        dataset.attrs["quantization_group_size"] = int8_group_size
        dataset.attrs["quantization_bits"] = 4
        dataset.attrs["original_dimension"] = representation.shape[1]
        scale_name = f"{name}__scale"
        layer_group.create_dataset(scale_name, data=scale.astype(np.float16), compression="lzf")
        dataset.attrs["scale_dataset"] = scale_name
    elif quantization in {"fp16", "float16"}:
        dataset = layer_group.create_dataset(name, data=representation.astype(np.float16), compression="lzf")
    else:
        raise ValueError("Supported quantization modes are int4, int8 and fp16")
    dataset.attrs["start"] = start
    dataset.attrs["end"] = end
    return cosine


def extract_esm2(
    dataset: str,
    csv_path: str | Path,
    model_path: str,
    layers: Sequence[int],
    output: str | Path,
    quantization: str,
    window_size: int,
    stride: int,
    device_name: str,
    resume: bool,
    maximum_samples: int | None,
    shard_index: int,
    num_shards: int,
    token_budget: int,
    max_windows: int,
    quantization_qc_samples: int,
    int8_group_size: int,
    id_list: str | Path | None,
    repair_id_list: str | Path | None = None,
    repair_audit_output: str | Path | None = None,
    source_commit: str | None = None,
) -> None:
    frame, id_column, sequence_column, _, _ = dataset_records(dataset, csv_path)
    frame, repair_ids = select_extraction_frame(
        frame,
        id_column,
        id_list,
        maximum_samples,
        shard_index,
        num_shards,
        repair_id_list,
    )
    if repair_ids is not None and not resume:
        raise ValueError("Repair mode requires --resume and an existing cache shard")
    if repair_ids is not None and repair_audit_output is None:
        raise ValueError("Repair mode requires --repair-audit-output")
    if repair_ids is not None and not source_commit:
        raise ValueError("Repair mode requires --source-commit for provenance")
    device = torch.device(device_name)
    model, alphabet, batch_converter = load_esm2(model_path, device)
    if device.type == "cuda":
        model.half()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if resume else "w"
    with h5py.File(output, mode) as handle:
        handle.attrs["format"] = "nuclear_sota_window_cache_v1"
        handle.attrs["dataset"] = dataset
        handle.attrs["encoder"] = Path(model_path).name
        handle.attrs["layers"] = json.dumps(list(layers))
        handle.attrs["window_size"] = window_size
        handle.attrs["stride"] = stride
        handle.attrs["quantization"] = quantization
        if quantization in {"int4", "int8"}:
            handle.attrs["quantization_scheme"] = f"groupwise_symmetric_{quantization}"
            handle.attrs["quantization_group_size"] = int8_group_size
        handle.attrs["shard_index"] = shard_index
        handle.attrs["num_shards"] = num_shards
        records = []
        remaining: dict[str, int] = {}
        repaired = []
        for _, row in frame.iterrows():
            identifier = str(row[id_column])
            sequence = str(row[sequence_column])
            digest = sequence_sha256(sequence)
            force_repair = repair_ids is not None and identifier in repair_ids
            if identifier in handle:
                group = handle[identifier]
                if str(group.attrs.get("sequence_sha256", "")) != digest:
                    raise ValueError(f"Sequence hash changed for existing cache entry {identifier}")
                previous_complete = bool(group.attrs.get("complete", False))
                if previous_complete and not force_repair:
                    continue
                del handle[identifier]
                handle.flush()
                if force_repair:
                    repaired.append(
                        {
                            "identifier": identifier,
                            "previous_complete": previous_complete,
                            "sequence_sha256": digest,
                        }
                    )
            elif force_repair:
                raise KeyError(f"Repair target {identifier} is absent from cache shard {shard_index}")
            group = handle.create_group(identifier)
            group.attrs["sequence_sha256"] = digest
            group.attrs["original_length"] = len(sequence)
            group.attrs["complete"] = False
            coordinates = window_coordinates(len(sequence), window_size, stride)
            remaining[identifier] = len(coordinates)
            for window_number, (start, end) in enumerate(coordinates):
                records.append(
                    {
                        "identifier": identifier,
                        "window_number": window_number,
                        "start": start,
                        "end": end,
                        "sequence": sequence[start:end],
                    }
                )
        records.sort(key=lambda record: (len(record["sequence"]), record["identifier"], record["start"]))
        rng = np.random.default_rng(20260830 + shard_index)
        qc_count = min(quantization_qc_samples, len(records)) if quantization in {"int4", "int8"} else 0
        qc_indices = set(rng.choice(len(records), size=qc_count, replace=False).tolist()) if qc_count else set()
        qc_values = []
        batches = []
        batch = []
        maximum_length = 0
        for record_index, record in enumerate(records):
            candidate_maximum = max(maximum_length, len(record["sequence"]))
            if batch and (
                candidate_maximum * (len(batch) + 1) > token_budget or len(batch) >= max_windows
            ):
                batches.append(batch)
                batch = []
                maximum_length = 0
            record["record_index"] = record_index
            batch.append(record)
            maximum_length = max(maximum_length, len(record["sequence"]))
        if batch:
            batches.append(batch)
        for batch_number, batch in enumerate(batches, 1):
            converted = [
                (f"{record['identifier']}:{record['start']}-{record['end']}", record["sequence"])
                for record in batch
            ]
            _, _, tokens = batch_converter(converted)
            tokens = tokens.to(device)
            with torch.inference_mode(), torch.autocast(
                device_type=device.type, dtype=torch.float16, enabled=device.type == "cuda"
            ):
                result = model(tokens, repr_layers=list(layers), return_contacts=False)
            for item_index, record in enumerate(batch):
                identifier = record["identifier"]
                for layer in layers:
                    representation = result["representations"][layer][
                        item_index, 1 : len(record["sequence"]) + 1
                    ]
                    layer_group = handle[identifier].require_group(f"layer_{layer}")
                    cosine = write_window(
                        layer_group,
                        f"window_{record['window_number']:05d}",
                        representation.float().cpu().numpy(),
                        quantization,
                        record["start"],
                        record["end"],
                        record["record_index"] in qc_indices,
                        int8_group_size,
                    )
                    if cosine is not None:
                        qc_values.append(cosine)
                remaining[identifier] -= 1
                if remaining[identifier] == 0:
                    handle[identifier].attrs["complete"] = True
            if batch_number == 1 or batch_number % 20 == 0 or batch_number == len(batches):
                handle.flush()
                print(
                    json.dumps(
                        {
                            "shard": shard_index,
                            "batch": batch_number,
                            "batches": len(batches),
                            "windows_complete": int(sum(len(item) for item in batches[:batch_number])),
                            "windows_total": len(records),
                        }
                    ),
                    flush=True,
                )
        if qc_values and repair_ids is None:
            handle.attrs["quantization_qc_cosine_count"] = len(qc_values)
            handle.attrs["quantization_qc_cosine_mean"] = float(np.mean(qc_values))
            handle.attrs["quantization_qc_cosine_minimum"] = float(np.min(qc_values))
            handle.flush()
        if repair_ids is not None:
            if len(repaired) != len(repair_ids):
                raise AssertionError("Not every requested cache group was repaired")
            repair_payload = {
                "format": "nuclear_sota_cache_repair_audit_v1",
                "timestamp_unix": time.time(),
                "cache": str(output.resolve()),
                "dataset": dataset,
                "encoder": Path(model_path).name,
                "source_commit": source_commit,
                "prepare_cache_sha256": sha256_file(Path(__file__)),
                "layers": list(layers),
                "shard_index": shard_index,
                "num_shards": num_shards,
                "quantization": quantization,
                "quantization_qc_cosine_count": len(qc_values),
                "quantization_qc_cosine_mean": float(np.mean(qc_values)) if qc_values else None,
                "quantization_qc_cosine_minimum": float(np.min(qc_values)) if qc_values else None,
                "repaired": repaired,
            }
            handle.attrs["last_repair_timestamp_unix"] = repair_payload["timestamp_unix"]
            handle.attrs["last_repair_count"] = len(repaired)
            handle.flush()
            atomic_json_dump(repair_payload, repair_audit_output)


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit or create full-coverage residue caches")
    parser.add_argument("--dataset", choices=("deeploc", "capsul"), required=True)
    parser.add_argument("--encoder", required=True)
    parser.add_argument("--layers", required=True)
    parser.add_argument("--window-policy", choices=("full1022_stride768",), required=True)
    parser.add_argument("--quantization", choices=("int4", "int8", "fp16"), required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--model-path")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--maximum-samples", type=int)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--num-shards", type=int, default=1)
    parser.add_argument("--token-budget", type=int, default=8192)
    parser.add_argument("--max-windows", type=int, default=8)
    parser.add_argument("--quantization-qc-samples", type=int, default=100)
    parser.add_argument("--int8-group-size", type=int, default=64)
    parser.add_argument("--id-list", help="Optional first-column CSV or NPZ ID subset")
    parser.add_argument(
        "--repair-id-list",
        help="Exact first-column CSV/NPZ IDs to replace in their original global shard",
    )
    parser.add_argument("--repair-audit-output", help="Required JSON audit path in repair mode")
    parser.add_argument("--source-commit", help="Git commit recorded in repair provenance")
    parser.add_argument(
        "--exclude-split",
        help="Development-only cache option; exclude one official partition without reading locked labels",
    )
    arguments = parser.parse_args()
    layers = parse_layers(arguments.layers)
    if arguments.audit_only:
        report = audit_dataset(arguments.dataset, arguments.csv, 1022, 768, arguments.output)
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    if not arguments.model_path:
        parser.error("--model-path is required unless --audit-only is set")
    if arguments.encoder.startswith("esmc"):
        if arguments.repair_id_list is not None or arguments.repair_audit_output is not None:
            raise NotImplementedError("Targeted repair is currently implemented only for ESM2 caches")
        from scripts.nuclear_sota.prepare_hf_cache import extract_hf_sequence

        extract_hf_sequence(
            arguments.dataset,
            arguments.csv,
            arguments.encoder,
            arguments.model_path,
            layers,
            arguments.output,
            arguments.quantization,
            1022,
            768,
            arguments.device,
            arguments.resume,
            arguments.maximum_samples,
            arguments.shard_index,
            arguments.num_shards,
            arguments.token_budget,
            arguments.max_windows,
            arguments.quantization_qc_samples,
            arguments.int8_group_size,
            arguments.exclude_split,
            arguments.id_list,
        )
        return
    if not arguments.encoder.startswith("esm2"):
        raise NotImplementedError(f"Unsupported sequence cache encoder: {arguments.encoder}")
    if arguments.exclude_split is not None:
        raise ValueError("--exclude-split is currently reserved for the Stage 2 ESMC development trial")
    extract_esm2(
        arguments.dataset,
        arguments.csv,
        arguments.model_path,
        layers,
        arguments.output,
        arguments.quantization,
        1022,
        768,
        arguments.device,
        arguments.resume,
        arguments.maximum_samples,
        arguments.shard_index,
        arguments.num_shards,
        arguments.token_budget,
        arguments.max_windows,
        arguments.quantization_qc_samples,
        arguments.int8_group_size,
        arguments.id_list,
        arguments.repair_id_list,
        arguments.repair_audit_output,
        arguments.source_commit,
    )


if __name__ == "__main__":
    main()
