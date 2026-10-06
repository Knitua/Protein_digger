from __future__ import annotations

import glob
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import h5py
import numpy as np
import torch

from scripts.nuclear_sota.prepare_cache import symmetric_int4_unpack


@dataclass
class ProteinWindows:
    arrays: list[np.ndarray]
    scales: list[np.ndarray | float | None]
    coordinates: list[tuple[int, int]]
    original_length: int


class WindowResidueStore:
    """Read nuclear_sota_window_cache_v1 shards without merging large files."""

    def __init__(
        self,
        paths: str | Sequence[str],
        accessions: Sequence[str],
        layer: int,
        repair_overlays: str | Sequence[str] | None = None,
    ):
        base_paths = self._resolve_paths(paths)
        overlay_paths = self._resolve_paths(repair_overlays) if repair_overlays else []
        if not base_paths:
            raise FileNotFoundError(f"No window cache shards match: {paths}")
        duplicate_paths = sorted(set(base_paths).intersection(overlay_paths))
        if duplicate_paths:
            raise ValueError(f"A cache path cannot be both base and repair overlay: {duplicate_paths[0]}")
        self.paths = overlay_paths + base_paths
        self.overlay_paths = overlay_paths
        self._overlay_count = len(overlay_paths)
        self.layer = int(layer)
        self._handles = [h5py.File(path, "r") for path in self.paths]
        self.location: dict[str, int] = {}
        self.costs: dict[str, int] = {}
        self.window_counts: dict[str, int] = {}
        self.original_lengths: dict[str, int] = {}
        self.shadowed_accessions: set[str] = set()
        self.dimension = 0
        requested = list(dict.fromkeys(map(str, accessions)))
        if not requested:
            raise ValueError("Window cache requires at least one accession")
        key_sets: list[set[str]] = []
        for shard_index, handle in enumerate(self._handles):
            if str(handle.attrs.get("format", "")) not in {
                "nuclear_sota_window_cache_v1",
                "nuclear_sota_structural_window_cache_v1",
            }:
                raise ValueError(f"Not a nuclear SOTA window cache: {self.paths[shard_index]}")
            advertised_layers = json.loads(str(handle.attrs["layers"]))
            if self.layer not in advertised_layers:
                raise ValueError(f"Layer {self.layer} is absent from {self.paths[shard_index]}")
            key_sets.append(set(map(str, handle.keys())))

        overlay_sets = key_sets[: self._overlay_count]
        base_sets = key_sets[self._overlay_count :]
        for role, sets in (("overlay", overlay_sets), ("base", base_sets)):
            seen: set[str] = set()
            for keys in sets:
                duplicate = sorted(seen.intersection(keys))
                if duplicate:
                    raise ValueError(
                        f"Duplicate accession across {role} cache shards: {duplicate[0]}"
                    )
                seen.update(keys)
        overlay_accessions = set().union(*overlay_sets) if overlay_sets else set()
        base_accessions = set().union(*base_sets) if base_sets else set()
        unshadowed = sorted(overlay_accessions.difference(base_accessions))
        if unshadowed:
            raise ValueError(
                "Repair overlay entries must shadow an accession in the base cache; "
                f"first unshadowed accession={unshadowed[0]}"
            )
        self.shadowed_accessions = overlay_accessions.copy()

        for accession in requested:
            overlay_locations = [
                index
                for index, keys in enumerate(overlay_sets)
                if accession in keys
            ]
            if overlay_locations:
                shard_index = overlay_locations[0]
            else:
                base_locations = [
                    self._overlay_count + index
                    for index, keys in enumerate(base_sets)
                    if accession in keys
                ]
                if not base_locations:
                    continue
                shard_index = base_locations[0]
            handle = self._handles[shard_index]
            group = handle[accession]
            if not bool(group.attrs.get("complete", False)):
                continue
            layer_group = group[f"layer_{self.layer}"]
            names = sorted(
                name for name in layer_group.keys() if not name.endswith("__scale")
            )
            if not names:
                raise ValueError(f"No windows for {accession} layer {self.layer}")
            coordinates = [
                (
                    int(layer_group[name].attrs["start"]),
                    int(layer_group[name].attrs["end"]),
                )
                for name in names
            ]
            original_length = int(group.attrs["original_length"])
            coverage = np.zeros(original_length, dtype=np.bool_)
            for start, end in coordinates:
                coverage[start:end] = True
            if not coverage.all():
                raise AssertionError(f"Cache windows do not fully cover {accession}")
            first = layer_group[names[0]]
            first_dimension = int(
                first.attrs.get("original_dimension", first.shape[1])
            )
            if self.dimension == 0:
                self.dimension = first_dimension
            elif self.dimension != first_dimension:
                raise ValueError("Embedding dimension differs across cache shards")
            self.location[accession] = shard_index
            self.costs[accession] = int(
                sum(end - start for start, end in coordinates)
            )
            self.window_counts[accession] = len(names)
            self.original_lengths[accession] = original_length
        missing = [accession for accession in requested if accession not in self.location]
        if missing:
            raise KeyError(f"Window cache lacks {len(missing)} requested proteins; first={missing[0]}")

    @staticmethod
    def _resolve_paths(paths: str | Sequence[str]) -> list[str]:
        if isinstance(paths, str):
            return sorted(glob.glob(paths))
        return [str(Path(path)) for path in paths]

    def close(self) -> None:
        for handle in self._handles:
            handle.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def read(self, accession: str) -> ProteinWindows:
        handle = self._handles[self.location[accession]]
        group = handle[accession]
        layer_group = group[f"layer_{self.layer}"]
        names = sorted(name for name in layer_group.keys() if not name.endswith("__scale"))
        arrays = []
        scales = []
        coordinates = []
        for name in names:
            dataset = layer_group[name]
            values = np.asarray(dataset)
            if "scale_dataset" in dataset.attrs:
                scale = np.asarray(layer_group[str(dataset.attrs["scale_dataset"])], dtype=np.float32)
            else:
                scale = np.asarray(dataset.attrs["scale"], dtype=np.float32) if "scale" in dataset.attrs else None
            if int(dataset.attrs.get("quantization_bits", 8)) == 4:
                values = symmetric_int4_unpack(values, int(dataset.attrs["original_dimension"]))
                values = dequantize(values, scale)
                scale = None
            arrays.append(values)
            scales.append(scale)
            coordinates.append((int(dataset.attrs["start"]), int(dataset.attrs["end"])))
        return ProteinWindows(
            arrays=arrays,
            scales=scales,
            coordinates=coordinates,
            original_length=int(group.attrs["original_length"]),
        )

    def window_records(
        self, accessions: Sequence[str]
    ) -> list[tuple[str, int, int]]:
        """Return coordinates in the exact flattened order used by make_window_batch."""
        records = []
        for accession in accessions:
            handle = self._handles[self.location[accession]]
            layer_group = handle[accession][f"layer_{self.layer}"]
            names = sorted(
                name for name in layer_group.keys() if not name.endswith("__scale")
            )
            for name in names:
                dataset = layer_group[name]
                records.append(
                    (
                        accession,
                        int(dataset.attrs["start"]),
                        int(dataset.attrs["end"]),
                    )
                )
        return records


def dequantize(values: np.ndarray, scale: np.ndarray | float | None) -> np.ndarray:
    if not np.issubdtype(values.dtype, np.integer):
        return values.astype(np.float16, copy=False)
    if scale is None:
        raise ValueError("Integer cache values require a quantization scale")
    scale_array = np.asarray(scale, dtype=np.float32)
    if scale_array.ndim == 0:
        return (values.astype(np.float32) * float(scale_array)).astype(np.float16)
    if scale_array.ndim == 2:
        if scale_array.shape[0] != values.shape[0] or values.shape[1] % scale_array.shape[1]:
            raise ValueError("Row-group scales do not align with the quantized embedding")
        group_size = values.shape[1] // scale_array.shape[1]
        grouped = values.astype(np.float32).reshape(
            values.shape[0], scale_array.shape[1], group_size
        )
        return (grouped * scale_array[:, :, None]).reshape(values.shape).astype(np.float16)
    if values.shape[1] % len(scale_array):
        raise ValueError("Groupwise quantization scale does not divide the embedding dimension")
    group_size = values.shape[1] // len(scale_array)
    grouped = values.astype(np.float32).reshape(values.shape[0], len(scale_array), group_size)
    return (grouped * scale_array[None, :, None]).reshape(values.shape).astype(np.float16)


def window_dynamic_batches(
    accessions: Sequence[str],
    costs: dict[str, int],
    token_budget: int,
    max_proteins: int,
    seed: int,
    epoch: int,
    shuffle: bool,
) -> list[list[str]]:
    ordered = sorted(accessions, key=lambda accession: (costs[accession], accession))
    batches: list[list[str]] = []
    batch: list[str] = []
    cost = 0
    for accession in ordered:
        accession_cost = costs[accession]
        if batch and (cost + accession_cost > token_budget or len(batch) >= max_proteins):
            batches.append(batch)
            batch = []
            cost = 0
        batch.append(accession)
        cost += accession_cost
    if batch:
        batches.append(batch)
    if shuffle:
        np.random.default_rng(seed + epoch).shuffle(batches)
    return batches


def make_window_batch(
    accessions: Sequence[str],
    store: WindowResidueStore,
    label_map: dict[str, np.ndarray],
    label_count: int,
    device: torch.device,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    proteins = [store.read(accession) for accession in accessions]
    maximum_window = max(array.shape[0] for protein in proteins for array in protein.arrays)
    window_count = sum(len(protein.arrays) for protein in proteins)
    quantized = all(
        np.issubdtype(array.dtype, np.integer)
        for protein in proteins
        for array in protein.arrays
    )
    features = np.zeros(
        (window_count, maximum_window, store.dimension),
        dtype=np.int8 if quantized else np.float16,
    )
    scale_count = 0
    if quantized:
        scale_lengths = {
            int(np.asarray(scale).size)
            for protein in proteins
            for scale in protein.scales
            if scale is not None
        }
        if len(scale_lengths) != 1:
            raise ValueError("All quantized windows in a batch must use the same scale grouping")
        scale_count = scale_lengths.pop()
        if store.dimension % scale_count:
            raise ValueError("Quantization groups do not divide the embedding dimension")
        group_scales = np.ones((window_count, scale_count), dtype=np.float32)
    mask = np.zeros((window_count, maximum_window), dtype=np.bool_)
    owner = np.zeros(window_count, dtype=np.int64)
    metadata = np.zeros((window_count, 5), dtype=np.float32)
    target = np.stack([label_map[accession] for accession in accessions]).astype(np.float32)
    window_index = 0
    for protein_index, protein in enumerate(proteins):
        length_denominator = max(1, protein.original_length - 1)
        for values, scale, (start, end) in zip(
            protein.arrays, protein.scales, protein.coordinates
        ):
            length = values.shape[0]
            if quantized:
                features[window_index, :length] = values
                group_scales[window_index] = np.asarray(scale, dtype=np.float32).reshape(-1)
            else:
                features[window_index, :length] = dequantize(values, scale)
            mask[window_index, :length] = True
            owner[window_index] = protein_index
            metadata[window_index] = (
                start / length_denominator,
                (end - 1) / length_denominator,
                float(start == 0),
                float(end == protein.original_length),
                np.log1p(protein.original_length) / np.log1p(40_000),
            )
            window_index += 1
    feature_tensor = torch.from_numpy(features).to(device, non_blocking=True)
    if quantized:
        group_size = store.dimension // scale_count
        feature_tensor = feature_tensor.to(torch.float32).reshape(
            window_count, maximum_window, scale_count, group_size
        )
        feature_tensor = feature_tensor * torch.from_numpy(group_scales).to(
            device, non_blocking=True
        )[:, None, :, None]
        feature_tensor = feature_tensor.reshape(
            window_count, maximum_window, store.dimension
        ).to(torch.float16)
    # CPU smoke/report paths have no CUDA autocast, so align cached fp16
    # features with the fp32 model parameters. CUDA deliberately retains fp16.
    if device.type != "cuda":
        feature_tensor = feature_tensor.float()
    return (
        feature_tensor,
        torch.from_numpy(mask).to(device, non_blocking=True),
        torch.from_numpy(owner).to(device, non_blocking=True),
        torch.from_numpy(metadata).to(device, non_blocking=True),
        torch.from_numpy(target).to(device, non_blocking=True),
    )


def make_window_availability(
    accessions: Sequence[str],
    store: WindowResidueStore,
    maximum_window: int,
    device: torch.device,
) -> torch.Tensor:
    """Return the explicit per-window structural availability aligned to batch order."""
    window_count = sum(store.window_counts[accession] for accession in accessions)
    availability = np.zeros((window_count, maximum_window), dtype=np.bool_)
    window_index = 0
    for accession in accessions:
        handle = store._handles[store.location[accession]]
        group = handle[accession]
        layer_group = group[f"layer_{store.layer}"]
        window_names = sorted(
            name for name in layer_group.keys() if not name.endswith("__scale")
        )
        structure_group = group.get("structure_mask")
        for name in window_names:
            start = int(layer_group[name].attrs["start"])
            end = int(layer_group[name].attrs["end"])
            length = end - start
            if structure_group is None:
                availability[window_index, :length] = True
            else:
                if name not in structure_group:
                    raise KeyError(f"Structure mask misses {accession}/{name}")
                structure_mask = np.asarray(structure_group[name], dtype=np.bool_)
                if structure_mask.shape != (length,):
                    raise ValueError(f"Structure mask length differs for {accession}/{name}")
                availability[window_index, :length] = structure_mask
            window_index += 1
    return torch.from_numpy(availability).to(device, non_blocking=True)
