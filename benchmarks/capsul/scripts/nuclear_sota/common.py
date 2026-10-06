from __future__ import annotations

import hashlib
import json
import math
import os
import random
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    matthews_corrcoef,
    precision_recall_fscore_support,
    roc_auc_score,
)


DEEPLOC_LABELS = [
    "Membrane",
    "Cytoplasm",
    "Nucleus",
    "Extracellular",
    "Cell membrane",
    "Mitochondrion",
    "Plastid",
    "Endoplasmic reticulum",
    "Lysosome/Vacuole",
    "Golgi apparatus",
    "Peroxisome",
]
CAPSUL_LABELS = [
    "Nucleus",
    "Nuclear Membrane",
    "Nucleoli",
    "Nucleoplasm",
    "Cytoplasm",
    "Cytosol",
    "Cytoskeleton",
    "Centrosome",
    "Mitochondria",
    "Endoplasmic Reticulum",
    "Golgi Apparatus",
    "Plasma Membrane/Cell Membrane",
    "Endosome",
    "Lipid droplet",
    "Lysosome/Vacuole",
    "Peroxisome",
    "Vesicle",
    "Primary Cilium",
    "Secreted Proteins",
    "Sperm",
]
DEEPLOC_EXCLUDED_ACCESSIONS = ("Q04656-5", "O43157", "Q9UPN3-2")
LOCKED_SPLITS = {"deeploc": {"hpa"}, "capsul": {"test"}}
CONFIRMATORY_REFERENCE_METRICS = (
    "micro_f1",
    "micro_mcc",
    "macro_f1",
    "macro_mcc",
    "nucleus_precision",
    "nucleus_recall",
    "nucleus_f1",
    "nucleus_mcc",
    "nucleus_roc_auc",
    "nucleus_pr_auc",
)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def capture_rng_state() -> dict[str, Any]:
    return {
        "python": random.getstate(),
        "numpy": np.random.get_state(),
        "torch": torch.get_rng_state(),
        "cuda": torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None,
    }


def restore_rng_state(state: Mapping[str, Any]) -> None:
    random.setstate(state["python"])
    np.random.set_state(state["numpy"])
    torch.set_rng_state(state["torch"])
    if torch.cuda.is_available() and state.get("cuda") is not None:
        torch.cuda.set_rng_state_all(state["cuda"])


def measure_inference(
    function: Callable[[], Any],
    protein_count: int,
    device: torch.device,
    split: str,
) -> tuple[Any, dict[str, Any]]:
    """Run one synchronized inference pass and return portable throughput telemetry."""
    if protein_count < 0:
        raise ValueError("protein_count must be nonnegative")
    if not split:
        raise ValueError("split must be nonempty")
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    started = time.perf_counter()
    result = function()
    if device.type == "cuda":
        torch.cuda.synchronize(device)
    elapsed = time.perf_counter() - started
    if elapsed <= 0:
        raise RuntimeError("Inference benchmark produced a nonpositive elapsed time")
    return result, {
        "split": split,
        "proteins": int(protein_count),
        "inference_seconds": float(elapsed),
        "proteins_per_second": float(protein_count / elapsed),
        "device_type": device.type,
        "cuda_synchronized": device.type == "cuda",
    }


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def sha256_path(path: str | Path) -> str:
    """Hash a file or a directory tree including stable relative filenames."""
    path = Path(path)
    if path.is_file():
        return sha256_file(path)
    if not path.is_dir():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    files = sorted(item for item in path.rglob("*") if item.is_file())
    if not files:
        raise ValueError(f"Cannot hash empty directory: {path}")
    for item in files:
        relative = item.relative_to(path).as_posix()
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(item).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def sequence_sha256(sequence: str) -> str:
    return sha256_bytes(str(sequence).strip().upper().encode("ascii"))


def stable_records_hash(records: Iterable[Sequence[Any]]) -> str:
    digest = hashlib.sha256()
    for record in sorted(tuple(map(str, row)) for row in records):
        digest.update("\t".join(record).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def git_commit(repo: str | Path = ".") -> str:
    root = Path(repo).resolve()
    command = ["git", "-C", str(root), "rev-parse", "--show-toplevel", "HEAD"]
    lines = subprocess.check_output(
        command,
        text=True,
        stderr=subprocess.DEVNULL,
    ).splitlines()
    if len(lines) != 2 or Path(lines[0]).resolve() != root:
        raise subprocess.CalledProcessError(128, command)
    return lines[1].strip()


def execution_commit_record(
    config: Mapping[str, Any], repo: str | Path = "."
) -> dict[str, Any]:
    """Resolve executed code separately from the commit that authored a reused config."""
    configured = str(config.get("source_git_commit", "")).strip() or None
    environment = os.environ.get("NUCLEAR_SOTA_SOURCE_COMMIT", "").strip() or None
    try:
        observed = git_commit(repo)
    except (OSError, subprocess.CalledProcessError):
        observed = None
    if environment and observed and environment != observed:
        raise PermissionError(
            "NUCLEAR_SOTA_SOURCE_COMMIT differs from the executing Git checkout"
        )
    execution = observed or environment or configured
    if not execution:
        raise PermissionError(
            "Cannot identify execution commit; use a Git checkout or set "
            "NUCLEAR_SOTA_SOURCE_COMMIT for a pinned source snapshot"
        )
    return {
        "execution_git_commit": execution,
        "execution_commit_source": (
            "git_checkout"
            if observed
            else "pinned_snapshot_environment"
            if environment
            else "config_fallback"
        ),
        "execution_commit_verified": bool(observed or environment),
        "config_authored_git_commit": configured,
    }


def atomic_json_dump(payload: Mapping[str, Any], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=False)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def load_deeploc(path: str | Path, strict: bool = True) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"ACC", "Partition", "Sequence", *DEEPLOC_LABELS}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"DeepLoc file is missing columns: {missing}")
    frame = frame.loc[
        ~frame["ACC"].astype(str).isin(DEEPLOC_EXCLUDED_ACCESSIONS)
    ].copy()
    frame["ACC"] = frame["ACC"].astype(str)
    frame["Partition"] = frame["Partition"].astype(int)
    frame = frame.reset_index(drop=True)
    if strict:
        if len(frame) != 28_300 or frame["ACC"].nunique() != 28_300:
            raise AssertionError(
                f"Expected 28,300 strict DeepLoc proteins, got rows={len(frame)}, "
                f"unique={frame['ACC'].nunique()}"
            )
        if sorted(frame["Partition"].unique().tolist()) != [0, 1, 2, 3, 4]:
            raise AssertionError("DeepLoc must contain exactly partitions 0..4")
    return frame


def _find_column(frame: pd.DataFrame, candidates: Sequence[str]) -> str:
    lowered = {column.lower(): column for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    raise ValueError(f"None of the expected columns are present: {candidates}")


def load_capsul(path: str | Path, include_locked_test: bool = True) -> tuple[pd.DataFrame, str, str]:
    frame = pd.read_csv(path)
    split_column = _find_column(frame, ("split", "dataset", "set"))
    sequence_column = _find_column(frame, ("Sequence", "sequence", "seq", "Fasta Sequence"))
    missing = sorted(set(CAPSUL_LABELS).difference(frame.columns))
    if missing:
        raise ValueError(f"CAPSUL file is missing labels: {missing}")
    frame = frame.copy()
    frame[split_column] = frame[split_column].astype(str).str.lower()
    expected = {"train": 14_126, "valid": 3_027, "test": 3_028}
    observed = frame[split_column].value_counts().to_dict()
    if include_locked_test:
        if any(int(observed.get(split, -1)) != size for split, size in expected.items()):
            raise AssertionError(f"Unexpected CAPSUL official splits: {observed}")
    else:
        if int(observed.get("train", -1)) != expected["train"] or int(observed.get("valid", -1)) != expected["valid"]:
            raise AssertionError(f"Unexpected CAPSUL development splits: {observed}")
        if int(observed.get("test", 0)) not in {0, expected["test"]}:
            raise AssertionError(f"Unexpected CAPSUL locked-test row count: {observed}")
    if not include_locked_test:
        frame = frame.loc[frame[split_column].isin(("train", "valid"))].copy()
        if (frame[split_column] == "test").any():
            raise AssertionError("CAPSUL test rows entered a development dataframe")
    frame[CAPSUL_LABELS] = (frame[CAPSUL_LABELS].to_numpy(dtype=float) > 0).astype(np.int8)
    return frame, split_column, sequence_column


def expected_calibration_error(
    target: np.ndarray, probability: np.ndarray, bins: int = 15
) -> float:
    target = np.asarray(target, dtype=float).reshape(-1)
    probability = np.asarray(probability, dtype=float).reshape(-1)
    edges = np.linspace(0.0, 1.0, bins + 1)
    indices = np.minimum(np.digitize(probability, edges[1:-1]), bins - 1)
    error = 0.0
    for index in range(bins):
        selected = indices == index
        if selected.any():
            error += selected.mean() * abs(target[selected].mean() - probability[selected].mean())
    return float(error)


def _safe_auc(metric, target: np.ndarray, probability: np.ndarray) -> float:
    if np.unique(target).size < 2:
        return float("nan")
    return float(metric(target, probability))


def binary_metrics(
    target: np.ndarray, prediction: np.ndarray, probability: np.ndarray
) -> dict[str, float | int]:
    target = np.asarray(target, dtype=np.int8).reshape(-1)
    prediction = np.asarray(prediction, dtype=np.int8).reshape(-1)
    probability = np.asarray(probability, dtype=float).reshape(-1)
    precision, recall, f1, _ = precision_recall_fscore_support(
        target, prediction, average="binary", zero_division=0
    )
    return {
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "mcc": float(matthews_corrcoef(target, prediction)),
        "roc_auc": _safe_auc(roc_auc_score, target, probability),
        "pr_auc": _safe_auc(average_precision_score, target, probability),
        "ece": expected_calibration_error(target, probability),
        "brier": float(brier_score_loss(target, probability)),
        "support": int(target.sum()),
        "predicted_positive": int(prediction.sum()),
    }


def multilabel_metrics_from_predictions(
    target: np.ndarray,
    prediction: np.ndarray,
    probability: np.ndarray,
    labels: Sequence[str],
    metric_label_indices: Sequence[int] | None = None,
    nucleus_index: int = 0,
) -> dict[str, Any]:
    """Compute metrics when predictions use sample-specific thresholds.

    Strict OOF and legacy fold aggregation can use a different threshold vector
    in every fold.  Keeping the already thresholded prediction avoids silently
    replacing those fold-local decisions with one pooled threshold vector.
    """

    target = np.asarray(target, dtype=np.int8)
    prediction = np.asarray(prediction, dtype=np.int8)
    probability = np.asarray(probability, dtype=float)
    if target.ndim != 2:
        raise ValueError("Target must be a two-dimensional multilabel array")
    if target.shape != prediction.shape or target.shape != probability.shape:
        raise ValueError("Target, prediction, and probability shapes do not match")
    if target.shape[1] != len(labels):
        raise ValueError("Target and label dimensions do not match")
    if not 0 <= nucleus_index < len(labels):
        raise ValueError("Nucleus index is outside the label order")
    if not np.isin(target, (0, 1)).all():
        raise ValueError("Targets must be binary")
    if not np.isin(prediction, (0, 1)).all():
        raise ValueError("Predictions must be binary")
    if not np.isfinite(probability).all() or np.any((probability < 0) | (probability > 1)):
        raise ValueError("Probabilities must be finite and lie in [0, 1]")

    indices = np.asarray(
        list(range(len(labels))) if metric_label_indices is None else metric_label_indices,
        dtype=int,
    )
    if indices.ndim != 1 or not len(indices):
        raise ValueError("At least one metric label index is required")
    if np.any(indices < 0) or np.any(indices >= len(labels)):
        raise ValueError("Metric label index is outside the label order")
    selected_target = target[:, indices]
    selected_prediction = prediction[:, indices]
    selected_probability = probability[:, indices]
    micro = binary_metrics(
        selected_target.reshape(-1), selected_prediction.reshape(-1), selected_probability.reshape(-1)
    )
    per_label = {
        labels[index]: binary_metrics(target[:, index], prediction[:, index], probability[:, index])
        for index in range(len(labels))
    }

    def finite_mean(values: Sequence[float]) -> float:
        array = np.asarray(values, dtype=float)
        finite_values = array[np.isfinite(array)]
        return float(finite_values.mean()) if len(finite_values) else float("nan")

    macro = {
        key: finite_mean([per_label[labels[index]][key] for index in indices])
        for key in ("precision", "recall", "f1", "mcc", "roc_auc", "pr_auc", "ece", "brier")
    }
    macro.update(
        {
            "support": int(selected_target.sum()),
            "predicted_positive": int(selected_prediction.sum()),
        }
    )
    return {
        "micro": micro,
        "macro": macro,
        "nucleus": per_label[labels[nucleus_index]],
        "per_label": per_label,
        "n_samples": int(target.shape[0]),
    }


def multilabel_metrics(
    target: np.ndarray,
    probability: np.ndarray,
    thresholds: np.ndarray,
    labels: Sequence[str],
    metric_label_indices: Sequence[int] | None = None,
    nucleus_index: int = 0,
) -> dict[str, Any]:
    target = np.asarray(target, dtype=np.int8)
    probability = np.asarray(probability, dtype=float)
    thresholds = np.asarray(thresholds, dtype=float)
    if target.shape != probability.shape or target.shape[1] != len(labels):
        raise ValueError("Target, probability, and label dimensions do not match")
    if thresholds.shape != (len(labels),):
        raise ValueError(f"Expected {len(labels)} thresholds, got {thresholds.shape}")
    prediction = (probability >= thresholds[None, :]).astype(np.int8)
    metrics = multilabel_metrics_from_predictions(
        target,
        prediction,
        probability,
        labels,
        metric_label_indices=metric_label_indices,
        nucleus_index=nucleus_index,
    )
    metrics["thresholds"] = {
        label: float(thresholds[index]) for index, label in enumerate(labels)
    }
    return metrics


def optimize_thresholds(
    target: np.ndarray,
    probability: np.ndarray,
    method: str,
    nucleus_index: int,
    precision_floor: float | None = None,
    grid: np.ndarray | None = None,
) -> np.ndarray:
    target = np.asarray(target, dtype=np.int8)
    probability = np.asarray(probability, dtype=float)
    grid = np.linspace(0.02, 0.98, 193) if grid is None else np.asarray(grid, dtype=float)
    thresholds = np.full(target.shape[1], 0.5, dtype=float)
    for label_index in range(target.shape[1]):
        best_score = -math.inf
        best_threshold = 0.5
        for threshold in grid:
            prediction = probability[:, label_index] >= threshold
            if method == "per_label_f1":
                _, _, score, _ = precision_recall_fscore_support(
                    target[:, label_index], prediction, average="binary", zero_division=0
                )
            elif method in {"per_label_mcc", "nucleus_mcc_precision_floor"}:
                score = matthews_corrcoef(target[:, label_index], prediction)
                if (
                    method == "nucleus_mcc_precision_floor"
                    and label_index == nucleus_index
                    and precision_floor is not None
                ):
                    precision, _, _, _ = precision_recall_fscore_support(
                        target[:, label_index], prediction, average="binary", zero_division=0
                    )
                    if precision < precision_floor:
                        continue
            else:
                raise ValueError(f"Unknown threshold method: {method}")
            if score > best_score or (score == best_score and abs(threshold - 0.5) < abs(best_threshold - 0.5)):
                best_score = float(score)
                best_threshold = float(threshold)
        thresholds[label_index] = best_threshold
    return thresholds


def promotion_score(metrics: Mapping[str, Any]) -> float:
    return float(
        0.25 * metrics["micro"]["f1"]
        + 0.25 * metrics["macro"]["f1"]
        + 0.20 * metrics["nucleus"]["f1"]
        + 0.20 * metrics["nucleus"]["mcc"]
        + 0.10 * metrics["nucleus"]["pr_auc"]
    )


def _normalize_experiment_field(value: Any) -> str:
    return "".join(
        character.lower() if character.isalnum() else "-"
        for character in str(value)
    ).strip("-")


def model_subject_id(
    stage: str,
    backbone: str,
    architecture: str,
    loss: str,
    threshold: str,
) -> str:
    """Return the dataset-independent identity shared by equivalent models."""
    return "-".join(
        _normalize_experiment_field(field)
        for field in (stage, backbone, architecture, loss, threshold)
    )


def experiment_id(
    dataset: str,
    stage: str,
    backbone: str,
    architecture: str,
    loss: str,
    threshold: str,
    seed: int,
    fold: int | str,
) -> str:
    subject = model_subject_id(stage, backbone, architecture, loss, threshold)
    return "-".join(
        (_normalize_experiment_field(dataset), subject, f"s{seed}", f"f{fold}")
    )


def validate_frozen_manifest(
    manifest: Mapping[str, Any], dataset: str, maximum_candidates: int = 3
) -> list[Mapping[str, Any]]:
    if manifest.get("protocol") != "locked_candidate_freeze_v3":
        raise PermissionError("Locked evaluation requires protocol='locked_candidate_freeze_v3'")
    if manifest.get("dataset") != dataset:
        raise PermissionError("Candidate manifest dataset does not match requested locked dataset")
    if manifest.get("status") != "frozen":
        raise PermissionError("Locked evaluation requires a manifest with status='frozen'")
    candidates = manifest.get("candidates")
    if not isinstance(candidates, list) or not 1 <= len(candidates) <= maximum_candidates:
        raise PermissionError(f"Locked evaluation allows 1..{maximum_candidates} frozen candidates")
    required = {
        "candidate_id",
        "config_path",
        "config_sha256",
        "oof_path",
        "oof_sha256",
        "locked_prediction_npz",
        "locked_prediction_sha256",
        "prediction_report_path",
        "prediction_report_sha256",
        "prediction_report_protocol",
        "provenance_summary",
    }
    candidate_ids = []
    for candidate in candidates:
        missing = required.difference(candidate)
        if missing:
            raise PermissionError(f"Frozen candidate lacks required fields: {sorted(missing)}")
        if candidate["prediction_report_protocol"] not in {
            "target_free_locked_cv_ensemble_prediction_v1",
            "target_free_locked_stacker_prediction_v1",
        }:
            raise PermissionError("Frozen candidate has an unsupported provenance protocol")
        if not isinstance(candidate["provenance_summary"], Mapping):
            raise PermissionError("Frozen candidate provenance_summary must be a mapping")
        candidate_id = str(candidate["candidate_id"])
        if (
            candidate_id in {"", ".", ".."}
            or Path(candidate_id).name != candidate_id
            or "\\" in candidate_id
        ):
            raise PermissionError(
                "Frozen candidate_id must be a single path-safe identifier"
            )
        candidate_ids.append(candidate_id)
    if len(set(candidate_ids)) != len(candidate_ids):
        raise PermissionError("Frozen candidate IDs must be unique")
    if not manifest.get("confirmation_batch"):
        raise PermissionError("Manifest must name a one-shot confirmation_batch")
    required_manifest = {
        "locked_target_path",
        "locked_target_sha256",
        "locked_id_order_sha256",
        "baseline_candidate_id",
        "bootstrap",
    }
    missing_manifest = required_manifest.difference(manifest)
    if missing_manifest:
        raise PermissionError(f"Frozen manifest lacks fields: {sorted(missing_manifest)}")
    if manifest["baseline_candidate_id"] not in candidate_ids:
        raise PermissionError("Frozen baseline_candidate_id is not one of the candidates")
    baseline_candidate_ids = manifest.get(
        "baseline_candidate_ids", [manifest["baseline_candidate_id"]]
    )
    if (
        not isinstance(baseline_candidate_ids, list)
        or not 1 <= len(baseline_candidate_ids) <= len(candidate_ids)
        or len(set(baseline_candidate_ids)) != len(baseline_candidate_ids)
        or any(item not in candidate_ids for item in baseline_candidate_ids)
    ):
        raise PermissionError(
            "Frozen baseline_candidate_ids must be unique members of the confirmation batch"
        )
    if baseline_candidate_ids[0] != manifest["baseline_candidate_id"]:
        raise PermissionError(
            "baseline_candidate_id must remain the first pre-registered reference candidate"
        )
    baseline_metric_map = manifest.get("baseline_metric_map")
    if len(baseline_candidate_ids) > 1 and not isinstance(baseline_metric_map, Mapping):
        raise PermissionError(
            "Multiple frozen baselines require a pre-registered baseline_metric_map"
        )
    if baseline_metric_map is not None:
        if not isinstance(baseline_metric_map, Mapping):
            raise PermissionError("baseline_metric_map must be a mapping")
        if (
            manifest.get("baseline_reference_protocol")
            != "pre_registered_per_metric_envelope_v1"
        ):
            raise PermissionError(
                "baseline_metric_map requires pre_registered_per_metric_envelope_v1"
            )
        expected_metrics = set(CONFIRMATORY_REFERENCE_METRICS)
        if set(baseline_metric_map) != expected_metrics:
            missing = sorted(expected_metrics.difference(baseline_metric_map))
            extra = sorted(set(baseline_metric_map).difference(expected_metrics))
            raise PermissionError(
                "baseline_metric_map must freeze every confirmatory reference metric; "
                f"missing={missing}, extra={extra}"
            )
        if any(value not in baseline_candidate_ids for value in baseline_metric_map.values()):
            raise PermissionError(
                "baseline_metric_map may reference only frozen baseline_candidate_ids"
            )
    bootstrap = manifest["bootstrap"]
    if (
        not isinstance(bootstrap, Mapping)
        or bootstrap.get("unit") != "protein"
        or not bootstrap.get("paired")
        or int(bootstrap.get("replicates", 0)) < 10_000
    ):
        raise PermissionError("Locked evaluation requires at least 10,000 paired protein bootstraps")
    return candidates
