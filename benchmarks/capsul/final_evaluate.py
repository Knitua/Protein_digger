#!/usr/bin/env python3
"""Unlock CAPSUL test targets only after all six checkpoints are frozen."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    matthews_corrcoef,
    precision_recall_fscore_support,
    roc_auc_score,
)


DATA_SHA = "50459272402c7e939108428298d444ee7c977d9a5630148665a98397ed0309f2"
LABELS = [
    "Nucleus", "Nuclear Membrane", "Nucleoli", "Nucleoplasm", "Cytoplasm",
    "Cytosol", "Cytoskeleton", "Centrosome", "Mitochondria",
    "Endoplasmic Reticulum", "Golgi Apparatus",
    "Plasma Membrane/Cell Membrane", "Endosome", "Lipid droplet",
    "Lysosome/Vacuole", "Peroxisome", "Vesicle", "Primary Cilium",
    "Secreted Proteins", "Sperm",
]
NUCLEAR = LABELS[:4]
PAPER = {
    "ESM-C 600M": {"Nucleus": 0.649, "micro_f1": 0.495, "macro_f1": 0.263},
    "Published best by nuclear label": {
        "Nuclear Membrane": 0.037, "Nucleoli": 0.203, "Nucleoplasm": 0.643,
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_auc(fn: Callable, target: np.ndarray, probability: np.ndarray) -> float:
    return float("nan") if np.unique(target).size < 2 else float(fn(target, probability))


def binary(target: np.ndarray, prediction: np.ndarray, probability: np.ndarray) -> dict[str, Any]:
    precision, recall, f1, _ = precision_recall_fscore_support(
        target, prediction, average="binary", zero_division=0
    )
    tp = int(((target == 1) & (prediction == 1)).sum())
    fp = int(((target == 0) & (prediction == 1)).sum())
    fn = int(((target == 1) & (prediction == 0)).sum())
    return {
        "precision": float(precision), "recall": float(recall), "f1": float(f1),
        "mcc": float(matthews_corrcoef(target, prediction)),
        "roc_auc": safe_auc(roc_auc_score, target, probability),
        "pr_auc": safe_auc(average_precision_score, target, probability),
        "tp": tp, "fp": fp, "fn": fn, "support": int(target.sum()),
        "predicted_positive": int(prediction.sum()),
    }


def metrics(target: np.ndarray, probability: np.ndarray) -> dict[str, Any]:
    prediction = (probability > 0.5).astype(np.int8)
    micro_p, micro_r, micro_f1, _ = precision_recall_fscore_support(
        target, prediction, average="micro", zero_division=0
    )
    macro_p, macro_r, macro_f1, _ = precision_recall_fscore_support(
        target, prediction, average="macro", zero_division=0
    )
    per_label = {
        label: binary(target[:, i], prediction[:, i], probability[:, i])
        for i, label in enumerate(LABELS)
    }
    nuclear_target, nuclear_prediction = target[:, :4], prediction[:, :4]
    n_micro = precision_recall_fscore_support(
        nuclear_target, nuclear_prediction, average="micro", zero_division=0
    )
    n_macro = precision_recall_fscore_support(
        nuclear_target, nuclear_prediction, average="macro", zero_division=0
    )
    sub_positive = prediction[:, 1:4].any(axis=1)
    inconsistent = sub_positive & (prediction[:, 0] == 0)
    return {
        "threshold": 0.5, "threshold_operator": ">",
        "overall": {
            "micro_precision": float(micro_p), "micro_recall": float(micro_r),
            "micro_f1": float(micro_f1), "macro_precision": float(macro_p),
            "macro_recall": float(macro_r), "macro_f1": float(macro_f1),
        },
        "per_label": per_label,
        "nuclear_family": {
            "micro_precision": float(n_micro[0]), "micro_recall": float(n_micro[1]),
            "micro_f1": float(n_micro[2]), "macro_precision": float(n_macro[0]),
            "macro_recall": float(n_macro[1]), "macro_f1": float(n_macro[2]),
            "hierarchy_inconsistent_proteins": int(inconsistent.sum()),
            "sub_label_positive_proteins": int(sub_positive.sum()),
            "hierarchy_inconsistency_rate": (
                float(inconsistent.sum() / sub_positive.sum()) if sub_positive.sum() else 0.0
            ),
        },
    }


def scalar_view(item: dict[str, Any]) -> dict[str, float]:
    result = {
        "micro_f1": item["overall"]["micro_f1"],
        "macro_f1": item["overall"]["macro_f1"],
        "nuclear_family_micro_f1": item["nuclear_family"]["micro_f1"],
        "nuclear_family_macro_f1": item["nuclear_family"]["macro_f1"],
    }
    for label in NUCLEAR:
        for metric in ("precision", "recall", "f1", "mcc", "roc_auc", "pr_auc"):
            result[f"{label}_{metric}"] = item["per_label"][label][metric]
    return result


def aggregate(seed_metrics: list[dict[str, Any]]) -> dict[str, Any]:
    views = [scalar_view(item) for item in seed_metrics]
    return {
        key: {
            "mean": float(np.mean([item[key] for item in views])),
            "sd": float(np.std([item[key] for item in views], ddof=1)),
        }
        for key in views[0]
    }


def bootstrap(
    target: np.ndarray, a: np.ndarray, b: np.ndarray, samples: int = 10000
) -> dict[str, Any]:
    rng = np.random.default_rng(20261006)
    functions = {
        "overall_micro_f1": lambda t, p: precision_recall_fscore_support(
            t, (p > 0.5), average="micro", zero_division=0
        )[2],
        "overall_macro_f1": lambda t, p: precision_recall_fscore_support(
            t, (p > 0.5), average="macro", zero_division=0
        )[2],
        "Nucleus_f1": lambda t, p: precision_recall_fscore_support(
            t[:, 0], (p[:, 0] > 0.5), average="binary", zero_division=0
        )[2],
        "nuclear_family_micro_f1": lambda t, p: precision_recall_fscore_support(
            t[:, :4], (p[:, :4] > 0.5), average="micro", zero_division=0
        )[2],
    }
    output = {}
    for name, fn in functions.items():
        observed = float(fn(target, b) - fn(target, a))
        values = np.empty(samples, dtype=np.float32)
        for i in range(samples):
            index = rng.integers(0, len(target), size=len(target))
            values[i] = fn(target[index], b[index]) - fn(target[index], a[index])
        output[name] = {
            "delta_fullwindow_minus_middleclip": observed,
            "ci95_low": float(np.quantile(values, 0.025)),
            "ci95_high": float(np.quantile(values, 0.975)),
            "samples": samples,
        }
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--union-csv", required=True)
    args = parser.parse_args()
    root = Path(args.run_root).resolve()
    union = Path(args.union_csv).resolve()
    if sha256_file(union) != DATA_SHA:
        raise AssertionError("CAPSUL data SHA mismatch")

    # Hard gate: no target is read until every checkpoint and target-free prediction is frozen.
    records: dict[str, dict[int, dict[str, Any]]] = {}
    prediction_archives: dict[str, dict[int, Any]] = {}
    for variant in ("middleclip1022", "fullwindow_mil"):
        records[variant], prediction_archives[variant] = {}, {}
        for seed in (42, 43, 44):
            run = root / "formal" / variant / f"seed{seed}"
            lock = run / "CHECKPOINT_LOCK.json"
            complete = run / "RUN_COMPLETE.json"
            prediction = run / "test_predictions.npz"
            manifest = run / "test_predictions.manifest.json"
            for path in (lock, complete, prediction, manifest):
                if not path.exists():
                    raise FileNotFoundError(path)
            lock_payload = json.loads(lock.read_text(encoding="utf-8"))
            manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
            if lock_payload.get("test_targets_used") is not False:
                raise AssertionError(f"Invalid checkpoint lock: {lock}")
            if manifest_payload.get("targets_read") is not False:
                raise AssertionError(f"Prediction scope violation: {manifest}")
            if sha256_file(prediction) != manifest_payload["prediction_sha256"]:
                raise AssertionError(f"Prediction hash mismatch: {prediction}")
            records[variant][seed] = {
                "checkpoint_lock": lock_payload,
                "run_complete": json.loads(complete.read_text(encoding="utf-8")),
            }
            prediction_archives[variant][seed] = np.load(prediction, allow_pickle=False)

    frame = pd.read_csv(union)
    frame["ID"] = frame["ID"].astype(str)
    frame["Dataset"] = frame["Dataset"].astype(str).str.lower()
    test = frame.loc[frame["Dataset"] == "test"].copy()
    if len(test) != 3028 or test["ID"].nunique() != 3028:
        raise AssertionError("Expected 3,028 unique test IDs")
    target = (test[LABELS].to_numpy(dtype=float) > 0).astype(np.int8)
    expected_ids = test["ID"].tolist()
    results: dict[str, Any] = {
        "protocol": "capsul_official_split_retrospective_validation_only_v1",
        "dataset_sha256": DATA_SHA, "test_count": len(test),
        "threshold": ">0.5", "paper_references": PAPER,
        "variants": {},
    }
    probabilities: dict[str, dict[int, np.ndarray]] = {}
    for variant in records:
        results["variants"][variant] = {"seeds": {}}
        probabilities[variant] = {}
        seed_items = []
        for seed in (42, 43, 44):
            archive = prediction_archives[variant][seed]
            ids = archive["ids"].astype(str).tolist()
            probability = archive["probability"].astype(np.float32)
            if ids != expected_ids:
                raise AssertionError(f"Test ID order mismatch for {variant}/seed{seed}")
            if probability.shape != (3028, 20) or not np.isfinite(probability).all():
                raise AssertionError(f"Invalid test probabilities for {variant}/seed{seed}")
            item = metrics(target, probability)
            item["checkpoint"] = records[variant][seed]["run_complete"]
            results["variants"][variant]["seeds"][str(seed)] = item
            seed_items.append(item)
            probabilities[variant][seed] = probability
        results["variants"][variant]["mean_sd"] = aggregate(seed_items)
    results["paired_bootstrap_seed42"] = bootstrap(
        target, probabilities["middleclip1022"][42], probabilities["fullwindow_mil"][42]
    )

    metrics_path = root / "metrics.json"
    metrics_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# CAPSUL official-split Nucleus specialist comparison", "",
        "This is a retrospective evaluation on the official CAPSUL split. Model and learning-rate selection used validation BCE only; the test set was evaluated once after all six checkpoints were frozen. The CAPSUL paper reports test-selected hyperparameters, so its values have a favorable selection bias relative to ours.", "",
        "All values below use the fixed decision rule `probability > 0.5`.", "",
        "| Method | Seed | Nucleus F1 | Nuclear Membrane F1 | Nucleoli F1 | Nucleoplasm F1 | Micro-F1 | Macro-F1 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
        "| ESM-C 600M (paper) | paper | 0.649 | – | 0.091 | 0.621 | 0.495 | 0.263 |",
        "| Published best by nuclear label | paper | – | 0.037 | 0.203 | 0.643 | – | – |",
    ]
    for variant in ("middleclip1022", "fullwindow_mil"):
        for seed in (42, 43, 44):
            item = results["variants"][variant]["seeds"][str(seed)]
            per = item["per_label"]
            lines.append(
                f"| {variant} | {seed} | {per['Nucleus']['f1']:.3f} | "
                f"{per['Nuclear Membrane']['f1']:.3f} | {per['Nucleoli']['f1']:.3f} | "
                f"{per['Nucleoplasm']['f1']:.3f} | {item['overall']['micro_f1']:.3f} | "
                f"{item['overall']['macro_f1']:.3f} |"
            )
        agg = results["variants"][variant]["mean_sd"]
        lines.append(
            f"| {variant} | mean±SD | {agg['Nucleus_f1']['mean']:.3f}±{agg['Nucleus_f1']['sd']:.3f} | "
            f"{agg['Nuclear Membrane_f1']['mean']:.3f}±{agg['Nuclear Membrane_f1']['sd']:.3f} | "
            f"{agg['Nucleoli_f1']['mean']:.3f}±{agg['Nucleoli_f1']['sd']:.3f} | "
            f"{agg['Nucleoplasm_f1']['mean']:.3f}±{agg['Nucleoplasm_f1']['sd']:.3f} | "
            f"{agg['micro_f1']['mean']:.3f}±{agg['micro_f1']['sd']:.3f} | "
            f"{agg['macro_f1']['mean']:.3f}±{agg['macro_f1']['sd']:.3f} |"
        )
    lines.extend(["", "## Paired bootstrap (seed 42)", ""])
    for name, item in results["paired_bootstrap_seed42"].items():
        lines.append(
            f"- {name}: fullwindow − middleclip = {item['delta_fullwindow_minus_middleclip']:.4f} "
            f"(95% CI {item['ci95_low']:.4f} to {item['ci95_high']:.4f})."
        )
    comparison = root / "CAPSUL_COMPARISON.md"
    comparison.write_text("\n".join(lines) + "\n", encoding="utf-8")
    complete = {
        "status": "complete", "metrics": str(metrics_path),
        "metrics_sha256": sha256_file(metrics_path), "comparison": str(comparison),
        "comparison_sha256": sha256_file(comparison),
        "statement": "Retrospective validation-only selection evaluation on the official CAPSUL split",
    }
    (root / "COMPLETE.json").write_text(
        json.dumps(complete, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(complete, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
