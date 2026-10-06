#!/usr/bin/env python3
"""Train or predict the CAPSUL Nucleus specialist without reading test targets."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
import time
from pathlib import Path
from typing import Any, Sequence

HERE = Path(__file__).resolve().parent
CODE_ROOT = Path(os.environ.get("NUCLEAR_SOTA_CODE_ROOT", HERE)).resolve()
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

from scripts.nuclear_sota.common import CAPSUL_LABELS, sha256_file
from scripts.nuclear_sota.models import (
    GeneralWithNucleusSpecialist,
    LabelQueryHead,
    build_window_mil_head,
)
from scripts.nuclear_sota.window_store import (
    WindowResidueStore,
    make_window_batch,
    window_dynamic_batches,
)


EXPECTED_LABELS = [
    "Nucleus", "Nuclear Membrane", "Nucleoli", "Nucleoplasm", "Cytoplasm",
    "Cytosol", "Cytoskeleton", "Centrosome", "Mitochondria",
    "Endoplasmic Reticulum", "Golgi Apparatus",
    "Plasma Membrane/Cell Membrane", "Endosome", "Lipid droplet",
    "Lysosome/Vacuole", "Peroxisome", "Vesicle", "Primary Cilium",
    "Secreted Proteins", "Sperm",
]


def atomic_json(payload: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def atomic_torch(payload: Any, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, temporary)
    os.replace(temporary, path)


def stable_hash(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def read_labeled_scope(path: Path, expected_split: str) -> tuple[list[str], dict[str, np.ndarray]]:
    frame = pd.read_csv(path)
    required = {"ID", "Dataset", *EXPECTED_LABELS}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Labeled scope misses columns: {missing}")
    frame["ID"] = frame["ID"].astype(str)
    frame["Dataset"] = frame["Dataset"].astype(str).str.lower()
    if set(frame["Dataset"]) != {expected_split}:
        raise AssertionError(f"{path} is not a pure {expected_split} scope")
    if frame["ID"].duplicated().any():
        raise AssertionError(f"Duplicate IDs in {path}")
    values = frame[EXPECTED_LABELS].to_numpy(dtype=np.float32)
    if not np.isin(values, (0.0, 1.0)).all():
        raise AssertionError("Training scopes must already be binary")
    ids = frame["ID"].tolist()
    return ids, {identifier: values[i] for i, identifier in enumerate(ids)}


def read_sequence_scope(path: Path, split: str) -> list[str]:
    frame = pd.read_csv(path, usecols=["ID", "Dataset"])
    frame["ID"] = frame["ID"].astype(str)
    frame["Dataset"] = frame["Dataset"].astype(str).str.lower()
    ids = frame.loc[frame["Dataset"] == split, "ID"].tolist()
    if len(ids) != len(set(ids)):
        raise AssertionError("Prediction IDs are not unique")
    return ids


def build_model(variant: str, input_dim: int) -> torch.nn.Module:
    if list(CAPSUL_LABELS) != EXPECTED_LABELS:
        raise AssertionError("Frozen source CAPSUL label order changed")
    if variant == "middleclip1022":
        general = LabelQueryHead(
            input_dim=input_dim, projection_dim=512, output_dim=20, heads=4,
            mode="label", dropout=0.15, smoothing_radius=0,
            dct_cutoff_fraction=0.25,
        )
        return GeneralWithNucleusSpecialist(
            general, input_dim=input_dim, nucleus_index=0,
            specialist_projection_dim=384, specialist_heads=4, dropout=0.15,
        )
    if variant == "fullwindow_mil":
        return build_window_mil_head(
            input_dim=input_dim, output_dim=20, nucleus_index=0,
            model_config={
                "projection_dim": 512, "heads": 4, "query_mode": "label",
                "dropout": 0.15, "smoothing_radius": 0,
                "dct_cutoff_fraction": 0.25, "mil_aggregation": "attention",
                "nucleus_specialist": True, "specialist_projection_dim": 384,
                "specialist_heads": 4,
            },
        )
    raise ValueError(f"Unknown variant: {variant}")


def forward_batch(
    model: torch.nn.Module,
    variant: str,
    batch: Sequence[str],
    store: WindowResidueStore,
    label_map: dict[str, np.ndarray],
    device: torch.device,
) -> tuple[torch.Tensor, torch.Tensor]:
    windows, mask, owner, metadata, target = make_window_batch(
        batch, store, label_map, 20, device
    )
    if variant == "middleclip1022":
        if windows.shape[0] != len(batch) or not torch.equal(
            owner, torch.arange(len(batch), device=device)
        ):
            raise AssertionError("middleclip cache must contain exactly one window per protein")
        output = model(windows, mask)
    else:
        output = model(windows, mask, owner, metadata, len(batch))
    return output.logits, target


@torch.inference_mode()
def evaluate_bce(
    model: torch.nn.Module, variant: str, ids: Sequence[str],
    store: WindowResidueStore, label_map: dict[str, np.ndarray],
    device: torch.device, token_budget: int, max_proteins: int, seed: int,
) -> float:
    model.eval()
    total = 0.0
    denominator = 0
    for batch in window_dynamic_batches(
        ids, store.costs, token_budget, max_proteins, seed, 0, False
    ):
        with torch.autocast(
            device_type=device.type, dtype=torch.float16,
            enabled=device.type == "cuda",
        ):
            logits, target = forward_batch(model, variant, batch, store, label_map, device)
            loss_sum = F.binary_cross_entropy_with_logits(logits, target, reduction="sum")
        total += float(loss_sum.cpu())
        denominator += int(target.numel())
    return total / denominator


@torch.inference_mode()
def predict(
    model: torch.nn.Module, variant: str, ids: Sequence[str],
    store: WindowResidueStore, device: torch.device,
    token_budget: int, max_proteins: int, seed: int,
) -> np.ndarray:
    model.eval()
    dummy = {identifier: np.zeros(20, dtype=np.float32) for identifier in ids}
    outputs = []
    ordered = []
    for batch in window_dynamic_batches(
        ids, store.costs, token_budget, max_proteins, seed, 0, False
    ):
        with torch.autocast(
            device_type=device.type, dtype=torch.float16,
            enabled=device.type == "cuda",
        ):
            logits, _ = forward_batch(model, variant, batch, store, dummy, device)
        outputs.append(torch.sigmoid(logits).float().cpu().numpy())
        ordered.extend(batch)
    matrix = np.concatenate(outputs, axis=0)
    by_id = {identifier: matrix[i] for i, identifier in enumerate(ordered)}
    return np.stack([by_id[identifier] for identifier in ids]).astype(np.float32)


def train(args: argparse.Namespace) -> None:
    out = Path(args.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    train_ids, train_map = read_labeled_scope(Path(args.train_csv), "train")
    valid_ids, valid_map = read_labeled_scope(Path(args.valid_csv), "valid")
    if args.smoke_limit:
        train_ids = train_ids[: args.smoke_limit]
        valid_ids = valid_ids[: args.smoke_limit]
        train_map = {k: train_map[k] for k in train_ids}
        valid_map = {k: valid_map[k] for k in valid_ids}
    if set(train_ids).intersection(valid_ids):
        raise AssertionError("Train/validation IDs overlap")
    label_map = {**train_map, **valid_map}
    set_seed(args.seed)
    device = torch.device(args.device)
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats()
    store = WindowResidueStore(args.cache, train_ids + valid_ids, layer=36)
    if args.variant == "middleclip1022" and any(
        store.window_counts[identifier] != 1 for identifier in train_ids + valid_ids
    ):
        raise AssertionError("middleclip cache contains a multi-window protein")
    model = build_model(args.variant, store.dimension).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.learning_rate, weight_decay=0.01
    )
    scaler = torch.cuda.amp.GradScaler(enabled=device.type == "cuda")
    config = {
        "variant": args.variant, "seed": args.seed,
        "learning_rate": args.learning_rate, "weight_decay": 0.01,
        "scheduler": None, "max_epochs": args.epochs, "patience": args.patience,
        "gradient_clip": 1.0, "loss": "unweighted_BCEWithLogitsLoss",
        "fixed_threshold": 0.5, "threshold_operator": ">",
        "train_csv": str(Path(args.train_csv).resolve()),
        "train_csv_sha256": sha256_file(args.train_csv),
        "valid_csv": str(Path(args.valid_csv).resolve()),
        "valid_csv_sha256": sha256_file(args.valid_csv),
        "cache_glob": args.cache, "cache_layer": 36,
        "test_targets_available_to_process": False,
        "smoke_limit": args.smoke_limit,
    }
    atomic_json(config, out / "config.json")
    config_hash = stable_hash(config)
    best_path = out / "best_checkpoint.pt"
    last_path = out / "last_checkpoint.pt"
    best_loss, best_epoch, stale, start_epoch = float("inf"), -1, 0, 0
    history: list[dict[str, Any]] = []
    if args.resume and last_path.exists():
        saved = torch.load(last_path, map_location=device, weights_only=False)
        if saved["config_hash"] != config_hash:
            raise AssertionError("Resume config differs from original")
        model.load_state_dict(saved["model"])
        optimizer.load_state_dict(saved["optimizer"])
        scaler.load_state_dict(saved["scaler"])
        best_loss, best_epoch = float(saved["best_loss"]), int(saved["best_epoch"])
        stale, start_epoch = int(saved["stale"]), int(saved["next_epoch"])
        history = list(saved["history"])
    started = time.time()
    for epoch in range(start_epoch, args.epochs):
        model.train()
        train_loss_sum, train_items = 0.0, 0
        batches = window_dynamic_batches(
            train_ids, store.costs, args.train_token_budget,
            args.max_proteins, args.seed, epoch, True,
        )
        for batch in batches:
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(
                device_type=device.type, dtype=torch.float16,
                enabled=device.type == "cuda",
            ):
                logits, target = forward_batch(
                    model, args.variant, batch, store, label_map, device
                )
                loss = F.binary_cross_entropy_with_logits(logits, target)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optimizer)
            scaler.update()
            train_loss_sum += float(loss.detach().cpu()) * int(target.numel())
            train_items += int(target.numel())
        valid_loss = evaluate_bce(
            model, args.variant, valid_ids, store, label_map, device,
            args.eval_token_budget, args.max_proteins, args.seed,
        )
        record = {
            "epoch": epoch + 1,
            "train_bce": train_loss_sum / train_items,
            "valid_bce": valid_loss,
        }
        history.append(record)
        atomic_json({"epochs": history}, out / "training_history.json")
        print(json.dumps(record, sort_keys=True), flush=True)
        if valid_loss < best_loss:
            best_loss, best_epoch, stale = valid_loss, epoch + 1, 0
            atomic_torch({
                "model": model.state_dict(), "variant": args.variant,
                "seed": args.seed, "learning_rate": args.learning_rate,
                "best_epoch": best_epoch, "best_valid_bce": best_loss,
                "config_hash": config_hash, "input_dim": store.dimension,
                "labels": EXPECTED_LABELS,
            }, best_path)
        else:
            stale += 1
        atomic_torch({
            "model": model.state_dict(), "optimizer": optimizer.state_dict(),
            "scaler": scaler.state_dict(), "best_loss": best_loss,
            "best_epoch": best_epoch, "stale": stale,
            "next_epoch": epoch + 1, "history": history,
            "config_hash": config_hash,
        }, last_path)
        if stale >= args.patience:
            break
    if not best_path.exists():
        raise RuntimeError("No best checkpoint was written")
    complete = {
        "protocol": "capsul_paper_protocol_training_v1", "status": "complete",
        "variant": args.variant, "seed": args.seed,
        "learning_rate": args.learning_rate, "best_epoch": best_epoch,
        "best_valid_bce": best_loss, "checkpoint": str(best_path),
        "checkpoint_sha256": sha256_file(best_path), "config_hash": config_hash,
        "train_count": len(train_ids), "valid_count": len(valid_ids),
        "test_targets_used": False, "test_predictions_generated": False,
        "elapsed_seconds": time.time() - started,
        "peak_gpu_memory_bytes": int(torch.cuda.max_memory_allocated()) if device.type == "cuda" else 0,
    }
    atomic_json(complete, out / "RUN_COMPLETE.json")
    atomic_json({
        "checkpoint": str(best_path), "checkpoint_sha256": complete["checkpoint_sha256"],
        "config_hash": config_hash, "validation_only_selection": True,
        "test_targets_used": False,
    }, out / "CHECKPOINT_LOCK.json")
    store.close()


def prediction(args: argparse.Namespace) -> None:
    checkpoint_path = Path(args.checkpoint).resolve()
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    if checkpoint["variant"] != args.variant:
        raise AssertionError("Checkpoint variant mismatch")
    ids = read_sequence_scope(Path(args.sequence_csv), args.predict_split)
    if args.smoke_limit:
        ids = ids[: args.smoke_limit]
    device = torch.device(args.device)
    store = WindowResidueStore(args.cache, ids, layer=36)
    model = build_model(args.variant, int(checkpoint["input_dim"])).to(device)
    model.load_state_dict(checkpoint["model"])
    probability = predict(
        model, args.variant, ids, store, device,
        args.eval_token_budget, args.max_proteins, int(checkpoint["seed"]),
    )
    if probability.shape != (len(ids), 20) or not np.isfinite(probability).all():
        raise AssertionError("Invalid prediction matrix")
    output = Path(args.prediction_output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output, ids=np.asarray(ids), probability=probability,
        labels=np.asarray(EXPECTED_LABELS), split=np.asarray([args.predict_split]),
    )
    atomic_json({
        "protocol": "capsul_unlabeled_prediction_v1", "split": args.predict_split,
        "count": len(ids), "unique_ids": len(set(ids)), "nan_count": int(np.isnan(probability).sum()),
        "checkpoint": str(checkpoint_path), "checkpoint_sha256": sha256_file(checkpoint_path),
        "sequence_scope": str(Path(args.sequence_csv).resolve()),
        "sequence_scope_sha256": sha256_file(args.sequence_csv),
        "targets_read": False, "prediction_file": str(output),
        "prediction_sha256": sha256_file(output),
    }, output.with_suffix(".manifest.json"))
    store.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("train", "predict"))
    parser.add_argument("--variant", choices=("middleclip1022", "fullwindow_mil"), required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--train-token-budget", type=int, default=8192)
    parser.add_argument("--eval-token-budget", type=int, default=12288)
    parser.add_argument("--max-proteins", type=int, default=8)
    parser.add_argument("--smoke-limit", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--train-csv")
    parser.add_argument("--valid-csv")
    parser.add_argument("--out-dir")
    parser.add_argument("--checkpoint")
    parser.add_argument("--sequence-csv")
    parser.add_argument("--predict-split", choices=("train", "valid", "test"), default="test")
    parser.add_argument("--prediction-output")
    args = parser.parse_args()
    if args.mode == "train":
        for field in ("train_csv", "valid_csv", "out_dir"):
            if getattr(args, field) is None:
                parser.error(f"--{field.replace('_', '-')} is required for train")
        train(args)
    else:
        for field in ("checkpoint", "sequence_csv", "prediction_output"):
            if getattr(args, field) is None:
                parser.error(f"--{field.replace('_', '-')} is required for predict")
        prediction(args)


if __name__ == "__main__":
    main()
