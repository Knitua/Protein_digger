#!/usr/bin/env python3
"""Freeze per-variant learning rates using validation BCE only."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


PRIORITY = {1e-4: 0, 4e-5: 1, 2e-4: 2}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = Path(args.screen_root).resolve()
    result = {
        "protocol": "capsul_validation_only_lr_selection_v1",
        "seed": 42, "candidates": [4e-5, 1e-4, 2e-4],
        "selection_metric": "minimum validation BCE",
        "tie_priority": [1e-4, 4e-5, 2e-4],
        "test_targets_used": False, "variants": {},
    }
    for variant in ("middleclip1022", "fullwindow_mil"):
        candidates = []
        for lr_name in ("4e-5", "1e-4", "2e-4"):
            path = root / variant / f"lr_{lr_name}" / "RUN_COMPLETE.json"
            if not path.exists():
                raise FileNotFoundError(path)
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("test_targets_used") is not False:
                raise AssertionError(f"Test-scope violation in {path}")
            candidates.append({
                "learning_rate": float(payload["learning_rate"]),
                "best_valid_bce": float(payload["best_valid_bce"]),
                "checkpoint": payload["checkpoint"],
                "checkpoint_sha256": payload["checkpoint_sha256"],
                "run_complete": str(path), "run_complete_sha256": sha256_file(path),
            })
        winner = min(
            candidates,
            key=lambda item: (round(item["best_valid_bce"], 12), PRIORITY[item["learning_rate"]]),
        )
        result["variants"][variant] = {
            "selected_learning_rate": winner["learning_rate"],
            "selected_best_valid_bce": winner["best_valid_bce"],
            "candidates": candidates,
        }
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
