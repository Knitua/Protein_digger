#!/usr/bin/env python3
"""Small deterministic protocol tests run before GPU work."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
from sklearn.metrics import f1_score, precision_recall_fscore_support

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
CODE_ROOT = Path(os.environ.get("NUCLEAR_SOTA_CODE_ROOT", HERE)).resolve()
sys.path.insert(0, str(CODE_ROOT))

from prepare_scopes import middleclip
from scripts.nuclear_sota.prepare_cache import window_coordinates


def run() -> None:
    sequence = "A" * 511 + "B" * 100 + "C" * 511
    clipped = middleclip(sequence)
    assert len(clipped) == 1022
    assert clipped == "A" * 511 + "C" * 511
    assert middleclip("ACD") == "ACD"

    for length in (1, 1022, 1023, 2046, 5000, 40000):
        coords = window_coordinates(length, 1022, 768)
        covered = np.zeros(length, dtype=bool)
        for start, end in coords:
            assert 0 <= start < end <= length
            covered[start:end] = True
        assert covered.all()

    raw = np.asarray([[0, 1, 2], [2, 0, 1]])
    assert np.array_equal((raw > 0).astype(np.int8), [[0, 1, 1], [1, 0, 1]])
    probability = np.asarray([0.49, 0.5, 0.5000001, 1.0])
    assert np.array_equal(probability > 0.5, [False, False, True, True])

    target = np.asarray([[1, 0], [1, 1], [0, 1]])
    prediction = np.asarray([[1, 0], [0, 1], [0, 1]])
    ours = precision_recall_fscore_support(
        target, prediction, average="micro", zero_division=0
    )[2]
    assert abs(ours - f1_score(target.reshape(-1), prediction.reshape(-1))) < 1e-12
    print("protocol unit checks passed")


if __name__ == "__main__":
    run()
