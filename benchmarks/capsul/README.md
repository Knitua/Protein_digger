# CAPSUL official-split benchmark

This directory contains a retrospective evaluation of two Nucleus Specialist variants on the official CAPSUL split: 14,126 train, 3,027 validation and 3,028 test proteins. Raw label values `1` and `2` are positive.

The supervised heads are freshly initialized. Only generic ESM2-3B pretrained representations are reused; no DeepLoc supervised checkpoint or external supervised labels are transferred. Learning rate and early stopping use validation BCE only. All six checkpoint/configuration locks are frozen before target-free test prediction, and final metrics use a fixed `probability > 0.5` rule.

Files:

- `PROTOCOL.md`: model, training, leakage guard and evaluation protocol;
- `CAPSUL_COMPARISON.md`: compact paper-style comparison and paired bootstrap results;
- `metrics.json`: complete machine-readable metrics for all six runs;
- `train_capsul_paper_protocol.py`: training implementation;
- `prepare_scopes.py`, `verify_cache.py`, `select_learning_rates.py`: data-scope, cache and validation-selection tooling;
- `final_evaluate.py`: guarded final test evaluation;
- `test_protocol.py`: unit checks for clipping, windows, label conversion, thresholds and test access.
- `scripts/nuclear_sota/`: the minimal frozen model, cache and window-store modules required by the benchmark scripts.

The CAPSUL dataset, representation caches and checkpoints are not redistributed.
The benchmark implementation retains its upstream CC BY-NC-SA 4.0 license in `LICENSE`.
