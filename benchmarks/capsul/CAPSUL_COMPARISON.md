# CAPSUL official-split Nucleus specialist comparison

This is a retrospective evaluation on the official CAPSUL split. Model and learning-rate selection used validation BCE only; the test set was evaluated once after all six checkpoints were frozen. The CAPSUL paper reports test-selected hyperparameters, so its values have a favorable selection bias relative to ours.

All values below use the fixed decision rule `probability > 0.5`.

| Method | Seed | Nucleus F1 | Nuclear Membrane F1 | Nucleoli F1 | Nucleoplasm F1 | Micro-F1 | Macro-F1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ESM-C 600M (paper) | paper | 0.649 | – | 0.091 | 0.621 | 0.495 | 0.263 |
| Published best by nuclear label | paper | – | 0.037 | 0.203 | 0.643 | – | – |
| middleclip1022 | 42 | 0.631 | 0.000 | 0.233 | 0.607 | 0.490 | 0.295 |
| middleclip1022 | 43 | 0.628 | 0.000 | 0.215 | 0.591 | 0.481 | 0.269 |
| middleclip1022 | 44 | 0.637 | 0.026 | 0.218 | 0.640 | 0.525 | 0.339 |
| middleclip1022 | mean±SD | 0.632±0.005 | 0.009±0.015 | 0.222±0.010 | 0.612±0.025 | 0.499±0.023 | 0.301±0.035 |
| fullwindow_mil | 42 | 0.642 | 0.000 | 0.242 | 0.604 | 0.509 | 0.307 |
| fullwindow_mil | 43 | 0.638 | 0.000 | 0.141 | 0.613 | 0.519 | 0.314 |
| fullwindow_mil | 44 | 0.631 | 0.000 | 0.127 | 0.610 | 0.503 | 0.292 |
| fullwindow_mil | mean±SD | 0.637±0.005 | 0.000±0.000 | 0.170±0.062 | 0.609±0.004 | 0.511±0.008 | 0.304±0.011 |

## Paired bootstrap (seed 42)

- overall_micro_f1: fullwindow − middleclip = 0.0195 (95% CI 0.0125 to 0.0265).
- overall_macro_f1: fullwindow − middleclip = 0.0120 (95% CI 0.0027 to 0.0215).
- Nucleus_f1: fullwindow − middleclip = 0.0111 (95% CI -0.0016 to 0.0240).
- nuclear_family_micro_f1: fullwindow − middleclip = 0.0056 (95% CI -0.0049 to 0.0158).
