# CAPSUL paper-protocol Nucleus specialist evaluation

This run evaluates two freshly initialized ESM2-3B layer-36 heads on the official
CAPSUL split (14,126 train / 3,027 validation / 3,028 test). Raw label values 1
and 2 are positive.

## Models

- `middleclip1022`: first 511 plus last 511 residues for sequences longer than
  1,022; 20-query LabelQueryHead (projection 512, four heads) plus the
  three-query Nucleus specialist (projection 384, four heads).
- `fullwindow_mil`: complete-sequence windows of 1,022 residues with stride 768;
  the same residue head and specialist, with label-specific attention over windows.

Only the Nucleus logit receives the specialist residual. Both variants start from
freshly initialized supervised heads and use only the generic ESM2-3B pretrained
encoder representations. No DeepLoc checkpoint or external supervision is used.

## Selection and training

The objective is unweighted `BCEWithLogitsLoss`. AdamW uses weight decay 0.01;
there is no scheduler. Training runs for at most 200 epochs with validation-BCE
early stopping (patience 5) and gradient clipping at 1.0. Seed 42 screens learning
rates 4e-5, 1e-4, and 2e-4 using validation BCE only. Ties prefer 1e-4, then 4e-5,
then 2e-4. The selected learning rate is frozen separately for each input variant,
after which seeds 42, 43, and 44 are trained.

The training process can read only labeled train and labeled validation scope files.
The all-split cache source contains ID, split, and sequence only. All six checkpoint
locks must exist before target-free test prediction begins. Test targets are unlocked
only by `final_evaluate.py` after all predictions are hashed.

## Evaluation

All decisions use the fixed rule `probability > 0.5`. Results include per-label
precision, recall, F1, MCC, ROC-AUC and PR-AUC, overall micro/macro metrics, nuclear
family metrics, counts, hierarchy inconsistency, three-seed mean±SD, and a paired
protein bootstrap comparison of the two seed-42 variants.

This is a retrospective, validation-only-selection evaluation on the official
CAPSUL split. It is not described as an organizationally blind test because the
CAPSUL test set has been used by earlier work. The paper comparison is footnoted:
the paper selected hyperparameters using test micro-F1, whereas this run does not.
