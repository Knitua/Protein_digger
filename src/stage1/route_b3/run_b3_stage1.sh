#!/usr/bin/env bash
set -euo pipefail

B3_ROOT=/root/autodl-tmp/Agent_analysis_v2/03_stage1_B/B3
PYTHON=/root/miniconda3/envs/nuclear_models/bin/python
MODEL_SCRIPT_ROOT=/root/autodl-tmp/Agent_analysis/04_Nuclear_Models/scripts
DEEPLOC_ROOT=/root/autodl-tmp/Agent_analysis/04_Nuclear_Models/DeepLoc-2.0
NLS_WEIGHTS=/root/autodl-tmp/Agent_analysis/04_Nuclear_Models/NLSExplorer_1.0/NLSEXplorer_A2KA/NLS_loc_modeltes
NLS_SCRIPT="$B3_ROOT/scripts/run_nlsexplorer_unpadded_head_inference.py"
FASTA="$B3_ROOT/00_input/b3_all_alternative_isoforms_22131.fasta"
DEEPLOC_OUT="$B3_ROOT/01_nuclear_evidence/02_model_scores/deeploc_full_22131"
NLS_OUT="$B3_ROOT/01_nuclear_evidence/02_model_scores/nlsexplorer_full_22131_unpadded_head"
LOG_DIR="$B3_ROOT/logs"

test -x "$PYTHON"
test -s "$FASTA"
test -s "$NLS_WEIGHTS"
test -s "$NLS_SCRIPT"
test ! -e "$DEEPLOC_OUT/deeploc_prott5_full_predictions.tsv"
test ! -e "$NLS_OUT/nlsexplorer_protein_scores.tsv"
mkdir -p "$DEEPLOC_OUT" "$NLS_OUT" "$LOG_DIR"

date -Iseconds > "$LOG_DIR/stage1_started_at.txt"

"$PYTHON" "$MODEL_SCRIPT_ROOT/run_deeploc_prott5_inference.py" \
  --fasta "$FASTA" \
  --deeploc-root "$DEEPLOC_ROOT" \
  --outdir "$DEEPLOC_OUT" \
  --clip-len 4000 \
  --max-tokens 8192 \
  --device cuda:0 \
  --progress-every 20 \
  > "$LOG_DIR/deeploc_full_22131.log" 2>&1 &
DEEPLOC_PID=$!

"$PYTHON" "$NLS_SCRIPT" \
  --fasta "$FASTA" \
  --nls-weights "$NLS_WEIGHTS" \
  --outdir "$NLS_OUT" \
  --window 1022 \
  --overlap 200 \
  --batch-size 4 \
  --device cuda:1 \
  --progress-every 1000 \
  > "$LOG_DIR/nlsexplorer_full_22131.log" 2>&1 &
NLS_PID=$!

wait "$DEEPLOC_PID"
wait "$NLS_PID"

"$PYTHON" "$B3_ROOT/scripts/build_b3_stage1_union.py" \
  --metadata "$B3_ROOT/00_input/b3_all_alternative_isoforms_22131.tsv" \
  --direct-uniprot "$B3_ROOT/01_nuclear_evidence/01_uniprot_direct/b3_uniprot_isoform_specific_nucleus_266.tsv" \
  --deeploc "$DEEPLOC_OUT/deeploc_prott5_full_predictions.tsv" \
  --nls "$NLS_OUT/nlsexplorer_protein_scores.tsv" \
  --out-root "$B3_ROOT" \
  > "$LOG_DIR/build_stage1_union.log" 2>&1

"$PYTHON" "$B3_ROOT/scripts/validate_b3_stage1.py" \
  --b3-root "$B3_ROOT" \
  > "$LOG_DIR/validate_stage1.log" 2>&1

date -Iseconds > "$LOG_DIR/stage1_finished_at.txt"
cat "$LOG_DIR/validate_stage1.log"
