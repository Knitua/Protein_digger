#!/usr/bin/env bash
set -euo pipefail

python3 /root/autodl-tmp/Agent_analysis_v2/03_stage0_AB_partition/scripts/validate_stage0_release.py \
  --project-root /root/autodl-tmp/Agent_analysis_v2
