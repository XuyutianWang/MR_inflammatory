#!/usr/bin/env bash
# Regenerate only MR-Egger intercept output from existing harmonised files.
# Set R_BIN if Rscript is not on PATH, for example:
# R_BIN=/home/df24679/miniconda3/envs/r_mr/bin/Rscript bash script/repair_all_egger_pleiotropy.sh
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$project_dir"

r_bin="${R_BIN:-Rscript}"
jobs=(
  "cis/harmonised_results_cis|cis/mr_results_cis"
  "cis/harmonised_results_cis_incidence|cis/mr_results_cis_incidence"
  "cis/harmonised_results_cis_sh|cis/mr_results_cis_sh"
  "replicate/decode/harmonised_F10|replicate/decode/mr_F10_survival"
  "replicate/decode/harmonised_F10_incidence|replicate/decode/mr_F10_incidence"
  "replicate/decode/harmonised_F10_sh|replicate/decode/mr_F10_sh"
  "replicate/fenland/harmonised_F10|replicate/fenland/mr_F10_survival"
  "replicate/fenland/harmonised_F10_incidence|replicate/fenland/mr_F10_incidence"
  "replicate/fenland/harmonised_F10_sh|replicate/fenland/mr_F10_sh"
)

for job in "${jobs[@]}"; do
  IFS='|' read -r harmonised_dir output_dir <<< "$job"
  "$r_bin" script/repair_egger_pleiotropy.R "$harmonised_dir" "$output_dir"
done
