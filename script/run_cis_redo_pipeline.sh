#!/usr/bin/env bash
# Label pQTLs, isolate cis-pQTLs, apply association/F filters, and LD-clump.
# Final lead-SNP tables are written to exposure_all/lead_snps_cis_redo.
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw_dir="$project_dir/exposure_all/UKB_1_22_all_rsid"
work_dir="$project_dir/exposure_all/cis_redo_work"
labelled_dir="$work_dir/01_labelled_cis_trans"
cis_dir="$work_dir/02_cis_only"
final_dir="$project_dir/exposure_all/lead_snps_cis_redo"

mkdir -p "$work_dir"

python3 "$project_dir/script/label_pqtl_cis_trans.py" \
  --input-dir "$raw_dir" \
  --labelled-dir "$labelled_dir" \
  --map-output "$work_dir/protein_gene_map.csv" \
  --audit-output "$work_dir/protein_gene_map_audit.csv"

Rscript "$project_dir/script/filter_cis_lead_snps.R" "$labelled_dir" "$cis_dir"

python3 "$project_dir/script/select_clump_cis_instruments.py" \
  --input-dir "$cis_dir" \
  --output-dir "$final_dir" \
  --reference "$project_dir/reference/EUR"
