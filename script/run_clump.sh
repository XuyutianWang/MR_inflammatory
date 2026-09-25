#!/usr/bin/env bash
# LD-clump cis pQTLs with PLINK 1.9.
#
# Optional overrides, for example:
#   PLINK_BIN=/path/to/plink script/run_clump.sh
#   EXPOSURE_DIR=/path/to/cis OUT_DIR=/path/to/output script/run_clump.sh
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLINK_BIN="${PLINK_BIN:-plink}"
REF_BFILE="${REF_BFILE:-$project_dir/reference/EUR}"
EXPOSURE_DIR="${EXPOSURE_DIR:-$project_dir/exposure_all/UKB_1_22_all_rsid_cis}"
OUT_DIR="${OUT_DIR:-$project_dir/exposure_all/UKB_1_22_all_rsid_cis_clumped}"

if ! command -v "$PLINK_BIN" >/dev/null 2>&1; then
    echo "Error: PLINK was not found: $PLINK_BIN" >&2
    echo "Install PLINK 1.9 or set PLINK_BIN=/absolute/path/to/plink." >&2
    exit 127
fi

if [[ ! -f "$REF_BFILE.bed" || ! -f "$REF_BFILE.bim" || ! -f "$REF_BFILE.fam" ]]; then
    echo "Error: PLINK reference files are incomplete for: $REF_BFILE" >&2
    exit 1
fi

shopt -s nullglob
txt_files=("$EXPOSURE_DIR"/*.txt)
if (( ${#txt_files[@]} == 0 )); then
    echo "Error: no .txt files found in: $EXPOSURE_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "PLINK clumping ${#txt_files[@]} exposure file(s)..."
for txt_file in "${txt_files[@]}"; do
    base_name="$(basename "$txt_file")"
    prefix="${base_name%.txt}"

    echo "Processing: $base_name"
    "$PLINK_BIN" \
        --bfile "$REF_BFILE" \
        --clump "$txt_file" \
        --clump-snp-field rsid \
        --clump-field P \
        --clump-p1 5e-8 \
        --clump-r2 0.001 \
        --clump-kb 10000 \
        --out "$OUT_DIR/$prefix"
done

echo "Done. PLINK reports are in: $OUT_DIR"
