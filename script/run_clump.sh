#!/bin/bash

REF_BFILE="MR/reference/EUR"

EXPOSURE_DIR="MR/exposure_all/UKB_1_22_all_rsid"

OUT_DIR="MR/exposure_all/clump_results_all"

mkdir -p "$OUT_DIR"

gz_files=("$EXPOSURE_DIR"/*.gz)
if [ ! -e "${gz_files[0]}" ]; then
    echo "Wrong: No .gz files found in the input directory: $EXPOSURE_DIR"
    exit 1
fi

echo "PLINK Clumping..."
echo "----------------------------------------"

for gz_file in "$EXPOSURE_DIR"/*.gz; do
    base_name=$(basename "$gz_file")
    
    prefix="${base_name%.gz}"
    prefix="${prefix%.txt}"
    prefix="${prefix%.tsv}"

    echo "Processing: $base_name ..."

    plink \
      --bfile "$REF_BFILE" \
      --clump "$gz_file" \
      --clump-snp-field rsid \
      --clump-p1 5e-8 \
      --clump-r2 0.001 \
      --clump-kb 10000 \
      --out "$OUT_DIR/$prefix"

    echo "finished: $prefix"
    echo "----------------------------------------"
done

echo "Done with all files."