#!/usr/bin/env python3
import pandas as pd
import gzip
import glob
import os

clump_dir = os.path.expanduser("MR/exposure_all/clump_results_all")              # .clumped files directory
gz_dir = os.path.expanduser("MR/exposure_all/UKB_1_22_all_rsid")  # original .gz files directory
out_dir = os.path.expanduser("MR/exposure_all/lead_snps_all")     # new output directory for lead SNPs

os.makedirs(out_dir, exist_ok=True)

clump_files = glob.glob(os.path.join(clump_dir, "*.clumped"))

for c_file in clump_files:
    prefix = os.path.basename(c_file).replace(".clumped", "")
    gz_file = os.path.join(gz_dir, f"{prefix}.gz")
    
    if not os.path.exists(gz_file):
        possible_gz = glob.glob(os.path.join(gz_dir, f"{prefix}*.gz"))
        if possible_gz:
            gz_file = possible_gz[0]
        else:
            print(f"cannot find corresponding .gz file: {prefix}")
            continue

    # 1. extract Lead SNP names from .clumped file (3rd column, i.e., index 2)
    lead_snps = set()
    with open(c_file, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                lead_snps.add(parts[2])  # SNP 列

    if not lead_snps:
        print(f"[{prefix}] no significant Lead SNPs found")
        continue

    # 2. extract Lead SNP rows from the original .gz file
    extracted_rows = []
    with gzip.open(gz_file, 'rt') as gz_in:
        header_line = gz_in.readline().strip().split('\t')
        
        if len(header_line) == 1:
            gz_in.seek(0)
            df_full = pd.read_csv(gz_file, sep=r'\s+')
        else:
            df_full = pd.read_csv(gz_file, sep='\t')

    # The source files use either SNP or rsid for the variant identifier.
    if 'SNP' not in df_full.columns:
        if 'rsid' in df_full.columns:
            df_full = df_full.rename(columns={'rsid': 'SNP'})
        else:
            raise ValueError(
                f"[{prefix}] input file has no SNP or rsid column; "
                f"found: {', '.join(map(str, df_full.columns))}"
            )

    # screen Lead SNP
    df_lead = df_full[df_full['SNP'].isin(lead_snps)].copy()

    # 3. calculate R2 and F-statistic
    if 'BETA' in df_lead.columns and 'SE' in df_lead.columns and 'N' in df_lead.columns:
        df_lead['R2'] = (df_lead['BETA'] ** 2) / ((df_lead['BETA'] ** 2) + (df_lead['N'] * (df_lead['SE'] ** 2)))
        df_lead['F_stat'] = (df_lead['R2'] * (df_lead['N'] - 2)) / (1 - df_lead['R2'])

    # 4. save the extracted Lead SNPs to a new .txt file
    out_path = os.path.join(out_dir, f"{prefix}_lead_snps.txt")
    df_lead.to_csv(out_path, sep='\t', index=False)
    print(f"[{prefix}] successfully extracted {len(df_lead)} lead SNPs -> {out_path}")

print("\nAll processing completed!")