#!/usr/bin/env python3"

import os
import gzip
import pandas as pd
from glob import glob

input_dir = "MR/exposure/UKB_1_22_rsid"
file_pattern = os.path.join(input_dir, "*.gz")

gz_files = glob(file_pattern)

for file_path in gz_files:
    print(f"processing: {os.path.basename(file_path)}")
    
    df = pd.read_csv(file_path, compression='gzip', sep='\t')
    
    if 'rsid' in df.columns:
        df.rename(columns={'rsid': 'SNP'}, inplace=True)
        
        df.to_csv(file_path, sep='\t', index=False, compression='gzip')
        print(f"  done: rsid → SNP")
    else:
        print(f"   'rsid' not found, current column names: {list(df.columns)[:5]}...")

print("Done！")