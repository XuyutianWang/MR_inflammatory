#!/usr/bin/env python3
import gzip
import glob
import os
import sys

def main():
    # Define input and output directory paths
    input_dir = os.path.expanduser("MR/exposure_all/UKB_1_22_all_pfiltered")
    output_dir = os.path.expanduser("MR/exposure_all/UKB_1_22_all_pfiltered_5108")

    # Check if the input directory exists
    if not os.path.exists(input_dir):
        print(f"Error: Input directory {input_dir} does not exist.")
        sys.exit(1)

    # Automatically create output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Find all .gz files
    gz_files = sorted(glob.glob(os.path.join(input_dir, "*.gz")))
    
    if not gz_files:
        print(f"No .gz files found in {input_dir}")
        sys.exit(0)

    print(f"Found {len(gz_files)} files to process.\n")

    for file_path in gz_files:
        file_name = os.path.basename(file_path)
        out_path = os.path.join(output_dir, file_name)
        print(f"Processing: {file_name} -> {out_path}")

        try:
            with gzip.open(file_path, "rt") as f_in, gzip.open(out_path, "wt") as f_out:
                header_line = f_in.readline()
                if not header_line:
                    continue
                
                header = header_line.strip().split("\t")

                # Locate the LOG10P column index
                if "LOG10P" in header:
                    log10p_idx = header.index("LOG10P")
                    header[log10p_idx] = "P"  # 将表头改为 P / Rename header to P
                else:
                    print(f"  [Warning] 'LOG10P' column not found in {file_name}, skipping file.")
                    continue

                # Write updated header
                f_out.write("\t".join(header) + "\n")

                # Convert values line by line
                for line in f_in:
                    parts = line.strip().split("\t")
                    if len(parts) > log10p_idx:
                        try:
                            log10p_val = float(parts[log10p_idx])
                            # Calculate p = 10^(-LOG10P)
                            p_val = 10 ** (-log10p_val)
                            parts[log10p_idx] = f"{p_val:.6e}"
                        except ValueError:
                            # Handle potential NA or invalid values gracefully
                            pass
                    
                    f_out.write("\t".join(parts) + "\n")

        except Exception as e:
            print(f"  [Error] Failed to process {file_name}: {e}")

    print("\nAll files processed successfully!")

if __name__ == "__main__":
    main()