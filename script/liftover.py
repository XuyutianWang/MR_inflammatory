import argparse
from pathlib import Path

import pandas as pd
from pyliftover import LiftOver


def norm_chr(c):
    return str(c).replace("chr", "")


def lift_variant(lo, chrom, pos):
    chrom = norm_chr(chrom)
    pos0 = int(pos) - 1

    for c in [chrom, f"chr{chrom}"]:
        try:
            hits = lo.convert_coordinate(c, pos0)
        except Exception:
            hits = None

        if hits:
            new_chrom, new_pos0, _, _ = hits[0]
            return str(new_chrom).replace("chr", ""), int(new_pos0) + 1

    return None, None


def process_file(input_file, output_file, lo):
    df = pd.read_csv(input_file, sep="\t")

    required = ["SNP", "beta", "se", "effect_allele", "other_allele", "chr", "pos", "pval"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{input_file} missing columns: {missing}")

    out_rows = []
    total_snps = len(df)
    success_count = 0

    for _, row in df.iterrows():
        lifted_chr, lifted_pos = lift_variant(lo, row["chr"], row["pos"])
        
        # Copy the original row data
        new_row = row.to_dict()
        
        if lifted_chr is not None:
            # Create new columns for GRCh38 coordinates
            new_row["chr_38"] = lifted_chr
            new_row["pos_38"] = lifted_pos
            new_row["SNP_38"] = f"{lifted_chr}_{lifted_pos}_{row['effect_allele']}_{row['other_allele']}"
            success_count += 1
        else:
            # Assign None (will become NA in pandas) if mapping fails
            new_row["chr_38"] = None
            new_row["pos_38"] = None
            new_row["SNP_38"] = None

        out_rows.append(new_row)

    out = pd.DataFrame(out_rows)
    # Ensure original column order is maintained, with new columns at the end
    out.to_csv(output_file, sep="\t", index=False)
    
    return total_snps, success_count


def main():
    parser = argparse.ArgumentParser(description="Lift MR incidence subtype files to GRCh38.")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--chain", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    lo = LiftOver(args.chain)

    for f in sorted(input_dir.iterdir()):
        if f.is_file() and f.suffix.lower() in {".tsv", ".txt", ".csv"}:
            out_file = output_dir / f.name
            total, success = process_file(f, out_file, lo)
            
            # Calculate and print the success rate
            success_rate = (success / total * 100) if total > 0 else 0
            print(f"Processed: {f.name} | Success Rate: {success_rate:.2f}% ({success}/{total})")


if __name__ == "__main__":
    main()