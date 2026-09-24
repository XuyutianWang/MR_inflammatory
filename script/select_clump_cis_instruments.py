#!/usr/bin/env python3
"""Select strong cis-pQTLs and retain LD-independent lead instruments.

This script expects the cis-only tables created by filter_cis_lead_snps.R.
It applies P <= 5e-8 and per-SNP F = (BETA / SE)^2 >= 10 before PLINK
clumping (EUR reference, r2 < 0.001 within 10 Mb).
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        return list(reader.fieldnames), list(reader)


def as_number(value: str | None) -> float:
    try:
        return float(value) if value is not None else math.nan
    except ValueError:
        return math.nan


def write_table(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_clumped_snps(path: Path) -> set[str]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        # PLINK aligns columns with variable whitespace; re-tokenise explicitly.
        lines = [line.split() for line in handle if line.strip()]
    if not lines:
        return set()
    header = lines[0]
    if "SNP" not in header:
        raise ValueError(f"PLINK output has no SNP column: {path}")
    index = header.index("SNP")
    return {fields[index] for fields in lines[1:] if len(fields) > index}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=PROJECT_DIR / "exposure_all/UKB_1_22_all_rsid_cis")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_DIR / "exposure_all/lead_snps_cis_redo")
    parser.add_argument("--reference", type=Path, default=PROJECT_DIR / "reference/EUR")
    parser.add_argument("--plink", default="plink")
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    parser.add_argument("--f-threshold", type=float, default=10.0)
    parser.add_argument("--r2", type=float, default=0.001)
    parser.add_argument("--kb", type=int, default=10_000)
    args = parser.parse_args()

    if not args.input_dir.is_dir():
        raise SystemExit(f"Cis input directory does not exist: {args.input_dir}")
    if not args.reference.with_suffix(".bed").is_file():
        raise SystemExit(f"PLINK reference is incomplete: {args.reference}.bed is missing")
    plink = shutil.which(args.plink)
    if plink is None:
        raise SystemExit(f"PLINK executable was not found: {args.plink}")

    filtered_dir = args.output_dir / "_filtered_before_clumping"
    plink_dir = args.output_dir / "_plink"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    filtered_dir.mkdir(exist_ok=True)
    plink_dir.mkdir(exist_ok=True)

    audit: list[dict[str, str | int]] = []
    input_files = sorted(args.input_dir.glob("*.txt"))
    if not input_files:
        raise SystemExit(f"No cis .txt files found in: {args.input_dir}")

    for input_file in input_files:
        exposure_id = input_file.stem
        fields, rows = read_table(input_file)
        required = {"BETA", "SE", "P"}
        missing = required.difference(fields)
        if missing:
            raise ValueError(f"{input_file} is missing: {', '.join(sorted(missing))}")
        snp_column = "SNP" if "SNP" in fields else "rsid" if "rsid" in fields else None
        if snp_column is None:
            raise ValueError(f"{input_file} has neither SNP nor rsid")

        output_fields = list(fields)
        if "SNP" not in output_fields:
            output_fields.append("SNP")
        if "F_stat" not in output_fields:
            output_fields.append("F_stat")

        retained: list[dict[str, str]] = []
        for row in rows:
            beta, se, p_value = (as_number(row.get(column)) for column in ("BETA", "SE", "P"))
            f_stat = (beta / se) ** 2 if math.isfinite(beta) and math.isfinite(se) and se != 0 else math.nan
            snp = (row.get(snp_column) or "").strip()
            if snp and math.isfinite(p_value) and p_value <= args.p_threshold and math.isfinite(f_stat) and f_stat >= args.f_threshold:
                row["SNP"] = snp
                row["F_stat"] = f"{f_stat:.8g}"
                retained.append(row)

        filtered_file = filtered_dir / f"{exposure_id}.txt"
        write_table(filtered_file, output_fields, retained)
        lead_file = args.output_dir / f"{exposure_id}_lead_snps.txt"
        clump_prefix = plink_dir / exposure_id
        status = "no_variants_after_association_and_F_filter"
        lead_rows: list[dict[str, str]] = []
        if retained:
            command = [
                plink, "--bfile", str(args.reference), "--clump", str(filtered_file),
                "--clump-snp-field", "SNP", "--clump-field", "P", "--clump-p1", "1",
                "--clump-p2", "1", "--clump-r2", str(args.r2), "--clump-kb", str(args.kb),
                "--out", str(clump_prefix),
            ]
            completed = subprocess.run(command, text=True, capture_output=True)
            (plink_dir / f"{exposure_id}.log").write_text(completed.stdout + completed.stderr, encoding="utf-8")
            if completed.returncode != 0:
                raise RuntimeError(f"PLINK failed for {exposure_id}; see {clump_prefix}.log")
            lead_snps = parse_clumped_snps(clump_prefix.with_suffix(".clumped"))
            lead_rows = [row for row in retained if row["SNP"] in lead_snps]
            status = "clumped" if lead_rows else "no_reference_matched_lead_snp"
        write_table(lead_file, output_fields, lead_rows)
        audit.append({"exposure_id": exposure_id, "cis_rows": len(rows), "retained_before_clumping": len(retained), "lead_snps": len(lead_rows), "status": status})

    with (args.output_dir / "instrument_selection_audit.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["exposure_id", "cis_rows", "retained_before_clumping", "lead_snps", "status"], delimiter="\t")
        writer.writeheader()
        writer.writerows(audit)
    print(f"Final LD-independent cis instruments: {args.output_dir}")


if __name__ == "__main__":
    main()
