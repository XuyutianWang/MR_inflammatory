#!/usr/bin/env python3
"""Combine MR-Egger intercept tests and flag directional-pleiotropy signals.

``nominal_p_lt_0_05`` is a screening flag.  ``qval_bh_all_tests`` controls
the false-discovery rate over every successful MR-Egger test in this project.
Neither result proves that pleiotropy is present; they indicate evidence from
the MR-Egger intercept test only.
"""
from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUTS = [
    ("cis_survival", ROOT / "cis/mr_results_cis/MR_Egger_pleiotropy.csv"),
    ("cis_incidence", ROOT / "cis/mr_results_cis_incidence/MR_Egger_pleiotropy.csv"),
    ("cis_slope_hunter", ROOT / "cis/mr_results_cis_sh/MR_Egger_pleiotropy.csv"),
    ("decode_survival", ROOT / "replicate/decode/mr_F10_survival/MR_Egger_pleiotropy.csv"),
    ("decode_incidence", ROOT / "replicate/decode/mr_F10_incidence/MR_Egger_pleiotropy.csv"),
    ("decode_slope_hunter", ROOT / "replicate/decode/mr_F10_sh/MR_Egger_pleiotropy.csv"),
    ("fenland_survival", ROOT / "replicate/fenland/mr_F10_survival/MR_Egger_pleiotropy.csv"),
    ("fenland_incidence", ROOT / "replicate/fenland/mr_F10_incidence/MR_Egger_pleiotropy.csv"),
    ("fenland_slope_hunter", ROOT / "replicate/fenland/mr_F10_sh/MR_Egger_pleiotropy.csv"),
]
OUTPUT_DIR = ROOT / "sensitivity"


def bh_qvalues(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted P values in the original row order."""
    order = sorted(range(len(pvalues)), key=pvalues.__getitem__)
    adjusted = [0.0] * len(pvalues)
    running = 1.0
    n = len(pvalues)
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, pvalues[index] * n / rank)
        adjusted[index] = running
    return adjusted


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    rows: list[dict[str, str]] = []
    for analysis, path in INPUTS:
        if not path.exists():
            raise FileNotFoundError(path)
        with path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                try:
                    pvalue = float(row["pval"])
                except (KeyError, TypeError, ValueError):
                    continue
                row["analysis_set"] = analysis
                row["pval"] = f"{pvalue:.12g}"
                rows.append(row)
    qvalues = bh_qvalues([float(row["pval"]) for row in rows])
    for row, qvalue in zip(rows, qvalues):
        pvalue = float(row["pval"])
        row["qval_bh_all_tests"] = f"{qvalue:.12g}"
        row["nominal_p_lt_0_05"] = str(pvalue < 0.05).lower()
        row["fdr_q_lt_0_05"] = str(qvalue < 0.05).lower()
    fields = ["analysis_set", "exposure_id_file", "outcome_id_file", "nsnp_retained",
              "egger_intercept", "se", "pval", "qval_bh_all_tests",
              "nominal_p_lt_0_05", "fdr_q_lt_0_05", "filename",
              "id.exposure", "id.outcome", "exposure", "outcome"]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT_DIR / "MR_Egger_pleiotropy_all_tested.csv", rows, fields)
    nominal = [row for row in rows if row["nominal_p_lt_0_05"] == "true"]
    fdr = [row for row in rows if row["fdr_q_lt_0_05"] == "true"]
    write_csv(OUTPUT_DIR / "MR_Egger_pleiotropy_nominal_p_lt_0_05.csv", nominal, fields)
    write_csv(OUTPUT_DIR / "MR_Egger_pleiotropy_FDR_q_lt_0_05.csv", fdr, fields)
    print(f"Successful MR-Egger tests: {len(rows):,}")
    print(f"Nominal P < 0.05: {len(nominal):,}")
    print(f"BH FDR q < 0.05: {len(fdr):,}")


if __name__ == "__main__":
    main()
