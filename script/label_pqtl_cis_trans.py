#!/usr/bin/env python3
"""Create a protein-to-gene map and label UKB-PPP lead pQTLs as cis or trans.

The input variant positions must be GRCh38. Cis is defined as a variant within
1 Mb of the strand-aware transcription start site (TSS) of the protein's
encoding gene. Original lead-SNP files are never modified.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path


ENSEMBL_URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"

# The source UniProt table contains several protein aliases rather than current
# HGNC symbols. These accession-specific overrides are deliberate and auditable.
# EBI3_IL27 is an IL-27 heterodimer, so cis includes either encoding gene.
OFFICIAL_GENE_OVERRIDES = {
    "P04141": ("CSF2",), "P09919": ("CSF3",), "P10145": ("CXCL8",),
    "Q16552": ("IL17A",), "P01374": ("LTA",), "P01137": ("TGFB1",),
    "P15692": ("VEGFA",), "P29965": ("CD40LG",), "Q8IV48": ("ERI1",),
    "P48023": ("FASLG",), "Q8IU54": ("IFNL1",), "P14778": ("IL1R1",),
    "P02776": ("PF4",), "P01009": ("SERPINA1",), "P05121": ("SERPINE1",),
    "P00995": ("SPINK1",), "P10451": ("SPP1",), "P50591": ("TNFSF10",),
    "O14788": ("TNFSF11",), "Q9Y275": ("TNFSF13B",), "O43557": ("TNFSF14",),
}
EXPOSURE_GENE_OVERRIDES = {
    "EBI3_IL27_Q14213_Q8NEV9_OID21389_v1_Oncology": ("IL27", "EBI3"),
}


def exposure_id_from_path(path: Path) -> str:
    suffix = "_lead_snps.txt"
    if not path.name.endswith(suffix):
        raise ValueError(f"Unexpected exposure filename: {path.name}")
    return path.name.removesuffix(suffix)


def uniprot_id_from_exposure(exposure_id: str) -> str:
    # UKB-PPP filenames begin: PROTEINNAME_UNIPROT_OID...
    fields = exposure_id.split("_")
    if len(fields) < 3 or not fields[1]:
        raise ValueError(f"Cannot extract UniProt accession from {exposure_id}")
    return fields[1]


def read_uniprot_map(path: Path) -> dict[str, list[dict[str, str]]]:
    by_accession: dict[str, list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            accession = (row.get("uniprot_id") or "").strip()
            gene = (row.get("matched_gene") or "").strip()
            if accession and gene:
                by_accession.setdefault(accession, []).append(row)
    return by_accession


def lookup_genes(symbols: list[str], retries: int = 4) -> dict[str, dict]:
    """Use Ensembl's batch endpoint to avoid one network request per protein."""
    request = urllib.request.Request(
        ENSEMBL_URL,
        data=json.dumps({"symbols": symbols, "expand": 0}).encode("utf-8"),
        headers={"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "cis-trans-pqtl-labeler/1.0"},
        method="POST",
    )
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return json.load(response)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
            last_error = error
            if isinstance(error, urllib.error.HTTPError) and error.code == 404:
                break
            time.sleep(1 + attempt)
    raise RuntimeError(f"Ensembl batch lookup failed: {last_error}")


def normalise_chromosome(value: str) -> str:
    return str(value).strip().removeprefix("chr").removeprefix("CHR")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("MR/exposure_all/lead_snps_all"))
    parser.add_argument("--uniprot-map", type=Path, default=Path("MR/exposure_all/download_UKBPPP/inflammatory_biomarkers_uniprot_ids.csv"))
    parser.add_argument("--map-output", type=Path, default=Path("MR/exposure_all/protein_gene_map.csv"))
    parser.add_argument("--labelled-dir", type=Path, default=Path("MR/exposure_all/lead_snps_cis_trans"))
    parser.add_argument("--audit-output", type=Path, default=Path("MR/exposure_all/protein_gene_map_audit.csv"))
    parser.add_argument("--cis-window", type=int, default=1_000_000)
    args = parser.parse_args()

    input_files = sorted(args.input_dir.glob("*_lead_snps.txt"))
    if not input_files:
        raise SystemExit(f"No *_lead_snps.txt files found in {args.input_dir}")
    by_accession = read_uniprot_map(args.uniprot_map)
    args.labelled_dir.mkdir(parents=True, exist_ok=True)

    exposure_symbols: dict[Path, tuple[str, str, tuple[str, ...]]] = {}
    audit_rows: list[dict[str, str]] = []
    for input_file in input_files:
        exposure_id = exposure_id_from_path(input_file)
        accession = uniprot_id_from_exposure(exposure_id)
        if exposure_id in EXPOSURE_GENE_OVERRIDES:
            exposure_symbols[input_file] = (exposure_id, accession, EXPOSURE_GENE_OVERRIDES[exposure_id])
            continue
        if accession in OFFICIAL_GENE_OVERRIDES:
            exposure_symbols[input_file] = (exposure_id, accession, OFFICIAL_GENE_OVERRIDES[accession])
            continue
        candidates = by_accession.get(accession, [])
        candidate_symbols = sorted({(x.get("matched_gene") or "").strip() for x in candidates if (x.get("matched_gene") or "").strip()})
        if len(candidate_symbols) != 1:
            audit_rows.append({"exposure_id": exposure_id, "uniprot_id": accession, "status": "unresolved_or_ambiguous_uniprot_mapping", "detail": "|".join(candidate_symbols)})
        else:
            exposure_symbols[input_file] = (exposure_id, accession, (candidate_symbols[0],))

    ensembl_records = lookup_genes(sorted({symbol for item in exposure_symbols.values() for symbol in item[2]}))
    resolved: dict[str, dict[str, str]] = {}
    for symbol, record in ensembl_records.items():
        if record and record.get("assembly_name") == "GRCh38" and record.get("seq_region_name"):
            resolved[symbol] = {
                "chromosome": normalise_chromosome(record["seq_region_name"]),
                "tss": str(record["start"] if record["strand"] == 1 else record["end"]),
            }

    map_rows: list[dict[str, str]] = []

    for input_file in input_files:
        if input_file not in exposure_symbols:
            continue
        exposure_id, accession, symbols = exposure_symbols[input_file]
        if any(symbol not in resolved for symbol in symbols):
            audit_rows.append({"exposure_id": exposure_id, "uniprot_id": accession, "status": "ensembl_lookup_failed", "detail": "|".join(symbols)})
            continue

        coords = [resolved[symbol] for symbol in symbols]
        map_rows.append({
            "exposure_id": exposure_id,
            "gene_symbol": "|".join(symbols),
            "chromosome": "|".join(coord["chromosome"] for coord in coords),
            "tss": "|".join(coord["tss"] for coord in coords),
        })

        with input_file.open(newline="", encoding="utf-8") as source:
            reader = csv.DictReader(source, delimiter="\t")
            required = {"CHROM", "GENPOS"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise RuntimeError(f"{input_file} must contain tab-delimited CHROM and GENPOS columns")
            fieldnames = list(reader.fieldnames)
            if "qtl_type" not in fieldnames:
                fieldnames.append("qtl_type")
            rows = []
            for row in reader:
                chrom = normalise_chromosome(row["CHROM"])
                try:
                    position = int(float(row["GENPOS"]))
                except ValueError:
                    row["qtl_type"] = "unclassified_invalid_position"
                else:
                    row["qtl_type"] = "cis" if any(
                        chrom == coord["chromosome"] and abs(position - int(coord["tss"])) <= args.cis_window
                        for coord in coords
                    ) else "trans"
                rows.append(row)

        output_file = args.labelled_dir / input_file.name
        with output_file.open("w", newline="", encoding="utf-8") as target:
            writer = csv.DictWriter(target, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    for output, rows, columns in [
        (args.map_output, map_rows, ["exposure_id", "gene_symbol", "chromosome", "tss"]),
        (args.audit_output, audit_rows, ["exposure_id", "uniprot_id", "status", "detail"]),
    ]:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)

    print(f"Mapped and labelled: {len(map_rows)} exposures")
    print(f"Requires review: {len(audit_rows)} exposures; see {args.audit_output}")
    print(f"Protein/gene map: {args.map_output}")
    print(f"Labelled files: {args.labelled_dir}")


if __name__ == "__main__":
    main()
