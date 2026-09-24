#!/usr/bin/env python3
"""Create a protein-to-gene map and label UKB-PPP lead pQTLs as cis or trans.

The input variant positions must be GRCh38. Cis is defined as a variant within
1 Mb of the strand-aware transcription start site (TSS) of the protein's
encoding gene. Original lead-SNP files are never modified.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import time
import urllib.error
import urllib.request
from pathlib import Path


# Ensembl's current GraphQL API and its GRCh38.p14 assembly accession.
ENSEMBL_GRAPHQL_URL = "https://www.ensembl.org/data/graphql"
GRCH38_ASSEMBLY_ACCESSION = "GCA_000001405.29"
GENOME_QUERY = """
query Genome($assembly: String!) {
  genomes(by_keyword: {assembly_accession_id: $assembly}) { genome_id }
}
"""
GENE_QUERY = """
query Gene($genome_id: String!, $symbol: String!) {
  genes(by_symbol: {genome_id: $genome_id, symbol: $symbol}) {
    slice {
      region { name }
      location { start end }
      strand { value }
    }
  }
}
"""

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
    """Return the UKB-PPP exposure ID from a raw or lead-SNP filename."""
    name = path.name.removesuffix(".gz")
    return name.removesuffix("_lead_snps.txt").removesuffix(".txt")


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


def lookup_genes(
    symbols: list[str],
    graphql_url: str = ENSEMBL_GRAPHQL_URL,
    retries: int = 4,
) -> dict[str, dict]:
    """Resolve symbols using Ensembl's GraphQL API on GRCh38.p14."""
    genome_id = lookup_grch38_genome_id(graphql_url, retries)
    records: dict[str, dict] = {}
    for symbol in symbols:
        record = lookup_gene(symbol, genome_id, graphql_url, retries)
        if record is not None:
            records[symbol] = record
    return records


def lookup_grch38_genome_id(graphql_url: str, retries: int) -> str:
    """Return the GraphQL genome ID for the current GRCh38 annotation."""
    data = graphql_request(
        GENOME_QUERY,
        {"assembly": GRCH38_ASSEMBLY_ACCESSION},
        graphql_url,
        retries,
    )
    genomes = data.get("genomes", [])
    if not genomes or not genomes[0].get("genome_id"):
        raise RuntimeError(f"Ensembl GraphQL did not return {GRCH38_ASSEMBLY_ACCESSION}")
    return genomes[0]["genome_id"]


def lookup_gene(
    symbol: str,
    genome_id: str,
    graphql_url: str,
    retries: int,
) -> dict | None:
    """Return one GRCh38 gene record, or None when the symbol is unknown."""
    data = graphql_request(
        GENE_QUERY,
        {"genome_id": genome_id, "symbol": symbol},
        graphql_url,
        retries,
    )
    genes = data.get("genes", [])
    if not genes:
        return None
    gene_slice = genes[0].get("slice") or {}
    location = gene_slice.get("location") or {}
    region = gene_slice.get("region") or {}
    strand = gene_slice.get("strand") or {}
    if not region.get("name") or not location.get("start") or not location.get("end"):
        return None
    return {
        "assembly_name": "GRCh38",
        "seq_region_name": str(region["name"]),
        "start": location["start"],
        "end": location["end"],
        "strand": strand.get("value"),
    }


def graphql_request(
    query: str,
    variables: dict[str, str],
    graphql_url: str,
    retries: int,
) -> dict:
    """Execute a GraphQL request and retry only transport failures."""
    request = urllib.request.Request(
        graphql_url,
        data=json.dumps({"query": query, "variables": variables}).encode("utf-8"),
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": "cis-trans-pqtl-labeler/2.0",
        },
        method="POST",
    )
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = json.load(response)
            if payload.get("errors"):
                messages = "; ".join(error.get("message", str(error)) for error in payload["errors"])
                raise RuntimeError(f"Ensembl GraphQL query failed: {messages}")
            if not isinstance(payload.get("data"), dict):
                raise RuntimeError("Ensembl GraphQL response did not contain a data object")
            return payload["data"]
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
            last_error = error
            if isinstance(error, urllib.error.HTTPError) and error.code not in {408, 429, 500, 502, 503, 504}:
                raise RuntimeError(f"Ensembl GraphQL request failed: {error}") from error
            retry_after = (
                (error.headers or {}).get("Retry-After")
                if isinstance(error, urllib.error.HTTPError)
                else None
            )
            try:
                delay = float(retry_after) if retry_after else 1 + attempt
            except ValueError:
                delay = 1 + attempt
            time.sleep(delay)
    raise RuntimeError(f"Ensembl GraphQL request failed after {retries} attempts: {last_error}")


def normalise_chromosome(value: str) -> str:
    return str(value).strip().removeprefix("chr").removeprefix("CHR")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("exposure_all/UKB_1_22_all_rsid"))
    parser.add_argument("--uniprot-map", type=Path, default=Path("exposure_all/download_UKBPPP/inflammatory_biomarkers_uniprot_ids.csv"))
    parser.add_argument("--map-output", type=Path, default=Path("exposure_all/protein_gene_map.csv"))
    parser.add_argument("--labelled-dir", type=Path, default=Path("exposure_all/UKB_1_22_all_rsid_cis_trans"))
    parser.add_argument("--audit-output", type=Path, default=Path("exposure_all/protein_gene_map_audit.csv"))
    parser.add_argument("--cis-window", type=int, default=1_000_000)
    parser.add_argument(
        "--ensembl-graphql-url",
        default=ENSEMBL_GRAPHQL_URL,
        help="Ensembl GraphQL endpoint URL",
    )
    args = parser.parse_args()

    input_files = sorted(args.input_dir.glob("*.gz"))
    if not input_files:
        raise SystemExit(f"No gzipped exposure files found in {args.input_dir}")
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

    ensembl_records = lookup_genes(
        sorted({symbol for item in exposure_symbols.values() for symbol in item[2]}),
        graphql_url=args.ensembl_graphql_url,
    )
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

        with gzip.open(input_file, mode="rt", newline="", encoding="utf-8") as source:
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

        output_file = args.labelled_dir / f"{exposure_id}.txt"
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
