#!/usr/bin/env python3
"""Concatenate UKB-PPP chromosome 1-22 files into one GWAS file per protein."""

import argparse
import gzip
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional


CHROMOSOMES = range(1, 23)
CHR_FILE_RE = re.compile(r"^discovery_chr(\d+)_.+\.gz$")


def chromosome_files(protein_dir: Path) -> List[Path]:
    """Return chromosome 1-22 files in numeric order, rejecting bad inputs."""
    by_chromosome = {}  # type: Dict[int, List[Path]]
    for chromosome in CHROMOSOMES:
        by_chromosome[chromosome] = []

    for path in protein_dir.iterdir():
        if not path.is_file():
            continue
        match = CHR_FILE_RE.match(path.name)
        if match:
            chromosome = int(match.group(1))
            if chromosome in by_chromosome:
                by_chromosome[chromosome].append(path)

    missing = [str(c) for c, paths in by_chromosome.items() if not paths]
    duplicated = [str(c) for c, paths in by_chromosome.items() if len(paths) > 1]
    if missing or duplicated:
        problems = []
        if missing:
            problems.append(f"missing chromosome(s): {', '.join(missing)}")
        if duplicated:
            problems.append(f"duplicate chromosome(s): {', '.join(duplicated)}")
        raise ValueError("; ".join(problems))

    return [by_chromosome[c][0] for c in CHROMOSOMES]


def concatenate_protein(protein_dir: Path, output_file: Path, overwrite: bool) -> None:
    """Stream chromosome files into a gzip file, retaining only the first header."""
    if output_file.exists() and not overwrite:
        print(f"Skipping existing file: {output_file}")
        return

    files = chromosome_files(protein_dir)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None  # type: Optional[str]

    try:
        with tempfile.NamedTemporaryFile(
            dir=output_file.parent,
            prefix=f".{output_file.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_name = temporary.name
            with gzip.GzipFile(fileobj=temporary, mode="wb") as destination:
                expected_header: bytes | None = None
                for chromosome, source_file in zip(CHROMOSOMES, files):
                    with gzip.open(source_file, "rb") as source:
                        header = source.readline()
                        if not header:
                            raise ValueError(f"empty chromosome file: {source_file}")
                        if expected_header is None:
                            expected_header = header
                            destination.write(header)
                        elif header != expected_header:
                            raise ValueError(
                                f"header in chromosome {chromosome} differs: {source_file}"
                            )
                        shutil.copyfileobj(source, destination, length=1024 * 1024)

        os.replace(temporary_name, output_file)
        temporary_name = None
        print(f"Created: {output_file}")
    finally:
        if temporary_name is not None:
            try:
                Path(temporary_name).unlink()
            except FileNotFoundError:
                pass


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=script_dir.parent / "exposure_all" / "UKB_all",
        help="Directory containing one subdirectory per protein (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir.parent / "exposure_all" / "UKB_1_22_all",
        help="Destination directory (default: %(default)s)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace output files that already exist",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.is_dir():
        print(f"Input directory does not exist: {input_dir}", file=sys.stderr)
        return 1

    protein_dirs = sorted(path for path in input_dir.iterdir() if path.is_dir())
    if not protein_dirs:
        print(f"No protein directories found in: {input_dir}", file=sys.stderr)
        return 1

    failures = 0
    for protein_dir in protein_dirs:
        output_file = output_dir / f"{protein_dir.name}.gz"
        try:
            concatenate_protein(protein_dir, output_file, args.overwrite)
        except (OSError, EOFError, ValueError) as error:
            failures += 1
            print(f"ERROR [{protein_dir.name}]: {error}", file=sys.stderr)

    if failures:
        print(f"Finished with {failures} failed protein(s).", file=sys.stderr)
        return 1

    print(f"Finished: {len(protein_dirs)} protein file(s) written to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
