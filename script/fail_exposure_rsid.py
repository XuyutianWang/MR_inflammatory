#!/usr/bin/env python3

import argparse
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_DIR = PROJECT_ROOT / "MR" / "exposure_all" / "UKB_1_22_all_pfiltered_5108"
DEFAULT_RSID_DIR = PROJECT_ROOT / "MR" / "exposure_all" / "UKB_1_22_all_rsid"


def gz_filenames(directory: Path) -> set[str]:
	return {
		file_path.name
		for file_path in directory.iterdir()
		if file_path.is_file() and file_path.suffix == ".gz"
	}


def find_failed_files(source_dir: Path, rsid_dir: Path) -> list[str]:
	return sorted(gz_filenames(source_dir) - gz_filenames(rsid_dir))


def main() -> None:
	parser = argparse.ArgumentParser(
		description="Find exposure files that do not have a matching rsid file."
	)
	parser.add_argument(
		"source_dir",
		nargs="?",
		type=Path,
		default=DEFAULT_SOURCE_DIR,
		help="Directory containing the original filtered exposure files.",
	)
	parser.add_argument(
		"rsid_dir",
		nargs="?",
		type=Path,
		default=DEFAULT_RSID_DIR,
		help="Directory containing the rsid-converted exposure files.",
	)
	args = parser.parse_args()

	for directory in (args.source_dir, args.rsid_dir):
		if not directory.is_dir():
			parser.error(f"Directory does not exist: {directory}")

	failed_files = find_failed_files(args.source_dir, args.rsid_dir)
	print(f"Found {len(failed_files)} files missing from the rsid directory:")
	print("\n".join(failed_files))


if __name__ == "__main__":
	main()
