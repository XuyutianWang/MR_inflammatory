#!/usr/bin/env python3

import argparse
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = PROJECT_ROOT / "MR" / "exposure_all" / "clump_results_all"
DEFAULT_OUTPUT_FILE = PROJECT_ROOT / "MR" / "exposure_all" / "biomarkers_without_clumped.txt"


def find_biomarkers_without_clumped(input_dir: Path) -> list[str]:
	result_stems = {
		file_path.stem
		for file_path in input_dir.iterdir()
		if file_path.is_file() and file_path.suffix in {".log", ".nosex"}
	}
	clumped_stems = {
		file_path.stem
		for file_path in input_dir.glob("*.clumped")
		if file_path.is_file()
	}

	return sorted(
		stem.split("_", 1)[0]
		for stem in result_stems - clumped_stems
	)


def main() -> None:
	parser = argparse.ArgumentParser(
		description="Find biomarkers without a matching .clumped file."
	)
	parser.add_argument(
		"input_dir",
		nargs="?",
		type=Path,
		default=DEFAULT_INPUT_DIR,
		help="Directory containing clump result files.",
	)
	parser.add_argument(
		"-o",
		"--output",
		type=Path,
		default=DEFAULT_OUTPUT_FILE,
		help="File to write the biomarker list to.",
	)
	args = parser.parse_args()

	if not args.input_dir.is_dir():
		parser.error(f"Input directory does not exist: {args.input_dir}")

	biomarkers = find_biomarkers_without_clumped(args.input_dir)
	args.output.parent.mkdir(parents=True, exist_ok=True)
	args.output.write_text("\n".join(biomarkers) + ("\n" if biomarkers else ""))

	print(f"Found {len(biomarkers)} biomarkers without .clumped files:")
	print("\n".join(biomarkers))
	print(f"\nList written to: {args.output}")


if __name__ == "__main__":
	main()
