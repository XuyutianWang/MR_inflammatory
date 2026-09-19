#!/usr/bin/env python3

import argparse
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EXPOSURE_DIR = PROJECT_ROOT / "MR" / "exposure_all" / "lead_snps_all"
DEFAULT_OUTCOME_DIR = PROJECT_ROOT / "MR" / "incidence" / "rsid_no_chr23"
DEFAULT_HARMONISED_DIR = PROJECT_ROOT / "MR" / "harmonised_results_incidence"


def exposure_stems(exposure_dir: Path) -> set[str]:
    return {
        path.name.removesuffix("_lead_snps.txt")
        for path in exposure_dir.glob("*_lead_snps.txt")
        if path.is_file()
    }


def harmonised_pairs(harmonised_dir: Path) -> set[tuple[str, str]]:
    pairs = set()
    for path in harmonised_dir.glob("*_harmonised.csv"):
        if not path.is_file():
            continue

        name = path.name.removesuffix("_harmonised.csv")
        if "__VS__" not in name:
            continue

        exposure, outcome = name.split("__VS__", maxsplit=1)
        pairs.add((exposure, outcome))
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find exposure/outcome pairs missing harmonised results."
    )
    parser.add_argument(
        "--exposure-dir",
        type=Path,
        default=DEFAULT_EXPOSURE_DIR,
        help="Directory containing *_lead_snps.txt exposure files.",
    )
    parser.add_argument(
        "--harmonised-dir",
        type=Path,
        default=DEFAULT_HARMONISED_DIR,
        help="Directory containing *_harmonised.csv files.",
    )
    args = parser.parse_args()

    for directory in (args.exposure_dir, args.harmonised_dir):
        if not directory.is_dir():
            parser.error(f"Directory does not exist: {directory}")

    exposures = exposure_stems(args.exposure_dir)
    observed_pairs = harmonised_pairs(args.harmonised_dir)
    outcomes = {outcome for _, outcome in observed_pairs}
    expected_pairs = {(exposure, outcome) for exposure in exposures for outcome in outcomes}
    missing_pairs = sorted(expected_pairs - observed_pairs)
    missing_exposures = sorted(
        exposure for exposure in exposures if not any(pair[0] == exposure for pair in observed_pairs)
    )

    print(f"Exposures: {len(exposures)}")
    print(f"Outcomes: {len(outcomes)}")
    print(f"Expected harmonised results: {len(expected_pairs)}")
    print(f"Existing harmonised results: {len(observed_pairs)}")
    print(f"Missing harmonised results: {len(missing_pairs)}")
    print(f"Exposures with no harmonised results: {len(missing_exposures)}")

    if missing_exposures:
        print("\nExposures with no harmonised results:")
        print("\n".join(missing_exposures))

    if missing_pairs:
        print("\nMissing exposure/outcome pairs:")
        for exposure, outcome in missing_pairs:
            print(f"{exposure}__VS__{outcome}")


if __name__ == "__main__":
    main()