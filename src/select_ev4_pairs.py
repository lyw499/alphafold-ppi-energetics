"""Select experimentally supported PPI candidates from Pedro Dataset EV4.

The EV4 workbook stores each score as a symmetric protein-by-protein matrix.
This script converts the matrices to a tidy table with one row per unique
protein pair. STRING experimental evidence is used for selection; AlphaFold3
scores are retained as annotations for later stratification and analysis.

Example
-------
python src/select_ev4_pairs.py \
    --input data/raw/pedro_dataset_EV4.xlsx \
    --output data/processed/string_supported_pairs.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SHEETS = {
    "string_experimental": "STRING_experimental channel",
    "raw_iptm": "raw ipTM",
    "size_corrected_iptm": "size_corrected ipTM",
    "corrected_iptm_correlation": "correlation of corrected ipTM",
    "iptm_ratio": "ipTM-ratio",
    "iptm_ratio_correlation": "correlation of ipTM-ratio",
}


def read_matrix(workbook: Path, sheet_name: str) -> pd.DataFrame:
    """Read one EV4 worksheet as a numeric matrix with protein IDs as labels."""
    matrix = pd.read_excel(workbook, sheet_name=sheet_name, index_col=0)
    matrix.index = matrix.index.astype(str).str.strip()
    matrix.columns = matrix.columns.astype(str).str.strip()
    matrix = matrix.apply(pd.to_numeric, errors="coerce")

    if matrix.index.has_duplicates or matrix.columns.has_duplicates:
        raise ValueError(f"Duplicate protein IDs found in sheet {sheet_name!r}.")
    return matrix


def validate_labels(matrices: dict[str, pd.DataFrame]) -> list[str]:
    """Confirm that all matrices contain the same proteins in both dimensions."""
    first_name, first = next(iter(matrices.items()))
    proteins = list(first.index)
    expected = set(proteins)

    if set(first.columns) != expected:
        raise ValueError(f"Rows and columns differ in sheet {SHEETS[first_name]!r}.")

    for name, matrix in matrices.items():
        if set(matrix.index) != expected or set(matrix.columns) != expected:
            raise ValueError(
                f"Protein labels in {SHEETS[name]!r} do not match "
                f"{SHEETS[first_name]!r}."
            )
    return proteins


def matrices_to_pairs(matrices: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Convert aligned symmetric matrices to one row per unordered pair."""
    proteins = validate_labels(matrices)
    aligned = {
        name: matrix.loc[proteins, proteins] for name, matrix in matrices.items()
    }

    row_indices, column_indices = np.triu_indices(len(proteins), k=1)
    pairs = pd.DataFrame(
        {
            "protein_1": np.asarray(proteins)[row_indices],
            "protein_2": np.asarray(proteins)[column_indices],
        }
    )

    for name, matrix in aligned.items():
        pairs[name] = matrix.to_numpy()[row_indices, column_indices]

    pairs["combined_score"] = (
        pairs["size_corrected_iptm"]
        + 0.2 * pairs["corrected_iptm_correlation"]
    )
    pairs["is_string_999"] = pairs["string_experimental"].eq(999)
    return pairs


def select_pairs(
    workbook: Path,
    string_threshold: float,
    include_equal: bool,
) -> pd.DataFrame:
    """Select pairs using only the STRING experimental-channel threshold."""
    matrices = {
        name: read_matrix(workbook, sheet_name)
        for name, sheet_name in SHEETS.items()
    }
    pairs = matrices_to_pairs(matrices)

    if include_equal:
        keep = pairs["string_experimental"] >= string_threshold
    else:
        keep = pairs["string_experimental"] > string_threshold

    selected = pairs.loc[keep].copy()
    return selected.sort_values(
        ["string_experimental", "combined_score", "size_corrected_iptm"],
        ascending=False,
        na_position="last",
    ).reset_index(drop=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="EV4 .xlsx file")
    parser.add_argument("--output", type=Path, required=True, help="Output .csv file")
    parser.add_argument(
        "--string-threshold",
        type=float,
        default=800,
        help="STRING experimental threshold (default: 800)",
    )
    parser.add_argument(
        "--include-equal",
        action="store_true",
        help="Use >= threshold instead of the paper's > threshold.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not args.input.exists():
        raise FileNotFoundError(f"Input workbook not found: {args.input}")
    if args.input.suffix.lower() != ".xlsx":
        raise ValueError("EV4 input must be an .xlsx workbook.")

    selected = select_pairs(
        args.input,
        string_threshold=args.string_threshold,
        include_equal=args.include_equal,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, index=False)

    operator = ">=" if args.include_equal else ">"
    strongest = int(selected["is_string_999"].sum())
    print(f"Selected {len(selected):,} unique protein pairs.")
    print(f"Criterion: STRING experimental {operator} {args.string_threshold:g}")
    print(f"Of these, {strongest:,} pairs have STRING experimental = 999.")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
