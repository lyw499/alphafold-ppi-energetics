"""Export every unique protein pair from Pedro Dataset EV4.

No STRING or AlphaFold confidence filtering is applied.
Self-pairs and symmetric duplicates are excluded.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def export_all_pairs(workbook: Path) -> pd.DataFrame:
    """Convert the EV4 raw-ipTM matrix to one row per unique protein pair."""

    matrix = pd.read_excel(
        workbook,
        sheet_name="raw ipTM",
        index_col=0,
    )

    matrix.index = matrix.index.astype(str).str.strip()
    matrix.columns = matrix.columns.astype(str).str.strip()

    if set(matrix.index) != set(matrix.columns):
        raise ValueError("Protein IDs in rows and columns do not match.")

    proteins = list(matrix.index)
    matrix = matrix.loc[proteins, proteins]

    # k=1 means use only the upper triangle:
    # no self-pairs and no symmetric duplicates.
    row_indices, column_indices = np.triu_indices(
        len(proteins),
        k=1,
    )

    pairs = pd.DataFrame(
        {
            "protein_1": np.asarray(proteins)[row_indices],
            "protein_2": np.asarray(proteins)[column_indices],
            "raw_iptm": matrix.to_numpy()[row_indices, column_indices],
        }
    )

    return pairs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")

    pairs = export_all_pairs(args.input)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.output, index=False)

    print(f"Proteins: {len(set(pairs['protein_1']) | set(pairs['protein_2'])):,}")
    print(f"All unique pairs: {len(pairs):,}")
    print("No STRING or confidence filtering was applied.")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()