"""
Analyze the 5 Å heavy-atom interface contact-count distribution.

Input:
    results/pair_level_contacts.csv

Outputs:
    results/contact_distribution_5A/
        summary_5A.csv
        small_interface_candidates_5A.csv
        contact_distribution_all_5A.png
        contact_distribution_nonzero_5A.png
        contact_vs_iptm_5A.png
        contact_distribution_qc_5A.png

The pair-level table has already been hierarchically aggregated:
    archive reruns -> job level -> protein-pair level.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INPUT = Path("results/pair_level_contacts.csv")
OUTPUT_DIR = Path("results/contact_distribution_5A")

# Definition used during parsing
CONTACT_CUTOFF = 5.0

# Severe steric overlap threshold used in previous QC
SEVERE_OVERLAP_CUTOFF = 0.5


def load_data():
    """Load pair-level contact data."""
    if not INPUT.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT}")

    df = pd.read_csv(INPUT)

    required = [
        "protein_pair",
        "contact_count",
        "mean_pair_iptm",
        "min_interchain_distance",
        "has_severe_overlap",
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    return df


def make_summary(df):
    """Create summary statistics for the contact distribution."""

    nonzero = df[df["contact_count"] > 0]

    qc = df[~df["has_severe_overlap"].astype(bool)]

    rows = []

    for name, x in [
        ("all_pairs", df),
        ("nonzero_pairs", nonzero),
        ("exclude_severe_overlap", qc),
    ]:
        contacts = x["contact_count"]

        rows.append(
            {
                "dataset": name,
                "n_pairs": len(x),
                "zero_contacts_n": int(
                    (contacts == 0).sum()
                ),
                "zero_contacts_pct": (
                    100 * (contacts == 0).mean()
                ),
                "nonzero_contacts_n": int(
                    (contacts > 0).sum()
                ),
                "nonzero_contacts_pct": (
                    100 * (contacts > 0).mean()
                ),
                "mean": contacts.mean(),
                "median": contacts.median(),
                "p75": contacts.quantile(0.75),
                "p90": contacts.quantile(0.90),
                "p95": contacts.quantile(0.95),
                "p99": contacts.quantile(0.99),
                "max": contacts.max(),
            }
        )

    return pd.DataFrame(rows)


def plot_all_pairs(df):
    """
    Plot all pairs.

    Limit the visible x-axis to the 99th percentile so that
    extreme outliers do not compress the main distribution.
    """

    contacts = df["contact_count"]

    upper = contacts.quantile(0.99)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(
        contacts[contacts <= upper],
        bins=60,
    )

    ax.set_xlabel(
        "Mean inter-chain heavy-atom contact count (5 Å cutoff)"
    )
    ax.set_ylabel("Number of protein pairs")

    ax.set_title(
        "Interface-size distribution: all protein pairs"
    )

    ax.text(
        0.97,
        0.95,
        (
            f"N = {len(df):,}\n"
            f"Zero contacts = "
            f"{100 * (contacts == 0).mean():.1f}%\n"
            f"Display ≤ P99 ({upper:.1f})"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / "contact_distribution_all_5A.png",
        dpi=300,
    )

    plt.close(fig)


def plot_nonzero(df):
    """Plot only protein pairs with non-zero mean contacts."""

    x = df.loc[
        df["contact_count"] > 0,
        "contact_count",
    ]

    upper = x.quantile(0.99)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(
        x[x <= upper],
        bins=60,
    )

    ax.set_xlabel(
        "Mean inter-chain heavy-atom contact count (5 Å cutoff)"
    )
    ax.set_ylabel("Number of protein pairs")

    ax.set_title(
        "Interface-size distribution: non-zero interfaces"
    )

    ax.text(
        0.97,
        0.95,
        (
            f"N = {len(x):,}\n"
            f"Median = {x.median():.1f}\n"
            f"P95 = {x.quantile(0.95):.1f}\n"
            f"Display ≤ P99 ({upper:.1f})"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR /
        "contact_distribution_nonzero_5A.png",
        dpi=300,
    )

    plt.close(fig)


def plot_contact_vs_iptm(df):
    """Plot interface contact count against mean pair ipTM."""

    # Keep the plot readable by displaying contacts through P99.
    upper = df["contact_count"].quantile(0.99)

    x = df[df["contact_count"] <= upper].copy()

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(
        x["mean_pair_iptm"],
        x["contact_count"],
        s=6,
        alpha=0.15,
    )

    ax.set_xlabel("Mean pair ipTM")

    ax.set_ylabel(
        "Mean inter-chain heavy-atom contact count"
    )

    ax.set_title(
        "Interface size versus AlphaFold3 confidence"
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / "contact_vs_iptm_5A.png",
        dpi=300,
    )

    plt.close(fig)


def plot_qc_comparison(df):
    """
    Compare the main distribution before and after removing
    pairs containing at least one severe-overlap observation.
    """

    raw = df["contact_count"]

    qc = df.loc[
        ~df["has_severe_overlap"].astype(bool),
        "contact_count",
    ]

    upper = raw.quantile(0.99)

    bins = np.linspace(0, upper, 61)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(
        raw[raw <= upper],
        bins=bins,
        density=True,
        histtype="step",
        linewidth=2,
        label="All pairs",
    )

    ax.hist(
        qc[qc <= upper],
        bins=bins,
        density=True,
        histtype="step",
        linewidth=2,
        label="Exclude severe overlap",
    )

    ax.set_xlabel(
        "Mean inter-chain heavy-atom contact count (5 Å cutoff)"
    )

    ax.set_ylabel("Density")

    ax.set_title(
        "Sensitivity to severe structural overlaps"
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        OUTPUT_DIR / "contact_distribution_qc_5A.png",
        dpi=300,
    )

    plt.close(fig)


def export_small_interfaces(df):
    """
    Export candidate small interfaces for later structural inspection.

    These are real protein pairs whose pair-level mean contact count
    is >0 and <=10.

    Severe-overlap cases are excluded here.
    """

    small = df[
        (df["contact_count"] > 0)
        & (df["contact_count"] <= 10)
        & (~df["has_severe_overlap"].astype(bool))
    ].copy()

    small = small.sort_values(
        [
            "contact_count",
            "mean_pair_iptm",
        ],
        ascending=[True, False],
    )

    preferred_columns = [
        "protein_pair",
        "protein_a",
        "protein_b",
        "annotation_a",
        "annotation_b",
        "length_a",
        "length_b",
        "contact_count",
        "interface_residues_a",
        "interface_residues_b",
        "min_interchain_distance",
        "mean_pair_iptm",
        "n_jobs",
        "n_raw_observations",
        "contact_count_range",
    ]

    columns = [
        c for c in preferred_columns
        if c in small.columns
    ]

    small[columns].to_csv(
        OUTPUT_DIR /
        "small_interface_candidates_5A.csv",
        index=False,
    )

    return small


def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    df = load_data()

    print("=" * 70)
    print("5 Å CONTACT DISTRIBUTION ANALYSIS")
    print("=" * 70)

    print("Protein pairs:", len(df))
    print(
        "Unique protein pairs:",
        df["protein_pair"].nunique(),
    )

    print()

    summary = make_summary(df)

    summary.to_csv(
        OUTPUT_DIR / "summary_5A.csv",
        index=False,
    )

    print("Summary:")
    print(summary.to_string(index=False))

    print()

    small = export_small_interfaces(df)

    print(
        "Small-interface candidates "
        "(0 < mean contacts <= 10, no severe overlap):",
        len(small),
    )

    print()

    plot_all_pairs(df)
    plot_nonzero(df)
    plot_contact_vs_iptm(df)
    plot_qc_comparison(df)

    print("Saved outputs:")
    for path in sorted(OUTPUT_DIR.iterdir()):
        print(" ", path)


if __name__ == "__main__":
    main()