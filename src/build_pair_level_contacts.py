"""
Build a pair-level protein interface table from parsed pooled AlphaFold3 outputs.

Aggregation is hierarchical:

1. Multiple archived outputs of the same job + protein pair are treated as
   reruns of the same pooled AF3 job and are averaged first.

2. If a protein pair appears in multiple different pooled jobs, the resulting
   job-level measurements are averaged across jobs.

This prevents duplicated archived outputs of one job from giving that job
extra weight.

The final output contains one row per M. genitalium protein pair.
"""

from pathlib import Path

import pandas as pd


INPUT_DIR = Path("results/pooled_parts")
OUTPUT_DIR = Path("results")

JOB_LEVEL_OUTPUT = OUTPUT_DIR / "job_level_contacts.csv"
PAIR_LEVEL_OUTPUT = OUTPUT_DIR / "pair_level_contacts.csv"

SEVERE_OVERLAP_CUTOFF = 0.5


def load_observations() -> pd.DataFrame:
    """Load all parsed pooled-AF3 CSV files."""

    files = sorted(INPUT_DIR.glob("*.csv"))

    if not files:
        raise FileNotFoundError(
            f"No CSV files found in {INPUT_DIR}"
        )

    frames = []

    for path in files:
        df = pd.read_csv(path)
        df["source_file"] = path.name
        frames.append(df)

    observations = pd.concat(frames, ignore_index=True)

    observations["is_nonstandard_job"] = (
        observations["job_id"]
        .astype(str)
        .str.endswith("_2")
    )

    observations["severe_overlap"] = (
        observations["min_interchain_distance"]
        < SEVERE_OVERLAP_CUTOFF
    )

    return observations


def build_job_level(observations: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse repeated archived outputs of the same job.

    One row = one job_id + protein_pair.
    """

    job_level = (
        observations
        .groupby(
            ["job_id", "protein_pair"],
            as_index=False
        )
        .agg(
            protein_a=("protein_a", "first"),
            protein_b=("protein_b", "first"),

            annotation_a=("annotation_a", "first"),
            annotation_b=("annotation_b", "first"),

            length_a=("length_a", "first"),
            length_b=("length_b", "first"),

            # Number of archived versions of this same job/pair
            n_reruns=("source_file", "nunique"),

            # Structural measurements averaged across reruns
            contact_count_mean=("contact_count", "mean"),
            interface_residues_a_mean=(
                "interface_residues_a", "mean"
            ),
            interface_residues_b_mean=(
                "interface_residues_b", "mean"
            ),

            # Keep structural variability as well
            contact_count_min=("contact_count", "min"),
            contact_count_max=("contact_count", "max"),

            min_interchain_distance=(
                "min_interchain_distance", "min"
            ),

            # AF3 confidence
            mean_pair_iptm=("mean_pair_iptm", "mean"),
            job_ranking_score=("job_ranking_score", "mean"),

            # QC
            n_severe_overlap=("severe_overlap", "sum"),
            any_nonstandard_job=("is_nonstandard_job", "max"),

            source_files=(
                "source_file",
                lambda x: " | ".join(sorted(set(x)))
            ),
        )
    )

    job_level["contact_count_range"] = (
        job_level["contact_count_max"]
        - job_level["contact_count_min"]
    )

    return job_level


def build_pair_level(job_level: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse different pooled jobs for the same protein pair.

    Each job receives equal weight.

    One row = one protein pair.
    """

    pair_level = (
        job_level
        .groupby("protein_pair", as_index=False)
        .agg(
            protein_a=("protein_a", "first"),
            protein_b=("protein_b", "first"),

            annotation_a=("annotation_a", "first"),
            annotation_b=("annotation_b", "first"),

            length_a=("length_a", "first"),
            length_b=("length_b", "first"),

            # Number of independent pooled jobs
            n_jobs=("job_id", "nunique"),

            # Total number of archived structural outputs
            n_raw_observations=("n_reruns", "sum"),

            # Main interface-size measurements
            contact_count=(
                "contact_count_mean", "mean"
            ),

            interface_residues_a=(
                "interface_residues_a_mean", "mean"
            ),

            interface_residues_b=(
                "interface_residues_b_mean", "mean"
            ),

            # Between-job variability
            contact_count_min=(
                "contact_count_mean", "min"
            ),

            contact_count_max=(
                "contact_count_mean", "max"
            ),

            # Closest inter-chain approach observed
            min_interchain_distance=(
                "min_interchain_distance", "min"
            ),

            # Confidence averaged across jobs
            mean_pair_iptm=(
                "mean_pair_iptm", "mean"
            ),

            # QC
            n_severe_overlap=(
                "n_severe_overlap", "sum"
            ),

            any_nonstandard_job=(
                "any_nonstandard_job", "max"
            ),

            job_ids=(
                "job_id",
                lambda x: " | ".join(sorted(set(x)))
            ),
        )
    )

    pair_level["contact_count_range"] = (
        pair_level["contact_count_max"]
        - pair_level["contact_count_min"]
    )

    pair_level["has_severe_overlap"] = (
        pair_level["n_severe_overlap"] > 0
    )

    return pair_level


def print_summary(
    observations: pd.DataFrame,
    job_level: pd.DataFrame,
    pair_level: pd.DataFrame,
) -> None:
    """Print basic integrity and distribution checks."""

    print("=" * 70)
    print("RAW OBSERVATIONS")
    print("=" * 70)

    print("Raw observations:", len(observations))
    print(
        "Unique protein pairs:",
        observations["protein_pair"].nunique(),
    )
    print(
        "Unique jobs:",
        observations["job_id"].nunique(),
    )

    print(
        "Non-standard _2 observations:",
        observations["is_nonstandard_job"].sum(),
    )

    print()

    print("=" * 70)
    print("JOB LEVEL")
    print("=" * 70)

    print("Job + pair rows:", len(job_level))

    print(
        "Job+pair entries with >1 archived rerun:",
        (job_level["n_reruns"] > 1).sum(),
    )

    print()

    print("=" * 70)
    print("PAIR LEVEL")
    print("=" * 70)

    print("Final protein pairs:", len(pair_level))

    print(
        "Pairs appearing in >1 pooled job:",
        (pair_level["n_jobs"] > 1).sum(),
    )

    print(
        "Pairs with any severe overlap:",
        pair_level["has_severe_overlap"].sum(),
    )

    print(
        "Pairs involving a non-standard _2 job:",
        pair_level["any_nonstandard_job"].sum(),
    )

    print()

    print("Contact-count distribution:")
    print(
        pair_level["contact_count"].describe(
            percentiles=[0.5, 0.75, 0.9, 0.95, 0.99]
        )
    )

    print()

    print(
        "Pairs with mean contact count > 0:",
        (pair_level["contact_count"] > 0).sum(),
    )

    print(
        "Percentage with mean contact count > 0:",
        100 * (pair_level["contact_count"] > 0).mean(),
    )


def main() -> None:

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    observations = load_observations()

    job_level = build_job_level(observations)

    pair_level = build_pair_level(job_level)

    job_level.to_csv(
        JOB_LEVEL_OUTPUT,
        index=False
    )

    pair_level.to_csv(
        PAIR_LEVEL_OUTPUT,
        index=False
    )

    print_summary(
        observations,
        job_level,
        pair_level
    )

    print()
    print("Saved:")
    print(JOB_LEVEL_OUTPUT)
    print(PAIR_LEVEL_OUTPUT)


if __name__ == "__main__":
    main()