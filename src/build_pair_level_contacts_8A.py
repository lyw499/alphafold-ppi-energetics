from pathlib import Path

import pandas as pd


INPUT_DIR = Path("results/pooled_parts_8A")
OUTPUT_DIR = Path("results")

JOB_LEVEL_OUTPUT = OUTPUT_DIR / "job_level_contacts_8A.csv"
PAIR_LEVEL_OUTPUT = OUTPUT_DIR / "pair_level_contacts_8A.csv"

SEVERE_OVERLAP_CUTOFF = 0.5


def load_observations():
    """Load all 8 Å pooled-AF3 parsed CSV files."""

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

    observations = pd.concat(
        frames,
        ignore_index=True,
    )

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


def join_unique_strings(values):
    """Join unique non-empty string values."""

    values = [
        str(x)
        for x in values
        if pd.notna(x)
        and str(x).strip()
    ]

    return " | ".join(
        sorted(set(values))
    )


def build_job_level(observations):
    """
    Collapse repeated archived outputs of the same
    job_id + protein_pair.

    One row = one pooled job observation for one protein pair.
    """

    job_level = (
        observations
        .groupby(
            ["job_id", "protein_pair"],
            as_index=False,
        )
        .agg(
            protein_a=("protein_a", "first"),
            protein_b=("protein_b", "first"),

            annotation_a=("annotation_a", "first"),
            annotation_b=("annotation_b", "first"),

            chain_a=("chain_a", "first"),
            chain_b=("chain_b", "first"),

            length_a=("length_a", "first"),
            length_b=("length_b", "first"),

            short_chain_id=("short_chain_id", "first"),
            long_chain_id=("long_chain_id", "first"),

            short_protein=("short_protein", "first"),
            long_protein=("long_protein", "first"),

            short_chain_length=(
                "short_chain_length",
                "first",
            ),

            long_chain_length=(
                "long_chain_length",
                "first",
            ),

            # Number of archived reruns of this job/pair
            n_reruns=("source_file", "nunique"),

            # Main interface measurements
            contact_count_mean=(
                "contact_count",
                "mean",
            ),

            contact_count_min=(
                "contact_count",
                "min",
            ),

            contact_count_max=(
                "contact_count",
                "max",
            ),

            interface_residues_a_mean=(
                "interface_residues_a",
                "mean",
            ),

            interface_residues_b_mean=(
                "interface_residues_b",
                "mean",
            ),

            interface_residues_short_mean=(
                "interface_residues_short",
                "mean",
            ),

            interface_residues_long_mean=(
                "interface_residues_long",
                "mean",
            ),

            fraction_interface_short_mean=(
                "fraction_interface_short",
                "mean",
            ),

            fraction_interface_long_mean=(
                "fraction_interface_long",
                "mean",
            ),

            contact_density_mean=(
                "contact_density",
                "mean",
            ),

            # Keep residue IDs from all archived runs
            interface_residue_ids_a=(
                "interface_residue_ids_a",
                join_unique_strings,
            ),

            interface_residue_ids_b=(
                "interface_residue_ids_b",
                join_unique_strings,
            ),

            interface_residue_ids_short=(
                "interface_residue_ids_short",
                join_unique_strings,
            ),

            interface_residue_ids_long=(
                "interface_residue_ids_long",
                join_unique_strings,
            ),

            # Geometry QC
            min_interchain_distance=(
                "min_interchain_distance",
                "min",
            ),

            # AF3 confidence
            mean_pair_iptm=(
                "mean_pair_iptm",
                "mean",
            ),

            job_ranking_score=(
                "job_ranking_score",
                "mean",
            ),

            # AF3 clash annotation
            any_af3_clash=(
                "has_clash",
                "max",
            ),

            # Our geometric overlap flag
            n_severe_overlap=(
                "severe_overlap",
                "sum",
            ),

            any_nonstandard_job=(
                "is_nonstandard_job",
                "max",
            ),

            source_files=(
                "source_file",
                join_unique_strings,
            ),
        )
    )

    job_level["contact_count_range"] = (
        job_level["contact_count_max"]
        - job_level["contact_count_min"]
    )

    job_level["has_severe_overlap"] = (
        job_level["n_severe_overlap"] > 0
    )

    return job_level


def build_pair_level(job_level):
    """
    Collapse different pooled jobs for the same protein pair.

    Each pooled job receives equal weight.

    One row = one protein pair.
    """

    pair_level = (
        job_level
        .groupby(
            "protein_pair",
            as_index=False,
        )
        .agg(
            protein_a=("protein_a", "first"),
            protein_b=("protein_b", "first"),

            annotation_a=("annotation_a", "first"),
            annotation_b=("annotation_b", "first"),

            length_a=("length_a", "first"),
            length_b=("length_b", "first"),

            short_protein=(
                "short_protein",
                "first",
            ),

            long_protein=(
                "long_protein",
                "first",
            ),

            short_chain_length=(
                "short_chain_length",
                "first",
            ),

            long_chain_length=(
                "long_chain_length",
                "first",
            ),

            # Number of pooled jobs
            n_jobs=(
                "job_id",
                "nunique",
            ),

            # Number of raw archived structural outputs
            n_raw_observations=(
                "n_reruns",
                "sum",
            ),

            # Main pair-level interface measurements
            contact_count=(
                "contact_count_mean",
                "mean",
            ),

            interface_residues_a=(
                "interface_residues_a_mean",
                "mean",
            ),

            interface_residues_b=(
                "interface_residues_b_mean",
                "mean",
            ),

            interface_residues_short=(
                "interface_residues_short_mean",
                "mean",
            ),

            interface_residues_long=(
                "interface_residues_long_mean",
                "mean",
            ),

            fraction_interface_short=(
                "fraction_interface_short_mean",
                "mean",
            ),

            fraction_interface_long=(
                "fraction_interface_long_mean",
                "mean",
            ),

            contact_density=(
                "contact_density_mean",
                "mean",
            ),

            # Variability across pooled jobs
            contact_count_min=(
                "contact_count_mean",
                "min",
            ),

            contact_count_max=(
                "contact_count_mean",
                "max",
            ),

            # Closest approach seen in any pooled job
            min_interchain_distance=(
                "min_interchain_distance",
                "min",
            ),

            # Confidence
            mean_pair_iptm=(
                "mean_pair_iptm",
                "mean",
            ),

            # QC
            n_severe_overlap=(
                "n_severe_overlap",
                "sum",
            ),

            any_af3_clash=(
                "any_af3_clash",
                "max",
            ),

            any_nonstandard_job=(
                "any_nonstandard_job",
                "max",
            ),

            job_ids=(
                "job_id",
                join_unique_strings,
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
    observations,
    job_level,
    pair_level,
):
    """Print basic integrity and QC statistics."""

    print("=" * 70)
    print("RAW 8 Å OBSERVATIONS")
    print("=" * 70)

    print(
        "Raw observations:",
        len(observations),
    )

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

    print(
        "Raw observations with severe overlap:",
        observations["severe_overlap"].sum(),
    )

    print()

    print("=" * 70)
    print("JOB LEVEL")
    print("=" * 70)

    print(
        "Job + pair rows:",
        len(job_level),
    )

    print(
        "Job+pair entries with >1 archived rerun:",
        (job_level["n_reruns"] > 1).sum(),
    )

    print(
        "Job+pair entries with AF3 clash flag:",
        (job_level["any_af3_clash"] > 0).sum(),
    )

    print()

    print("=" * 70)
    print("PAIR LEVEL")
    print("=" * 70)

    print(
        "Final protein pairs:",
        len(pair_level),
    )

    print(
        "Pairs appearing in >1 pooled job:",
        (pair_level["n_jobs"] > 1).sum(),
    )

    print(
        "Pairs with any severe overlap:",
        pair_level["has_severe_overlap"].sum(),
    )

    print(
        "Pairs with any AF3 clash flag:",
        (pair_level["any_af3_clash"] > 0).sum(),
    )

    print(
        "Pairs involving a non-standard _2 job:",
        pair_level["any_nonstandard_job"].sum(),
    )

    print()

    print("8 Å contact-count distribution:")

    print(
        pair_level["contact_count"].describe(
            percentiles=[
                0.5,
                0.75,
                0.9,
                0.95,
                0.99,
            ]
        )
    )

    print()

    print(
        "Pairs with mean contact count > 0:",
        (pair_level["contact_count"] > 0).sum(),
    )

    print(
        "Percentage with mean contact count > 0:",
        100
        * (pair_level["contact_count"] > 0).mean(),
    )

    print()

    print("Short-chain length summary:")

    print(
        pair_level[
            "short_chain_length"
        ].describe(
            percentiles=[
                0.1,
                0.25,
                0.5,
                0.75,
                0.9,
            ]
        )
    )


def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    observations = load_observations()

    job_level = build_job_level(
        observations
    )

    pair_level = build_pair_level(
        job_level
    )

    job_level.to_csv(
        JOB_LEVEL_OUTPUT,
        index=False,
    )

    pair_level.to_csv(
        PAIR_LEVEL_OUTPUT,
        index=False,
    )

    print_summary(
        observations,
        job_level,
        pair_level,
    )

    print()
    print("Saved:")
    print(JOB_LEVEL_OUTPUT)
    print(PAIR_LEVEL_OUTPUT)


if __name__ == "__main__":
    main()