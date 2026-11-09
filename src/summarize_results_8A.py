from pathlib import Path

import pandas as pd


INPUT = Path("results/pair_level_contacts_8A.csv")
OUTPUT_DIR = Path("results/summary_8A")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# --------------------------------------------------
# Load data
# --------------------------------------------------

df = pd.read_csv(INPUT)

print("Protein pairs:", len(df))

contacts = df["contact_count"]


# ==================================================
# 1. Overall summary table
# ==================================================

nonzero = df[
    df["contact_count"] > 0
].copy()

no_severe_overlap = df[
    ~df["has_severe_overlap"].astype(bool)
].copy()

summary_rows = []


def summarize(name, x):

    c = x["contact_count"]

    return {
        "dataset": name,

        "n_pairs": len(x),

        "zero_contacts_n":
            int((c == 0).sum()),

        "zero_contacts_pct":
            100 * (c == 0).mean(),

        "nonzero_contacts_n":
            int((c > 0).sum()),

        "nonzero_contacts_pct":
            100 * (c > 0).mean(),

        "mean_contacts":
            c.mean(),

        "median_contacts":
            c.median(),

        "p75_contacts":
            c.quantile(0.75),

        "p90_contacts":
            c.quantile(0.90),

        "p95_contacts":
            c.quantile(0.95),

        "p99_contacts":
            c.quantile(0.99),

        "max_contacts":
            c.max(),
    }


summary_rows.append(
    summarize(
        "all_pairs",
        df,
    )
)

summary_rows.append(
    summarize(
        "nonzero_pairs",
        nonzero,
    )
)

summary_rows.append(
    summarize(
        "exclude_severe_overlap",
        no_severe_overlap,
    )
)

summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    OUTPUT_DIR / "summary_8A.csv",
    index=False,
)

print()
print("Summary:")
print(summary.to_string(index=False))


# ==================================================
# 2. Contact-count bins
# ==================================================

bins = [
    -0.001,
    0,
    10,
    25,
    50,
    100,
    200,
    500,
    float("inf"),
]

labels = [
    "0",
    "0–10",
    "10–25",
    "25–50",
    "50–100",
    "100–200",
    "200–500",
    ">500",
]

df["contact_bin"] = pd.cut(
    df["contact_count"],
    bins=bins,
    labels=labels,
    include_lowest=True,
    right=True,
)

bin_counts = (
    df["contact_bin"]
    .value_counts(sort=False)
    .rename("n_pairs")
    .reset_index()
)

bin_counts["percentage"] = (
    100
    * bin_counts["n_pairs"]
    / len(df)
)

bin_counts.to_csv(
    OUTPUT_DIR / "contact_bins_8A.csv",
    index=False,
)

print()
print("Contact-count bins:")
print(
    bin_counts.to_string(index=False)
)


# ==================================================
# 3. Largest interfaces
# ==================================================

large_cols = [
    "protein_pair",
    "protein_a",
    "protein_b",
    "annotation_a",
    "annotation_b",

    "short_protein",
    "long_protein",

    "short_chain_length",
    "long_chain_length",

    "contact_count",

    "interface_residues_short",
    "interface_residues_long",

    "fraction_interface_short",
    "fraction_interface_long",

    "contact_density",

    "mean_pair_iptm",
    "min_interchain_distance",

    "n_jobs",
    "n_raw_observations",

    "has_severe_overlap",
    "any_af3_clash",
]

large_cols = [
    c for c in large_cols
    if c in df.columns
]

largest = (
    df
    .sort_values(
        "contact_count",
        ascending=False,
    )
    [large_cols]
    .head(100)
)

largest.to_csv(
    OUTPUT_DIR /
    "top_large_interfaces_8A.csv",
    index=False,
)


# ==================================================
# 4. Candidate small interfaces
# ==================================================

# Exploratory definition:
# short chain <= 100 aa
# long chain >= 150 aa
# real interface
# no severe overlap
#
# This is NOT yet a final peptide definition.
# It is only a candidate list for inspection.

small = df[
    (df["short_chain_length"] <= 100)
    & (df["long_chain_length"] >= 150)
    & (df["contact_count"] > 0)
    & (~df["has_severe_overlap"].astype(bool))
].copy()

# Useful additional metric:
# proportion of the short chain involved in the interface.

small = small.sort_values(
    [
        "contact_count",
        "mean_pair_iptm",
    ],
    ascending=[
        True,
        False,
    ],
)

small_cols = [
    "protein_pair",

    "short_protein",
    "long_protein",

    "annotation_a",
    "annotation_b",

    "short_chain_length",
    "long_chain_length",

    "contact_count",

    "interface_residues_short",
    "interface_residues_long",

    "fraction_interface_short",
    "fraction_interface_long",

    "contact_density",

    "mean_pair_iptm",

    "min_interchain_distance",

    "n_jobs",
    "n_raw_observations",

    "contact_count_range",

    "any_af3_clash",
]

small_cols = [
    c for c in small_cols
    if c in small.columns
]

small[
    small_cols
].to_csv(
    OUTPUT_DIR /
    "small_interface_candidates_8A.csv",
    index=False,
)

print()
print(
    "Short-chain / long-chain "
    "candidate pairs:",
    len(small),
)


# ==================================================
# 5. High-confidence small-interface candidates
# ==================================================

# Do NOT treat this as a biological truth cutoff.
# This is just a convenient shortlist.
#
# We keep mean_pair_iptm >= 0.3 as an exploratory
# confidence threshold.

small_confident = small[
    small["mean_pair_iptm"] >= 0.3
].copy()

small_confident[
    small_cols
].to_csv(
    OUTPUT_DIR /
    "small_interface_candidates_high_conf_8A.csv",
    index=False,
)

print(
    "Candidates with mean pair ipTM >= 0.3:",
    len(small_confident),
)


# ==================================================
# 6. QC table
# ==================================================

qc = df[
    df["has_severe_overlap"].astype(bool)
    | (df["any_af3_clash"] > 0)
].copy()

qc_cols = [
    "protein_pair",
    "contact_count",
    "mean_pair_iptm",
    "min_interchain_distance",
    "has_severe_overlap",
    "any_af3_clash",
    "n_severe_overlap",
    "short_chain_length",
    "long_chain_length",
]

qc_cols = [
    c for c in qc_cols
    if c in qc.columns
]

qc[
    qc_cols
].sort_values(
    "min_interchain_distance"
).to_csv(
    OUTPUT_DIR /
    "structural_qc_flags_8A.csv",
    index=False,
)

print()
print(
    "Pairs with structural QC flags:",
    len(qc),
)


print()
print("Saved files:")

for p in sorted(
    OUTPUT_DIR.glob("*.csv")
):
    print(p)