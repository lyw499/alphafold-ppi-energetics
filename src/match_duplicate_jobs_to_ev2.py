from pathlib import Path
import pandas as pd
import numpy as np

PARTS_DIR = Path("results/pooled_parts")
EV2_PATH = Path("data/metadata/EV2.xlsx")
OUT_PATH = Path("results/duplicate_job_ev2_match.csv")


# ---------------------------------------------------------
# 1. Load all parsed structural observations
# ---------------------------------------------------------
frames = []

for f in sorted(PARTS_DIR.glob("*.csv")):
    d = pd.read_csv(f)
    d["source_file"] = f.name
    frames.append(d)

df = pd.concat(frames, ignore_index=True)

# Keep only official jobs: allbyall_1 ... allbyall_2027
df["job_number"] = (
    df["job_id"]
    .str.extract(r"allbyall_(\d+)$")[0]
    .astype("Int64")
)

df = df[
    df["job_number"].between(1, 2027)
].copy()

print("Official raw observations:", len(df))
print("Official jobs:", df["job_id"].nunique())


# ---------------------------------------------------------
# 2. Find jobs occurring in multiple source files
# ---------------------------------------------------------
job_sources = (
    df.groupby("job_id")["source_file"]
    .nunique()
)

duplicate_jobs = job_sources[
    job_sources > 1
].index

dup = df[
    df["job_id"].isin(duplicate_jobs)
].copy()

print("Jobs in multiple source files:", len(duplicate_jobs))
print("Rows belonging to these jobs:", len(dup))


# ---------------------------------------------------------
# 3. Load EV2 published uncorrected ipTM
# ---------------------------------------------------------
ev2 = pd.read_excel(
    EV2_PATH,
    sheet_name="uncorrected_ipTMS"
)

# EV2 values can contain multiple ipTMs separated in a cell.
# Convert every cell into a list of numeric values.
def parse_iptm_values(x):
    if pd.isna(x):
        return []

    if isinstance(x, (int, float, np.number)):
        return [float(x)]

    text = str(x)

    # Handle common separators
    for sep in [",", ";", "|"]:
        text = text.replace(sep, " ")

    values = []

    for token in text.split():
        try:
            values.append(float(token))
        except ValueError:
            pass

    return values


ev2["published_values"] = (
    ev2["uncorrected_ipTM"]
    .apply(parse_iptm_values)
)

ev2_map = dict(
    zip(
        ev2["protein_pair"],
        ev2["published_values"]
    )
)


# ---------------------------------------------------------
# 4. Compare each source observation with EV2
# ---------------------------------------------------------
def closest_ev2_difference(row):
    values = ev2_map.get(row["protein_pair"], [])

    if not values:
        return np.nan

    observed = row["mean_pair_iptm"]

    return min(
        abs(observed - x)
        for x in values
    )


dup["ev2_abs_diff"] = dup.apply(
    closest_ev2_difference,
    axis=1
)

dup["ev2_values"] = dup["protein_pair"].map(
    lambda p: ev2_map.get(p, [])
)


# ---------------------------------------------------------
# 5. Summarize source quality within each duplicated job
# ---------------------------------------------------------
summary = (
    dup.groupby(["job_id", "source_file"])
    .agg(
        n_pairs=("protein_pair", "size"),
        mean_ev2_abs_diff=("ev2_abs_diff", "mean"),
        median_ev2_abs_diff=("ev2_abs_diff", "median"),
        exact_matches=(
            "ev2_abs_diff",
            lambda x: (x < 1e-6).sum()
        ),
        mean_contact_count=("contact_count", "mean"),
        mean_min_distance=("min_interchain_distance", "mean"),
    )
    .reset_index()
)

summary["exact_match_fraction"] = (
    summary["exact_matches"] /
    summary["n_pairs"]
)

summary = summary.sort_values(
    ["job_id", "mean_ev2_abs_diff"]
)

summary.to_csv(
    OUT_PATH,
    index=False
)


# ---------------------------------------------------------
# 6. Determine best source for each duplicated job
# ---------------------------------------------------------
best = (
    summary
    .sort_values(
        ["job_id", "mean_ev2_abs_diff"]
    )
    .groupby("job_id")
    .first()
    .reset_index()
)

print("\nBest source per duplicated job:")
print(
    best[
        [
            "job_id",
            "source_file",
            "n_pairs",
            "mean_ev2_abs_diff",
            "exact_match_fraction",
        ]
    ]
    .head(30)
    .to_string(index=False)
)


print("\nMean EV2 difference across selected sources:")
print(best["mean_ev2_abs_diff"].describe())

print("\nExact-match fraction across selected sources:")
print(best["exact_match_fraction"].describe())

print("\nSaved:")
print(OUT_PATH)