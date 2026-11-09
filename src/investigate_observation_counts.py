from pathlib import Path
import pandas as pd

PARTS_DIR = Path("results/pooled_parts")

# ----------------------------
# Load all parsed CSV files
# ----------------------------
frames = []

for f in sorted(PARTS_DIR.glob("*.csv")):
    df = pd.read_csv(f)
    df["source_file"] = f.name
    frames.append(df)

df = pd.concat(frames, ignore_index=True)

print("Raw observations:", len(df))
print("Unique protein pairs:", df["protein_pair"].nunique())

# ----------------------------
# One row per job + protein pair
# ----------------------------
job_pair = df.drop_duplicates(
    subset=["job_id", "protein_pair"]
).copy()

print("\nUnique job + pair:", len(job_pair))
print("Paper reports:       160772")
print("Difference:          ", len(job_pair) - 160772)

# ----------------------------
# Count observations per job
# ----------------------------
job_counts = (
    job_pair.groupby("job_id")
    .size()
    .rename("n_pairs")
    .reset_index()
)

print("\nUnique jobs:", job_counts["job_id"].nunique())

print("\nPairs per job:")
print(job_counts["n_pairs"].describe())

# ----------------------------
# Which source files contribute?
# ----------------------------
source_counts = (
    job_pair.groupby("source_file")
    .agg(
        n_observations=("protein_pair", "size"),
        n_jobs=("job_id", "nunique"),
    )
    .sort_values("n_observations", ascending=False)
)

print("\nObservations by source file:")
print(source_counts.to_string())

# ----------------------------
# Jobs appearing in >1 source
# ----------------------------
job_sources = (
    df.groupby("job_id")["source_file"]
    .agg(lambda x: sorted(set(x)))
    .reset_index()
)

job_sources["n_sources"] = job_sources["source_file"].str.len()

multi_source_jobs = job_sources[
    job_sources["n_sources"] > 1
].copy()

print("\nJobs appearing in multiple source files:",
      len(multi_source_jobs))

print("\nExamples:")
print(multi_source_jobs.head(30).to_string(index=False))

# ----------------------------
# Look at job-number range
# ----------------------------
job_pair["job_number"] = (
    job_pair["job_id"]
    .str.extract(r"allbyall_(\d+)$")[0]
    .astype("Int64")
)

print("\nJob number range:")
print(job_pair["job_number"].describe())

print("\nHighest job numbers:")
print(
    job_pair[
        ["job_id", "job_number", "source_file"]
    ]
    .drop_duplicates()
    .sort_values("job_number", ascending=False)
    .head(30)
    .to_string(index=False)
)

# ----------------------------
# Check expected 2027 jobs
# ----------------------------
observed_numbers = set(
    job_pair["job_number"].dropna().astype(int)
)

expected_numbers = set(range(1, 2028))

missing = sorted(expected_numbers - observed_numbers)
outside = sorted(observed_numbers - expected_numbers)

print("\nMissing job numbers from 1..2027:")
print(missing)

print("\nJob numbers outside 1..2027:")
print(outside)

# ----------------------------
# Save useful tables
# ----------------------------
source_counts.to_csv(
    "results/observation_counts_by_source.csv"
)

multi_source_jobs.to_csv(
    "results/jobs_in_multiple_sources.csv",
    index=False
)

job_counts.to_csv(
    "results/pair_counts_per_job.csv",
    index=False
)