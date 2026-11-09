from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


INPUT = Path(
    "results/pair_level_contacts_8A_with_best_iptm.csv"
)

OUTPUT_DIR = Path(
    "results/figures_8A"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# --------------------------------------------------
# Load pair-level data
# --------------------------------------------------

df = pd.read_csv(INPUT)

required_columns = [
    "contact_count",
    "best_model_pair_iptm",
    "mean_pair_iptm",
]

missing = [
    c for c in required_columns
    if c not in df.columns
]

if missing:
    raise ValueError(
        f"Missing required columns: {missing}"
    )


# --------------------------------------------------
# Keep non-zero interfaces
# --------------------------------------------------

plot_df = df[
    (df["contact_count"] > 0)
    & df["best_model_pair_iptm"].notna()
    & df["mean_pair_iptm"].notna()
].copy()

print(
    f"Non-zero interfaces: "
    f"{len(plot_df):,}"
)


# ==================================================
# 1. Best-model pair ipTM correlations
# ==================================================

pearson_best = plot_df[
    [
        "contact_count",
        "best_model_pair_iptm",
    ]
].corr(
    method="pearson"
).iloc[0, 1]

spearman_best = plot_df[
    [
        "contact_count",
        "best_model_pair_iptm",
    ]
].corr(
    method="spearman"
).iloc[0, 1]


# ==================================================
# 2. Mean-model pair ipTM correlations
# ==================================================

pearson_mean = plot_df[
    [
        "contact_count",
        "mean_pair_iptm",
    ]
].corr(
    method="pearson"
).iloc[0, 1]

spearman_mean = plot_df[
    [
        "contact_count",
        "mean_pair_iptm",
    ]
].corr(
    method="spearman"
).iloc[0, 1]


print()
print(
    "Contact count vs BEST-MODEL pair ipTM"
)
print(
    f"Pearson r = {pearson_best:.4f}"
)
print(
    f"Spearman rho = {spearman_best:.4f}"
)

print()
print(
    "Contact count vs MEAN pair ipTM"
)
print(
    f"Pearson r = {pearson_mean:.4f}"
)
print(
    f"Spearman rho = {spearman_mean:.4f}"
)


# --------------------------------------------------
# P99 display cutoff
#
# Important:
# correlations above use ALL non-zero interfaces.
# P99 only controls visualization.
# --------------------------------------------------

p99 = plot_df[
    "contact_count"
].quantile(
    0.99
)

display_df = plot_df[
    plot_df["contact_count"] <= p99
].copy()

print()
print(
    f"P99 contact count = {p99:.1f}"
)

print(
    f"Displayed pairs = "
    f"{len(display_df):,}"
)


# ==================================================
# Figure 1:
# Contact count vs BEST-MODEL pair ipTM
# ==================================================

plt.figure(
    figsize=(9, 6)
)

plt.scatter(
    display_df["contact_count"],
    display_df["best_model_pair_iptm"],
    s=8,
    alpha=0.25,
)

plt.xlabel(
    "Mean inter-chain heavy-atom contact count "
    "(8 Å cutoff)"
)

plt.ylabel(
    "Best-model pair ipTM"
)

plt.title(
    "Interface size vs best-model pair ipTM"
)

plt.text(
    0.97,
    0.95,
    f"N = {len(plot_df):,}\n"
    f"Pearson r = {pearson_best:.3f}\n"
    f"Spearman ρ = {spearman_best:.3f}\n"
    f"Display ≤ P99 ({p99:.1f})",
    transform=plt.gca().transAxes,
    ha="right",
    va="top",
)

plt.tight_layout()

best_output = (
    OUTPUT_DIR
    / "contact_count_vs_best_model_iptm_8A.png"
)

plt.savefig(
    best_output,
    dpi=300,
)

plt.close()


# ==================================================
# Figure 2:
# Contact count vs MEAN pair ipTM
# sensitivity / comparison
# ==================================================

plt.figure(
    figsize=(9, 6)
)

plt.scatter(
    display_df["contact_count"],
    display_df["mean_pair_iptm"],
    s=8,
    alpha=0.25,
)

plt.xlabel(
    "Mean inter-chain heavy-atom contact count "
    "(8 Å cutoff)"
)

plt.ylabel(
    "Mean pair ipTM across AF3 models"
)

plt.title(
    "Interface size vs mean pair ipTM"
)

plt.text(
    0.97,
    0.95,
    f"N = {len(plot_df):,}\n"
    f"Pearson r = {pearson_mean:.3f}\n"
    f"Spearman ρ = {spearman_mean:.3f}\n"
    f"Display ≤ P99 ({p99:.1f})",
    transform=plt.gca().transAxes,
    ha="right",
    va="top",
)

plt.tight_layout()

mean_output = (
    OUTPUT_DIR
    / "contact_count_vs_mean_iptm_8A.png"
)

plt.savefig(
    mean_output,
    dpi=300,
)

plt.close()


# ==================================================
# Save correlation summary
# ==================================================

summary = pd.DataFrame(
    [
        {
            "iptm_metric":
                "best_model_pair_iptm",
            "n_nonzero_pairs":
                len(plot_df),
            "pearson_r":
                pearson_best,
            "spearman_rho":
                spearman_best,
            "p99_contact_count":
                p99,
        },
        {
            "iptm_metric":
                "mean_pair_iptm",
            "n_nonzero_pairs":
                len(plot_df),
            "pearson_r":
                pearson_mean,
            "spearman_rho":
                spearman_mean,
            "p99_contact_count":
                p99,
        },
    ]
)

summary_output = (
    OUTPUT_DIR
    / "contact_iptm_correlations_8A.csv"
)

summary.to_csv(
    summary_output,
    index=False,
)


print()
print("Saved:")
print(best_output)
print(mean_output)
print(summary_output)