from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


INPUT = Path("results/pair_level_contacts_8A.csv")
OUTPUT_DIR = Path("results/figures_8A")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# --------------------------------------------------
# Load pair-level data
# --------------------------------------------------

df = pd.read_csv(INPUT)

contacts = df["contact_count"]

print(f"Protein pairs: {len(df):,}")
print(f"Zero contacts: {(contacts == 0).sum():,}")
print(f"Non-zero contacts: {(contacts > 0).sum():,}")
print(
    f"Non-zero percentage: "
    f"{100 * (contacts > 0).mean():.2f}%"
)

print()
print("Contact-count summary:")
print(
    contacts.describe(
        percentiles=[
            0.50,
            0.75,
            0.90,
            0.95,
            0.99,
        ]
    )
)


# ==================================================
# Figure 1
# Full distribution
# ==================================================

plt.figure(figsize=(8, 5))

plt.hist(
    contacts,
    bins=100,
)

plt.xlabel(
    "Residue–residue contact count"
)
plt.ylabel(
    "Number of protein pairs"
)

plt.title(
    "Protein–protein interface size distribution\n"
    "8 Å heavy-atom contact cutoff"
)

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR
    / "contact_count_distribution_8A_full.png",
    dpi=300,
)

plt.close()


# ==================================================
# Figure 2
# Non-zero interfaces only
# ==================================================

nonzero = contacts[
    contacts > 0
]

plt.figure(figsize=(8, 5))

plt.hist(
    nonzero,
    bins=100,
)

plt.xlabel(
    "Residue–residue contact count"
)
plt.ylabel(
    "Number of protein pairs"
)

plt.title(
    "Non-zero protein–protein interfaces\n"
    "8 Å heavy-atom contact cutoff"
)

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR
    / "contact_count_distribution_8A_nonzero.png",
    dpi=300,
)

plt.close()


# ==================================================
# Figure 3
# Non-zero interfaces, zoomed to 99th percentile
# ==================================================

p99 = nonzero.quantile(0.99)

zoomed = nonzero[
    nonzero <= p99
]

plt.figure(figsize=(8, 5))

plt.hist(
    zoomed,
    bins=80,
)

plt.xlabel(
    "Residue–residue contact count"
)
plt.ylabel(
    "Number of protein pairs"
)

plt.title(
    "Non-zero protein–protein interfaces\n"
    "8 Å cutoff, zoomed to 99th percentile"
)

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR
    / "contact_count_distribution_8A_nonzero_p99.png",
    dpi=300,
)

plt.close()


print()
print("Saved figures:")

for path in sorted(
    OUTPUT_DIR.glob("*.png")
):
    print(path)