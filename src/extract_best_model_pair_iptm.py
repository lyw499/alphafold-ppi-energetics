import csv
import io
import json
import re
import shutil
import sys
import zipfile
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# Paths
# ============================================================

OUTER_ZIP = Path(
    "data/AF3_completed_zips.zip"
)

EV1_PATH = Path(
    "data/metadata/EV1.xlsx"
)

EV2_PATH = Path(
    "data/metadata/EV2.xlsx"
)

TMP_DIR = Path(
    "tmp/best_model_iptm"
)

OUTPUT_DIR = Path(
    "results"
)

RAW_OUTPUT = OUTPUT_DIR / (
    "best_model_pair_iptm_raw.csv"
)

JOB_OUTPUT = OUTPUT_DIR / (
    "best_model_pair_iptm_job_level.csv"
)

PAIR_OUTPUT = OUTPUT_DIR / (
    "best_model_pair_iptm_pair_level.csv"
)

JOB_CONTACT_INPUT = OUTPUT_DIR / (
    "job_level_contacts_8A.csv"
)

PAIR_CONTACT_INPUT = OUTPUT_DIR / (
    "pair_level_contacts_8A.csv"
)

JOB_CONTACT_OUTPUT = OUTPUT_DIR / (
    "job_level_contacts_8A_with_best_iptm.csv"
)

PAIR_CONTACT_OUTPUT = OUTPUT_DIR / (
    "pair_level_contacts_8A_with_best_iptm.csv"
)


# ============================================================
# Metadata
# ============================================================

def load_metadata():
    """
    Load EV1 sequence mappings and EV2 pooled-job definitions.
    """

    ev1 = pd.read_excel(
        EV1_PATH,
        sheet_name="mgen_protein_sequences",
    )

    ev2 = pd.read_excel(
        EV2_PATH,
        sheet_name="pool_descriptions",
    )

    # Protein sequence -> MG locus tag
    seq_to_mg = dict(
        zip(
            ev1["aa_sequence"],
            ev1["locus_tag"],
        )
    )

    # Job title -> ordered protein list
    job_to_proteins = {}

    for _, row in ev2.iterrows():

        job_to_proteins[
            row["job_title"]
        ] = (
            row["proteins_tested"]
            .split("::::")
        )

    return (
        seq_to_mg,
        job_to_proteins,
    )


# ============================================================
# ZIP utilities
# ============================================================

def safe_name(member_name):
    """
    Convert an outer-ZIP member name into a safe filename.

    Example:
    AF3_completed_zips/folds_2025_03_18_20_18.zip

    ->
    folds_2025_03_18_20_18
    """

    name = Path(member_name).name

    if name.lower().endswith(".zip"):
        name = name[:-4]

    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        name,
    )


def find_inner_zips(outer_zip):
    """
    Find real nested ZIP files, ignoring macOS metadata.
    """

    members = []

    for info in outer_zip.infolist():

        name = info.filename

        if info.is_dir():
            continue

        if not name.lower().endswith(".zip"):
            continue

        if "__MACOSX/" in name:
            continue

        if Path(name).name.startswith("._"):
            continue

        members.append(info)

    return members


def find_jobs(names):
    """
    Find AF3 jobs and their summary_confidences files.

    Returns:
        {
            job_dir: {
                model_number: summary_filename
            }
        }
    """

    pattern = re.compile(
        r"^([^/]+)/.*_summary_confidences_(\d+)\.json$"
    )

    jobs = {}

    for name in names:

        match = pattern.match(name)

        if not match:
            continue

        job_dir = match.group(1)
        model_num = int(
            match.group(2)
        )

        jobs.setdefault(
            job_dir,
            {},
        )[model_num] = name

    return jobs


def find_job_request(
    zf,
    names,
    job_dir,
):
    """
    Read one AF3 job_request.json.
    """

    candidates = [
        name
        for name in names
        if name.startswith(
            job_dir + "/"
        )
        and name.endswith(
            "_job_request.json"
        )
    ]

    if len(candidates) != 1:

        raise ValueError(
            f"{job_dir}: expected one "
            f"job_request.json, found "
            f"{len(candidates)}"
        )

    request = json.loads(
        zf.read(candidates[0])
    )

    return request[0]


# ============================================================
# Protein mapping
# ============================================================

def validate_and_map_proteins(
    request,
    job_dir,
    seq_to_mg,
    job_to_proteins,
):
    """
    Convert AF3 sequence order to MG protein IDs.

    The sequence-derived order is validated against EV2.
    """

    sequences = [
        item["proteinChain"][
            "sequence"
        ]
        for item in request[
            "sequences"
        ]
    ]

    mapped = []

    for sequence in sequences:

        mg_id = seq_to_mg.get(
            sequence
        )

        if mg_id is None:

            raise ValueError(
                f"{job_dir}: AF3 sequence "
                "not found in EV1"
            )

        mapped.append(
            mg_id
        )

    # Handle rerun names such as:
    # 250310_mgen_allbyall_196_2
    metadata_job = job_dir

    if metadata_job not in job_to_proteins:

        metadata_job = re.sub(
            r"_\d+$",
            "",
            job_dir,
        )

    if metadata_job not in job_to_proteins:

        raise ValueError(
            f"{job_dir}: not found in "
            "EV2 pool_descriptions"
        )

    expected = (
        job_to_proteins[
            metadata_job
        ]
    )

    if mapped != expected:

        raise ValueError(
            f"{job_dir}: protein "
            "order mismatch"
        )

    return mapped


# ============================================================
# Parse one inner ZIP
# ============================================================

def parse_inner_zip(
    inner_zip_path,
    source_file,
    seq_to_mg,
    job_to_proteins,
):
    """
    Extract best-model pairwise ipTM values from one inner ZIP.

    No CIF files are read.
    """

    rows = []
    errors = []

    with zipfile.ZipFile(
        inner_zip_path
    ) as zf:

        names = zf.namelist()

        jobs = find_jobs(
            names
        )

        print(
            f"  Jobs found: {len(jobs)}"
        )

        for job_index, job_dir in enumerate(
            sorted(jobs),
            start=1,
        ):

            try:

                # ---------------------------------------------
                # Protein order
                # ---------------------------------------------

                request = find_job_request(
                    zf,
                    names,
                    job_dir,
                )

                proteins = (
                    validate_and_map_proteins(
                        request,
                        job_dir,
                        seq_to_mg,
                        job_to_proteins,
                    )
                )

                n_chains = len(
                    proteins
                )

                # ---------------------------------------------
                # Load summaries for all AF3 models
                # ---------------------------------------------

                summaries = {}

                for (
                    model_num,
                    summary_name,
                ) in jobs[
                    job_dir
                ].items():

                    summaries[
                        model_num
                    ] = json.loads(
                        zf.read(
                            summary_name
                        )
                    )

                if not summaries:

                    raise ValueError(
                        f"{job_dir}: no "
                        "summary files"
                    )

                # ---------------------------------------------
                # Select same best model used for structure
                # ---------------------------------------------

                best_model = max(
                    summaries,
                    key=lambda model_num:
                    summaries[
                        model_num
                    ]["ranking_score"],
                )

                best_summary = (
                    summaries[
                        best_model
                    ]
                )

                matrix = np.asarray(
                    best_summary[
                        "chain_pair_iptm"
                    ],
                    dtype=float,
                )

                if matrix.shape != (
                    n_chains,
                    n_chains,
                ):

                    raise ValueError(
                        f"{job_dir}: best-model "
                        f"pair ipTM matrix shape "
                        f"{matrix.shape}; expected "
                        f"{(n_chains, n_chains)}"
                    )

                # ---------------------------------------------
                # One row per protein pair
                # ---------------------------------------------

                for i, j in combinations(
                    range(n_chains),
                    2,
                ):

                    protein_a = (
                        proteins[i]
                    )

                    protein_b = (
                        proteins[j]
                    )

                    protein_pair = (
                        "::::".join(
                            sorted(
                                [
                                    protein_a,
                                    protein_b,
                                ]
                            )
                        )
                    )

                    rows.append(
                        {
                            "source_file":
                                source_file,

                            "job_id":
                                job_dir,

                            "protein_pair":
                                protein_pair,

                            "protein_a":
                                protein_a,

                            "protein_b":
                                protein_b,

                            "chain_index_a":
                                i,

                            "chain_index_b":
                                j,

                            "best_model":
                                best_model,

                            "best_model_ranking_score":
                                float(
                                    best_summary[
                                        "ranking_score"
                                    ]
                                ),

                            "best_model_pair_iptm":
                                float(
                                    matrix[
                                        i,
                                        j,
                                    ]
                                ),
                        }
                    )

            except Exception as exc:

                errors.append(
                    (
                        job_dir,
                        str(exc),
                    )
                )

                print(
                    f"    ERROR: "
                    f"{job_dir}: {exc}",
                    file=sys.stderr,
                )

    return (
        rows,
        errors,
    )


# ============================================================
# Hierarchical aggregation
# ============================================================

def build_job_level(
    raw,
):
    """
    First aggregation level:

    same job_id + protein_pair may occur in multiple archive
    source files.

    Average these reruns so duplicated archive outputs do not
    give one pooled job extra weight.
    """

    job = (
        raw
        .groupby(
            [
                "job_id",
                "protein_pair",
            ],
            as_index=False,
        )
        .agg(
            protein_a=(
                "protein_a",
                "first",
            ),

            protein_b=(
                "protein_b",
                "first",
            ),

            best_model_pair_iptm=(
                "best_model_pair_iptm",
                "mean",
            ),

            best_model_pair_iptm_min=(
                "best_model_pair_iptm",
                "min",
            ),

            best_model_pair_iptm_max=(
                "best_model_pair_iptm",
                "max",
            ),

            n_reruns=(
                "source_file",
                "nunique",
            ),

            source_files=(
                "source_file",
                lambda x:
                " | ".join(
                    sorted(
                        set(x)
                    )
                ),
            ),
        )
    )

    job[
        "best_model_pair_iptm_range"
    ] = (
        job[
            "best_model_pair_iptm_max"
        ]
        -
        job[
            "best_model_pair_iptm_min"
        ]
    )

    return job


def build_pair_level(
    job,
):
    """
    Second aggregation level:

    same protein pair may occur in different pooled jobs.

    Each pooled job receives equal weight.
    """

    pair = (
        job
        .groupby(
            "protein_pair",
            as_index=False,
        )
        .agg(
            best_model_pair_iptm=(
                "best_model_pair_iptm",
                "mean",
            ),

            best_model_pair_iptm_min=(
                "best_model_pair_iptm",
                "min",
            ),

            best_model_pair_iptm_max=(
                "best_model_pair_iptm",
                "max",
            ),

            n_best_iptm_jobs=(
                "job_id",
                "nunique",
            ),
        )
    )

    pair[
        "best_model_pair_iptm_range"
    ] = (
        pair[
            "best_model_pair_iptm_max"
        ]
        -
        pair[
            "best_model_pair_iptm_min"
        ]
    )

    return pair


# ============================================================
# Merge with 8 Å contact tables
# ============================================================

def merge_contact_tables(
    job_iptm,
    pair_iptm,
):
    """
    Add best-model pair ipTM to existing 8 Å contact tables.
    """

    if not JOB_CONTACT_INPUT.exists():

        raise FileNotFoundError(
            JOB_CONTACT_INPUT
        )

    if not PAIR_CONTACT_INPUT.exists():

        raise FileNotFoundError(
            PAIR_CONTACT_INPUT
        )

    # --------------------------------------------------------
    # Job level
    # --------------------------------------------------------

    contacts_job = pd.read_csv(
        JOB_CONTACT_INPUT
    )

    job_cols = [
        "job_id",
        "protein_pair",
        "best_model_pair_iptm",
        "best_model_pair_iptm_min",
        "best_model_pair_iptm_max",
        "best_model_pair_iptm_range",
    ]

    merged_job = (
        contacts_job
        .merge(
            job_iptm[
                job_cols
            ],
            on=[
                "job_id",
                "protein_pair",
            ],
            how="left",
            validate="one_to_one",
        )
    )

    merged_job.to_csv(
        JOB_CONTACT_OUTPUT,
        index=False,
    )

    # --------------------------------------------------------
    # Pair level
    # --------------------------------------------------------

    contacts_pair = pd.read_csv(
        PAIR_CONTACT_INPUT
    )

    pair_cols = [
        "protein_pair",
        "best_model_pair_iptm",
        "best_model_pair_iptm_min",
        "best_model_pair_iptm_max",
        "best_model_pair_iptm_range",
        "n_best_iptm_jobs",
    ]

    merged_pair = (
        contacts_pair
        .merge(
            pair_iptm[
                pair_cols
            ],
            on="protein_pair",
            how="left",
            validate="one_to_one",
        )
    )

    merged_pair.to_csv(
        PAIR_CONTACT_OUTPUT,
        index=False,
    )

    return (
        merged_job,
        merged_pair,
    )


# ============================================================
# Main
# ============================================================

def main():

    if not OUTER_ZIP.exists():

        raise FileNotFoundError(
            OUTER_ZIP
        )

    TMP_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        seq_to_mg,
        job_to_proteins,
    ) = load_metadata()

    all_rows = []
    all_errors = []

    # --------------------------------------------------------
    # Read the outer ZIP
    # --------------------------------------------------------

    print(
        f"Opening outer ZIP: "
        f"{OUTER_ZIP}"
    )

    with zipfile.ZipFile(
        OUTER_ZIP
    ) as outer:

        inner_zips = find_inner_zips(
            outer
        )

        print(
            f"Inner ZIP files: "
            f"{len(inner_zips)}"
        )

        for index, info in enumerate(
            inner_zips,
            start=1,
        ):

            member_name = (
                info.filename
            )

            stem = safe_name(
                member_name
            )

            temp_zip = (
                TMP_DIR
                / f"{stem}.zip"
            )

            print()
            print(
                "=" * 70
            )

            print(
                f"[{index}/{len(inner_zips)}] "
                f"{member_name}"
            )

            # ---------------------------------------------
            # Remove stale temp ZIP
            # ---------------------------------------------

            if temp_zip.exists():
                temp_zip.unlink()

            # ---------------------------------------------
            # Extract one inner ZIP
            # ---------------------------------------------

            print(
                "  Extracting..."
            )

            with outer.open(
                info
            ) as source:

                with open(
                    temp_zip,
                    "wb",
                ) as destination:

                    shutil.copyfileobj(
                        source,
                        destination,
                        length=(
                            16
                            * 1024
                            * 1024
                        ),
                    )

            # ---------------------------------------------
            # Parse summaries only
            # ---------------------------------------------

            print(
                "  Reading confidence summaries..."
            )

            rows, errors = (
                parse_inner_zip(
                    temp_zip,
                    Path(
                        member_name
                    ).name,
                    seq_to_mg,
                    job_to_proteins,
                )
            )

            all_rows.extend(
                rows
            )

            all_errors.extend(
                [
                    (
                        member_name,
                        job,
                        message,
                    )
                    for (
                        job,
                        message,
                    ) in errors
                ]
            )

            print(
                f"  Extracted pair rows: "
                f"{len(rows)}"
            )

            # ---------------------------------------------
            # Delete temp inner ZIP
            # ---------------------------------------------

            if temp_zip.exists():
                temp_zip.unlink()

    # --------------------------------------------------------
    # Raw table
    # --------------------------------------------------------

    raw = pd.DataFrame(
        all_rows
    )

    raw.to_csv(
        RAW_OUTPUT,
        index=False,
    )

    print()
    print(
        "=" * 70
    )

    print("RAW EXTRACTION")

    print(
        "Rows:",
        len(raw),
    )

    print(
        "Unique jobs:",
        raw[
            "job_id"
        ].nunique(),
    )

    print(
        "Unique pairs:",
        raw[
            "protein_pair"
        ].nunique(),
    )

    # --------------------------------------------------------
    # Hierarchical aggregation
    # --------------------------------------------------------

    job = build_job_level(
        raw
    )

    pair = build_pair_level(
        job
    )

    job.to_csv(
        JOB_OUTPUT,
        index=False,
    )

    pair.to_csv(
        PAIR_OUTPUT,
        index=False,
    )

    print()
    print("JOB LEVEL")

    print(
        "Job + pair rows:",
        len(job),
    )

    print(
        "Repeated job+pair entries:",
        (
            job["n_reruns"] > 1
        ).sum(),
    )

    print()
    print("PAIR LEVEL")

    print(
        "Protein pairs:",
        len(pair),
    )

    print(
        "Pairs with >1 pooled job:",
        (
            pair[
                "n_best_iptm_jobs"
            ] > 1
        ).sum(),
    )

    # --------------------------------------------------------
    # Merge with contact tables
    # --------------------------------------------------------

    merged_job, merged_pair = (
        merge_contact_tables(
            job,
            pair,
        )
    )

    print()
    print(
        "Merged job rows:",
        len(merged_job),
    )

    print(
        "Missing best-model ipTM "
        "at job level:",
        merged_job[
            "best_model_pair_iptm"
        ].isna().sum(),
    )

    print()

    print(
        "Merged pair rows:",
        len(merged_pair),
    )

    print(
        "Missing best-model ipTM "
        "at pair level:",
        merged_pair[
            "best_model_pair_iptm"
        ].isna().sum(),
    )

    # --------------------------------------------------------
    # Compare best-model versus current mean pair ipTM
    # --------------------------------------------------------

    comparison = (
        merged_pair[
            [
                "best_model_pair_iptm",
                "mean_pair_iptm",
            ]
        ]
        .dropna()
    )

    pearson = (
        comparison.corr(
            method="pearson"
        )
        .iloc[0, 1]
    )

    spearman = (
        comparison.corr(
            method="spearman"
        )
        .iloc[0, 1]
    )

    diff = (
        comparison[
            "best_model_pair_iptm"
        ]
        -
        comparison[
            "mean_pair_iptm"
        ]
    )

    print()
    print(
        "Best-model vs mean-model pair ipTM:"
    )

    print(
        f"Pearson r: "
        f"{pearson:.4f}"
    )

    print(
        f"Spearman rho: "
        f"{spearman:.4f}"
    )

    print()

    print(
        "Difference "
        "(best - mean):"
    )

    print(
        diff.describe(
            percentiles=[
                0.05,
                0.25,
                0.5,
                0.75,
                0.95,
            ]
        )
    )

    # --------------------------------------------------------
    # Errors
    # --------------------------------------------------------

    if all_errors:

        error_path = (
            OUTPUT_DIR
            / "best_model_pair_iptm_errors.tsv"
        )

        with open(
            error_path,
            "w",
        ) as handle:

            for (
                inner_zip,
                job,
                message,
            ) in all_errors:

                handle.write(
                    f"{inner_zip}\t"
                    f"{job}\t"
                    f"{message}\n"
                )

        print()
        print(
            f"WARNING: "
            f"{len(all_errors)} "
            f"job errors."
        )

        print(
            f"See: {error_path}"
        )

    print()
    print("Saved:")

    for path in [
        RAW_OUTPUT,
        JOB_OUTPUT,
        PAIR_OUTPUT,
        JOB_CONTACT_OUTPUT,
        PAIR_CONTACT_OUTPUT,
    ]:

        print(
            f"  {path}"
        )


if __name__ == "__main__":
    main()