import csv
from scipy.spatial import cKDTree
import io
import json
import re
import sys
import zipfile
from itertools import combinations

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch


CUTOFF = 5.0

EV1_PATH = "data/metadata/EV1.xlsx"
EV2_PATH = "data/metadata/EV2.xlsx"


def load_metadata():
    """Load protein metadata and AF3 pool descriptions."""

    ev1 = pd.read_excel(
        EV1_PATH,
        sheet_name="mgen_protein_sequences",
    )

    ev2 = pd.read_excel(
        EV2_PATH,
        sheet_name="pool_descriptions",
    )

    # sequence -> MG locus tag
    seq_to_mg = dict(
        zip(
            ev1["aa_sequence"],
            ev1["locus_tag"],
        )
    )

    # MG locus tag -> annotation
    mg_to_annotation = dict(
        zip(
            ev1["locus_tag"],
            ev1["annotation"],
        )
    )

    # job title -> ordered list of MG proteins
    job_to_proteins = {}

    for _, row in ev2.iterrows():

        job_to_proteins[row["job_title"]] = (
            row["proteins_tested"].split("::::")
        )

    return (
        seq_to_mg,
        mg_to_annotation,
        job_to_proteins,
    )


def load_structure(cif_bytes):
    """Load an AF3 mmCIF structure directly from ZIP bytes."""

    cif_text = cif_bytes.decode("utf-8")
    handle = io.StringIO(cif_text)

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(
        "af3",
        handle,
    )

    return structure[0]


def get_residues(chain):
    """Return standard amino-acid residues."""

    return [
        residue
        for residue in chain
        if residue.id[0] == " "
    ]


def get_heavy_atoms(residues):
    """Return non-hydrogen atoms."""

    return [
        atom
        for residue in residues
        for atom in residue
        if atom.element != "H"
    ]


def count_pair_contacts(
    chain_a,
    chain_b,
    cutoff=CUTOFF,
):
    """
    Count unique residue-residue contacts.

    Two residues are considered in contact if at least one
    pair of non-hydrogen atoms is <= cutoff Angstrom apart.
    """

    residues_a = get_residues(chain_a)
    residues_b = get_residues(chain_b)

    atoms_b = get_heavy_atoms(residues_b)

    if not atoms_b:
        return {
            "length_a": len(residues_a),
            "length_b": len(residues_b),
            "contact_count": 0,
            "interface_residues_a": 0,
            "interface_residues_b": 0,
        }

    neighbor_search = NeighborSearch(atoms_b)

    contacts = set()

    for residue_a in residues_a:

        for atom_a in residue_a:

            if atom_a.element == "H":
                continue

            nearby_atoms = neighbor_search.search(
                atom_a.coord,
                cutoff,
                level="A",
            )

            for atom_b in nearby_atoms:

                residue_b = atom_b.get_parent()

                contacts.add(
                    (
                        residue_a.id[1],
                        residue_b.id[1],
                    )
                )

    interface_a = {
        a for a, b in contacts
    }

    interface_b = {
        b for a, b in contacts
    }

    return {
        "length_a": len(residues_a),
        "length_b": len(residues_b),
        "contact_count": len(contacts),
        "interface_residues_a": len(interface_a),
        "interface_residues_b": len(interface_b),
    }


def find_jobs(names):
    """Find pooled AF3 job directories."""

    pattern = re.compile(
        r"^([^/]+)/.*_summary_confidences_(\d+)\.json$"
    )

    jobs = {}

    for name in names:

        match = pattern.match(name)

        if not match:
            continue

        job_dir = match.group(1)
        model_num = int(match.group(2))

        jobs.setdefault(
            job_dir,
            {},
        )[model_num] = name

    return jobs


def find_job_request(zf, names, job_dir):
    """Read the AF3 job request for one pooled job."""

    candidates = [
        name
        for name in names
        if name.startswith(job_dir + "/")
        and name.endswith("_job_request.json")
    ]

    if len(candidates) != 1:
        raise ValueError(
            f"Expected one job request, "
            f"found {len(candidates)}"
        )

    request = json.loads(
        zf.read(candidates[0])
    )

    return request[0]


def validate_and_map_proteins(
    request,
    job_dir,
    seq_to_mg,
    job_to_proteins,
):
    """
    Map AF3 sequence order to MG locus tags and verify
    against EV2 pool_descriptions.
    """

    sequences = [
        item["proteinChain"]["sequence"]
        for item in request["sequences"]
    ]

    mapped_proteins = []

    for sequence in sequences:

        mg_id = seq_to_mg.get(sequence)

        if mg_id is None:
            raise ValueError(
                "AF3 sequence not found in EV1"
            )

        mapped_proteins.append(mg_id)

    if job_dir not in job_to_proteins:
        raise ValueError(
            f"{job_dir} not found in EV2 "
            f"pool_descriptions"
        )

    expected_proteins = (
        job_to_proteins[job_dir]
    )

    if mapped_proteins != expected_proteins:

        raise ValueError(
            f"Protein order mismatch for {job_dir}"
        )

    return mapped_proteins


def main(zip_path, output_path):

    (
        seq_to_mg,
        mg_to_annotation,
        job_to_proteins,
    ) = load_metadata()

    results = []

    with zipfile.ZipFile(zip_path) as zf:

        names = zf.namelist()
        jobs = find_jobs(names)

        print(
            f"Found {len(jobs)} pooled AF3 jobs."
        )

        for job_index, job_dir in enumerate(
            sorted(jobs),
            start=1,
        ):

            print()
            print(
                f"[Job {job_index}/{len(jobs)}] "
                f"{job_dir}"
            )

            try:

                # -------------------------------
                # Read + validate protein mapping
                # -------------------------------

                request = find_job_request(
                    zf,
                    names,
                    job_dir,
                )

                proteins = validate_and_map_proteins(
                    request,
                    job_dir,
                    seq_to_mg,
                    job_to_proteins,
                )

                print(
                    f"  Protein mapping: OK "
                    f"({len(proteins)} proteins)"
                )

                # -------------------------------
                # Read confidence summaries
                # -------------------------------

                summaries = {}

                for (
                    model_num,
                    summary_name,
                ) in jobs[job_dir].items():

                    summaries[model_num] = (
                        json.loads(
                            zf.read(summary_name)
                        )
                    )

                print(
                    "  Models found: "
                    f"{sorted(summaries)}"
                )

                # -------------------------------
                # Select best structural model
                # -------------------------------

                best_model = max(
                    summaries,
                    key=lambda model_num:
                    summaries[model_num][
                        "ranking_score"
                    ],
                )

                best_summary = (
                    summaries[best_model]
                )

                print(
                    f"  Best model: {best_model} "
                    f"(ranking_score="
                    f"{best_summary['ranking_score']})"
                )

                # -------------------------------
                # Load corresponding CIF
                # -------------------------------

                cif_suffix = (
                    f"_model_{best_model}.cif"
                )

                cif_candidates = [
                    name
                    for name in names
                    if name.startswith(
                        job_dir + "/"
                    )
                    and name.endswith(
                        cif_suffix
                    )
                ]

                if len(cif_candidates) != 1:
                    raise ValueError(
                        "Expected one CIF for "
                        f"model {best_model}, found "
                        f"{len(cif_candidates)}"
                    )

                model = load_structure(
                    zf.read(
                        cif_candidates[0]
                    )
                )

                chains = list(model)
                n_chains = len(chains)

                if n_chains != len(proteins):
                    raise ValueError(
                        f"CIF has {n_chains} chains "
                        f"but metadata has "
                        f"{len(proteins)} proteins"
                    )

                print(
                    f"  Chains: {n_chains}"
                )

                # -------------------------------
                # Average pairwise ipTM
                # across AF3 models
                # -------------------------------

                pair_iptm_matrices = []

                for model_num in sorted(
                    summaries
                ):

                    matrix = np.asarray(
                        summaries[model_num][
                            "chain_pair_iptm"
                        ],
                        dtype=float,
                    )

                    if matrix.shape != (
                        n_chains,
                        n_chains,
                    ):
                        raise ValueError(
                            f"Model {model_num}: "
                            f"ipTM matrix shape "
                            f"{matrix.shape}"
                        )

                    pair_iptm_matrices.append(
                        matrix
                    )

                mean_pair_iptm = np.mean(
                    pair_iptm_matrices,
                    axis=0,
                )

                # -------------------------------
                # Calculate every chain pair
                # -------------------------------

                expected_pairs = (
                    n_chains
                    * (n_chains - 1)
                    // 2
                )

                print(
                    "  Pairwise interfaces: "
                    f"{expected_pairs}"
                )

                job_rows = 0

                for i, j in combinations(
                    range(n_chains),
                    2,
                ):

                    chain_a = chains[i]
                    chain_b = chains[j]

                    protein_a = proteins[i]
                    protein_b = proteins[j]

                    contact_data = (
                        count_pair_contacts(
                            chain_a,
                            chain_b,
                        )
                    )

                    # Canonical protein-pair ID.
                    # This lets us merge duplicate
                    # observations across pools later.
                    protein_pair = "::::".join(
                        sorted(
                            [
                                protein_a,
                                protein_b,
                            ]
                        )
                    )

                    results.append({
                        "job_id": job_dir,

                        "protein_pair":
                            protein_pair,

                        "protein_a":
                            protein_a,

                        "protein_b":
                            protein_b,

                        "annotation_a":
                            mg_to_annotation[
                                protein_a
                            ],

                        "annotation_b":
                            mg_to_annotation[
                                protein_b
                            ],

                        "chain_a":
                            chain_a.id,

                        "chain_b":
                            chain_b.id,

                        "best_model":
                            best_model,

                        "job_ranking_score":
                            best_summary[
                                "ranking_score"
                            ],

                        "length_a":
                            contact_data[
                                "length_a"
                            ],

                        "length_b":
                            contact_data[
                                "length_b"
                            ],

                        "contact_count":
                            contact_data[
                                "contact_count"
                            ],

                        "interface_residues_a":
                            contact_data[
                                "interface_residues_a"
                            ],

                        "interface_residues_b":
                            contact_data[
                                "interface_residues_b"
                            ],

                        "mean_pair_iptm":
                            float(
                                mean_pair_iptm[
                                    i,
                                    j,
                                ]
                            ),
                    })

                    job_rows += 1

                print(
                    f"  Finished: "
                    f"{job_rows} pairs"
                )

            except Exception as exc:

                print(
                    f"  ERROR: {exc}",
                    file=sys.stderr,
                )

    fieldnames = [
        "job_id",
        "protein_pair",
        "protein_a",
        "protein_b",
        "annotation_a",
        "annotation_b",
        "chain_a",
        "chain_b",
        "best_model",
        "job_ranking_score",
        "length_a",
        "length_b",
        "contact_count",
        "interface_residues_a",
        "interface_residues_b",
        "mean_pair_iptm",
    ]

    with open(
        output_path,
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(results)

    print()
    print("=" * 60)
    print(
        f"Successfully wrote "
        f"{len(results)} chain-pair rows."
    )
    print(
        f"Results: {output_path}"
    )


if __name__ == "__main__":

    if len(sys.argv) != 3:

        print(
            "Usage: python "
            "parse_pooled_inner_zip.py "
            "<input.zip> <output.csv>"
        )

        sys.exit(1)

    main(
        sys.argv[1],
        sys.argv[2],
    )
