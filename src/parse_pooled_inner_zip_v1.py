import csv
import io
import json
import re
import sys
import zipfile
from itertools import combinations

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch


CUTOFF = 5.0


def load_structure(cif_bytes):
    """Load an AF3 mmCIF structure directly from bytes."""
    cif_text = cif_bytes.decode("utf-8")
    handle = io.StringIO(cif_text)

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("af3", handle)

    return structure[0]


def get_residues(chain):
    """Return standard amino-acid residues from a chain."""
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


def count_pair_contacts(chain_a, chain_b, cutoff=CUTOFF):
    """
    Count unique residue-residue contacts between two chains.

    A residue pair is a contact if at least one pair of heavy atoms
    is within the specified distance cutoff.
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

    interface_a = {a for a, b in contacts}
    interface_b = {b for a, b in contacts}

    return {
        "length_a": len(residues_a),
        "length_b": len(residues_b),
        "contact_count": len(contacts),
        "interface_residues_a": len(interface_a),
        "interface_residues_b": len(interface_b),
    }


def find_jobs(names):
    """Find pooled AF3 job directories in the ZIP."""

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

        jobs.setdefault(job_dir, {})[model_num] = name

    return jobs


def main(zip_path, output_path):

    results = []

    with zipfile.ZipFile(zip_path) as zf:

        names = zf.namelist()
        jobs = find_jobs(names)

        print(f"Found {len(jobs)} pooled AF3 jobs.")

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
                # -----------------------------------------
                # 1. Read summary confidence files
                # -----------------------------------------

                summaries = {}

                for model_num, summary_name in jobs[job_dir].items():

                    summaries[model_num] = json.loads(
                        zf.read(summary_name)
                    )

                print(
                    f"  Models found: "
                    f"{sorted(summaries)}"
                )

                # -----------------------------------------
                # 2. Select top-ranking structural model
                # -----------------------------------------

                best_model = max(
                    summaries,
                    key=lambda model_num:
                    summaries[model_num]["ranking_score"],
                )

                best_summary = summaries[best_model]

                print(
                    f"  Best model: {best_model} "
                    f"(ranking_score="
                    f"{best_summary['ranking_score']})"
                )

                # -----------------------------------------
                # 3. Find corresponding CIF
                # -----------------------------------------

                cif_suffix = f"_model_{best_model}.cif"

                cif_candidates = [
                    name
                    for name in names
                    if name.startswith(job_dir + "/")
                    and name.endswith(cif_suffix)
                ]

                if len(cif_candidates) != 1:
                    raise ValueError(
                        f"Expected one CIF for model "
                        f"{best_model}, found "
                        f"{len(cif_candidates)}"
                    )

                model = load_structure(
                    zf.read(cif_candidates[0])
                )

                chains = list(model)
                n_chains = len(chains)

                print(
                    f"  Chains: {n_chains}"
                )

                # -----------------------------------------
                # 4. Validate pairwise ipTM matrices
                # -----------------------------------------

                pair_iptm_matrices = []

                for model_num in sorted(summaries):

                    matrix = np.asarray(
                        summaries[model_num]["chain_pair_iptm"],
                        dtype=float,
                    )

                    if matrix.shape != (
                        n_chains,
                        n_chains,
                    ):
                        raise ValueError(
                            f"Model {model_num}: "
                            f"chain_pair_iptm shape "
                            f"{matrix.shape}, expected "
                            f"({n_chains}, {n_chains})"
                        )

                    pair_iptm_matrices.append(matrix)

                mean_pair_iptm = np.mean(
                    pair_iptm_matrices,
                    axis=0,
                )

                # -----------------------------------------
                # 5. Calculate every unique chain pair
                # -----------------------------------------

                expected_pairs = (
                    n_chains * (n_chains - 1) // 2
                )

                print(
                    f"  Pairwise interfaces: "
                    f"{expected_pairs}"
                )

                job_rows = 0

                for i, j in combinations(
                    range(n_chains),
                    2,
                ):

                    chain_a = chains[i]
                    chain_b = chains[j]

                    contact_data = count_pair_contacts(
                        chain_a,
                        chain_b,
                    )

                    results.append({
                        "job_id": job_dir,
                        "best_model": best_model,
                        "job_ranking_score":
                            best_summary["ranking_score"],

                        "chain_a": chain_a.id,
                        "chain_b": chain_b.id,

                        "chain_index_a": i,
                        "chain_index_b": j,

                        "length_a":
                            contact_data["length_a"],
                        "length_b":
                            contact_data["length_b"],

                        "contact_count":
                            contact_data["contact_count"],

                        "interface_residues_a":
                            contact_data[
                                "interface_residues_a"
                            ],

                        "interface_residues_b":
                            contact_data[
                                "interface_residues_b"
                            ],

                        "mean_pair_iptm":
                            float(mean_pair_iptm[i, j]),
                    })

                    job_rows += 1

                print(
                    f"  Finished: {job_rows} pairs"
                )

            except Exception as exc:

                print(
                    f"  ERROR: {exc}",
                    file=sys.stderr,
                )

    # -----------------------------------------
    # Write CSV
    # -----------------------------------------

    fieldnames = [
        "job_id",
        "best_model",
        "job_ranking_score",
        "chain_a",
        "chain_b",
        "chain_index_a",
        "chain_index_b",
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
        f"Successfully wrote {len(results)} "
        f"chain-pair rows."
    )
    print(
        f"Results: {output_path}"
    )


if __name__ == "__main__":

    if len(sys.argv) != 3:

        print(
            "Usage: python parse_pooled_inner_zip.py "
            "<input.zip> <output.csv>"
        )

        sys.exit(1)

    main(
        sys.argv[1],
        sys.argv[2],
    )
