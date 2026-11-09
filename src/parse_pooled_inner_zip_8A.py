import csv
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
from scipy.spatial import cKDTree


# Main interface-contact definition:
# two residues are in contact if any pair of non-H atoms
# from the two residues is <= 8 Å apart.
CUTOFF = 8.0

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


def residue_id_string(residue):
    """
    Convert a Bio.PDB residue ID to a compact string.

    For normal AF3 protein residues this will usually just be
    the residue number, e.g. "42".

    If an insertion code exists, it is retained.
    """

    seq_id = residue.id[1]
    insertion_code = residue.id[2].strip()

    if insertion_code:
        return f"{seq_id}{insertion_code}"

    return str(seq_id)


def count_pair_contacts(
    chain_a,
    chain_b,
    cutoff=CUTOFF,
):
    """
    Count unique residue-residue contacts and calculate the
    minimum heavy-atom distance between two chains.

    Two residues are considered in contact if at least one
    pair of non-hydrogen atoms is <= cutoff Angstrom apart.

    Also returns the residue IDs participating in the interface.
    """

    residues_a = get_residues(chain_a)
    residues_b = get_residues(chain_b)

    atoms_a = get_heavy_atoms(residues_a)
    atoms_b = get_heavy_atoms(residues_b)

    if not atoms_a or not atoms_b:
        return {
            "length_a": len(residues_a),
            "length_b": len(residues_b),

            "contact_count": 0,

            "interface_residues_a": 0,
            "interface_residues_b": 0,

            "interface_residue_ids_a": "",
            "interface_residue_ids_b": "",

            "min_interchain_distance": None,
        }

    # --------------------------------------------------
    # Minimum heavy-atom distance
    # --------------------------------------------------

    coords_a = np.asarray(
        [atom.coord for atom in atoms_a],
        dtype=float,
    )

    coords_b = np.asarray(
        [atom.coord for atom in atoms_b],
        dtype=float,
    )

    tree_b = cKDTree(coords_b)

    nearest_distances, _ = tree_b.query(
        coords_a,
        k=1,
    )

    min_interchain_distance = float(
        np.min(nearest_distances)
    )

    # --------------------------------------------------
    # Residue-residue contacts
    # --------------------------------------------------

    neighbor_search = NeighborSearch(atoms_b)

    # Store actual residue objects as contact keys.
    # This avoids ambiguity during interface-residue extraction.
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
                        residue_a,
                        residue_b,
                    )
                )

    interface_a = {
        residue_a
        for residue_a, residue_b in contacts
    }

    interface_b = {
        residue_b
        for residue_a, residue_b in contacts
    }

    # Sort by residue sequence number for readable output.
    interface_a_sorted = sorted(
        interface_a,
        key=lambda r: (
            r.id[1],
            r.id[2],
        ),
    )

    interface_b_sorted = sorted(
        interface_b,
        key=lambda r: (
            r.id[1],
            r.id[2],
        ),
    )

    interface_ids_a = ";".join(
        residue_id_string(residue)
        for residue in interface_a_sorted
    )

    interface_ids_b = ";".join(
        residue_id_string(residue)
        for residue in interface_b_sorted
    )

    return {
        "length_a": len(residues_a),
        "length_b": len(residues_b),

        "contact_count": len(contacts),

        "interface_residues_a": len(interface_a),
        "interface_residues_b": len(interface_b),

        "interface_residue_ids_a":
            interface_ids_a,

        "interface_residue_ids_b":
            interface_ids_b,

        "min_interchain_distance":
            min_interchain_distance,
    }


def make_interface_orientation_fields(
    chain_a,
    chain_b,
    protein_a,
    protein_b,
    contact_data,
):
    """
    Convert A/B-oriented interface information into
    short-chain / long-chain-oriented information.

    If the chains have identical lengths, chain A is treated
    as the short chain deterministically.
    """

    length_a = contact_data["length_a"]
    length_b = contact_data["length_b"]

    interface_a = (
        contact_data["interface_residues_a"]
    )

    interface_b = (
        contact_data["interface_residues_b"]
    )

    ids_a = (
        contact_data["interface_residue_ids_a"]
    )

    ids_b = (
        contact_data["interface_residue_ids_b"]
    )

    if length_a <= length_b:

        short_chain_id = chain_a.id
        long_chain_id = chain_b.id

        short_protein = protein_a
        long_protein = protein_b

        short_chain_length = length_a
        long_chain_length = length_b

        interface_residues_short = interface_a
        interface_residues_long = interface_b

        interface_residue_ids_short = ids_a
        interface_residue_ids_long = ids_b

    else:

        short_chain_id = chain_b.id
        long_chain_id = chain_a.id

        short_protein = protein_b
        long_protein = protein_a

        short_chain_length = length_b
        long_chain_length = length_a

        interface_residues_short = interface_b
        interface_residues_long = interface_a

        interface_residue_ids_short = ids_b
        interface_residue_ids_long = ids_a

    if short_chain_length > 0:
        fraction_interface_short = (
            interface_residues_short
            / short_chain_length
        )
    else:
        fraction_interface_short = None

    if long_chain_length > 0:
        fraction_interface_long = (
            interface_residues_long
            / long_chain_length
        )
    else:
        fraction_interface_long = None

    total_interface_residues = (
        interface_residues_short
        + interface_residues_long
    )

    if total_interface_residues > 0:
        contact_density = (
            contact_data["contact_count"]
            / total_interface_residues
        )
    else:
        contact_density = 0.0

    return {
        "short_chain_id":
            short_chain_id,

        "long_chain_id":
            long_chain_id,

        "short_protein":
            short_protein,

        "long_protein":
            long_protein,

        "short_chain_length":
            short_chain_length,

        "long_chain_length":
            long_chain_length,

        "interface_residues_short":
            interface_residues_short,

        "interface_residues_long":
            interface_residues_long,

        "fraction_interface_short":
            fraction_interface_short,

        "fraction_interface_long":
            fraction_interface_long,

        "interface_residue_ids_short":
            interface_residue_ids_short,

        "interface_residue_ids_long":
            interface_residue_ids_long,

        "contact_density":
            contact_density,
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


def find_job_request(
    zf,
    names,
    job_dir,
):
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

    metadata_job = job_dir

    if metadata_job not in job_to_proteins:
        metadata_job = re.sub(
            r"_\d+$",
            "",
            job_dir,
        )

    if metadata_job not in job_to_proteins:
        raise ValueError(
            f"{job_dir} not found in EV2 "
            f"pool_descriptions"
        )

    expected_proteins = (
        job_to_proteins[metadata_job]
    )

    if mapped_proteins != expected_proteins:
        raise ValueError(
            f"Protein order mismatch for {job_dir}"
        )

    return mapped_proteins


def main(
    zip_path,
    output_path,
):

    (
        seq_to_mg,
        mg_to_annotation,
        job_to_proteins,
    ) = load_metadata()

    results = []
    errors = []

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

                # --------------------------------------------------
                # Read + validate protein mapping
                # --------------------------------------------------

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

                # --------------------------------------------------
                # Read confidence summaries
                # --------------------------------------------------

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

                # --------------------------------------------------
                # Select top-ranking structural model
                # --------------------------------------------------

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

                # Save AF3 clash annotation if present.
                has_clash = best_summary.get(
                    "has_clash"
                )

                # --------------------------------------------------
                # Load corresponding CIF
                # --------------------------------------------------

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

                # --------------------------------------------------
                # Average pairwise ipTM across AF3 models
                # --------------------------------------------------

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

                # --------------------------------------------------
                # Calculate every unique chain pair
                # --------------------------------------------------

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

                    orientation_data = (
                        make_interface_orientation_fields(
                            chain_a,
                            chain_b,
                            protein_a,
                            protein_b,
                            contact_data,
                        )
                    )

                    # Canonical protein pair ID
                    protein_pair = "::::".join(
                        sorted(
                            [
                                protein_a,
                                protein_b,
                            ]
                        )
                    )

                    results.append({
                        "job_id":
                            job_dir,

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

                        "has_clash":
                            has_clash,

                        "contact_cutoff_A":
                            CUTOFF,

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

                        "interface_residue_ids_a":
                            contact_data[
                                "interface_residue_ids_a"
                            ],

                        "interface_residue_ids_b":
                            contact_data[
                                "interface_residue_ids_b"
                            ],

                        "min_interchain_distance":
                            contact_data[
                                "min_interchain_distance"
                            ],

                        "mean_pair_iptm":
                            float(
                                mean_pair_iptm[
                                    i,
                                    j,
                                ]
                            ),

                        # ------------------------------------------
                        # Short-chain / long-chain fields
                        # ------------------------------------------

                        "short_chain_id":
                            orientation_data[
                                "short_chain_id"
                            ],

                        "long_chain_id":
                            orientation_data[
                                "long_chain_id"
                            ],

                        "short_protein":
                            orientation_data[
                                "short_protein"
                            ],

                        "long_protein":
                            orientation_data[
                                "long_protein"
                            ],

                        "short_chain_length":
                            orientation_data[
                                "short_chain_length"
                            ],

                        "long_chain_length":
                            orientation_data[
                                "long_chain_length"
                            ],

                        "interface_residues_short":
                            orientation_data[
                                "interface_residues_short"
                            ],

                        "interface_residues_long":
                            orientation_data[
                                "interface_residues_long"
                            ],

                        "fraction_interface_short":
                            orientation_data[
                                "fraction_interface_short"
                            ],

                        "fraction_interface_long":
                            orientation_data[
                                "fraction_interface_long"
                            ],

                        "interface_residue_ids_short":
                            orientation_data[
                                "interface_residue_ids_short"
                            ],

                        "interface_residue_ids_long":
                            orientation_data[
                                "interface_residue_ids_long"
                            ],

                        "contact_density":
                            orientation_data[
                                "contact_density"
                            ],
                    })

                    job_rows += 1

                print(
                    f"  Finished: "
                    f"{job_rows} pairs"
                )

            except Exception as exc:

                errors.append(
                    (job_dir, str(exc))
                )

                print(
                    f"  ERROR: {exc}",
                    file=sys.stderr,
                )

    # --------------------------------------------------
    # Write output CSV
    # --------------------------------------------------

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
        "has_clash",

        "contact_cutoff_A",

        "length_a",
        "length_b",

        "contact_count",

        "interface_residues_a",
        "interface_residues_b",

        "interface_residue_ids_a",
        "interface_residue_ids_b",

        "min_interchain_distance",

        "mean_pair_iptm",

        "short_chain_id",
        "long_chain_id",

        "short_protein",
        "long_protein",

        "short_chain_length",
        "long_chain_length",

        "interface_residues_short",
        "interface_residues_long",

        "fraction_interface_short",
        "fraction_interface_long",

        "interface_residue_ids_short",
        "interface_residue_ids_long",

        "contact_density",
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

    if errors:

        print()

        print(
            f"WARNING: {len(errors)} jobs failed.",
            file=sys.stderr,
        )

        for job_dir, message in errors:

            print(
                f"  {job_dir}: {message}",
                file=sys.stderr,
            )

        sys.exit(2)


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