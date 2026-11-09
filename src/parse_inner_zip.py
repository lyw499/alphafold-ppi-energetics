import csv
import io
import json
import re
import sys
import zipfile

from Bio.PDB import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch


CUTOFF = 5.0


def count_contacts(cif_bytes):
    """Count unique inter-chain residue-residue contacts."""

    # MMCIFParser expects a text file-like object
    cif_text = cif_bytes.decode("utf-8")
    handle = io.StringIO(cif_text)

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("af3", handle)
    model = structure[0]

    chains = list(model)

    if len(chains) != 2:
        raise ValueError(f"Expected 2 chains, found {len(chains)}")

    chain_a, chain_b = chains

    residues_a = [r for r in chain_a if r.id[0] == " "]
    residues_b = [r for r in chain_b if r.id[0] == " "]

    atoms_b = [
        atom
        for residue in residues_b
        for atom in residue
        if atom.element != "H"
    ]

    neighbor_search = NeighborSearch(atoms_b)

    contacts = set()

    for residue_a in residues_a:
        for atom_a in residue_a:

            if atom_a.element == "H":
                continue

            nearby_atoms = neighbor_search.search(
                atom_a.coord,
                CUTOFF,
                level="A",
            )

            for atom_b in nearby_atoms:
                residue_b = atom_b.get_parent()

                contacts.add(
                    (residue_a.id[1], residue_b.id[1])
                )

    interface_a = {a for a, b in contacts}
    interface_b = {b for a, b in contacts}

    return {
        "chain_a_length": len(residues_a),
        "chain_b_length": len(residues_b),
        "contact_count": len(contacts),
        "interface_residues_a": len(interface_a),
        "interface_residues_b": len(interface_b),
    }


def main(zip_path, output_path):

    results = []

    with zipfile.ZipFile(zip_path) as zf:

        names = zf.namelist()

        # Find all summary confidence files.
        summary_pattern = re.compile(
            r"^(mgen_pairs_v1_\d+)/.*_summary_confidences_(\d+)\.json$"
        )

        pairs = {}

        for name in names:
            match = summary_pattern.match(name)

            if match:
                pair_dir = match.group(1)
                model_num = int(match.group(2))

                pairs.setdefault(pair_dir, []).append(
                    (model_num, name)
                )

        print(f"Found {len(pairs)} protein pairs.")

        for index, pair_dir in enumerate(sorted(pairs), start=1):

            try:
                models = []

                # Read confidence for every available model.
                for model_num, summary_name in pairs[pair_dir]:

                    summary = json.loads(
                        zf.read(summary_name)
                    )

                    models.append({
                        "model": model_num,
                        "summary": summary,
                    })

                # Select the highest-ranking AF3 model.
                best = max(
                    models,
                    key=lambda x: x["summary"]["ranking_score"],
                )

                best_model = best["model"]
                summary = best["summary"]

                # Locate corresponding CIF.
                cif_suffix = f"_model_{best_model}.cif"

                cif_candidates = [
                    name
                    for name in names
                    if name.startswith(pair_dir + "/")
                    and name.endswith(cif_suffix)
                ]

                if len(cif_candidates) != 1:
                    raise ValueError(
                        f"Expected one CIF for model {best_model}, "
                        f"found {len(cif_candidates)}"
                    )

                cif_bytes = zf.read(cif_candidates[0])

                contact_data = count_contacts(cif_bytes)

                result = {
                    "pair_id": pair_dir,
                    "best_model": best_model,
                    "ranking_score": summary.get("ranking_score"),
                    "iptm": summary.get("iptm"),
                    "ptm": summary.get("ptm"),
                    **contact_data,
                }

                results.append(result)

                print(
                    f"[{index}/{len(pairs)}] "
                    f"{pair_dir}: "
                    f"{contact_data['contact_count']} contacts, "
                    f"ipTM={summary.get('iptm')}"
                )

            except Exception as exc:
                print(
                    f"[{index}/{len(pairs)}] "
                    f"ERROR {pair_dir}: {exc}",
                    file=sys.stderr,
                )

    fieldnames = [
        "pair_id",
        "best_model",
        "ranking_score",
        "iptm",
        "ptm",
        "chain_a_length",
        "chain_b_length",
        "contact_count",
        "interface_residues_a",
        "interface_residues_b",
    ]

    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(results)

    print()
    print(f"Successfully parsed {len(results)} pairs.")
    print(f"Results written to: {output_path}")


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print(
            "Usage: python parse_inner_zip.py "
            "<input.zip> <output.csv>"
        )
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
