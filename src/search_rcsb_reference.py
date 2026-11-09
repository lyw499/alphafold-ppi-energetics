"""
Search RCSB PDB for candidate experimental heteromeric PPI references.

Selection:
1. Requested taxonomic lineage
2. Experimental structures only
3. Biological assembly contains >= 2 protein chains
4. Biological assembly contains >= 2 distinct protein entities

Outputs:
    results/rcsb_reference_counts.csv
    results/rcsb_reference_assemblies.csv
"""

from pathlib import Path
import csv
import requests


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

RESULTS_DIR = Path("results")
COUNT_OUTPUT = RESULTS_DIR / "rcsb_reference_counts.csv"
ASSEMBLY_OUTPUT = RESULTS_DIR / "rcsb_reference_assemblies.csv"


TAXONOMIC_SCOPES = {
    "M_genitalium": "Mycoplasma genitalium",
    "Mycoplasma": "Mycoplasma",
    "Mollicutes": "Mollicutes",
    "Bacteria": "Bacteria",
}


def terminal(attribute, operator, value):
    """Create one RCSB attribute-search node."""
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": attribute,
            "operator": operator,
            "value": value,
        },
    }


def build_query(taxon_name):
    """
    Experimental biological assemblies containing at least
    two protein chains AND at least two distinct protein entities.
    """

    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [

                # Taxonomic lineage
                terminal(
                    "rcsb_entity_source_organism.taxonomy_lineage.name",
                    "exact_match",
                    taxon_name,
                ),

                # Experimental structures only
                terminal(
                    "rcsb_entry_info.structure_determination_methodology",
                    "exact_match",
                    "experimental",
                ),

                # At least two protein chains
                terminal(
                    "rcsb_assembly_info.polymer_entity_instance_count_protein",
                    "greater_or_equal",
                    2,
                ),

                # At least two DISTINCT protein entities
                terminal(
                    "rcsb_assembly_info.polymer_entity_count_protein",
                    "greater_or_equal",
                    2,
                ),
            ],
        },

        "return_type": "assembly",

        "request_options": {
            "return_all_hits": True,
        },
    }


def search_rcsb(taxon_name):
    """Run one RCSB Search API query."""

    response = requests.post(
        SEARCH_URL,
        json=build_query(taxon_name),
        timeout=120,
    )

    response.raise_for_status()

    data = response.json()

    return [
        result["identifier"]
        for result in data.get("result_set", [])
    ]


def split_assembly_id(identifier):
    """
    RCSB assembly identifiers normally look like:

        1ABC-1

    Returns:
        pdb_id, assembly_id
    """

    if "-" in identifier:
        pdb_id, assembly_id = identifier.split("-", 1)
    else:
        pdb_id = identifier
        assembly_id = ""

    return pdb_id, assembly_id


def main():

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    count_rows = []
    assembly_rows = []

    for label, taxon in TAXONOMIC_SCOPES.items():

        print("\n" + "=" * 60)
        print(f"Searching: {label}")
        print(f"Taxon:     {taxon}")

        try:
            assemblies = search_rcsb(taxon)

        except requests.RequestException as exc:
            print(f"Search failed: {exc}")
            continue

        pdb_ids = set()

        for identifier in assemblies:

            pdb_id, assembly_id = split_assembly_id(identifier)

            pdb_ids.add(pdb_id)

            assembly_rows.append({
                "scope": label,
                "taxon": taxon,
                "pdb_id": pdb_id,
                "assembly_id": assembly_id,
                "rcsb_assembly_identifier": identifier,
            })

        print(
            f"Heteromeric candidate assemblies: "
            f"{len(assemblies):,}"
        )

        print(
            f"Unique PDB entries:              "
            f"{len(pdb_ids):,}"
        )

        count_rows.append({
            "scope": label,
            "taxon": taxon,
            "heteromeric_assembly_count": len(assemblies),
            "unique_pdb_count": len(pdb_ids),
        })

    # ---------------------------------------------------------
    # Write summary counts
    # ---------------------------------------------------------

    with COUNT_OUTPUT.open("w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "scope",
                "taxon",
                "heteromeric_assembly_count",
                "unique_pdb_count",
            ],
        )

        writer.writeheader()
        writer.writerows(count_rows)

    # ---------------------------------------------------------
    # Write every candidate assembly
    # ---------------------------------------------------------

    with ASSEMBLY_OUTPUT.open("w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "scope",
                "taxon",
                "pdb_id",
                "assembly_id",
                "rcsb_assembly_identifier",
            ],
        )

        writer.writeheader()
        writer.writerows(assembly_rows)

    print("\n" + "=" * 60)
    print("Finished.")
    print(f"Counts:     {COUNT_OUTPUT}")
    print(f"Assemblies: {ASSEMBLY_OUTPUT}")


if __name__ == "__main__":
    main()