"""Map EV4 MG locus tags to UniProt and find shared PDB entries.

Input is the CSV produced by ``select_ev4_pairs.py``. A shared PDB entry is a
candidate experimental complex because both UniProt accessions are cross-
referenced to that entry. Biological-assembly correctness must still be checked
manually before using a structure as an energetic ground truth.

Example
-------
python src/query_pdb_coverage.py \
    --input data/processed/string_supported_pairs.csv \
    --output-dir data/processed/pdb_coverage
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import io
import json
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
USER_AGENT = "alphafold-ppi-energetics/1.0 (academic research)"


def get_text(url: str, timeout: int = 60, attempts: int = 3) -> str:
    """Download text with a short retry policy."""
    request = Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(1, attempts + 1):
        try:
            with urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8")
        except (HTTPError, URLError, TimeoutError) as error:
            if attempt == attempts:
                raise RuntimeError(f"Could not download {url}: {error}") from error
            time.sleep(2 ** (attempt - 1))
    raise AssertionError("Unreachable")


def normalise_locus_tag(value: object) -> str:
    """Make EV4 'MG_191' and UniProt 'MG191' comparable."""
    return re.sub(r"[^A-Z0-9]", "", str(value).upper().strip())


def split_pdb_ids(value: object) -> list[str]:
    """Parse UniProt's semicolon-separated PDB cross-references."""
    if pd.isna(value):
        return []
    return sorted({item.strip().upper() for item in str(value).split(";") if item.strip()})


def download_uniprot_proteome(organism_id: int) -> pd.DataFrame:
    """Download locus-tag, accession, protein, length and PDB mappings."""
    parameters = urlencode(
        {
            "query": f"organism_id:{organism_id}",
            "format": "tsv",
            "fields": (
                "accession,reviewed,gene_primary,gene_oln,"
                "protein_name,length,xref_pdb"
            ),
        }
    )
    table = pd.read_csv(io.StringIO(get_text(f"{UNIPROT_STREAM_URL}?{parameters}")), sep="\t")
    table = table.rename(
        columns={
            "Entry": "uniprot_accession",
            "Reviewed": "reviewed",
            "Gene Names (primary)": "gene_primary",
            "Gene Names (ordered locus)": "ordered_locus",
            "Protein names": "protein_name",
            "Length": "sequence_length",
            "PDB": "pdb_ids",
        }
    )
    required = {
        "uniprot_accession",
        "ordered_locus",
        "protein_name",
        "sequence_length",
        "pdb_ids",
    }
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Unexpected UniProt response; missing columns: {sorted(missing)}")

    table["locus_key"] = table["ordered_locus"].map(normalise_locus_tag)
    table["pdb_id_list"] = table["pdb_ids"].map(split_pdb_ids)
    table["pdb_ids"] = table["pdb_id_list"].map(lambda values: ";".join(values))
    table["pdb_entry_count"] = table["pdb_id_list"].map(len)
    return table


def build_locus_mapping(uniprot: pd.DataFrame) -> dict[str, dict[str, object]]:
    """Build one mapping record per ordered locus tag."""
    valid = uniprot.loc[uniprot["locus_key"].str.match(r"^MG\d+$", na=False)].copy()
    duplicates = valid.loc[valid["locus_key"].duplicated(keep=False), "locus_key"]
    if not duplicates.empty:
        names = ", ".join(sorted(duplicates.unique()))
        raise ValueError(f"Multiple UniProt records found for locus tags: {names}")
    return valid.set_index("locus_key").to_dict(orient="index")


def annotate_pairs(pairs: pd.DataFrame, mapping: dict[str, dict[str, object]]) -> pd.DataFrame:
    """Add UniProt accessions and shared PDB entry IDs to every input pair."""
    required = {"protein_1", "protein_2"}
    missing = required.difference(pairs.columns)
    if missing:
        raise ValueError(f"Input pair table is missing columns: {sorted(missing)}")

    rows: list[dict[str, object]] = []
    for record in pairs.to_dict(orient="records"):
        key_1 = normalise_locus_tag(record["protein_1"])
        key_2 = normalise_locus_tag(record["protein_2"])
        map_1 = mapping.get(key_1, {})
        map_2 = mapping.get(key_2, {})
        pdb_1 = set(map_1.get("pdb_id_list", []))
        pdb_2 = set(map_2.get("pdb_id_list", []))
        shared = sorted(pdb_1.intersection(pdb_2))

        record.update(
            {
                "uniprot_1": map_1.get("uniprot_accession", pd.NA),
                "uniprot_2": map_2.get("uniprot_accession", pd.NA),
                "protein_name_1": map_1.get("protein_name", pd.NA),
                "protein_name_2": map_2.get("protein_name", pd.NA),
                "sequence_length_1": map_1.get("sequence_length", pd.NA),
                "sequence_length_2": map_2.get("sequence_length", pd.NA),
                "uniprot_mapping_complete": bool(map_1 and map_2),
                "pdb_entries_protein_1": ";".join(sorted(pdb_1)),
                "pdb_entries_protein_2": ";".join(sorted(pdb_2)),
                "shared_pdb_ids": ";".join(shared),
                "shared_pdb_count": len(shared),
                "has_shared_pdb_entry": bool(shared),
            }
        )
        rows.append(record)
    return pd.DataFrame(rows)


def fetch_rcsb_metadata(pdb_id: str) -> dict[str, object]:
    """Retrieve entry-level method, resolution, title and assembly count."""
    data = json.loads(get_text(RCSB_ENTRY_URL.format(pdb_id=pdb_id)))
    methods = sorted(
        {item.get("method", "") for item in data.get("exptl", []) if item.get("method")}
    )
    info = data.get("rcsb_entry_info", {})
    resolutions = info.get("resolution_combined") or []
    return {
        "pdb_id": pdb_id,
        "experimental_method": ";".join(methods),
        "resolution_angstrom": min(resolutions) if resolutions else pd.NA,
        "assembly_count": info.get("assembly_count", pd.NA),
        "structure_title": data.get("struct", {}).get("title", pd.NA),
        "rcsb_url": f"https://www.rcsb.org/structure/{pdb_id}",
    }


def build_complex_table(coverage: pd.DataFrame) -> pd.DataFrame:
    """Create one row for every pair-by-shared-PDB combination."""
    covered = coverage.loc[coverage["has_shared_pdb_entry"]].copy()
    if covered.empty:
        return pd.DataFrame(
            columns=[
                "protein_1", "protein_2", "uniprot_1", "uniprot_2", "pdb_id",
                "experimental_method", "resolution_angstrom", "assembly_count",
                "structure_title", "rcsb_url", "manual_biological_assembly_check",
            ]
        )

    pdb_ids = sorted(
        {pdb_id for value in covered["shared_pdb_ids"] for pdb_id in value.split(";")}
    )
    with ThreadPoolExecutor(max_workers=min(8, len(pdb_ids))) as executor:
        metadata = dict(zip(pdb_ids, executor.map(fetch_rcsb_metadata, pdb_ids)))
    rows: list[dict[str, object]] = []
    carry = [
        "protein_1", "protein_2", "uniprot_1", "uniprot_2",
        "string_experimental", "size_corrected_iptm", "combined_score",
    ]
    for record in covered.to_dict(orient="records"):
        for pdb_id in record["shared_pdb_ids"].split(";"):
            row = {column: record.get(column, pd.NA) for column in carry}
            row.update(metadata[pdb_id])
            row["manual_biological_assembly_check"] = "pending"
            rows.append(row)
    return pd.DataFrame(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Selected PPI CSV")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--organism-id",
        type=int,
        default=243273,
        help="UniProt organism taxonomy ID (default: 243273)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not args.input.exists():
        raise FileNotFoundError(f"Input pair table not found: {args.input}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(args.input)
    uniprot = download_uniprot_proteome(args.organism_id)
    mapping = build_locus_mapping(uniprot)
    coverage = annotate_pairs(pairs, mapping)
    complexes = build_complex_table(coverage)

    mapping_output = args.output_dir / "uniprot_mapping.csv"
    coverage_output = args.output_dir / "pdb_pair_coverage.csv"
    complex_output = args.output_dir / "pdb_complex_candidates.csv"

    uniprot.drop(columns=["pdb_id_list"]).to_csv(mapping_output, index=False)
    coverage.to_csv(coverage_output, index=False)
    complexes.to_csv(complex_output, index=False)

    print(f"Input protein pairs: {len(coverage):,}")
    print(f"Complete UniProt mappings: {int(coverage['uniprot_mapping_complete'].sum()):,}")
    print(f"Pairs sharing at least one PDB entry: {int(coverage['has_shared_pdb_entry'].sum()):,}")
    print(f"Pair-PDB candidates requiring assembly review: {len(complexes):,}")
    print(f"Saved outputs in: {args.output_dir}")


if __name__ == "__main__":
    main()
