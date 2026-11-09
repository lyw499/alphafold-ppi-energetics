"""Annotate candidate chain sources without changing counts or applying filters."""
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/contact_comparison_atom8A'


def main():
    path = OUT / 'mgen_reference_experimental.csv'
    data = pd.read_csv(path)
    sources, metadata_rows = {}, []
    for pdb in sorted(data.pdb_id.unique()):
        header = ROOT / 'data/mgen_reference_atom8A/metadata' / (pdb + '.cif')
        d = MMCIF2Dict(str(header))
        source = defaultdict(set)
        for category, name in [('_entity_src_gen', 'pdbx_gene_src_scientific_name'), ('_entity_src_nat', 'pdbx_organism_scientific')]:
            for entity, organism in zip(d.get(category + '.entity_id', []), d.get(category + '.' + name, [])):
                if organism not in ('.', '?'):
                    source[entity].add(organism)
        for entity, organisms in source.items():
            value = ' | '.join(sorted(organisms))
            sources[(pdb, entity)] = value
            metadata_rows.append(dict(pdb_id=pdb, entity_id=entity, gene_source_organisms=value,
                                      derived_from_mgen=any('genitalium' in o.lower() for o in organisms),
                                      source_metadata_path=str(header.relative_to(ROOT))))
    for i, row in data.iterrows():
        for chain in (1, 2):
            entity = row.get(f'chain{chain}_entity_id')
            entity = str(int(entity)) if pd.notna(entity) else ''
            value = sources.get((row.pdb_id, entity), '')
            data.at[i, f'chain{chain}_source_organisms'] = value
            data.at[i, f'chain{chain}_derived_from_mgen'] = 'genitalium' in value.lower()
        a, b = [bool(data.at[i, f'chain{c}_derived_from_mgen']) for c in (1, 2)]
        known = all(bool(data.at[i, f'chain{c}_source_organisms']) for c in (1, 2))
        scope = ('both_Mgen' if a and b else 'one_Mgen_partner' if a or b else 'neither_Mgen') if known else 'source_unknown'
        data.at[i, 'reference_taxonomy_scope'] = scope
        flags = [f for f in str(row.get('qc_flag', '')).split(';') if f and f != 'nan' and not f.startswith('reference_taxonomy_')]
        if scope != 'both_Mgen':
            flags.append('reference_taxonomy_' + scope)
        data.at[i, 'qc_flag'] = ';'.join(flags)
    data.to_csv(path, index=False)
    pd.DataFrame(metadata_rows).to_csv(OUT / 'reference_entity_source_metadata.csv', index=False)
    print(data.reference_taxonomy_scope.value_counts().to_string())


if __name__ == '__main__':
    main()
