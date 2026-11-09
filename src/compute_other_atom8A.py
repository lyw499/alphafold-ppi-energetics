"""Count SAFES, downloaded reference assemblies, and existing FoldBench selections."""
import gzip
import io
import json
from itertools import combinations
from pathlib import Path

import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from atom_contact_common import atom_contacts, atom_tree, finish_record, parse_structure, protein_chains

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/contact_comparison_atom8A'


def write_rows(name, rows):
    frame = pd.DataFrame(rows)
    path = OUT / (name + '.csv')
    frame.to_csv(path, index=False)
    print(name, 'total', len(frame), 'calculable', frame.raw_included.sum(), 'qc', frame.qc_included.sum(), flush=True)


def all_pairs(path, dataset, structure_type, structure_id, metadata=None):
    metadata = metadata or {}
    base = dict(dataset=dataset, structure_type=structure_type, structure_id=structure_id,
                source_path=str(path.relative_to(ROOT)), pair_scope='all_unordered_protein_or_peptide_chain_pairs', **metadata)
    try:
        model, info = parse_structure(path)
        chains = protein_chains(model)
    except Exception as e:
        return [finish_record(dict(**base, chain_pair='', count_error=str(e)))], []
    chain_inventory, trees = [], {}
    for c in chains:
        try:
            trees[c.id] = atom_tree(c)
            chain_inventory.append(dict(dataset=dataset, structure_id=structure_id, chain_id=c.id,
                                        chain_length=trees[c.id][1], heavy_atoms=trees[c.id][0].n))
        except Exception as e:
            chain_inventory.append(dict(dataset=dataset, structure_id=structure_id, chain_id=c.id, error=str(e)))
    if len(chains) < 2:
        return [finish_record(dict(**base, **info, chain_pair='', count_error='fewer_than_two_protein_chains'))], chain_inventory
    rows = []
    for a, b in combinations(chains, 2):
        row = dict(**base, **info, chain1=a.id, chain2=b.id, chain_pair=a.id + ':' + b.id,
                   qc_flag='chain_pair_scope_not_curated_binding_partner_list')
        try:
            row.update(atom_contacts(trees[a.id], trees[b.id]))
        except Exception as e:
            row['count_error'] = str(e)
        rows.append(finish_record(row))
    return rows, chain_inventory


def safes():
    inventory = []
    for kind, suffix in [('experimental', '.pdb'), ('AF3', '_af3.pdb')]:
        rows = []
        paths = sorted((ROOT / 'data/SAFES_PDBs').glob('*' + suffix))
        if kind == 'experimental':
            paths = [p for p in paths if not p.stem.endswith('_af3')]
        for path in paths:
            new, chains = all_pairs(path, 'SAFES', kind, path.stem, dict(pdb_id=path.stem.split('_')[0]))
            rows.extend(new)
            inventory.extend(chains)
        write_rows('safes_' + kind.lower(), rows)
    pd.DataFrame(inventory).to_csv(OUT / 'safes_chain_inventory.csv', index=False)


def references():
    manifest = pd.read_csv(ROOT / 'results/rcsb_reference/rcsb_reference_assemblies.csv')
    manifest = manifest[manifest.scope == 'M_genitalium']
    rows, chains, assemblies = [], [], []
    for _, source in manifest.iterrows():
        target = source.pdb_id.lower() + '-assembly' + str(source.assembly_id)
        path = ROOT / 'data/mgen_reference_atom8A' / (target + '.cif.gz')
        meta = dict(pdb_id=source.pdb_id, assembly_id=source.assembly_id, source_candidate_identifier=source.rcsb_assembly_identifier)
        try:
            cif = MMCIF2Dict(io.StringIO(gzip.open(path, 'rt').read()))
            taxa = sorted(set(cif.get('_entity_src_gen.pdbx_gene_src_scientific_name', []) + cif.get('_entity_src_nat.pdbx_organism_scientific', [])))
            meta['source_organisms'] = ' | '.join(taxa)
            meta['experimental_method'] = ' | '.join(cif.get('_exptl.method', []))
            chain_to_entity = dict(zip(cif['_atom_site.label_asym_id'], cif['_atom_site.label_entity_id']))
        except Exception:
            chain_to_entity = {}
        new, inventory = all_pairs(path, 'Mgen_reference', 'experimental', target, meta)
        for row in new:
            row['chain1_entity_id'] = chain_to_entity.get(row.get('chain1'), '')
            row['chain2_entity_id'] = chain_to_entity.get(row.get('chain2'), '')
            if ' | ' in meta.get('source_organisms', ''):
                row['qc_flag'] += ';multiple_source_organisms_in_candidate_assembly'
        rows.extend(new)
        chains.extend(inventory)
        assemblies.append(dict(**meta, structure_id=target, n_protein_or_peptide_chains=len(inventory),
                               n_enumerated_pairs=len(new), n_calculable=sum(r['raw_included'] for r in new)))
    write_rows('mgen_reference_experimental', rows)
    pd.DataFrame(chains).to_csv(OUT / 'reference_chain_inventory.csv', index=False)
    pd.DataFrame(assemblies).to_csv(OUT / 'reference_assembly_inventory.csv', index=False)


def foldbench():
    original = pd.read_csv(ROOT / 'results/foldbench/paired_interfaces_8A.csv')
    output = {'experimental': [], 'AF3': []}
    for i, (target, group) in enumerate(original.groupby('target_id', sort=False), 1):
        for kind, prefix in [('experimental', 'experimental'), ('AF3', 'af3')]:
            path_col = 'experimental_path' if kind == 'experimental' else 'af3_path'
            value = group.iloc[0][path_col]
            model, trees, info, error = None, {}, {}, ''
            try:
                if not isinstance(value, str) or not value:
                    raise FileNotFoundError('no_AF3_for_this_FoldBench_target')
                model, info = parse_structure(ROOT / value)
                for c in protein_chains(model):
                    trees[c.id] = atom_tree(c)
            except Exception as e:
                error = str(e)
            for _, old in group.iterrows():
                a, b = (old.chain1, old.chain2) if kind == 'experimental' else (old.af3_chain1, old.af3_chain2)
                row = dict(dataset='FoldBench', structure_type=kind, structure_id=target,
                           interface_id=old.interface_id, source_path=value if isinstance(value, str) else '',
                           chain1=a, chain2=b, chain_pair=f'{a}:{b}',
                           benchmark_chain1=old.chain1, benchmark_chain2=old.chain2,
                           pair_scope='official_FoldBench_PPI_interface_row',
                           legacy_residue_contact_count=old.get(prefix + '_contact_count'),
                           legacy_qc_flags=old.qc_flags,
                           selected_seed=old.af3_selected_seed if kind == 'AF3' else None,
                           selected_sample=old.af3_selected_sample if kind == 'AF3' else None,
                           count_error=error, **info)
                if model is not None:
                    try:
                        row.update(atom_contacts(trees[a], trees[b]))
                        if kind == 'AF3' and isinstance(old.get('tied_interface_chain_pairs_json'), str):
                            alternatives = json.loads(old.tied_interface_chain_pairs_json)
                            counts = [atom_contacts(trees[x], trees[y])['contact_count'] for x, y in alternatives]
                            row['mapping_tie_atom_count_min'], row['mapping_tie_atom_count_max'] = min(counts), max(counts)
                            row['qc_flag'] = 'legacy_geometric_chain_mapping_tie'
                            if min(counts) != max(counts):
                                row['qc_flag'] += ';chain_mapping_tie_changes_atom_count'
                    except Exception as e:
                        row['count_error'] = str(e)
                output[kind].append(finish_record(row))
        if i % 40 == 0:
            print('FoldBench', i, 'assemblies', flush=True)
    for kind, rows in output.items():
        write_rows('foldbench_' + kind.lower(), rows)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    safes()
    references()
    from annotate_reference_atom8A import main as annotate_sources
    annotate_sources()
    foldbench()
