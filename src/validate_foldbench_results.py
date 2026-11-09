"""Validate output identities, paired arithmetic, selection, and plot populations."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/foldbench'
OUT = ROOT / 'results/foldbench'


def main():
    definition = pd.read_csv(DATA / 'metadata/interface_protein_protein.csv')
    result = pd.read_csv(OUT / 'paired_interfaces_8A.csv')
    selection = pd.read_csv(OUT / 'af3_model_selection.csv')
    assert len(result) == len(definition)
    assert result.interface_id.tolist() == list(range(1, len(definition) + 1))
    for out, original in [('target_id', 'pdb_id'), ('chain1', 'interface_chain_id_1'), ('chain2', 'interface_chain_id_2')]:
        assert result[out].tolist() == definition[original].tolist(), out
    assert selection.target_id.is_unique
    for _, row in selection.iterrows():
        summary = pd.read_csv(DATA / 'af3/summaries' / (row.target_id + '_results_summary.csv'))
        assert row.ranking_score == summary.ranking_score.max()
        top = summary[summary.ranking_score == summary.ranking_score.max()]
        numbers = top['seed-sample'].str.extract(r'seed_(\d+)_sample_(\d+)').astype(int)
        expected = min(map(tuple, numbers.to_numpy()))
        assert (row.selected_seed, row.selected_sample) == expected
        assert (ROOT / row.selected_pdb_path).is_file()
    assert result.paired_success.equals(result.experimental_contact_count.notna() & result.af3_contact_count.notna())
    for label in ('experimental', 'af3'):
        sub = result[result[label + '_contact_count'].notna()]
        contacts = sub[label + '_contact_count']
        assert np.all(contacts == contacts.astype(int))
        ia, ib = [sub[f'{label}_interface_residues_chain{i}'] for i in (1, 2)]
        assert np.all(sub[label + '_total_interface_residues'] == ia + ib)
        assert np.all((contacts >= np.maximum(ia, ib)) & (contacts <= ia * ib))
        for i in (1, 2):
            length = sub[f'{label}_chain{i}_length']
            interface = sub[f'{label}_interface_residues_chain{i}']
            assert np.all((length > 0) & (interface >= 0) & (interface <= length))
            assert np.allclose(sub[f'{label}_fraction_interface_chain{i}'], interface / length)
        zero = contacts == 0
        flagged_zero = sub.qc_flags.fillna('').str.contains(label + '_zero_contacts')
        assert zero.equals(flagged_zero)
    paired = result[result.paired_success]
    for prop in ('contact_count', 'total_interface_residues'):
        assert np.allclose(paired['delta_' + prop], paired['af3_' + prop] - paired['experimental_' + prop])
    maps = pd.read_csv(OUT / 'chain_mapping.csv')
    for _, group in maps.groupby('target_id'):
        assert group.experimental_chain.is_unique and group.af3_chain.is_unique
    bins = pd.read_csv(OUT / 'histogram_bins.csv')
    for (prop, _), group in bins.groupby(['property', 'dataset']):
        assert group['count'].sum() == len(paired) * (2 if prop == 'fraction' else 1)
    coverage = pd.read_csv(OUT / 'target_coverage.csv')
    # This particular run resolves every available AF3 interface; missing records remain.
    assert len(paired) == coverage.af3_selected_present.sum()
    report = dict(status='passed', original_interfaces_preserved=len(result),
                  highest_score_models_verified=len(selection), paired_interfaces_verified=len(paired),
                  checks=['original interface identity and ordering', 'numeric ranking-score tie break',
                          'integer residue-pair counts and bounds', 'residue totals and fractions',
                          'zero-contact preservation and flags', 'paired delta arithmetic',
                          'injective chain mappings', 'identical paired plot populations',
                          'all archive-available interfaces paired'])
    (OUT / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
