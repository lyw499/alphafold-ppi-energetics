"""FoldBench paired analysis using the unchanged Pedro 8-A residue-pair function.

Experimental chains use label_asym_id, as in the official FoldBench parser.
Sequence-compatible chain assignments are resolved by whole-complex C-alpha RMSD,
never by contact count, interface size, or an interface-specific score.
"""
import argparse
import hashlib
import json
import os
import platform
from collections import Counter
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/foldbench-matplotlib')
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1

# Reuse the exact existing implementation; do not fork the contact definition.
from parse_pooled_inner_zip_8A import count_pair_contacts, get_residues

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/foldbench'
OUT = ROOT / 'results/foldbench'
MAX_MAPPINGS = 50000
MIN_IDENTITY = 0.95
MIN_OBSERVED_COVERAGE = 0.90
RMSD_TIE_TOLERANCE = 1e-4


def protein_residues(chain):
    """Mapping uses the same coordinate-residue selection as contact counting."""
    residues = get_residues(chain)
    sequence = ''.join(seq1(r.resname, undef_code='X') for r in residues)
    if not sequence or sequence.count('X') / len(sequence) > 0.2:
        return [], ''
    return residues, sequence


def load_model(path, experimental=False):
    # Keep Bio.PDB's default auth residue numbering, altloc selection, and model 0.
    parser = MMCIFParser(QUIET=True, auth_chains=False) if experimental else PDBParser(QUIET=True)
    structure = parser.get_structure('structure', str(path))
    return structure[0], len(structure)


def alignment_info(exp_res, exp_seq, pred_res, pred_seq):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    # Missing terminal residues in experimental coordinates are not penalized.
    aligner.end_gap_score = 0
    aln = aligner.align(exp_seq, pred_seq)[0]
    pairs = []
    matched = 0
    aligned = 0
    for a, b in zip(*aln.aligned):
        for i, j in zip(range(*a), range(*b)):
            aligned += 1
            if exp_seq[i] == pred_seq[j] and exp_seq[i] != 'X':
                matched += 1
                if 'CA' in exp_res[i] and 'CA' in pred_res[j]:
                    pairs.append((exp_res[i]['CA'].coord, pred_res[j]['CA'].coord))
    identity = matched / aligned if aligned else 0
    coverage = aligned / len(exp_seq) if exp_seq else 0
    if identity < MIN_IDENTITY or coverage < MIN_OBSERVED_COVERAGE or len(pairs) < 3:
        return None
    return dict(identity=identity, experimental_observed_coverage=coverage,
                af3_observed_coverage=aligned / len(pred_seq),
                n_ca_pairs=len(pairs), xyz=np.asarray(pairs, dtype=float))


def fitted_rmsd(x, y):
    x = x - x.mean(axis=0)
    y = y - y.mean(axis=0)
    u, _, vt = np.linalg.svd(y.T @ x)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    return float(np.sqrt(np.mean(np.sum((y @ rotation - x) ** 2, axis=1))))


def map_chains(experimental, prediction, required):
    exp = {c.id: protein_residues(c) for c in experimental}
    pred = {c.id: protein_residues(c) for c in prediction}
    exp = {k: v for k, v in exp.items() if v[1]}
    pred = {k: v for k, v in pred.items() if v[1]}
    absent = set(required) - set(exp)
    if absent:
        raise ValueError('missing_experimental_protein_chains:' + ','.join(sorted(absent)))
    compatible = {}
    for a, (ra, sa) in exp.items():
        compatible[a] = {}
        for b, (rb, sb) in pred.items():
            # Experimental coordinates can have gaps but cannot greatly exceed AF3.
            if len(sa) > len(sb) / MIN_OBSERVED_COVERAGE:
                continue
            info = alignment_info(ra, sa, rb, sb)
            if info is not None:
                compatible[a][b] = info
    missing = [a for a in required if not compatible.get(a)]
    if missing:
        raise ValueError('no_sequence_compatible_af3_chain:' + ','.join(sorted(missing)))
    # Use all matchable protein chains for a single assembly-wide assignment.
    order = sorted((a for a in exp if compatible[a]), key=lambda a: (len(compatible[a]), a))
    if len(order) > len(pred):
        raise ValueError('fewer_af3_protein_chains_than_matchable_experimental_chains')
    candidates = []

    def visit(index, mapping, used):
        if len(candidates) > MAX_MAPPINGS:
            raise ValueError('sequence_mapping_search_exceeds_limit')
        if index == len(order):
            xyz = np.concatenate([compatible[a][b]['xyz'] for a, b in mapping.items()])
            rmsd = fitted_rmsd(xyz[:, 0], xyz[:, 1])
            candidates.append((rmsd, dict(mapping)))
            return
        a = order[index]
        for b in sorted(compatible[a]):
            if b not in used:
                mapping[a] = b
                visit(index + 1, mapping, used | {b})
                mapping.pop(a)

    visit(0, {}, set())
    if not candidates:
        raise ValueError('no_injective_sequence_compatible_chain_mapping')
    candidates.sort(key=lambda item: (item[0], tuple(sorted(item[1].items()))))
    best_rmsd = candidates[0][0]
    ties = [item for item in candidates if item[0] <= best_rmsd + RMSD_TIE_TOLERANCE]
    # For geometric ties, prefer matching labels, then lexicographic ordering.
    _, best = min(ties, key=lambda item: (-sum(a == b for a, b in item[1].items()), tuple(sorted(item[1].items()))))
    details = []
    for a, b in sorted(best.items()):
        info = {k: v for k, v in compatible[a][b].items() if k != 'xyz'}
        details.append(dict(experimental_chain=a, af3_chain=b, **info,
                            same_chain_id=a == b, n_sequence_candidates=len(compatible[a]),
                            n_assembly_mappings=len(candidates), n_rmsd_tied_mappings=len(ties),
                            assembly_ca_rmsd_A=best_rmsd,
                            method='sequence_unique' if len(candidates) == 1 else 'sequence_plus_global_CA_RMSD'))
    return best, details, [mapping for _, mapping in ties], best_rmsd


def metrics(chain1, chain2):
    m = count_pair_contacts(chain1, chain2, cutoff=8.0)
    return dict(contact_count=m['contact_count'],
                interface_residues_chain1=m['interface_residues_a'],
                interface_residues_chain2=m['interface_residues_b'],
                total_interface_residues=m['interface_residues_a'] + m['interface_residues_b'],
                chain1_length=m['length_a'], chain2_length=m['length_b'],
                fraction_interface_chain1=m['interface_residues_a'] / m['length_a'] if m['length_a'] else np.nan,
                fraction_interface_chain2=m['interface_residues_b'] / m['length_b'] if m['length_b'] else np.nan,
                interface_residue_ids_chain1=m['interface_residue_ids_a'],
                interface_residue_ids_chain2=m['interface_residue_ids_b'],
                min_interchain_distance_A=m['min_interchain_distance'])


def analysis(limit=None):
    definitions = pd.read_csv(DATA / 'metadata/interface_protein_protein.csv')
    definitions.insert(0, 'interface_id', np.arange(1, len(definitions) + 1))
    selections = pd.read_csv(OUT / 'af3_model_selection.csv').set_index('target_id')
    rows, chain_rows = [], []
    groups = list(definitions.groupby('pdb_id', sort=False))
    if limit:
        groups = groups[:limit]
    for idx, (target, group) in enumerate(groups, 1):
        exp_path = DATA / 'experimental/structures' / (target + '.cif')
        selected = selections.loc[target].to_dict() if target in selections.index else {}
        pred_path = ROOT / selected['selected_pdb_path'] if selected else None
        exp, pred, mapping, details, tied, rmsd = None, None, {}, [], [], np.nan
        exp_error, pred_error, mapping_error = '', '', ''
        exp_models = pred_models = 0
        try:
            if not exp_path.exists():
                raise FileNotFoundError('missing_experimental_structure')
            exp, exp_models = load_model(exp_path, experimental=True)
        except Exception as e:
            exp_error = str(e)
        try:
            if pred_path is None or not pred_path.exists():
                raise FileNotFoundError('missing_af3_selected_structure')
            pred, pred_models = load_model(pred_path)
        except Exception as e:
            pred_error = str(e)
        if exp is not None and pred is not None:
            try:
                required = set(group.interface_chain_id_1) | set(group.interface_chain_id_2)
                mapping, details, tied, rmsd = map_chains(exp, pred, required)
                chain_rows.extend(dict(target_id=target, **d) for d in details)
            except Exception as e:
                mapping_error = str(e)
        for _, definition in group.iterrows():
            a, b = definition.interface_chain_id_1, definition.interface_chain_id_2
            flags = []
            row = dict(interface_id=definition.interface_id, target_id=target, chain1=a, chain2=b,
                       contact_definition='unique_residue_pairs_any_nonH_atoms_le_8A', cutoff_A=8.0,
                       experimental_path=str(exp_path.relative_to(ROOT)),
                       af3_path=str(pred_path.relative_to(ROOT)) if pred_path else '',
                       af3_chain1=mapping.get(a, ''), af3_chain2=mapping.get(b, ''),
                       af3_selected_seed=selected.get('selected_seed'), af3_selected_sample=selected.get('selected_sample'),
                       af3_ranking_score=selected.get('ranking_score'), af3_ptm=selected.get('ptm'), af3_iptm=selected.get('iptm'),
                       af3_n_top_score_ties=selected.get('n_top_score_ties'),
                       assembly_ca_rmsd_A=rmsd, n_rmsd_tied_mappings=len(tied),
                       experimental_n_models=exp_models, af3_n_models=pred_models)
            if exp_error:
                flags.append('experimental_error:' + exp_error)
            if pred_error:
                flags.append('af3_error:' + pred_error)
            if mapping_error:
                flags.append('mapping_error:' + mapping_error)
            for label, model, x, y in [('experimental', exp, a, b), ('af3', pred, mapping.get(a), mapping.get(b))]:
                if model is None:
                    continue
                if x not in model or y not in model:
                    flags.append(label + '_missing_or_unmapped_chain')
                    continue
                try:
                    values = metrics(model[x], model[y])
                    row.update({label + '_' + k: v for k, v in values.items()})
                    if values['contact_count'] == 0:
                        flags.append(label + '_zero_contacts')
                    distance = values['min_interchain_distance_A']
                    if distance is not None and distance < 0.5:
                        flags.append(label + '_severe_overlap_lt_0.5A')
                    if min(values['chain1_length'], values['chain2_length']) == 0:
                        flags.append(label + '_empty_chain')
                except Exception as e:
                    flags.append(label + '_contact_error:' + str(e))
            if mapping:
                if rmsd > 10:
                    flags.append('assembly_CA_RMSD_gt_10A')
                pair_maps = sorted({(t[a], t[b]) for t in tied})
                row['n_tied_interface_mappings'] = len(pair_maps)
                if len(pair_maps) > 1:
                    flags.append('geometrically_tied_chain_mapping')
                    alternatives = [metrics(pred[x], pred[y]) for x, y in pair_maps]
                    changed = False
                    for prop in ('contact_count', 'total_interface_residues', 'fraction_interface_chain1', 'fraction_interface_chain2'):
                        values = [m[prop] for m in alternatives]
                        row['mapping_tie_' + prop + '_min'] = min(values)
                        row['mapping_tie_' + prop + '_max'] = max(values)
                        changed |= max(values) != min(values)
                    if row['mapping_tie_contact_count_max'] != row['mapping_tie_contact_count_min']:
                        flags.append('mapping_tie_changes_contact_count')
                    if changed:
                        flags.append('mapping_tie_changes_interface_metrics')
                    row['tied_interface_chain_pairs_json'] = json.dumps(pair_maps)
                if a != mapping[a] or b != mapping[b]:
                    flags.append('chain_ids_remapped')
            paired = 'experimental_contact_count' in row and 'af3_contact_count' in row
            row['paired_success'] = paired
            if paired:
                for field in ('contact_count', 'total_interface_residues'):
                    row['delta_' + field] = row['af3_' + field] - row['experimental_' + field]
                if any(row['experimental_chain' + str(c) + '_length'] != row['af3_chain' + str(c) + '_length'] for c in (1, 2)):
                    flags.append('coordinate_chain_lengths_differ')
            if selected.get('n_top_score_ties', 0) > 1:
                flags.append('rounded_ranking_score_tie')
            if exp_models > 1 or pred_models > 1:
                flags.append('multiple_coordinate_models_first_used')
            row['qc_flags'] = ';'.join(flags)
            rows.append(row)
        if idx % 10 == 0 or idx == len(groups):
            print(f'{idx}/{len(groups)} assemblies; {len(rows)} interfaces; {sum(r["paired_success"] for r in rows)} paired', flush=True)
    df = pd.DataFrame(rows).sort_values('interface_id')
    # Flag unusually large counts against one common threshold; retain every row.
    all_counts = pd.concat([df.get('experimental_contact_count', pd.Series(dtype=float)), df.get('af3_contact_count', pd.Series(dtype=float))]).dropna()
    q1, q3 = all_counts.quantile([0.25, 0.75])
    extreme_threshold = float(q3 + 3 * (q3 - q1))
    for label in ('experimental', 'af3'):
        col = label + '_contact_count'
        if col in df:
            mask = df[col] > extreme_threshold
            df.loc[mask, 'qc_flags'] = df.loc[mask, 'qc_flags'].apply(lambda s: (s + ';' if s else '') + label + '_extreme_contact_count')
    df['extreme_contact_threshold'] = extreme_threshold
    prefix = 'pilot_' if limit else ''
    df.to_csv(OUT / (prefix + 'paired_interfaces_8A.csv'), index=False)
    pd.DataFrame(chain_rows).to_csv(OUT / (prefix + 'chain_mapping.csv'), index=False)
    qc_cols = [c for c in df.columns if c in ['interface_id', 'target_id', 'chain1', 'chain2', 'af3_chain1', 'af3_chain2', 'paired_success', 'qc_flags', 'assembly_ca_rmsd_A', 'n_rmsd_tied_mappings', 'n_tied_interface_mappings', 'experimental_n_models', 'af3_n_models', 'af3_n_top_score_ties', 'extreme_contact_threshold', 'tied_interface_chain_pairs_json'] or c.startswith('mapping_tie_') or c.endswith('min_interchain_distance_A') or c.endswith('_contact_count')]
    df[qc_cols].to_csv(OUT / (prefix + 'qc_interfaces.csv'), index=False)
    return df


def summarize_and_plot(df):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    paired = df[df.paired_success].copy()
    if paired.empty:
        raise ValueError('No paired interfaces available for summary/plots')
    summary = []
    def add_summary(dataset, prop, values):
        values = pd.Series(values).dropna()
        summary.append(dict(dataset=dataset, property=prop, n=len(values), median=values.median(),
                            q1=values.quantile(.25), q3=values.quantile(.75),
                            iqr=values.quantile(.75) - values.quantile(.25)))
    fractions = {}
    for label in ('experimental', 'af3'):
        for prop in ('contact_count', 'total_interface_residues'):
            add_summary(label, prop, paired[label + '_' + prop])
        # Both specified chains contribute once per benchmark interface.
        fractions[label] = np.concatenate([paired[label + '_fraction_interface_chain1'], paired[label + '_fraction_interface_chain2']])
        add_summary(label, 'interface_fraction_pooled_chains', fractions[label])
    for prop in ('contact_count', 'total_interface_residues'):
        add_summary('af3_minus_experimental', prop, paired['delta_' + prop])
    pd.DataFrame(summary).to_csv(OUT / 'distribution_summary.csv', index=False)
    primary_summary = summary
    summary = []
    stable = paired[~paired.qc_flags.fillna('').str.contains('mapping_tie_changes_interface_metrics')]
    for label in ('experimental', 'af3'):
        for prop in ('contact_count', 'total_interface_residues'):
            add_summary(label, prop, stable[label + '_' + prop])
        add_summary(label, 'interface_fraction_pooled_chains', np.concatenate([
            stable[label + '_fraction_interface_chain1'], stable[label + '_fraction_interface_chain2']]))
    for prop in ('contact_count', 'total_interface_residues'):
        add_summary('af3_minus_experimental', prop, stable['delta_' + prop])
    pd.DataFrame(summary).to_csv(OUT / 'mapping_sensitivity_summary.csv', index=False)
    summary = primary_summary
    missing = df[df.af3_path.fillna('') == ''].groupby('target_id').size().reset_index(name='n_missing_interface_rows')
    missing.to_csv(OUT / 'missing_af3_targets.csv', index=False)
    df[df.qc_flags.fillna('').str.contains('mapping_tie_changes_interface_metrics')].to_csv(OUT / 'mapping_sensitive_interfaces.csv', index=False)
    # Long table makes pooled-chain weighting explicit and preserves the pair identity.
    long_rows = []
    for _, row in paired.iterrows():
        for dataset in ('experimental', 'af3'):
            for chain in (1, 2):
                long_rows.append(dict(interface_id=row.interface_id, target_id=row.target_id, dataset=dataset,
                                      benchmark_chain=row['chain' + str(chain)], chain_role='chain' + str(chain),
                                      chain_length=row[f'{dataset}_chain{chain}_length'],
                                      interface_residues=row[f'{dataset}_interface_residues_chain{chain}'],
                                      interface_fraction=row[f'{dataset}_fraction_interface_chain{chain}']))
    pd.DataFrame(long_rows).to_csv(OUT / 'chain_level_fractions.csv', index=False)
    flag_counts = Counter(flag for text in df.qc_flags.fillna('') for flag in text.split(';') if flag)
    report = dict(n_interfaces=len(df), n_successfully_paired=len(paired),
                  n_pairs_without_tie_sensitive_interface_metrics=len(stable),
                  n_unpaired=len(df) - len(paired), n_paired_assemblies=paired.target_id.nunique(),
                  zero_contact_interfaces={label: int((df[label + '_contact_count'] == 0).sum()) for label in ('experimental', 'af3')},
                  zero_contact_interfaces_paired={label: int((paired[label + '_contact_count'] == 0).sum()) for label in ('experimental', 'af3')},
                  qc_flag_counts=dict(flag_counts),
                  unpaired_targets=sorted(df.loc[~df.paired_success, 'target_id'].unique().tolist()),
                  note='Descriptive statistics use the same successfully paired interfaces; zero and flagged contacts retained. Fractions pool both chains per interface; lengths count observed coordinate residues, not full sequences.')
    (OUT / 'analysis_summary.json').write_text(json.dumps(report, indent=2) + '\n')
    figdir = OUT / 'figures'
    figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 180})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    bin_rows = []
    labels = ['Residue-pair contacts (8 Å)', 'Total interface residues', 'Interface fraction (both chains pooled)']
    for ax, prop, title in zip(axes, ('contact_count', 'total_interface_residues', 'fraction'), labels):
        values = fractions if prop == 'fraction' else {s: paired[s + '_' + prop].to_numpy() for s in ('experimental', 'af3')}
        combined = np.concatenate(list(values.values()))
        bins = np.linspace(0, 1, 21) if prop == 'fraction' else np.linspace(0, max(float(combined.max()) + 1, 1), 31)
        for label, color in [('experimental', '#2878B5'), ('af3', '#D66B27')]:
            counts, _ = np.histogram(values[label], bins=bins)
            ax.hist(values[label], bins=bins, weights=np.full(len(values[label]), 100 / len(values[label])),
                    histtype='step', linewidth=2, color=color, label=label.upper() if label == 'af3' else 'Experimental')
            for left, right, count in zip(bins[:-1], bins[1:], counts):
                bin_rows.append(dict(property=prop, dataset=label, bin_left=left, bin_right=right, count=int(count)))
        ax.set_xlabel(title)
        ax.set_ylabel('Observations (%)')
        ax.set_xlim(bins[0], bins[-1])
        ax.legend(frameon=False)
    fig.suptitle(f'FoldBench paired interfaces (n={len(paired)}); zeros and QC-flagged pairs retained')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(figdir / ('interface_distributions_8A.' + ext))
    plt.close(fig)
    pd.DataFrame(bin_rows).to_csv(OUT / 'histogram_bins.csv', index=False)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, prop, title in zip(axes[:2], ('contact_count', 'total_interface_residues'), labels[:2]):
        x, y = paired['experimental_' + prop], paired['af3_' + prop]
        high = max(x.max(), y.max()) * 1.04
        ax.scatter(x, y, s=19, alpha=.6, color='#2878B5')
        ax.plot([0, high], [0, high], '--', color='0.4', label='y = x')
        ax.set(xlim=(0, high), ylim=(0, high), xlabel='Experimental: ' + title, ylabel='AF3: ' + title)
        ax.set_aspect('equal', adjustable='box')
        ax.legend(frameon=False)
    delta = paired.delta_contact_count
    axes[2].hist(delta, bins=30, color='#6B8095', edgecolor='white')
    axes[2].axvline(0, color='0.3', linestyle='--', label='Zero difference')
    axes[2].axvline(delta.median(), color='#D66B27', label=f'Median = {delta.median():g}')
    axes[2].set(xlabel='AF3 − experimental residue-pair contacts', ylabel='Interfaces')
    axes[2].legend(frameon=False)
    fig.suptitle(f'FoldBench paired diagnostics (n={len(paired)})')
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(figdir / ('paired_diagnostics_8A.' + ext))
    plt.close(fig)
    print(json.dumps(report, indent=2), flush=True)


def write_provenance():
    from importlib.metadata import version
    paths = [ROOT / 'src' / name for name in (
        'prepare_foldbench.py', 'foldbench_interfaces.py', 'test_foldbench_interfaces.py',
        'validate_foldbench_results.py', 'parse_pooled_inner_zip_8A.py', 'build_pair_level_contacts_8A.py')]
    paths += list((DATA / 'metadata').glob('*'))
    paths += list((DATA / 'experimental/structures').glob('*.cif'))
    paths += list((DATA / 'af3/selected').glob('*.pdb'))
    hashes = {}
    for path in sorted(paths):
        if path.is_file():
            hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    report = dict(python=platform.python_version(),
                  packages={name: version(name) for name in ('biopython', 'numpy', 'scipy', 'pandas', 'matplotlib')},
                  input_and_code_sha256=hashes,
                  contact_function='parse_pooled_inner_zip_8A.count_pair_contacts',
                  cutoff_A=8.0, min_mapping_identity=MIN_IDENTITY,
                  min_mapping_experimental_observed_coverage=MIN_OBSERVED_COVERAGE,
                  max_chain_mappings=MAX_MAPPINGS, rmsd_tie_tolerance_A=RMSD_TIE_TOLERANCE)
    (OUT / 'provenance.json').write_text(json.dumps(report, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, help='Pilot only; separate pilot output names, no final plots')
    parser.add_argument('--plots-only', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.plots_only:
        df = pd.read_csv(OUT / 'paired_interfaces_8A.csv')
    else:
        df = analysis(args.limit)
    if not args.limit:
        summarize_and_plot(df)
        write_provenance()


if __name__ == '__main__':
    main()
