"""One checkpointed Pedro ZIP per Slurm array task; no existing files are edited."""
import argparse
import io
import json
import os
import shutil
import tempfile
import time
import zipfile
from pathlib import Path

import pandas as pd

from atom_contact_common import atom_contacts, atom_tree, finish_record, parse_structure

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/contact_comparison_atom8A'


def run(index):
    spec = json.loads((OUT / 'pedro_manifest.json').read_text())[index]
    dest = OUT / 'pedro_parts' / spec['source_file']
    marker = dest.with_suffix('.done.json')
    if dest.exists() and marker.exists():
        print('Already completed:', dest, flush=True)
        return
    started = time.monotonic()
    observations = pd.read_csv(ROOT / 'results/pooled_parts_8A' / spec['source_file'])
    rows = []
    # TemporaryDirectory cleans only this task's own temporary archive on completion.
    scratch_parent = os.environ.get('SLURM_TMPDIR') or str(ROOT / 'tmp')
    with tempfile.TemporaryDirectory(prefix='pedro_atom8A_', dir=scratch_parent) as scratch:
        inner_path = Path(scratch) / 'input.zip'
        print('Copying', spec['member'], flush=True)
        with zipfile.ZipFile(ROOT / 'data/AF3_completed_zips.zip') as outer:
            with outer.open(spec['member']) as src, inner_path.open('wb') as dst:
                shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
        print('Inner archive ready', flush=True)
        with zipfile.ZipFile(inner_path) as inner:
            names = inner.namelist()
            for j, ((job, best), group) in enumerate(observations.groupby(['job_id', 'best_model'], sort=False), 1):
                suffix = f'_model_{int(best)}.cif'
                candidates = [n for n in names if n.startswith(str(job) + '/') and n.endswith(suffix)]
                model, trees, parse_info, error = None, {}, {}, ''
                try:
                    if len(candidates) != 1:
                        raise ValueError(f'expected_one_CIF_found_{len(candidates)}')
                    model, parse_info = parse_structure(cif_text=inner.read(candidates[0]).decode(), author_chains=True)
                    for chain in model:
                        try:
                            trees[chain.id] = atom_tree(chain)
                        except ValueError:
                            pass
                except Exception as exc:
                    error = str(exc)
                for _, old in group.iterrows():
                    row = dict(dataset='Pedro', structure_type='AF3',
                               structure_id=f'{Path(spec["source_file"]).stem}/{job}/model_{int(best)}',
                               source_file=spec['source_file'], archive_member=spec['member'],
                               cif_member=candidates[0] if len(candidates) == 1 else '',
                               job_id=job, protein_pair=old.protein_pair,
                               chain1=old.chain_a, chain2=old.chain_b,
                               chain_pair=str(old.chain_a) + ':' + str(old.chain_b),
                               best_model=int(best), legacy_residue_contact_count=old.contact_count,
                               legacy_min_distance=old.min_interchain_distance,
                               legacy_has_clash=old.has_clash,
                               count_error=error, **parse_info)
                    if model is not None:
                        try:
                            # Recount every pair, including legacy zeros; no zero shortcut.
                            row.update(atom_contacts(trees[old.chain_a], trees[old.chain_b]))
                            if old.has_clash == 1:
                                row['qc_flag'] = 'legacy_AF3_has_clash_annotation'
                            if (old.contact_count == 0) != (row['contact_count'] == 0):
                                row['qc_flag'] = row.get('qc_flag', '') + ';legacy_zero_status_changed'
                        except Exception as exc:
                            row['count_error'] = str(exc)
                    rows.append(finish_record(row, inherited_severe=old.min_interchain_distance < 0.5))
                if j % 5 == 0:
                    print(f'{j} jobs; {len(rows)} chain-pair observations', flush=True)
    frame = pd.DataFrame(rows)
    assert len(frame) == len(observations)
    partial = dest.with_suffix('.partial.csv')
    frame.to_csv(partial, index=False)
    partial.replace(dest)
    report = dict(source_file=spec['source_file'], n_rows=len(frame), n_calculable=int(frame.raw_included.sum()),
                  n_legacy_severe=int((observations.min_interchain_distance < 0.5).sum()),
                  elapsed_seconds=time.monotonic() - started)
    marker.write_text(json.dumps(report, indent=2) + '\n')
    print(report, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index', required=True, type=int)
    run(p.parse_args().index)
