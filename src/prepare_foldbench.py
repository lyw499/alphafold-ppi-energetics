"""Inventory official FoldBench archives and extract only required structures.

No existing Pedro/SAFES outputs are touched. Run from any working directory.
"""
import argparse
import csv
import hashlib
import io
import json
import re
import shutil
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/foldbench'
OUT = ROOT / 'results/foldbench'


def write_csv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def copy_member(archive, member, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    # Never use extractall: names in third-party archives are not trusted paths.
    with archive.extractfile(member) as src, path.open('wb') as dst:
        shutil.copyfileobj(src, dst)


def experimental(targets):
    found = set()
    rows = []
    archive = DATA / 'experimental/ground_truth_1522.tar'
    with tarfile.open(archive, 'r:') as t:
        for m in t:
            if not m.isfile():
                continue
            rows.append({'member': m.name, 'bytes': m.size})
            target = Path(m.name).stem
            if target in targets and m.name.endswith('.cif'):
                if target in found:
                    raise ValueError(f'Duplicate experimental target: {target}')
                copy_member(t, m, DATA / 'experimental/structures' / (target + '.cif'))
                found.add(target)
    write_csv(OUT / 'experimental_archive_inventory.csv', rows)
    print(f'Experimental structures: {len(found)}/{len(targets)}', flush=True)
    return found


def select_summary(text, target):
    rows = list(csv.DictReader(io.StringIO(text)))
    candidates = []
    for row in rows:
        match = re.fullmatch(r'seed_(\d+)_sample_(\d+)', row['seed-sample'])
        if not match or row['target'] != target:
            raise ValueError(f'Unexpected summary identity: {target}: {row}')
        score = float(row['ranking_score'])
        if not (-float('inf') < score < float('inf')):
            raise ValueError(f'Nonfinite ranking score: {target}')
        candidates.append((score, int(match[1]), int(match[2]), row))
    # Scores are rounded by the depositor: record ties; deterministic numeric tie break.
    score, seed, sample, row = sorted(candidates, key=lambda x: (-x[0], x[1], x[2]))[0]
    return dict(target_id=target, selected_seed=seed, selected_sample=sample,
                ranking_score=score, ptm=float(row['ptm']), iptm=float(row['iptm']),
                n_predictions=len(rows), n_top_score_ties=sum(x[0] == score for x in candidates),
                tie_break='lowest_numeric_seed_then_sample',
                selected_pdb_path=str((DATA / 'af3/selected' / f'{target}_seed_{seed}_sample_{sample}.pdb').relative_to(ROOT)))


def af3(targets):
    path = DATA / 'af3/af3_foldbench.tar.gz'
    meta = json.loads((DATA / 'metadata/zenodo_record.json').read_text())
    spec = next(f for f in meta['files'] if f['key'] == path.name)
    if path.stat().st_size != spec['size']:
        raise ValueError('AF3 archive is incomplete; wait for the download to finish.')
    digest = hashlib.md5()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            digest.update(block)
    if 'md5:' + digest.hexdigest() != spec['checksum']:
        raise ValueError('Zenodo checksum mismatch')
    print('AF3 archive MD5 verified', flush=True)
    summaries, selected, extracted, pdb_counts = {}, {}, set(), {}
    required_members = {}
    with (OUT / 'af3_archive_inventory.csv').open('w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['member', 'bytes'])
        with tarfile.open(path, 'r|gz') as t:
            for i, m in enumerate(t):
                if not m.isfile():
                    continue
                writer.writerow([m.name, m.size])
                name = Path(m.name).name
                if name in ('MODIFICATIONS.md', 'OUTPUT_TERMS_OF_USE.md', 'Legally_Binding_Terms_of_Use.txt'):
                    copy_member(t, m, DATA / 'af3' / name)
                if name.endswith('_results_summary.csv'):
                    target = name.removesuffix('_results_summary.csv')
                    text = t.extractfile(m).read().decode()
                    if target in summaries:
                        raise ValueError(f'Duplicate AF3 summary: {target}')
                    summaries[target] = m.name
                    summary_path = DATA / 'af3/summaries' / name
                    summary_path.parent.mkdir(parents=True, exist_ok=True)
                    summary_path.write_text(text)
                    if target in targets:
                        row = select_summary(text, target)
                        wanted = str(Path(m.name).parent / Path(row['selected_pdb_path']).name)
                        row['archive_member'] = wanted
                        selected[target] = row
                        required_members[wanted] = target
                if name.endswith('.pdb'):
                    target = name.split('_seed_')[0]
                    pdb_counts[target] = pdb_counts.get(target, 0) + 1
                if m.name in required_members:
                    target = required_members[m.name]
                    copy_member(t, m, ROOT / selected[target]['selected_pdb_path'])
                    extracted.add(target)
                if i % 10000 == 0:
                    print(f'AF3 archive: {i} members, {len(summaries)} summaries, {len(extracted)} selected', flush=True)
    # General fallback for archives whose summaries follow their PDBs.
    remaining = {r['archive_member']: target for target, r in selected.items() if target not in extracted}
    if remaining:
        with tarfile.open(path, 'r|gz') as t:
            for m in t:
                if m.name in remaining:
                    target = remaining.pop(m.name)
                    copy_member(t, m, ROOT / selected[target]['selected_pdb_path'])
                    extracted.add(target)
                if not remaining:
                    break
    for target, row in selected.items():
        row['selected_pdb_present'] = target in extracted
        row['archive_pdb_count'] = pdb_counts.get(target, 0)
    write_csv(OUT / 'af3_model_selection.csv', [selected[t] for t in sorted(selected)])
    write_csv(OUT / 'af3_target_inventory.csv', [dict(target_id=t, summary_member=summaries[t], n_pdb_files=pdb_counts.get(t, 0)) for t in sorted(summaries)])
    return summaries, extracted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experimental-only', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((DATA / 'metadata/interface_protein_protein.csv').open()))
    targets = {r['pdb_id'] for r in rows}
    exp = experimental(targets)
    if args.experimental_only:
        return
    inventory, present = af3(targets)
    coverage = []
    for i, r in enumerate(rows, 1):
        target = r['pdb_id']
        coverage.append(dict(interface_id=i, target_id=target,
                             chain1=r['interface_chain_id_1'], chain2=r['interface_chain_id_2'],
                             experimental_present=target in exp, af3_summary_present=target in inventory,
                             af3_selected_present=target in present))
    write_csv(OUT / 'target_coverage.csv', coverage)
    summary = dict(n_ppi_interface_rows=len(rows), n_unique_ppi_assemblies=len(targets),
                   n_af3_archive_targets=len(inventory), n_ppi_assemblies_in_af3=len(targets & set(inventory)),
                   n_ppi_assemblies_with_selected_pdb=len(present),
                   n_ppi_rows_with_selected_pdb=sum(r['pdb_id'] in present for r in rows),
                   missing_af3_targets=sorted(targets - present), missing_experimental_targets=sorted(targets - exp),
                   af3_archive_bytes=path_size(DATA / 'af3/af3_foldbench.tar.gz'))
    (OUT / 'coverage_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)


def path_size(path):
    return path.stat().st_size


if __name__ == '__main__':
    main()
