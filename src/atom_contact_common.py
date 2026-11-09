"""8-A inter-chain atom-pair counts; separate from legacy residue-pair outputs."""
import gzip
import io
import json
import re
import warnings
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, is_aa
from scipy.spatial import cKDTree

from parse_pooled_inner_zip_8A import get_residues

CUTOFF = 8.0
DEFINITION = 'interchain_nonH_nonD_atom_pairs_distance_le_8A'


def parse_structure(path=None, cif_text=None, author_chains=False):
    """Keep warnings as evidence; do not discard calculable malformed structures."""
    if cif_text is not None:
        text = cif_text
        handle = io.StringIO(text)
        is_cif = True
    else:
        path = Path(path)
        text = gzip.open(path, 'rt').read() if path.suffix == '.gz' else path.read_text()
        handle = io.StringIO(text)
        is_cif = '.cif' in path.name
    parser = MMCIFParser(QUIET=False, auth_chains=author_chains) if is_cif else PDBParser(QUIET=False)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        structure = parser.get_structure('structure', handle)
    messages = sorted({str(w.message) for w in caught})
    # A warning alone is not proof of malformed eligible protein atoms: duplicated
    # ions and legal alternative residue conformers are irrelevant to our count.
    seen, malformed_chains = set(), set()
    if is_cif:
        d = parser._mmcif_dict
        n = len(d['_atom_site.id'])
        chains = d['_atom_site.auth_asym_id' if author_chains else '_atom_site.label_asym_id']
        models = d.get('_atom_site.pdbx_PDB_model_num', ['1'] * n)
        icodes = d.get('_atom_site.pdbx_PDB_ins_code', ['?'] * n)
        atoms = zip(d['_atom_site.group_PDB'], d['_atom_site.type_symbol'], models, chains,
                    d.get('_atom_site.auth_seq_id', d['_atom_site.label_seq_id']), icodes,
                    d['_atom_site.label_atom_id'], d['_atom_site.label_alt_id'])
    else:
        atoms = []
        current_model = '1'
        for line in text.splitlines():
            if line.startswith('MODEL '):
                current_model = line[10:14].strip()
            if line.startswith('ATOM  '):
                atoms.append(('ATOM', line[76:78].strip(), current_model, line[21:22],
                              line[22:26].strip(), line[26:27], line[12:16].strip(), line[16:17]))
    first_model = None
    for group, element, model_number, chain_id, number, icode, name, altloc in atoms:
        if group != 'ATOM' or element in ('H', 'D'):
            continue
        if first_model is None:
            first_model = model_number
        if model_number != first_model:
            continue
        key = (chain_id, number, icode if icode not in ('.', '?') else ' ', name,
               altloc if altloc not in ('.', '?') else ' ')
        if key in seen:
            malformed_chains.add(chain_id)
        seen.add(key)
    info = dict(parser_warnings=' | '.join(messages), malformed_duplicate_identity=bool(malformed_chains),
                malformed_chain_ids=json.dumps(sorted(malformed_chains)),
                n_coordinate_models=len(structure))
    return structure[0], info


def protein_chains(model):
    # Include short peptides; no arbitrary chain-length cutoff.
    return [c for c in model if any(is_aa(r, standard=False) for r in get_residues(c))]


def atom_tree(chain):
    residues = get_residues(chain)
    atoms = [a for r in residues for a in r if a.element not in ('H', 'D')]
    coords = np.asarray([a.coord for a in atoms], dtype=float).reshape((-1, 3))
    if not np.isfinite(coords).all():
        raise ValueError('nonfinite_atom_coordinates')
    if not len(coords):
        raise ValueError('no_eligible_heavy_atoms_in_required_chain')
    return cKDTree(coords), len(residues)


def atom_contacts(first, second):
    a, la = first
    b, lb = second
    # Distinct trees: each cross-chain pair is counted once (no division by 2).
    count = int(a.count_neighbors(b, CUTOFF))
    distance = float(np.min(b.query(a.data, k=1)[0]))
    return dict(contact_count=count, min_interchain_distance=distance,
                chain1_length=la, chain2_length=lb, chain1_heavy_atoms=a.n, chain2_heavy_atoms=b.n)


def finish_record(row, inherited_severe=False):
    row['contact_definition'] = DEFINITION
    row['cutoff_angstrom'] = CUTOFF
    flags = [x for x in str(row.get('qc_flag', '')).split(';') if x]
    reasons = []
    count = row.get('contact_count')
    calculable = count is not None and np.isfinite(float(count))
    distance = row.get('min_interchain_distance')
    severe = distance is not None and np.isfinite(float(distance)) and float(distance) < 0.5
    if severe or inherited_severe:
        flags.append('severe_overlap_lt_0.5A')
        reasons.append('severe_overlap_lt_0.5A')
    malformed = bool(row.get('malformed_duplicate_identity', False))
    if isinstance(row.get('malformed_chain_ids'), str):
        malformed = bool(set(json.loads(row['malformed_chain_ids'])) & {row.get('chain1'), row.get('chain2')})
    row['malformed_duplicate_identity'] = malformed
    if malformed:
        flags.append('malformed_duplicate_identity')
        reasons.append('malformed_duplicate_identity')
    if row.get('parser_warnings'):
        flags.append('parser_warning_present')
    if row.get('n_coordinate_models', 1) > 1:
        flags.append('multiple_models_first_used')
    if not calculable:
        flags.append('not_calculable')
    row['raw_included'] = calculable
    row['qc_included'] = calculable and not reasons
    row['qc_filter_reason'] = ';'.join(sorted(set(reasons)))
    row['qc_flag'] = ';'.join(sorted(set(flags)))
    row['included_main_distribution'] = calculable and float(count) > 0
    row['included_qc_nonzero_distribution'] = row['qc_included'] and float(count) > 0
    return row
