"""Small synthetic checks; no benchmark archive or network access required."""
import unittest

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue

from foldbench_interfaces import map_chains, metrics
from prepare_foldbench import select_summary


def chain(name, entries):
    result = Chain(name)
    for number, insertion, atoms in entries:
        residue = Residue((' ', number, insertion), 'ALA', '')
        for serial, (atom_name, element, xyz) in enumerate(atoms, 1):
            residue.add(Atom(atom_name, np.array(xyz, dtype=float), 0, 1, ' ', atom_name, serial, element=element))
        result.add(residue)
    return result


class InterfaceTests(unittest.TestCase):
    def test_residue_pairs_not_atom_pairs_and_hydrogen_exclusion(self):
        a = chain('A', [(1, ' ', [('CA', 'C', (0, 0, 0)), ('CB', 'C', (0, 1, 0))]),
                        (2, ' ', [('CA', 'C', (100, 0, 0)), ('H', 'H', (8, 0, 0))])])
        b = chain('B', [(1, ' ', [('CA', 'C', (7, 0, 0)), ('CB', 'C', (7, 1, 0))]),
                        (2, ' ', [('CA', 'C', (200, 0, 0)), ('H', 'H', (0, 0, 0))])])
        m = metrics(a, b)
        self.assertEqual(m['contact_count'], 1)
        self.assertEqual(m['total_interface_residues'], 2)
        self.assertEqual(m['fraction_interface_chain1'], .5)
        self.assertEqual(m['fraction_interface_chain2'], .5)

    def test_boundary_insertion_codes_zero_and_symmetry(self):
        a = chain('A', [(1, ' ', [('CA', 'C', (0, 0, 0))]),
                        (1, 'A', [('CA', 'C', (0, 0, 0))])])
        b = chain('B', [(1, ' ', [('CA', 'C', (8, 0, 0))])])
        self.assertEqual(metrics(a, b)['contact_count'], 2)
        self.assertEqual(metrics(b, a)['contact_count'], 2)
        b[(' ', 1, ' ')]['CA'].coord[0] = 8.001
        m = metrics(a, b)
        self.assertEqual(m['contact_count'], 0)
        self.assertEqual(m['fraction_interface_chain1'], 0)

    def test_seed_selection_and_numeric_tie_break(self):
        text = ('target,seed-sample,ranking_score,ptm,iptm\n'
                'example,seed_10_sample_1,0.9,0.8,0.9\n'
                'example,seed_2_sample_3,0.9,0.8,0.9\n'
                'example,seed_1_sample_0,0.8,0.7,0.8\n')
        r = select_summary(text, 'example')
        self.assertEqual((r['selected_seed'], r['selected_sample']), (2, 3))
        self.assertEqual(r['n_top_score_ties'], 2)

    def test_sequence_mapping_checks_names_instead_of_assuming(self):
        exp, pred = Model(0), Model(0)
        a = chain('A', [(i, ' ', [('CA', 'C', (i, i % 2, 0))]) for i in range(1, 5)])
        b = chain('B', [(i, ' ', [('CA', 'C', (i, i % 2, 5))]) for i in range(1, 5)])
        for r in b:
            r.resname = 'GLY'
        exp.add(a)
        exp.add(b)
        pa, pb = b.copy(), a.copy()
        pa.id, pb.id = 'A', 'B'
        pred.add(pa)
        pred.add(pb)
        mapping, _, _, rmsd = map_chains(exp, pred, {'A', 'B'})
        self.assertEqual(mapping, {'A': 'B', 'B': 'A'})
        self.assertLess(rmsd, 1e-8)

    def test_global_geometry_resolves_identical_sequences_after_rotation(self):
        exp, pred = Model(0), Model(0)
        rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        expected = {'A': 'Z', 'B': 'X', 'C': 'Y'}
        local = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 2]])
        for label, center in zip('ABC', [[0, 0, 0], [15, 2, 0], [2, 21, 9]]):
            xyz = local + np.array(center)
            exp.add(chain(label, [(i, ' ', [('CA', 'C', p)]) for i, p in enumerate(xyz, 1)]))
            moved = xyz @ rotation + np.array([40, -20, 10])
            pred.add(chain(expected[label], [(i, ' ', [('CA', 'C', p)]) for i, p in enumerate(moved, 1)]))
        mapping, _, ties, rmsd = map_chains(exp, pred, set('ABC'))
        self.assertEqual(mapping, expected)
        self.assertEqual(len(ties), 1)
        self.assertLess(rmsd, 1e-8)


if __name__ == '__main__':
    unittest.main()
