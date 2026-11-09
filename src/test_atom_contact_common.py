import unittest
import tempfile
from pathlib import Path

from atom_contact_common import atom_contacts, atom_tree, finish_record, parse_structure
from test_foldbench_interfaces import chain


class AtomContactTests(unittest.TestCase):
    def test_atom_pairs_not_residue_pairs(self):
        a = chain('A', [(1, ' ', [('CA', 'C', (0, 0, 0)), ('CB', 'C', (0, 1, 0))])])
        b = chain('B', [(1, ' ', [('CA', 'C', (7, 0, 0)), ('CB', 'C', (7, 1, 0))])])
        self.assertEqual(atom_contacts(atom_tree(a), atom_tree(b))['contact_count'], 4)
        self.assertEqual(atom_contacts(atom_tree(b), atom_tree(a))['contact_count'], 4)

    def test_hydrogen_exclusion_and_inclusive_boundary(self):
        a = chain('A', [(1, ' ', [('CA', 'C', (0, 0, 0)), ('H', 'H', (16, 0, 0))])])
        b = chain('B', [(1, ' ', [('CA', 'C', (8, 0, 0)), ('H', 'H', (0, 0, 0))])])
        self.assertEqual(atom_contacts(atom_tree(a), atom_tree(b))['contact_count'], 1)
        b[(' ', 1, ' ')]['CA'].coord[0] = 8.001
        self.assertEqual(atom_contacts(atom_tree(a), atom_tree(b))['contact_count'], 0)

    def test_raw_retains_zero_and_severe_overlap(self):
        zero = finish_record(dict(contact_count=0, min_interchain_distance=20))
        self.assertTrue(zero['raw_included'] and zero['qc_included'])
        self.assertFalse(zero['included_main_distribution'])
        clash = finish_record(dict(contact_count=9999, min_interchain_distance=.1))
        self.assertTrue(clash['raw_included'] and clash['included_main_distribution'])
        self.assertFalse(clash['qc_included'])

    def test_large_count_and_annotation_alone_do_not_filter(self):
        row = finish_record(dict(contact_count=9999999, min_interchain_distance=2,
                                 qc_flag='low_confidence;mapping_ambiguity;legacy_AF3_has_clash_annotation'))
        self.assertTrue(row['raw_included'] and row['qc_included'])

    def test_unreadable_not_replaced_by_zero(self):
        row = finish_record(dict(count_error='unreadable'))
        self.assertFalse(row['raw_included'] or row['qc_included'])
        self.assertNotIn('contact_count', row)

    def test_duplicate_ion_warning_does_not_remove_protein_interface(self):
        text = ('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n'
                'ATOM      2  CA  ALA B   1       5.000   0.000   0.000  1.00 20.00           C  \n'
                'HETATM    3 CL    CL A   2      30.000   0.000   0.000  1.00 20.00          CL  \n'
                'HETATM    4 CL    CL A   2      40.000   0.000   0.000  1.00 20.00          CL  \n')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.pdb'
            path.write_text(text)
            _, info = parse_structure(path)
        self.assertTrue(info['parser_warnings'])
        row = finish_record(dict(chain1='A', chain2='B', contact_count=1, min_interchain_distance=5, **info))
        self.assertTrue(row['raw_included'] and row['qc_included'])

    def test_duplicate_eligible_atom_is_pair_scoped_qc(self):
        text = ('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n'
                'ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C  \n'
                'ATOM      3  CA  ALA B   1       5.000   0.000   0.000  1.00 20.00           C  \n')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.pdb'
            path.write_text(text)
            _, info = parse_structure(path)
        bad = finish_record(dict(chain1='A', chain2='B', contact_count=1, min_interchain_distance=5, **info))
        good = finish_record(dict(chain1='B', chain2='C', contact_count=1, min_interchain_distance=5, **info))
        self.assertTrue(bad['raw_included'])
        self.assertFalse(bad['qc_included'])
        self.assertTrue(good['qc_included'])


if __name__ == '__main__':
    unittest.main()
