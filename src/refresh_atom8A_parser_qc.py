"""Refresh warning-related QC only; never recalculate contacts or select models."""
from pathlib import Path

import pandas as pd

from atom_contact_common import finish_record, parse_structure

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/contact_comparison_atom8A'


def main():
    names = ['safes_experimental', 'safes_af3', 'mgen_reference_experimental', 'foldbench_experimental', 'foldbench_af3']
    for name in names:
        path = OUT / (name + '.csv')
        data = pd.read_csv(path).astype(object)
        warned = data[data.parser_warnings.fillna('') != '']
        for source, group in warned.groupby('source_path'):
            _, info = parse_structure(ROOT / source)
            for index in group.index:
                row = data.loc[index].to_dict()
                flags = [f for f in str(row.get('qc_flag', '')).split(';') if f not in ('malformed_duplicate_identity', 'parser_warning_present', 'nan')]
                row['qc_flag'] = ';'.join(flags)
                row.update(info)
                row = finish_record(row)
                for key, value in row.items():
                    if key not in data:
                        data[key] = None
                    data.at[index, key] = value
        data.to_csv(path, index=False)
        print(name, 'raw', data.raw_included.sum(), 'qc', data.qc_included.sum())


if __name__ == '__main__':
    main()
