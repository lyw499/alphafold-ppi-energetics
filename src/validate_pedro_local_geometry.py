"""Redirect benchmark brute-force test outputs into the new Pedro validation directory."""
import json
from pathlib import Path
import parse_pedro_local_geometry as parser
import test_local_geometry as brute
import pandas as pd
OUT=parser.OUT
class ValidationRoot:
 def __truediv__(self,relative):
  return OUT/'validation/brute_force_validation.json' if str(relative)=='results/local_interface_geometry/algorithm_validation.json' else parser.ROOT/relative
brute.ROOT=ValidationRoot()
brute.main()
chosen=pd.read_csv(OUT/'validation/selected_jobs.csv')
for idx,g in chosen.groupby('archive_index'):parser.run(int(idx),set(g.job_key))
table=pd.concat([pd.read_csv(OUT/f'audit_parts/{k}.csv') for k in chosen.job_key],ignore_index=True)
table.rename(columns={'selected_model':'model'}).to_csv(OUT/'pedro_parser_validation.csv',index=False)
assert len(table)==int(chosen.n_interfaces.sum()) and table.status.eq('passed').all() and table.difference.eq(0).all()
parser.atomic_json(OUT/'validation/status.json',dict(status='passed',interfaces=len(table),jobs=len(chosen),zero_interfaces=int(table.new_contact_count.eq(0).sum()),contacts=int(table.new_contact_count.sum()),max_contacts=int(table.new_contact_count.max()),deuterium_atoms=int(table.n_deuterium.sum()),protocol=parser.fingerprint()))
print((OUT/'validation/status.json').read_text(),flush=True)
