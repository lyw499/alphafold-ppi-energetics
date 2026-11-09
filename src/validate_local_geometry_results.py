"""Independent output checks, targeted brute force spot checks, and provenance."""
from pathlib import Path
import sys,json,hashlib,platform,datetime
import numpy as np
import pandas as pd
import scipy,Bio
from scipy.spatial.distance import cdist
from atom_contact_common import parse_structure
from parse_pooled_inner_zip_8A import get_residues
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/local_interface_geometry'
def residue(chain,identifier):
 return next(r for r in get_residues(chain) if str(r.id[1])+r.id[2].strip()==str(identifier))
def coords(r):return np.asarray([a.coord for a in r if a.element not in ('H','D')],dtype=float)
def main():
 checks=[];cache={}
 for dataset in ['foldbench','safes']:
  base=OUT/dataset;contacts=pd.read_csv(base/'residue_contacts.csv',dtype={'interface_id':str,'residue_a_id':str,'residue_b_id':str})
  summaries=pd.read_csv(base/'interface_summaries.csv',dtype={'interface_id':str});audit=pd.read_csv(base/'chain_pair_audit.csv',dtype={'interface_id':str})
  assert len(summaries)==2*audit.included_primary.sum()
  for kind,g in contacts.groupby('structure_type'):
   sample=pd.concat([g.nsmallest(1,'d_min'),g.sample(3,random_state=421)])
   for _,row in sample.iterrows():
    path=ROOT/row.source_path
    if path not in cache:cache[path]=parse_structure(path)[0]
    model=cache[path];a=residue(model[row.chain_a],row.residue_a_id);b=residue(model[row.chain_b],row.residue_b_id)
    brute=float(cdist(coords(a),coords(b)).min());atomdist=float(np.linalg.norm(a[row.atom_a_name].coord.astype(float)-b[row.atom_b_name].coord.astype(float)))
    assert abs(brute-row.d_min)<1e-10 and abs(atomdist-row.d_min)<1e-10
    checks.append(dict(dataset=dataset,structure_type=kind,interface_id=row.interface_id,structure_id=row.structure_id,reference_minimum=brute,saved_minimum=row.d_min,nearest_atom_distance=atomdist))
  if dataset=='foldbench':
   old=pd.read_csv(ROOT/'results/foldbench/paired_interfaces_8A.csv').set_index('interface_id')
   for _,r in summaries.iterrows():assert r.n_contacts==old.loc[int(r.interface_id),('experimental' if r.structure_type=='Experimental' else 'af3')+'_contact_count']
  else:
   old=pd.read_csv(ROOT/'results/contact_comparison_residue8A/safes_provisional_two_chain_paired_metrics.csv').set_index('interface_id')
   for _,r in summaries[summaries.structure_id.isin(old.index)].iterrows():assert r.n_contacts==old.loc[r.structure_id,('experimental' if r.structure_type=='Experimental' else 'af3')+'_contact_count']
  paired=pd.read_csv(base/'paired_residue_distances.csv',dtype={'interface_id':str})
  for state in ['both','lost_in_AF3','gained_in_AF3']:
   group=paired[paired.contact_state==state]
   if group.empty:continue
   row=group.sample(1,random_state=213).iloc[0]
   for prefix,kind in [('experimental','Experimental'),('af3','AF3')]:
    src=summaries[(summaries.interface_id==row.interface_id)&(summaries.structure_type==kind)].iloc[0]
    path=ROOT/src.source_path
    if path not in cache:cache[path]=parse_structure(path)[0]
    model=cache[path];a=residue(model[row[prefix+'_chain_a']],row[prefix+'_residue_a']);b=residue(model[row[prefix+'_chain_b']],row[prefix+'_residue_b'])
    brute=float(cdist(coords(a),coords(b)).min());assert abs(brute-row[prefix+'_d_min'])<1e-10
    checks.append(dict(dataset=dataset,structure_type=kind,interface_id=row.interface_id,contact_state=state,reference_minimum=brute,saved_minimum=row[prefix+'_d_min']))
 pd.DataFrame(checks).to_csv(OUT/'independent_geometry_spot_checks.csv',index=False)
 v=json.loads((OUT/'validation.json').read_text());v['independent_output_spot_checks']=len(checks);v['all_FoldBench_contact_counts_equal_legacy']=True;v['SAFES_previous_10_counts_equal_legacy']=True
 previous=json.loads((OUT/'computation_status.json').read_text())['legacy_table_sha256']
 assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in previous.items())
 v['historical_input_table_hashes_unchanged']=True;(OUT/'validation.json').write_text(json.dumps(v,indent=2)+'\n')
 scripts=['local_geometry_common.py','test_local_geometry.py','compute_local_geometry.py','report_local_geometry.py','audit_local_geometry_context.py','validate_local_geometry_results.py','local_geometry.sbatch']
 provenance=dict(recorded_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),HPC_compute_job=65790588,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,biopython=Bio.__version__,source_code_sha256={name:hashlib.sha256((ROOT/'src'/name).read_bytes()).hexdigest() for name in scripts},historical_input_table_sha256=previous,Pedro_processed=False)
 (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 print('Passed:',len(checks),'independent distance/nearest-atom checks, all legacy counts and source hashes.')
if __name__=='__main__':main()
