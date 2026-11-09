"""Explain SAFES mapping failures and short-distance chemical context, without exclusions."""
from pathlib import Path
import itertools
import numpy as np
import json
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from foldbench_interfaces import load_model,protein_residues
from atom_contact_common import parse_structure
from parse_pooled_inner_zip_8A import get_residues
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/local_interface_geometry'
def main():
 rows=[]
 for sid in ['1OGA','5E9D']:
  e,_=load_model(ROOT/'data/SAFES_PDBs'/f'{sid}.pdb');p,_=load_model(ROOT/'data/SAFES_PDBs'/f'{sid}_af3.pdb')
  for ec in e:
   er,es=protein_residues(ec)
   if not es:continue
   for pc in p:
    pr,ps=protein_residues(pc)
    if not ps:continue
    al=PairwiseAligner();al.mode='global';al.match_score=2;al.mismatch_score=-3;al.open_gap_score=-5;al.extend_gap_score=-.5;al.end_gap_score=0
    match=aligned=nca=0;aln=al.align(es,ps)[0]
    for ea,pa in zip(*aln.aligned):
     for i,j in zip(range(*ea),range(*pa)):
      aligned+=1
      if es[i]==ps[j] and es[i]!='X':match+=1;nca+=int('CA' in er[i] and 'CA' in pr[j])
    identity=match/aligned if aligned else 0;coverage=aligned/len(es)
    reasons=[]
    if len(es)>len(ps)/.9:reasons.append('AF3_coordinate_chain_too_short_for_90pct_experimental_coverage')
    if identity<.95:reasons.append('identity_below_0.95')
    if coverage<.9:reasons.append('experimental_observed_coverage_below_0.90')
    if nca<3:reasons.append('fewer_than_3_matched_CA')
    rows.append(dict(structure_id=sid,experimental_chain=ec.id,af3_chain=pc.id,experimental_coordinate_residues=len(es),af3_coordinate_residues=len(ps),identity=identity,experimental_observed_coverage=coverage,n_matched_CA=nca,passes_existing_candidate_rules=not reasons,failure_reasons=';'.join(reasons)))
 data=pd.DataFrame(rows);data.to_csv(OUT/'safes/failed_target_all_chain_sequence_candidates.csv',index=False)
 best=data.sort_values(['structure_id','experimental_chain','n_matched_CA','experimental_observed_coverage','identity'],ascending=[True,True,False,False,False]).groupby(['structure_id','experimental_chain'],sort=False).head(1)
 best.to_csv(OUT/'safes/failed_target_best_chain_candidates.csv',index=False)
 short=pd.read_csv(OUT/'qc/extreme_short_distance_cases.csv',dtype={'residue_a_id':str,'residue_b_id':str})
 short=short[short.d_min<2].copy()
 # Distinguish merely SG-SG-like geometry from bonds annotated in the experimental CIF.
 cache={};model_cache={};alternate_rows=[]
 for idx,row in short.iterrows():
  status='not_available_in_supplied_coordinate_only_PDB';evidence=[]
  path=ROOT/row.source_path
  if path not in model_cache:model_cache[path]=parse_structure(path)[0]
  model=model_cache[path]
  for side in ['a','b']:
   atom=get_residues(model[row['chain_'+side]])[int(row['residue_index_'+side])][row['atom_'+side+'_name']]
   short.at[idx,'selected_atom_'+side+'_occupancy']=atom.get_occupancy()
   short.at[idx,'selected_atom_'+side+'_altloc']=atom.get_altloc().strip()
  atom_a=get_residues(model[row.chain_a])[int(row.residue_index_a)][row.atom_a_name]
  atom_b=get_residues(model[row.chain_b])[int(row.residue_index_b)][row.atom_b_name]
  alternatives_a=list(atom_a.child_dict.values()) if atom_a.is_disordered()==2 else [atom_a]
  alternatives_b=list(atom_b.child_dict.values()) if atom_b.is_disordered()==2 else [atom_b]
  if len(alternatives_a)>1 or len(alternatives_b)>1:
   for aa,bb in itertools.product(alternatives_a,alternatives_b):
    alternate_rows.append(dict(dataset=row.dataset,structure_type=row.structure_type,structure_id=row.structure_id,interface_id=row.interface_id,chain_a=row.chain_a,chain_b=row.chain_b,residue_a_id=row.residue_a_id,residue_b_id=row.residue_b_id,atom_a_name=aa.name,atom_b_name=bb.name,altloc_a=aa.altloc,altloc_b=bb.altloc,occupancy_a=aa.occupancy,occupancy_b=bb.occupancy,atom_distance_A=float(np.linalg.norm(aa.coord.astype(float)-bb.coord.astype(float))),note='Alternative positions of the originally closest atoms only; not a replacement residue minimum or coordinated conformer reconstruction'))
  if path.suffix=='.cif':
   if path not in cache:cache[path]=MMCIF2Dict(str(path))
   d=cache[path];keys=[k for k in d if k.startswith('_struct_conn.')];status='no_matching_declared_connection'
   for values in zip(*(d[k] for k in keys)):
    conn=dict(zip([k.split('.',1)[1] for k in keys],values))
    atom1=(conn.get('ptnr1_label_asym_id'),conn.get('ptnr1_auth_seq_id'),conn.get('ptnr1_label_atom_id'))
    atom2=(conn.get('ptnr2_label_asym_id'),conn.get('ptnr2_auth_seq_id'),conn.get('ptnr2_label_atom_id'))
    wanted={(row.chain_a,row.residue_a_id,row.atom_a_name),(row.chain_b,row.residue_b_id,row.atom_b_name)}
    if {atom1,atom2}==wanted:evidence.append(conn.get('conn_type_id','unknown'))
   if evidence:status='declared_connection:'+','.join(evidence)
  short.at[idx,'source_bond_annotation_status']=status
 short.to_csv(OUT/'qc/short_distance_chemical_context.csv',index=False)
 pd.DataFrame(alternate_rows).to_csv(OUT/'qc/short_distance_altloc_diagnostics.csv',index=False)
 counts=[]
 for (dataset,kind),g in short.groupby(['dataset','structure_type']):
  sg=(g.atom_a_name=='SG')&(g.atom_b_name=='SG')
  counts.append(dict(dataset=dataset,structure_type=kind,n_lt2=len(g),n_SG_SG_like=int(sg.sum()),n_non_SG_SG_short=int((~sg).sum()),n_lt1_5=int(g.d_min.lt(1.5).sum())))
 pd.DataFrame(counts).to_csv(OUT/'qc/short_distance_chemical_context_counts.csv',index=False)
 print(best.to_string(index=False));print(pd.DataFrame(counts).to_string(index=False))
if __name__=='__main__':main()
