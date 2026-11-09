"""Paired interface local geometry. All writes confined to local_interface_geometry."""
import json,itertools,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from atom_contact_common import parse_structure
from foldbench_interfaces import map_chains,protein_residues
from local_geometry_common import pack_chain,min_contacts,exact_targeted_distance,summarize,residue_map
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/local_interface_geometry'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run_dataset(dataset):
 audits=[];summaries=[];contacts=[];paired=[];residue_maps=[];tie_rows=[];sources=[]
 if dataset=='foldbench':
  table=pd.read_csv(ROOT/'results/foldbench/paired_interfaces_8A.csv');groups=list(table.groupby('target_id',sort=False))
 else:
  groups=[(p.stem,p) for p in sorted((ROOT/'data/SAFES_PDBs').glob('*.pdb')) if not p.stem.endswith('_af3')]
 for num,(sid,g) in enumerate(groups,1):
  if dataset=='foldbench':
   if not g.paired_success.any():
    for _,r in g.iterrows():audits.append(dict(dataset=dataset,interface_id=str(int(r.interface_id)),structure_id=sid,experimental_chain_a=r.chain1,experimental_chain_b=r.chain2,status='unavailable_AF3',failure_reason='target_absent_from_provided_AF3_archive',included_primary=False,legacy_mapping_stable=False))
    continue
   ep=ROOT/g.iloc[0].experimental_path;pp=ROOT/g.iloc[0].af3_path
  else:ep=g;pp=g.with_name(g.stem+'_af3.pdb')
  try:em,ei=parse_structure(ep);pm,pi=parse_structure(pp)
  except Exception as exc:
   audits.append(dict(dataset=dataset,structure_id=sid,status='parsing_error',failure_reason=str(exc),included_primary=False));continue
  sources.extend([dict(dataset=dataset,structure_id=sid,structure_type=kind,path=str(path.relative_to(ROOT)),sha256=digest(path),**info) for kind,path,info in [('Experimental',ep,ei),('AF3',pp,pi)]])
  trees_e={c.id:pack_chain(c) for c in em if protein_residues(c)[1]};trees_p={c.id:pack_chain(c) for c in pm if protein_residues(c)[1]}
  specifications=[]
  if dataset=='foldbench':
   for _,r in g.iterrows():
    pairs=json.loads(r.tied_interface_chain_pairs_json) if isinstance(r.tied_interface_chain_pairs_json,str) else [[r.af3_chain1,r.af3_chain2]]
    sensitive=any(pd.notna(r[f'mapping_tie_{m}_min']) and r[f'mapping_tie_{m}_min']!=r[f'mapping_tie_{m}_max'] for m in ['contact_count','total_interface_residues','fraction_interface_chain1','fraction_interface_chain2'])
    specifications.append(dict(id=str(int(r.interface_id)),ea=r.chain1,eb=r.chain2,pa=r.af3_chain1,pb=r.af3_chain2,alternatives=pairs,reliable=True,ordered_unique=len(pairs)==1,legacy_mapping_stable=not sensitive,legacy_flags=str(r.qc_flags) if pd.notna(r.qc_flags) else '',expected_exp=int(r.experimental_contact_count),expected_af3=int(r.af3_contact_count)))
  else:
   mapping={};ties=[];error=''
   try:
    # Allow reliably matchable chains even if some chains have no sequence candidate.
    mapping,details,ties,rmsd=map_chains(em,pm,[])
    residue_maps.extend(dict(dataset=dataset,structure_id=sid,record_type='chain_mapping',**d) for d in details)
   except Exception as exc:error=str(exc)
   for a,b in itertools.combinations(trees_e,2):
    alternatives=sorted({(t[a],t[b]) for t in ties if a in t and b in t})
    chain_pair_unique=len({frozenset(v) for v in alternatives})==1
    reliable=bool(alternatives) and chain_pair_unique
    specifications.append(dict(id=f'{sid}:{a}:{b}',ea=a,eb=b,pa=mapping.get(a,''),pb=mapping.get(b,''),alternatives=alternatives,reliable=reliable,ordered_unique=len(alternatives)==1,legacy_mapping_stable=reliable,legacy_flags='all_corresponding_chain_pairs_not_curated_binding_partners',mapping_error=error or ('unmapped_required_chain' if not alternatives else 'ambiguous_chain_pair' if not reliable else '')))
  exp_cache={};pred_cache={};rm_cache={}
  for spec in specifications:
   iid=spec['id'];ea,eb,pa,pb=[spec[k] for k in ['ea','eb','pa','pb']]
   status='mapped_unique_ordered_pair' if spec['ordered_unique'] else 'mapped_chain_copy_ambiguity' if spec['reliable'] else 'unreliable_mapping'
   audit=dict(dataset=dataset,interface_id=iid,pair_identifier=f'{dataset}:{iid}',structure_id=sid,experimental_chain_a=ea,experimental_chain_b=eb,af3_chain_a=pa,af3_chain_b=pb,status=status,included_primary=spec['reliable'],chain_mapping_unambiguous=spec['ordered_unique'],legacy_mapping_stable=spec['legacy_mapping_stable'],failure_reason=spec.get('mapping_error',''),candidate_af3_pairs=json.dumps(spec['alternatives']),existing_qc_flags=spec['legacy_flags'])
   erows=exp_cache.setdefault((ea,eb),min_contacts(trees_e[ea],trees_e[eb]))
   prows=[]
   if spec['reliable']:
    if (pa,pb) not in pred_cache:pred_cache[pa,pb]=min_contacts(trees_p[pa],trees_p[pb])
    prows=pred_cache[pa,pb]
   audit['experimental_n_contacts']=len(erows);audit['af3_n_contacts']=len(prows) if spec['reliable'] else np.nan
   if dataset=='foldbench':
    assert len(erows)==spec['expected_exp'],(sid,iid,'experimental residue count differs from legacy')
    assert len(prows)==spec['expected_af3'],(sid,iid,'AF3 residue count differs from legacy')
   local_stable=True
   for ta,tb in spec['alternatives']:
    if not spec['reliable']:continue
    if (ta,tb) not in pred_cache:pred_cache[ta,tb]=min_contacts(trees_p[ta],trees_p[tb])
    alt=pred_cache[ta,tb];d1=sorted(r['d_min'] for r in prows);d2=sorted(r['d_min'] for r in alt)
    same=len(d1)==len(d2) and np.allclose(d1,d2,atol=1e-6,rtol=0)
    local_stable &= bool(same)
    if len(spec['alternatives'])>1:tie_rows.append(dict(dataset=dataset,interface_id=iid,structure_id=sid,chain_a=ta,chain_b=tb,distance_spectrum_equal_to_selected=same,**summarize(alt)))
   audit['local_distance_spectrum_stable']=local_stable if spec['reliable'] else False
   # Do not pool an unmatched experimental chain pair into a paired benchmark.
   if not spec['reliable']:
    audits.append(audit);continue
   for kind,path,info,a,b,rows,trees in [('Experimental',ep,ei,ea,eb,erows,trees_e),('AF3',pp,pi,pa,pb,prows,trees_p)]:
    malformed=bool({a,b}&set(json.loads(info['malformed_chain_ids'])))
    flags=spec['legacy_flags']+(';' if spec['legacy_flags'] else '')+('mapping_ambiguity;' if not spec['ordered_unique'] else '')+('parser_warning;' if info['parser_warnings'] else '')+('malformed_duplicate_identity;' if malformed else '')
    meta=dict(dataset=dataset,interface_id=iid,pair_identifier=f'{dataset}:{iid}',structure_type=kind,structure_id=sid,chain_a=a,chain_b=b,source_path=str(path.relative_to(ROOT)),existing_qc_flags=flags.rstrip(';'),mapping_status=status,legacy_mapping_stable=spec['legacy_mapping_stable'],chain_mapping_unambiguous=spec['ordered_unique'],local_distance_spectrum_stable=local_stable)
    summary=summarize(rows)
    summaries.append(dict(**meta,chain_a_length=len(trees[a]['res']),chain_b_length=len(trees[b]['res']),n_deuterium_excluded=trees[a]['n_deuterium']+trees[b]['n_deuterium'],parser_warnings=info['parser_warnings'],malformed_duplicate_identity=malformed,**summary))
    contacts.extend(dict(**meta,**r) for r in rows)
   # Residue mapping is deliberately stricter than the interface-level selection.
   if not spec['ordered_unique']:
    audit['residue_pair_mapping_status']='skipped_ambiguous_chain_correspondence';audits.append(audit);continue
   maps=[]
   for ex,pr in [(ea,pa),(eb,pb)]:
    if (ex,pr) not in rm_cache:
     m,details=residue_map(em[ex],pm[pr]);rm_cache[ex,pr]=(m,details)
     residue_maps.append(dict(dataset=dataset,structure_id=sid,record_type='residue_mapping_summary',experimental_chain=ex,af3_chain=pr,**details))
    maps.append(rm_cache[ex,pr][0])
   ma,mb=maps;ia={v:k for k,v in ma.items()};ib={v:k for k,v in mb.items()}
   ed={(r['residue_index_a'],r['residue_index_b']):r['d_min'] for r in erows};pdct={(r['residue_index_a'],r['residue_index_b']):r['d_min'] for r in prows}
   union={(i,j) for i,j in ed if i in ma and j in mb}|{(ia[i],ib[j]) for i,j in pdct if i in ia and j in ib}
   audit.update(residue_pair_mapping_status='stable_identical_sequence_positions_only',n_experimental_contacts_residue_mappable=sum(i in ma and j in mb for i,j in ed),n_af3_contacts_residue_mappable=sum(i in ia and j in ib for i,j in pdct),n_union_residue_pairs=len(union))
   for i,j in sorted(union):
    u,v=ma[i],mb[j]
    ex=ed.get((i,j));pr=pdct.get((u,v))
    if ex is None:ex=exact_targeted_distance(trees_e[ea],trees_e[eb],i,j)
    if pr is None:pr=exact_targeted_distance(trees_p[pa],trees_p[pb],u,v)
    er1,er2=trees_e[ea]['res'][i],trees_e[eb]['res'][j];ar1,ar2=trees_p[pa]['res'][u],trees_p[pb]['res'][v]
    from local_geometry_common import residue_id
    paired.append(dict(dataset=dataset,interface_id=iid,structure_id=sid,experimental_chain_a=ea,experimental_chain_b=eb,af3_chain_a=pa,af3_chain_b=pb,experimental_residue_a=residue_id(er1),experimental_residue_b=residue_id(er2),af3_residue_a=residue_id(ar1),af3_residue_b=residue_id(ar2),residue_a_type=er1.resname,residue_b_type=er2.resname,experimental_d_min=ex,af3_d_min=pr,delta_d_min=pr-ex,experimental_contact=ex<=8,af3_contact=pr<=8,contact_state='both' if ex<=8 and pr<=8 else 'lost_in_AF3' if ex<=8 else 'gained_in_AF3',mapping='same_amino_acids_at_positions_common_to_all_optimal_alignments',existing_qc_flags=spec['legacy_flags']))
   audits.append(audit)
  if num%10==0:print(dataset,num,'/',len(groups),'structures processed',flush=True)
 dest=OUT/dataset;dest.mkdir(exist_ok=True,parents=True)
 for name,rows in [('chain_pair_audit',audits),('interface_summaries',summaries),('residue_contacts',contacts),('paired_residue_distances',paired),('residue_mapping_audit',residue_maps),('mapping_tie_geometry_sensitivity',tie_rows),('structure_source_audit',sources)]:
  frame=pd.DataFrame(rows)
  if frame.empty:frame=pd.DataFrame(columns=['dataset','interface_id'])
  frame.to_csv(dest/f'{name}.csv',index=False)
 print(dataset,'DONE',len(summaries),'structure-interface observations',len(contacts),'residue contacts',len(paired),'mapped union residue pairs',flush=True)

def main():
 assert json.loads((OUT/'algorithm_validation.json').read_text())['status']=='passed','Validate optimized algorithm before full run'
 inputs=[ROOT/'results/foldbench/paired_interfaces_8A.csv',ROOT/'results/contact_comparison_residue8A/foldbench_matching_audit.csv',ROOT/'results/contact_comparison_residue8A/safes_provisional_two_chain_paired_metrics.csv']
 before={str(p.relative_to(ROOT)):digest(p) for p in inputs}
 run_dataset('foldbench');run_dataset('safes')
 assert all(digest(ROOT/p)==h for p,h in before.items())
 (OUT/'computation_status.json').write_text(json.dumps(dict(status='complete',datasets=['foldbench','safes'],Pedro_processed=False,legacy_table_sha256=before),indent=2)+'\n')
if __name__=='__main__':main()
