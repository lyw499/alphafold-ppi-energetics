"""Reuse validated residue counts; preserve zeros, source scope, and conservative QC."""
import os,json
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/residue8A-mpl')
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr
from scipy.stats import ks_2samp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from foldbench_interfaces import metrics,load_model,map_chains,protein_residues
from atom_contact_common import parse_structure
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/contact_comparison_residue8A'
METRICS=['contact_count','total_interface_residues','interface_fraction']
LABEL={'contact_count':'Residue-pair contacts (8 Å)','total_interface_residues':'Total interface residues','interface_fraction':'Interface fraction (pooled coordinate residues)'}

def standard(row):
 row['interface_fraction']=row['total_interface_residues']/(row['chain1_length']+row['chain2_length'])
 row['qc_included']=not row.get('severe_overlap',False)
 row['nonzero']=row['contact_count']>0
 row['contact_definition']='unique_interchain_residue_pairs_any_nonH_atom_distance_le_8A'
 return row

def desc(x):
 return dict(n=len(x),mean=x.mean(),median=x.median(),std=x.std(ddof=1),q1=x.quantile(.25),q3=x.quantile(.75),min=x.min(),max=x.max())

def savefig(fig,folder,name):
 p=OUT/'figures'/folder;p.mkdir(parents=True,exist_ok=True)
 fig.tight_layout();fig.savefig(p/(name+'.png'),dpi=160);fig.savefig(p/(name+'.pdf'));plt.close(fig)

def distribution(a,b,labels,folder,stem):
 rows=[]
 for metric in METRICS:
  x=a.loc[a.contact_count>0,metric].dropna();y=b.loc[b.contact_count>0,metric].dropna()
  if not len(x) or not len(y):continue
  # Shared bins in transformed coordinates, normalized separately; no clipping.
  transform=(lambda v:np.log10(v)) if metric!='interface_fraction' else (lambda v:v)
  xx,yy=transform(x),transform(y);lo=min(xx.min(),yy.min());hi=max(xx.max(),yy.max())
  bins=np.linspace(lo if lo!=hi else lo-.5,hi if lo!=hi else hi+.5,31)
  fig,ax=plt.subplots(figsize=(7,4))
  for values,label in [(xx,labels[0]),(yy,labels[1])]:
   ax.hist(values,bins=bins,density=True,histtype='step',linewidth=1.7,label=f'{label} (n={len(values):,})')
  ax.set(xlabel=('log10 ' if metric!='interface_fraction' else '')+LABEL[metric],ylabel='Normalized density',title='Non-zero interface distributions');ax.legend(fontsize=8)
  savefig(fig,folder,f'{stem}_{metric}_distribution')
  fig,ax=plt.subplots(figsize=(7,4))
  for values,label in [(x,labels[0]),(y,labels[1])]:
   v=np.sort(values);ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',label=f'{label} (n={len(v):,})')
  if metric!='interface_fraction':ax.set_xscale('log')
  ax.set(xlabel=LABEL[metric],ylabel='ECDF',title='Non-zero interface distributions');ax.legend(fontsize=8)
  savefig(fig,folder,f'{stem}_{metric}_ecdf')
  for values,label in [(x,labels[0]),(y,labels[1])]:rows.append(dict(comparison=stem,group=label,metric=metric,**desc(values)))
 return rows

def paired_analysis(exp,pred,name):
 e=exp.set_index('interface_id');p=pred.set_index('interface_id');ids=e.index.intersection(p.index);e=e.loc[ids];p=p.loc[ids]
 assert e.index.is_unique and p.index.is_unique
 wide=e.add_prefix('experimental_').join(p.add_prefix('af3_'));wide.to_csv(OUT/f'{name}_paired_metrics.csv')
 summaries=distribution(e,p,['Experimental','AF3'],name,f'{name}_exp_vs_af3') or []
 stats=[]
 for metric in METRICS:
  x=e[metric];y=p[metric];delta=y-x
  valid=x.notna()&y.notna();x=x[valid];y=y[valid];delta=delta[valid]
  pear=pearsonr(x,y).statistic if len(x)>1 and x.nunique()>1 and y.nunique()>1 else np.nan
  spear=spearmanr(x,y).statistic if len(x)>1 and x.nunique()>1 and y.nunique()>1 else np.nan
  stats.append(dict(comparison=name,metric=metric,n=len(x),pearson=pear,spearman=spear,median_difference=delta.median(),mean_difference=delta.mean(),q1_difference=delta.quantile(.25),q3_difference=delta.quantile(.75),std_difference=delta.std(),mae=delta.abs().mean(),rmse=np.sqrt((delta**2).mean()),n_paired_with_any_zero=int(((e.contact_count==0)|(p.contact_count==0)).sum())))
  fig,ax=plt.subplots(figsize=(5,5));ax.scatter(x,y,s=17,alpha=.6);upper=max(x.max(),y.max())*1.03;ax.plot([0,upper],[0,upper],color='black',linestyle='--',linewidth=1);ax.set(xlabel='Experimental: '+LABEL[metric],ylabel='AF3: '+LABEL[metric],xlim=(0,upper),ylim=(0,upper),title=f'{name}: valid paired cases, zeros retained\nn={len(x)}, Pearson={pear:.3f}, Spearman={spear:.3f}')
  savefig(fig,name,f'{name}_exp_vs_af3_{metric}_scatter')
  wide[f'difference_{metric}']=delta
 wide.to_csv(OUT/f'{name}_paired_metrics.csv');pd.DataFrame(stats).to_csv(OUT/f'{name}_paired_statistics.csv',index=False)
 return summaries

def foldbench():
 audit=pd.read_csv(OUT/'foldbench_matching_audit.csv');assert len(audit)==279 and audit.af3_file_available.sum()==132
 old=pd.read_csv(ROOT/'results/foldbench/paired_interfaces_8A.csv');rows=[]
 for _,r in old.iterrows():
  for kind,prefix in [('experimental','experimental'),('AF3','af3')]:
   if pd.isna(r[prefix+'_contact_count']):continue
   row=dict(dataset='FoldBench',structure_type=kind,interface_id=r.interface_id,structure_id=r.target_id,chain1=r.chain1 if kind=='experimental' else r.af3_chain1,chain2=r.chain2 if kind=='experimental' else r.af3_chain2,paired_success=r.paired_success,source='results/foldbench/paired_interfaces_8A.csv',qc_flag=r.qc_flags)
   for col in ['contact_count','total_interface_residues','chain1_length','chain2_length','interface_residues_chain1','interface_residues_chain2','fraction_interface_chain1','fraction_interface_chain2']:row[col]=r[prefix+'_'+col]
   row['severe_overlap']=r[prefix+'_min_interchain_distance_A']<.5
   row['mapping_metric_sensitive']=any(pd.notna(r[f'mapping_tie_{m}_min']) and r[f'mapping_tie_{m}_min']!=r[f'mapping_tie_{m}_max'] for m in ['contact_count','total_interface_residues','fraction_interface_chain1','fraction_interface_chain2'])
   rows.append(standard(row))
 data=pd.DataFrame(rows);data.to_csv(OUT/'foldbench_all_computable.csv',index=False)
 return data

def pedro():
 frames=[]
 cols=['job_id','protein_pair','chain_a','chain_b','best_model','length_a','length_b','contact_count','interface_residues_a','interface_residues_b','min_interchain_distance','has_clash']
 for path in sorted((ROOT/'results/pooled_parts_8A').glob('*.csv')):
  d=pd.read_csv(path,usecols=cols);d['source_file']=path.name;frames.append(d)
 d=pd.concat(frames,ignore_index=True);assert len(d)==170644
 d=d.rename(columns={'length_a':'chain1_length','length_b':'chain2_length','chain_a':'chain1','chain_b':'chain2','interface_residues_a':'interface_residues_chain1','interface_residues_b':'interface_residues_chain2'})
 d['dataset']='Pedro';d['structure_type']='AF3';d['structure_id']=d.source_file+':'+d.job_id+':'+d.best_model.astype(str)
 d['total_interface_residues']=d.interface_residues_chain1+d.interface_residues_chain2;d['interface_fraction']=d.total_interface_residues/(d.chain1_length+d.chain2_length)
 d['severe_overlap']=d.min_interchain_distance<.5;d['qc_included']=~d.severe_overlap;d['nonzero']=d.contact_count>0;d['qc_flag']=np.where(d.severe_overlap,'severe_overlap_min_distance_lt_0.5A','');d['contact_definition']='unique_interchain_residue_pairs_any_nonH_atom_distance_le_8A'
 assert not d[METRICS].isna().any().any()
 for version,x in [('raw',d),('qc',d[d.qc_included])]:
  x.to_csv(OUT/f'pedro_{version}.csv',index=False)
  jobs=x.groupby(['job_id','protein_pair'])[METRICS].mean();pairs=jobs.groupby('protein_pair')[METRICS].mean();pairs.to_csv(OUT/f'pedro_pair_level_{version}.csv')
 return d

def reference():
 source=pd.read_csv(ROOT/'results/contact_comparison_atom8A/mgen_reference_experimental.csv');rows=[];audit=[]
 for sid,g in source.groupby('structure_id',sort=False):
  model,info=parse_structure(ROOT/g.iloc[0].source_path);header=MMCIF2Dict(str(ROOT/'data/mgen_reference_atom8A/metadata'/f'{g.iloc[0].pdb_id}.cif'))
  descs=dict(zip(header['_entity.id'],header['_entity.pdbx_description']));assembly=dict(zip(header['_pdbx_struct_assembly.id'],header['_pdbx_struct_assembly.details']))
  for _,r in g.iterrows():
   m=metrics(model[r.chain1],model[r.chain2]);m={k:v for k,v in m.items() if 'residue_ids' not in k}
   row=standard(dict(dataset='Reference_candidates',structure_type='experimental',interface_id=sid+':'+r.chain_pair,structure_id=sid,chain1=r.chain1,chain2=r.chain2,source=r.source_path,reference_taxonomy_scope=r.reference_taxonomy_scope,severe_overlap=m['min_interchain_distance_A']<.5,**m))
   rows.append(row)
   audit.append(dict(interface_id=row['interface_id'],structure_id=sid,chain1=r.chain1,chain2=r.chain2,chain1_source=r.chain1_source_organisms,chain2_source=r.chain2_source_organisms,reference_taxonomy_scope=r.reference_taxonomy_scope,assembly_evidence=assembly.get(str(r.assembly_id),''),inclusion_reason='Previously enumerated protein-chain pair within deposited biological assembly; retained as candidate audit scope, not proof of pair-specific biological interface.',intended_biological_interface='not_independently_established_for_each_pair',main_comparison_included=r.reference_taxonomy_scope!='neither_Mgen',scope_exclusion_reason='neither chain derives from Mgen' if r.reference_taxonomy_scope=='neither_Mgen' else '',contact_count=m['contact_count'],zero_contact=m['contact_count']==0))
 d=pd.DataFrame(rows);d.to_csv(OUT/'reference_all_29_candidates.csv',index=False);pd.DataFrame(audit).to_csv(OUT/'reference_interface_selection_audit.csv',index=False)
 return d

def safes():
 # No inferred multi-chain binding-partner assignment: preserve an explicit unresolved inventory.
 rows=[];audit=[];mapping_rows=[]
 for path in sorted((ROOT/'data/SAFES_PDBs').glob('*.pdb')):
  if path.stem.endswith('_af3'):continue
  exp,_=load_model(path);pred,_=load_model(path.with_name(path.stem+'_af3.pdb'));chains=[c.id for c in exp if protein_residues(c)[1]]
  row=dict(pdb_id=path.stem,experimental_chains=';'.join(chains),af3_chains=';'.join(c.id for c in pred),n_experimental_protein_chains=len(chains),intended_interface_metadata='not_available_in_local_files',status='pending_binding_partner_metadata',notes='No automatic all-chain-pair enumeration; multichain cases await intended partner chain groups.')
  if len(chains)==2:
   try:
    mapping,details,ties,rmsd=map_chains(exp,pred,chains)
    mapping_rows.extend(dict(pdb_id=path.stem,**v) for v in details)
    row.update(status='provisional_two_chain_complex',mapped_af3_pair=';'.join(mapping[c] for c in chains),notes='Only one protein chain pair in supplied complex; sequence/geometry verified. Intended SAFES binding-partner metadata still unverified.',n_tied_mappings=len(ties))
    for kind,model,a,b in [('experimental',exp,*chains),('AF3',pred,*(mapping[c] for c in chains))]:
     m=metrics(model[a],model[b]);m={k:v for k,v in m.items() if 'residue_ids' not in k}
     rows.append(standard(dict(dataset='SAFES_provisional_two_chain',structure_type=kind,interface_id=path.stem,structure_id=path.stem,chain1=a,chain2=b,source=str(path.relative_to(ROOT)),severe_overlap=m['min_interchain_distance_A']<.5,**m)))
   except Exception as exc:row.update(status='mapping_failed',notes=str(exc))
  audit.append(row)
 pd.DataFrame(audit).to_csv(OUT/'safes_interface_selection_audit.csv',index=False);pd.DataFrame(mapping_rows).to_csv(OUT/'safes_chain_mapping.csv',index=False)
 d=pd.DataFrame(rows);d.to_csv(OUT/'safes_provisional_two_chain_metrics.csv',index=False);return d

def main():
 assert (OUT/'FOLDBENCH_MATCHING_AUDIT.md').exists(),'Phase 1 must complete first'
 f=foldbench();p=pedro();r=reference();s=safes();summaries=[]
 paired=f[f.paired_success];summaries+=paired_analysis(paired[paired.structure_type=='experimental'],paired[paired.structure_type=='AF3'],'foldbench')
 stable=paired[~paired.mapping_metric_sensitive]
 summaries+=paired_analysis(stable[stable.structure_type=='experimental'],stable[stable.structure_type=='AF3'],'foldbench_mapping_stable')
 if len(s):summaries+=paired_analysis(s[s.structure_type=='experimental'],s[s.structure_type=='AF3'],'safes_provisional_two_chain')
 refmain=r[r.reference_taxonomy_scope!='neither_Mgen']
 summaries+=distribution(p,refmain,['Pedro AF3 raw','Mgen-containing reference candidates'],'pedro_reference','pedro_vs_reference')
 summaries+=distribution(p,p[p.qc_included],['Pedro raw','Pedro severe-overlap QC'],'qc','pedro_raw_vs_qc')
 summaries+=distribution(p,r,['Pedro AF3 raw','All 29 reference candidates'],'pedro_reference','pedro_vs_all29_reference_sensitivity')
 pd.DataFrame(summaries).to_csv(OUT/'nonzero_distribution_summary.csv',index=False)
 populations={'FoldBench experimental all':f[f.structure_type=='experimental'],'FoldBench experimental paired':paired[paired.structure_type=='experimental'],'FoldBench AF3 paired':paired[paired.structure_type=='AF3'],'Pedro raw':p,'Pedro QC':p[p.qc_included],'Reference all 29':r,'Reference Mgen-containing':refmain,'SAFES provisional experimental':s[s.structure_type=='experimental'],'SAFES provisional AF3':s[s.structure_type=='AF3']}
 zeros=[];allstats=[]
 for name,d in populations.items():
  zeros.append(dict(dataset=name,n_computable=len(d),n_structures=d.structure_id.nunique(),n_zero=int((d.contact_count==0).sum()),pct_zero=100*(d.contact_count==0).mean(),n_nonzero=int((d.contact_count>0).sum()),n_qc_flagged=int((~d.qc_included).sum())))
  for metric in METRICS:
   for scope,x in [('all_including_zeros',d[metric]),('nonzero_interfaces',d.loc[d.contact_count>0,metric])]:allstats.append(dict(dataset=name,metric=metric,scope=scope,**desc(x)))
 pd.DataFrame(zeros).to_csv(OUT/'zero_contact_summary.csv',index=False);pd.DataFrame(allstats).to_csv(OUT/'metric_summary.csv',index=False)
 impacts=[]
 for metric in METRICS:
  a=p.loc[p.nonzero,metric];b=p.loc[p.nonzero&p.qc_included,metric]
  impacts.append(dict(metric=metric,n_raw=len(a),n_qc=len(b),n_removed=len(a)-len(b),median_raw=a.median(),median_qc=b.median(),mean_raw=a.mean(),mean_qc=b.mean(),max_raw=a.max(),max_qc=b.max(),ecdf_max_gap=ks_2samp(a,b).statistic))
 pd.DataFrame(impacts).to_csv(OUT/'pedro_qc_effect.csv',index=False)
 # Save raw/QC subsets explicitly even if equal; avoid duplicate equal figures.
 for name,d in [('foldbench',f),('reference_candidates',r),('safes_provisional_two_chain',s)]:
  d.to_csv(OUT/f'{name}_raw.csv',index=False);d[d.qc_included].to_csv(OUT/f'{name}_qc.csv',index=False)
 assert (p.interface_fraction.between(0,1)).all() and (f.interface_fraction.between(0,1)).all()
 assert (p.total_interface_residues<=p.chain1_length+p.chain2_length).all()
 assert (f[f.structure_type=='AF3'].contact_count==0).sum()==1
 (OUT/'validation.json').write_text(json.dumps(dict(status='passed',foldbench_audited=279,foldbench_paired=132,pedro_reused_rows=len(p),reference_candidate_pairs=len(r),safes_provisional_paired=int(len(s)/2),safes_full_intended_benchmark_complete=False,checks=['missing AF3 never imputed zero','same FoldBench matched population on both sides','valid paired zero retained','zero counts separate from QC','existing residue results reused','weighted interface fractions in [0,1]']),indent=2)+'\n')
 print(pd.DataFrame(zeros).to_string(index=False));print(pd.DataFrame(impacts).to_string(index=False))
if __name__=='__main__':main()
