"""Exact closest-atom BB/SC decomposition from existing contact tables only."""
import json,collections
from pathlib import Path
from audit_contact_type_data import ROOT,BASE,OUT,sources
import pandas as pd
import numpy as np
import pyarrow.parquet as pq
from scipy.stats import spearmanr,pearsonr
BB={'N','CA','C','O'};TYPES=['BB-BB','BB-SC','SC-SC'];BINS=['lt2','2to3','3to4','4to5','5to6','6to8'];EDGES=np.arange(0,8.00001,.25)
def classify(a,b):
 aa=a.isin(BB);bb=b.isin(BB);return np.where(aa&bb,'BB-BB',np.where(aa|bb,'BB-SC','SC-SC'))
def fullframes():
 for p,label in sources():
  d=pq.read_table(p).to_pandas() if p.suffix=='.parquet' else pd.read_csv(p,low_memory=False,dtype={'interface_id':str})
  d['interface_id']=d.interface_id.astype(str)
  if p.suffix=='.csv':
   for typ,g in d.groupby('structure_type'):yield ('FoldBench' if label=='foldbench' else 'SAFES')+' '+typ,g
  else:yield label,d

def metadata():
 frames=[]
 for ds,label in [('foldbench','FoldBench'),('safes','SAFES')]:
  d=pd.read_csv(BASE/ds/'interface_summaries.csv',dtype={'interface_id':str})
  for typ,g in d.groupby('structure_type'):
   g=g.copy();g['dataset']=label+' '+typ;frames.append(g[['dataset','interface_id','structure_id','n_contacts','fraction_3_4','existing_qc_flags','chain_mapping_unambiguous','legacy_mapping_stable']])
 d=pd.read_csv(BASE/'mgen_same_species_control/mgen_experimental_interface_summary.csv',dtype={'interface_id':str});d['dataset']='Mgen Experimental';frames.append(d[['dataset','interface_id','structure_id','n_contacts','fraction_3_4','existing_qc_flags','pdb_id']])
 for name,excluded in [('Mgen excl 8PBZ 8PC1',['8PBZ','8PC1']),('Mgen excl 8PBZ 8PC0 8PC1',['8PBZ','8PC0','8PC1'])]:
  g=d[~d.pdb_id.isin(excluded)].copy();g['dataset']=name;frames.append(g[['dataset','interface_id','structure_id','n_contacts','fraction_3_4','existing_qc_flags','pdb_id']])
 d=pd.read_parquet(BASE/'pedro/pedro_interface_geometry_summary.parquet',columns=['interface_id','job_key','protein_pair','n_contacts','fraction_3_4','existing_severe_overlap']);d['dataset']='Pedro AF3';d['structure_id']=d.job_key;frames.append(d)
 return pd.concat(frames,ignore_index=True)
def main():
 gate=json.loads((OUT/'data_audit.json').read_text());assert gate['status']=='passed' and gate['missing_required_cells']==0
 a=pd.Series(['N','CA','C','O','CB','CG','CD','CE','NZ','OD1','OE1','SG','OXT']);b=pd.Series(['CA']*4+['SG']*9)
 assert list(classify(a,b))==['BB-BB']*4+['SC-SC']*9
 assert classify(pd.Series(['CA','SG']),pd.Series(['SG','CA'])).tolist()==['BB-SC','BB-SC']
 counts=collections.Counter();oxt_deltas=collections.Counter();oxt_counts=collections.Counter();hist={};rows=[];samples=[]
 for num,(dataset,d) in enumerate(fullframes(),1):
  assert not d[['interface_id','residue_index_a','residue_index_b']].duplicated().any() if 'residue_index_a' in d else True
  d=d.copy();d['contact_type']=classify(d.atom_a_name,d.atom_b_name);v=d.d_min.to_numpy();d['distance_bin']=np.select([v<2,v<=3,v<=4,v<=5,v<=6],BINS[:5],default='6to8')
  versions=[(dataset,d)]
  if dataset=='Mgen Experimental':
   versions += [('Mgen excl 8PBZ 8PC1',d[~d.pdb_id.isin(['8PBZ','8PC1'])]),('Mgen excl 8PBZ 8PC0 8PC1',d[~d.pdb_id.isin(['8PBZ','8PC0','8PC1'])])]
  for label,g in versions:
   ox=g[g.atom_a_name.eq('OXT')|g.atom_b_name.eq('OXT')].copy();oxt_counts[label]+=len(ox)
   if len(ox):
    old=ox.contact_type;aa=ox.atom_a_name.isin(BB|{'OXT'});bb=ox.atom_b_name.isin(BB|{'OXT'});ox['alternative']=np.where(aa&bb,'BB-BB',np.where(aa|bb,'BB-SC','SC-SC'))
    for (typ,binname),n in ox.groupby(['contact_type','distance_bin']).size().items():oxt_deltas[label,typ,binname]-=int(n)
    for (typ,binname),n in ox.groupby(['alternative','distance_bin']).size().items():oxt_deltas[label,typ,binname]+=int(n)
   for (typ,binname),n in g.groupby(['contact_type','distance_bin']).size().items():counts[label,typ,binname]+=int(n)
   for typ,z in g.groupby('contact_type'):
    h=np.histogram(z.d_min,EDGES)[0];hist[label,typ]=hist.get((label,typ),np.zeros(len(EDGES)-1,dtype=np.int64))+h
   grouped=g.groupby(['interface_id','contact_type'],sort=False)
   stats=grouped.d_min.agg(n_type='size',median_d_min='median',mean_d_min='mean')
   bins=g.groupby(['interface_id','contact_type','distance_bin']).size().unstack('distance_bin',fill_value=0).reindex(columns=BINS,fill_value=0).rename(columns={b:'count_'+b for b in BINS})
   part=stats.join(bins).reset_index();part['dataset']=label;rows.append(part)
  if len(samples)<120:samples.extend(d.groupby('contact_type').head(2).assign(dataset=dataset).to_dict('records'))
  if num%200==0:print('Classified input frames',num,flush=True)
 meta=metadata();assert not meta.duplicated(['dataset','interface_id']).any()
 skeleton=meta.merge(pd.DataFrame({'contact_type':TYPES}),how='cross');parts=pd.concat(rows,ignore_index=True)
 assert not parts.duplicated(['dataset','interface_id','contact_type']).any()
 interface=skeleton.merge(parts,on=['dataset','interface_id','contact_type'],how='left',validate='one_to_one')
 for col in ['n_type']+['count_'+b for b in BINS]:interface[col]=interface[col].fillna(0).astype(int)
 for b in BINS:
  interface['fraction_all_'+b]=interface['count_'+b]/interface.n_contacts.replace(0,np.nan)
  interface['fraction_within_'+b]=interface['count_'+b]/interface.n_type.replace(0,np.nan)
 interface['type_fraction']=interface.n_type/interface.n_contacts.replace(0,np.nan)
 assert np.array_equal(interface[['count_'+b for b in BINS]].sum(axis=1),interface.n_type)
 check=interface.groupby(['dataset','interface_id']).agg(n=('n_type','sum'),old=('n_contacts','first'),close=('count_3to4','sum'),old_close_fraction=('fraction_3_4','first'))
 assert check.n.eq(check.old).all()
 nz=check[check.n>0];assert np.allclose(nz.close/nz.n,nz.old_close_fraction,atol=1e-12,rtol=0)
 interface.to_csv(OUT/'interface_contact_type_summary.csv',index=False)
 interface.to_parquet(OUT/'interface_contact_type_summary.parquet',index=False,compression='zstd')
 pd.DataFrame(samples).to_csv(OUT/'classified_contact_samples.csv',index=False)
 datasets=meta.dataset.unique();countrows=[];binrows=[];histrows=[]
 for dataset in datasets:
  total=int(meta.loc[meta.dataset==dataset,'n_contacts'].sum())
  for typ in TYPES:
   nt=sum(counts[dataset,typ,b] for b in BINS);countrows.append(dict(dataset=dataset,contact_type=typ,n_contacts=nt,total_contacts=total,fraction=nt/total,unclassified=0))
   for b in BINS:binrows.append(dict(dataset=dataset,contact_type=typ,distance_bin=b,count=counts[dataset,typ,b],fraction_all=counts[dataset,typ,b]/total,fraction_within=counts[dataset,typ,b]/nt if nt else np.nan))
   for lo,hi,n in zip(EDGES[:-1],EDGES[1:],hist.get((dataset,typ),np.zeros(32))):histrows.append(dict(dataset=dataset,contact_type=typ,left=lo,right=hi,count=int(n),within_type_density=n/nt/.25 if nt else np.nan))
 oxtrows=[]
 for (dataset,typ,b),n in counts.items():
  total=int(meta.loc[meta.dataset==dataset,'n_contacts'].sum());delta=oxt_deltas[dataset,typ,b];oxtrows.append(dict(dataset=dataset,contact_type=typ,distance_bin=b,primary_count=n,OXT_backbone_count=n+delta,change_all_fraction=delta/total,n_contacts_with_OXT=oxt_counts[dataset]))
 pd.DataFrame(oxtrows).to_csv(OUT/'OXT_classification_sensitivity.csv',index=False)
 ct=pd.DataFrame(countrows);bt=pd.DataFrame(binrows);ct.to_csv(OUT/'contact_type_counts.csv',index=False);bt.to_csv(OUT/'contact_type_distance_bins.csv',index=False);pd.DataFrame(histrows).to_csv(OUT/'contact_type_histogram.csv',index=False)
 distribution=[]
 for (dataset,typ),g in interface.groupby(['dataset','contact_type']):
  for metric in ['type_fraction','median_d_min']+['fraction_all_'+b for b in BINS]+['fraction_within_'+b for b in BINS]:
   x=g[metric].dropna();distribution.append(dict(dataset=dataset,contact_type=typ,metric=metric,n_defined=len(x),mean=x.mean(),median=x.median(),q25=x.quantile(.25),q75=x.quantile(.75)))
 pd.DataFrame(distribution).to_csv(OUT/'interface_type_distribution_summary.csv',index=False)
 # Residual comparisons use interface-level summaries, not atom observations.
 ped=interface[(interface.dataset=='Pedro AF3')&(interface.n_contacts>0)].copy();unique=ped[['interface_id','n_contacts']].drop_duplicates();unique['size_quartile'],edges=pd.qcut(unique.n_contacts,4,labels=False,duplicates='drop',retbins=True);ped=ped.merge(unique[['interface_id','size_quartile']],on='interface_id',validate='many_to_one')
 assoc=[];quartiles=[]
 for version,g in [('raw',ped),('existing_QC',ped[~ped.existing_severe_overlap.fillna(False).astype(bool)])]:
  for typ,z in g.groupby('contact_type'):
   assoc.append(dict(version=version,contact_type=typ,n_interfaces=len(z),spearman=spearmanr(z.n_contacts,z.fraction_all_3to4)[0],pearson=pearsonr(z.n_contacts,z.fraction_all_3to4)[0]))
   for (q,typ),z in g.groupby(['size_quartile','contact_type']):quartiles.append(dict(version=version,size_quartile=int(q)+1,contact_type=typ,n_interfaces=len(z),contact_min=z.n_contacts.min(),contact_max=z.n_contacts.max(),median=z.fraction_all_3to4.median(),mean=z.fraction_all_3to4.mean(),q25=z.fraction_all_3to4.quantile(.25),q75=z.fraction_all_3to4.quantile(.75),pooled=z.count_3to4.sum()/z.n_contacts.sum()))
 pd.DataFrame(assoc).to_csv(OUT/'pedro_contact_type_size_association.csv',index=False);pd.DataFrame(quartiles).drop_duplicates().to_csv(OUT/'pedro_size_quantile_summary.csv',index=False)
 pairs=[]
 for ds in ['FoldBench','SAFES']:
  for typ in TYPES:
   ex=interface[(interface.dataset==ds+' Experimental')&(interface.contact_type==typ)];af=interface[(interface.dataset==ds+' AF3')&(interface.contact_type==typ)];joined=ex.merge(af,on='interface_id',suffixes=('_exp','_af3'),validate='one_to_one')
   for scope in ['full_matched','unambiguous_mapping']:
    z=joined if scope=='full_matched' else joined[joined.chain_mapping_unambiguous_exp.fillna(False).astype(bool)]
    for metric in ['fraction_all_3to4','median_d_min']:
     x=z[[metric+'_exp',metric+'_af3']].dropna();delta=x.iloc[:,1]-x.iloc[:,0];valid=len(x)>2 and x.iloc[:,0].nunique()>1 and x.iloc[:,1].nunique()>1
     pairs.append(dict(dataset=ds,scope=scope,contact_type=typ,metric=metric,n_pairs=len(x),median_delta=delta.median(),q25_delta=delta.quantile(.25),q75_delta=delta.quantile(.75),spearman=spearmanr(x.iloc[:,0],x.iloc[:,1])[0] if valid else np.nan,pearson=pearsonr(x.iloc[:,0],x.iloc[:,1])[0] if valid else np.nan))
 pd.DataFrame(pairs).to_csv(OUT/'paired_benchmark_contact_type_summary.csv',index=False)
 # Contribution differences and symmetric composition/within-type decomposition: exact algebra, no causality.
 diffs=[]
 for ref in ['FoldBench Experimental','SAFES Experimental','Mgen Experimental','Mgen excl 8PBZ 8PC1','Mgen excl 8PBZ 8PC0 8PC1']:
  for typ in TYPES:
   rp=ct[(ct.dataset==ref)&(ct.contact_type==typ)].iloc[0].fraction;pp=ct[(ct.dataset=='Pedro AF3')&(ct.contact_type==typ)].iloc[0].fraction
   for b in BINS:
    r=bt[(bt.dataset==ref)&(bt.contact_type==typ)&(bt.distance_bin==b)].iloc[0];p=bt[(bt.dataset=='Pedro AF3')&(bt.contact_type==typ)&(bt.distance_bin==b)].iloc[0]
    diffs.append(dict(reference=ref,contact_type=typ,distance_bin=b,pedro_minus_reference_all_fraction=p.fraction_all-r.fraction_all,pedro_minus_reference_within_fraction=p.fraction_within-r.fraction_within,composition_component=(pp-rp)*(p.fraction_within+r.fraction_within)/2,within_type_component=(p.fraction_within-r.fraction_within)*(pp+rp)/2))
 pd.DataFrame(diffs).to_csv(OUT/'pedro_experimental_decomposition.csv',index=False)
 (OUT/'classification_validation.json').write_text(json.dumps(dict(status='passed',all_contacts_classified=True,unclassified=0,interface_counts_reconciled=True,close_bins_reconciled=True,synthetic_atom_rules_passed=True,size_quantile_edges=edges.tolist(),structural_reparse=False),indent=2))
 print(ct.to_string(index=False));print(bt[bt.distance_bin=='3to4'].to_string(index=False));print(pd.DataFrame(assoc).to_string(index=False),flush=True)
if __name__=='__main__':main()
