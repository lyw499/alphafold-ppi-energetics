"""Streaming contact audit and contextual benchmark report; new Pedro outputs only."""
import sys,json,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/local_interface_geometry/pedro';sys.path.insert(0,str(OUT/'_deps'))
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import pearsonr,spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from local_geometry_common import BIN_NAMES,CHEM_CLASSES
from parse_pedro_local_geometry import atomic_json,fingerprint
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
EDGES=np.arange(0,8.0001,.25);THRESH=[.5,1.,1.5,2.]
LABELS=['FoldBench Experimental','FoldBench AF3','SAFES Experimental','SAFES AF3','Pedro AF3','Pedro AF3 existing-QC']
def save(fig,name):
 fig.tight_layout();fig.savefig(OUT/'figures'/name,dpi=180);plt.close(fig)
def diststats(x,label):
 return dict(group=label,n_contacts=len(x),mean=float(np.mean(x)),median=float(np.median(x)),q25=float(np.quantile(x,.25)),q75=float(np.quantile(x,.75)),**{f'n_lt_{t}':int(np.sum(x<t)) for t in THRESH},**{f'fraction_lt_{t}':float(np.mean(x<t)) for t in THRESH})
def main():
 expected=pd.read_csv(OUT/'job_manifest.csv');keys=expected.job_key.tolist()
 missing=[k for k in keys if not (OUT/f'checkpoints/{k}.done.json').exists()]
 if missing:raise RuntimeError(f'Incomplete parse: {len(missing)} jobs; no scientific report produced')
 for k in keys:assert json.loads((OUT/f'checkpoints/{k}.done.json').read_text())['protocol']==fingerprint()
 s=pd.concat([pd.read_parquet(OUT/f'summary_parts/{k}.parquet') for k in keys],ignore_index=True)
 assert len(s)==170644 and s.interface_id.is_unique and s.status.eq('passed').all()
 assert s.n_contacts.sum()==6671664 and s.n_contacts.eq(0).sum()==116535 and s.existing_severe_overlap.sum()==464
 s['fraction_lt_2']=np.where(s.n_contacts>0,s.n_lt_2_0/s.n_contacts,np.nan)
 s.to_parquet(OUT/'pedro_interface_geometry_summary.parquet',index=False,compression='zstd');s.to_csv(OUT/'pedro_interface_geometry_summary.csv',index=False)
 audits=pd.concat([pd.read_csv(OUT/f'audit_parts/{k}.csv') for k in keys],ignore_index=True);audits.to_csv(OUT/'pedro_processing_audit.csv',index=False)
 n=int(s.n_contacts.sum());nq=int(s.loc[~s.existing_severe_overlap,'n_contacts'].sum())
 raw=np.lib.format.open_memmap(OUT/'qc/pooled_distances.npy',mode='w+',dtype='float64',shape=(n,))
 qc=np.lib.format.open_memmap(OUT/'qc/qc_pooled_distances.npy',mode='w+',dtype='float64',shape=(nq,))
 ir=iq=0;extreme=None;tailtypes=collections.Counter();tailclass=collections.Counter()
 for k in keys:
  f=OUT/f'contacts/contacts_{k}.parquet'
  if not f.exists():continue
  for batch in pq.ParquetFile(f).iter_batches(batch_size=50000):
   d=batch.to_pandas();v=d.d_min.to_numpy();raw[ir:ir+len(v)]=v;ir+=len(v)
   w=v[~d.existing_severe_overlap.to_numpy()];qc[iq:iq+len(w)]=w;iq+=len(w)
   short=d[d.d_min<2].copy()
   if len(short):
    extreme=short if extreme is None else pd.concat([extreme,short],ignore_index=True)
    extreme=extreme.nsmallest(10000,'d_min')
    for row in short.itertuples():
     atoms=tuple(sorted([f'{row.residue_a_type}:{row.atom_a_name}',f'{row.residue_b_type}:{row.atom_b_name}']))
     for t in THRESH:
      if row.d_min<t:
       tailtypes[t,atoms[0],atoms[1]]+=1;tailclass[t,row.contact_class]+=1
 assert ir==n and iq==nq
 if extreme is not None:extreme.sort_values('d_min').to_csv(OUT/'pedro_short_contact_audit.csv',index=False)
 pd.DataFrame([dict(threshold_A=t,atom_residue_a=a,atom_residue_b=b,n_contacts=c) for (t,a,b),c in tailtypes.items()]).sort_values(['threshold_A','n_contacts'],ascending=[True,False]).to_csv(OUT/'qc/short_contact_atom_residue_types.csv',index=False)
 pd.DataFrame([dict(threshold_A=t,contact_class=a,n_contacts=c) for (t,a),c in tailclass.items()]).to_csv(OUT/'qc/short_contact_chemical_classes.csv',index=False)
 values={};interfaces={}
 for ds,display in [('foldbench','FoldBench'),('safes','SAFES')]:
  contacts=pd.read_csv(OUT.parent/ds/'residue_contacts.csv',usecols=['structure_type','d_min','contact_class'])
  summary=pd.read_csv(OUT.parent/ds/'interface_summaries.csv')
  for typ in ['Experimental','AF3']:
   label=f'{display} {typ}';values[label]=contacts.loc[contacts.structure_type==typ,'d_min'].to_numpy();interfaces[label]=summary[summary.structure_type==typ].copy()
 values['Pedro AF3']=raw;values['Pedro AF3 existing-QC']=qc
 interfaces['Pedro AF3']=s;interfaces['Pedro AF3 existing-QC']=s[~s.existing_severe_overlap].copy()
 pooled=pd.DataFrame([diststats(values[k],k) for k in LABELS]);pooled.to_csv(OUT/'pooled_distance_comparison.csv',index=False)
 rows=[];chem=[];qcrows=[];bins=[]
 for label,g in interfaces.items():
  nz=g[g.n_contacts>0]
  rows.append(dict(group=label,n_interfaces=len(g),n_zero=int(g.n_contacts.eq(0).sum()),zero_percent=float(g.n_contacts.eq(0).mean()*100),nonzero_interfaces=len(nz),median_interface_median=float(nz.median_d_min.median()),mean_interface_median=float(nz.median_d_min.mean()),mean_interface_fraction_lt_2=float((nz.n_lt_2_0/nz.n_contacts).mean())))
  for cls in CHEM_CLASSES:
   count=float((g['chem_'+cls].fillna(0)*g.n_contacts).sum());chem.append(dict(group=label,contact_class=cls,n_contacts=int(round(count)),pooled_fraction=count/g.n_contacts.sum(),mean_interface_fraction=float(nz['chem_'+cls].mean())))
  for b in BIN_NAMES:bins.append(dict(group=label,bin=b,pooled_fraction=float((g[b].fillna(0)*g.n_contacts).sum()/g.n_contacts.sum()),mean_interface_fraction=float(nz[b].mean())))
  for t in THRESH:
   col='n_lt_'+str(t).replace('.','_');affected=g[col]>0
   qcrows.append(dict(group=label,threshold_A=t,n_contacts=int(g[col].sum()),fraction_contacts=float(g[col].sum()/g.n_contacts.sum()),n_interfaces=int(affected.sum()),fraction_all_interfaces=float(affected.mean()),fraction_nonzero_interfaces=float(affected.sum()/len(nz)),affected_existing_severe=int((affected & g.existing_severe_overlap).sum()) if 'existing_severe_overlap' in g else None))
 interface_summary=pd.DataFrame(rows);interface_summary.to_csv(OUT/'interface_level_comparison.csv',index=False)
 chemistry=pd.DataFrame(chem);chemistry.to_csv(OUT/'chemistry_comparison.csv',index=False)
 bins=pd.DataFrame(bins);bins.to_csv(OUT/'distance_bin_comparison.csv',index=False)
 tail=pd.DataFrame(qcrows);tail.to_csv(OUT/'qc/short_distance_qc.csv',index=False)
 concentration=s.groupby(['source_file','job_id','best_model'],dropna=False).agg(n_interfaces=('interface_id','size'),contacts=('n_contacts','sum'),n_lt_0_5=('n_lt_0_5','sum'),n_lt_1_0=('n_lt_1_0','sum'),n_lt_1_5=('n_lt_1_5','sum'),n_lt_2_0=('n_lt_2_0','sum'),severe_interfaces=('existing_severe_overlap','sum')).reset_index().sort_values('n_lt_2_0',ascending=False)
 concentration.to_csv(OUT/'qc/short_contact_structure_concentration.csv',index=False)
 s[s.n_lt_2_0>0].sort_values(['min_d_min','n_lt_2_0'],ascending=[True,False]).to_csv(OUT/'qc/short_distance_interfaces.csv',index=False)
 correlations=[]
 for version,g in [('raw',s),('existing_QC',s[~s.existing_severe_overlap])]:
  for x in ['n_contacts','total_interface_residues','interface_fraction']:
   for y in ['median_d_min','fraction_lt_2']+BIN_NAMES:
    z=g[[x,y]].dropna();valid=len(z)>2 and z[x].nunique()>1 and z[y].nunique()>1
    correlations.append(dict(version=version,x=x,y=y,n_interfaces=len(z),pearson=float(pearsonr(z[x],z[y])[0]) if valid else np.nan,spearman=float(spearmanr(z[x],z[y])[0]) if valid else np.nan))
 corr=pd.DataFrame(correlations);corr.to_csv(OUT/'interface_size_geometry_correlations.csv',index=False)
 # Plots: pooled contacts descriptive, no residue-level inferential tests.
 def hist(ax,labels):
  for label in labels:
   counts,_=np.histogram(values[label],bins=EDGES);ax.stairs(counts/(len(values[label])*.25),EDGES,label=label,linewidth=1.6)
  ax.set(xlabel='Residue-pair minimum heavy-atom distance (Å)',ylabel='Density (Å⁻¹)',xlim=(0,8));ax.legend(fontsize=8)
 def ecdf(ax,labels):
  for label in labels:
   v=values[label]
   # In-place memmap sort uses bounded RAM; plot selected exact empirical ranks.
   v.sort();ii=np.unique(np.linspace(0,len(v)-1,min(10000,len(v)),dtype=int));ax.plot(v[ii],(ii+1)/len(v),label=label,lw=1.5)
  ax.set(xlabel='Residue-pair minimum heavy-atom distance (Å)',ylabel='Empirical cumulative fraction',xlim=(0,8),ylim=(0,1));ax.legend(fontsize=8)
 fig,ax=plt.subplots(figsize=(7,4));hist(ax,['Pedro AF3']);save(fig,'pedro_min_distance_distribution.png')
 fig,ax=plt.subplots(figsize=(7,4));ecdf(ax,['Pedro AF3']);save(fig,'pedro_min_distance_ecdf.png')
 fig,axes=plt.subplots(1,2,figsize=(13,4),sharey=True)
 for ax,ls in zip(axes,[LABELS[:2]+['Pedro AF3'],LABELS[2:4]+['Pedro AF3']]):hist(ax,ls)
 save(fig,'benchmark_min_distance_comparison.png')
 fig,axes=plt.subplots(1,2,figsize=(13,4),sharey=True)
 for ax,ls in zip(axes,[LABELS[:2]+['Pedro AF3'],LABELS[2:4]+['Pedro AF3']]):ecdf(ax,ls)
 save(fig,'benchmark_min_distance_ecdf.png')
 fig,axes=plt.subplots(1,2,figsize=(13,4),sharey=True)
 for ax,ls in zip(axes,[LABELS[:2]+['Pedro AF3'],LABELS[2:4]+['Pedro AF3']]):
  for label in ls:
   v=interfaces[label].median_d_min.dropna();ax.hist(v,bins=EDGES,density=True,histtype='step',label=label,lw=1.6)
  ax.set(xlabel='Within-interface median distance (Å)',ylabel='Density (Å⁻¹)',xlim=(0,8));ax.legend(fontsize=8)
 save(fig,'interface_median_distance_comparison.png')
 fig,ax=plt.subplots(figsize=(7,4));hist(ax,['Pedro AF3','Pedro AF3 existing-QC']);save(fig,'pedro_raw_vs_existing_qc_distance.png')
 fig,ax=plt.subplots(figsize=(8,4))
 for label in LABELS:
  g=tail[tail.group==label];ax.plot(g.threshold_A,g.fraction_contacts*100,marker='o',label=label)
 ax.set(xlabel='Strict distance threshold (Å)',ylabel='Residue contacts below threshold (%)',yscale='symlog');ax.legend(fontsize=7);save(fig,'extreme_short_contact_comparison.png')
 for y,name in [('median_d_min','contact_count_vs_median_distance.png'),('fraction_lt_2','contact_count_vs_short_fraction.png')]:
  fig,axes=plt.subplots(1,2,figsize=(12,4),sharey=True)
  for ax,(label,g) in zip(axes,[('Raw',s),('Existing QC',s[~s.existing_severe_overlap])]):
   g=g[g.n_contacts>0];h=ax.hexbin(np.log10(g.n_contacts),g[y],gridsize=50,mincnt=1,bins='log');fig.colorbar(h,ax=ax,label='Interface count');ax.set(xlabel='log10 residue contact count',ylabel=y,title=label)
  save(fig,name)
 for frame,col,file in [(chemistry,'contact_class','chemical_contact_class_comparison.png'),(bins,'bin','distance_bin_fraction_comparison.png')]:
  pivot=frame.pivot(index=col,columns='group',values='pooled_fraction').reindex(columns=LABELS)
  fig,ax=plt.subplots(figsize=(12,5));pivot.plot.bar(ax=ax);ax.set(ylabel='Pooled residue-contact fraction',xlabel='');ax.tick_params(axis='x',labelrotation=30);ax.legend(fontsize=7);save(fig,file)
 rawrow=pooled[pooled.group=='Pedro AF3'].iloc[0];qcrow=pooled[pooled.group=='Pedro AF3 existing-QC'].iloc[0]
 def table(frame):return '```\n'+frame.to_string(index=False)+'\n```\n'
 md='# Pedro local interface geometry\n\nAll expected observations were reconciled before interpretation. No old result was overwritten.\n\n'
 md+=f'Expected/successful interfaces: {len(s):,}/{len(s):,}; failed: 0; skipped scientifically: 0; zero: {int(s.n_contacts.eq(0).sum()):,}; nonzero: {int(s.n_contacts.gt(0).sum()):,}. Archived job/model observations: {len(keys)}. Residue contacts: {n:,}. Existing severe-overlap interfaces: {int(s.existing_severe_overlap.sum())}.\n\n'
 md+='## Pooled descriptive results\n'+table(pooled)+'\n## Interface-level results\n'+table(interface_summary)+'\n## Short-distance diagnostics\n'+table(tail)
 md+='\n## Size versus geometry (no significance tests)\n'+table(corr)
 md+='\n## Explicit interpretation\n'
 md+=f'1. Pedro pooled median is {rawrow["median"]:.4f} Å; benchmark medians and quartiles above quantify similarity. Compare full density and ECDF shapes, not only medians.\n'
 md+='2. The short-distance table gives exact Pedro and benchmark contact fractions below each threshold; these are descriptive comparisons, not independent-sample tests.\n'
 for t in THRESH:
  row=tail[(tail.group=='Pedro AF3') & (tail.threshold_A==t)].iloc[0]
  md+=f'   Below {t:g} Å: {int(row.n_interfaces)} Pedro interfaces ({100*row.fraction_all_interfaces:.4f}% of all interfaces); {int(row.affected_existing_severe)} of these carry the existing severe-overlap flag.\n'
 md+=f'3. Structure concentration: the ten archived jobs with most <2 Å contacts contribute {100*concentration.head(10).n_lt_2_0.sum()/max(1,concentration.n_lt_2_0.sum()):.2f}% of the short tail. Full per-job and atom/residue tables are in qc/.\n'
 md+=f'4. Existing QC changes the pooled median from {rawrow["median"]:.4f} to {qcrow["median"]:.4f} Å and the <2 Å fraction from {rawrow["fraction_lt_2.0"]:.6f} to {qcrow["fraction_lt_2.0"]:.6f}. The raw/QC figure retains identical bins.\n'
 md+='5. Broad chemistry frequencies, including equal-interface-weighted frequencies, are in chemistry_comparison.csv and the chemistry figure.\n'
 md+='6. Large-interface geometry is evaluated by raw/QC hexbin plots and correlations above; contact count is not a substitute for a local-distance metric. Zero interfaces have no defined local distribution and are excluded only from these geometry associations.\n'
 md+='7. Size and geometry measure distinct properties; empirical associations above quantify how strongly they co-vary in this collection. Neither independence nor biological validity is established by a similar pooled spectrum.\n'
 md+='8. Limitations: Pedro is not experimentally paired here; its modeled complexes, repeated protein pairs, sequence composition, chain sizes and selection differ from benchmarks. Contacts cluster within interfaces and structures. Local plausibility cannot demonstrate a correct PPI, correct chain assignment, biological assembly, or model confidence calibration. Short contacts are diagnostics, not automatically verified steric clashes.\n'
 md+='\n## Output conventions\nZero interfaces remain in summaries, with undefined contact-distance statistics NA. QC excludes only the pre-existing <0.5 Å severe-overlap classification. The ranked contact audit retains the 10,000 smallest distances below 2 Å (or all if fewer); aggregate tail tables count every contact. ECDF curves display up to 10,000 exact empirical ranks per distribution. Histograms use 0.25 Å bins. No residue-level p-values are reported. Chemistry is exploratory.\n'
 (OUT/'PEDRO_SUMMARY.md').write_text(md)
 atomic_json(OUT/'report_complete.json',dict(status='complete',interfaces=len(s),contacts=n,zero=int(s.n_contacts.eq(0).sum()),jobs=len(keys)))
 print(pooled.to_string(index=False));print(interface_summary.to_string(index=False));print('Report complete',flush=True)
if __name__=='__main__':main()
