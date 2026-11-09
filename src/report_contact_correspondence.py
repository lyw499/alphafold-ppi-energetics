"""Contact network correspondence from saved paired distances only. No structures opened."""
import json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from local_geometry_common import chemistry,contact_class,CHEM_CLASSES
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/contact_correspondence';SOURCE=ROOT/'results/local_interface_geometry'
STATES=['shared','Experimental_only','AF3_only'];BINS=['(8,9]','(9,10]','(10,12]','(12,15]','>15']
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
def save(fig,name):
 fig.tight_layout();fig.savefig(OUT/'figures'/name,dpi=180);plt.close(fig)
def ratio(a,b):return a/b if b else np.nan
def metrics(d,t):
 e=d.experimental_d_min.le(t);a=d.af3_d_min.le(t);shared=int((e&a).sum());lost=int((e&~a).sum());gain=int((~e&a).sum())
 return dict(experimental_contacts=int(e.sum()),af3_contacts=int(a.sum()),shared=shared,Experimental_only=lost,AF3_only=gain,contact_count_difference=gain-lost,recovery=ratio(shared,shared+lost),precision=ratio(shared,shared+gain),jaccard=ratio(shared,shared+lost+gain))
def texttable(d):return '\n```\n'+d.to_string(index=False)+'\n```\n'
def main():
 hashes={};pooled=[];perstats=[];tails=[];deltas=[];correlations=[];chemrows=[];coverage=[];cutoffdist=[]
 for ds in ['foldbench','safes']:
  dest=OUT/ds
  inputs=['paired_residue_distances.csv','chain_pair_audit.csv','interface_summaries.csv','residue_mapping_audit.csv']
  for name in inputs:
   p=SOURCE/ds/name;hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
  d=pd.read_csv(SOURCE/ds/inputs[0],dtype={'interface_id':str},low_memory=False)
  audit=pd.read_csv(SOURCE/ds/inputs[1],dtype={'interface_id':str})
  summary=pd.read_csv(SOURCE/ds/inputs[2],dtype={'interface_id':str})
  reliable=audit[audit.residue_pair_mapping_status.eq('stable_identical_sequence_positions_only')].copy()
  assert reliable.interface_id.is_unique and set(d.interface_id)<=set(reliable.interface_id)
  assert np.isfinite(d[['experimental_d_min','af3_d_min']]).all().all()
  assert ((d.experimental_d_min<=8)|(d.af3_d_min<=8)).all()
  keys=['interface_id','experimental_residue_a','experimental_residue_b','af3_residue_a','af3_residue_b'];assert not d.duplicated(keys).any()
  d['classification']=np.select([(d.experimental_d_min<=8)&(d.af3_d_min<=8),d.experimental_d_min<=8],['shared','Experimental_only'],default='AF3_only')
  assert d.contact_state.map({'both':'shared','lost_in_AF3':'Experimental_only','gained_in_AF3':'AF3_only'}).eq(d.classification).all()
  d['delta_d']=d.af3_d_min-d.experimental_d_min
  assert np.allclose(d.delta_d,d.delta_d_min)
  d['contact_class']=[contact_class(chemistry(a),chemistry(b)) for a,b in zip(d.residue_a_type,d.residue_b_type)]
  d.to_csv(dest/'paired_contact_classification.csv',index=False)
  audit['correspondence_included']=audit.interface_id.isin(reliable.interface_id)
  audit.to_csv(dest/'mapping_inclusion_audit.csv',index=False)
  groups=dict(tuple(d.groupby('interface_id')));empty=d.iloc[:0];rows=[]
  for r in reliable.to_dict('records'):
   g=groups.get(r['interface_id'],empty)
   m=metrics(g,8)
   assert m['experimental_contacts']==r['n_experimental_contacts_residue_mappable'] and m['af3_contacts']==r['n_af3_contacts_residue_mappable'] and len(g)==r['n_union_residue_pairs']
   for cutoff in [5,6,8]:
    rows.append(dict(dataset=ds,interface_id=r['interface_id'],structure_id=r['structure_id'],cutoff_A=cutoff,**metrics(g,cutoff),full_experimental_contacts_8A=r['experimental_n_contacts'],full_af3_contacts_8A=r['af3_n_contacts'],legacy_mapping_stable=r['legacy_mapping_stable'],local_distance_spectrum_stable=r['local_distance_spectrum_stable'],existing_qc_flags=r.get('existing_qc_flags',''),experimental_mapping_coverage=ratio(m['experimental_contacts'],r['experimental_n_contacts']),af3_mapping_coverage=ratio(m['af3_contacts'],r['af3_n_contacts'])))
  per=pd.DataFrame(rows);per.to_csv(dest/'cutoff_sensitivity_per_interface.csv',index=False);p=per[per.cutoff_A==8].copy();p.to_csv(dest/'interface_contact_overlap.csv',index=False)
  for cutoff in [5,6,8]:
   pooled.append(dict(dataset=ds,cutoff_A=cutoff,n_interfaces=len(reliable),**metrics(d,cutoff)))
   z=per[per.cutoff_A==cutoff]
   cutoffdist.append(dict(dataset=ds,cutoff_A=cutoff,n_af3_more=int(z.contact_count_difference.gt(0).sum()),n_equal=int(z.contact_count_difference.eq(0).sum()),n_af3_fewer=int(z.contact_count_difference.lt(0).sum()),median_count_difference=z.contact_count_difference.median()))
   for col in ['recovery','precision','jaccard','contact_count_difference']:
    v=z[col].dropna();perstats.append(dict(dataset=ds,cutoff_A=cutoff,metric=col,n_defined=len(v),mean=v.mean(),median=v.median(),q25=v.quantile(.25),q75=v.quantile(.75)))
  coverage.append(dict(dataset=ds,matched_interfaces=len(summary)//2,reliable_interfaces=len(reliable),zero_union_interfaces=int((p.shared+p.Experimental_only+p.AF3_only).eq(0).sum()),full_exp_contacts_selected=int(p.full_experimental_contacts_8A.sum()),mappable_exp_contacts=int(p.experimental_contacts.sum()),full_af3_contacts_selected=int(p.full_af3_contacts_8A.sum()),mappable_af3_contacts=int(p.af3_contacts.sum()),unavailable_or_failed=int((~audit.included_primary).sum()),ambiguous_excluded=int(audit.residue_pair_mapping_status.eq('skipped_ambiguous_chain_correspondence').sum())))
  for state,other in [('AF3_only','experimental_d_min'),('Experimental_only','af3_d_min')]:
   v=d.loc[d.classification.eq(state),other];counts=np.bincount(np.searchsorted([9,10,12,15],v,side='left'),minlength=5)
   tails.append(dict(dataset=ds,state=state,n=len(v),median=v.median(),q25=v.quantile(.25),q75=v.quantile(.75),**{f'fraction_{b}':ratio(int(c),len(v)) for b,c in zip(BINS,counts)}))
   for iid,g in d[d.classification.eq(state)].groupby('interface_id'):
    vv=g[other];cc=np.bincount(np.searchsorted([9,10,12,15],vv,side='left'),minlength=5)
    # Interfaces without discordant contacts have undefined tail fractions and are absent here.
    row=dict(interface_id=iid,state=state,n=len(vv),median=vv.median(),q25=vv.quantile(.25),q75=vv.quantile(.75),**{f'fraction_{b}':c/len(vv) for b,c in zip(BINS,cc)})
    tail_interface_rows.append(dict(dataset=ds,**row))
  for state in STATES:
   g=d[d.classification.eq(state)];v=g.delta_d
   deltas.append(dict(dataset=ds,state=state,n=len(v),mean=v.mean(),median=v.median(),q25=v.quantile(.25),q75=v.quantile(.75),fraction_shorter=(v<0).mean()))
   gc=g.groupby('interface_id').delta_d.agg(['size','mean','median']).reset_index();gc['state']=state;gc['dataset']=ds;delta_interface_frames.append(gc)
   for cls in CHEM_CLASSES:
    perfrac=g.assign(hit=g.contact_class.eq(cls)).groupby('interface_id').hit.mean()
    chemrows.append(dict(dataset=ds,state=state,contact_class=cls,n_contacts=int(g.contact_class.eq(cls).sum()),pooled_fraction=g.contact_class.eq(cls).mean(),mean_interface_fraction=perfrac.mean()))
  for x in ['AF3_only','jaccard','experimental_contacts','full_experimental_contacts_8A']:
   z=p[[x,'contact_count_difference']].dropna();valid=len(z)>2 and z[x].nunique()>1 and z.contact_count_difference.nunique()>1
   correlations.append(dict(dataset=ds,x=x,n=len(z),pearson=pearsonr(z[x],z.contact_count_difference)[0] if valid else np.nan,spearman=spearmanr(z[x],z.contact_count_difference)[0] if valid else np.nan))
  # QC and mapping are descriptive strata, not new filters.
  for col in ['existing_qc_flags','legacy_mapping_stable','local_distance_spectrum_stable']:
   agg=p.groupby(col,dropna=False).agg(n_interfaces=('interface_id','size'),median_excess=('contact_count_difference','median'),mean_excess=('contact_count_difference','mean'),median_jaccard=('jaccard','median')).reset_index()
   agg.to_csv(dest/f'stratified_{col}.csv',index=False)
  fig,ax=plt.subplots(figsize=(7,4))
  for col in ['recovery','precision','jaccard']:
   v=np.sort(p[col].dropna());ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',label=col)
  ax.set(xlabel='Per-interface contact overlap metric',ylabel='ECDF',xlim=(0,1),title=ds);ax.legend();save(fig,f'{ds}_interface_overlap_ecdf.png')
  fig,axes=plt.subplots(1,2,figsize=(11,4))
  for state,other in [('AF3_only','experimental_d_min'),('Experimental_only','af3_d_min')]:
   v=np.sort(d.loc[d.classification.eq(state),other]);axes[0].step(v,np.arange(1,len(v)+1)/len(v),where='post',label=state)
   c=np.bincount(np.searchsorted([9,10,12,15],v,side='left'),minlength=5);axes[1].plot(BINS,c/len(v),marker='o',label=state)
  axes[0].set(xlabel='Distance in non-contact structure (Å)',ylabel='ECDF',xscale='log');axes[1].set(xlabel='Distance in non-contact structure (Å)',ylabel='Fraction');axes[0].legend();axes[1].legend();fig.suptitle(ds);save(fig,f'{ds}_discordant_contact_distances.png')
  fig,ax=plt.subplots(figsize=(7,4))
  for state in STATES:
   v=np.sort(d.loc[d.classification.eq(state),'delta_d']);ax.step(v,np.arange(1,len(v)+1)/len(v),where='post',label=state)
  ax.set(xlabel='AF3 − Experimental distance (Å)',ylabel='ECDF',xscale='symlog');ax.axvline(0,color='gray',ls='--');ax.legend();save(fig,f'{ds}_distance_change_ecdf.png')
  fig,axes=plt.subplots(1,3,figsize=(13,4))
  for ax,x in zip(axes,['AF3_only','jaccard','experimental_contacts']):
   ax.scatter(p[x],p.contact_count_difference,s=15,alpha=.65);ax.axhline(0,color='gray',ls='--');ax.set(xlabel=x,ylabel='AF3 − Experimental contact count')
  save(fig,f'{ds}_contact_excess_associations.png')
 pd.DataFrame(pooled).to_csv(OUT/'pooled_overlap_cutoff_sensitivity.csv',index=False)
 pd.DataFrame(perstats).to_csv(OUT/'interface_overlap_distribution_summary.csv',index=False)
 pd.DataFrame(cutoffdist).to_csv(OUT/'cutoff_interface_direction_counts.csv',index=False)
 pd.DataFrame(tails).to_csv(OUT/'discordant_distance_summary.csv',index=False)
 ti=pd.DataFrame(tail_interface_rows);ti.to_csv(OUT/'discordant_distance_per_interface.csv',index=False)
 pd.DataFrame(deltas).to_csv(OUT/'distance_change_summary.csv',index=False)
 di=pd.concat(delta_interface_frames,ignore_index=True);di.to_csv(OUT/'distance_change_per_interface.csv',index=False)
 pd.DataFrame(correlations).to_csv(OUT/'contact_excess_correlations.csv',index=False)
 pd.DataFrame(coverage).to_csv(OUT/'mapping_coverage.csv',index=False)
 chem=pd.DataFrame(chemrows);chem.to_csv(OUT/'chemistry_by_contact_state.csv',index=False)
 for ds in ['foldbench','safes']:
  fig,axes=plt.subplots(1,2,figsize=(13,5),sharey=True)
  for ax,col in zip(axes,['pooled_fraction','mean_interface_fraction']):
   chem[chem.dataset==ds].pivot(index='contact_class',columns='state',values=col).plot.bar(ax=ax);ax.set(ylabel='Fraction',xlabel='',title=col);ax.tick_params(axis='x',rotation=35)
  save(fig,f'{ds}_contact_state_chemistry.png')
 for path,h in hashes.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
 (OUT/'provenance.json').write_text(json.dumps(dict(input_sha256=hashes,source_coordinates_read=False,mapping_recomputed=False,Pedro_accessed=False,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
 md='# 接触网络对应分析：FoldBench 与 SAFES\n\n仅复用已保存的配对残基距离和已验证映射，未读取任何结构或 Pedro 归档；旧结果未改动。以任一结构在 8 Å 内接触的可对应残基对之并集为分析全集。5/6 Å 汇总直接复用这组距离，故不会遗漏任何更小截断下的接触。\n'
 md+='\n## 覆盖范围与分母\n'+texttable(pd.DataFrame(coverage))
 md+='FoldBench 的链对应有歧义界面不强制配对；不可用 AF3 不当作零。SAFES 的可靠但双侧零接触链对保留，recovery/precision/Jaccard 分母为零时记 NA。残基映射仅纳入所有最优序列比对共同支持的相同氨基酸位置。全界面接触数与可对应子集接触数分列保存；不能把不可对应接触当作丢失。\n'
 md+='\n## 接触恢复、精确率与截断敏感性\n'+texttable(pd.DataFrame(pooled))+texttable(pd.DataFrame(perstats))+texttable(pd.DataFrame(cutoffdist))
 md+='\n## 新增与丢失接触：另一侧距离\n分箱严格为 (8,9]、(9,10]、(10,12]、(12,15]、>15 Å。\n'+texttable(pd.DataFrame(tails))
 md+='\n逐界面汇总（仅该类接触存在的界面；无该类接触的界面没有定义的条件比例）：\n'+texttable(ti.groupby(['dataset','state']).agg(n_interfaces=('interface_id','size'),median_of_medians=('median','median'),mean_fraction_within9=('fraction_(8,9]','mean'),mean_fraction_beyond15=('fraction_>15','mean')).reset_index())
 md+='\n## 距离变化\n'+texttable(pd.DataFrame(deltas))+'\n共享接触的逐界面距离变化：\n'+texttable(di[di.state=='shared'].groupby('dataset').agg(n_interfaces=('interface_id','size'),median_of_interface_median_delta=('median','median'),mean_of_interface_mean_delta=('mean','mean')).reset_index())
 md+='\n## 接触数过量关联\n'+texttable(pd.DataFrame(correlations))+'AF3 过量 = AF3-only − Experimental-only，因此它与 AF3-only 的相关性部分来自代数关系，不能解释为独立机制。QC 标记和既有 mapping-stable 标记的分层表保留在各数据集目录；不新增过滤。\n'
 md+='\n## 化学组成（探索性）\n'+texttable(chem)+'\n分类直接复用基准的宽泛化学分类。既给出池化比例，也给出等权界面均值；不把接触数当作独立样本量。\n'
 md+='\n## 结论与限制\n局部距离直方图不记录“哪些残基彼此接触”，因此距离谱相似与接触网络改变可以同时成立。新增与丢失接触都要计入，不能只看 AF3-only。用逐界面指标作为主要解释，池化指标仅用于描述；没有残基级显著性检验。界面仍嵌套在结构内，且同源/重复复合物并非独立样本。映射可靠子集可能不代表全部原始界面。\n'
 (OUT/'CONTACT_CORRESPONDENCE_SUMMARY.md').write_text(md)
 print(pd.DataFrame(pooled).to_string(index=False));print(pd.DataFrame(tails).to_string(index=False));print(pd.DataFrame(coverage).to_string(index=False))
tail_interface_rows=[];delta_interface_frames=[]
if __name__=='__main__':main()
