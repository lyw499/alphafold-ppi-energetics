"""Interface-level descriptive confidence analysis from audited saved tables only."""
import json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/local_interface_geometry';OUT=BASE/'pedro_confidence_analysis'
PRIMARY='best_model_pair_iptm';METRICS=[PRIMARY,'mean_pair_iptm','job_ranking_score'];OUTCOMES=['fraction_3_4','median_d_min','fraction_lt_2']
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
def save(fig,name):fig.tight_layout();fig.savefig(OUT/'figures'/name,dpi=180);plt.close(fig)
def strata(values,q):
 codes,edges=pd.qcut(values,q,labels=False,retbins=True,duplicates='drop');return codes.add(1).astype('Int64'),edges.tolist()
def correlations(g,x,y,**context):
 z=g[[x,y]].dropna();valid=len(z)>2 and z[x].nunique()>1 and z[y].nunique()>1
 return dict(**context,x=x,y=y,n_interfaces=len(z),spearman=float(spearmanr(z[x],z[y])[0]) if valid else np.nan,pearson=float(pearsonr(z[x],z[y])[0]) if valid else np.nan)
def stats(g,**context):
 nz=g[g.n_contacts>0];n=int(nz.n_contacts.sum());v=nz.fraction_3_4
 return dict(**context,n_interfaces=len(g),n_nonzero=len(nz),n_zero=int(g.n_contacts.eq(0).sum()),n_contacts=n,confidence_min=g[PRIMARY].min(),confidence_max=g[PRIMARY].max(),confidence_median=g[PRIMARY].median(),median_contact_count=g.n_contacts.median(),median_interface_fraction=g.interface_fraction.median(),pooled_close_fraction=nz.count_3_4.sum()/n if n else np.nan,interface_close_mean=v.mean(),interface_close_median=v.median(),interface_close_q25=v.quantile(.25),interface_close_q75=v.quantile(.75),median_interface_median_d_min=nz.median_d_min.median(),pooled_short_fraction=nz.n_lt_2_0.sum()/n if n else np.nan,interface_short_mean=nz.fraction_lt_2.mean(),interface_short_median=nz.fraction_lt_2.median())
def main():
 assert json.loads((OUT/'join_audit.json').read_text())['status']=='passed'
 s=pd.read_csv(OUT/'pedro_confidence_interface_table.csv',low_memory=False);nz=s[s.n_contacts>0].copy()
 bins={};nz['confidence_decile'],bins['nonzero_confidence_deciles']=strata(nz[PRIMARY],10);nz['confidence_tercile'],bins['nonzero_confidence_terciles']=strata(nz[PRIMARY],3);nz['size_quartile'],bins['nonzero_contact_count_quartiles']=strata(nz.n_contacts,4)
 s['all_interface_confidence_decile'],bins['all_interface_confidence_deciles']=strata(s[PRIMARY],10)
 # Quantile strata never split tied scores; actual numerical boundaries and occupancies are saved.
 assignment=nz[['interface_id','confidence_decile','confidence_tercile','size_quartile']]
 s=s.merge(assignment,on='interface_id',how='left',validate='one_to_one');s.to_csv(OUT/'pedro_confidence_interface_table.csv',index=False)
 (OUT/'stratum_boundaries.json').write_text(json.dumps(bins,indent=2))
 rows=[]
 for version,g in [('raw',nz),('existing_QC_passing',nz[~nz.existing_severe_overlap]),('severe_overlap',nz[nz.existing_severe_overlap])]:
  for metric in METRICS:
   for y in OUTCOMES:rows.append(correlations(g,metric,y,version=version))
 corr=pd.DataFrame(rows);corr.to_csv(OUT/'confidence_correlations.csv',index=False)
 dec=pd.DataFrame([stats(g,decile=int(k)) for k,g in nz.groupby('confidence_decile')]);dec.to_csv(OUT/'pedro_confidence_decile_summary.csv',index=False)
 terc=pd.DataFrame([stats(g,stratum=['low','middle','high'][int(k)-1]) for k,g in nz.groupby('confidence_tercile')]);terc.to_csv(OUT/'pedro_confidence_tercile_summary.csv',index=False)
 qc=pd.DataFrame([stats(g,qc_group=label) for label,g in [('all',s),('existing_QC_passing',s[~s.existing_severe_overlap]),('severe_overlap',s[s.existing_severe_overlap])]]);qc.to_csv(OUT/'pedro_qc_confidence_summary.csv',index=False)
 sizecorr=pd.DataFrame([correlations(nz,x,'fraction_3_4') for x in ['n_contacts','total_interface_residues','interface_fraction','total_chain_length']]);sizecorr.to_csv(OUT/'interface_size_correlations.csv',index=False)
 size=[];sizebins=[]
 for version,g in [('raw',nz),('existing_QC_passing',nz[~nz.existing_severe_overlap])]:
  for k,z in g.groupby('size_quartile'):
   size.append({**stats(z,version=version,size_quartile=int(k)),**correlations(z,PRIMARY,'fraction_3_4'), 'min_contacts':int(z.n_contacts.min()),'max_contacts':int(z.n_contacts.max())})
   for decile,sub in z.groupby('confidence_decile'):sizebins.append(stats(sub,version=version,size_quartile=int(k),confidence_decile=int(decile)))
 pd.DataFrame(size).to_csv(OUT/'pedro_confidence_size_stratified_summary.csv',index=False);sz=pd.DataFrame(sizebins);sz.to_csv(OUT/'confidence_deciles_within_size_quartiles.csv',index=False)
 zero=[]
 for k,g in s.groupby('all_interface_confidence_decile'):
  zero.append(dict(decile=int(k),n_interfaces=len(g),n_zero=int(g.n_contacts.eq(0).sum()),zero_fraction=g.n_contacts.eq(0).mean(),confidence_min=g[PRIMARY].min(),confidence_max=g[PRIMARY].max(),confidence_median=g[PRIMARY].median()))
 zero=pd.DataFrame(zero);zero.to_csv(OUT/'pedro_zero_contact_confidence_summary.csv',index=False)
 zd=[]
 for metric in METRICS:
  for group,g in [('zero',s[s.n_contacts==0]),('nonzero',nz)]:
   v=g[metric].dropna();zd.append(dict(metric=metric,group=group,n_interfaces=len(v),mean=v.mean(),median=v.median(),q25=v.quantile(.25),q75=v.quantile(.75),min=v.min(),max=v.max()))
 pd.DataFrame(zd).to_csv(OUT/'zero_nonzero_confidence_distributions.csv',index=False)
 # Experimental contexts are read as existing summary outputs, never recomputed from coordinates.
 comparison=[];context_files=[BASE/'mgen_same_species_control/mgen_same_species_comparison_summary.csv',BASE/'mgen_same_species_control/mgen_reference_sensitivity_summary.csv',BASE/'mgen_same_species_control/interface_level_distributions.csv',BASE/'mgen_same_species_control/mgen_experimental_interface_summary.csv']
 context_hash={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in context_files}
 c=pd.read_csv(context_files[0]);ints=pd.read_csv(context_files[2]);ms=pd.read_csv(context_files[3])
 for name in ['FoldBench Experimental','SAFES Experimental','Mgen Experimental']:
  r=c[c.dataset==name].iloc[0];t=ints[(ints.dataset==name)&(ints.metric=='fraction_3_4')].iloc[0]
  comparison.append(dict(group=name,n_interfaces=int(r.n_interfaces),pooled_close_fraction=r.fraction_3_4,interface_close_median=t['median'],interface_close_q25=t.q25,interface_close_q75=t.q75))
 sens=pd.read_csv(context_files[1])
 for name,exclude in [('exclude_8PBZ_8PC1',['8PBZ','8PC1']),('exclude_8PBZ_8PC0_8PC1',['8PBZ','8PC0','8PC1'])]:
  r=sens[(sens.view==name)&sens.weighting.eq('pooled')].iloc[0];v=ms.loc[~ms.pdb_id.isin(exclude),'fraction_3_4'].dropna()
  comparison.append(dict(group='Mgen '+name,n_interfaces=int(r.n_pairs),pooled_close_fraction=r.fraction_3_4,interface_close_median=v.median(),interface_close_q25=v.quantile(.25),interface_close_q75=v.quantile(.75)))
 comparison.append(stats(nz,group='Pedro overall nonzero'))
 for r in terc.to_dict('records'):comparison.append(dict(group='Pedro '+r['stratum']+' tercile',**r))
 top=nz[nz.confidence_decile==nz.confidence_decile.max()];comparison.append(stats(top,group='Pedro top decile'))
 comp=pd.DataFrame(comparison);comp.to_csv(OUT/'high_confidence_benchmark_comparison.csv',index=False)
 for y,file in [('fraction_3_4','01_confidence_vs_close_contact_fraction.png'),('median_d_min','03_confidence_vs_median_distance.png'),('fraction_lt_2','04_confidence_vs_short_contact_fraction.png')]:
  fig,ax=plt.subplots(figsize=(8,5));h=ax.hexbin(nz[PRIMARY],nz[y],gridsize=60,mincnt=1,bins='log');fig.colorbar(h,ax=ax,label='Interface count (log color)')
  trend=nz.groupby('confidence_decile').agg(x=(PRIMARY,'median'),median=(y,'median'),q25=(y,lambda v:v.quantile(.25)),q75=(y,lambda v:v.quantile(.75)))
  ax.plot(trend.x,trend['median'],color='red',marker='o',label='Decile median');ax.fill_between(trend.x,trend.q25,trend.q75,color='red',alpha=.18,label='Decile IQR');ax.set(xlabel='Selected-model pair ipTM',ylabel=y,title='Nonzero interfaces only');ax.legend();save(fig,file)
 fig,ax=plt.subplots(figsize=(8,4));ax.plot(dec.decile,100*dec.pooled_close_fraction,'o-',label='Pooled contacts (descriptive)');ax.plot(dec.decile,100*dec.interface_close_median,'o-',label='Interface median');ax.fill_between(dec.decile,100*dec.interface_close_q25,100*dec.interface_close_q75,alpha=.2,label='Interface IQR');ax.set(xlabel='Pair-ipTM decile among nonzero interfaces',ylabel='Contacts in (3,4] Å (%)');ax.legend();save(fig,'02_confidence_decile_close_contact_fraction.png')
 fig,ax=plt.subplots(figsize=(10,6));y=np.arange(len(comp));ax.scatter(100*comp.pooled_close_fraction,y,label='Pooled contacts',marker='s');ax.errorbar(100*comp.interface_close_median,y+.15,xerr=[100*(comp.interface_close_median-comp.interface_close_q25),100*(comp.interface_close_q75-comp.interface_close_median)],fmt='o',label='Interface median / IQR');ax.set_yticks(y);ax.set_yticklabels(comp.group);ax.invert_yaxis();ax.set(xlabel='Contacts in (3,4] Å (%)');ax.legend();save(fig,'05_high_confidence_vs_experimental_benchmarks.png')
 fig,axes=plt.subplots(2,2,figsize=(11,8))
 for ax,x in zip(axes.flat,['n_contacts','total_interface_residues','interface_fraction','total_chain_length']):
  v=nz[x] if x=='interface_fraction' else np.log10(nz[x]);h=ax.hexbin(v,nz.fraction_3_4,gridsize=45,mincnt=1,bins='log');ax.set(xlabel=x if x=='interface_fraction' else 'log10 '+x,ylabel='Fraction (3,4] Å');fig.colorbar(h,ax=ax,label='Interface count')
 save(fig,'06_interface_size_vs_close_contact_fraction.png')
 fig,axes=plt.subplots(1,2,figsize=(12,4),sharey=True)
 for ax,version in zip(axes,['raw','existing_QC_passing']):
  for k,g in sz[sz.version==version].groupby('size_quartile'):ax.plot(g.confidence_median,100*g.interface_close_median,'o-',label=f'Size quartile {k}')
  ax.set(xlabel='Selected-model pair ipTM (bin median)',ylabel='Interface median (3,4] fraction (%)',title=version);ax.legend()
 save(fig,'07_confidence_size_stratified.png')
 fig,ax=plt.subplots(figsize=(7,4));y=np.arange(len(qc));ax.bar(y-.17,100*qc.pooled_close_fraction,width=.34,label='Pooled');ax.bar(y+.17,100*qc.interface_close_median,width=.34,label='Interface median');ax.set_xticks(y);ax.set_xticklabels(qc.qc_group);ax.set(ylabel='Contacts in (3,4] Å (%)');ax.legend();save(fig,'08_qc_stratified_close_contact_fraction.png')
 fig,axes=plt.subplots(1,2,figsize=(11,4));axes[0].plot(zero.confidence_median,zero.zero_fraction,'o-');axes[0].set(xlabel='Pair ipTM: all-interface decile median',ylabel='Zero-contact fraction',ylim=(0,1))
 for label,g in [('Zero',s[s.n_contacts==0]),('Nonzero',nz)]:
  v=np.sort(g[PRIMARY]);ix=np.unique(np.linspace(0,len(v)-1,min(5000,len(v)),dtype=int));axes[1].plot(v[ix],(ix+1)/len(v),label=label)
 axes[1].set(xlabel='Selected-model pair ipTM',ylabel='ECDF');axes[1].legend();save(fig,'09_zero_contact_fraction_by_confidence.png')
 # Optional model selection audit is limited to existing scores/indices.
 s.groupby('selected_model').agg(n_interfaces=('interface_id','size'),n_nonzero=('geometry_defined','sum')).to_csv(OUT/'selected_model_inventory.csv')
 for path,h in context_hash.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
 (OUT/'benchmark_context_sha256.json').write_text(json.dumps(context_hash,indent=2))
 md='# Pedro 置信度与局部接触几何\n\n'
 md+='仅使用已保存的界面汇总和既有置信度/QC 元数据；未读取 CIF/ZIP、未重跑任何基准或同种对照。\n\n'
 md+='## 连接与定义\n'+(OUT/'join_audit.json').read_text()+'\n\n主置信度 best_model_pair_iptm 是所选模型的链对 ipTM；mean_pair_iptm 是原流程跨模型平均链对 ipTM；job_ranking_score 是整个预测复合物的分数，同一作业多个界面共享，不能视为独立界面置信度。现有输入表未包含 pTM，未从结构或归档提取新字段。\n\n'
 md+='主结局为 3 < d_min ≤4 Å 的残基对比例，直接复用 count_3_4 / n_contacts。零界面不参与条件几何分析；单独用于接触形成概率的描述。分位组按非零界面定义，零接触分析重新使用全部界面分位组；两者边界不可混用。使用 qcut，保持相同分数不跨组，重复边界合并，实际边界和样本量已保存。高置信度表示相对分位数，不是外部校准的可信 PPI 阈值。\n\n'
 md+='几何分层的 median_d_min 是“界面内距离中位数的跨界面中位数”，不是池化距离中位数。池化接触比例从已有精确分箱计数加总，无需再次读取接触明细。Pearson/Spearman 均为描述性关联，不报告 p 值。界面仍嵌套于复合物，重复蛋白对和作业不独立，不推断因果。\n\n'
 for title,frame in [('置信度字段',pd.read_csv(OUT/'confidence_field_inventory.csv')),('连续关联',corr),('置信度十分位',dec),('三分位',terc),('大小关联',sizecorr),('同尺寸分层',pd.DataFrame(size)),('QC 对比',qc),('零与非零置信度',pd.DataFrame(zd)),('零概率十分位',zero),('实验上下文',comp)]:md+='## '+title+'\n\n```\n'+frame.to_string(index=False)+'\n```\n\n'
 (OUT/'PEDRO_CONFIDENCE_ANALYSIS.md').write_text(md)
 (OUT/'completion.json').write_text(json.dumps(dict(status='complete',all_interfaces=len(s),nonzero_interfaces=len(nz),primary=PRIMARY,geometry_deciles=len(dec),zero_deciles=len(zero),figures=9),indent=2))
 print(corr[corr.version=='raw'].to_string(index=False));print(dec.to_string(index=False));print(comp[['group','pooled_close_fraction','interface_close_median']].to_string(index=False))
if __name__=='__main__':main()
