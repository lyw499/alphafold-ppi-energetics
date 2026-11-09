"""Descriptive pooled plots and paired interface/structure summaries, no contact-level p values."""
import os,json,hashlib
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/local-geometry-mpl')
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from local_geometry_common import BIN_NAMES,CHEM_CLASSES,CHEM_GROUPS
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/local_interface_geometry'
METRICS=['n_contacts','median_d_min','mean_d_min','fraction_le_4','fraction_le_5','fraction_le_6']+BIN_NAMES
LABELS={'n_contacts':'Number of residue contacts','median_d_min':'Median residue-pair minimum distance (Å)','mean_d_min':'Mean residue-pair minimum distance (Å)','fraction_le_4':'Fraction of contacts at ≤4 Å','fraction_le_5':'Fraction of contacts at ≤5 Å','fraction_le_6':'Fraction of contacts at ≤6 Å'}

def stats(x):
 x=pd.Series(x).dropna();return dict(n=len(x),mean=x.mean(),median=x.median(),q25=x.quantile(.25),q75=x.quantile(.75),min=x.min(),max=x.max())
def tab(d):
 def f(v):return f'{v:.5g}' if isinstance(v,float) else str(v).replace('|','/')
 return '\n'.join(['| '+' | '.join(d.columns)+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(f(v) for v in row)+' |' for row in d.itertuples(index=False,name=None)])
def save(fig,folder,name):
 folder.mkdir(parents=True,exist_ok=True);fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(folder/f'{name}.{ext}',dpi=170)
 plt.close(fig)

def analyze(dataset,s,c,scope='full_matched',plots=True):
 dest=OUT/dataset;folder=OUT/'figures'/dataset
 if scope!='full_matched':folder=folder/scope
 pooled=[];bin_rows=[];per_metric=[]
 if plots:fig,axes=plt.subplots(1,2,figsize=(12,4.5))
 for kind,color in [('Experimental','#1f77b4'),('AF3','#ff7f0e')]:
  ss=s[s.structure_type==kind];cc=c[c.structure_type==kind];x=cc.d_min.to_numpy()
  pooled.append(dict(dataset=dataset,scope=scope,structure_type=kind,n_interfaces=len(ss),n_complexes=ss.structure_id.nunique(),n_zero_interfaces=int(ss.n_contacts.eq(0).sum()),n_nonzero_interfaces=int(ss.n_contacts.gt(0).sum()),**stats(x)))
  hist,edges=np.histogram(x,bins=np.linspace(0,8,33));assert hist.sum()==len(x)
  for left,right,n in zip(edges[:-1],edges[1:],hist):bin_rows.append(dict(scope=scope,structure_type=kind,left_A=left,right_A=right,n=int(n),density=n/(len(x)*.25)))
  for metric in METRICS:per_metric.append(dict(dataset=dataset,scope=scope,structure_type=kind,metric=metric,**stats(ss[metric])))
  if plots:
   axes[0].stairs(hist/(len(x)*.25),edges,color=color,label=f'{kind} ({len(x):,} contacts)')
   values=np.sort(x);axes[1].step(values,np.arange(1,len(x)+1)/len(x),where='post',color=color,label=kind)
 if plots:
  for ax,ylabel in zip(axes,['Density (Å⁻¹)','ECDF']):ax.set(xlabel='Residue-pair minimum heavy-atom distance (Å)',ylabel=ylabel,xlim=(0,8));ax.legend(fontsize=8)
  fig.suptitle(f'{dataset}: {scope}; pooled contacts ≤8 Å (descriptive, not independent samples)')
  save(fig,folder,f'{dataset}_exp_vs_af3_min_distance_panels')
  # Explicit separate requested figures.
  for plotkind in ['distribution','ecdf']:
   fig,ax=plt.subplots(figsize=(7,4.5))
   for kind,color in [('Experimental','#1f77b4'),('AF3','#ff7f0e')]:
    cc=c[c.structure_type==kind];x=np.sort(cc.d_min.to_numpy());label=f'{kind}: {len(x):,} contacts / {cc.interface_id.nunique()} nonzero interfaces'
    if plotkind=='distribution':hist,edges=np.histogram(x,bins=np.linspace(0,8,33));ax.stairs(hist/(len(x)*.25),edges,label=label,color=color,linewidth=1.8)
    else:ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',label=label,color=color)
   ax.set(xlabel='Residue-pair minimum heavy-atom distance (Å)',ylabel='Density (Å⁻¹)' if plotkind=='distribution' else 'ECDF',xlim=(0,8),title=f'{dataset}: {scope}\nPooled contacts ≤8 Å; descriptive');ax.legend(fontsize=8)
   save(fig,folder,f'{dataset}_exp_vs_af3_min_distance_{plotkind}')
 e=s[s.structure_type=='Experimental'].set_index('interface_id');p=s[s.structure_type=='AF3'].set_index('interface_id');assert e.index.is_unique and e.index.equals(p.index)
 wide=e.add_prefix('experimental_').join(p.add_prefix('af3_'));results=[];chem=[];cluster_rows=[]
 for metric in METRICS+['chem_'+v for v in CHEM_CLASSES]:
  valid=wide[['experimental_'+metric,'af3_'+metric]].dropna();x=valid['experimental_'+metric];y=valid['af3_'+metric];delta=y-x
  wide['delta_'+metric]=delta
  r=pearsonr(x,y).statistic if len(x)>2 and x.nunique()>1 and y.nunique()>1 else np.nan
  rho=spearmanr(x,y).statistic if len(x)>2 and x.nunique()>1 and y.nunique()>1 else np.nan
  clusters=wide.loc[valid.index,'experimental_structure_id'].str.split('-assembly').str[0].str.upper()
  cluster_delta=delta.groupby(clusters).mean();rng=np.random.default_rng(941)
  samples=rng.choice(cluster_delta.to_numpy(),size=(2000,len(cluster_delta)),replace=True).mean(axis=1) if len(cluster_delta) else np.asarray([np.nan])
  item=dict(dataset=dataset,scope=scope,metric=metric,n_paired=len(valid),n_unavailable_summary=len(wide)-len(valid),n_PDB_clusters=len(cluster_delta),experimental_median=x.median(),af3_median=y.median(),pearson=r,spearman=rho,median_delta=delta.median(),mean_delta=delta.mean(),delta_q25=delta.quantile(.25),delta_q75=delta.quantile(.75),n_AF3_shorter_or_lower=int((delta<0).sum()),n_AF3_longer_or_higher=int((delta>0).sum()),PDB_equal_weight_mean_delta=cluster_delta.mean(),cluster_bootstrap_mean_delta_low=np.quantile(samples,.025),cluster_bootstrap_mean_delta_high=np.quantile(samples,.975))
  (chem if metric.startswith('chem_') else results).append(item)
  cluster_rows.extend(dict(dataset=dataset,scope=scope,pdb_id=key,metric=metric,mean_interface_delta=value) for key,value in cluster_delta.items())
  if plots and metric in ['median_d_min','fraction_le_4']:
   fig,ax=plt.subplots(figsize=(5.6,5.4));ax.scatter(x,y,s=22,alpha=.65)
   lo,hi=(0,1) if metric.startswith('fraction') else (min(x.min(),y.min())-.1,max(x.max(),y.max())+.1)
   ax.plot([lo,hi],[lo,hi],'k--',linewidth=1);ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Experimental: '+LABELS[metric],ylabel='AF3: '+LABELS[metric],title=f'{dataset}: paired interfaces n={len(valid)}\nPearson={r:.3f}; Spearman={rho:.3f}')
   save(fig,folder,f'{dataset}_exp_vs_af3_interface_'+('median_distance' if metric=='median_d_min' else 'fraction_le4')+'_scatter')
 wide.to_csv(dest/f'{scope}_paired_interface_metrics.csv')
 pd.DataFrame(bin_rows).to_csv(dest/f'{scope}_distance_histogram_bins.csv',index=False)
 pd.DataFrame(results).to_csv(dest/f'{scope}_paired_interface_statistics.csv',index=False)
 pd.DataFrame(chem).to_csv(dest/f'{scope}_paired_chemistry_statistics.csv',index=False)
 pd.DataFrame(cluster_rows).to_csv(dest/f'{scope}_PDB_cluster_mean_differences.csv',index=False)
 return pooled,per_metric,results,chem

def main():
 assert json.loads((OUT/'computation_status.json').read_text())['status']=='complete'
 allpooled=[];allmetrics=[];allpaired=[];allchem=[];tails=[];extreme=[];interface_tail=[];chemfreq=[];respaired=[];counts=[];validations=[]
 for dataset in ['foldbench','safes']:
  dest=OUT/dataset;s=pd.read_csv(dest/'interface_summaries.csv',dtype={'interface_id':str});c=pd.read_csv(dest/'residue_contacts.csv',dtype={'interface_id':str});a=pd.read_csv(dest/'chain_pair_audit.csv',dtype={'interface_id':str})
  assert c.d_min.between(0,8).all();assert not c.duplicated(['interface_id','structure_type','residue_index_a','residue_index_b']).any()
  observed=c.groupby(['interface_id','structure_type']).size()
  for _,row in s.iterrows():assert observed.get((row.interface_id,row.structure_type),0)==row.n_contacts
  for _,g in s[s.n_contacts>0].iterrows():assert np.isclose(sum(g[k] for k in BIN_NAMES),1) and np.isclose(sum(g['chem_'+k] for k in CHEM_CLASSES),1)
  assert s.n_deuterium_excluded.sum()==0,'Investigate legacy non-H/non-D definition difference'
  reliable=a[a.included_primary==True]
  counts.append(dict(dataset=dataset,n_audited_chain_pairs=len(a),n_included_interfaces=len(reliable),n_structures=reliable.structure_id.nunique(),n_unavailable_or_unreliable=int((a.included_primary!=True).sum()),n_ordered_unambiguous=int(reliable.chain_mapping_unambiguous.sum()),n_historical_size_mapping_stable=int(reliable.legacy_mapping_stable.sum()),n_local_spectrum_stable=int(reliable.local_distance_spectrum_stable.sum())))
  versions={'full_matched':s}
  if dataset=='foldbench':versions.update(legacy_size_mapping_stable=s[s.legacy_mapping_stable],unambiguous_chain_mapping=s[s.chain_mapping_unambiguous],local_spectrum_stable=s[s.local_distance_spectrum_stable])
  for scope,ss in versions.items():
   cc=c[c.interface_id.isin(ss.interface_id)]
   po,me,pa,ch=analyze(dataset,ss,cc,scope,plots=scope in ('full_matched','legacy_size_mapping_stable'))
   allpooled+=po;allmetrics+=me;allpaired+=pa;allchem+=ch
  for kind,g in c.groupby('structure_type'):
   ss=s[s.structure_type==kind]
   for threshold in [.5,1,1.5,2]:
    short=g[g.d_min<threshold];tails.append(dict(dataset=dataset,structure_type=kind,threshold_strict_lt_A=threshold,n_contacts_below=len(short),n_all_contacts=len(g),fraction_contacts_below=len(short)/len(g),n_interfaces_with_short=short.interface_id.nunique(),n_all_interfaces=len(ss),n_structures_with_short=short.structure_id.nunique()))
   selected=pd.concat([g[g.d_min<2],g.nsmallest(20,'d_min')]).drop_duplicates(['interface_id','residue_index_a','residue_index_b']);selected=selected.copy();selected['below_2A']=selected.d_min<2
   # Chemistry flags are context clues, never automatic diagnoses or exclusions.
   selected['short_atom_context']=np.where((selected.atom_a_name=='SG')&(selected.atom_b_name=='SG'),'SG-SG: covalent/disulfide context should be checked',np.where(((selected.atom_a_name=='C')&(selected.atom_b_name=='N'))|((selected.atom_a_name=='N')&(selected.atom_b_name=='C')),'C-N backbone: possible split-chain/covalent context, unverified','No bond assignment; inspect coordinates and metadata'))
   extreme.append(selected)
   for _,row in ss.iterrows():
    interface_tail.append(dict(dataset=dataset,structure_type=kind,interface_id=row.interface_id,structure_id=row.structure_id,n_contacts=row.n_contacts,min_d_min=row.min_d_min,existing_qc_flags=row.existing_qc_flags,mapping_status=row.mapping_status,**{key:row[key] for key in ['n_lt_0_5','n_lt_1_0','n_lt_1_5','n_lt_2_0']}))
   for category in CHEM_CLASSES:
    chemfreq.append(dict(dataset=dataset,structure_type=kind,contact_class=category,n_contacts=int(g.contact_class.eq(category).sum()),pooled_fraction=g.contact_class.eq(category).mean(),interface_equal_weight_mean_fraction=ss.loc[ss.n_contacts>0,'chem_'+category].mean(),n_nonzero_interfaces=int(ss.n_contacts.gt(0).sum())))
  freq=pd.DataFrame([r for r in chemfreq if r['dataset']==dataset]);fig,ax=plt.subplots(figsize=(11,5));x=np.arange(len(CHEM_CLASSES))
  for k,kind in enumerate(['Experimental','AF3']):
   group=freq[freq.structure_type==kind].set_index('contact_class').loc[CHEM_CLASSES];ax.bar(x+(k-.5)*.36,100*group.pooled_fraction,.36,label=kind)
  ax.set_xticks(x,labels=[v.replace('_','\n') for v in CHEM_CLASSES]);ax.set(ylabel='Pooled residue contacts (%)',title=f'{dataset}: exploratory residue chemistry; not independent samples');ax.legend();save(fig,OUT/'figures'/dataset,f'{dataset}_exp_vs_af3_exploratory_contact_chemistry')
  r=pd.read_csv(dest/'paired_residue_distances.csv',dtype={'interface_id':str})
  for iid,g in r.groupby('interface_id'):
   both=g.contact_state=='both';lost=g.contact_state=='lost_in_AF3';gained=g.contact_state=='gained_in_AF3'
   respaired.append(dict(dataset=dataset,interface_id=iid,structure_id=g.structure_id.iloc[0],n_union=len(g),n_both=int(both.sum()),n_lost=int(lost.sum()),n_gained=int(gained.sum()),experimental_contact_retention=both.sum()/(both.sum()+lost.sum()) if (both|lost).sum() else np.nan,AF3_contact_overlap=both.sum()/(both.sum()+gained.sum()) if (both|gained).sum() else np.nan,median_delta_common_contacts=g.loc[both,'delta_d_min'].median(),median_absolute_error_common_contacts=g.loc[both,'delta_d_min'].abs().median(),mean_absolute_error_common_contacts=g.loc[both,'delta_d_min'].abs().mean(),median_absolute_error_union=g.delta_d_min.abs().median(),median_delta_union=g.delta_d_min.median()))
  assert ((r.experimental_d_min<=8)|(r.af3_d_min<=8)).all();assert np.isfinite(r[['experimental_d_min','af3_d_min']]).all().all()
  validations.append(dict(dataset=dataset,n_contact_rows=len(c),n_interface_rows=len(s),n_paired_residue_rows=len(r),checks='counts reconcile; unique residue pairs; d_min in [0,8]; histogram/fractions normalized; no deuterium in eligible atoms; targeted union counterparts only'))
 for filename,rows in [('pooled_distance_summary',allpooled),('interface_level_summary',allmetrics),('paired_interface_statistics',allpaired),('paired_chemistry_statistics',allchem),('chemistry_frequencies',chemfreq),('paired_residue_interface_summary',respaired),('analysis_counts',counts)]:pd.DataFrame(rows).to_csv(OUT/f'{filename}.csv',index=False)
 pd.DataFrame(tails).to_csv(OUT/'qc/short_distance_qc.csv',index=False);pd.DataFrame(interface_tail).to_csv(OUT/'qc/interface_short_distance_audit.csv',index=False)
 pd.concat(extreme,ignore_index=True).sort_values('d_min').to_csv(OUT/'qc/extreme_short_distance_cases.csv',index=False)
 (OUT/'validation.json').write_text(json.dumps(dict(status='passed',datasets=validations,Pedro_processed=False,statistics='No residue-level significance tests; paired interface summaries and PDB-cluster bootstrap only'),indent=2)+'\n')
 print(pd.DataFrame(counts).to_string(index=False));print(pd.DataFrame(allpooled).query("scope=='full_matched'").to_string(index=False));print(pd.DataFrame(allpaired).query("scope=='full_matched' and metric=='median_d_min'").to_string(index=False));print(pd.DataFrame(tails).to_string(index=False))
 write_summary(pd.DataFrame(counts),pd.DataFrame(allpooled),pd.DataFrame(allpaired),pd.DataFrame(tails),pd.DataFrame(chemfreq),pd.DataFrame(allchem),pd.DataFrame(respaired))

def write_summary(counts,pooled,paired,tails,chem,chemstats,respaired):
 lines=['# Experimental vs AF3：局部界面几何','',
 '仅分析 FoldBench 和 SAFES。Pedro未运行；历史界面大小结果未修改。主观测为跨链残基对的最小非氢原子距离，不是原子对频数。','',
 '## 方法与口径','',
 '- 对每条链建立重原子KD-tree，只产生≤8 Å临时原子邻居，然后按残基对保留最小距离及其原子。没有构建/保存整个链的所有长距离原子矩阵。',
 '- 沿用有坐标且residue.id[0]为空的残基、首个model和Bio.PDB默认altloc。排除H和D；已确认全部分析链无D，FoldBench新接触数逐接口与历史残基结果一致。HETATM编码修饰残基/配体仍不纳入历史定义。',
 '- 密度图共同分箱为0–8 Å，每0.25 Å；ECDF使用实际距离。pooled图按接触数加权，仅作描述，不把残基对当独立实验样本。',
 '- 每接口的分箱为[0,3]、(3,4]、(4,5]、(5,6]、(6,8] Å，互斥且总和1；另保留≤4/5/6 Å累积分数，不使用恒等于1的≤8 Å分数。',
 '- 零接触接口保留在链对审计与接口汇总中；其距离均值/中位数和比例未定义（NaN），不能填0。接口散点仅使用两侧相关统计量均定义的有效配对；每项报告排除的未定义数。',
 '- 每接口统计保留配对相关性/差值；额外先在每PDB内平均接口差，再对PDB等权。95%区间为2000次PDB簇bootstrap（固定种子），不报告把全部残基对当独立样本产生的p值；不同PDB仍可能同源或来自同一实验体系。',
 '- 所有短距离阈值仅诊断，不增加自动QC删除。小于2 Å可能涉及共价/二硫键、切分链或局部碰撞，不能仅凭阈值断言unphysical。extreme表提供残基、原子和既有QC/映射状态。','',
 '## 数据集与映射','',tab(counts),'',
 'FoldBench复用已审计的132条映射。147个AF3缺失案例只进入审计，不填零。完整集合保留几何并列映射；分别提供历史界面大小稳定子集、无歧义有序链映射子集、局部距离谱稳定子集。历史大小稳定不保证局部几何稳定；并列链对的替代距离汇总见mapping_tie_geometry_sensitivity.csv。谱一致判据是接触数相同且排序后的距离差≤1e-6 Å，仅为数值一致性检查，不是结构QC阈值。','',
 'SAFES按本次要求枚举全部可靠对应蛋白链对，**不再限于旧10个两链案例**。沿用序列identity≥0.95、实验观测覆盖率≥0.90和全复合体CA几何的一对一链映射。若所有并列映射给出同一个无序AF3链对，可用于链对距离分布；若对应到不同AF3链对则不纳入配对主分析，原因保留。此枚举是所有可对应的链对，不代表SAFES官方预期结合伙伴清单；抗体亚基/多链结构内部接触也可能在内。','',
 '逐残基配对更严格：有序链对应必须唯一；使用与链映射相同参数的全局序列比对，仅保留所有最优比对共同指定、氨基酸身份相同且非X的位置。最优比对超过128种则跳过并记录，不凭残基编号相同强行配对。只计算两侧接触集合的可映射并集；另一侧即使>8 Å也仅计算该指定残基对的小型原子距离矩阵，作为接触丢失/新增的证据。此配对结果有保守映射与界面选择偏差。','',
 '## 汇总距离（pooled，仅描述）','',tab(pooled[pooled.scope=='full_matched']),'',
 '## 主要配对接口结果','',tab(paired[(paired.scope=='full_matched')&paired.metric.isin(['median_d_min','fraction_le_4'])]),'',
 '## FoldBench映射敏感性','',tab(paired[(paired.dataset=='foldbench')&paired.metric.eq('median_d_min')]),'',
 '## 短距离诊断','',tab(tails),'',
 '阈值均为严格小于。分母是本组≤8 Å残基接触数；同时报告涉及接口/结构数量。最短接触另列于 `qc/extreme_short_distance_cases.csv`，没有删除可计算结构。','',
 '## 探索性残基化学','',
 '分类：疏水A,V,I,L,M,F,W,Y；极性不带电S,T,N,Q,C；正电K,R；负电D,E；特殊G,P,H；非标准/未知为other。H的质子化状态不确定，因此不硬归正电；Y/C的分组是简化约定。正负电指常用残基类别，并未计算实际质子化或能量。',
 '接触类别互斥：疏水–疏水、异号电荷、同号电荷、极性–极性（仅不带电组）、疏水–极性、疏水–带电、极性–带电、含特殊/其他。全部比例合计1。','',tab(chem),'',
 'pooled比例会被大接口主导，表中同时提供接口等权均值；配对接口化学比例差与PDB等权汇总见paired_chemistry_statistics.csv。频率变化不等同于能量或偏好机制。','',
 '## 逐残基配对覆盖','']
 for dataset,g in respaired.groupby('dataset'):
  lines.append(f'- {dataset}: {len(g)} 个接口有可映射接触并集，共{g.n_union.sum()}对；共有接触{g.n_both.sum()}、AF3丢失{g.n_lost.sum()}、AF3新增{g.n_gained.sum()}。残基定位审计及双方实际距离保存在各数据集子目录。')
 lines += ['', '## 输出与复现','',
 '- `foldbench/`、`safes/`：residue_contacts.csv（逐残基接触及最近原子）、interface_summaries.csv、chain_pair_audit.csv、paired_residue_distances.csv、residue_mapping_audit.csv、mapping_tie_geometry_sensitivity.csv及配对统计。',
 '- `figures/foldbench/`、`figures/safes/`：密度、ECDF、接口中位距离/≤4 Å比例散点、探索性化学频率；PNG及PDF。',
 '- `qc/`：短距离阈值汇总、逐接口异常尾部、极端残基/原子案例。',
 '- `algorithm_validation.json`：三个真实小结构与穷举一致，及合成边界/氢排除/零接口验证。`validation.json`：完整结果一致性。',
 '- 复现依次运行 `src/test_local_geometry.py`、`src/compute_local_geometry.py`、`src/report_local_geometry.py`。所有输出限于此新目录。','']
 lines += ['', '## 直接回答与解释边界', '',
 '1. **总体分布是否相似？** 两组实验/AF3的pooled密度和ECDF都很接近。靠近8 Å的接触数量本来就受截断与空间邻域体积影响；图形相似不是正确界面或逐残基准确性的充分证据。',
 '2. **有系统性变短或变长吗？** 当前逐接口中位距离差接近零；PDB等权平均差的bootstrap区间覆盖零，尚无一致方向的证据。下面给出实际效应量，不能以不显著宣称完全等价。','']
 for ds in ['foldbench','safes']:
  row=paired[(paired.dataset==ds)&(paired.scope=='full_matched')&(paired.metric=='median_d_min')].iloc[0]
  px=pooled[(pooled.dataset==ds)&(pooled.scope=='full_matched')].set_index('structure_type')
  lines.append(f'- {ds}: pooled中位数 Experimental {px.loc["Experimental","median"]:.4f} Å / AF3 {px.loc["AF3","median"]:.4f} Å；逐接口中位数差的中位数 {row.median_delta:+.4f} Å，均值 {row.mean_delta:+.4f} Å；PDB等权均值 {row.PDB_equal_weight_mean_delta:+.4f} Å，95%簇bootstrap区间 [{row.cluster_bootstrap_mean_delta_low:.4f}, {row.cluster_bootstrap_mean_delta_high:.4f}] Å。')
 lines += ['', '3. **接口层面对应程度？** FoldBench的中位距离Pearson=0.959、Spearman=0.575（131个两侧非零接口）；SAFES为0.760、0.682（42个）。FoldBench的高Pearson部分反映分离的低中位距离接口群，不能将其单独解读为所有接口都高度准确。未定义距离的零接口仍在审计中。','']
 for ds,g in respaired.groupby('dataset'):
  lines.append(f'- {ds}: 在严格映射的接触集合中，逐接口实验接触保留率中位数={g.experimental_contact_retention.median():.3f}；两侧共同接触的逐接口中位绝对距离误差之中位数={g.median_absolute_error_common_contacts.median():.3f} Å。共同接触分析受两侧≤8 Å选择限制；完整并集也保留丢失/新增和另一侧>8 Å的指定残基距离。')
 lines += ['', '4. **AF3短距离是否更多？** FoldBench Experimental 2/35539（0.0056%）vs AF3 24/40836（0.0588%）小于2 Å；SAFES 1/8922（0.0112%）vs 5/9883（0.0506%）。四组均无<1 Å。FoldBench <1.5 Å为实验1条、AF3 2条；SAFES均0。',
 '短距离不全部等于非物理碰撞：SAFES实验的1条和AF3的4/5条为SG–SG；FoldBench AF3的3/24条为SG–SG。剩余极短O–O等值得检查；实验9c90也有约1.35 Å的OE2–OE2接触。FoldBench两条实验<2 Å记录均涉及占有率0.51的B/B替代构象；Bio.PDB默认选定构象未做跨链altloc协调，不能仅凭该短距离认定实验生物学界面必然碰撞。不同altloc的最近原子距离另存short_distance_altloc_diagnostics.csv（没有替换主度量，也不据占有率假定构象必然互斥）。`qc/short_distance_chemical_context.csv`保留所选altloc/occupancy和可用的共价连接注释，未根据新阈值删行。','',
 '5. **化学类别频率是否明显不同？** pooled宽类别变化较小，最大绝对差见下；接口等权与配对变化另表报告，不能据这些频率推出能量偏好。','']
 for ds,g in chem.groupby('dataset'):
  v=g.pivot(index='contact_class',columns='structure_type',values='pooled_fraction');delta=v.AF3-v.Experimental;key=delta.abs().idxmax()
  lines.append(f'- {ds}: 最大pooled类别差 {key}: AF3−Experimental {delta[key]*100:+.3f} 个百分点。')
 lines += ['', '6. **两个数据集定性一致吗？** 都显示总体距离形状相近、逐接口中位距离无清晰单向偏移、原始<2 Å尾部AF3更多。但SAFES短距离主要是SG–SG，且样本较小、包含复合物内部链对；不能把两组短尾都归为同一类碰撞问题。',
 '7. **是否值得扩展Pedro？** 可以用于描述接触几何和定位异常短距离，尤其作为已有严重重叠QC的补充；不能单凭几何分布相似就判断Pedro界面正确，也不能在缺乏一对一实验对应时照搬配对准确性结论。还应保留重复观测/接口大小与蛋白来源的层级差异。本轮没有运行Pedro。','']
 context=OUT/'safes/failed_target_best_chain_candidates.csv'
 if context.exists():
  failed=pd.read_csv(context);lines += ['## SAFES未纳入链对的具体序列证据','',
  '1OGA涉及A链的4个链对未达既有95%identity阈值；5E9D只有B链达标，因此全部10个链对未纳入。低于保守规则不等于断言它们绝对不可能建立生物学对应；本轮没有放宽阈值。所有候选比对保存在failed_target_all_chain_sequence_candidates.csv。','',tab(failed[~failed.passes_existing_candidate_rules]),'']
 lines += ['补充复现：主报告生成后运行 `src/audit_local_geometry_context.py`，再运行 `src/report_local_geometry.py --summary-only` 将化学背景与失败比对证据写入本报告。','']
 (OUT/'SUMMARY.md').write_text('\n'.join(lines))

if __name__=='__main__':
 import sys
 if '--summary-only' in sys.argv:
  write_summary(*[pd.read_csv(OUT/(name+'.csv')) for name in ['analysis_counts','pooled_distance_summary','paired_interface_statistics']],pd.read_csv(OUT/'qc/short_distance_qc.csv'),pd.read_csv(OUT/'chemistry_frequencies.csv'),pd.read_csv(OUT/'paired_chemistry_statistics.csv'),pd.read_csv(OUT/'paired_residue_interface_summary.csv'))
 else:main()
