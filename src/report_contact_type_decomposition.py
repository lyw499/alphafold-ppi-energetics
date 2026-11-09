"""Plots and readable reporting from new decomposition tables only."""
from audit_contact_type_data import ROOT,BASE,OUT
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
TYPES=['BB-BB','BB-SC','SC-SC'];COLORS=['#4c78a8','#f2a541','#54a24b'];BINS=['lt2','2to3','3to4','4to5','5to6','6to8'];LABELS=['<2','[2,3]','(3,4]','(4,5]','(5,6]','(6,8]']
MAJOR=['FoldBench Experimental','FoldBench AF3','SAFES Experimental','SAFES AF3','Mgen Experimental','Pedro AF3'];ALL=MAJOR[:-1]+['Mgen excl 8PBZ 8PC1','Mgen excl 8PBZ 8PC0 8PC1']+['Pedro AF3']
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':9})
def save(fig,name):fig.tight_layout();fig.savefig(OUT/'figures'/name,dpi=180);plt.close(fig)
def table(d):return '\n```\n'+d.to_string(index=False)+'\n```\n'
def main():
 assert json.loads((OUT/'classification_validation.json').read_text())['status']=='passed'
 c=pd.read_csv(OUT/'contact_type_counts.csv');b=pd.read_csv(OUT/'contact_type_distance_bins.csv');h=pd.read_csv(OUT/'contact_type_histogram.csv');i=pd.read_csv(OUT/'interface_type_distribution_summary.csv');dif=pd.read_csv(OUT/'pedro_experimental_decomposition.csv');a=pd.read_csv(OUT/'pedro_contact_type_size_association.csv');q=pd.read_csv(OUT/'pedro_size_quantile_summary.csv');paired=pd.read_csv(OUT/'paired_benchmark_contact_type_summary.csv')
 close=b[b.distance_bin=='3to4'];diffclose=dif[dif.distance_bin=='3to4'];driver=diffclose.groupby('contact_type').pedro_minus_reference_all_fraction.mean().idxmin()
 (OUT/'primary_driver.json').write_text(json.dumps(dict(contact_type=driver,selection='Most negative mean Pedro minus experimental contribution across the five explicitly reported reference views; inspect each reference separately, not a significance test'),indent=2))
 fig,ax=plt.subplots(figsize=(10,5));p=c.pivot(index='dataset',columns='contact_type',values='fraction').reindex(ALL)[TYPES];(p*100).plot.bar(stacked=True,ax=ax,color=COLORS);ax.set(ylabel='All residue contacts (%)',xlabel='');ax.tick_params(axis='x',rotation=30);save(fig,'01_contact_type_composition.png')
 fig,axes=plt.subplots(3,3,figsize=(15,11),sharex=True,sharey=True)
 panels=[['FoldBench Experimental','FoldBench AF3'],['SAFES Experimental','SAFES AF3'],['Mgen Experimental','Mgen excl 8PBZ 8PC0 8PC1','Pedro AF3']]
 for row,typ in enumerate(TYPES):
  for ax,names in zip(axes[row],panels):
   for name in names:
    x=h[(h.dataset==name)&(h.contact_type==typ)];ax.stairs(x.within_type_density,np.r_[x.left,x.right.iloc[-1]],label=name)
   ax.set(title=typ,xlabel='Minimum heavy-atom distance (Å)',ylabel='Within-type density (Å⁻¹)',xlim=(0,8));ax.legend(fontsize=7)
 fig.suptitle('Closest-atom types; pooled residue contacts are descriptive, not independent samples');save(fig,'02_contact_type_distance_distributions.png')
 fig,axes=plt.subplots(1,3,figsize=(17,6),sharey=True)
 for ax,typ in zip(axes,TYPES):
  p=b[b.contact_type==typ].pivot(index='dataset',columns='distance_bin',values='fraction_within').reindex(ALL)[BINS];p.columns=LABELS;(p*100).plot.bar(stacked=True,ax=ax,colormap='viridis');ax.set(title=typ,ylabel='Within-type contacts (%)',xlabel='');ax.tick_params(axis='x',rotation=50);ax.legend(fontsize=7)
 save(fig,'03_distance_bin_by_contact_type.png')
 fig,ax=plt.subplots(figsize=(12,6));p=close.pivot(index='dataset',columns='contact_type',values='fraction_all').reindex(ALL)[TYPES];(p*100).plot.bar(stacked=True,ax=ax,color=COLORS);ax.set(ylabel='All contacts in (3,4] Å (%)',xlabel='',title='Contribution to close-contact fraction = BB-BB + BB-SC + SC-SC');ax.tick_params(axis='x',rotation=30);save(fig,'04_close_contact_3to4_decomposition.png')
 fig,axes=plt.subplots(1,3,figsize=(15,4),sharey=True)
 for ax,typ in zip(axes,TYPES):
  for ref in ['FoldBench Experimental','SAFES Experimental','Mgen excl 8PBZ 8PC0 8PC1']:
   x=dif[(dif.reference==ref)&(dif.contact_type==typ)].set_index('distance_bin').reindex(BINS);ax.plot(LABELS,x.pedro_minus_reference_all_fraction*100,'o-',label=ref)
  ax.axhline(0,color='gray',ls='--');ax.set(title=typ,xlabel='Distance (Å)',ylabel='Pedro − reference (percentage points of ALL contacts)');ax.legend(fontsize=7)
 save(fig,'05_pedro_vs_experimental_contact_type_shift.png')
 fig,axes=plt.subplots(1,3,figsize=(15,6),sharey=True)
 for ax,typ in zip(axes,TYPES):
  x=i[(i.contact_type==typ)&(i.metric=='fraction_all_3to4')].set_index('dataset').reindex(ALL);pos=np.arange(len(x));ax.errorbar(pos,100*x['median'],yerr=[100*(x['median']-x.q25),100*(x.q75-x['median'])],fmt='o');ax.set_xticks(pos);ax.set_xticklabels(ALL,rotation=55,ha='right');ax.set(title=typ,ylabel='Interface contribution to (3,4] Å (%)')
 save(fig,'06_interface_level_contact_type_comparison.png')
 fig,axes=plt.subplots(1,2,figsize=(11,4),sharey=True)
 for ax,version in zip(axes,['raw','existing_QC']):
  for typ,color in zip(TYPES,COLORS):
   z=q[(q.version==version)&(q.contact_type==typ)];ax.plot(z.size_quartile,100*z['mean'],'o-',label=typ,color=color)
  ax.set(title=version,xlabel='Pedro contact-count quartile',ylabel='Mean interface contribution to (3,4] Å (%)');ax.legend()
 save(fig,'07_interface_size_vs_contact_type_close_fraction.png')
 # Difference between pooled weight and conditional type spectrum is explicit.
 md='# 最近原子接触类型分解：Pedro近接触缺口\n\n'
 md+='## 为什么做这一步\n\n此前知道Pedro的(3,4] Å比例偏低，但不知道是主链接触还是侧链接触贡献了差异。本轮将所有≤8 Å残基接触分类，而非只看仍留在3–4 Å的那部分。比较三层：类型在全部接触中的份额；类型内部的距离组成；每个界面的类型贡献。三者不同。\n\n'
 md+='## 输入审计：无需重解析\n'+(OUT/'data_audit.json').read_text()+'\n\n逐文件路径、字节数、行数、全部列名见input_file_audit.csv；缺失值逐列核对见missing_field_audit.csv。所有数据来自已有接触CSV/Parquet和界面摘要，未读取任何CIF/ZIP或重新执行基准/参考结构分析。\n'
 md+='## 定义、边界与验证\n\nN、CA、C、O为主链BB，其余已保存蛋白残基重原子名称为SC；BB-SC不分链方向。甘氨酸CA仍为BB。标签仅描述实现d_min的那一对原子；同一残基对仍可能有稍远的其他类型相互作用。已有算法对相同最小距离使用确定性原子顺序选择，本轮不重新解析坐标或枚举平局替代标签。\n\nOXT若存在，按本次明确给定的四原子主链规则记SC，但它是末端原子，不能在解释中硬说为生物学侧链；所有出现的原子名/残基组合及特殊标记列在atom_name_inventory.csv。非标准残基沿同一名称规则分类并保留审计，不凭空指定化学作用。\n\n距离箱采用 <2、[2,3]、(3,4]、(4,5]、(5,6]、(6,8] Å，恰好2 Å进入第二箱，确保原[0,3]区间完整拆分。<2是严格阈值，8包含边界。标准名程序测试覆盖N/CA/C/O、CB/CG/CD/CE/NZ/OD1/OE1/SG及两侧方向，样本表保留原字段供检查。\n\n'
 md+='fraction_all = 某类型某距离箱计数 / 全部接触数；fraction_within = 同一计数 / 该类型全部接触数。例如前者低，可能是类型少，也可能是类型内部近距离少。三种fraction_all在3–4 Å的和必须等于旧总近接触比例。界面零接触时所有比例NA；非零界面缺少某类型时贡献为0，但该类型内部比例和距离中位数为NA。\n'
 md+='\n## 类型总量\n'+table(c)+'\n每组分类总数与已有总接触数一致；unclassified为0。绝对总数受数据集大小影响，跨数据集“减少”指比例，不能把百万与万的原始计数直接比较。\n'
 md+='\n## 主结果：3–4 Å贡献\n'+table(close)+'\nPedro减实验的贡献差（负号为Pedro较少）：\n'+table(diffclose)
 md+='\ncomposition_component与within_type_component把贡献p×q的差作对称代数分解：(pP−pE)(qP+qE)/2 与 (qP−qE)(pP+pE)/2。p为该类型占全部接触的比例，q为类型内落入该箱的比例。这是构成分解，不是因果效应；对每个参考单独解释，不把多个参考当作独立复制。\n'
 md+='\n## 类型内部是否偏向更长距离\n'+table(b)+'\n所有距离箱都保留。较长区间富集可支持“距离组成偏向更长的最近原子距离”，不能说同一批残基接触移动了；Pedro与这些实验集合没有直接残基对应。\n'
 md+='\n## 界面层面与大小\n'+table(i[i.metric.isin(['fraction_all_3to4','median_d_min'])])+table(a)+table(q)
 md+='\n这里相关性的分母是整个界面的接触数，因此三类型贡献相加为总近接触比例；各类型相关系数不能相加分摊0.556。尺寸分位层中的贡献均值可显示增量主要体现在哪一类，仍不能推出界面大小导致几何变化。\n'
 md+='\n## 配对基准检查\n'+table(paired)+'\n配对贡献包含已形成界面中缺少某类型的零贡献；类型中位距离只比较双方确有该类型的界面。FoldBench另外保留有序链映射明确子集。这里是配对界面统计，不声称三类型中的每一个残基对都相互对应。\n'
 md+='\n## 同种敏感性与解释边界\n\nMgen原始、排除8PBZ/8PC1、再排除8PC0各自保留，主原始结果未替换。样本集中于少数相关装配、低分辨率模型和黏附体系，不能推广为物种全局规律。本轮未删除原始接触；所有接触型比较为描述性，不做残基级p值。最近原子侧链几何不等于完整侧链堆积能、氢键网络、盐桥能、范德华能或结合自由能。\n'
 (OUT/'CONTACT_TYPE_DECOMPOSITION_CN.md').write_text(md)
 print('Driver for focused follow-up:',driver);print(diffclose.to_string(index=False))
if __name__=='__main__':main()
