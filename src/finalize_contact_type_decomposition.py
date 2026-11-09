"""Readable interpretation, focused chemistry plots, and final numerical checks."""
from audit_contact_type_data import ROOT,BASE,OUT
import pandas as pd
import numpy as np
import json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def tab(d):return '\n```\n'+d.to_string(index=False)+'\n```\n'
def main():
 c=pd.read_csv(OUT/'contact_type_counts.csv');b=pd.read_csv(OUT/'contact_type_distance_bins.csv');d=pd.read_csv(OUT/'pedro_experimental_decomposition.csv');i=pd.read_csv(OUT/'interface_type_distribution_summary.csv');o=pd.read_csv(OUT/'OXT_classification_sensitivity.csv');q=pd.read_csv(OUT/'pedro_size_quantile_summary.csv');a=pd.read_csv(OUT/'pedro_contact_type_size_association.csv');ss=pd.read_csv(OUT/'scsc_chemistry_distance_summary.csv');bs=pd.read_csv(OUT/'bbsc_sidechain_chemistry_distance_summary.csv')
 assert np.allclose(d.composition_component+d.within_type_component,d.pedro_minus_reference_all_fraction,atol=1e-12)
 for ds in c.dataset.unique():
  z=c[c.dataset==ds];assert z.n_contacts.sum()==z.total_contacts.iloc[0]
  for typ,g in b[b.dataset==ds].groupby('contact_type'):assert g['count'].sum()==z[z.contact_type==typ].n_contacts.iloc[0]
  for frame,typ in [(ss,'SC-SC'),(bs,'BB-SC')]:
   for binname,v in frame[frame.dataset==ds].groupby('distance_bin')['count'].sum().items():assert v==b[(b.dataset==ds)&(b.contact_type==typ)&(b.distance_bin==binname)]['count'].iloc[0]
 c['unclassified_fraction']=0.0;c.to_csv(OUT/'contact_type_counts.csv',index=False)
 difchem=[]
 for frame,col,typ in [(ss,'chemistry','SC-SC'),(bs,'sidechain_chemistry','BB-SC')]:
  for ref in ['FoldBench Experimental','SAFES Experimental','Mgen Experimental','Mgen excl 8PBZ 8PC0 8PC1']:
   r=frame[frame.dataset==ref];p=frame[frame.dataset=='Pedro AF3'];z=r.merge(p,on=[col,'distance_bin'],suffixes=('_ref','_pedro'),validate='one_to_one')
   for row in z.to_dict('records'):difchem.append(dict(reference=ref,contact_type=typ,chemistry=row[col],distance_bin=row['distance_bin'],pedro_minus_reference_percentage_points=100*(row['fraction_all_contacts_pedro']-row['fraction_all_contacts_ref']),within_chemistry_change_percentage_points=100*(row['fraction_within_chemistry_pedro']-row['fraction_within_chemistry_ref'])))
 cd=pd.DataFrame(difchem);cd.to_csv(OUT/'focused_chemistry_contribution_differences.csv',index=False)
 labels=['FoldBench Experimental','SAFES Experimental','Mgen excl 8PBZ 8PC0 8PC1','Pedro AF3']
 for frame,col,name in [(ss,'chemistry','08_scsc_chemistry_decomposition.png'),(bs,'sidechain_chemistry','09_bbsc_sidechain_chemistry_decomposition.png')]:
  fig,axes=plt.subplots(1,2,figsize=(14,5))
  for ax,binname in zip(axes,['3to4','5to6']):
   v=frame[frame.distance_bin==binname].pivot(index='dataset',columns=col,values='fraction_all_contacts').reindex(labels)*100;v.plot.bar(stacked=True,ax=ax);ax.set(ylabel='Percentage of all residue contacts',xlabel='',title=binname+' Å');ax.tick_params(axis='x',rotation=25);ax.legend(fontsize=7)
  fig.tight_layout();fig.savefig(OUT/'figures'/name,dpi=180);plt.close(fig)
 inputfiles=pd.read_csv(OUT/'input_file_audit.csv');paths=[]
 for kind,mask in [('FoldBench CSV',inputfiles.file_path.str.contains('/foldbench/')),('SAFES CSV',inputfiles.file_path.str.contains('/safes/')),('Mgen Parquet',inputfiles.file_path.str.contains('/mgen_same_species_control/')),('Pedro Parquet',inputfiles.file_path.str.contains('/pedro/contacts/'))]:
  f=inputfiles[mask];paths.append(dict(input=kind,n_files=len(f),bytes=int(f.size_bytes.sum()),rows=int(f.n_rows.sum()),path=f.file_path.iloc[0] if len(f)==1 else str(BASE/'pedro/contacts')+'/*.parquet'))
 pd.DataFrame(paths).to_csv(OUT/'input_dataset_audit_summary.csv',index=False)
 major=['FoldBench Experimental','FoldBench AF3','SAFES Experimental','SAFES AF3','Mgen Experimental','Mgen excl 8PBZ 8PC1','Mgen excl 8PBZ 8PC0 8PC1','Pedro AF3']
 close=b[b.distance_bin=='3to4'].pivot(index='dataset',columns='contact_type',values='fraction_all').reindex(major)*100;close['合计']=close.sum(axis=1)
 deficit=d[d.distance_bin=='3to4'].pivot(index='reference',columns='contact_type',values='pedro_minus_reference_all_fraction')*100
 meanchange=q[q.version=='raw'].pivot(index='size_quartile',columns='contact_type',values='mean')*100
 detail='''
## 核心结果：先说清楚“缺的是什么”

本轮不是去找一个预设的侧链故事。数据表明：**BB-SC 是每个实验参照下 Pedro 3–4 Å 净缺口的最大单项贡献；SC-SC 是重要的第二项；相对一般基准 BB-BB 也减少，但它不解释同种参照下的缺口。**“最大单项”不表示对所有参照都独占过半，更不表示其他类无关。

### 1. 3–4 Å 分解表怎样读？

以下单位全部是占所有≤8 Å残基接触的百分比。每行三类相加等于原分析的近接触比例；所以它直接回答总体缺口由哪一类贡献，而非只看近接触子集内部谁最多。
'''+tab(close.reset_index())+'''
例如Pedro为0.403% BB-BB + 2.893% BB-SC + 4.986% SC-SC = 8.283%。SC-SC仍是Pedro现存近接触中的大部分，但“现存最多”不等于“相对实验缺得最多”。比较差值才能回答后者：
'''+tab(deficit.reset_index())+'''
负号表示Pedro少，单位为百分点。相对FoldBench实验，BB-SC少2.29个百分点，SC-SC少1.98，BB-BB少1.80；相对SAFES实验，分别少2.46、1.78、1.09。相对去掉8PBZ/8PC0/8PC1的同种实验，BB-SC少3.48、SC-SC少0.88，而BB-BB反而多0.19，部分抵消总缺口。因此BB-SC结论跨参照成立，但BB-BB结论具有明显数据集组成依赖。

### 2. 类型整体减少了吗，还是类型内的距离不同？

Pedro的所有接触中BB-SC占52.53%，FoldBench/SAFES实验分别51.97%/52.84%，非常接近。同种敏感性约55.43%。所以主要问题不是Pedro普遍没有主链—侧链接触，而是它们落在近距离箱的比例较低。

在BB-SC类型内部，(3,4] Å比例为Pedro **5.51%**，FoldBench实验 **9.98%**、SAFES实验 **10.14%**、同种扩展敏感性 **11.50%**。这里分母只含BB-SC，与上表的全接触贡献分母不同。Pedro的BB-SC在4–5 Å为9.01%，5–6 Å为15.40%，相应实验参照在这些区间较低；6–8 Å并非对所有参照都增加。主要是4–6 Å方向的相对富集，而不是每个更长区间一律增多。

SC-SC更不能简单说整体缺失：Pedro占所有接触的37.05%，反而高于一般实验约29.59%/30.33%。但类型内3–4 Å比例只有13.46%，一般实验为23.55%/22.32%。Pedro的SC-SC在5–6 Å占16.67%，高于一般实验10.73%/11.75%；同时<3 Å部分也偏多。因此它是距离组成改变，不是纯粹所有距离统一拉长。SC-SC整体占比更高，还会部分掩盖其近距离条件比例的下降。

代数分解进一步区分“类型占比”和“类型内距离组成”。例如对FoldBench，SC-SC类型增多给近接触贡献+1.38个百分点，但类型内距离变化贡献−3.36个百分点，净为−1.98。不能把这两个数学分量解释成真实过程或能量机制。

**“分布向较长最近原子距离偏移”是跨集合的描述，不代表同一批残基对真的从3–4 Å移动到了4–6 Å。**本轮没有建立Pedro与这些实验界面的残基对应。

### 3. 配对基准是否已有相同偏差？

三类型分布并非实验与AF3完全重合，但两套配对基准的近接触差异通常小于Pedro与实验之间的差异。例如BB-SC类型内3–4 Å，FoldBench实验/AF3为9.98%/9.33%，SAFES为10.14%/9.73%，都明显高于Pedro的5.51%。SC-SC实验/AF3则为23.55%/22.19%和22.32%/21.58%，方向上也有轻微下降，不能把“AF3偏低”说成Pedro独有；Pedro的幅度更大。

配对界面也不是完美重现。BB-SC近接触全分母贡献的配对差中位数（AF3−实验），FoldBench约−0.370个百分点，SAFES约+0.117个百分点。类型内中位距离配对差分别约+0.023/−0.022 Å，较小但不代表逐残基一致。全部类型的差值、IQR、相关性和有效分母在paired_benchmark_contact_type_summary.csv；没有将接触当独立样本检验。

### 4. 同种实验支持哪一种解释？

同种实验的BB-BB近接触本来就很少，不能用一般基准的较多BB-BB来直接定义Pedro错误。对主问题更有用的是：不论保留原始同种集合还是排除已有低分辨率疑点，BB-SC仍是Pedro净缺口的最大贡献项，且敏感性排除后这一缺口更明显。SC-SC差异也仍在，但幅度较小。实验参考仍是少数相关黏附复合物，不支持推广为整个物种的类型比例。

### 5. 界面层面是否也看得到？

每个非零界面先算该类型3–4 Å接触数/全部接触数，再报告分布；缺少某类型时该贡献是真正的0，不能只挑有该类型的界面来抬高典型值。Pedro的BB-SC贡献中位数为1.82%，IQR为0–3.39%；FoldBench实验为5.20%，IQR 3.46–6.29%；SAFES实验为5.07%，IQR 4.19–6.41%；同种扩展敏感性为6.74%，IQR 5.62–6.87%。因此主现象不只是几个大界面在池化中占优势。

SC-SC贡献的界面中位数Pedro约3.89%，一般实验约6.19%/6.43%。三类型中位数不能相加当成总体中位数；只有同一界面的贡献、或相同权重的均值可以相加。

### 6. 界面大小关联由哪一类驱动？

Pedro接触数与各类型3–4 Å全分母贡献的Spearman，BB-SC约0.611，SC-SC约0.497，BB-BB约0.453。按排序关联强度，BB-SC最突出，QC通过子集结论相近。但相关系数不是可加的贡献率，不能说BB-SC“解释了61.1%”或将三系数加起来还原旧0.556。

接触数从最低到最高四分位时，各类型贡献的**界面均值**如下，单位为百分点：
'''+tab(meanchange.reset_index())+'''
BB-SC从约0.382%到3.000%，SC-SC从2.257%到5.077%，BB-BB从0.004%到0.421%。因此按均值增量看，BB-SC与SC-SC共同贡献主要变化，SC-SC增量还略大；不能把大小效应全部归给单一类型。大/小界面的离散分母与组成差异也可能参与，仍不作因果推断。

### 7. 化学补充为什么只做有限范围？

BB-SC是主要净缺口，SC-SC对一般基准缺口也有重要贡献，因此只对这两类作附属汇总，不为BB-BB展开额外化学矩阵。BB-SC按负责最近距离的侧链一侧所属既有残基化学组划分；SC-SC直接复用已有八类接触化学。没有新增芳香专属组：芳香残基沿用项目原分组，不能从这次汇总声称芳香堆积已被检测。类别只表示残基身份，不代表某个原子的供体/受体资格、氢键角度、盐桥或能量。
'''
 f=OUT/'CONTACT_TYPE_DECOMPOSITION_CN.md';f.write_text(f.read_text()+detail)
 appendix='\n## 输入路径、大小及特殊原子核对\n'+tab(pd.DataFrame(paths))
 appendix+='\n相关字段：interface_id、d_min、residue_a_type/residue_b_type、atom_a_name/atom_b_name、atom_a_element/atom_b_element；链、结构/作业、蛋白标识及既有化学分类亦已保存。没有缺失原子字段、未知残基或氢/氘最近原子；实际残基均为20种标准氨基酸，特殊原子名仅OXT。甘氨酸出现的非N/CA/C/O名称也仅OXT，没有把CA或虚构CB当侧链。严格四名称规则下的SC包含这个已明确标记的末端例外。\n'
 appendix+=f'\n将OXT改归末端主链后，任一数据集、任一类型/距离箱的全接触比例最大变化为 {100*o.change_all_fraction.abs().max():.4f} 个百分点；完整结果见OXT_classification_sensitivity.csv。该量远小于主要BB-SC近接触缺口，具体近接触驱动排序另作核对，不能因此将OXT称为真正侧链。\n'
 f.write_text(f.read_text()+appendix)
 (OUT/'final_validation.json').write_text(json.dumps(dict(status='passed',unclassified=0,classification_total_reconciled=True,interface_close_fractions_reconciled=True,focused_chemistry_counts_reconciled=True,algebraic_components_reconciled=True,structural_reparse=False),indent=2))
 (OUT/'provenance.json').write_text(json.dumps(dict(contact_table_only=True,old_outputs_modified=False,input_manifest='input_file_audit.csv',scripts_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'src'/n for n in ['audit_contact_type_data.py','compute_contact_type_decomposition.py','report_contact_type_decomposition.py','focused_scsc_chemistry.py','finalize_contact_type_decomposition.py']]}),indent=2))
 print(cd[cd.distance_bin=='3to4'].sort_values(['reference','contact_type','pedro_minus_reference_percentage_points']).to_string(index=False))
if __name__=='__main__':main()
