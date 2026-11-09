"""Transparent upper-quantile sensitivity; reads only the new audited interface table."""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1]/'results/local_interface_geometry/pedro_confidence_analysis'
s=pd.read_csv(P/'pedro_confidence_interface_table.csv',usecols=['n_contacts','best_model_pair_iptm','fraction_3_4','count_3_4','existing_severe_overlap','job_key','protein_pair']);s=s[s.n_contacts>0];rows=[]
for q in [.9,.95,.99,.995,.999]:
 t=s.best_model_pair_iptm.quantile(q);g=s[s.best_model_pair_iptm>=t]
 rows.append(dict(quantile=q,minimum_score=t,n=len(g),actual_fraction=len(g)/len(s),median_score=g.best_model_pair_iptm.median(),pooled_close=g.count_3_4.sum()/g.n_contacts.sum(),median_close=g.fraction_3_4.median(),q25=g.fraction_3_4.quantile(.25),q75=g.fraction_3_4.quantile(.75),jobs=g.job_key.nunique(),protein_pairs=g.protein_pair.nunique(),median_contacts=g.n_contacts.median(),severe=int(g.existing_severe_overlap.sum())))
u=pd.DataFrame(rows);u.to_csv(P/'upper_tail_confidence_sensitivity.csv',index=False)
c=pd.read_csv(P/'high_confidence_benchmark_comparison.csv');c=c[~c.group.str.startswith('Pedro upper quantile')]
for q in [.99,.999]:
 r=u[u['quantile']==q].iloc[0];c=pd.concat([c,pd.DataFrame([dict(group=f'Pedro upper quantile {q:g} (n={int(r.n)})',n_interfaces=r.n,pooled_close_fraction=r.pooled_close,interface_close_median=r.median_close,interface_close_q25=r.q25,interface_close_q75=r.q75)])],ignore_index=True)
c.to_csv(P/'high_confidence_benchmark_comparison.csv',index=False)
fig,ax=plt.subplots(figsize=(11,7));y=np.arange(len(c));ax.scatter(100*c.pooled_close_fraction,y,label='Pooled contacts',marker='s');ax.errorbar(100*c.interface_close_median,y+.15,xerr=[100*(c.interface_close_median-c.interface_close_q25),100*(c.interface_close_q75-c.interface_close_median)],fmt='o',label='Interface median / IQR');ax.set_yticks(y);ax.set_yticklabels(c.group);ax.invert_yaxis();ax.set(xlabel='Contacts in (3,4] Å (%)');ax.legend();fig.tight_layout();fig.savefig(P/'figures/05_high_confidence_vs_experimental_benchmarks.png',dpi=180);plt.close(fig)
f=P/'PEDRO_CONFIDENCE_ANALYSIS.md'
text='''
## 模式判断：更支持 B，但需区分整体趋势和极高分尾部

整体是较弱的正关联（ρ=0.160），并非 C 所指的完全无关联；最高十分位池化比例仅 8.90%，仍明显低于实验 12–14%，不支持将总体缺口主要归因为低置信度并选择 A。更符合 **B：置信度只解释部分现象，存在明显剩余差异**。这里的“解释”仅是描述，不是回归解释方差或因果解释。

接触数与近接触比例的 ρ=0.556，总界面残基数为 0.532，界面覆盖比例为 0.534，总链长约为 −0.006。按接触数四分位分层后，置信度关系减弱至 ρ≈−0.028 至 +0.014，说明总体关系与尺寸组成密切相关；不能断言尺寸造成缺口或已经消除所有混杂。较大尺寸层 Pearson 仍为正，也提示可能存在非线性高分尾部，不能说任何尺寸下都没有置信度效应。

### 对“最高置信度”的补充检查

最高十分位中位 ipTM 只有 0.19，故另检查上 5%、1%、0.5%、0.1% 的嵌套分位阈值。阈值保留同分数全部界面，实际样本比例略有偏离；这只是透明的分位敏感性，不是新增生物学置信度标准。上尾 >=90% 阈值包含边界同分数，和原 qcut 最后一组的右闭分箱不是完全同一个集合。

```
'''+u.to_string(index=False)+'''
```

最顶端约 0.1% 只有 55 个界面（52 个归档作业、39 个蛋白对），池化 11.67%、界面中位数 11.57%，确实更接近 Mgen 敏感性 12.11–12.45%，但仍低于一般实验 13.62–14.36%。这些界面更大（接触数中位数 445），且样本很小、选择特殊，所以不能用它们证明整个 Pedro 缺口只是低置信度预测驱动。结论应保留这一尾部趋近现象，而不是声称所有高分预测均远低于实验。

### 零接触与极短接触需分别理解

全部界面的分数存在大量并列值，qcut 原请求十组，重复边界合并后实际为 **9 组**；非零几何分析实际为 **10 组**，各组样本数并不相等。零接触率从最低组 88.16% 降到最高组 24.49%；零/非零 ipTM 中位数分别为 0.07/0.10。这是接触形成关联，不是条件几何证据。

<2 Å 比例对置信度的 Spearman 为 +0.147，但 Pearson 约 +0.002，许多界面该比例为零，分位池化短接触率没有明显单调下降；不能声称高 ipTM 自动消除了短距离疑点。QC 通过子集的相关性相近。既有严重重叠排除后近接触池化比例仅从 8.283% 变为 8.270%，说明此类 QC 失败不是主要缺口来源。
'''
f.write_text(f.read_text()+text)
