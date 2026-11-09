"""Add interpretation and explicit histogram diagnostics from derived tables only."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1]/'results/contact_correspondence'
lines=['\n## 七个核心问题的明确回答\n']
lines+=['1. **实验接触恢复率**：8 Å 下，逐界面中位数 FoldBench 为 87.7%，SAFES 为 93.0%；池化为 64.3% 和 87.5%。分布异质性很大，不能仅引用池化结果代表典型界面。',
'2. **AF3 接触有多少与实验共享**：逐界面精确率中位数为 85.0% 和 90.5%；池化为 61.2% 和 85.1%。逐界面 Jaccard 中位数为 77.3% 和 85.4%。',
'3. **是否主要是刚跨过 8 Å**：FoldBench AF3-only 的实验距离落在 (8,9] Å 的池化比例仅 12.9%，SAFES 为 44.6%；放宽到 (8,10] Å 为 17.1% 和 56.0%。逐界面等权的 (8,9] 比例均值则是 45.3% 和 68.9%，说明 SAFES 的典型界面更符合阈值附近微调。',
'4. **是否有远距离残基形成新接触**：有。AF3-only 的实验距离 >15 Å 比例为 74.1% 和 22.8%，实验距离池化中位数为 29.42 Å 和 9.47 Å；逐界面等权 >15 Å 比例均值为 28.8% 和 8.4%。这些是已映射残基接触网络的明显改变，不是仅跨过 8 Å。但本分析未区分结构构象改变、相对链位置改变、装配差异或预测错误等原因。',
'5. **为何接触更多但距离谱相似**：距离谱忽略残基身份。可靠对应子集中，FoldBench 新增 11,680、丢失 10,226，净增 1,454；SAFES 新增 1,369、丢失 1,113，净增 256。接触网络替换可保持相近距离谱。共享接触逐界面距离变化中位数的中位数仅为 −0.029 Å 和 −0.040 Å，表现为轻微缩短，不能把全部新增接触归因为统一压缩。',
'6. **两组是否一致**：方向一致（净增且共享接触仅轻微缩短），机制比例不一致：FoldBench 的远距离网络变化明显更多，SAFES 典型界面更常表现为阈值附近变化。化学频率变化也不能视为统一机制。',
'7. **截断是否稳健**：5/6/8 Å 下，FoldBench 可对应接触数均净增（480/733/1,454，约 +6.3%/+6.4%/+5.1%）；SAFES 均净增（42/121/256，约 +1.8%/+3.5%/+2.9%）。但 SAFES 在 5 Å 的逐界面净增中位数为 0，不能称每个界面都增多。重叠度的两组差异也在三个阈值下保持。此处只做接触汇总敏感性，不声称所有机制解释都已在每个阈值重新验证。']
lines+=['\n### 映射覆盖率对原始“接触过量”的限制\nFoldBench 109 个已选界面的全接触净增为 4,005；其中可对应子集净增 1,454，未对应部分净差为 2,551（AF3 2,973 个、实验 422 个接触无法纳入对应比较）。SAFES 全接触净增 961，其中可对应子集净增 256，未对应部分净差 705（AF3 714 个、实验 9 个）。因此本分析解释的是可靠对应的接触子集，不能把所有原始接触过量归因于阈值或网络变化。未对应不等于假接触；可能受序列覆盖、缺失残基或比对歧义影响，本轮没有重做映射。']
c=pd.read_csv(P/'chemistry_by_contact_state.csv');con=[]
for ds in ['foldbench','safes']:
 d=pd.read_csv(P/ds/'paired_contact_classification.csv');g=d[d.classification.eq('AF3_only')].copy();g['far']=g.experimental_d_min>15
 counts=g.groupby(['structure_id','interface_id']).agg(AF3_only=('far','size'),far_contacts=('far','sum')).reset_index().sort_values('far_contacts',ascending=False)
 counts.to_csv(P/ds/'far_AF3_only_contact_concentration.csv',index=False)
 con.append(f'{ds}：贡献远距离 AF3-only 接触最多的 5 个界面，占全部该类 >15 Å 接触的 {100*counts.head(5).far_contacts.sum()/counts.far_contacts.sum():.1f}%。')
 x=c[c.dataset==ds].pivot(index='contact_class',columns='state',values='pooled_fraction');change=100*(x.AF3_only-x.shared)
 lines.append(f'\n{ds} 化学探索：AF3-only 相对共享接触，polar–charged 比例变化 {change["polar_charged"]:+.2f} 个百分点，hydrophobic–hydrophobic 变化 {change["hydrophobic_hydrophobic"]:+.2f} 个百分点。池化和界面等权比例均已保存；这是组成描述，未做富集显著性断言。')
 fig,axes=plt.subplots(1,2,figsize=(11,4))
 for ax,(state,col) in zip(axes,[('AF3_only','experimental_d_min'),('Experimental_only','af3_d_min')]):
  v=d.loc[d.classification.eq(state),col];counts=np.bincount(np.searchsorted([9,10,12,15],v,side='left'),minlength=5)
  ax.bar(['(8,9]','(9,10]','(10,12]','(12,15]','>15'],100*counts/len(v));ax.set(xlabel='Distance in non-contact structure (Å)',ylabel='Contacts (%)',title=state)
 fig.suptitle(ds);fig.tight_layout();fig.savefig(P/'figures'/f'{ds}_noncontact_distance_histogram.png',dpi=180);plt.close(fig)
lines+=['\n### 远距离接触集中程度\n'+'\n'.join(con)]
f=P/'CONTACT_CORRESPONDENCE_SUMMARY.md';f.write_text(f.read_text()+'\n\n'.join(lines)+'\n')
