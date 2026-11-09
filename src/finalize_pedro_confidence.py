"""Numerical interpretation and reconciliation of the saved confidence analysis."""
from pathlib import Path
import pandas as pd
import numpy as np
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'results/local_interface_geometry/pedro_confidence_analysis'
s=pd.read_csv(P/'pedro_confidence_interface_table.csv',low_memory=False);nz=s[s.n_contacts>0]
d=pd.read_csv(P/'pedro_confidence_decile_summary.csv');c=pd.read_csv(P/'confidence_correlations.csv');q=pd.read_csv(P/'pedro_qc_confidence_summary.csv');z=pd.read_csv(P/'pedro_zero_contact_confidence_summary.csv');sz=pd.read_csv(P/'pedro_confidence_size_stratified_summary.csv');context=pd.read_csv(P/'high_confidence_benchmark_comparison.csv')
assert len(s)==170644 and len(nz)==54109 and s.interface_id.is_unique
assert d.n_interfaces.sum()==len(nz) and d.n_contacts.sum()==nz.n_contacts.sum() and z.n_interfaces.sum()==len(s) and z.n_zero.sum()==116535
assert np.allclose(nz.fraction_3_4,nz.count_3_4/nz.n_contacts)
assert s.loc[s.n_contacts==0,['fraction_3_4','median_d_min','fraction_lt_2']].isna().all().all()
# Explicitly separate all-interface confidence from conditional nonzero confidence in the QC audit.
q['nonzero_confidence_median']=[nz.best_model_pair_iptm.median(),nz.loc[~nz.existing_severe_overlap,'best_model_pair_iptm'].median(),nz.loc[nz.existing_severe_overlap,'best_model_pair_iptm'].median()]
q.to_csv(P/'pedro_qc_confidence_summary.csv',index=False)
r=c[(c.version=='raw')&(c.x=='best_model_pair_iptm')&(c.y=='fraction_3_4')].iloc[0];short=c[(c.version=='raw')&(c.x=='best_model_pair_iptm')&(c.y=='fraction_lt_2')].iloc[0];low=d.iloc[0];high=d.iloc[-1];validsz=sz[sz.version=='raw']
lines=['\n## 八个问题：结果解释\n']
lines.append(f'1. **置信度与近接触比例的关联**：在 {int(r.n_interfaces):,} 个非零界面中，所选模型链对 ipTM 与 (3,4] Å 比例的 Spearman ρ={r.spearman:.4f}，Pearson r={r.pearson:.4f}。这是界面层面的关联，不能解释为置信度造成几何变化。')
lines.append(f'2. **效应大小**：最低可用分位组 ipTM 范围 {low.confidence_min:g}–{low.confidence_max:g}，最高组 {high.confidence_min:g}–{high.confidence_max:g}。近接触比例逐界面中位数由 {100*low.interface_close_median:.2f}% 到 {100*high.interface_close_median:.2f}%（变化 {100*(high.interface_close_median-low.interface_close_median):+.2f} 个百分点）；池化比例由 {100*low.pooled_close_fraction:.2f}% 到 {100*high.pooled_close_fraction:.2f}%。应结合全部分位组判断是否单调。')
lines.append(f'3. **高置信度与实验比较**：最高分位组池化比例 {100*high.pooled_close_fraction:.2f}%，逐界面中位数 {100*high.interface_close_median:.2f}%。实验池化参考：FoldBench 14.36%、SAFES 13.62%、Mgen 原始 10.75%、Mgen 敏感性 12.11–12.45%。最高组是相对分位数，范围较宽，不能称作全部经校准的高可信结构。')
lines.append('4. **控制界面大小后的描述性关联**：相同接触数四分位层内的 Spearman ρ 为 '+', '.join(f'Q{int(v.size_quartile)} ({int(v.min_contacts)}–{int(v.max_contacts)} contacts): {v.spearman:.4f}' for v in validsz.itertuples())+'。尺寸层仍有剩余跨度；这是分层检查，不是因果调整。QC 通过子集的同一分层结果也已保存。')
lines.append(f'5. **是否由严重重叠主导**：原始与既有 QC 通过数据的池化近接触比例分别为 {100*q.iloc[0].pooled_close_fraction:.3f}% 和 {100*q.iloc[1].pooled_close_fraction:.3f}%；严重重叠组为 {100*q.iloc[2].pooled_close_fraction:.3f}%。QC 前后总体差异很小，不支持严重重叠是整体近接触缺口的主要来源。高 (3,4] Å 比例没有被定义为 QC 失败。')
lines.append(f'6. **极短距离与置信度**：ipTM 与 <2 Å 接触比例的 Spearman ρ={short.spearman:.4f}，Pearson r={short.pearson:.4f}。需结合尺寸、QC 分层以及大量零短接触比例解释；没有自动将 <2 Å 定义为碰撞。')
lines.append(f'7. **零接触与置信度**：按全体界面定义的最低/最高置信度分位组，零接触率分别为 {100*z.iloc[0].zero_fraction:.2f}%/{100*z.iloc[-1].zero_fraction:.2f}%。这些分位边界与非零几何分析不同。零接触的条件几何始终为 NA，没有置零参与距离分析。')
lines.append('8. **解释边界**：最终模式判断应同时参考关联方向、幅度、尺寸分层和高分位相对实验的缺口。尺寸与置信度可能共同关联于接触形成及局部几何；本分析没有识别因果机制。多个界面来自同一预测作业，且存在重复蛋白对；不把接触数量或界面数量视为独立重复的实验样本数。')
lines.append('\n## 模型选择及评分限制\n全部既有记录 selected_model=0；仅凭这些结果无法比较不同所选模型索引或判断换模型是否改变几何。ranking score 有负值，沿用原始数据而不截断或重定义；它属于复合物级指标，不作为主链对置信度。mean_pair_iptm 使用跨模型平均，而本几何只属于所选模型，其相关性只能作为辅助稳健性检查。')
lines.append('\n## QC 条件置信度补充\nconfidence_median 对 QC 表中所有界面计算，新增 nonzero_confidence_median 则仅针对几何有定义的非零界面，避免两种分母混淆。\n\n```\n'+q.to_string(index=False)+'\n```')
f=P/'PEDRO_CONFIDENCE_ANALYSIS.md';f.write_text(f.read_text()+'\n\n'.join(lines)+'\n')
for file,h in json.loads((P/'input_sha256.json').read_text()).items():assert hashlib.sha256((ROOT/file).read_bytes()).hexdigest()==h
(P/'validation.json').write_text(json.dumps(dict(status='passed',all_interfaces=170644,nonzero=54109,zero=116535,old_inputs_unchanged=True,decile_counts_reconciled=True,close_fractions_reconciled=True,zero_geometry_NA=True),indent=2))
print('\n\n'.join(lines[:9]))
