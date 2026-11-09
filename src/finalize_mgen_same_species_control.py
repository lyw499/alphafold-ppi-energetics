"""Interpret new same-species outputs; no coordinate access."""
from pathlib import Path
import pandas as pd
P=Path(__file__).resolve().parents[1]/'results/local_interface_geometry/mgen_same_species_control'
a=pd.read_csv(P/'mgen_same_species_comparison_summary.csv');b=pd.read_csv(P/'mgen_reference_sensitivity_summary.csv');i=pd.read_csv(P/'interface_level_distributions.csv');s=pd.read_csv(P/'mgen_experimental_interface_summary.csv')
text='\n## 八个问题的直接回答\n\n'
text+='1. **同种实验是否类似 Pedro？** 两者中位数接近，但局部距离分箱不同。不能将相近中位数解释为相同几何。下表同时给出峰区和较远区比例。\n\n'
cols=['dataset','median','fraction_3_4','fraction_4_5','fraction_5_6','fraction_6_8'];t=a[cols].copy();t.iloc[:,2:]=100*t.iloc[:,2:]
text+='比例单位为 %；中位数单位为 Å。\n\n```\n'+t.to_string(index=False)+'\n```\n\n'
text+='2. **是否也有弱 3–4 Å 峰？** 原始 Mgen 实验 (3,4] Å 比例 10.75%，低于 FoldBench/SAFES 实验的 14.36%/13.62%，但高于 Pedro 的 8.28%。因此有部分相似性，程度并不相同。\n\n'
text+='3. **是否更像一般实验基准？** 针对 (3,4] Å：去掉预先标记的 8PBZ/8PC1 后为 12.11%，再去掉 8PC0 后为 12.45%，接近一般实验而非 Pedro。但不同分箱、完整曲线可能给出不同相似程度，不应据此把参考集定义为通用标准。\n\n'
text+='4. **重复观测的影响？** 以下同时列出装配等权与 PDB 等权。等权消除一个装配仅因链对/接触数多而主导，但仍保留同一复合物跨 PDB 的相关性；不等于获得独立生物学样本。\n\n'
t=b[['view','weighting','n_pairs','n_nonzero','n_pdb','n_assemblies','fraction_3_4','fraction_4_5','fraction_5_6','fraction_6_8']].copy();t.iloc[:,6:]*=100
text+='```\n'+t.to_string(index=False)+'\n```\n\n'
text+='5. **低分辨率是否影响答案？** 排除既有低分辨率/几何疑点结构后，3–4 Å 比例升高，原始“弱峰”部分由这些结构贡献。因此不能将弱峰直接解释为物种特性。主结果没有删除这些结构。\n\n'
text+='6. **界面层面如何？** 下表用等权界面统计；零接触界面不虚构距离或条件比例。界面大小、总界面残基数和覆盖比例另见 interface_level_distributions.csv。\n\n'
t=i[i.dataset.isin(['Mgen Experimental','Pedro AF3']) & i.metric.isin(['median_d_min','fraction_3_4','fraction_4_5','fraction_5_6','fraction_6_8','fraction_lt_2_0'])]
text+='```\n'+t.to_string(index=False)+'\n```\n\n'
text+='7. **极短接触是否主要来自已有问题结构？** 同种实验 10,528 个接触中有 26 个 <2 Å（0.247%）；25 个在 8PBZ，另 1 个在 8PBY。排除 8PBZ/8PC1 后只剩 1/6,491（0.0154%），再排除 8PC0 后为 1/6,054（0.0165%）。Pedro 原始为 0.5034%，已有 QC 后仍为 0.4337%；FoldBench/SAFES AF3 为 0.0588%/0.0506%。不能把所有 <2 Å 一律称作碰撞；接触层面的比例不是独立样本显著性证据。\n\n'
text+='8. **是否支持物种解释？** 证据不足。参考主集只有 7 个 PDB、8 个装配，23 个链对中 2 个零接触，而且基本集中于同一黏附复合物及相关研究；分辨率、建模方式、蛋白组成和装配重复混杂。弱峰随敏感性处理而减弱，不能归因为 M. genitalium 的普遍几何特征，也不能定义该物种真实全局 PPI 分布。\n\n'
text+='## 化学与短距离补充\n\n'
chem=pd.read_csv(P/'chemistry_comparison.csv');text+='```\n'+chem[chem.dataset.isin(['Mgen Experimental','Pedro AF3'])].to_string(index=False)+'\n```\n'
text+='化学分类只是组成描述；参考集以少数黏附复合物为主，与全蛋白组预测不匹配，不应把比例差异直接归因于物种或 AF3 的化学偏好。\n'
f=P/'MGEN_SAME_SPECIES_CONTROL.md';f.write_text(f.read_text()+text)
print(t.to_string(index=False))
