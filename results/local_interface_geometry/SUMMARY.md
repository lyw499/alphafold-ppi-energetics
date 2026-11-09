# Experimental vs AF3：局部界面几何

仅分析 FoldBench 和 SAFES。Pedro未运行；历史界面大小结果未修改。主观测为跨链残基对的最小非氢原子距离，不是原子对频数。

## 方法与口径

- 对每条链建立重原子KD-tree，只产生≤8 Å临时原子邻居，然后按残基对保留最小距离及其原子。没有构建/保存整个链的所有长距离原子矩阵。
- 沿用有坐标且residue.id[0]为空的残基、首个model和Bio.PDB默认altloc。排除H和D；已确认全部分析链无D，FoldBench新接触数逐接口与历史残基结果一致。HETATM编码修饰残基/配体仍不纳入历史定义。
- 密度图共同分箱为0–8 Å，每0.25 Å；ECDF使用实际距离。pooled图按接触数加权，仅作描述，不把残基对当独立实验样本。
- 每接口的分箱为[0,3]、(3,4]、(4,5]、(5,6]、(6,8] Å，互斥且总和1；另保留≤4/5/6 Å累积分数，不使用恒等于1的≤8 Å分数。
- 零接触接口保留在链对审计与接口汇总中；其距离均值/中位数和比例未定义（NaN），不能填0。接口散点仅使用两侧相关统计量均定义的有效配对；每项报告排除的未定义数。
- 每接口统计保留配对相关性/差值；额外先在每PDB内平均接口差，再对PDB等权。95%区间为2000次PDB簇bootstrap（固定种子），不报告把全部残基对当独立样本产生的p值；不同PDB仍可能同源或来自同一实验体系。
- 所有短距离阈值仅诊断，不增加自动QC删除。小于2 Å可能涉及共价/二硫键、切分链或局部碰撞，不能仅凭阈值断言unphysical。extreme表提供残基、原子和既有QC/映射状态。

## 数据集与映射

| dataset | n_audited_chain_pairs | n_included_interfaces | n_structures | n_unavailable_or_unreliable | n_ordered_unambiguous | n_historical_size_mapping_stable | n_local_spectrum_stable |
| --- | --- | --- | --- | --- | --- | --- | --- |
| foldbench | 279 | 132 | 108 | 147 | 109 | 118 | 126 |
| safes | 64 | 50 | 19 | 14 | 50 | 50 | 50 |

FoldBench复用已审计的132条映射。147个AF3缺失案例只进入审计，不填零。完整集合保留几何并列映射；分别提供历史界面大小稳定子集、无歧义有序链映射子集、局部距离谱稳定子集。历史大小稳定不保证局部几何稳定；并列链对的替代距离汇总见mapping_tie_geometry_sensitivity.csv。谱一致判据是接触数相同且排序后的距离差≤1e-6 Å，仅为数值一致性检查，不是结构QC阈值。

SAFES按本次要求枚举全部可靠对应蛋白链对，**不再限于旧10个两链案例**。沿用序列identity≥0.95、实验观测覆盖率≥0.90和全复合体CA几何的一对一链映射。若所有并列映射给出同一个无序AF3链对，可用于链对距离分布；若对应到不同AF3链对则不纳入配对主分析，原因保留。此枚举是所有可对应的链对，不代表SAFES官方预期结合伙伴清单；抗体亚基/多链结构内部接触也可能在内。

逐残基配对更严格：有序链对应必须唯一；使用与链映射相同参数的全局序列比对，仅保留所有最优比对共同指定、氨基酸身份相同且非X的位置。最优比对超过128种则跳过并记录，不凭残基编号相同强行配对。只计算两侧接触集合的可映射并集；另一侧即使>8 Å也仅计算该指定残基对的小型原子距离矩阵，作为接触丢失/新增的证据。此配对结果有保守映射与界面选择偏差。

## 汇总距离（pooled，仅描述）

| dataset | scope | structure_type | n_interfaces | n_complexes | n_zero_interfaces | n_nonzero_interfaces | n | mean | median | q25 | q75 | min | max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| foldbench | full_matched | Experimental | 132 | 108 | 0 | 132 | 35539 | 6.0228 | 6.4772 | 4.7989 | 7.283 | 1.3524 | 8 |
| foldbench | full_matched | AF3 | 132 | 108 | 1 | 131 | 40836 | 5.9988 | 6.4451 | 4.7802 | 7.2835 | 1.3996 | 7.9999 |
| safes | full_matched | Experimental | 50 | 19 | 8 | 42 | 8922 | 6.0633 | 6.5257 | 4.8945 | 7.3166 | 1.9853 | 7.9999 |
| safes | full_matched | AF3 | 50 | 19 | 7 | 43 | 9883 | 6.0595 | 6.5257 | 4.9199 | 7.3131 | 1.5859 | 7.9991 |

## 主要配对接口结果

| dataset | scope | metric | n_paired | n_unavailable_summary | n_PDB_clusters | experimental_median | af3_median | pearson | spearman | median_delta | mean_delta | delta_q25 | delta_q75 | n_AF3_shorter_or_lower | n_AF3_longer_or_higher | PDB_equal_weight_mean_delta | cluster_bootstrap_mean_delta_low | cluster_bootstrap_mean_delta_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| foldbench | full_matched | median_d_min | 131 | 1 | 107 | 6.5187 | 6.5099 | 0.95875 | 0.57504 | -0.00067855 | 0.014758 | -0.077412 | 0.098085 | 66 | 65 | -0.0031861 | -0.033004 | 0.028537 |
| foldbench | full_matched | fraction_le_4 | 131 | 1 | 107 | 0.15789 | 0.15847 | 0.89608 | 0.69755 | -0.0025705 | -0.0082132 | -0.024668 | 0.01302 | 74 | 57 | -0.0041224 | -0.010531 | 0.0024932 |
| safes | full_matched | median_d_min | 42 | 8 | 19 | 6.5369 | 6.5471 | 0.75969 | 0.6822 | -0.0051017 | 0.011691 | -0.083104 | 0.055023 | 21 | 21 | 0.010228 | -0.054694 | 0.093247 |
| safes | full_matched | fraction_le_4 | 42 | 8 | 19 | 0.16062 | 0.15569 | 0.75217 | 0.74658 | -0.00097033 | -0.0029873 | -0.012215 | 0.012574 | 21 | 21 | -0.0064431 | -0.016844 | 0.0013346 |

## FoldBench映射敏感性

| dataset | scope | metric | n_paired | n_unavailable_summary | n_PDB_clusters | experimental_median | af3_median | pearson | spearman | median_delta | mean_delta | delta_q25 | delta_q75 | n_AF3_shorter_or_lower | n_AF3_longer_or_higher | PDB_equal_weight_mean_delta | cluster_bootstrap_mean_delta_low | cluster_bootstrap_mean_delta_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| foldbench | full_matched | median_d_min | 131 | 1 | 107 | 6.5187 | 6.5099 | 0.95875 | 0.57504 | -0.00067855 | 0.014758 | -0.077412 | 0.098085 | 66 | 65 | -0.0031861 | -0.033004 | 0.028537 |
| foldbench | legacy_size_mapping_stable | median_d_min | 117 | 1 | 96 | 6.5128 | 6.5038 | 0.95956 | 0.59058 | -0.0051427 | 0.019456 | -0.077315 | 0.10456 | 60 | 57 | -0.00090072 | -0.033141 | 0.033227 |
| foldbench | unambiguous_chain_mapping | median_d_min | 109 | 0 | 88 | 6.5078 | 6.5026 | 0.96094 | 0.58138 | 0.00018385 | 0.022782 | -0.077315 | 0.10827 | 54 | 55 | 0.0013688 | -0.031517 | 0.035674 |
| foldbench | local_spectrum_stable | median_d_min | 125 | 1 | 104 | 6.5162 | 6.5038 | 0.95861 | 0.57731 | -0.0051427 | 0.014715 | -0.079605 | 0.10436 | 64 | 61 | -0.0050329 | -0.040614 | 0.027597 |

## 短距离诊断

| dataset | structure_type | threshold_strict_lt_A | n_contacts_below | n_all_contacts | fraction_contacts_below | n_interfaces_with_short | n_all_interfaces | n_structures_with_short |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| foldbench | AF3 | 0.5 | 0 | 40836 | 0 | 0 | 132 | 0 |
| foldbench | AF3 | 1 | 0 | 40836 | 0 | 0 | 132 | 0 |
| foldbench | AF3 | 1.5 | 2 | 40836 | 4.8976e-05 | 2 | 132 | 2 |
| foldbench | AF3 | 2 | 24 | 40836 | 0.00058772 | 15 | 132 | 15 |
| foldbench | Experimental | 0.5 | 0 | 35539 | 0 | 0 | 132 | 0 |
| foldbench | Experimental | 1 | 0 | 35539 | 0 | 0 | 132 | 0 |
| foldbench | Experimental | 1.5 | 1 | 35539 | 2.8138e-05 | 1 | 132 | 1 |
| foldbench | Experimental | 2 | 2 | 35539 | 5.6276e-05 | 2 | 132 | 2 |
| safes | AF3 | 0.5 | 0 | 9883 | 0 | 0 | 50 | 0 |
| safes | AF3 | 1 | 0 | 9883 | 0 | 0 | 50 | 0 |
| safes | AF3 | 1.5 | 0 | 9883 | 0 | 0 | 50 | 0 |
| safes | AF3 | 2 | 5 | 9883 | 0.00050592 | 5 | 50 | 4 |
| safes | Experimental | 0.5 | 0 | 8922 | 0 | 0 | 50 | 0 |
| safes | Experimental | 1 | 0 | 8922 | 0 | 0 | 50 | 0 |
| safes | Experimental | 1.5 | 0 | 8922 | 0 | 0 | 50 | 0 |
| safes | Experimental | 2 | 1 | 8922 | 0.00011208 | 1 | 50 | 1 |

阈值均为严格小于。分母是本组≤8 Å残基接触数；同时报告涉及接口/结构数量。最短接触另列于 `qc/extreme_short_distance_cases.csv`，没有删除可计算结构。

## 探索性残基化学

分类：疏水A,V,I,L,M,F,W,Y；极性不带电S,T,N,Q,C；正电K,R；负电D,E；特殊G,P,H；非标准/未知为other。H的质子化状态不确定，因此不硬归正电；Y/C的分组是简化约定。正负电指常用残基类别，并未计算实际质子化或能量。
接触类别互斥：疏水–疏水、异号电荷、同号电荷、极性–极性（仅不带电组）、疏水–极性、疏水–带电、极性–带电、含特殊/其他。全部比例合计1。

| dataset | structure_type | contact_class | n_contacts | pooled_fraction | interface_equal_weight_mean_fraction | n_nonzero_interfaces |
| --- | --- | --- | --- | --- | --- | --- |
| foldbench | AF3 | hydrophobic_hydrophobic | 7939 | 0.19441 | 0.17541 | 131 |
| foldbench | AF3 | opposite_charge | 1586 | 0.038838 | 0.047861 | 131 |
| foldbench | AF3 | same_charge | 1333 | 0.032643 | 0.037871 | 131 |
| foldbench | AF3 | polar_polar | 1925 | 0.04714 | 0.047198 | 131 |
| foldbench | AF3 | hydrophobic_polar | 6096 | 0.14928 | 0.14157 | 131 |
| foldbench | AF3 | hydrophobic_charged | 7474 | 0.18302 | 0.17965 | 131 |
| foldbench | AF3 | polar_charged | 3707 | 0.090778 | 0.090018 | 131 |
| foldbench | AF3 | special_or_other | 10776 | 0.26388 | 0.28042 | 131 |
| foldbench | Experimental | hydrophobic_hydrophobic | 6687 | 0.18816 | 0.16995 | 132 |
| foldbench | Experimental | opposite_charge | 1413 | 0.039759 | 0.045689 | 132 |
| foldbench | Experimental | same_charge | 1266 | 0.035623 | 0.040982 | 132 |
| foldbench | Experimental | polar_polar | 1646 | 0.046315 | 0.049579 | 132 |
| foldbench | Experimental | hydrophobic_polar | 5509 | 0.15501 | 0.1458 | 132 |
| foldbench | Experimental | hydrophobic_charged | 6489 | 0.18259 | 0.17623 | 132 |
| foldbench | Experimental | polar_charged | 3399 | 0.095641 | 0.093954 | 132 |
| foldbench | Experimental | special_or_other | 9130 | 0.2569 | 0.27782 | 132 |
| safes | AF3 | hydrophobic_hydrophobic | 1510 | 0.15279 | 0.1379 | 43 |
| safes | AF3 | opposite_charge | 251 | 0.025397 | 0.024322 | 43 |
| safes | AF3 | same_charge | 155 | 0.015683 | 0.023461 | 43 |
| safes | AF3 | polar_polar | 779 | 0.078822 | 0.069511 | 43 |
| safes | AF3 | hydrophobic_polar | 2103 | 0.21279 | 0.20579 | 43 |
| safes | AF3 | hydrophobic_charged | 1331 | 0.13468 | 0.13261 | 43 |
| safes | AF3 | polar_charged | 986 | 0.099767 | 0.10724 | 43 |
| safes | AF3 | special_or_other | 2768 | 0.28008 | 0.29916 | 43 |
| safes | Experimental | hydrophobic_hydrophobic | 1414 | 0.15848 | 0.14347 | 42 |
| safes | Experimental | opposite_charge | 199 | 0.022304 | 0.021166 | 42 |
| safes | Experimental | same_charge | 130 | 0.014571 | 0.014329 | 42 |
| safes | Experimental | polar_polar | 702 | 0.078682 | 0.073295 | 42 |
| safes | Experimental | hydrophobic_polar | 1891 | 0.21195 | 0.21137 | 42 |
| safes | Experimental | hydrophobic_charged | 1162 | 0.13024 | 0.13127 | 42 |
| safes | Experimental | polar_charged | 873 | 0.097848 | 0.098343 | 42 |
| safes | Experimental | special_or_other | 2551 | 0.28592 | 0.30676 | 42 |

pooled比例会被大接口主导，表中同时提供接口等权均值；配对接口化学比例差与PDB等权汇总见paired_chemistry_statistics.csv。频率变化不等同于能量或偏好机制。

## 逐残基配对覆盖

- foldbench: 109 个接口有可映射接触并集，共40348对；共有接触18442、AF3丢失10226、AF3新增11680。残基定位审计及双方实际距离保存在各数据集子目录。
- safes: 43 个接口有可映射接触并集，共10282对；共有接触7800、AF3丢失1113、AF3新增1369。残基定位审计及双方实际距离保存在各数据集子目录。

## 输出与复现

- `foldbench/`、`safes/`：residue_contacts.csv（逐残基接触及最近原子）、interface_summaries.csv、chain_pair_audit.csv、paired_residue_distances.csv、residue_mapping_audit.csv、mapping_tie_geometry_sensitivity.csv及配对统计。
- `figures/foldbench/`、`figures/safes/`：密度、ECDF、接口中位距离/≤4 Å比例散点、探索性化学频率；PNG及PDF。
- `qc/`：短距离阈值汇总、逐接口异常尾部、极端残基/原子案例。
- `algorithm_validation.json`：三个真实小结构与穷举一致，及合成边界/氢排除/零接口验证。`validation.json`：完整结果一致性。
- 复现依次运行 `src/test_local_geometry.py`、`src/compute_local_geometry.py`、`src/report_local_geometry.py`。所有输出限于此新目录。


## 直接回答与解释边界

1. **总体分布是否相似？** 两组实验/AF3的pooled密度和ECDF都很接近。靠近8 Å的接触数量本来就受截断与空间邻域体积影响；图形相似不是正确界面或逐残基准确性的充分证据。
2. **有系统性变短或变长吗？** 当前逐接口中位距离差接近零；PDB等权平均差的bootstrap区间覆盖零，尚无一致方向的证据。下面给出实际效应量，不能以不显著宣称完全等价。

- foldbench: pooled中位数 Experimental 6.4772 Å / AF3 6.4451 Å；逐接口中位数差的中位数 -0.0007 Å，均值 +0.0148 Å；PDB等权均值 -0.0032 Å，95%簇bootstrap区间 [-0.0330, 0.0285] Å。
- safes: pooled中位数 Experimental 6.5257 Å / AF3 6.5257 Å；逐接口中位数差的中位数 -0.0051 Å，均值 +0.0117 Å；PDB等权均值 +0.0102 Å，95%簇bootstrap区间 [-0.0547, 0.0932] Å。

3. **接口层面对应程度？** FoldBench的中位距离Pearson=0.959、Spearman=0.575（131个两侧非零接口）；SAFES为0.760、0.682（42个）。FoldBench的高Pearson部分反映分离的低中位距离接口群，不能将其单独解读为所有接口都高度准确。未定义距离的零接口仍在审计中。

- foldbench: 在严格映射的接触集合中，逐接口实验接触保留率中位数=0.877；两侧共同接触的逐接口中位绝对距离误差之中位数=0.244 Å。共同接触分析受两侧≤8 Å选择限制；完整并集也保留丢失/新增和另一侧>8 Å的指定残基距离。
- safes: 在严格映射的接触集合中，逐接口实验接触保留率中位数=0.930；两侧共同接触的逐接口中位绝对距离误差之中位数=0.214 Å。共同接触分析受两侧≤8 Å选择限制；完整并集也保留丢失/新增和另一侧>8 Å的指定残基距离。

4. **AF3短距离是否更多？** FoldBench Experimental 2/35539（0.0056%）vs AF3 24/40836（0.0588%）小于2 Å；SAFES 1/8922（0.0112%）vs 5/9883（0.0506%）。四组均无<1 Å。FoldBench <1.5 Å为实验1条、AF3 2条；SAFES均0。
短距离不全部等于非物理碰撞：SAFES实验的1条和AF3的4/5条为SG–SG；FoldBench AF3的3/24条为SG–SG。剩余极短O–O等值得检查；实验9c90也有约1.35 Å的OE2–OE2接触。FoldBench两条实验<2 Å记录均涉及占有率0.51的B/B替代构象；Bio.PDB默认选定构象未做跨链altloc协调，不能仅凭该短距离认定实验生物学界面必然碰撞。不同altloc的最近原子距离另存short_distance_altloc_diagnostics.csv（没有替换主度量，也不据占有率假定构象必然互斥）。`qc/short_distance_chemical_context.csv`保留所选altloc/occupancy和可用的共价连接注释，未根据新阈值删行。

5. **化学类别频率是否明显不同？** pooled宽类别变化较小，最大绝对差见下；接口等权与配对变化另表报告，不能据这些频率推出能量偏好。

- foldbench: 最大pooled类别差 special_or_other: AF3−Experimental +0.698 个百分点。
- safes: 最大pooled类别差 special_or_other: AF3−Experimental -0.585 个百分点。

6. **两个数据集定性一致吗？** 都显示总体距离形状相近、逐接口中位距离无清晰单向偏移、原始<2 Å尾部AF3更多。但SAFES短距离主要是SG–SG，且样本较小、包含复合物内部链对；不能把两组短尾都归为同一类碰撞问题。
7. **是否值得扩展Pedro？** 可以用于描述接触几何和定位异常短距离，尤其作为已有严重重叠QC的补充；不能单凭几何分布相似就判断Pedro界面正确，也不能在缺乏一对一实验对应时照搬配对准确性结论。还应保留重复观测/接口大小与蛋白来源的层级差异。本轮没有运行Pedro。

## SAFES未纳入链对的具体序列证据

1OGA涉及A链的4个链对未达既有95%identity阈值；5E9D只有B链达标，因此全部10个链对未纳入。低于保守规则不等于断言它们绝对不可能建立生物学对应；本轮没有放宽阈值。所有候选比对保存在failed_target_all_chain_sequence_candidates.csv。

| structure_id | experimental_chain | af3_chain | experimental_coordinate_residues | af3_coordinate_residues | identity | experimental_observed_coverage | n_matched_CA | passes_existing_candidate_rules | failure_reasons |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1OGA | A | A | 276 | 275 | 0.93091 | 0.99638 | 256 | False | identity_below_0.95 |
| 5E9D | A | A | 266 | 275 | 0.92857 | 1 | 247 | False | identity_below_0.95 |
| 5E9D | C | C | 10 | 10 | 0.9 | 1 | 9 | False | identity_below_0.95 |
| 5E9D | D | D | 110 | 91 | 0.97753 | 0.80909 | 87 | False | AF3_coordinate_chain_too_short_for_90pct_experimental_coverage;experimental_observed_coverage_below_0.90 |
| 5E9D | E | E | 115 | 94 | 0.94565 | 0.8 | 87 | False | AF3_coordinate_chain_too_short_for_90pct_experimental_coverage;identity_below_0.95;experimental_observed_coverage_below_0.90 |

补充复现：主报告生成后运行 `src/audit_local_geometry_context.py`，再运行 `src/report_local_geometry.py --summary-only` 将化学背景与失败比对证据写入本报告。
