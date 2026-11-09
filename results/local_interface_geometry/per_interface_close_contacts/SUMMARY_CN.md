# 每界面 (3,4] Å 接触数量：结果与统计解释

## 审计与范围

复用已验证汇总与残基距离表；没有结构重算或覆盖旧输出。边界为3 < d_min ≤4 Å，分母接触为≤8 Å，一对残基仅计一次；无口径变化。FoldBench主比较132个完整已配对界面，严格配对控制109个；SAFES50个链对。

{
  "pedro_archive_job_records": 2156,
  "pedro_job_names": 2033,
  "pedro_archives": 70,
  "pedro_archive_components": 60,
  "proteins": 476,
  "unique_protein_pairs": 113050,
  "extra_protein_pair_representations": 57594,
  "possible_repeat_keys": 7283,
  "possible_repeat_extra_rows": 9574,
  "possible_repeat_different_count": 2251
}

Pedro跨预测重复蛋白对与跨归档job均存在；同一源文件/链对无重复。可能重复表示的详细表及只取首个表示的敏感性已保存，未据此删除原始数据。

## 1. 描述统计（分位数为逆经验CDF）

| version | label | scope | n | n_clusters | zero_8A_percent | zero_close_percent | close_count_median | close_count_q25 | close_count_q75 | close_count_mean | pooled_close_fraction | close_fraction_median |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| raw | FoldBench AF3 | zero_inclusive | 132 | 108 | 0.75758 | 1.5152 | 37 | 23 | 50 | 42.144 | 0.13623 | 0.12632 |
| raw | FoldBench AF3 | nonzero_only | 131 | 107 | 0 | 0.76336 | 37 | 24 | 50 | 42.466 | 0.13623 | 0.12632 |
| raw | FoldBench Experimental | zero_inclusive | 132 | 108 | 0 | 0 | 32 | 20 | 42 | 38.667 | 0.14362 | 0.13077 |
| raw | FoldBench Experimental | nonzero_only | 132 | 108 | 0 | 0 | 32 | 20 | 42 | 38.667 | 0.14362 | 0.13077 |
| raw | Pedro AF3 | zero_inclusive | 170644 | 2033 | 68.291 | 75.719 | 0 | 0 | 0 | 3.2384 | 0.082829 | 0.068966 |
| raw | Pedro AF3 | nonzero_only | 54109 | 2033 | 0 | 23.425 | 5 | 1 | 13 | 10.213 | 0.082829 | 0.068966 |
| raw | SAFES AF3 | zero_inclusive | 50 | 19 | 14 | 16 | 16 | 7 | 34 | 25.88 | 0.13093 | 0.12281 |
| raw | SAFES AF3 | nonzero_only | 43 | 19 | 0 | 2.3256 | 25 | 11 | 38 | 30.093 | 0.13093 | 0.12281 |
| raw | SAFES Experimental | zero_inclusive | 50 | 19 | 16 | 16 | 16 | 6 | 35 | 24.3 | 0.13618 | 0.12766 |
| raw | SAFES Experimental | nonzero_only | 42 | 19 | 0 | 0 | 20 | 12 | 38 | 28.929 | 0.13618 | 0.12766 |
| existing_qc | Pedro AF3 | zero_inclusive | 170180 | 2033 | 68.477 | 75.925 | 0 | 0 | 0 | 3.151 | 0.082702 | 0.068966 |
| existing_qc | Pedro AF3 | nonzero_only | 53645 | 2033 | 0 | 23.628 | 5 | 1 | 13 | 9.996 | 0.082702 | 0.068966 |

## 2. 聚类效应与检验

方向均为Pedro−实验。主要p值检验平均log(1+数量)差异，不是中位数或所有分布形状的联合检验；4项Holm校正。效应区间是边际95%区间。KS仅报告D，未用接触级/朴素界面级p值。

| reference | scope | median_count_difference | median_ci_low | median_ci_high | cliffs_delta | cliffs_ci_low | cliffs_ci_high | mean_log1p_difference | logdiff_ci_low | logdiff_ci_high | p_cluster_approx | p_holm | ks_D_descriptive |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| FoldBench Experimental | zero_inclusive | -32 | -36 | -30 | -0.91871 | -0.93889 | -0.8978 | -2.8541 | -3.0206 | -2.6956 | 1.2145e-58 | 4.8581e-58 | 0.79378 |
| SAFES Experimental | zero_inclusive | -16 | -29 | -12 | -0.70158 | -0.86925 | -0.58297 | -2.0551 | -2.586 | -1.7022 | 3.0351e-08 | 3.0351e-08 | 0.64195 |
| FoldBench Experimental | nonzero_only | -27 | -31 | -25 | -0.74365 | -0.80764 | -0.67812 | -1.7058 | -1.8752 | -1.5367 | 4.8861e-38 | 1.4658e-37 | 0.63048 |
| SAFES Experimental | nonzero_only | -15 | -27 | -11 | -0.62619 | -0.71304 | -0.56821 | -1.3998 | -1.5664 | -1.2823 | 1.1163e-13 | 2.2327e-13 | 0.51054 |

QC及更粗归档聚类完整结果见 cluster_aware_comparisons.csv：

| version | grouping | reference | scope | clusters_pedro | median_count_difference | cliffs_delta | logdiff_ci_low | logdiff_ci_high | p_holm |
|---|---|---|---|---|---|---|---|---|---|
| raw | archive_component | FoldBench Experimental | zero_inclusive | 60 | -32 | -0.91871 | -3.0192 | -2.695 | 4.7434e-58 |
| raw | archive_component | SAFES Experimental | zero_inclusive | 60 | -16 | -0.70158 | -2.5867 | -1.7005 | 3.0342e-08 |
| raw | archive_component | FoldBench Experimental | nonzero_only | 60 | -27 | -0.74365 | -1.8735 | -1.5365 | 1.5502e-37 |
| raw | archive_component | SAFES Experimental | nonzero_only | 60 | -15 | -0.62619 | -1.5622 | -1.2828 | 2.3716e-13 |
| existing_qc | job | FoldBench Experimental | zero_inclusive | 2033 | -32 | -0.92043 | -3.025 | -2.6999 | 3.8354e-58 |
| existing_qc | job | SAFES Experimental | zero_inclusive | 2033 | -16 | -0.70386 | -2.5908 | -1.7076 | 2.8848e-08 |
| existing_qc | job | FoldBench Experimental | nonzero_only | 2033 | -27 | -0.74757 | -1.8842 | -1.5526 | 7.9682e-38 |
| existing_qc | job | SAFES Experimental | nonzero_only | 2033 | -15 | -0.63134 | -1.5751 | -1.2931 | 1.9133e-13 |

这些结果依赖job/PDB聚类独立性的近似；476种蛋白跨job重复仍有交叉依赖，且实验样本为非随机选择。应称条件性/探索性统计证据，不宜宣称已解决所有相关性或证明AF3错误。零膨胀差异尤其受数据集选择影响。

## 3. 面积、总接触数与早期结果

Pedro仅有选定小样本的ΔSASA/2，无法对全量Pedro完成Å²校正。现有面积可用数量见 grouping_and_area_audit.csv。图7仅作配对基准面积背景，不构成面积调整。总接触数匹配不是面积匹配，也不提供独立因果贡献。

| version | reference | n | experimental_mean_count | matched_pedro_mean_count | mean_count_difference | median_count_difference | within_support |
|---|---|---|---|---|---|---|---|
| existing_qc | FoldBench Experimental | 132 | 38.667 | 23.104 | -15.563 | -11.613 | 132 |
| existing_qc | SAFES Experimental | 42 | 28.929 | 17.646 | -11.282 | -9.452 | 42 |
| raw | FoldBench Experimental | 132 | 38.667 | 23.099 | -15.568 | -11.613 | 132 |
| raw | SAFES Experimental | 42 | 28.929 | 17.62 | -11.309 | -9.5346 | 42 |

原8.28%/14.36%/13.62%是接触加权pooled比例；本次数量分布是界面绝对接触数，零接触界面也贡献观测。界面比例及总接触数匹配分别回答局部比例和条件大小问题，不能互相替代。旧匹配结果路径：./results/local_interface_geometry/size_matched_packing/size_matched_effect_differences.csv。

## 4. 独立报告的严格配对控制

| dataset | n | n_pdb | median_difference | mean_difference | ci_low | ci_high |
|---|---|---|---|---|---|---|
| FoldBench | 109 | 88 | 1 | 3.2202 | 0 | 5.025 |
| SAFES | 50 | 19 | 0 | 1.58 | 0 | 1 |

此处方向为AF3−Experimental，置信区间按PDB重采样配对差，不能与Pedro跨数据集差当作同一问题。

## 5. 逐图说明

- `./results/local_interface_geometry/per_interface_close_contacts/figures/01_close_count_ecdf.png`：主要数量分布：左含零，右仅总8 Å接触>0；对称log轴保留0，原始与旧QC同时显示。
- `./results/local_interface_geometry/per_interface_close_contacts/figures/02_close_fraction_ecdf.png`：次要界面比例ECDF，仅比例有定义界面；不把零接触的比例当0。
- `./results/local_interface_geometry/per_interface_close_contacts/figures/03_close_count_vs_total_contacts.png`：全部界面数量—总接触数散点；长尾坐标含0。重叠点较多，散点密度不代表独立生物重复。
- `./results/local_interface_geometry/per_interface_close_contacts/figures/04_cluster_effect_sizes.png`：主要中位数量差与Cliff效应及job/PDB bootstrap区间，区分统计差异幅度和p值。
- `./results/local_interface_geometry/per_interface_close_contacts/figures/05_paired_count_controls.png`：FoldBench109/SAFES50严格配对控制，y=x显示相同蛋白链对的计数对应。
- `./results/local_interface_geometry/per_interface_close_contacts/figures/06_count_size_matching_descriptive.png`：非零界面总接触数匹配后的数量差，描述性辅助；不做重复匹配样本的朴素配对检验。
- `./results/local_interface_geometry/per_interface_close_contacts/figures/07_area_context_paired_only.png`：仅已有完整面积的配对控制中，接触数与ΔSASA/2的关系，不能代替Pedro面积校正。

## 6. 可向Amelie报告及限制

- 对 FoldBench Experimental，zero_inclusive：数量中位差 -32，95%区间 [-36, -30]；Cliff δ=-0.919；平均log(1+数量)差的近似聚类检验Holm p=4.86e-58。
- 对 SAFES Experimental，zero_inclusive：数量中位差 -16，95%区间 [-29, -12]；Cliff δ=-0.702；平均log(1+数量)差的近似聚类检验Holm p=3.04e-08。
- 对 FoldBench Experimental，nonzero_only：数量中位差 -27，95%区间 [-31, -25]；Cliff δ=-0.744；平均log(1+数量)差的近似聚类检验Holm p=1.47e-37。
- 对 SAFES Experimental，nonzero_only：数量中位差 -15，95%区间 [-27, -11]；Cliff δ=-0.626；平均log(1+数量)差的近似聚类检验Holm p=2.23e-13。

可以讨论当前所采样数据中每界面近接触数量分布的差异，并明确包含零和非零条件结果。不能把跨数据集显著差异直接解释为AF3精度缺陷、物种因果效应或结合强度差异。尚未完成全量Å²校正，也未完全解决跨job共享蛋白的相关性。主要分析以效应和不确定性为先，配对基准作为不同问题的控制。

所有输入路径及SHA256在INPUT_MANIFEST.json；脚本、标准表、聚类比较、重复表示审计与图均在本目录。

## 7. 补充分布检验及最终核验

除平均log(1+数量)检验外，另以整个数量ECDF的最大差异做中心化聚类bootstrap KS检验（2,000次，独立种子20261025）；它不是假定所有界面独立的标准KS检验。详见DISTRIBUTIONAL_TEST_METHOD.md。

- FoldBench Experimental / zero_inclusive：KS D=0.794，bootstrap超越次数 0/2000，Monte Carlo p=0.00050，补充4项Holm校正p=0.00200。
- SAFES Experimental / zero_inclusive：KS D=0.642，bootstrap超越次数 0/2000，Monte Carlo p=0.00050，补充4项Holm校正p=0.00200。
- FoldBench Experimental / nonzero_only：KS D=0.630，bootstrap超越次数 0/2000，Monte Carlo p=0.00050，补充4项Holm校正p=0.00200。
- SAFES Experimental / nonzero_only：KS D=0.511，bootstrap超越次数 0/2000，Monte Carlo p=0.00050，补充4项Holm校正p=0.00200。

0次超越表示分辨率受2000次模拟限制，不能将p写成0，也不应把显示值当作精确真实p。分布检验仍依赖job/PDB独立性近似，不能消除共享蛋白和非随机选择的限制。

最终检查：主表171,008行（170,644 Pedro + 264 FoldBench + 100 SAFES）；7张PNG和7张PDF已人工检查；所有实际输入的SHA256前后相同；配对计数、区间边界、零值/NA、数值算法检查通过。主任务66900648；补充分布任务66900657。旧结果未覆盖。
