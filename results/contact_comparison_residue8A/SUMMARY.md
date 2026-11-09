# 8 Å 残基界面比较

状态：FoldBench 审计与配对分析完成；Pedro 既有残基结果及实验参考候选的分布比较完成。**完整 SAFES 预期界面基准尚待结合伙伴元数据确认**；仅提供明确标记的两蛋白链结构暂定子集，不能冒充全部20个案例的正式配对结果。

## 统一定义与统计单位

**8 Å heavy-atom cutoff 定义两个残基是否接触；主 contact_count 是不同链的唯一残基对数，不是重原子对数。** 每对残基只计一次。直接复用 `src/parse_pooled_inner_zip_8A.py::count_pair_contacts`（经 `foldbench_interfaces.metrics` 调用）。原函数选择 `residue.id[0] == " "` 的有坐标残基，排除元素 H；保留 Bio.PDB 默认 altloc 选择及首个 model。HETATM 编码的修饰残基/配体不在该历史统计范围。该历史函数未单独排除元素 D，不能把它写成已验证的 non-H/non-D 实现；这是继承历史计算时的明确技术边界。

- 总界面残基 = 两条链各自参与接触的唯一残基数之和。
- 接口级 interface fraction = 总界面残基 / (两条链有坐标的残基数之和)。这是长度加权比例，不是两条链比例的算术均值；旧链级比例仍保留在 FoldBench 表中。
- 主分布和 ECDF 仅使用 contact_count > 0；散点、误差和配对相关性保留所有有效映射零值。未匹配 AF3 没有伪装为零。
- RAW 和 QC 分开保留。QC 沿用最小链间非氢原子距离 <0.5 Å 的严重重叠规则；零值、大接触数、低预测置信度及几何偏差都不自动排除。其他数据集中当前明确规则没有触发时，raw/QC 表保留，但不复制相同的 QC 图。
- Pedro 主单位为结构内链对观测，包含重复任务/重跑；另提供 `pedro_pair_level_raw.csv` 与 `pedro_pair_level_qc.csv`，先对同 job/蛋白对平均，再对不同 job 平均。不能将170644条结构观测当作独立蛋白对。

## 1. SAFES：配对基准，预期界面待确认

本地40个PDB（20实验+20 AF3）仅含坐标；本地 EV1/EV2 属于 Pedro，AF-assessment 压缩包中未发现 SAFES 接口定义。没有足够依据把多链复合物所有链对当作 SAFES 的目标结合伙伴。未使用同链 ID 直接配对，也未为增加数量枚举多链组合。

暂定子集仅使用所提供实验结构恰有两条蛋白链的案例，重新验证序列匹配并作一对一映射。它是结构层面的唯一链对，仍不能替代 SAFES 官方目标定义；多链案例逐条保留 pending 原因。见 `safes_interface_selection_audit.csv` 和 `safes_chain_mapping.csv`。
暂定配对数：10；图表均使用 `safes_provisional_two_chain` 名称。完整 SAFES 仍需要每个案例的 receptor/ligand 链组，尤其是否应合并多链结合伙伴。

## 2. FoldBench：配对基准

审计结果：279实验接口 / 239 assembly；132可映射接口 /108 assembly，147缺失接口 /131 assembly。147全部因为提供的AF3归档中缺少目标，同PDB其他assembly或命名变体也未命中。没有找到可修复后恢复这些案例的代码错误，没有重跑已有接触计算。详见 `FOLDBENCH_MATCHING_AUDIT.md` 和逐接口 CSV。
主配对集合为132，其中23有几何并列映射，14会影响至少一项指标；另给118个映射稳定接口的敏感性结果。映射不确定性不是结构质量错误，未静默删除。AF3 模型仍沿用历史最大 ranking_score 选择，未按接触数/参考一致性挑选。

## 3. Pedro vs 实验参考：仅分布级比较

没有建立逐一生物学对应关系，因此没有 Pedro-reference 配对散点或配对相关性。Pedro 复用70份既有8 Å残基计数文件，不重新解析94 GB结构归档。

参考原始范围29条候选链对完整保留；逐条审计其assembly来源、源物种和范围理由，见 `reference_interface_selection_audit.csv`。其中23条两条链均来自M. genitalium，4条仅一条来自该物种，2条两条均不来自该物种。主探索性图采用至少一条Mgen链的27条候选；全部29条另有敏感性图。这是数据集范围选择，不是QC，2条被排除主物种范围的记录未删除。
参考候选来自已提交的biological assemblies，但assembly成员身份本身不证明任意链对都是独立、预期生物学界面。表中将此项标为未逐对独立确认，不能将候选集合称为已验证的独立阳性基准。存在混合物种、重复条目/构象及同一assembly内相关链对，分布比较是探索性的。

## 零值与样本量

| dataset | n_computable | n_structures | n_zero | pct_zero | n_nonzero | n_qc_flagged |
| --- | --- | --- | --- | --- | --- | --- |
| FoldBench experimental all | 279 | 239 | 0 | 0.0 | 279 | 0 |
| FoldBench experimental paired | 132 | 108 | 0 | 0.0 | 132 | 0 |
| FoldBench AF3 paired | 132 | 108 | 1 | 0.75758 | 131 | 0 |
| Pedro raw | 170644 | 2156 | 116535 | 68.2913 | 54109 | 464 |
| Pedro QC | 170180 | 2156 | 116535 | 68.47749 | 53645 | 0 |
| Reference all 29 | 29 | 10 | 2 | 6.89655 | 27 | 0 |
| Reference Mgen-containing | 27 | 10 | 2 | 7.40741 | 25 | 0 |
| SAFES provisional experimental | 10 | 10 | 0 | 0.0 | 10 | 0 |
| SAFES provisional AF3 | 10 | 10 | 0 | 0.0 | 10 | 0 |

## 非零接口的中位数与 IQR

| dataset | metric | n | median | q1 | q3 |
| --- | --- | --- | --- | --- | --- |
| FoldBench experimental all | contact_count | 279 | 238.0 | 126.5 | 359.0 |
| FoldBench experimental all | total_interface_residues | 279 | 88.0 | 53.5 | 124.0 |
| FoldBench experimental all | interface_fraction | 279 | 0.32 | 0.20704 | 0.55903 |
| FoldBench experimental paired | contact_count | 132 | 231.5 | 128.0 | 299.0 |
| FoldBench experimental paired | total_interface_residues | 132 | 87.0 | 57.0 | 111.25 |
| FoldBench experimental paired | interface_fraction | 132 | 0.30935 | 0.1926 | 0.51136 |
| FoldBench AF3 paired | contact_count | 131 | 254.0 | 178.5 | 365.0 |
| FoldBench AF3 paired | total_interface_residues | 131 | 94.0 | 68.5 | 126.0 |
| FoldBench AF3 paired | interface_fraction | 131 | 0.31461 | 0.21014 | 0.56441 |
| Pedro raw | contact_count | 54109 | 69.0 | 16.0 | 165.0 |
| Pedro raw | total_interface_residues | 54109 | 41.0 | 15.0 | 81.0 |
| Pedro raw | interface_fraction | 54109 | 0.05556 | 0.02012 | 0.11024 |
| Pedro QC | contact_count | 53645 | 68.0 | 16.0 | 163.0 |
| Pedro QC | total_interface_residues | 53645 | 41.0 | 15.0 | 80.0 |
| Pedro QC | interface_fraction | 53645 | 0.05493 | 0.01988 | 0.10945 |
| Reference all 29 | contact_count | 27 | 494.0 | 93.5 | 759.0 |
| Reference all 29 | total_interface_residues | 27 | 193.0 | 59.5 | 292.0 |
| Reference all 29 | interface_fraction | 27 | 0.13625 | 0.02536 | 0.15834 |
| Reference Mgen-containing | contact_count | 25 | 480.0 | 91.0 | 760.0 |
| Reference Mgen-containing | total_interface_residues | 25 | 249.0 | 57.0 | 296.0 |
| Reference Mgen-containing | interface_fraction | 25 | 0.12505 | 0.02318 | 0.14906 |
| SAFES provisional experimental | contact_count | 10 | 230.0 | 201.0 | 247.5 |
| SAFES provisional experimental | total_interface_residues | 10 | 81.5 | 77.0 | 87.5 |
| SAFES provisional experimental | interface_fraction | 10 | 0.29373 | 0.27132 | 0.30582 |
| SAFES provisional AF3 | contact_count | 10 | 225.0 | 200.75 | 272.75 |
| SAFES provisional AF3 | total_interface_residues | 10 | 78.0 | 75.25 | 92.5 |
| SAFES provisional AF3 | interface_fraction | 10 | 0.28113 | 0.26896 | 0.31261 |

## 配对统计（保留有效零值）

### foldbench

| comparison | metric | n | pearson | spearman | median_difference | mean_difference | q1_difference | q3_difference | std_difference | mae | rmse | n_paired_with_any_zero |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| foldbench | contact_count | 132 | 0.66268 | 0.79617 | 17.0 | 40.12879 | -3.75 | 61.0 | 177.79031 | 75.97727 | 181.60465 | 1 |
| foldbench | total_interface_residues | 132 | 0.72431 | 0.85606 | 5.5 | 11.92424 | 0.0 | 16.0 | 46.65344 | 20.84848 | 47.98169 | 1 |
| foldbench | interface_fraction | 132 | 0.92975 | 0.89621 | 0.0 | 0.01421 | -0.01062 | 0.0189 | 0.11478 | 0.04961 | 0.11523 | 1 |

### foldbench_mapping_stable

| comparison | metric | n | pearson | spearman | median_difference | mean_difference | q1_difference | q3_difference | std_difference | mae | rmse | n_paired_with_any_zero |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| foldbench_mapping_stable | contact_count | 118 | 0.65556 | 0.79515 | 15.0 | 35.15254 | -2.5 | 60.5 | 176.1295 | 73.59322 | 178.8698 | 1 |
| foldbench_mapping_stable | total_interface_residues | 118 | 0.72242 | 0.8655 | 5.0 | 10.48305 | 0.0 | 15.0 | 45.81862 | 20.12712 | 46.81292 | 1 |
| foldbench_mapping_stable | interface_fraction | 118 | 0.937 | 0.91129 | 0.0 | 0.01113 | -0.01091 | 0.01689 | 0.11327 | 0.04856 | 0.11334 | 1 |

### safes_provisional_two_chain

| comparison | metric | n | pearson | spearman | median_difference | mean_difference | q1_difference | q3_difference | std_difference | mae | rmse | n_paired_with_any_zero |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| safes_provisional_two_chain | contact_count | 10 | 0.94127 | 0.97576 | 2.5 | 10.9 | -2.5 | 10.75 | 25.94203 | 15.9 | 26.91654 | 0 |
| safes_provisional_two_chain | total_interface_residues | 10 | 0.96311 | 0.96951 | -0.5 | 2.3 | -1.75 | 2.5 | 7.74668 | 4.5 | 7.70065 | 0 |
| safes_provisional_two_chain | interface_fraction | 10 | 0.99012 | 0.89091 | -0.00637 | -0.00025 | -0.01665 | 0.00472 | 0.03443 | 0.02176 | 0.03267 | 0 |

## Pedro 严重重叠 QC 对非零分布的影响

| metric | n_raw | n_qc | n_removed | median_raw | median_qc | mean_raw | mean_qc | max_raw | max_qc | ecdf_max_gap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| contact_count | 54109 | 53645 | 464 | 69.0 | 68.0 | 123.30045 | 120.8681 | 18018.0 | 4315.0 | 0.00442 |
| total_interface_residues | 54109 | 53645 | 464 | 41.0 | 41.0 | 59.20568 | 58.48142 | 1645.0 | 1420.0 | 0.00411 |
| interface_fraction | 54109 | 53645 | 464 | 0.05556 | 0.05493 | 0.07981 | 0.07921 | 0.94952 | 0.83505 | 0.00366 |

Pedro 非零残基接触中位数69→68，总界面残基中位数41→41，三个指标ECDF最大纵向差均小于0.0045。因此严重重叠过滤对主体分布影响较小；接触数最大值18018→4315，说明极端尾部明显受影响。这是对实际效应量的描述，不是新增过滤阈值。FoldBench、参考候选和SAFES暂定子集均无当前规则的排除项，raw与QC结果完全相同。

## 文件与图

- `metric_summary.csv`：含零/非零两套三指标均值、中位数、标准差、Q1/Q3、最小最大值。
- `zero_contact_summary.csv`：零值数量/比例，与QC分开。
- `foldbench_matching_audit.csv`：全部279接口的可用性、序列证据、链映射与失败原因。
- `foldbench_paired_metrics.csv` / `foldbench_paired_statistics.csv`：配对三指标、AF3−experimental差、Pearson和Spearman。
- `foldbench_raw.csv` / `foldbench_qc.csv`：所有可计算实验/AF3行；只有 `paired_success` 集合进入配对图。
- `pedro_raw.csv` / `pedro_qc.csv`、`pedro_qc_effect.csv`：逐观测及QC影响。
- `reference_all_29_candidates.csv`、`reference_candidates_raw.csv` / `reference_candidates_qc.csv`：全部候选原始/QC数据。
- `safes_provisional_two_chain_*`：暂定子集表；完整SAFES未完成。
- `figures/foldbench/`、`figures/foldbench_mapping_stable/`：三指标实验/AF3重叠归一化直方图、ECDF、保留零值散点图。
- `figures/safes_provisional_two_chain/`：暂定SAFES子集的相同图。
- `figures/pedro_reference/`：三指标重叠分布与ECDF，仅非零；27条主候选和29条范围敏感性分别命名。
- `figures/qc/`：Pedro raw/QC 非零三指标分布与ECDF。
- 所有图保留完整范围，计数指标采用log10直方图/log轴ECDF，比例使用0–1内原值；没有裁掉大值。PNG及PDF均保留。
- `validation.json`：验证及明确的SAFES未完成状态。旧原子对/残基对结果均未覆盖。

## 复现

依次运行 `src/audit_foldbench_residue8A.py`，检查审计结论，再运行 `src/compare_residue8A.py` 和 `src/summarize_residue8A.py`。新输出仅写入本目录；复现会更新本轮生成文件，不改历史结果目录。
