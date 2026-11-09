# Mgen 同物种参考：BSA 控制后的 3–4 Å 接触数量

原始参考共 29 条链对；主组只纳入 Mgen–Mgen 23 条，其中 21 条非零接触，来自 7 个 PDB、8 个 biological assemblies。mixed-species 和非 Mgen 链对均保留在完整表，但不进入主模型。

## 数据、面积和零值
BSA=(SASA_A+SASA_B−SASA_AB)/2，Å²。直接导入已冻结 Pedro 引擎：Bio.PDB ShrakeRupley 对应的100点球面、probe 1.4 Å、相同原子半径和 H/D 排除，第一坐标model及默认altloc；只纳入 residue.id[0] 为空的实际坐标残基。不补缺失原子或残基，排除HETATM修饰残基/配体。按已验证的biological assembly与label chain选择，仅两条指定链计算遮挡。三个参考链对与直接Bio.PDB SASA计算交叉验证，完整记录见direct_sasa_validation.csv。

接触数量复用逐残基最小距离，严格3<d_min≤4 Å；不计独立原子对。零接触保留且close fraction为NA；BSA=0时密度为NA。Pedro原始170644行与既有severe-overlap QC分开，原结果不覆盖。

## 主模型与不确定性
参考样本很小，预先采用NB2 + log1p(BSA)线性项 + 数据集交互（两组可有不同斜率），不用原4结点样条。标准化到所纳入实验接口的实际面积分布；主分析为raw Pedro、正BSA、共同面积[min,max]范围。共同范围仅避免端点外推，并不保证范围内部密度充分。参考PDB是主聚类单位，已覆盖其嵌套assembly依赖；assembly聚类仅敏感性，不能当作更多独立实验。Pedro沿用2033个预测簇，仍有共享蛋白相关性。

95% CI是cluster sandwich + delta method + t(df=min两组簇数−1)的近似，不是已证明有限样本有效的区间。仅约7个参考PDB时，即使使用t修正仍可能低估不确定性；不把巨大Pedro样本视为弥补参考不足。未报告接触级独立检验或用少簇bootstrap制造精确性。p仅探索性，优先看效应、逐PDB影响与模型诊断。一个预设新主比较，未与旧两项Holm家族合并；敏感性p不作额外确证。

- primary：Pedro n=37,069，参考 n=21，7 PDB/8 assemblies；共同面积范围 101.83–3257.30 Å²，参考scope排除2、support排除0。调整比 0.748，近似95% CI 0.545–1.025，探索性p=0.0652。
- pedro_existing_qc：Pedro n=36,634，参考 n=21，7 PDB/8 assemblies；共同面积范围 101.83–3257.30 Å²，参考scope排除2、support排除0。拟合失败：Unstable Hessian；不填补推断。
- zero_inclusive_common：Pedro n=170,096，参考 n=23，7 PDB/8 assemblies；共同面积范围 0.00–3257.30 Å²，参考scope排除0、support排除0。调整比 0.772，近似95% CI 0.563–1.060，探索性p=0.0927。
- zero_inclusive_full_descriptive_model：Pedro n=170,644，参考 n=23，7 PDB/8 assemblies；共同面积范围 0.00–3257.30 Å²，参考scope排除0、support排除0。调整比 0.758，近似95% CI 0.552–1.042，探索性p=0.0773。
- pdb_equal：Pedro n=37,069，参考 n=21，7 PDB/8 assemblies；共同面积范围 101.83–3257.30 Å²，参考scope排除2、support排除0。调整比 0.767，近似95% CI 0.641–0.916，探索性p=0.0108。
- target_pdb_equal：Pedro n=37,069，参考 n=21，7 PDB/8 assemblies；共同面积范围 101.83–3257.30 Å²，参考scope排除2、support排除0。调整比 0.750，近似95% CI 0.550–1.023，探索性p=0.0641。
- common_area_slope：Pedro n=37,069，参考 n=21，7 PDB/8 assemblies；共同面积范围 101.83–3257.30 Å²，参考scope排除2、support排除0。调整比 0.747，近似95% CI 0.540–1.033，探索性p=0.0701。
- exclude_8PBZ_8PC1：Pedro n=22,784，参考 n=14，5 PDB/6 assemblies；共同面积范围 377.71–3216.92 Å²，参考scope排除2、support排除0。调整比 0.640，近似95% CI 0.515–0.796，探索性p=0.00475。
- exclude_8PBZ_8PC1_8PC0：Pedro n=22,784，参考 n=13，4 PDB/5 assemblies；共同面积范围 377.71–3216.92 Å²，参考scope排除2、support排除0。调整比 0.622，近似95% CI 0.508–0.761，探索性p=0.00492。

逐PDB删除得到的点估计范围 0.653–0.869。每次重新计算共同支持和标准化目标，因此变化包含组成/面积支持变化，不能只归因于该PDB的残差。

PDB等权敏感性：参考每个PDB总权重相同，归一化后参考总权重仍等于参考接口数；同时在拟合及标准化中应用。target_pdb_equal仅改变标准化权重。Pedro保持接口权重，未宣称整个两组都变为生物实体等权。

## 与原比较的关系
- FoldBench Experimental：0.548 [0.529, 0.567]；previous flexible-spline benchmark。
- SAFES Experimental：0.528 [0.496, 0.562]；previous flexible-spline benchmark。
- Mgen–Mgen Experimental：0.748 [0.545, 1.025]；new simple-model exploratory sensitivity。

三组结果只描述比较，参考分布、共同支持、标准化目标和主模型复杂度不同；另提供同一简化模型拟合FoldBench/SAFES的敏感性行，但仍非同一目标人群。未做正式效应差异检验，不能因一组显著另一组不显著宣称差异不同。三组共享Pedro，区间不能当作独立研究区间。

## 各图如何阅读
1. 图01：全部接口与正面积接口的BSA ECDF；检查样本选择和面积覆盖，零面积保留在左图。
2. 图02：每个接口的近接触数量与面积；参考按PDB着色，可看到同一复合物重复贡献与面积趋势。
3. 图03：左侧为正面积近接触密度，右侧为含零数量ECDF；密度不是主调整模型，也不证明严格面积比例关系。
4. 图04：主模型、QC、PDB等权和共同斜率的调整比及近似区间；1表示相同标准化均值，跨1不证明等效。
5. 图05：旧FoldBench/SAFES主结果与新Mgen探索性结果并列，不能据此直接推断物种的因果贡献。
6. 图06：逐个删去PDB后重拟合，查看结论是否依赖单个条目；并不是7个新独立重复实验。
7. 图07：左图面积分层的中位数/IQR和样本簇数；右图观察计数与拟合趋势，用于识别模型不充分和局部稀疏。

## 限制
Mgen参考是同物种敏感性集合，不是与FoldBench/SAFES等价的独立第三基准。少数PDB、同一蛋白体系、重复assembly/对称拷贝、低分辨率（尤其8PBZ/8PC0/8PC1）、缺失侧链和非随机选择限制推断。PDB聚类不能消除跨PDB同源体系依赖。BSA匹配无法恢复未观察的原子，也不等于功能性界面或亲和力匹配。低分辨率排除按既有指定条目进行，不新增任意阈值。

共 24 个模型/敏感性，1 个失败；完整失败原因和样本数见model_results.csv。
## 补充诊断与可报告结论

下面给出描述性数字，避免将模型结果代替原始分布。
- Pedro raw / all：n=170,644；零接触116,535（68.29%）；近接触数中位数0，IQR 0–0；BSA中位数0.00 Å²。
- Pedro raw / nonzero_contacts：n=54,109；零接触0（0.00%）；近接触数中位数5，IQR 1–13；BSA中位数291.48 Å²。
- Pedro QC / all：n=170,180；零接触116,535（68.48%）；近接触数中位数0，IQR 0–0；BSA中位数0.00 Å²。
- Pedro QC / nonzero_contacts：n=53,645；零接触0（0.00%）；近接触数中位数5，IQR 1–13；BSA中位数287.20 Å²。
- Mgen–Mgen Experimental / all：n=23；零接触2（8.70%）；近接触数中位数51，IQR 11–91.5；BSA中位数2606.54 Å²。
- Mgen–Mgen Experimental / nonzero_contacts：n=21；零接触0（0.00%）；近接触数中位数55，IQR 12–93；BSA中位数2623.73 Å²。

主模型参考组的平均Pearson残差平方/模型方差为 3.19（仅描述性诊断，未作自由度校正），明显高于1。共同NB离散参数主要受巨大Pedro组影响，参考异质性未被充分解释；不应将模型视为已验证的真实生成机制。cluster sandwich不能在只有7个PDB时保证区间覆盖，当前p值只可作为探索性近似。
共同支持充分性：21条正面积参考全部处于Pedro面积端点范围内，主分析保留37,069条Pedro。但参考面积分布不连续，四个分层只有2、3、3、4个参考PDB；2623.73–2754.12 Å²层仅134条Pedro，最高层351条。因此有端点重叠不等于局部密度和独立样本充分，不能把37,069理解为各面积区域都有强支持。
QC拟合的Unstable Hessian是诊断失败而非无差异；不使用该拟合的比值/CI/p。raw/QC原始计数中位数和密度描述很接近，但不能因此补写面积调整QC模型的结论。图04明确标记失败。

同一简化模型的附加比较：
- same_simple_model_FoldBench Experimental：调整比0.546，近似95% CI 0.527–0.566。
- same_simple_model_SAFES Experimental：调整比0.529，近似95% CI 0.492–0.567。

Mgen主点估计0.748描述性上比旧0.548/0.528更接近1，简化模型附加检查也保留这个方向；这不证明效应差异具有统计显著性，也不能把缩小的差距定量归因于物种。Mgen估计约低25%，但区间同时容许约45%更低及约3%更高，不足以证明两组相同。
排除8PBZ/8PC1/8PC0后估计约0.622，参考仅13个接口、4个PDB；PDB等权估计约0.767。这些变化表明参考组成和权重有影响，不能选取达到显著性的子集替代主结果。

最终交付包括全部29条参考BSA、23条同物种主参考表、24项模型结果（1项明确失败）、7组PNG/PDF图。未重算Pedro，输入SHA256复核保持不变。