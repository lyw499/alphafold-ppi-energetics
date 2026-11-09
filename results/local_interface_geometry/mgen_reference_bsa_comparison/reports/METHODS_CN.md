# 预设方法与执行记录

主组为既有审计中Mgen-Mgen，不混入其他物种。参考BSA读取既有assembly文件，复用pedro_full_bsa/scripts/engine.py及其冻结PROTOCOL.json；结构阶段不打开Pedro归档。接触数量复用all_reference_contact_geometry.parquet并核对旧计数与fraction。原结果保持不变。

主模型：NB2，以log1p(BSA)线性项和dataset交互，4个均值参数加一个共同离散参数，替代旧4结点样条，避免小参考组过拟合。正BSA、raw Pedro、共同面积范围；标准化到参考面积经验分布。实验按PDB聚类（包含嵌套assembly），Pedro沿用已有prediction_cluster。cluster sandwich、delta method、t(df=min两组簇数-1)仅近似探索性区间和p；少簇风险须披露。不以不显著证明等效。

敏感性：原QC、非零接触、含零/全支持、PDB等权拟合与标准化、仅标准化等权、共同斜率、assembly聚类、指定低分辨率条目单独/联合排除、leave-one-PDB-out。额外用相同简化模型比较已有FoldBench/SAFES表，仅作描述，不替换旧结果。

不运行正式三组效应对比：共享Pedro、不同参考支持和目标使独立区间比较不充分。保存失败拟合与诊断，不新增为显著性服务的过滤。

小规模模型验证：单位权重与现有NB框架在无样条（log1p线性）情形的估计及CI一致。参考阶段三个链对与直接Bio.PDB SASA交叉核验。所有结构计算和全量拟合经Slurm后台执行。
