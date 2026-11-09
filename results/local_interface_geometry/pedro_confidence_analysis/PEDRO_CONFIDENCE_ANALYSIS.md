# Pedro 置信度与局部接触几何

仅使用已保存的界面汇总和既有置信度/QC 元数据；未读取 CIF/ZIP、未重跑任何基准或同种对照。

## 连接与定义
{
  "status": "passed",
  "join_keys": [
    "source_file",
    "job_id",
    "protein_pair",
    "best_model"
  ],
  "interfaces_before": 170644,
  "confidence_matched": 170644,
  "confidence_unmatched": 0,
  "summary_duplicate_keys": 0,
  "confidence_duplicate_keys": 0,
  "qc_matched": 170644,
  "qc_unmatched": 0,
  "interface_id_duplicates": 0,
  "nonzero": 54109,
  "zero": 116535,
  "severe": 464,
  "confidence_values_exactly_reconciled": true,
  "qc_flags_reconciled": true,
  "primary_metric": "best_model_pair_iptm",
  "available_metrics": [
    "best_model_pair_iptm",
    "mean_pair_iptm",
    "job_ranking_score"
  ],
  "ptm_available": false,
  "selected_model_counts": {
    "0": 170644
  }
}

主置信度 best_model_pair_iptm 是所选模型的链对 ipTM；mean_pair_iptm 是原流程跨模型平均链对 ipTM；job_ranking_score 是整个预测复合物的分数，同一作业多个界面共享，不能视为独立界面置信度。现有输入表未包含 pTM，未从结构或归档提取新字段。

主结局为 3 < d_min ≤4 Å 的残基对比例，直接复用 count_3_4 / n_contacts。零界面不参与条件几何分析；单独用于接触形成概率的描述。分位组按非零界面定义，零接触分析重新使用全部界面分位组；两者边界不可混用。使用 qcut，保持相同分数不跨组，重复边界合并，实际边界和样本量已保存。高置信度表示相对分位数，不是外部校准的可信 PPI 阈值。

几何分层的 median_d_min 是“界面内距离中位数的跨界面中位数”，不是池化距离中位数。池化接触比例从已有精确分箱计数加总，无需再次读取接触明细。Pearson/Spearman 均为描述性关联，不报告 p 值。界面仍嵌套于复合物，重复蛋白对和作业不独立，不推断因果。

## 置信度字段

```
              metric  n_present  n_missing  n_unique     min  max  median
best_model_pair_iptm     170644          0        80   0.010 0.89    0.08
      mean_pair_iptm     170644          0       439   0.016 0.89    0.08
   job_ranking_score     170644          0        25 -99.660 0.56    0.33
```

## 连续关联

```
            version                    x             y  n_interfaces  spearman   pearson
                raw best_model_pair_iptm  fraction_3_4         54109  0.160189  0.151444
                raw best_model_pair_iptm  median_d_min         54109 -0.061037 -0.066800
                raw best_model_pair_iptm fraction_lt_2         54109  0.147342  0.001917
                raw       mean_pair_iptm  fraction_3_4         54109  0.129080  0.133554
                raw       mean_pair_iptm  median_d_min         54109 -0.039835 -0.053474
                raw       mean_pair_iptm fraction_lt_2         54109  0.121801 -0.005606
                raw    job_ranking_score  fraction_3_4         54109  0.031496  0.004271
                raw    job_ranking_score  median_d_min         54109 -0.018278  0.000174
                raw    job_ranking_score fraction_lt_2         54109 -0.005833 -0.021413
existing_QC_passing best_model_pair_iptm  fraction_3_4         53645  0.160452  0.151588
existing_QC_passing best_model_pair_iptm  median_d_min         53645 -0.060669 -0.066831
existing_QC_passing best_model_pair_iptm fraction_lt_2         53645  0.146021  0.001202
existing_QC_passing       mean_pair_iptm  fraction_3_4         53645  0.129268  0.133684
existing_QC_passing       mean_pair_iptm  median_d_min         53645 -0.039361 -0.053487
existing_QC_passing       mean_pair_iptm fraction_lt_2         53645  0.120369 -0.006595
existing_QC_passing    job_ranking_score  fraction_3_4         53645  0.031361  0.004899
existing_QC_passing    job_ranking_score  median_d_min         53645 -0.018318 -0.001837
existing_QC_passing    job_ranking_score fraction_lt_2         53645 -0.006391 -0.013977
     severe_overlap best_model_pair_iptm  fraction_3_4           464 -0.063914 -0.046834
     severe_overlap best_model_pair_iptm  median_d_min           464  0.223821  0.201687
     severe_overlap best_model_pair_iptm fraction_lt_2           464 -0.244754 -0.211844
     severe_overlap       mean_pair_iptm  fraction_3_4           464 -0.068384 -0.046021
     severe_overlap       mean_pair_iptm  median_d_min           464  0.215653  0.187351
     severe_overlap       mean_pair_iptm fraction_lt_2           464 -0.226570 -0.195101
     severe_overlap    job_ranking_score  fraction_3_4           464  0.028777 -0.048006
     severe_overlap    job_ranking_score  median_d_min           464  0.025234  0.170280
     severe_overlap    job_ranking_score fraction_lt_2           464 -0.059149 -0.173553
```

## 置信度十分位

```
 decile  n_interfaces  n_nonzero  n_zero  n_contacts  confidence_min  confidence_max  confidence_median  median_contact_count  median_interface_fraction  pooled_close_fraction  interface_close_mean  interface_close_median  interface_close_q25  interface_close_q75  median_interface_median_d_min  pooled_short_fraction  interface_short_mean  interface_short_median
      1          7409       7409       0      395285            0.02            0.06               0.06                  27.0                   0.055882               0.080165              0.056427                0.056180             0.000000             0.090909                       6.588784               0.004751              0.002679                     0.0
      2          5634       5634       0      413268            0.07            0.07               0.07                  41.0                   0.053743               0.081402              0.059668                0.062500             0.000000             0.090909                       6.572349               0.005195              0.003103                     0.0
      3          6461       6461       0      546418            0.08            0.08               0.08                  48.0                   0.051370               0.081284              0.059679                0.063830             0.000000             0.090909                       6.579246               0.005053              0.003182                     0.0
      4          6440       6440       0      615025            0.09            0.09               0.09                  56.0                   0.048813               0.080228              0.060375                0.064297             0.000000             0.089058                       6.575258               0.005058              0.003332                     0.0
      5          5817       5817       0      622727            0.10            0.10               0.10                  66.0                   0.049199               0.080612              0.060296                0.066265             0.014493             0.089286                       6.567496               0.005367              0.003382                     0.0
      6          4896       4896       0      613159            0.11            0.11               0.11                  77.5                   0.049468               0.080575              0.062511                0.068648             0.030303             0.089552                       6.566099               0.005364              0.003549                     0.0
      7          3889       3889       0      535738            0.12            0.12               0.12                  86.0                   0.049195               0.081607              0.064255                0.069767             0.037037             0.090517                       6.562778               0.005066              0.003528                     0.0
      8          3013       3013       0      452661            0.13            0.13               0.13                 106.0                   0.053452               0.081529              0.066513                0.071429             0.045455             0.090452                       6.544608               0.004794              0.003456                     0.0
      9          5802       5802       0     1117723            0.14            0.16               0.15                 137.0                   0.057908               0.082517              0.069618                0.075000             0.051039             0.092551                       6.535868               0.004899              0.003712                     0.0
     10          4748       4748       0     1359660            0.17            0.89               0.19                 207.0                   0.090352               0.089039              0.083430                0.084424             0.068182             0.100354                       6.527106               0.004925              0.003158                     0.0
```

## 三分位

```
stratum  n_interfaces  n_nonzero  n_zero  n_contacts  confidence_min  confidence_max  confidence_median  median_contact_count  median_interface_fraction  pooled_close_fraction  interface_close_mean  interface_close_median  interface_close_q25  interface_close_q75  median_interface_median_d_min  pooled_short_fraction  interface_short_mean  interface_short_median
    low         19504      19504       0     1354971            0.02            0.08               0.07                  37.0                   0.053630               0.080994              0.058440                0.061224             0.000000             0.090909                       6.579612               0.005008              0.002968                     0.0
 middle         17153      17153       0     1850911            0.09            0.11               0.10                  65.0                   0.049118               0.080472              0.060958                0.066667             0.016949             0.089286                       6.569645               0.005263              0.003411                     0.0
   high         17452      17452       0     3465782            0.12            0.89               0.14                 136.0                   0.063172               0.084806              0.071645                0.076923             0.053254             0.094241                       6.539066               0.004921              0.003476                     0.0
```

## 大小关联

```
                       x            y  n_interfaces  spearman  pearson
              n_contacts fraction_3_4         54109  0.556201 0.250303
total_interface_residues fraction_3_4         54109  0.532271 0.292758
      interface_fraction fraction_3_4         54109  0.534044 0.321068
      total_chain_length fraction_3_4         54109 -0.005585 0.006106
```

## 同尺寸分层

```
            version  size_quartile  n_interfaces  n_nonzero  n_zero  n_contacts  confidence_min  confidence_max  confidence_median  median_contact_count  median_interface_fraction  pooled_close_fraction  interface_close_mean  interface_close_median  interface_close_q25  interface_close_q75  median_interface_median_d_min  pooled_short_fraction  interface_short_mean  interface_short_median                    x            y  spearman   pearson  min_contacts  max_contacts
                raw              1         13600      13600       0       78889            0.02            0.38               0.08                   4.0                   0.008065               0.038193              0.026432                0.000000             0.000000             0.000000                       6.927774               0.000735              0.000362                0.000000 best_model_pair_iptm fraction_3_4 -0.024594 -0.019886             1            16
                raw              2         13523      13523       0      543380            0.02            0.54               0.09                  39.0                   0.038208               0.067345              0.065130                0.060606             0.037037             0.090909                       6.609411               0.003316              0.003020                0.000000 best_model_pair_iptm fraction_3_4 -0.028185  0.014194            17            69
                raw              3         13465      13465       0     1513821            0.03            0.79               0.10                 110.0                   0.079661               0.078246              0.077729                0.076923             0.059524             0.094737                       6.527712               0.004738              0.004668                0.000000 best_model_pair_iptm fraction_3_4 -0.013787  0.085385            70           165
                raw              4         13521      13521       0     4535574            0.03            0.89               0.12                 263.0                   0.149079               0.086991              0.084973                0.084656             0.072100             0.096923                       6.493105               0.005413              0.005062                0.003390 best_model_pair_iptm fraction_3_4  0.013827  0.155845           166         18018
existing_QC_passing              1         13598      13598       0       78859            0.02            0.38               0.08                   4.0                   0.008065               0.038182              0.026426                0.000000             0.000000             0.000000                       6.927774               0.000710              0.000352                0.000000 best_model_pair_iptm fraction_3_4 -0.024787 -0.019973             1            16
existing_QC_passing              2         13503      13503       0      542363            0.02            0.54               0.09                  39.0                   0.038186               0.067324              0.065109                0.060606             0.037037             0.090909                       6.610012               0.003254              0.002966                0.000000 best_model_pair_iptm fraction_3_4 -0.027373  0.014829            17            69
existing_QC_passing              3         13358      13358       0     1500556            0.03            0.79               0.10                 110.0                   0.079661               0.078231              0.077717                0.076923             0.059524             0.094737                       6.528872               0.004552              0.004497                0.000000 best_model_pair_iptm fraction_3_4 -0.013907  0.085674            70           165
existing_QC_passing              4         13186      13186       0     4362191            0.03            0.89               0.12                 262.0                   0.149127               0.086956              0.085011                0.084728             0.072109             0.097004                       6.495289               0.004463              0.004624                0.003247 best_model_pair_iptm fraction_3_4  0.017087  0.159408           166          4315
```

## QC 对比

```
           qc_group  n_interfaces  n_nonzero  n_zero  n_contacts  confidence_min  confidence_max  confidence_median  median_contact_count  median_interface_fraction  pooled_close_fraction  interface_close_mean  interface_close_median  interface_close_q25  interface_close_q75  median_interface_median_d_min  pooled_short_fraction  interface_short_mean  interface_short_median
                all        170644      54109  116535     6671664            0.01            0.89               0.08                   0.0                   0.000000               0.082829              0.063497                0.068966             0.026316             0.091463                       6.559880               0.005034              0.003272                 0.00000
existing_QC_passing        170180      53645  116535     6483969            0.01            0.89               0.08                   0.0                   0.000000               0.082702              0.063335                0.068966             0.025000             0.091398                       6.561726               0.004337              0.003092                 0.00000
     severe_overlap           464        464       0      187695            0.04            0.30               0.11                 243.0                   0.125184               0.087237              0.082248                0.081725             0.069672             0.095415                       6.371963               0.029111              0.024094                 0.01889
```

## 零与非零置信度

```
              metric   group  n_interfaces     mean  median   q25   q75     min   max
best_model_pair_iptm    zero        116535 0.075847   0.070 0.060 0.090   0.010 0.440
best_model_pair_iptm nonzero         54109 0.107580   0.100 0.080 0.130   0.020 0.890
      mean_pair_iptm    zero        116535 0.076734   0.072 0.060 0.090   0.016 0.462
      mean_pair_iptm nonzero         54109 0.104919   0.094 0.074 0.122   0.016 0.890
   job_ranking_score    zero        116535 0.282246   0.330 0.320 0.340 -99.660 0.560
   job_ranking_score nonzero         54109 0.280227   0.330 0.320 0.340 -99.660 0.560
```

## 零概率十分位

```
 decile  n_interfaces  n_zero  zero_fraction  confidence_min  confidence_max  confidence_median
      1         28612   25223       0.881553            0.01            0.05               0.05
      2         22336   18316       0.820021            0.06            0.06               0.06
      3         25007   19373       0.774703            0.07            0.07               0.07
      4         23235   16774       0.721928            0.08            0.08               0.08
      5         19371   12931       0.667544            0.09            0.09               0.09
      6         14865    9048       0.608678            0.10            0.10               0.10
      7         10673    5777       0.541272            0.11            0.11               0.11
      8         12574    5672       0.451090            0.12            0.13               0.12
      9         13971    3421       0.244864            0.14            0.89               0.16
```

## 实验上下文

```
                      group  n_interfaces  pooled_close_fraction  interface_close_median  interface_close_q25  interface_close_q75  n_nonzero  n_zero  n_contacts  confidence_min  confidence_max  confidence_median  median_contact_count  median_interface_fraction  interface_close_mean  median_interface_median_d_min  pooled_short_fraction  interface_short_mean  interface_short_median stratum
     FoldBench Experimental           132               0.143617                0.131010             0.118228             0.158104        NaN     NaN         NaN             NaN             NaN                NaN                   NaN                        NaN                   NaN                            NaN                    NaN                   NaN                     NaN     NaN
         SAFES Experimental            50               0.136180                0.128346             0.115320             0.154315        NaN     NaN         NaN             NaN             NaN                NaN                   NaN                        NaN                   NaN                            NaN                    NaN                   NaN                     NaN     NaN
          Mgen Experimental            23               0.107523                0.105882             0.087391             0.132895        NaN     NaN         NaN             NaN             NaN                NaN                   NaN                        NaN                   NaN                            NaN                    NaN                   NaN                     NaN     NaN
     Mgen exclude_8PBZ_8PC1            16               0.121091                0.131926             0.108058             0.135826        NaN     NaN         NaN             NaN             NaN                NaN                   NaN                        NaN                   NaN                            NaN                    NaN                   NaN                     NaN     NaN
Mgen exclude_8PBZ_8PC0_8PC1            15               0.124546                0.131926             0.114583             0.136187        NaN     NaN         NaN             NaN             NaN                NaN                   NaN                        NaN                   NaN                            NaN                    NaN                   NaN                     NaN     NaN
      Pedro overall nonzero         54109               0.082829                0.068966             0.026316             0.091463    54109.0     0.0   6671664.0            0.02            0.89               0.10                  69.0                   0.055556              0.063497                       6.559880               0.005034              0.003272                     0.0     NaN
          Pedro low tercile         19504               0.080994                0.061224             0.000000             0.090909    19504.0     0.0   1354971.0            0.02            0.08               0.07                  37.0                   0.053630              0.058440                       6.579612               0.005008              0.002968                     0.0     low
       Pedro middle tercile         17153               0.080472                0.066667             0.016949             0.089286    17153.0     0.0   1850911.0            0.09            0.11               0.10                  65.0                   0.049118              0.060958                       6.569645               0.005263              0.003411                     0.0  middle
         Pedro high tercile         17452               0.084806                0.076923             0.053254             0.094241    17452.0     0.0   3465782.0            0.12            0.89               0.14                 136.0                   0.063172              0.071645                       6.539066               0.004921              0.003476                     0.0    high
           Pedro top decile          4748               0.089039                0.084424             0.068182             0.100354     4748.0     0.0   1359660.0            0.17            0.89               0.19                 207.0                   0.090352              0.083430                       6.527106               0.004925              0.003158                     0.0     NaN
```


## 八个问题：结果解释


1. **置信度与近接触比例的关联**：在 54,109 个非零界面中，所选模型链对 ipTM 与 (3,4] Å 比例的 Spearman ρ=0.1602，Pearson r=0.1514。这是界面层面的关联，不能解释为置信度造成几何变化。

2. **效应大小**：最低可用分位组 ipTM 范围 0.02–0.06，最高组 0.17–0.89。近接触比例逐界面中位数由 5.62% 到 8.44%（变化 +2.82 个百分点）；池化比例由 8.02% 到 8.90%。应结合全部分位组判断是否单调。

3. **高置信度与实验比较**：最高分位组池化比例 8.90%，逐界面中位数 8.44%。实验池化参考：FoldBench 14.36%、SAFES 13.62%、Mgen 原始 10.75%、Mgen 敏感性 12.11–12.45%。最高组是相对分位数，范围较宽，不能称作全部经校准的高可信结构。

4. **控制界面大小后的描述性关联**：相同接触数四分位层内的 Spearman ρ 为 Q1 (1–16 contacts): -0.0246, Q2 (17–69 contacts): -0.0282, Q3 (70–165 contacts): -0.0138, Q4 (166–18018 contacts): 0.0138。尺寸层仍有剩余跨度；这是分层检查，不是因果调整。QC 通过子集的同一分层结果也已保存。

5. **是否由严重重叠主导**：原始与既有 QC 通过数据的池化近接触比例分别为 8.283% 和 8.270%；严重重叠组为 8.724%。QC 前后总体差异很小，不支持严重重叠是整体近接触缺口的主要来源。高 (3,4] Å 比例没有被定义为 QC 失败。

6. **极短距离与置信度**：ipTM 与 <2 Å 接触比例的 Spearman ρ=0.1473，Pearson r=0.0019。需结合尺寸、QC 分层以及大量零短接触比例解释；没有自动将 <2 Å 定义为碰撞。

7. **零接触与置信度**：按全体界面定义的最低/最高置信度分位组，零接触率分别为 88.16%/24.49%。这些分位边界与非零几何分析不同。零接触的条件几何始终为 NA，没有置零参与距离分析。

8. **解释边界**：最终模式判断应同时参考关联方向、幅度、尺寸分层和高分位相对实验的缺口。尺寸与置信度可能共同关联于接触形成及局部几何；本分析没有识别因果机制。多个界面来自同一预测作业，且存在重复蛋白对；不把接触数量或界面数量视为独立重复的实验样本数。


## 模型选择及评分限制
全部既有记录 selected_model=0；仅凭这些结果无法比较不同所选模型索引或判断换模型是否改变几何。ranking score 有负值，沿用原始数据而不截断或重定义；它属于复合物级指标，不作为主链对置信度。mean_pair_iptm 使用跨模型平均，而本几何只属于所选模型，其相关性只能作为辅助稳健性检查。


## QC 条件置信度补充
confidence_median 对 QC 表中所有界面计算，新增 nonzero_confidence_median 则仅针对几何有定义的非零界面，避免两种分母混淆。

```
           qc_group  n_interfaces  n_nonzero  n_zero  n_contacts  confidence_min  confidence_max  confidence_median  median_contact_count  median_interface_fraction  pooled_close_fraction  interface_close_mean  interface_close_median  interface_close_q25  interface_close_q75  median_interface_median_d_min  pooled_short_fraction  interface_short_mean  interface_short_median  nonzero_confidence_median
                all        170644      54109  116535     6671664            0.01            0.89               0.08                   0.0                   0.000000               0.082829              0.063497                0.068966             0.026316             0.091463                       6.559880               0.005034              0.003272                 0.00000                       0.10
existing_QC_passing        170180      53645  116535     6483969            0.01            0.89               0.08                   0.0                   0.000000               0.082702              0.063335                0.068966             0.025000             0.091398                       6.561726               0.004337              0.003092                 0.00000                       0.10
     severe_overlap           464        464       0      187695            0.04            0.30               0.11                 243.0                   0.125184               0.087237              0.082248                0.081725             0.069672             0.095415                       6.371963               0.029111              0.024094                 0.01889                       0.11
```

## 模式判断：更支持 B，但需区分整体趋势和极高分尾部

整体是较弱的正关联（ρ=0.160），并非 C 所指的完全无关联；最高十分位池化比例仅 8.90%，仍明显低于实验 12–14%，不支持将总体缺口主要归因为低置信度并选择 A。更符合 **B：置信度只解释部分现象，存在明显剩余差异**。这里的“解释”仅是描述，不是回归解释方差或因果解释。

接触数与近接触比例的 ρ=0.556，总界面残基数为 0.532，界面覆盖比例为 0.534，总链长约为 −0.006。按接触数四分位分层后，置信度关系减弱至 ρ≈−0.028 至 +0.014，说明总体关系与尺寸组成密切相关；不能断言尺寸造成缺口或已经消除所有混杂。较大尺寸层 Pearson 仍为正，也提示可能存在非线性高分尾部，不能说任何尺寸下都没有置信度效应。

### 对“最高置信度”的补充检查

最高十分位中位 ipTM 只有 0.19，故另检查上 5%、1%、0.5%、0.1% 的嵌套分位阈值。阈值保留同分数全部界面，实际样本比例略有偏离；这只是透明的分位敏感性，不是新增生物学置信度标准。上尾 >=90% 阈值包含边界同分数，和原 qcut 最后一组的右闭分箱不是完全同一个集合。

```
 quantile  minimum_score    n  actual_fraction  median_score  pooled_close  median_close      q25      q75  jobs  protein_pairs  median_contacts  severe
    0.900        0.16000 6239         0.115304          0.18      0.087875      0.082609 0.065382 0.098784  1923           5185            192.0      85
    0.950        0.19000 2788         0.051526          0.23      0.091794      0.088608 0.072840 0.105618  1484           2347            223.5      27
    0.990        0.32000  581         0.010738          0.42      0.106834      0.107143 0.090090 0.120301   479            460            238.0       0
    0.995        0.43000  276         0.005101          0.52      0.111249      0.111318 0.097410 0.124338   248            201            279.5       0
    0.999        0.61892   55         0.001016          0.68      0.116703      0.115745 0.106675 0.129451    52             39            445.0       0
```

最顶端约 0.1% 只有 55 个界面（52 个归档作业、39 个蛋白对），池化 11.67%、界面中位数 11.57%，确实更接近 Mgen 敏感性 12.11–12.45%，但仍低于一般实验 13.62–14.36%。这些界面更大（接触数中位数 445），且样本很小、选择特殊，所以不能用它们证明整个 Pedro 缺口只是低置信度预测驱动。结论应保留这一尾部趋近现象，而不是声称所有高分预测均远低于实验。

### 零接触与极短接触需分别理解

全部界面的分数存在大量并列值，qcut 原请求十组，重复边界合并后实际为 **9 组**；非零几何分析实际为 **10 组**，各组样本数并不相等。零接触率从最低组 88.16% 降到最高组 24.49%；零/非零 ipTM 中位数分别为 0.07/0.10。这是接触形成关联，不是条件几何证据。

<2 Å 比例对置信度的 Spearman 为 +0.147，但 Pearson 约 +0.002，许多界面该比例为零，分位池化短接触率没有明显单调下降；不能声称高 ipTM 自动消除了短距离疑点。QC 通过子集的相关性相近。既有严重重叠排除后近接触池化比例仅从 8.283% 变为 8.270%，说明此类 QC 失败不是主要缺口来源。
