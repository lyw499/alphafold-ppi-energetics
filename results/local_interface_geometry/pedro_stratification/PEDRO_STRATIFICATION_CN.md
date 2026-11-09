# Pedro AF3 界面形成、局部接触及大小分层

本报告只分析已保存的 Pedro 界面和类型汇总。8 Å 定义为跨链残基的最近非氢原子距离≤8 Å，每个残基对计一次；近接触为3<d_min≤4 Å。零接触界面完整保留在形成分析中，其近接触比例未定义。原始群体与去除既有 severe-overlap 标记的结果并列展示，后者不替代原始结果。

## 覆盖率、来源和层级

```text
   all  nonzero   zero  severe  unique_pairs  repeated_pairs  max_models_per_pair  jobs selected_models  both_annotation_coverage
170644    54109 116535     464        113050           43190                    8  2156   {'0': 170644}                    1.0000
```
所有选中模型均为model 0；这不表示蛋白对唯一。pair_key由原蛋白ID无向排序得到；一个蛋白对的多作业/多模型记录在pair_multiplicity.csv中。模型的不确定性按pair_key聚类，另作pair等权敏感性。一个蛋白可出现在多个pair中，因此pair聚类仍不能完全消除蛋白级依赖。

原始记录包含蛋白ID、链长、界面残基比例、pair ipTM、既有severe-overlap、源archive/CIF成员和原始蛋白注释。注释可作为候选名称来源，但有描述不等于功能经实验证实；本轮不按ID猜测功能。pTM、独立实验结合标签及可直接访问的解压CIF路径不可用。当前archive/CIF成员只是既有来源索引，代表列表供后续检查。

结构QC字段核查如下。has_clash是原预测附带标记，与本项目既有severe-overlap定义不同；本轮敏感性只按后者排除，其他标记保留并单独报告。parser warning为空不等于实验质量已验证。

## 变量与分母

形成率：所有有效预测中n_contacts>0的比例。局部packing：仅n_contacts>0，分子count_3_4，分母n_contacts。类型组成share_type=n_type/n_contacts；类型内部近接触within_type=count_3to4_type/n_type，类型缺失时为NA；加性贡献contribution_type=count_3to4_type/n_contacts，三类贡献相加等于总近接触比例。BB只包括N、CA、C、O；OXT按既有主分类落在非BB并保留终端例外提示，既有OXT敏感性见contact_type_decomposition。最近原子类型不等于具体键或能量。

大小层沿用Pedro非零接触数四分位：1–16、17–69、70–165、≥166。各层样本量见表。小界面的近接触比例离散且常为零；全体和层内、池化与界面均值不可混用。

```text
  scope      n  severe_overlap  has_clash  parser_warning_present  status_nonpassed  selected_model_count  best_model_count
    all 170644             464         91                       0                 0                     1                 1
nonzero  54109             464         30                       0                 0                     1                 1
   zero 116535               0         61                       0                 0                     1                 1
```

### raw：形成及大小

```text
best_model_pair_iptm     n  zero  median_length  zero_rate
      (-0.001, 0.05] 28612 25223       341.0000     0.8816
        (0.05, 0.08] 70578 54463       566.0000     0.7717
        (0.08, 0.12] 52383 31341       891.0000     0.5983
         (0.12, 0.2] 17211  5494      1442.0000     0.3192
          (0.2, 1.0]  1860    14       911.0000     0.0075
```

```text
size_stratum     n  median_contacts  median_close  mean_close  pooled_close  contacts  median_iptm  median_length  median_fraction
        1-16 13600           4.0000        0.0000      0.0264        0.0382     78889       0.0800       702.0000           0.0081
       17-69 13523          39.0000        0.0606      0.0651        0.0673    543380       0.0900       703.0000           0.0382
      70-165 13465         110.0000        0.0769      0.0777        0.0782   1513821       0.1000       732.0000           0.0797
        166+ 13521         263.0000        0.0847      0.0850        0.0870   4535574       0.1200       854.0000           0.1491
```

### qc：形成及大小

```text
best_model_pair_iptm     n  zero  median_length  zero_rate
      (-0.001, 0.05] 28600 25223       341.0000     0.8819
        (0.05, 0.08] 70470 54463       566.0000     0.7729
        (0.08, 0.12] 52216 31341       891.0000     0.6002
         (0.12, 0.2] 17050  5494      1442.0000     0.3222
          (0.2, 1.0]  1844    14       903.0000     0.0076
```

```text
size_stratum     n  median_contacts  median_close  mean_close  pooled_close  contacts  median_iptm  median_length  median_fraction
        1-16 13598           4.0000        0.0000      0.0264        0.0382     78859       0.0800       702.0000           0.0081
       17-69 13503          39.0000        0.0606      0.0651        0.0673    542363       0.0900       704.0000           0.0382
      70-165 13358         110.0000        0.0769      0.0777        0.0782   1500556       0.1000       732.0000           0.0797
        166+ 13186         262.0000        0.0847      0.0850        0.0870   4362191       0.1200       851.0000           0.1491
```

## 单变量与大小层内关联

### raw

```text
scope stratum                  feature     n  spearman  pearson
  raw     all               n_contacts 54109    0.5562   0.2503
  raw    1-16               n_contacts 13600    0.4213   0.2025
  raw   17-69               n_contacts 13523    0.1668   0.1413
  raw  70-165               n_contacts 13465    0.0931   0.0815
  raw    166+               n_contacts 13521    0.1850   0.1296
  raw     all total_interface_residues 54109    0.5323   0.2928
  raw    1-16 total_interface_residues 13600    0.3941   0.1783
  raw   17-69 total_interface_residues 13523    0.0502   0.0213
  raw  70-165 total_interface_residues 13465   -0.0685  -0.0903
  raw    166+ total_interface_residues 13521    0.0797   0.0922
  raw     all       interface_fraction 54109    0.5340   0.3211
  raw    1-16       interface_fraction 13600    0.3434   0.1486
  raw   17-69       interface_fraction 13523    0.1025   0.0851
  raw  70-165       interface_fraction 13465    0.0935   0.1135
  raw    166+       interface_fraction 13521    0.2146   0.2220
  raw     all       total_chain_length 54109   -0.0056   0.0061
  raw    1-16       total_chain_length 13600   -0.0307  -0.0238
  raw   17-69       total_chain_length 13523   -0.0825  -0.0677
  raw  70-165       total_chain_length 13465   -0.1294  -0.1126
  raw    166+       total_chain_length 13521   -0.1507  -0.1363
  raw     all       short_chain_length 54109   -0.0154  -0.0064
  raw    1-16       short_chain_length 13600   -0.0179  -0.0157
  raw   17-69       short_chain_length 13523   -0.0831  -0.0760
  raw  70-165       short_chain_length 13465   -0.1264  -0.1220
  raw    166+       short_chain_length 13521   -0.1677  -0.1534
  raw     all        long_chain_length 54109   -0.0025   0.0106
  raw    1-16        long_chain_length 13600   -0.0307  -0.0228
  raw   17-69        long_chain_length 13523   -0.0721  -0.0509
  raw  70-165        long_chain_length 13465   -0.1111  -0.0865
  raw    166+        long_chain_length 13521   -0.1251  -0.1022
  raw     all             length_ratio 54109   -0.0127  -0.0047
  raw    1-16             length_ratio 13600    0.0077   0.0129
  raw   17-69             length_ratio 13523   -0.0114  -0.0125
  raw  70-165             length_ratio 13465   -0.0167  -0.0150
  raw    166+             length_ratio 13521   -0.0305  -0.0259
  raw     all     best_model_pair_iptm 54109    0.1602   0.1514
  raw    1-16     best_model_pair_iptm 13600   -0.0246  -0.0199
  raw   17-69     best_model_pair_iptm 13523   -0.0282   0.0142
  raw  70-165     best_model_pair_iptm 13465   -0.0138   0.0854
  raw    166+     best_model_pair_iptm 13521    0.0138   0.1558
  raw     all              share_BB-BB 54109    0.4060   0.2360
  raw    1-16              share_BB-BB 13600    0.1366   0.0258
  raw   17-69              share_BB-BB 13523    0.0913   0.0847
  raw  70-165              share_BB-BB 13465    0.0693   0.0848
  raw    166+              share_BB-BB 13521    0.0592   0.0907
  raw     all              share_BB-SC 54109    0.3724   0.3430
  raw    1-16              share_BB-SC 13600    0.2620   0.1122
  raw   17-69              share_BB-SC 13523    0.1697   0.1756
  raw  70-165              share_BB-SC 13465    0.1315   0.1262
  raw    166+              share_BB-SC 13521    0.1354   0.1201
  raw     all              share_SC-SC 54109   -0.4184  -0.3800
  raw    1-16              share_SC-SC 13600   -0.2634  -0.1129
  raw   17-69              share_SC-SC 13523   -0.1919  -0.1952
  raw  70-165              share_SC-SC 13465   -0.1671  -0.1748
  raw    166+              share_SC-SC 13521   -0.1731  -0.1856
```

### qc

```text
scope stratum                  feature     n  spearman  pearson
   qc     all               n_contacts 53645    0.5580   0.2773
   qc    1-16               n_contacts 13598    0.4210   0.2024
   qc   17-69               n_contacts 13503    0.1667   0.1412
   qc  70-165               n_contacts 13358    0.0927   0.0811
   qc    166+               n_contacts 13186    0.1884   0.1549
   qc     all total_interface_residues 53645    0.5342   0.2969
   qc    1-16 total_interface_residues 13598    0.3939   0.1782
   qc   17-69 total_interface_residues 13503    0.0502   0.0213
   qc  70-165 total_interface_residues 13358   -0.0681  -0.0901
   qc    166+ total_interface_residues 13186    0.0823   0.0940
   qc     all       interface_fraction 53645    0.5354   0.3210
   qc    1-16       interface_fraction 13598    0.3433   0.1486
   qc   17-69       interface_fraction 13503    0.1019   0.0842
   qc  70-165       interface_fraction 13358    0.0943   0.1142
   qc    166+       interface_fraction 13186    0.2153   0.2223
   qc     all       total_chain_length 53645   -0.0063   0.0052
   qc    1-16       total_chain_length 13598   -0.0309  -0.0238
   qc   17-69       total_chain_length 13503   -0.0818  -0.0671
   qc  70-165       total_chain_length 13358   -0.1299  -0.1132
   qc    166+       total_chain_length 13186   -0.1510  -0.1366
   qc     all       short_chain_length 53645   -0.0159  -0.0072
   qc    1-16       short_chain_length 13598   -0.0180  -0.0158
   qc   17-69       short_chain_length 13503   -0.0826  -0.0756
   qc  70-165       short_chain_length 13358   -0.1275  -0.1229
   qc    166+       short_chain_length 13186   -0.1677  -0.1534
   qc     all        long_chain_length 53645   -0.0032   0.0099
   qc    1-16        long_chain_length 13598   -0.0308  -0.0228
   qc   17-69        long_chain_length 13503   -0.0713  -0.0504
   qc  70-165        long_chain_length 13358   -0.1112  -0.0867
   qc    166+        long_chain_length 13186   -0.1257  -0.1025
   qc     all             length_ratio 53645   -0.0126  -0.0046
   qc    1-16             length_ratio 13598    0.0076   0.0129
   qc   17-69             length_ratio 13503   -0.0113  -0.0125
   qc  70-165             length_ratio 13358   -0.0179  -0.0163
   qc    166+             length_ratio 13186   -0.0294  -0.0249
   qc     all     best_model_pair_iptm 53645    0.1605   0.1516
   qc    1-16     best_model_pair_iptm 13598   -0.0248  -0.0200
   qc   17-69     best_model_pair_iptm 13503   -0.0274   0.0148
   qc  70-165     best_model_pair_iptm 13358   -0.0139   0.0857
   qc    166+     best_model_pair_iptm 13186    0.0171   0.1594
   qc     all              share_BB-BB 53645    0.4064   0.2345
   qc    1-16              share_BB-BB 13598    0.1363   0.0257
   qc   17-69              share_BB-BB 13503    0.0908   0.0844
   qc  70-165              share_BB-BB 13358    0.0669   0.0824
   qc    166+              share_BB-BB 13186    0.0599   0.0914
   qc     all              share_BB-SC 53645    0.3739   0.3425
   qc    1-16              share_BB-SC 13598    0.2617   0.1122
   qc   17-69              share_BB-SC 13503    0.1695   0.1755
   qc  70-165              share_BB-SC 13358    0.1317   0.1263
   qc    166+              share_BB-SC 13186    0.1383   0.1232
   qc     all              share_SC-SC 53645   -0.4189  -0.3790
   qc    1-16              share_SC-SC 13598   -0.2632  -0.1128
   qc   17-69              share_SC-SC 13503   -0.1915  -0.1950
   qc  70-165              share_SC-SC 13358   -0.1647  -0.1727
   qc    166+              share_SC-SC 13186   -0.1759  -0.1882
```

Spearman与Pearson均为描述性；n_contacts与界面残基数有定义上的联系，类型组成总和为1。层内相关不能单独证明大小调整后的因果关系。

## 多变量描述模型

以近接触数为分子、总接触数为二项分母，使用logit均值模型；该似然按界面接触数加权，与前述界面等权均值的目标量不同。大小项为log1p(contact_count)及其以样本中位数中心化的平方。初始模型中接触数与原始界面残基比例的VIF偏高，原拟合另存`original_collinear_*`供审计；主模型先将界面比例对大小及总链长回归，使用残差`interface_spread_excess`表示“相对于相近大小和链长的额外界面覆盖”。残差化系数及尺度在interface_spread_residualization.csv；这是描述性条件变量，不等同于原始界面比例。structural模型包含pair ipTM、log总链长、log长短比和该残差；composition模型另加BB-BB和BB-SC份额，以SC-SC为参照，避免三份额全入模。所有自变量按该模型样本的标准差标准化。`odds_ratio_per_sd`仅为条件均值描述，不是独立接触事件的优势比；报告按蛋白对聚类的sandwich区间，并用Pearson离散度检查二项方差假设。分位预测差固定其他变量于中位数，单位百分点；不是干预效应。强共线性以VIF诊断，须谨慎解释个别系数。

```text
   scope       model                    term  beta_std  se_cluster  odds_ratio_per_sd  ci_low  ci_high     n  n_pairs  pearson_dispersion
     raw  structural               intercept   -2.6851      0.0073             0.0682  0.0672   0.0692 54109    44075              0.9289
     raw  structural               log_count    0.2819      0.0065             1.3256  1.3088   1.3426 54109    44075              0.9289
     raw  structural            log_count_sq   -0.0626      0.0077             0.9394  0.9253   0.9537 54109    44075              0.9289
     raw  structural    best_model_pair_iptm    0.0329      0.0013             1.0334  1.0308   1.0361 54109    44075              0.9289
     raw  structural        log_total_length   -0.0855      0.0024             0.9181  0.9137   0.9225 54109    44075              0.9289
     raw  structural           log_asymmetry    0.0183      0.0017             1.0185  1.0151   1.0219 54109    44075              0.9289
     raw  structural interface_spread_excess   -0.0097      0.0013             0.9904  0.9878   0.9930 54109    44075              0.9289
     raw composition               intercept   -2.7148      0.0050             0.0662  0.0656   0.0669 54109    44075              0.9232
     raw composition               log_count    0.2358      0.0049             1.2659  1.2538   1.2781 54109    44075              0.9232
     raw composition            log_count_sq   -0.0548      0.0039             0.9467  0.9395   0.9539 54109    44075              0.9232
     raw composition    best_model_pair_iptm    0.0313      0.0013             1.0318  1.0292   1.0344 54109    44075              0.9232
     raw composition        log_total_length   -0.0704      0.0023             0.9320  0.9278   0.9363 54109    44075              0.9232
     raw composition           log_asymmetry    0.0159      0.0017             1.0160  1.0127   1.0194 54109    44075              0.9232
     raw composition interface_spread_excess   -0.0038      0.0013             0.9962  0.9937   0.9988 54109    44075              0.9232
     raw composition             share_BB-BB    0.0447      0.0023             1.0457  1.0409   1.0505 54109    44075              0.9232
     raw composition             share_BB-SC    0.1378      0.0053             1.1477  1.1358   1.1597 54109    44075              0.9232
      qc  structural               intercept   -2.6998      0.0051             0.0672  0.0666   0.0679 53645    43773              0.9318
      qc  structural               log_count    0.2934      0.0047             1.3410  1.3288   1.3533 53645    43773              0.9318
      qc  structural            log_count_sq   -0.0766      0.0044             0.9262  0.9182   0.9343 53645    43773              0.9318
      qc  structural    best_model_pair_iptm    0.0333      0.0013             1.0339  1.0312   1.0366 53645    43773              0.9318
      qc  structural        log_total_length   -0.0840      0.0023             0.9195  0.9153   0.9237 53645    43773              0.9318
      qc  structural           log_asymmetry    0.0187      0.0017             1.0189  1.0154   1.0223 53645    43773              0.9318
      qc  structural interface_spread_excess   -0.0094      0.0014             0.9907  0.9880   0.9933 53645    43773              0.9318
      qc composition               intercept   -2.7230      0.0051             0.0657  0.0650   0.0663 53645    43773              0.9240
      qc composition               log_count    0.2423      0.0049             1.2741  1.2619   1.2865 53645    43773              0.9240
      qc composition            log_count_sq   -0.0606      0.0042             0.9412  0.9335   0.9489 53645    43773              0.9240
      qc composition    best_model_pair_iptm    0.0315      0.0013             1.0320  1.0293   1.0346 53645    43773              0.9240
      qc composition        log_total_length   -0.0701      0.0023             0.9323  0.9281   0.9366 53645    43773              0.9240
      qc composition           log_asymmetry    0.0160      0.0017             1.0161  1.0127   1.0195 53645    43773              0.9240
      qc composition interface_spread_excess   -0.0041      0.0014             0.9959  0.9932   0.9985 53645    43773              0.9240
      qc composition             share_BB-BB    0.0435      0.0024             1.0445  1.0396   1.0494 53645    43773              0.9240
      qc composition             share_BB-SC    0.1381      0.0054             1.1481  1.1361   1.1602 53645    43773              0.9240
pair_one  structural               intercept   -2.7000      0.0082             0.0672  0.0661   0.0683 44075    44075              0.9226
pair_one  structural               log_count    0.2872      0.0073             1.3327  1.3138   1.3518 44075    44075              0.9226
pair_one  structural            log_count_sq   -0.0613      0.0084             0.9406  0.9252   0.9562 44075    44075              0.9226
pair_one  structural    best_model_pair_iptm    0.0304      0.0014             1.0308  1.0280   1.0337 44075    44075              0.9226
pair_one  structural        log_total_length   -0.0827      0.0026             0.9206  0.9158   0.9254 44075    44075              0.9226
pair_one  structural           log_asymmetry    0.0172      0.0019             1.0173  1.0136   1.0210 44075    44075              0.9226
pair_one  structural interface_spread_excess   -0.0093      0.0014             0.9907  0.9880   0.9934 44075    44075              0.9226
pair_one composition               intercept   -2.7321      0.0055             0.0651  0.0644   0.0658 44075    44075              0.9176
pair_one composition               log_count    0.2390      0.0053             1.2700  1.2567   1.2833 44075    44075              0.9176
pair_one composition            log_count_sq   -0.0541      0.0042             0.9473  0.9396   0.9551 44075    44075              0.9176
pair_one composition    best_model_pair_iptm    0.0289      0.0014             1.0293  1.0265   1.0322 44075    44075              0.9176
pair_one composition        log_total_length   -0.0672      0.0025             0.9350  0.9305   0.9396 44075    44075              0.9176
pair_one composition           log_asymmetry    0.0147      0.0018             1.0148  1.0111   1.0185 44075    44075              0.9176
pair_one composition interface_spread_excess   -0.0035      0.0014             0.9965  0.9939   0.9991 44075    44075              0.9176
pair_one composition             share_BB-BB    0.0456      0.0025             1.0467  1.0416   1.0518 44075    44075              0.9176
pair_one composition             share_BB-SC    0.1448      0.0057             1.1558  1.1431   1.1688 44075    44075              0.9176
```

```text
   scope       model                    term    vif
     raw  structural               log_count 2.0276
     raw  structural            log_count_sq 1.8253
     raw  structural    best_model_pair_iptm 1.5374
     raw  structural        log_total_length 1.4269
     raw  structural           log_asymmetry 1.0801
     raw  structural interface_spread_excess 1.0023
     raw composition               log_count 2.4600
     raw composition            log_count_sq 2.0558
     raw composition    best_model_pair_iptm 1.5421
     raw composition        log_total_length 1.4477
     raw composition           log_asymmetry 1.0818
     raw composition interface_spread_excess 1.0246
     raw composition             share_BB-BB 1.2800
     raw composition             share_BB-SC 1.7399
      qc  structural               log_count 2.0341
      qc  structural            log_count_sq 1.8312
      qc  structural    best_model_pair_iptm 1.5361
      qc  structural        log_total_length 1.4229
      qc  structural           log_asymmetry 1.0804
      qc  structural interface_spread_excess 1.0022
      qc composition               log_count 2.4604
      qc composition            log_count_sq 2.0616
      qc composition    best_model_pair_iptm 1.5413
      qc composition        log_total_length 1.4441
      qc composition           log_asymmetry 1.0821
      qc composition interface_spread_excess 1.0238
      qc composition             share_BB-BB 1.2767
      qc composition             share_BB-SC 1.7392
pair_one  structural               log_count 2.0069
pair_one  structural            log_count_sq 1.8020
pair_one  structural    best_model_pair_iptm 1.6149
pair_one  structural        log_total_length 1.4858
pair_one  structural           log_asymmetry 1.0792
pair_one  structural interface_spread_excess 1.0030
pair_one composition               log_count 2.4583
pair_one composition            log_count_sq 2.0142
pair_one composition    best_model_pair_iptm 1.6186
pair_one composition        log_total_length 1.5057
pair_one composition           log_asymmetry 1.0808
pair_one composition interface_spread_excess 1.0243
pair_one composition             share_BB-BB 1.2827
pair_one composition             share_BB-SC 1.7183
```

```text
   scope       model                 feature     q25    q75  predicted_probability_q25  predicted_probability_q75  difference_pp
     raw  structural    best_model_pair_iptm  0.0800 0.1300                     0.0680                     0.0699         0.1914
     raw  structural        log_total_length  6.2823 6.9470                     0.0724                     0.0652        -0.7249
     raw  structural           log_asymmetry  0.2991 1.1022                     0.0680                     0.0697         0.1629
     raw  structural interface_spread_excess -0.0201 0.0133                     0.0690                     0.0684        -0.0566
     raw composition    best_model_pair_iptm  0.0800 0.1300                     0.0670                     0.0688         0.1795
     raw composition        log_total_length  6.2823 6.9470                     0.0707                     0.0648        -0.5886
     raw composition           log_asymmetry  0.2991 1.1022                     0.0671                     0.0685         0.1395
     raw composition interface_spread_excess -0.0201 0.0133                     0.0678                     0.0676        -0.0218
      qc  structural    best_model_pair_iptm  0.0800 0.1200                     0.0676                     0.0691         0.1537
      qc  structural        log_total_length  6.2804 6.9451                     0.0719                     0.0648        -0.7088
      qc  structural           log_asymmetry  0.2991 1.1022                     0.0676                     0.0693         0.1655
      qc  structural interface_spread_excess -0.0201 0.0132                     0.0686                     0.0680        -0.0545
      qc composition    best_model_pair_iptm  0.0800 0.1200                     0.0667                     0.0682         0.1433
      qc composition        log_total_length  6.2804 6.9451                     0.0704                     0.0646        -0.5846
      qc composition           log_asymmetry  0.2991 1.1022                     0.0669                     0.0683         0.1395
      qc composition interface_spread_excess -0.0201 0.0132                     0.0676                     0.0673        -0.0237
pair_one  structural    best_model_pair_iptm  0.0700 0.1200                     0.0669                     0.0687         0.1821
pair_one  structural        log_total_length  6.2757 6.9402                     0.0715                     0.0646        -0.6962
pair_one  structural           log_asymmetry  0.3009 1.1012                     0.0673                     0.0688         0.1509
pair_one  structural interface_spread_excess -0.0192 0.0127                     0.0682                     0.0677        -0.0537
pair_one composition    best_model_pair_iptm  0.0700 0.1200                     0.0659                     0.0677         0.1710
pair_one composition        log_total_length  6.2757 6.9402                     0.0698                     0.0642        -0.5574
pair_one composition           log_asymmetry  0.3009 1.1012                     0.0664                     0.0677         0.1273
pair_one composition interface_spread_excess -0.0192 0.0127                     0.0670                     0.0668        -0.0199
```

形成模型使用所有有效预测、结局n_contacts>0；同样报告pair聚类标准误，并不把零界面近接触比例设为0。此简单模型的Pearson离散度极高，系数不应用作校准概率或外部准确率，形成趋势主要看实际分组率。

```text
                term  beta_std  se_cluster  odds_ratio_per_sd   ci_low  ci_high      n  n_pairs  pearson_dispersion
           intercept   -0.5379      0.0079             0.5840   0.5750   0.5931 170644   113050            372.8502
best_model_pair_iptm    4.8021      0.0442           121.7612 111.6622 132.7735 170644   113050            372.8502
    log_total_length   -2.7815      0.0295             0.0619   0.0585   0.0656 170644   113050            372.8502
       log_asymmetry   -0.2356      0.0073             0.7901   0.7789   0.8015 170644   113050            372.8502
```

## 类型组成与类型内部距离

### raw

```text
scope stratum  type  n_interfaces  n_with_type  share_interface_mean  within_mean  contribution_mean  pooled_share  pooled_within
  raw    1-16 BB-BB         13600          959                0.0144       0.0018             0.0000        0.0233         0.0043
  raw    1-16 BB-SC         13600         7873                0.2599       0.0116             0.0038        0.3486         0.0211
  raw    1-16 SC-SC         13600        13091                0.7257       0.0408             0.0226        0.6280         0.0489
  raw   17-69 BB-BB         13523         8742                0.0780       0.0110             0.0019        0.0840         0.0263
  raw   17-69 BB-SC         13523        13516                0.4821       0.0424             0.0218        0.4891         0.0471
  raw   17-69 SC-SC         13523        13520                0.4398       0.1006             0.0415        0.4268         0.0986
  raw  70-165 BB-BB         13465        12822                0.1035       0.0196             0.0034        0.1040         0.0333
  raw  70-165 BB-SC         13465        13465                0.5167       0.0539             0.0281        0.5179         0.0545
  raw  70-165 SC-SC         13465        13465                0.3798       0.1247             0.0462        0.3781         0.1231
  raw    166+ BB-BB         13521        13469                0.1080       0.0289             0.0042        0.1081         0.0417
  raw    166+ BB-SC         13521        13521                0.5316       0.0564             0.0300        0.5352         0.0565
  raw    166+ SC-SC         13521        13521                0.3604       0.1423             0.0508        0.3567         0.1464
```

### qc

```text
scope stratum  type  n_interfaces  n_with_type  share_interface_mean  within_mean  contribution_mean  pooled_share  pooled_within
   qc    1-16 BB-BB         13598          958                0.0144       0.0018             0.0000        0.0233         0.0044
   qc    1-16 BB-SC         13598         7871                0.2599       0.0116             0.0038        0.3485         0.0211
   qc    1-16 SC-SC         13598        13089                0.7257       0.0408             0.0226        0.6281         0.0489
   qc   17-69 BB-BB         13503         8725                0.0779       0.0109             0.0019        0.0839         0.0261
   qc   17-69 BB-SC         13503        13496                0.4821       0.0423             0.0217        0.4891         0.0471
   qc   17-69 SC-SC         13503        13500                0.4400       0.1006             0.0415        0.4271         0.0986
   qc  70-165 BB-BB         13358        12717                0.1031       0.0195             0.0034        0.1036         0.0331
   qc  70-165 BB-SC         13358        13358                0.5167       0.0538             0.0280        0.5179         0.0544
   qc  70-165 SC-SC         13358        13358                0.3802       0.1248             0.0463        0.3785         0.1232
   qc    166+ BB-BB         13186        13135                0.1069       0.0284             0.0041        0.1058         0.0396
   qc    166+ BB-SC         13186        13186                0.5318       0.0562             0.0299        0.5356         0.0561
   qc    166+ SC-SC         13186        13186                0.3613       0.1429             0.0511        0.3586         0.1470
```

查看share、within、contribution三列才能判断大小趋势来自类型比例、类型内部几何或二者。`pooled_share×pooled_within`等于该类型的池化贡献；等权界面均值的share与within乘积一般不等于贡献均值。OXT按既有规则处理；本目录另附复用的Pedro OXT敏感性表。

## 大小层内低/高近接触类别及代表界面

主低/高组按完全相同的整数接触数分别抽取互不重叠的低/高近接触分位成员，每个接触数两组人数相同。每个精确接触数至少20个可用界面；若20%和80%阈值同为零，高组改用>0的候选；无法分开两组的接触数不参与。支持度见exact_size_match_support_*.csv。原宽层分位结果因最小层阈值并列而重叠，另存broad_stratum_low_high_diagnostic_*.csv，不用于主比较。精确匹配组再按四个展示层汇总链长、长度比、界面残基比例、ipTM、类型份额及类型内部比例。代表例先排除已有severe-overlap，再按距候选组四个数值指标中位数的标准化平方距离选取，每组最多3个不同蛋白对；这仅定义未来结构检查的可复现抽样，不是验证其正确性。

OXT分类敏感性（复用既有审计）：

```text
  dataset contact_type distance_bin  primary_count  OXT_backbone_count  change_all_fraction  n_contacts_with_OXT
Pedro AF3        BB-BB         2to3           9157                9357               0.0000                19787
Pedro AF3        BB-BB         3to4          26908               27374               0.0001                19787
Pedro AF3        BB-BB         4to5          49117               49900               0.0001                19787
Pedro AF3        BB-BB         5to6          96967               98274               0.0002                19787
Pedro AF3        BB-BB         6to8         510207              515378               0.0008                19787
Pedro AF3        BB-SC         2to3          71688               72252               0.0001                19787
Pedro AF3        BB-SC         3to4         193034              193802               0.0001                19787
Pedro AF3        BB-SC         4to5         315831              316735               0.0001                19787
Pedro AF3        BB-SC         5to6         539830              540616               0.0001                19787
Pedro AF3        BB-SC         6to8        2368549             2369135               0.0001                19787
Pedro AF3        BB-SC          lt2          15827               16024               0.0000                19787
Pedro AF3        SC-SC         2to3         138275              137511              -0.0001                19787
Pedro AF3        SC-SC         3to4         332667              331433              -0.0002                19787
Pedro AF3        SC-SC         4to5         364931              363244              -0.0003                19787
Pedro AF3        SC-SC         5to6         411949              409856              -0.0003                19787
Pedro AF3        SC-SC         6to8        1208970             1203213              -0.0009                19787
Pedro AF3        SC-SC          lt2          14998               14732              -0.0000                19787
Pedro AF3        BB-BB          lt2           2759                2828               0.0000                19787
```

### raw

```text
scope stratum group              feature    n      mean   median      q25       q75  n_exact_contact_counts
  raw    1-16  high           n_contacts 1482    8.8644   9.0000   6.0000   12.0000                      16
  raw    1-16  high       close_fraction 1482    0.1890   0.1538   0.1250    0.2000                      16
  raw    1-16  high   total_chain_length 1482  756.2328 676.0000 489.2500  943.7500                      16
  raw    1-16  high         length_ratio 1482    0.5429   0.5341   0.3382    0.7550                      16
  raw    1-16  high   interface_fraction 1482    0.0158   0.0131   0.0085    0.0200                      16
  raw    1-16  high best_model_pair_iptm 1482    0.0868   0.0800   0.0600    0.1000                      16
  raw    1-16  high          share_BB-BB 1482    0.0212   0.0000   0.0000    0.0000                      16
  raw    1-16  high          share_BB-SC 1482    0.3755   0.4000   0.2000    0.5556                      16
  raw    1-16  high          share_SC-SC 1482    0.6034   0.6000   0.4286    0.8000                      16
  raw    1-16  high         within_BB-BB 1482    0.0065   0.0000   0.0000    0.0000                      16
  raw    1-16  high         within_BB-SC 1482    0.0502   0.0000   0.0000    0.0000                      16
  raw    1-16  high         within_SC-SC 1482    0.2772   0.2500   0.1667    0.3333                      16
  raw    1-16   low           n_contacts 1482    8.8644   9.0000   6.0000   12.0000                      16
  raw    1-16   low       close_fraction 1482    0.0000   0.0000   0.0000    0.0000                      16
  raw    1-16   low   total_chain_length 1482  793.2314 716.0000 504.0000  984.2500                      16
  raw    1-16   low         length_ratio 1482    0.5282   0.5284   0.3152    0.7311                      16
  raw    1-16   low   interface_fraction 1482    0.0158   0.0129   0.0084    0.0203                      16
  raw    1-16   low best_model_pair_iptm 1482    0.0899   0.0900   0.0700    0.1100                      16
  raw    1-16   low          share_BB-BB 1482    0.0196   0.0000   0.0000    0.0000                      16
  raw    1-16   low          share_BB-SC 1482    0.3402   0.3333   0.1667    0.5000                      16
  raw    1-16   low          share_SC-SC 1482    0.6402   0.6429   0.4667    0.8333                      16
  raw    1-16   low         within_BB-BB 1482    0.0000   0.0000   0.0000    0.0000                      16
  raw    1-16   low         within_BB-SC 1482    0.0000   0.0000   0.0000    0.0000                      16
  raw    1-16   low         within_SC-SC 1482    0.0000   0.0000   0.0000    0.0000                      16
  raw    166+  high           n_contacts 1933  241.4801 228.0000 193.0000  280.0000                     203
  raw    166+  high       close_fraction 1933    0.1116   0.1095   0.1031    0.1176                     203
  raw    166+  high   total_chain_length 1933  810.5753 696.0000 514.0000  996.0000                     203
  raw    166+  high         length_ratio 1933    0.5392   0.5235   0.3389    0.7419                     203
  raw    166+  high   interface_fraction 1933    0.1607   0.1433   0.1029    0.2009                     203
  raw    166+  high best_model_pair_iptm 1933    0.1481   0.1100   0.0900    0.1600                     203
  raw    166+  high          share_BB-BB 1933    0.1190   0.1114   0.0654    0.1636                     203
  raw    166+  high          share_BB-SC 1933    0.5353   0.5385   0.5000    0.5721                     203
  raw    166+  high          share_SC-SC 1933    0.3458   0.3432   0.3027    0.3857                     203
  raw    166+  high         within_BB-BB 1933    0.0445   0.0312   0.0000    0.0769                     203
  raw    166+  high         within_BB-SC 1933    0.0726   0.0714   0.0544    0.0896                     203
  raw    166+  high         within_SC-SC 1933    0.1904   0.1871   0.1613    0.2182                     203
  raw    166+   low           n_contacts 1933  241.4801 228.0000 193.0000  280.0000                     203
  raw    166+   low       close_fraction 1933    0.0574   0.0588   0.0518    0.0637                     203
  raw    166+   low   total_chain_length 1933 1016.5535 923.0000 687.0000 1269.0000                     203
  raw    166+   low         length_ratio 1933    0.5521   0.5473   0.3629    0.7539                     203
  raw    166+   low   interface_fraction 1933    0.1348   0.1223   0.0911    0.1636                     203
  raw    166+   low best_model_pair_iptm 1933    0.1267   0.1200   0.1000    0.1500                     203
  raw    166+   low          share_BB-BB 1933    0.1048   0.0981   0.0653    0.1378                     203
  raw    166+   low          share_BB-SC 1933    0.5202   0.5208   0.4865    0.5533                     203
  raw    166+   low          share_SC-SC 1933    0.3750   0.3750   0.3307    0.4198                     203
  raw    166+   low         within_BB-BB 1933    0.0165   0.0000   0.0000    0.0312                     203
  raw    166+   low         within_BB-SC 1933    0.0405   0.0400   0.0283    0.0518                     203
  raw    166+   low         within_SC-SC 1933    0.0914   0.0909   0.0727    0.1111                     203
  raw   17-69  high           n_contacts 2684   40.1673  39.0000  27.0000   53.0000                      53
  raw   17-69  high       close_fraction 2684    0.1242   0.1176   0.1061    0.1373                      53
  raw   17-69  high   total_chain_length 2684  746.6021 665.5000 463.7500  929.2500                      53
  raw   17-69  high         length_ratio 2684    0.5320   0.5179   0.3225    0.7396                      53
  raw   17-69  high   interface_fraction 2684    0.0470   0.0389   0.0263    0.0591                      53
  raw   17-69  high best_model_pair_iptm 2684    0.0954   0.0900   0.0700    0.1100                      53
  raw   17-69  high          share_BB-BB 2684    0.0891   0.0476   0.0000    0.1381                      53
  raw   17-69  high          share_BB-SC 2684    0.5087   0.5135   0.4340    0.5882                      53
  raw   17-69  high          share_SC-SC 2684    0.4023   0.3943   0.3000    0.5000                      53
  raw   17-69  high         within_BB-BB 2684    0.0237   0.0000   0.0000    0.0000                      53
  raw   17-69  high         within_BB-SC 2684    0.0859   0.0800   0.0357    0.1282                      53
  raw   17-69  high         within_SC-SC 2684    0.1892   0.1818   0.1333    0.2353                      53
  raw   17-69   low           n_contacts 2684   40.1673  39.0000  27.0000   53.0000                      53
  raw   17-69   low       close_fraction 2684    0.0159   0.0175   0.0000    0.0303                      53
  raw   17-69   low   total_chain_length 2684  832.7925 747.5000 544.7500 1033.2500                      53
  raw   17-69   low         length_ratio 2684    0.5423   0.5293   0.3431    0.7486                      53
  raw   17-69   low   interface_fraction 2684    0.0438   0.0373   0.0261    0.0540                      53
  raw   17-69   low best_model_pair_iptm 2684    0.0962   0.0900   0.0700    0.1100                      53
  raw   17-69   low          share_BB-BB 2684    0.0717   0.0385   0.0000    0.1064                      53
  raw   17-69   low          share_BB-SC 2684    0.4605   0.4667   0.3871    0.5417                      53
  raw   17-69   low          share_SC-SC 2684    0.4678   0.4661   0.3636    0.5607                      53
  raw   17-69   low         within_BB-BB 2684    0.0034   0.0000   0.0000    0.0000                      53
  raw   17-69   low         within_BB-SC 2684    0.0101   0.0000   0.0000    0.0000                      53
  raw   17-69   low         within_SC-SC 2684    0.0244   0.0000   0.0000    0.0476                      53
  raw  70-165  high           n_contacts 2655  112.3341 110.0000  88.0000  135.0000                      96
  raw  70-165  high       close_fraction 2655    0.1153   0.1121   0.1048    0.1225                      96
  raw  70-165  high   total_chain_length 2655  761.0821 657.0000 482.0000  923.5000                      96
  raw  70-165  high         length_ratio 2655    0.5330   0.5219   0.3292    0.7397                      96
  raw  70-165  high   interface_fraction 2655    0.0951   0.0836   0.0587    0.1158                      96
  raw  70-165  high best_model_pair_iptm 2655    0.1181   0.1000   0.0800    0.1300                      96
  raw  70-165  high          share_BB-BB 2655    0.1112   0.0938   0.0436    0.1627                      96
  raw  70-165  high          share_BB-SC 2655    0.5295   0.5324   0.4846    0.5769                      96
  raw  70-165  high          share_SC-SC 2655    0.3593   0.3583   0.2985    0.4140                      96
  raw  70-165  high         within_BB-BB 2655    0.0340   0.0000   0.0000    0.0625                      96
  raw  70-165  high         within_BB-SC 2655    0.0796   0.0781   0.0545    0.1042                      96
  raw  70-165  high         within_SC-SC 2655    0.1892   0.1852   0.1500    0.2222                      96
  raw  70-165   low           n_contacts 2655  112.3341 110.0000  88.0000  135.0000                      96
  raw  70-165   low       close_fraction 2655    0.0436   0.0448   0.0370    0.0523                      96
  raw  70-165   low   total_chain_length 2655  899.8105 800.0000 592.0000 1112.5000                      96
  raw  70-165   low         length_ratio 2655    0.5410   0.5249   0.3359    0.7375                      96
  raw  70-165   low   interface_fraction 2655    0.0835   0.0754   0.0551    0.1026                      96
  raw  70-165   low best_model_pair_iptm 2655    0.1092   0.1000   0.0800    0.1300                      96
  raw  70-165   low          share_BB-BB 2655    0.0947   0.0803   0.0385    0.1328                      96
  raw  70-165   low          share_BB-SC 2655    0.5049   0.5057   0.4574    0.5549                      96
  raw  70-165   low          share_SC-SC 2655    0.4004   0.3986   0.3333    0.4643                      96
  raw  70-165   low         within_BB-BB 2655    0.0097   0.0000   0.0000    0.0000                      96
  raw  70-165   low         within_BB-SC 2655    0.0299   0.0278   0.0154    0.0441                      96
  raw  70-165   low         within_SC-SC 2655    0.0682   0.0682   0.0465    0.0889                      96
```

### qc

```text
scope stratum group              feature    n      mean   median      q25       q75  n_exact_contact_counts
   qc    1-16  high           n_contacts 1481    8.8609   9.0000   6.0000   12.0000                      16
   qc    1-16  high       close_fraction 1481    0.1890   0.1538   0.1250    0.2000                      16
   qc    1-16  high   total_chain_length 1481  755.9892 676.0000 489.0000  943.0000                      16
   qc    1-16  high         length_ratio 1481    0.5426   0.5335   0.3381    0.7545                      16
   qc    1-16  high   interface_fraction 1481    0.0158   0.0131   0.0085    0.0200                      16
   qc    1-16  high best_model_pair_iptm 1481    0.0868   0.0800   0.0600    0.1000                      16
   qc    1-16  high          share_BB-BB 1481    0.0212   0.0000   0.0000    0.0000                      16
   qc    1-16  high          share_BB-SC 1481    0.3752   0.4000   0.2000    0.5556                      16
   qc    1-16  high          share_SC-SC 1481    0.6036   0.6000   0.4286    0.8000                      16
   qc    1-16  high         within_BB-BB 1481    0.0065   0.0000   0.0000    0.0000                      16
   qc    1-16  high         within_BB-SC 1481    0.0502   0.0000   0.0000    0.0000                      16
   qc    1-16  high         within_SC-SC 1481    0.2774   0.2500   0.1667    0.3333                      16
   qc    1-16   low           n_contacts 1481    8.8609   9.0000   6.0000   12.0000                      16
   qc    1-16   low       close_fraction 1481    0.0000   0.0000   0.0000    0.0000                      16
   qc    1-16   low   total_chain_length 1481  793.4922 716.0000 504.0000  985.0000                      16
   qc    1-16   low         length_ratio 1481    0.5283   0.5286   0.3151    0.7311                      16
   qc    1-16   low   interface_fraction 1481    0.0158   0.0129   0.0084    0.0203                      16
   qc    1-16   low best_model_pair_iptm 1481    0.0899   0.0900   0.0700    0.1100                      16
   qc    1-16   low          share_BB-BB 1481    0.0196   0.0000   0.0000    0.0000                      16
   qc    1-16   low          share_BB-SC 1481    0.3403   0.3333   0.1667    0.5000                      16
   qc    1-16   low          share_SC-SC 1481    0.6401   0.6429   0.4667    0.8333                      16
   qc    1-16   low         within_BB-BB 1481    0.0000   0.0000   0.0000    0.0000                      16
   qc    1-16   low         within_BB-SC 1481    0.0000   0.0000   0.0000    0.0000                      16
   qc    1-16   low         within_SC-SC 1481    0.0000   0.0000   0.0000    0.0000                      16
   qc    166+  high           n_contacts 1888  240.8978 228.0000 193.0000  279.0000                     200
   qc    166+  high       close_fraction 1888    0.1116   0.1095   0.1030    0.1176                     200
   qc    166+  high   total_chain_length 1888  809.1674 694.5000 513.0000  993.5000                     200
   qc    166+  high         length_ratio 1888    0.5386   0.5231   0.3382    0.7401                     200
   qc    166+  high   interface_fraction 1888    0.1609   0.1433   0.1035    0.2009                     200
   qc    166+  high best_model_pair_iptm 1888    0.1487   0.1100   0.0900    0.1600                     200
   qc    166+  high          share_BB-BB 1888    0.1180   0.1106   0.0647    0.1628                     200
   qc    166+  high          share_BB-SC 1888    0.5355   0.5385   0.5016    0.5721                     200
   qc    166+  high          share_SC-SC 1888    0.3464   0.3448   0.3033    0.3860                     200
   qc    166+  high         within_BB-BB 1888    0.0438   0.0303   0.0000    0.0755                     200
   qc    166+  high         within_BB-SC 1888    0.0721   0.0714   0.0541    0.0885                     200
   qc    166+  high         within_SC-SC 1888    0.1912   0.1877   0.1616    0.2188                     200
   qc    166+   low           n_contacts 1888  240.8978 228.0000 193.0000  279.0000                     200
   qc    166+   low       close_fraction 1888    0.0573   0.0588   0.0518    0.0637                     200
   qc    166+   low   total_chain_length 1888 1013.6854 918.0000 686.7500 1263.0000                     200
   qc    166+   low         length_ratio 1888    0.5509   0.5476   0.3625    0.7530                     200
   qc    166+   low   interface_fraction 1888    0.1351   0.1227   0.0915    0.1641                     200
   qc    166+   low best_model_pair_iptm 1888    0.1264   0.1200   0.1000    0.1500                     200
   qc    166+   low          share_BB-BB 1888    0.1039   0.0976   0.0646    0.1366                     200
   qc    166+   low          share_BB-SC 1888    0.5204   0.5210   0.4867    0.5534                     200
   qc    166+   low          share_SC-SC 1888    0.3757   0.3755   0.3313    0.4205                     200
   qc    166+   low         within_BB-BB 1888    0.0163   0.0000   0.0000    0.0303                     200
   qc    166+   low         within_BB-SC 1888    0.0403   0.0400   0.0280    0.0517                     200
   qc    166+   low         within_SC-SC 1888    0.0916   0.0917   0.0729    0.1111                     200
   qc   17-69  high           n_contacts 2682   40.1506  39.0000  27.0000   53.0000                      53
   qc   17-69  high       close_fraction 2682    0.1242   0.1176   0.1061    0.1373                      53
   qc   17-69  high   total_chain_length 2682  747.6868 666.0000 465.0000  930.7500                      53
   qc   17-69  high         length_ratio 2682    0.5324   0.5195   0.3231    0.7406                      53
   qc   17-69  high   interface_fraction 2682    0.0469   0.0388   0.0262    0.0590                      53
   qc   17-69  high best_model_pair_iptm 2682    0.0954   0.0900   0.0700    0.1100                      53
   qc   17-69  high          share_BB-BB 2682    0.0889   0.0476   0.0000    0.1379                      53
   qc   17-69  high          share_BB-SC 2682    0.5087   0.5135   0.4340    0.5882                      53
   qc   17-69  high          share_SC-SC 2682    0.4024   0.3943   0.3004    0.5000                      53
   qc   17-69  high         within_BB-BB 2682    0.0234   0.0000   0.0000    0.0000                      53
   qc   17-69  high         within_BB-SC 2682    0.0858   0.0800   0.0357    0.1274                      53
   qc   17-69  high         within_SC-SC 2682    0.1894   0.1818   0.1333    0.2353                      53
   qc   17-69   low           n_contacts 2682   40.1506  39.0000  27.0000   53.0000                      53
   qc   17-69   low       close_fraction 2682    0.0159   0.0175   0.0000    0.0303                      53
   qc   17-69   low   total_chain_length 2682  833.1044 748.0000 546.0000 1033.7500                      53
   qc   17-69   low         length_ratio 2682    0.5422   0.5293   0.3431    0.7484                      53
   qc   17-69   low   interface_fraction 2682    0.0438   0.0373   0.0261    0.0540                      53
   qc   17-69   low best_model_pair_iptm 2682    0.0962   0.0900   0.0700    0.1100                      53
   qc   17-69   low          share_BB-BB 2682    0.0716   0.0385   0.0000    0.1064                      53
   qc   17-69   low          share_BB-SC 2682    0.4605   0.4667   0.3871    0.5417                      53
   qc   17-69   low          share_SC-SC 2682    0.4679   0.4667   0.3636    0.5609                      53
   qc   17-69   low         within_BB-BB 2682    0.0034   0.0000   0.0000    0.0000                      53
   qc   17-69   low         within_BB-SC 2682    0.0101   0.0000   0.0000    0.0000                      53
   qc   17-69   low         within_SC-SC 2682    0.0245   0.0000   0.0000    0.0476                      53
   qc  70-165  high           n_contacts 2636  112.2394 110.0000  88.0000  135.0000                      96
   qc  70-165  high       close_fraction 2636    0.1153   0.1121   0.1048    0.1226                      96
   qc  70-165  high   total_chain_length 2636  762.4310 660.0000 483.0000  925.2500                      96
   qc  70-165  high         length_ratio 2636    0.5317   0.5216   0.3272    0.7382                      96
   qc  70-165  high   interface_fraction 2636    0.0950   0.0833   0.0587    0.1157                      96
   qc  70-165  high best_model_pair_iptm 2636    0.1184   0.1000   0.0800    0.1300                      96
   qc  70-165  high          share_BB-BB 2636    0.1106   0.0933   0.0431    0.1613                      96
   qc  70-165  high          share_BB-SC 2636    0.5291   0.5314   0.4840    0.5766                      96
   qc  70-165  high          share_SC-SC 2636    0.3603   0.3595   0.2991    0.4156                      96
   qc  70-165  high         within_BB-BB 2636    0.0335   0.0000   0.0000    0.0625                      96
   qc  70-165  high         within_BB-SC 2636    0.0792   0.0779   0.0541    0.1034                      96
   qc  70-165  high         within_SC-SC 2636    0.1897   0.1852   0.1515    0.2222                      96
   qc  70-165   low           n_contacts 2636  112.2394 110.0000  88.0000  135.0000                      96
   qc  70-165   low       close_fraction 2636    0.0436   0.0448   0.0370    0.0523                      96
   qc  70-165   low   total_chain_length 2636  902.7272 803.0000 593.0000 1115.2500                      96
   qc  70-165   low         length_ratio 2636    0.5407   0.5252   0.3345    0.7377                      96
   qc  70-165   low   interface_fraction 2636    0.0834   0.0751   0.0550    0.1023                      96
   qc  70-165   low best_model_pair_iptm 2636    0.1093   0.1000   0.0800    0.1300                      96
   qc  70-165   low          share_BB-BB 2636    0.0947   0.0800   0.0381    0.1330                      96
   qc  70-165   low          share_BB-SC 2636    0.5046   0.5052   0.4571    0.5545                      96
   qc  70-165   low          share_SC-SC 2636    0.4008   0.4000   0.3354    0.4648                      96
   qc  70-165   low         within_BB-BB 2636    0.0096   0.0000   0.0000    0.0000                      96
   qc  70-165   low         within_BB-SC 2636    0.0299   0.0278   0.0154    0.0441                      96
   qc  70-165   low         within_SC-SC 2636    0.0682   0.0682   0.0465    0.0889                      96
```

代表界面及来源：

```text
size_stratum class_label  rank     interface_id      pair_key    job_key protein_a protein_b                                                       protein_name_a                                                    protein_name_b  contact_count  close_fraction   iptm                                       archive                                                           cif_member  representative_distance
        1-16        high     1  a027_j0001:r168 MG_054|MG_303 a027_j0001    MG_303    MG_054             metal ion ABC transporter, ATP-binding protein, putative              transcription antitermination protein NusG, putative              7          0.1429 0.0800 AF3_completed_zips/folds_2025_03_25_17_17.zip   250310_mgen_allbyall_652/fold_250310_mgen_allbyall_652_model_0.cif                   0.1326
        1-16        high     2 a022_j0029:r2381 MG_332|MG_461 a022_j0029    MG_461    MG_332                                                    HD domain protein                             expressed protein of unknown function              7          0.1429 0.0800 AF3_completed_zips/folds_2025_03_21_22_23.zip   250310_mgen_allbyall_460/fold_250310_mgen_allbyall_460_model_0.cif                   0.1333
        1-16        high     3 a040_j0034:r2635 MG_123|MG_248 a040_j0034    MG_123    MG_248                                       conserved hypothetical protein                                    conserved hypothetical protein              7          0.1429 0.0800 AF3_completed_zips/folds_2025_03_31_00_55.zip   250310_mgen_allbyall_745/fold_250310_mgen_allbyall_745_model_0.cif                   0.1334
        1-16         low     1  a015_j0003:r364 MG_008|MG_108 a015_j0003    MG_108    MG_008                                     protein phosphatase 2C, putative                                     tRNA modification GTPase TrmE              9          0.0000 0.0900 AF3_completed_zips/folds_2025_03_18_20_17.zip 250310_mgen_allbyall_1830/fold_250310_mgen_allbyall_1830_model_0.cif                   0.0008
        1-16         low     2  a010_j0008:r700 MG_107|MG_249 a010_j0008    MG_249    MG_107                                     RNA polymerase sigma factor RpoD                                                  guanylate kinase              9          0.0000 0.0900 AF3_completed_zips/folds_2025_03_17_19_48.zip   250310_mgen_allbyall_149/fold_250310_mgen_allbyall_149_model_0.cif                   0.0039
        1-16         low     3  a010_j0003:r257 MG_281|MG_296 a010_j0003    MG_296    MG_281                                       conserved hypothetical protein                                    conserved hypothetical protein              9          0.0000 0.0900 AF3_completed_zips/folds_2025_03_17_19_48.zip   250310_mgen_allbyall_144/fold_250310_mgen_allbyall_144_model_0.cif                   0.0042
        166+        high     1 a036_j0038:r3154 MG_094|MG_106 a036_j0038    MG_094    MG_106                                             replicative DNA helicase                                               peptide deformylase            218          0.1101 0.1200 AF3_completed_zips/folds_2025_03_26_17_22.zip   250310_mgen_allbyall_518/fold_250310_mgen_allbyall_518_model_0.cif                   0.0189
        166+        high     2 a057_j0022:r1763 MG_130|MG_524 a057_j0022    MG_524    MG_130                                       conserved hypothetical protein                                       uncharacterized domain HDIG            212          0.1085 0.1200 AF3_completed_zips/folds_2025_04_04_17_15.zip 250310_mgen_allbyall_1723/fold_250310_mgen_allbyall_1723_model_0.cif                   0.0569
        166+        high     3 a031_j0035:r2602 MG_011|MG_297 a031_j0035    MG_297    MG_011                     signal recognition particle-docking protein FtsY                                    conserved hypothetical protein            223          0.1076 0.1100 AF3_completed_zips/folds_2025_03_25_23_16.zip   250310_mgen_allbyall_626/fold_250310_mgen_allbyall_626_model_0.cif                   0.0609
        166+         low     1 a024_j0020:r1644 MG_138|MG_410 a024_j0020    MG_138    MG_410                                             GTP-binding protein LepA                    phosphate ABC transporter, ATP-binding protein            239          0.0586 0.1200 AF3_completed_zips/folds_2025_03_24_17_32.zip   250310_mgen_allbyall_611/fold_250310_mgen_allbyall_611_model_0.cif                   0.0165
        166+         low     2  a062_j0008:r669 MG_043|MG_185 a062_j0008    MG_185    MG_043                                                lipoprotein, putative spermidine/putrescine ABC transporter, permease protein, putative            235          0.0596 0.1200 AF3_completed_zips/folds_2025_04_06_03_25.zip 250310_mgen_allbyall_1639/fold_250310_mgen_allbyall_1639_model_0.cif                   0.0229
        166+         low     3 a052_j0013:r1007 MG_194|MG_305 a052_j0013    MG_305    MG_194                                               chaperone protein DnaK                       phenylalanyl-tRNA synthetase, alpha subunit            214          0.0607 0.1200 AF3_completed_zips/folds_2025_04_02_19_28.zip 250310_mgen_allbyall_1394/fold_250310_mgen_allbyall_1394_model_0.cif                   0.0525
       17-69        high     1  a022_j0003:r292 MG_070|MG_181 a022_j0003    MG_181    MG_070                          metal ion ABC transporter, permease protein                                              ribosomal protein S2             34          0.1176 0.0900 AF3_completed_zips/folds_2025_03_21_22_23.zip   250310_mgen_allbyall_424/fold_250310_mgen_allbyall_424_model_0.cif                   0.0437
       17-69        high     2  a018_j0009:r759 MG_223|MG_273 a018_j0009    MG_223    MG_273                                       conserved hypothetical protein                 pyruvate dehydrogenase component E1, beta subunit             43          0.1163 0.0900 AF3_completed_zips/folds_2025_03_21_17_04.zip   250310_mgen_allbyall_370/fold_250310_mgen_allbyall_370_model_0.cif                   0.0489
       17-69        high     3  a065_j0002:r176 MG_035|MG_118 a065_j0002    MG_035    MG_118                                             histidyl-tRNA synthetase                                           UDP-glucose 4-epimerase             42          0.1190 0.0900 AF3_completed_zips/folds_2025_04_07_17_29.zip 250310_mgen_allbyall_1323/fold_250310_mgen_allbyall_1323_model_0.cif                   0.0511
       17-69         low     1  a064_j0005:r556 MG_249|MG_252 a064_j0005    MG_249    MG_252                                     RNA polymerase sigma factor RpoD                       RNA methyltransferase, TrmH family, group 3             43          0.0233 0.0900 AF3_completed_zips/folds_2025_04_07_16_56.zip 250310_mgen_allbyall_1446/fold_250310_mgen_allbyall_1446_model_0.cif                   0.0595
       17-69         low     2 a046_j0020:r1545 MG_227|MG_329 a046_j0020    MG_329    MG_227                                             GTP-binding protein engA                                              thymidylate synthase             43          0.0233 0.0900 AF3_completed_zips/folds_2025_03_31_18_09.zip   250310_mgen_allbyall_886/fold_250310_mgen_allbyall_886_model_0.cif                   0.0599
       17-69         low     3 a065_j0037:r2877 MG_200|MG_482 a065_j0037    MG_482    MG_200                                 holo-(acyl-carrier-protein) synthase                                               DnaJ domain protein             42          0.0238 0.0900 AF3_completed_zips/folds_2025_04_07_17_29.zip 250310_mgen_allbyall_1358/fold_250310_mgen_allbyall_1358_model_0.cif                   0.0606
      70-165        high     1 a046_j0039:r2884 MG_327|MG_372 a046_j0039    MG_372    MG_327                 thiamine biosynthesis/tRNA modification protein ThiI                                 hydrolase, alpha/beta fold family            109          0.1101 0.0900 AF3_completed_zips/folds_2025_03_31_18_09.zip   250310_mgen_allbyall_910/fold_250310_mgen_allbyall_910_model_0.cif                   0.0539
      70-165        high     2  a008_j0009:r717 MG_057|MG_329 a008_j0009    MG_329    MG_057                                             GTP-binding protein engA                                        small primase-like protein            109          0.1101 0.0900 AF3_completed_zips/folds_2025_03_14_22_45.zip 250310_mgen_allbyall_1790/fold_250310_mgen_allbyall_1790_model_0.cif                   0.0585
      70-165        high     3 a010_j0028:r2514 MG_138|MG_233 a010_j0028    MG_138    MG_233                                             GTP-binding protein LepA                                    conserved hypothetical protein            116          0.1121 0.0900 AF3_completed_zips/folds_2025_03_17_19_48.zip   250310_mgen_allbyall_189/fold_250310_mgen_allbyall_189_model_0.cif                   0.0650
      70-165         low     1 a040_j0036:r2783 MG_042|MG_304 a040_j0036    MG_042    MG_304 spermidine/putrescine ABC transporter, ATP-binding protein, putative          metal ion ABC transporter, ATP binding protein, putative            112          0.0446 0.1000 AF3_completed_zips/folds_2025_03_31_00_55.zip   250310_mgen_allbyall_747/fold_250310_mgen_allbyall_747_model_0.cif                   0.0048
      70-165         low     2  a026_j0004:r337 MG_231|MG_446 a026_j0004    MG_446    MG_231                                                ribosomal protein S16                 ribonucleoside-diphosphate reductase, alpha chain            109          0.0459 0.1000 AF3_completed_zips/folds_2025_03_24_20_04.zip   250310_mgen_allbyall_695/fold_250310_mgen_allbyall_695_model_0.cif                   0.0064
      70-165         low     3 a028_j0019:r1281 MG_108|MG_200 a028_j0019    MG_200    MG_108                                                  DnaJ domain protein                                  protein phosphatase 2C, putative            110          0.0455 0.1000 AF3_completed_zips/folds_2025_03_25_18_34.zip   250310_mgen_allbyall_790/fold_250310_mgen_allbyall_790_model_0.cif                   0.0153
```

## QC、重复观察及极端大小稳健性

```text
        scope      unit              feature     n  mean_close  spearman  pearson
          raw interface           n_contacts 54109      0.0635    0.5562   0.2503
          raw interface best_model_pair_iptm 54109      0.0635    0.1602   0.1514
          raw interface   total_chain_length 54109      0.0635   -0.0056   0.0061
          raw interface   interface_fraction 54109      0.0635    0.5340   0.3211
          raw      pair           n_contacts 44075      0.0625    0.5573   0.2646
          raw      pair best_model_pair_iptm 44075      0.0625    0.1579   0.1474
          raw      pair   total_chain_length 44075      0.0625    0.0012   0.0123
          raw      pair   interface_fraction 44075      0.0625    0.5355   0.3240
           qc interface           n_contacts 53645      0.0633    0.5580   0.2773
           qc interface best_model_pair_iptm 53645      0.0633    0.1605   0.1516
           qc interface   total_chain_length 53645      0.0633   -0.0063   0.0052
           qc interface   interface_fraction 53645      0.0633    0.5354   0.3210
           qc      pair           n_contacts 43773      0.0623    0.5590   0.2796
           qc      pair best_model_pair_iptm 43773      0.0623    0.1576   0.1473
           qc      pair   total_chain_length 43773      0.0623   -0.0000   0.0110
           qc      pair   interface_fraction 43773      0.0623    0.5370   0.3241
contacts_ge20 interface           n_contacts 39394      0.0766    0.2883   0.1850
contacts_ge20 interface best_model_pair_iptm 39394      0.0766    0.0726   0.1365
contacts_ge20 interface   total_chain_length 39394      0.0766   -0.0697  -0.0514
contacts_ge20 interface   interface_fraction 39394      0.0766    0.2652   0.2338
contacts_ge20      pair           n_contacts 32622      0.0761    0.2818   0.1938
contacts_ge20      pair best_model_pair_iptm 32622      0.0761    0.0672   0.1260
contacts_ge20      pair   total_chain_length 32622      0.0761   -0.0671  -0.0466
contacts_ge20      pair   interface_fraction 32622      0.0761    0.2600   0.2311
    trim_5_95 interface           n_contacts 51411      0.0621    0.5537   0.3359
    trim_5_95 interface best_model_pair_iptm 51411      0.0621    0.1360   0.1352
    trim_5_95 interface   total_chain_length 51411      0.0621   -0.0229  -0.0106
    trim_5_95 interface   interface_fraction 51411      0.0621    0.5284   0.3228
    trim_5_95      pair           n_contacts 42453      0.0613    0.5552   0.3392
    trim_5_95      pair best_model_pair_iptm 42453      0.0613    0.1340   0.1319
    trim_5_95      pair   total_chain_length 42453      0.0613   -0.0170  -0.0049
    trim_5_95      pair   interface_fraction 42453      0.0613    0.5298   0.3254
```

raw是所有有效观察；qc只排除既有severe-overlap；contacts_ge20和trim_5_95用于检查小分母/巨大界面影响。pair单位先在相同蛋白对的模型/作业间平均，再计算相关。结果若随单位或QC变化，不能称独立PPI的普遍性质。蛋白参与频率见protein_frequency.csv，重复pair完整分布见pair_multiplicity.csv。没有新增结构排除阈值用于主结果。

## 图件

图01为ipTM十分位零接触率；02为大小层等权近接触比例；03–05为大小层内ipTM、链长、界面残基比例的分箱描述，04b补充长短链比（点为分箱界面均值，图例显示层样本量）；06为三类组成；07为类型内部近接触比例；08为相同接触数低/高类别结构特征；09为原始、QC和大小端点敏感性。PNG/PDF成对保存；分箱图是描述，不作逐点显著性检验。

## 解释边界与下一步

这些是AF3预测结构内的描述性关联；没有实测结合标签、独立实验结构或能量评分。低近接触比例不等于错误界面，高ipTM也不保证真实相互作用。类型标签只代表最近原子，不能直接断言氢键或结合自由能。注释可能来自自动数据库文本，功能类别如要做推断应先核对覆盖、歧义和外部实验证据。下一步可优先检查代表列表的几何组织及严重重叠个案，再将相同蛋白对的独立实验结构纳入配对分析。

## 复现与验证

`run.py`只读原Parquet/CSV汇总，输入SHA256见input_manifest.csv。程序核对170644总界面、54109非零、116535零、464个既有severe-overlap；所有界面的三类接触数及3–4 Å数精确加和为旧总数。作业日志保存在logs/，输出文件及哈希在OUTPUT_FILES.csv。没有读取Pedro ZIP/CIF，没有重新运行AF3，既有结果未改动。

## 八个问题的具体解读

### 1. 哪些结构特征关联最强？

在54,109个非零界面中，近接触比例与接触数的Spearman相关为0.556，与界面残基比例为0.534；两者在定义和几何上强相关，不能算独立发现。pair ipTM的边际相关为0.160。总链长的全体边际相关约−0.006，掩盖层内负关联：在四个接触数层内依次约−0.031、−0.083、−0.129、−0.151。极小界面（1–16接触）近接触比例常为零，整数分母的离散性很强。

全体预测的零接触率为116,535/170,644=68.29%。既有ipTM区间中，≤0.05组零率88.16%，>0.20组为0.75%；这描述AF3预测的形成状态与其自身置信度同向变化，不是外部结合准确率。形成logit模型的Pearson离散度约372.9，提示极端预测/校准不佳；其系数和优势比不作为可靠定量结论。零界面的近接触比例始终保持NA。

### 2. 大小控制后还剩什么？

接触数与近接触比例的四层界面均值依次为2.64%、6.51%、7.77%、8.50%，说明明显非线性。层内pair ipTM关联远弱于边际关联；多变量描述模型中，将ipTM从第25到第75百分位（0.08→0.13）变化、其他变量固定中位数，预测近接触比例仅约+0.19个百分点。按蛋白对聚类的模型区间很窄，但这是预测数据的条件关联，不是独立实验PPI层面的准确度或生物学效应。

大小之外，总链长呈较稳定负向条件关联：其log长度四分位变化对应模型约−0.72个百分点。长短链比的条件关联更小，约+0.16个百分点。原界面残基比例与大小共线性偏高；将其对大小和总链长残差化后，额外覆盖比例四分位变化对应约−0.057个百分点。此残差结果不能解释成提高原始界面比例会降低近接触，也不能作为因果模型。主模型VIF最高2.46；初始共线模型及其VIF另存诊断。

### 3. 大小趋势来自类型组成、类型内部几何，还是两者？

**两者都有。** 从1–16到≥166接触层，界面等权的SC-SC份额约72.6%→36.0%，BB-SC约26.0%→53.2%，BB-BB约1.4%→10.8%；类型组成变化很大。同时在存在该类型的界面内，BB-SC的3–4 Å比例均值约1.16%→5.64%，SC-SC约4.08%→14.23%，BB-BB约0.18%→2.89%。类型内部近距离几何也明显随大小变动，不能只用份额解释。

最大层的三类对总近接触比例的界面等权贡献约0.42、3.00、5.08个百分点；份额高的BB-SC并非近接触贡献唯一来源。非常小界面的类型内部比例有很多零分母/整数边界，解释要更保守。OXT改归BB的既有敏感性最大改变约0.086个百分点（任一类型/距离箱的全接触份额），不足以解释主要层间差别。

### 4. 低近接触界面是否集中于某些结构类别？

主类别比较在**完全相同的整数接触数**上抽取互不重叠的低/高组，共有368个可比较接触数，每组8,754个界面；上限441接触，故不能代表更大的极端界面。1–16层因分位并列采用低组为零、高组为正的规则，属于特殊离散比较。四层各组接触数中位数严格相同，但组内仍可能存在蛋白功能、状态和模型质量差异。

在≥166层，低/高组总链长中位数约923/696 aa，界面残基比例约12.23%/14.33%，pair ipTM中位数约0.12/0.11；较低近接触组并非简单低ipTM组。其BB-SC份额约52.1%/53.8%、SC-SC份额37.5%/34.3%，而类型内部BB-SC近接触比例约4.0%/7.1%、SC-SC约9.1%/18.7%。类型内部差别较明显，但分组变量本身是这些类型近接触贡献之和，不能把这种联系误说成独立预测因子。已有注释字段覆盖全部记录，不过文本是自动/来源描述，未验证足以定义可靠的跨蛋白功能类别；目前先保留原注释和代表例名称，不根据ID造功能标签。

### 5. QC与重复蛋白对会改变结论吗？

464个既有severe-overlap均落在非零界面；排除后仍有53,645个非零界面。原始与QC下的界面平均近接触比例分别6.350%和6.334%；接触数Spearman相关0.556/0.558，pair ipTM相关0.160/0.160。主要方向变化很小。去掉接触数<20后，近接触比例均值升至7.658%，ipTM的Spearman相关降至0.073；显示小分母和大小组成对全体边际关联的重要性。

170,644条观察只对应113,050个无向蛋白对，其中43,190个对重复出现，单对最多8条记录。非零界面54,109条对应44,075个蛋白对。pair等权的主要Spearman相关与界面等权接近：接触数0.557、ipTM0.158、界面残基比例0.536。但同一蛋白可出现在许多pair中，pair等权仍不等于独立生物学样本。蛋白频率及重复数文件可用于后续更严格聚类。

### 6. ipTM在大小之外有多少信息？

ipTM对**是否形成预测界面**有很强描述性区分，但对**已形成界面中3–4 Å比例**的额外信息较弱。层内相关与条件模型效应明显小于边际相关，排除severe-overlap及pair聚合后仍如此。ipTM来自同一AF3预测流程，不能作为外部结构正确性的验证标签，也不能用本轮模型建立真实结合概率。

### 7. 哪些界面值得下一步看？

`representative_interfaces.csv`每个大小层的低/高组最多3个非severe例，按预先固定的数值中位距离规则选出，并提供protein pair、job、原注释和archive/CIF成员。可先看最大层：低组MG_138|MG_410、MG_043|MG_185、MG_194|MG_305；高组MG_094|MG_106、MG_130|MG_524、MG_011|MG_297。它们分别位于低/高类中，不能视为配对实验对照或质量裁决。结构检查应确认链组织、端点、局部遮挡和是否由重复模型驱动；本轮未重新读取归档。

### 8. 关于物理packing能推断什么？

可可靠描述的是：Pedro预测界面的大小、类型组成与最近原子距离分布协同变化；不同大小/链长类别具有不同的3–4 Å比例。无法由这些内部关联证明界面正确或错误，也不能推出普遍AF3 packing缺陷、氢键或结合能偏差。此前配对FoldBench/SAFES显示AF3可有明显接触身份重排，却未出现一致的大幅净近接触缺口；Pedro的组间差异仍可能涉及蛋白集合、预测状态和选择机制。最有价值的下一步是为代表蛋白对寻找可信实验对应，并先核对预测装配，再研究具体局部相互作用。
