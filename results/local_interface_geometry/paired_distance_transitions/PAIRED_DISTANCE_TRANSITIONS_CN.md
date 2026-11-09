# 配对 Experimental–AF3 接触距离转移分析

**核心结果：近接触身份会变化，但新增接触大致补偿流出。FoldBench/SAFES实验(3,4] Å接触的池化保留率为54.0%/75.5%；近接触绝对数净增67/6对。界面等权近接触比例差分别为−0.389/+0.185个百分点，不支持两组共同存在大幅净缺口。远距流失与整体装配偏差关联，不能全部归因为局部packing。**

全部计算、输入不变及输出一致性验证通过；16组图同时提供PNG/PDF。本轮没有运行Pedro或覆盖既有分析。详细的六问结论和建议位于报告末尾。

## 方法、数据及映射审计

沿用既有链和残基对应，只分析任一侧d_min≤8 Å的残基对并集。d_min为跨链残基间非H、非D原子的最小欧氏距离，每个残基对一次。所有原始有限距离保留，包括>8 Å；未进行全链/全结构长距离原子矩阵计算。直接复用已保存距离，并对并集内每对残基的小原子矩阵补齐双侧最近原子名称和元素。

FoldBench链映射歧义的23个界面不强制配对，147个不可用AF3不当零。SAFES沿用序列/元数据支持的全部可对应蛋白链对，不等同于人工确认的全部生物学结合伙伴。残基仅保留所有最优序列比对共同支持的相同氨基酸位置。不对应的接触既不当丢失，也不当错误。

距离箱严格为<2、[2,3]、(3,4]、(4,5]、(5,6]、(6,8]、(8,10]、(10,15]、>15 Å，边界单元检查通过。近接触指(3,4]。逐界面近接触比例分母分别为该结构可对应的≤8 Å接触，不是并集大小，也不是原全界面接触数。保留双侧零界面，分母为零时记NA。

主分析对界面等权；池化用于描述，未进行接触级显著性检验。每PDB等权和逐PDB留一另报，不给未经支持的独立样本p值或置信区间。没有删除QC界面。

```text
  dataset  n_interfaces  n_zero_union  n_pdb  n_union  full_exp  full_af3  mapped_exp  mapped_af3  new_coordinate_exclusions  unmapped_exp  unmapped_af3
foldbench           109             0     88    40348     29090     33095       28668       30122                          0           422          2973
    safes            50             7     19    10282      8922      9883        8913        9169                          0             9           714
```
full/mapped为全界面/可对应子集的接触数；unmapped为本次身份转移无法覆盖的接触，不能把这部分解释为AF3丢失。new_coordinate_exclusions为本次新增坐标排除。完整既有纳入/排除和原因在各数据集_mapping_audit.csv。

foldbench：legacy_mapping_stable={True: 109}；local_distance_spectrum_stable={True: 109}；可靠子集中这些标记没有变化时，不能人为构造不稳定对照。原几何并列/失败记录仍留在mapping_audit中。

safes：legacy_mapping_stable={True: 50}；local_distance_spectrum_stable={True: 50}；可靠子集中这些标记没有变化时，不能人为构造不稳定对照。原几何并列/失败记录仍留在mapping_audit中。

## 主要数值：保留、离开、进入与净变化

```text
  dataset  exp_contacts  af3_contacts  exp_close  af3_close  preserved  leaving  entering  net_close  lost_beyond8  shifted_other_le8  retention  new_close_fraction  exp_close_fraction  af3_close_fraction
foldbench         28668         30122       4090       4157       2210     1880      1947         67          1116                764     0.5403              0.4684              0.1427              0.1380
    safes          8913          9169       1214       1220        917      297       303          6            66                231     0.7554              0.2484              0.1362              0.1331
```
上述比例为0–1。leaving=实验近接触−保留，entering=AF3近接触−保留，net=entering−leaving；“新近接触”包含实验4 Å以外以及≤3 Å的残基对，不能全部叫从非接触形成。lost_beyond8才表示失去8 Å接触。

### foldbench 的条件去向和来源

```text
  dataset direction other_distance_bin  pooled_fraction  n_interfaces  interface_mean  interface_median    q25    q75
foldbench   forward                 <2           0.0002           109          0.0003            0.0000 0.0000 0.0000
foldbench   forward              [2,3]           0.0443           109          0.0480            0.0278 0.0000 0.0714
foldbench   forward              (3,4]           0.5403           109          0.5174            0.6667 0.0103 0.8235
foldbench   forward              (4,5]           0.0939           109          0.0879            0.0789 0.0000 0.1429
foldbench   forward              (5,6]           0.0218           109          0.0309            0.0000 0.0000 0.0333
foldbench   forward              (6,8]           0.0267           109          0.0272            0.0000 0.0000 0.0238
foldbench   forward             (8,10]           0.0144           109          0.0137            0.0000 0.0000 0.0000
foldbench   forward            (10,15]           0.0313           109          0.0272            0.0000 0.0000 0.0000
foldbench   forward                >15           0.2271           109          0.2475            0.0000 0.0000 0.4444
foldbench   reverse                 <2           0.0000           108          0.0000            0.0000 0.0000 0.0000
foldbench   reverse              [2,3]           0.0354           108          0.0368            0.0177 0.0000 0.0670
foldbench   reverse              (3,4]           0.5316           108          0.5213            0.6667 0.0297 0.8021
foldbench   reverse              (4,5]           0.1058           108          0.1069            0.0870 0.0000 0.1603
foldbench   reverse              (5,6]           0.0207           108          0.0229            0.0000 0.0000 0.0288
foldbench   reverse              (6,8]           0.0224           108          0.0269            0.0000 0.0000 0.0391
foldbench   reverse             (8,10]           0.0101           108          0.0167            0.0000 0.0000 0.0027
foldbench   reverse            (10,15]           0.0183           108          0.0182            0.0000 0.0000 0.0000
foldbench   reverse                >15           0.2557           108          0.2503            0.0000 0.0000 0.5606
```

### safes 的条件去向和来源

```text
dataset direction other_distance_bin  pooled_fraction  n_interfaces  interface_mean  interface_median    q25    q75
  safes   forward                 <2           0.0000            42          0.0000            0.0000 0.0000 0.0000
  safes   forward              [2,3]           0.0535            42          0.0586            0.0324 0.0000 0.0902
  safes   forward              (3,4]           0.7554            42          0.6944            0.7916 0.6580 0.8712
  safes   forward              (4,5]           0.1005            42          0.0889            0.0845 0.0482 0.1250
  safes   forward              (5,6]           0.0132            42          0.0130            0.0000 0.0000 0.0162
  safes   forward              (6,8]           0.0231            42          0.0266            0.0000 0.0000 0.0109
  safes   forward             (8,10]           0.0115            42          0.0181            0.0000 0.0000 0.0000
  safes   forward            (10,15]           0.0173            42          0.0354            0.0000 0.0000 0.0000
  safes   forward                >15           0.0255            42          0.0649            0.0000 0.0000 0.0000
  safes   reverse                 <2           0.0000            42          0.0000            0.0000 0.0000 0.0000
  safes   reverse              [2,3]           0.0459            42          0.0570            0.0373 0.0000 0.0735
  safes   reverse              (3,4]           0.7516            42          0.6595            0.7721 0.5885 0.8235
  safes   reverse              (4,5]           0.1139            42          0.1200            0.1176 0.0345 0.1555
  safes   reverse              (5,6]           0.0221            42          0.0298            0.0000 0.0000 0.0378
  safes   reverse              (6,8]           0.0197            42          0.0266            0.0000 0.0000 0.0055
  safes   reverse             (8,10]           0.0107            42          0.0087            0.0000 0.0000 0.0000
  safes   reverse            (10,15]           0.0139            42          0.0241            0.0000 0.0000 0.0000
  safes   reverse                >15           0.0221            42          0.0743            0.0000 0.0000 0.0000
```

## 界面级配对比较与复合物稳健性

### foldbench

```text
  dataset scope                 metric   n    mean  median     q25     q75
foldbench   all              net_close 109  0.6147  0.0000 -5.0000  5.0000
foldbench   all   delta_close_fraction 109 -0.0039 -0.0016 -0.0213  0.0136
foldbench   all              retention 109  0.5174  0.6667  0.0103  0.8235
foldbench   all     new_close_fraction 108  0.4787  0.3333  0.1979  0.9703
foldbench   all          lost_fraction 109  0.2884  0.0000  0.0000  0.8519
foldbench   all    shift_4to6_fraction 109  0.1188  0.1127  0.0000  0.1724
foldbench   all       outward_fraction 109  0.6659  0.5946  0.5000  1.0000
foldbench   all   exp_close_mean_delta 109 10.1432  0.2058  0.0376 10.0520
foldbench   all exp_close_median_delta 109 10.0857  0.0575  0.0006 10.9935
```

可定义双侧近接触比例的界面n=109；AF3−实验均值-0.389个百分点，中位数-0.164，IQR [-2.129, +1.360]。负/零/正变化界面数分别为59/1/49。

```text
  dataset               metric  n_pdb  equal_pdb_mean  equal_interface_mean
foldbench delta_close_fraction     88         -0.0032               -0.0039
foldbench            retention     88          0.4886                0.5174
foldbench   new_close_fraction     88          0.5101                0.4787
foldbench        lost_fraction     88          0.3335                0.2884
foldbench exp_close_mean_delta     88         11.8682               10.1432
```

逐PDB留一后的界面等权比例差范围为-0.544至-0.245个百分点；PDB等权范围-0.413至-0.228。这是影响诊断范围，不是置信区间。

exp_close前5个PDB占比18.51%；排序为8SBQ:245, 9BT4:155, 8KBH:130, 8FMG:114, 8IQ5:113。

lost_beyond8前5个PDB占比47.40%；排序为8SBQ:245, 8IXW:90, 8Z4F:88, 8TW0:60, 8I82:46。

### safes

```text
dataset scope                 metric  n   mean  median     q25    q75
  safes   all              net_close 50 0.1200  0.0000 -1.0000 1.0000
  safes   all   delta_close_fraction 42 0.0019 -0.0013 -0.0092 0.0126
  safes   all              retention 42 0.6944  0.7916  0.6580 0.8712
  safes   all     new_close_fraction 42 0.3405  0.2279  0.1765 0.4115
  safes   all          lost_fraction 42 0.1184  0.0000  0.0000 0.0000
  safes   all    shift_4to6_fraction 42 0.1019  0.0969  0.0534 0.1570
  safes   all       outward_fraction 42 0.5739  0.5375  0.4402 0.6667
  safes   all   exp_close_mean_delta 42 1.7371  0.1079  0.0073 0.3968
  safes   all exp_close_median_delta 42 1.6915  0.0196 -0.0420 0.2050
```

可定义双侧近接触比例的界面n=42；AF3−实验均值+0.185个百分点，中位数-0.128，IQR [-0.920, +1.255]。负/零/正变化界面数分别为22/0/20。

```text
dataset               metric  n_pdb  equal_pdb_mean  equal_interface_mean
  safes delta_close_fraction     19          0.0016                0.0019
  safes            retention     19          0.7275                0.6944
  safes   new_close_fraction     19          0.2923                0.3405
  safes        lost_fraction     19          0.0866                0.1184
  safes exp_close_mean_delta     19          1.2934                1.7371
```

逐PDB留一后的界面等权比例差范围为-0.093至+0.369个百分点；PDB等权范围-0.044至+0.295。这是影响诊断范围，不是置信区间。

exp_close前5个PDB占比51.24%；排序为1CHO:193, 1AO7:138, 3BN9:102, 1MHP:95, 1OGA:94。

lost_beyond8前5个PDB占比98.48%；排序为1VFB:24, 2JEL:23, 3BN9:16, 1AO7:1, 1MHP:1。

```text
  dataset   n  pearson  spearman
foldbench 109   0.8325    0.5917
    safes  42   0.7298    0.6239
```

## 最近原子类型与OXT敏感性

BB严格为N、CA、C、O；其余为SC，但OXT参与的最近原子对主分析单列为OXT-involved。敏感性把OXT并入BB。该分类针对使d_min最小的原子，不是氢键、盐桥或能量打分；超过8 Å的类别也只表示最近原子身份。

下表以实验近接触类型为条件；pooled为池化，interface_mean/median仅纳入该类实验近接触存在的界面。不要把“实验类型条件下AF3的变化”误认为AF3边际类型组成。

```text
  dataset        definition         type  n_exp_close  n_interfaces  pooled_retention  pooled_shift4to6  pooled_lost  pooled_outward  pooled_type_changed  interface_mean_retention  interface_median_retention  interface_mean_shift_4to6_fraction  interface_median_shift_4to6_fraction  interface_mean_lost_fraction  interface_median_lost_fraction  interface_mean_outward_fraction  interface_median_outward_fraction  interface_mean_exp_close_mean_delta  interface_median_exp_close_mean_delta
foldbench              type        BB-BB          582            73            0.6237            0.0893       0.2113          0.5704               0.1821                    0.6346                      0.8421                              0.1025                                0.0000                        0.1828                          0.0000                           0.5683                             0.5333                               7.3239                                 0.0581
foldbench              type        BB-SC         1497           104            0.5471            0.1029       0.2806          0.6533               0.2371                    0.5247                      0.6667                              0.1206                                0.0769                        0.2942                          0.0000                           0.6869                             0.6795                               9.7314                                 0.2419
foldbench              type        SC-SC         2011           104            0.5112            0.1328       0.2849          0.6748               0.3197                    0.4732                      0.5314                              0.1284                                0.0770                        0.3115                          0.0000                           0.6886                             0.6667                              10.6295                                 0.3598
foldbench type_oxt_backbone        BB-BB          582            73            0.6237            0.0893       0.2113          0.5704               0.1787                    0.6346                      0.8421                              0.1025                                0.0000                        0.1828                          0.0000                           0.5683                             0.5333                               7.3239                                 0.0581
foldbench type_oxt_backbone        BB-SC         1497           104            0.5471            0.1029       0.2806          0.6533               0.2371                    0.5247                      0.6667                              0.1206                                0.0769                        0.2942                          0.0000                           0.6869                             0.6795                               9.7314                                 0.2419
foldbench type_oxt_backbone        SC-SC         2011           104            0.5112            0.1328       0.2849          0.6748               0.3197                    0.4732                      0.5314                              0.1284                                0.0770                        0.3115                          0.0000                           0.6886                             0.6667                              10.6295                                 0.3598
    safes              type        BB-BB          133            26            0.9023            0.0150       0.0226          0.4962               0.1128                    0.9043                      1.0000                              0.0141                                0.0000                        0.0385                          0.0000                           0.5427                             0.5000                               0.7442                                 0.0064
    safes              type        BB-SC          478            42            0.7594            0.1067       0.0649          0.5397               0.1527                    0.7028                      0.8297                              0.0903                                0.0690                        0.1351                          0.0000                           0.5661                             0.5000                               1.9033                                 0.0811
    safes              type OXT-involved            7             5            0.4286            0.1429       0.1429          0.5714               1.0000                    0.4000                      0.5000                              0.2000                                0.0000                        0.2000                          0.0000                           0.7000                             1.0000                               4.0994                                 0.0560
    safes              type        SC-SC          596            40            0.7232            0.1409       0.0520          0.5789               0.1577                    0.6595                      0.7596                              0.1337                                0.1372                        0.1115                          0.0000                           0.6044                             0.5505                               1.7225                                 0.1223
    safes type_oxt_backbone        BB-BB          133            26            0.9023            0.0150       0.0226          0.4962               0.1128                    0.9043                      1.0000                              0.0141                                0.0000                        0.0385                          0.0000                           0.5427                             0.5000                               0.7442                                 0.0064
    safes type_oxt_backbone        BB-SC          485            42            0.7546            0.1072       0.0660          0.5402               0.1546                    0.7015                      0.8297                              0.0905                                0.0667                        0.1351                          0.0000                           0.5670                             0.5086                               1.8859                                 0.0811
    safes type_oxt_backbone        SC-SC          596            40            0.7232            0.1409       0.0520          0.5789               0.1577                    0.6595                      0.7596                              0.1337                                0.1372                        0.1115                          0.0000                           0.6044                             0.5505                               1.7225                                 0.1223
```

```text
  dataset        definition    structure         type  n_close  n_contacts
foldbench              type experimental        BB-BB      582        5085
foldbench              type experimental        BB-SC     1497       15069
foldbench              type experimental        SC-SC     2011        8507
foldbench              type experimental OXT-involved        0           7
foldbench              type          af3        BB-BB      648        5646
foldbench              type          af3        BB-SC     1499       15594
foldbench              type          af3        SC-SC     2009        8858
foldbench              type          af3 OXT-involved        1          24
foldbench type_oxt_backbone experimental        BB-BB      582        5085
foldbench type_oxt_backbone experimental        BB-SC     1497       15076
foldbench type_oxt_backbone experimental        SC-SC     2011        8507
foldbench type_oxt_backbone experimental OXT-involved        0           0
foldbench type_oxt_backbone          af3        BB-BB      648        5657
foldbench type_oxt_backbone          af3        BB-SC     1500       15607
foldbench type_oxt_backbone          af3        SC-SC     2009        8858
foldbench type_oxt_backbone          af3 OXT-involved        0           0
    safes              type experimental        BB-BB      133        1502
    safes              type experimental        BB-SC      478        4706
    safes              type experimental        SC-SC      596        2686
    safes              type experimental OXT-involved        7          19
    safes              type          af3        BB-BB      145        1571
    safes              type          af3        BB-SC      488        4893
    safes              type          af3        SC-SC      586        2691
    safes              type          af3 OXT-involved        1          14
    safes type_oxt_backbone experimental        BB-BB      133        1505
    safes type_oxt_backbone experimental        BB-SC      485        4722
    safes type_oxt_backbone experimental        SC-SC      596        2686
    safes type_oxt_backbone experimental OXT-involved        0           0
    safes type_oxt_backbone          af3        BB-BB      145        1578
    safes type_oxt_backbone          af3        BB-SC      489        4900
    safes type_oxt_backbone          af3        SC-SC      586        2691
    safes type_oxt_backbone          af3 OXT-involved        0           0
```

按装配偏差分层的类型条件转移（不删除任何类别）：

```text
  dataset          scope         type  n_close  n_interfaces  pooled_retention  pooled_shift4to6  pooled_lost  mean_interface_retention  mean_interface_shift4to6  mean_interface_lost
foldbench            all        BB-BB      582            73            0.6237            0.0893       0.2113                    0.6346                    0.1025               0.1828
foldbench            all        BB-SC     1497           104            0.5471            0.1029       0.2806                    0.5247                    0.1206               0.2942
foldbench            all        SC-SC     2011           104            0.5112            0.1328       0.2849                    0.4732                    0.1284               0.3115
foldbench assembly_le10A        BB-BB      375            54            0.8320            0.0747       0.0187                    0.7901                    0.1025               0.0185
foldbench assembly_le10A        BB-SC     1096            72            0.7217            0.1369       0.0529                    0.7128                    0.1700               0.0475
foldbench assembly_le10A        SC-SC     1405            70            0.6804            0.1765       0.0541                    0.6582                    0.1794               0.0467
foldbench assembly_gt10A        BB-BB      207            19            0.2464            0.1159       0.5604                    0.1927                    0.1025               0.6495
foldbench assembly_gt10A        BB-SC      401            32            0.0698            0.0100       0.9027                    0.1015                    0.0096               0.8494
foldbench assembly_gt10A        SC-SC      606            34            0.1188            0.0314       0.8201                    0.0926                    0.0235               0.8568
    safes            all        BB-BB      133            26            0.9023            0.0150       0.0226                    0.9043                    0.0141               0.0385
    safes            all        BB-SC      478            42            0.7594            0.1067       0.0649                    0.7028                    0.0903               0.1351
    safes            all OXT-involved        7             5            0.4286            0.1429       0.1429                    0.4000                    0.2000               0.2000
    safes            all        SC-SC      596            40            0.7232            0.1409       0.0520                    0.6595                    0.1337               0.1115
    safes assembly_le10A        BB-BB      130            25            0.9231            0.0154       0.0000                    0.9405                    0.0147               0.0000
    safes assembly_le10A        BB-SC      454            39            0.7687            0.1101       0.0485                    0.7329                    0.0955               0.0942
    safes assembly_le10A OXT-involved        6             4            0.5000            0.1667       0.0000                    0.5000                    0.2500               0.0000
    safes assembly_le10A        SC-SC      565            37            0.7345            0.1434       0.0354                    0.6913                    0.1405               0.0665
    safes assembly_gt10A        BB-BB        3             1            0.0000            0.0000       1.0000                    0.0000                    0.0000               1.0000
    safes assembly_gt10A        BB-SC       24             3            0.5833            0.0417       0.3750                    0.3111                    0.0222               0.6667
    safes assembly_gt10A OXT-involved        1             1            0.0000            0.0000       1.0000                    0.0000                    0.0000               1.0000
    safes assembly_gt10A        SC-SC       31             3            0.5161            0.0968       0.3548                    0.2667                    0.0500               0.6667
```

foldbench实验近接触的最近原子类型转移，行百分比：

```text
experimental_type   BB-BB   BB-SC  OXT-involved   SC-SC
            BB-BB 81.7869 13.7457        0.3436  4.1237
            BB-SC 13.2933 76.2859        0.4676  9.9532
            SC-SC 10.6415 21.1338        0.1989 68.0259
```

safes实验近接触的最近原子类型转移，行百分比：

```text
experimental_type   BB-BB   BB-SC   SC-SC
            BB-BB 88.7218 10.5263  0.7519
            BB-SC  6.9038 84.7280  8.3682
     OXT-involved  0.0000 71.4286 28.5714
            SC-SC  0.6711 15.1007 84.2282
```

类型转移另存union、shared_contacts和experimental_close三个范围，包含逐界面原始计数。OXT敏感性不会改变距离、接触身份或总近接触缺口，只改变类别归属。

## 整体装配偏差的关联

复用已保存的全装配共同CA配准RMSD，FoldBench来自results/foldbench/chain_mapping.csv，SAFES来自local_interface_geometry/safes/residue_mapping_audit.csv中的完整chain_mapping记录。展示连续关联，并沿用FoldBench已有10 Å阈值；同一阈值在SAFES仅为透明对照分层，不新增质量排除，不据结果调参。

该RMSD混合链内部构象和相对放置差异，不等于独立测得的相对链方向误差；≤10 Å不能保证局部结构正确。SAFES缺少可信装配ID，本指标只描述输入复合物的配准，不补造生物学装配。不能把每个丢失接触都归因为局部侧链packing缺陷。

### foldbench

```text
  dataset               metric   n  pearson  spearman
foldbench delta_close_fraction 109   0.0605    0.0599
foldbench        lost_fraction 109   0.8261    0.7498
foldbench            retention 109  -0.7534   -0.6512
foldbench exp_close_mean_delta 109   0.8386    0.6656
```

```text
  dataset          scope               metric  n    mean  median     q25     q75
foldbench assembly_le10A delta_close_fraction 75 -0.0046 -0.0036 -0.0187  0.0096
foldbench assembly_le10A            retention 75  0.7009  0.7524  0.6323  0.8473
foldbench assembly_le10A        lost_fraction 75  0.0467  0.0000  0.0000  0.0000
foldbench assembly_le10A exp_close_mean_delta 75  0.7322  0.1180  0.0141  0.3007
foldbench assembly_gt10A delta_close_fraction 34 -0.0024  0.0025 -0.0343  0.0357
foldbench assembly_gt10A            retention 34  0.1126  0.0000  0.0000  0.0077
foldbench assembly_gt10A        lost_fraction 34  0.8216  1.0000  0.9126  1.0000
foldbench assembly_gt10A exp_close_mean_delta 34 30.9028 30.3110 16.5493 42.7141
```

### safes

```text
dataset               metric  n  pearson  spearman
  safes delta_close_fraction 42  -0.3531   -0.2849
  safes        lost_fraction 42   0.7536    0.5991
  safes            retention 42  -0.7053   -0.5271
  safes exp_close_mean_delta 42   0.7287    0.5641
```

```text
dataset          scope               metric  n    mean  median     q25     q75
  safes assembly_le10A delta_close_fraction 39  0.0037  0.0023 -0.0071  0.0135
  safes assembly_le10A            retention 39  0.7259  0.8000  0.6750  0.8819
  safes assembly_le10A        lost_fraction 39  0.0762  0.0000  0.0000  0.0000
  safes assembly_le10A exp_close_mean_delta 39  0.8448  0.0737 -0.0027  0.2916
  safes assembly_gt10A delta_close_fraction  3 -0.0221 -0.0149 -0.0316 -0.0089
  safes assembly_gt10A            retention  3  0.2857  0.0000  0.0000  0.4286
  safes assembly_gt10A        lost_fraction  3  0.6667  1.0000  0.5000  1.0000
  safes assembly_gt10A exp_close_mean_delta  3 13.3367 18.7776  9.4540 19.9399
```

## 图件阅读

每组均有PNG和PDF：

1. `*_distance_heatmap`：完整有限距离范围、双轴symlog，色阶为对数计数；红虚线为8 Å，灰线y=x。8 Å以外的横/纵臂表现接触身份改变。所有并集点保留，没有把非接触截断到8 Å。
2. `*_transition_matrix`：按实验距离箱行归一化。>8 Å行只包含AF3形成接触的并集成员，不能解释为所有远距非接触残基对的成接触概率。
3. `*_close_conditional_transitions`：池化实验近接触去向与AF3近接触来源；`*_interface_conditional_transitions`补充界面等权均值及中位数。
4. `*_close_transition_summary`：前三柱划分实验近接触，第四柱为AF3新近接触，第四柱不是前三柱分母的一部分。
5. `*_interface_close_fraction`：逐界面比例散点、y=x与配对差值分布，颜色为全装配CA RMSD。零分母仍在审计中，但不画为零比例。
6. `*_assembly_association`：装配RMSD与丢失比例/净比例差的连续关系；阈值线只为既有规则敏感性。
7. `*_closest_atom_types`：池化与界面等权的类型条件转移。展示三项并非完备分割（未展示所有短距去向），不可把三柱硬加为100%。

## 限制与Pedro意义

本次是可对应残基子集，不能无条件外推到原全界面或不可用/歧义病例。缺失残基、序列覆盖和最优比对歧义导致的未映射接触可能有选择性。实验方法/分辨率元数据仍不完整，现有QC不等于质量校正。界面嵌套在PDB中，PDB间也可能存在同源/重复蛋白；PDB加权不保证生物学独立。

配对结果可揭示同一复合物中的AF3相关几何变化，但Pedro没有这些直接实验对应。即使存在残余向外移动，也不能由这些基准证明Pedro差异由同一机制引起；也不能将最近原子类型直接等同于具体相互作用或能量缺陷。没有运行Pedro、没有下载或重跑AF3，既有结果保持原样。

## 验证和可复现性

分箱边界和OXT规则先做单元检查；两组各4个代表结构先验证，其中包含最大并集、最远AF3距离和典型结构。正式计算逐并集残基对复核已有距离、所有≤8 Å接触原子标签及残基类型；检查源结构SHA-256与既有审计一致。另用显式欧氏距离循环抽检。新排除为0时，不代表旧未映射部分已解决。

`run.py --validate`执行抽查；`run.py`补齐原子并生成主结果；`summarize.py`生成类型及条件摘要；`write_report.py`生成本报告数据主体。`interpretation.md`记录结合结果的六问解读并附入报告。`run.sbatch`用于计算节点；logs/保留标准输出/错误，checkpoints/按结构保存补算表。input_manifest.csv记录本次读取数据哈希。目录内结果可重新生成，不改动旧分析。

```text
  dataset  n_distances  n_contact_atoms  max_abs_error  n_ties
foldbench        80696            58790         0.0000       0
    safes        20564            18082         0.0000       0
```

## 六个科学问题的结论

### 1. AF3是否系统性地把实验3–4 Å接触向外移动？

**存在向外移动的条件趋势，但不支持“全部局部接触统一膨胀”。** 实验近接触中，AF3距离变长的池化比例为FoldBench 65.2%、SAFES 55.4%；条件距离差中位数分别+0.160和+0.039 Å。逐界面“条件距离差中位数”的中位数仅+0.058和+0.020 Å，说明典型局部变化远小于尾部造成的均值。

FoldBench全体实验近接触条件平均距离增加9.07 Å，主要包含装配相关的远距离流失，不能将9 Å称为局部packing松弛。RMSD≤10 Å子集逐界面条件距离差中位数的中位数为+0.016 Å，SAFES对应+0.018 Å。多数局部保留接触的改变较小，同时存在整体接触身份替换。

此外，先按实验距离筛选(3,4]本身会产生条件选择/向邻近箱回归的效应；一个方向的流出不能证明净缺口或因果偏差，必须同时看AF3近接触来源及总体净变化。

### 2. 流失是否被新增近接触补偿？

**在总数上基本补偿，而且略有净增。** FoldBench实验近接触4090对，AF3为4157对；保留2210，离开1880，进入1947，净增67。SAFES为1214→1220；保留917，离开297，进入303，净增6。

池化保留率为54.0%/75.5%，逐界面保留率中位数为66.7%/79.2%。AF3近接触中新形成的池化比例为46.8%/24.8%，逐界面中位数为33.3%/22.8%。因此总体分布近似可以掩盖明显的残基身份交换。

接触总数分母也变化：可对应≤8 Å接触为28668→30122和8913→9169。因此池化近接触比例小幅下降0.466和0.315个百分点，不能把比例下降读成近接触绝对数下降。主要的界面等权配对比例差则为FoldBench **−0.389个百分点**、SAFES **+0.185个百分点**；中位数分别−0.164、−0.128个百分点，两组不存在一致的大幅净下降。

每PDB等权差为−0.320和+0.164个百分点。FoldBench逐PDB留一的界面均值差始终为负（−0.544至−0.245个百分点），SAFES跨零（−0.093至+0.369），后者方向不能称稳健。这里使用的是可对应子集，与上一轮全界面比例不同；特别是SAFES的AF3有714个接触未能对应、实验仅9个，不能混用两种分母并声称结果冲突或将映射外接触都视为假接触。

### 3. 是否集中于特定最近原子类型？

涉及多类，侧链参与类别的保留率通常低于BB-BB。FoldBench池化BB-BB/BB-SC/SC-SC保留率为62.4%/54.7%/51.1%；SAFES为90.2%/75.9%/72.3%。逐界面条件均值对应63.5%/52.5%/47.3%和90.4%/70.3%/65.9%。

BB-SC实验近接触转到(4,6]的池化比例为10.3%/10.7%，界面等权均值为12.1%/9.0%。SC-SC相应池化为13.3%/14.1%，并不支持只有BB-SC产生变化。即使在RMSD≤10 Å的FoldBench子集，BB-BB/BB-SC/SC-SC池化保留率为83.2%/72.2%/68.0%，侧链参与类型仍有较多交换，但这不是侧链能量错误的证明。

实验近接触的最近原子类别改变池化比例，FoldBench三类分别18.2%/23.7%/32.0%，SAFES为11.3%/15.3%/15.8%。分类身份本身会改变，不能拿实验类型标签代替AF3最近原子类型。SAFES实验近接触仅7对涉及OXT，AF3为1对；FoldBench为0/1。OXT改归BB不改变总体距离或转移计数，类型结论变化很小；小样本OXT单列不作强解释。

### 4. 大变化是否主要关联装配偏差？

**远距离接触身份丢失与装配偏差明显关联，但净近接触比例差并非同一件事。** FoldBench 109个界面中34个RMSD>10 Å，实验近接触流失到>8 Å的界面平均比例为82.2%，其余75个为4.7%；保留率均值为11.3%与70.1%。RMSD与流失比例Spearman相关为0.750，而与净近接触比例差仅0.060。装配明显改变时，新的接触也可能补偿旧接触，所以低身份恢复不必产生大净缺口。

SAFES的50个链对全部有完整既有RMSD，其中3个>10 Å且来自同一个PDB 1VFB；这不是3个独立实验复合物。SAFES RMSD与流失比例Spearman为0.599，但高偏差分层仅单一复合物，证据有限。其66对流失近接触中，1VFB、2JEL、3BN9贡献63对（95.5%），说明少数复杂体系影响远距尾部，即使没有越过10 Å阈值也可能存在局部链对重排。

FoldBench最大贡献PDB 8SBQ贡献245/1116对流失（22.0%），前5个合计47.4%。这说明远距尾部有集中性，但并非仅一个PDB。RMSD混合链内构象与相对链放置，现有证据只能称“装配/整体构象相关”，不能确诊所有观察为装配错误或晶体接触伪影。

### 5. FoldBench与SAFES一致吗？

定性一致之处是：有身份交换、侧链参与类别保留率较低、流入近接触能补偿流出，而且整体偏差与远距流失相关。差异是FoldBench远距重排更强：实验近接触流失到>8 Å的池化比例为27.3%，SAFES为5.4%。界面近接触比例Experimental–AF3相关，Pearson为0.833/0.730，Spearman为0.592/0.624；分布/比例对应尚可不等于同一接触身份被完整保留。

两组净近接触比例的界面均值方向并不一致，均远小于先前Pedro相对实验参考的约5–7个百分点残余缺口。不能用池化两组均下降来替代界面级结果。

### 6. 能解释Pedro的近接触不足吗？

这些配对结果为“局部侧链参与接触交换、部分近接触向外移动、整体链放置改变”提供了可检验机制线索，**但没有证明它们导致Pedro的残余缺口**。一般基准的新增接触几乎补偿流出，且SAFES界面均值不下降，因此不能把本轮概括为AF3普遍缺少3–4 Å接触，更不能声称发现统一能量缺陷。

观察到的结论是基准内接触身份重排及其装配关联；合理假说是这些过程可能与Pedro低近接触比例有关；未获支持的因果断言是“同一过程已经解释Pedro分布差异”。Pedro仍存在不同蛋白、状态、样本选择和质量组成的混杂。

## 推荐下一步（本轮未执行）

1. 优先为已有可靠对应基准补充分链配准后的partner位移/方向指标，将链内构象变化与相对放置变化更直接分开；不要把现有整体RMSD当作完全分离。
2. 审计不可对应的AF3接触是否集中于缺失实验坐标、序列覆盖或比对歧义，尤其SAFES，保留原全界面与可对应子集并列。
3. 补齐实验方法/分辨率及复合物家族冗余，继续采用界面/PDB层级比较。若研究Pedro机制，优先寻找同一蛋白对的可信实验对应，再设计配对验证。
4. 当前不启动氢键、盐桥或能量缺陷归因，不因尾部或类别自动删结构。
