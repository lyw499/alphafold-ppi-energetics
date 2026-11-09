# 控制界面大小后，Pedro 是否仍有局部接触几何差异？

## 本轮回答什么

此前池化残基接触中，Pedro的(3,4] Å比例较低。本轮把界面作为主要观察单位：每个实验界面先得到自身比例，再与接触数可比的Pedro界面比较。目标是判断大小组成能否充分解释差异，而不是建立因果匹配研究或评估能量函数。

## 输入审计与零界面

主输入是上一轮interface_contact_type_summary.parquet，包含每界面的三类计数、类型总占比和全部六个距离箱。Pedro/Mgen的两项次要大小指标直接读既有界面表；FoldBench/SAFES按既有接触CSV中的残基索引去重得到界面残基数，再除以既有两链长度之和。没有重新计算距离或读取结构。


```
                 dataset  n_interfaces  n_zero  n_nonzero  zero_percent
           FoldBench AF3           132       1        131        0.7576
  FoldBench Experimental           132       0        132        0.0000
       Mgen Experimental            23       2         21        8.6957
Mgen excl 8PBZ 8PC0 8PC1            15       2         13       13.3333
     Mgen excl 8PBZ 8PC1            16       2         14       12.5000
               Pedro AF3        170644  116535      54109       68.2913
               SAFES AF3            50       7         43       14.0000
      SAFES Experimental            50       8         42       16.0000
```

所有必要计数、总体近接触比例及三项大小指标齐全。某非零界面没有某接触类型时，该类型的全接触贡献为0，类型内部比例为NA；这是定义域缺失而非坏数据。零界面保留在analysis_interface_table.parquet和审计表，比例分析只使用非零界面。未把FoldBench不可用AF3案例当作零。

输入文件及哈希清单：


```
                                                                                                                                               path    bytes
./results/local_interface_geometry/contact_type_decomposition/interface_contact_type_summary.parquet  9445486
                                ./results/local_interface_geometry/foldbench/interface_summaries.csv   171469
                                   ./results/local_interface_geometry/foldbench/residue_contacts.csv 24467719
                                    ./results/local_interface_geometry/safes/interface_summaries.csv    74498
                                       ./results/local_interface_geometry/safes/residue_contacts.csv  5330342
                   ./results/local_interface_geometry/pedro/pedro_interface_geometry_summary.parquet 19904528
./results/local_interface_geometry/mgen_same_species_control/mgen_experimental_interface_summary.csv    29537
```

## 分母、权重和精确方法

- 总近接触比例：界面内(3,4] Å残基对数/全部≤8 Å残基对数。
- 类型贡献：界面内该类型(3,4] Å残基对数/全部≤8 Å残基对数；三类贡献逐界面相加为总比例。
- 类型内部比例：界面内该类型某距离箱的计数/该类型所有接触数；只在类型存在时定义。
- 主结果对界面等权；匹配Pedro使用下述权重。堆叠图采用均值，因此贡献可相加。中位数/IQR另存表，不能把三类中位数相加。
- 全距离范围沿用<2、[2,3]、(3,4]、(4,5]、(5,6]、(6,8] Å六箱。图07/08显示全范围分箱概率，而非假定箱内连续密度或精细峰形。
- BB名称仍为N、CA、C、O；OXT沿上一轮主定义作为明确的末端例外。没有改变类型或距离定义。

每套实验参考分别匹配：对每个实验界面，在Pedro非零界面中寻找log1p(contact_count)绝对距离最小的候选；所有并列最近候选平均分配该目标的1份权重，可被不同实验界面重复使用。这样不根据几何结果挑选候选，也避免在相同接触数界面中随意挑一个。总体权重和等于支持范围内实验界面数。k=20敏感性扩大到第20近距离以内的所有并列候选。

超出Pedro大小范围的实验界面单独记录，不作边界强行外推。每套参考重新匹配。权重有效样本数ESS=(Σw)²/Σw²仅描述权重集中度，不是独立生物学样本量；界面仍嵌套在结构/复合物内。分位数统一采用加权经验CDF反函数，故偶数样本的中位数可能不同于相邻两点平均的惯例。

## 一、共同绝对大小分层及支持度

Pedro非零接触数四分位边界为16、69、165，统一应用为1–16、17–69、70–165、≥166。最大层较宽，所以分层仅用于直观检查，主要大小控制依据更精确的最近邻匹配。n<5只标为低支持，不删除观察，不作为统计显著性门槛。


```
                 dataset    stratum  n_interfaces  median_contact_count  median    q25    q75
           FoldBench AF3    (0, 16]             2                8.5000  0.0000 0.0000 0.0833
           FoldBench AF3   (16, 69]             6               63.0000  0.1212 0.1000 0.1562
           FoldBench AF3  (69, 165]            14              134.0000  0.1146 0.0915 0.1863
           FoldBench AF3 (165, inf]           109              282.0000  0.1289 0.1106 0.1536
  FoldBench Experimental    (0, 16]             3               10.0000  0.1000 0.0833 0.1667
  FoldBench Experimental   (16, 69]             7               41.0000  0.1176 0.0732 0.1250
  FoldBench Experimental  (69, 165]            37              121.0000  0.1519 0.1221 0.2735
  FoldBench Experimental (165, inf]            85              278.0000  0.1303 0.1208 0.1481
               SAFES AF3    (0, 16]             1                5.0000  0.0000 0.0000 0.0000
               SAFES AF3   (16, 69]             7               48.0000  0.1364 0.1020 0.2292
               SAFES AF3  (69, 165]            13              125.0000  0.1151 0.0909 0.1368
               SAFES AF3 (165, inf]            22              289.5000  0.1228 0.1198 0.1452
      SAFES Experimental    (0, 16]             0                   NaN     NaN    NaN    NaN
      SAFES Experimental   (16, 69]             7               43.0000  0.1277 0.0794 0.2000
      SAFES Experimental  (69, 165]            14              120.0000  0.1190 0.1111 0.1538
      SAFES Experimental (165, inf]            21              248.0000  0.1303 0.1232 0.1520
       Mgen Experimental    (0, 16]             0                   NaN     NaN    NaN    NaN
       Mgen Experimental   (16, 69]             1               50.0000  0.0000 0.0000 0.0000
       Mgen Experimental  (69, 165]             7               88.0000  0.1250 0.1059 0.1573
       Mgen Experimental (165, inf]            13              760.0000  0.1016 0.0874 0.1319
     Mgen excl 8PBZ 8PC1    (0, 16]             0                   NaN     NaN    NaN    NaN
     Mgen excl 8PBZ 8PC1   (16, 69]             0                   NaN     NaN    NaN    NaN
     Mgen excl 8PBZ 8PC1  (69, 165]             6               86.5000  0.1250 0.1146 0.1573
     Mgen excl 8PBZ 8PC1 (165, inf]             8              759.0000  0.1319 0.0986 0.1329
Mgen excl 8PBZ 8PC0 8PC1    (0, 16]             0                   NaN     NaN    NaN    NaN
Mgen excl 8PBZ 8PC0 8PC1   (16, 69]             0                   NaN     NaN    NaN    NaN
Mgen excl 8PBZ 8PC0 8PC1  (69, 165]             6               86.5000  0.1250 0.1146 0.1573
Mgen excl 8PBZ 8PC0 8PC1 (165, inf]             7              760.0000  0.1319 0.1040 0.1347
               Pedro AF3    (0, 16]         13600                4.0000  0.0000 0.0000 0.0000
               Pedro AF3   (16, 69]         13523               39.0000  0.0606 0.0370 0.0909
               Pedro AF3  (69, 165]         13465              110.0000  0.0769 0.0595 0.0947
               Pedro AF3 (165, inf]         13521              263.0000  0.0847 0.0721 0.0969
```

上表比例为0–1；图02–04显示百分比。FoldBench最小层仅3个实验界面，SAFES最小层没有实验界面，Mgen的17–69层仅1个。它们不能提供强分层证据。完整三类及类型内部指标均保存在size_stratified_contact_type_summary.csv。

## 二、匹配是否真的改善了大小平衡

先看大小，后看几何。下表n是有正权重的独特Pedro界面数或实验界面数；effective_n是上述权重ESS。标准化差用log1p接触数的均值差除以匹配前实验/Pedro合并标准差，前后使用同一分母。weighted_ecdf_max_gap为两条加权经验CDF的最大差，不是检验p值。


```
               reference                  phase     n  effective_n   median      q25      q75  standardized_mean_difference  weighted_ecdf_max_gap  n_experimental_outside_support
  FoldBench Experimental       Experimental_all   132     132.0000 230.0000 128.0000 299.0000                        0.0000                 0.0000                               0
  FoldBench Experimental Experimental_supported   132     132.0000 230.0000 128.0000 299.0000                        0.0000                 0.0000                               0
  FoldBench Experimental           Pedro_before 54109   54109.0000  69.0000  16.0000 165.0000                       -1.1128                 0.5157                               0
  FoldBench Experimental            Pedro_after 10189    1511.5077 230.0000 128.0000 299.0000                        0.0000                 0.0076                               0
      SAFES Experimental       Experimental_all    42      42.0000 162.0000  82.0000 248.0000                        0.0000                 0.0000                               0
      SAFES Experimental Experimental_supported    42      42.0000 162.0000  82.0000 248.0000                        0.0000                 0.0000                               0
      SAFES Experimental           Pedro_before 54109   54109.0000  69.0000  16.0000 165.0000                       -0.9054                 0.3561                               0
      SAFES Experimental            Pedro_after  3879    1083.2863 162.0000  82.0000 248.0000                        0.0000                 0.0000                               0
       Mgen Experimental       Experimental_all    21      21.0000 699.0000  89.0000 771.0000                        0.0000                 0.0000                               0
       Mgen Experimental Experimental_supported    21      21.0000 699.0000  89.0000 771.0000                        0.0000                 0.0000                               0
       Mgen Experimental           Pedro_before 54109   54109.0000  69.0000  16.0000 165.0000                       -1.3816                 0.5730                               0
       Mgen Experimental            Pedro_after  1391      47.9663 699.0000  89.0000 772.0000                        0.0001                 0.0952                               0
     Mgen excl 8PBZ 8PC1       Experimental_all    14      14.0000 437.0000  88.0000 760.0000                        0.0000                 0.0000                               0
     Mgen excl 8PBZ 8PC1 Experimental_supported    14      14.0000 437.0000  88.0000 760.0000                        0.0000                 0.0000                               0
     Mgen excl 8PBZ 8PC1           Pedro_before 54109   54109.0000  69.0000  16.0000 165.0000                       -1.3208                 0.5416                               0
     Mgen excl 8PBZ 8PC1            Pedro_after   990      36.0741 437.0000  88.0000 760.0000                        0.0001                 0.1429                               0
Mgen excl 8PBZ 8PC0 8PC1       Experimental_all    13      13.0000 757.0000  88.0000 760.0000                        0.0000                 0.0000                               0
Mgen excl 8PBZ 8PC0 8PC1 Experimental_supported    13      13.0000 757.0000  88.0000 760.0000                        0.0000                 0.0000                               0
Mgen excl 8PBZ 8PC0 8PC1           Pedro_before 54109   54109.0000  69.0000  16.0000 165.0000                       -1.2850                 0.5416                               0
Mgen excl 8PBZ 8PC0 8PC1            Pedro_after   974      31.4666 756.0000  88.0000 760.0000                        0.0001                 0.1538                               0
```

## 三、匹配前后的总近接触差异

下表均为界面均值的百分比/百分点；负差表示Pedro较低。gap_attenuation_percent=(1−匹配后差/匹配前差)×100，只是描述性差值衰减，不能解释为大小的因果解释率。


```
               reference  experimental_mean  pedro_before_mean  pedro_matched_mean  before_difference  matched_difference  gap_attenuation_percent
  FoldBench Experimental            14.8516             6.3497              8.1665            -8.5019             -6.6851                  21.3691
      SAFES Experimental            13.5644             6.3497              7.9125            -7.2147             -5.6518                  21.6617
       Mgen Experimental            10.6543             6.3497              8.6416            -4.3046             -2.0127                  53.2436
     Mgen excl 8PBZ 8PC1            12.5831             6.3497              8.5454            -6.2334             -4.0377                  35.2244
Mgen excl 8PBZ 8PC0 8PC1            12.9878             6.3497              8.5712            -6.6380             -4.4165                  33.4665
```

这些数值不应直接与上一轮“全部残基接触池化比例”相减：池化隐含按接触数给大界面更大权重，本轮首先改成界面等权，再在同一估计口径内比较匹配前后。

## 四、匹配后哪一种类型贡献最多

Pedro减实验的三类界面均值差，单位百分点：


```
               reference  BB-BB_close  BB-SC_close  SC-SC_close
  FoldBench Experimental      -3.0815      -2.1222      -1.4814
       Mgen Experimental       0.2401      -2.5317       0.2789
Mgen excl 8PBZ 8PC0 8PC1       0.2787      -4.5513      -0.1439
     Mgen excl 8PBZ 8PC1       0.2567      -4.1935      -0.1009
      SAFES Experimental      -1.6175      -2.6106      -1.4238
```

逐参考最大的单项净缺口为：FoldBench Experimental：BB-BB；SAFES Experimental：BB-SC；Mgen Experimental：BB-SC；Mgen excl 8PBZ 8PC1：BB-SC；Mgen excl 8PBZ 8PC0 8PC1：BB-SC。

必须允许主导类别随参考及统计权重改变，不能预设BB-SC一定最大。此处是匹配后的界面均值贡献，上一轮主结论来自池化接触；两者区别应通过同口径的匹配前后表判断，不能把所有排序变化归因于大小匹配。

## 五、类型数量还是类型内部距离组成

share是类型占界面全部接触的比例；within是该类型内部距离比例。以下列出实验、匹配Pedro及两者差，全部换成百分比/百分点。


```
               reference            metric  experimental_mean  pedro_matched_mean  matched_difference
  FoldBench Experimental       BB-SC_share            49.6090             52.4271              2.8181
  FoldBench Experimental       SC-SC_share            29.1344             37.4055              8.2711
  FoldBench Experimental BB-SC_within_3to4            12.0222              5.4553             -6.5669
  FoldBench Experimental BB-SC_within_4to5            10.8588              8.9390             -1.9198
  FoldBench Experimental BB-SC_within_5to6            12.1831             15.3733              3.1902
  FoldBench Experimental BB-SC_within_6to8            62.9118             67.7262              4.8143
  FoldBench Experimental SC-SC_within_3to4            23.9773             13.5238            -10.4535
  FoldBench Experimental SC-SC_within_4to5            14.2610             14.6812              0.4202
  FoldBench Experimental SC-SC_within_5to6            11.0991             16.7011              5.6021
  FoldBench Experimental SC-SC_within_6to8            46.4522             48.6320              2.1798
      SAFES Experimental       BB-SC_share            51.4022             52.2363              0.8341
      SAFES Experimental       SC-SC_share            29.4838             37.7298              8.2460
      SAFES Experimental BB-SC_within_3to4            10.6951              5.2935             -5.4016
      SAFES Experimental BB-SC_within_4to5             8.4464              8.8132              0.3668
      SAFES Experimental BB-SC_within_5to6            13.3772             15.4154              2.0382
      SAFES Experimental BB-SC_within_6to8            65.4932             68.0875              2.5943
      SAFES Experimental SC-SC_within_3to4            21.2687             13.0608             -8.2079
      SAFES Experimental SC-SC_within_4to5            13.4871             14.5446              1.0575
      SAFES Experimental SC-SC_within_5to6            11.2862             16.7559              5.4697
      SAFES Experimental SC-SC_within_6to8            50.0046             49.5370             -0.4676
       Mgen Experimental       BB-SC_share            53.6839             53.8735              0.1896
       Mgen Experimental       SC-SC_share            37.2273             36.0121             -1.2152
       Mgen Experimental BB-SC_within_3to4             9.5301              5.5480             -3.9821
       Mgen Experimental BB-SC_within_4to5             7.8228              8.5232              0.7004
       Mgen Experimental BB-SC_within_5to6            10.8560             14.9558              4.0998
       Mgen Experimental BB-SC_within_6to8            69.9780             68.6847             -1.2934
       Mgen Experimental SC-SC_within_3to4            14.9714             14.8096             -0.1617
       Mgen Experimental SC-SC_within_4to5            15.0403             16.0183              0.9780
       Mgen Experimental SC-SC_within_5to6            16.4655             16.1831             -0.2824
       Mgen Experimental SC-SC_within_6to8            50.2578             46.6156             -3.6422
     Mgen excl 8PBZ 8PC1       BB-SC_share            57.6418             54.0719             -3.5700
     Mgen excl 8PBZ 8PC1       SC-SC_share            31.6740             36.2994              4.6254
     Mgen excl 8PBZ 8PC1 BB-SC_within_3to4            12.2457              5.3900             -6.8557
     Mgen excl 8PBZ 8PC1 BB-SC_within_4to5             8.6320              8.4530             -0.1790
     Mgen excl 8PBZ 8PC1 BB-SC_within_5to6            11.0379             14.9994              3.9615
     Mgen excl 8PBZ 8PC1 BB-SC_within_6to8            66.2454             68.9172              2.6719
     Mgen excl 8PBZ 8PC1 SC-SC_within_3to4            16.6715             14.7311             -1.9404
     Mgen excl 8PBZ 8PC1 SC-SC_within_4to5            15.9711             16.0663              0.0953
     Mgen excl 8PBZ 8PC1 SC-SC_within_5to6            16.0303             16.4870              0.4567
     Mgen excl 8PBZ 8PC1 SC-SC_within_6to8            48.4272             46.0782             -2.3490
Mgen excl 8PBZ 8PC0 8PC1       BB-SC_share            58.0448             53.9263             -4.1185
Mgen excl 8PBZ 8PC0 8PC1       SC-SC_share            31.0124             36.5085              5.4961
Mgen excl 8PBZ 8PC0 8PC1 BB-SC_within_3to4            12.8517              5.4256             -7.4261
Mgen excl 8PBZ 8PC0 8PC1 BB-SC_within_4to5             8.8594              8.4489             -0.4105
Mgen excl 8PBZ 8PC0 8PC1 BB-SC_within_5to6            10.5769             14.9913              4.4144
Mgen excl 8PBZ 8PC0 8PC1 BB-SC_within_6to8            65.9666             68.9319              2.9652
Mgen excl 8PBZ 8PC0 8PC1 SC-SC_within_3to4            17.0798             14.6912             -2.3886
Mgen excl 8PBZ 8PC0 8PC1 SC-SC_within_4to5            16.1069             15.9759             -0.1311
Mgen excl 8PBZ 8PC0 8PC1 SC-SC_within_5to6            15.7337             16.4542              0.7205
Mgen excl 8PBZ 8PC0 8PC1 SC-SC_within_6to8            48.0440             46.3030             -1.7410
```

图07/08及size_matched_distance_bin_summary.csv保留0–8 Å六箱，含短于3 Å部分。主列mean是等界面的类型内部比例均值；pooled_fraction另外给出加权计数池化比例，分母不同，不混用。类型不存在的界面不进入within均值，相关n及ESS逐指标保存。任何“较长区间富集”均是跨集合组成描述，不表示对应残基对真实移动。

## 六、连续大小关联和既有QC

Spearman在全部界面上计算，log1p为单调变换，故排序相关与原接触数相同；图09的分箱仅用于显示。


```
    version size_variable            metric     n  spearman
        raw    n_contacts       total_close 54109    0.5562
        raw    n_contacts       BB-BB_close 54109    0.4580
        raw    n_contacts       BB-SC_close 54109    0.6107
        raw    n_contacts       SC-SC_close 54109    0.4969
        raw    n_contacts BB-BB_within_3to4 35992    0.3505
        raw    n_contacts BB-SC_within_3to4 48375    0.4850
        raw    n_contacts SC-SC_within_3to4 53597    0.5644
existing_QC    n_contacts       total_close 53645    0.5580
existing_QC    n_contacts       BB-BB_close 53645    0.4535
existing_QC    n_contacts       BB-SC_close 53645    0.6119
existing_QC    n_contacts       SC-SC_close 53645    0.5004
existing_QC    n_contacts BB-BB_within_3to4 35535    0.3464
existing_QC    n_contacts BB-SC_within_3to4 47911    0.4860
existing_QC    n_contacts SC-SC_within_3to4 53133    0.5672
```

相关系数不可相加，也不是解释方差。BB-SC/SC-SC各自的全分母贡献与条件比例都报告；不将一个较大的相关系数解释成唯一机制。

## 七、敏感性：不同大小指标、最近邻数及原有QC

各次敏感性重新匹配，原始结果始终保留。总界面残基数使用log1p，界面比例直接使用原比例距离。单独匹配一个次要指标不等于同时平衡三个变量；不能替代接触数主分析。QC仅沿用existing_severe_overlap，无新阈值。


```
               reference            size_variable pedro_version  nearest_k  before_difference  matched_difference
  FoldBench Experimental               n_contacts           raw          1            -0.0850             -0.0669
      SAFES Experimental               n_contacts           raw          1            -0.0721             -0.0565
       Mgen Experimental               n_contacts           raw          1            -0.0430             -0.0201
     Mgen excl 8PBZ 8PC1               n_contacts           raw          1            -0.0623             -0.0404
Mgen excl 8PBZ 8PC0 8PC1               n_contacts           raw          1            -0.0664             -0.0442
  FoldBench Experimental total_interface_residues           raw          1            -0.0850             -0.0695
  FoldBench Experimental       interface_fraction           raw          1            -0.0632             -0.0371
      SAFES Experimental total_interface_residues           raw          1            -0.0721             -0.0587
      SAFES Experimental       interface_fraction           raw          1            -0.0721             -0.0395
       Mgen Experimental total_interface_residues           raw          1            -0.0430             -0.0215
       Mgen Experimental       interface_fraction           raw          1            -0.0430             -0.0320
     Mgen excl 8PBZ 8PC1 total_interface_residues           raw          1            -0.0623             -0.0430
     Mgen excl 8PBZ 8PC1       interface_fraction           raw          1            -0.0623             -0.0492
Mgen excl 8PBZ 8PC0 8PC1 total_interface_residues           raw          1            -0.0664             -0.0477
Mgen excl 8PBZ 8PC0 8PC1       interface_fraction           raw          1            -0.0664             -0.0540
  FoldBench Experimental               n_contacts   existing_QC          1            -0.0852             -0.0669
      SAFES Experimental               n_contacts   existing_QC          1            -0.0723             -0.0565
       Mgen Experimental               n_contacts   existing_QC          1            -0.0432             -0.0195
     Mgen excl 8PBZ 8PC1               n_contacts   existing_QC          1            -0.0625             -0.0403
Mgen excl 8PBZ 8PC0 8PC1               n_contacts   existing_QC          1            -0.0665             -0.0441
  FoldBench Experimental               n_contacts           raw         20            -0.0850             -0.0669
      SAFES Experimental               n_contacts           raw         20            -0.0721             -0.0557
       Mgen Experimental               n_contacts           raw         20            -0.0430             -0.0201
     Mgen excl 8PBZ 8PC1               n_contacts           raw         20            -0.0623             -0.0398
Mgen excl 8PBZ 8PC0 8PC1               n_contacts           raw         20            -0.0664             -0.0438
```

上表差值为0–1比例；乘100得百分点。次要指标的平衡、域外名单、样本数同样见matching_balance_summary.csv和outside_size_support.csv。界面比例是近连续变量，最近邻可能只有一个，需结合ESS看待，不能与接触数大量精确并列匹配具有同等稳定性。

## 八、FoldBench和SAFES保留配对的控制

使用相同interface_id配对，双方非零且该指标可定义才计算差值；零接触仍在配对ID审计。保留全映射和明确链映射子集。为了描述配对大小变化的影响，拟合Δfraction=截距+斜率×Δlog1p(contact_count)，截距代表该线性模型在大小不变时的预测差，不能称为实验控制或因果效应；不报告p值。


```
  dataset       scope      metric   n    mean  median     q25    q75  size_adjusted_intercept  spearman
FoldBench        full total_close 131 -0.0098 -0.0085 -0.0278 0.0097                  -0.0126    0.5486
FoldBench        full BB-SC_close 131 -0.0042 -0.0037 -0.0158 0.0072                  -0.0048    0.4871
FoldBench        full SC-SC_close 131 -0.0036 -0.0035 -0.0153 0.0073                  -0.0053    0.6961
FoldBench unambiguous total_close 109 -0.0114 -0.0064 -0.0298 0.0092                  -0.0138    0.5758
FoldBench unambiguous BB-SC_close 109 -0.0040 -0.0035 -0.0149 0.0071                  -0.0046    0.5084
FoldBench unambiguous SC-SC_close 109 -0.0048 -0.0035 -0.0153 0.0058                  -0.0063    0.7444
    SAFES        full total_close  42 -0.0018 -0.0017 -0.0118 0.0108                  -0.0003    0.6697
    SAFES        full BB-SC_close  42 -0.0001  0.0000 -0.0101 0.0097                   0.0030    0.5763
    SAFES        full SC-SC_close  42 -0.0007 -0.0002 -0.0184 0.0105                  -0.0023    0.2630
    SAFES unambiguous total_close  42 -0.0018 -0.0017 -0.0118 0.0108                  -0.0003    0.6697
    SAFES unambiguous BB-SC_close  42 -0.0001  0.0000 -0.0101 0.0097                   0.0030    0.5763
    SAFES unambiguous SC-SC_close  42 -0.0007 -0.0002 -0.0184 0.0105                  -0.0023    0.2630
```

另外对FoldBench AF3/SAFES AF3各自的接触数分布匹配Pedro，结果保存在主结果表中，作为补充而非替代原始配对。所有控制差值均按界面定义，不把残基接触当独立样本。

## 九、Mgen敏感性与限制

沿用全23界面（21非零）、排除8PBZ/8PC1的16界面（14非零）、再排除8PC0的15界面（13非零）。各自独立匹配，没有新增过滤门槛。Mgen包含同一装配多个链对及相关黏附复合物，因此只作描述性支持，不能把13–21个界面理解成同样数量的独立实验复合物。

只匹配接触数不平衡蛋白种类、构象状态、长度、实验方法、序列、复合物选择或可信度。接触数本身也是几何派生量；观察到的“匹配后差异”不是证明某种物理机制。最近原子BB-SC几何更不等同于氢键几何。

## 文件与验证

输入计数与旧比例逐界面对账、每目标匹配权重和、零界面NA、输入文件前后SHA-256均通过。没有读取CIF/ZIP、重新计算原子距离、下载数据或进行新的化学相互作用分析。结构文件读取数为0。

主表：size_stratified_contact_type_summary.csv、matching_balance_summary.csv、size_matched_close_contact_summary.csv、size_matched_distance_bin_summary.csv、pedro_size_contact_type_correlations.csv。补充：size_matched_effect_differences.csv、matched_interface_ids.csv、outside_size_support.csv、benchmark_paired_size_context.csv、benchmark_pair_ids.csv、data_audit.csv、analysis_interface_table.parquet、input_manifest.csv、methods.json、validation.json。图01–09位于figures/。

## 匹配支持度的进一步核对

接触数主分析没有实验界面超出Pedro支持范围。FoldBench 132个实验界面中131个可精确匹配接触数，余下1个最大只差5个接触（0.318%）；SAFES 42个全部精确匹配。Mgen全参考21个中14个精确匹配，剩余最大只差1个（0.132%）；两套敏感性最大也只差1个。

因此，Mgen的小样本经验CDF最大差可达0.095–0.154，不能脱离距离尺度误解为严重不平衡：一些重复接触数的观测整体相差1，就会形成较高的CDF台阶差。逐目标距离误差及log尺度均值差均很小。完整记录见matching_pair_size_gaps.csv；Mgen匹配权重ESS约48、36、31，仍不能改变实验端只有21、14、13个非独立界面的限制。

**次要界面比例匹配存在支持限制。**Pedro非零界面比例最大约0.9495；FoldBench有21个实验界面超过此范围。它们只在该项敏感性比较中单列为域外，保留于主接触数比较和原始数据。此敏感性只代表111个支持范围内实验界面，不能把它与全部132界面的主分析差异直接归因于调整变量。逐个ID见outside_size_support.csv。其他主/次匹配没有域外观察。

最终审计修正了3个加权分位数单元格的浮点累积边界误差（CDF恰好0.25/0.5/0.75时的数值精度问题），均值、差值、权重、图及科学结论均未改变，详见weighted_quantile_numerical_audit.json。分位数搜索统一采用1e−12数值容差。

## 十、九个问题的明确回答

### 1. 界面大小分布有多大不同？

差别很明显。非零Pedro接触数中位数为69，IQR 16–165；FoldBench实验为230（128–299），SAFES为162（82–248），Mgen全参考为699（89–771）。这里统一用经验CDF分位数。匹配前log接触数标准化均值差约−1.11、−0.91、−1.38；匹配后绝对值均不超过0.00014。Pedro包含很多更小的界面，因此直接比较原始界面比例会混入明显的大小组成差异。

### 2. 相近大小内，总3–4 Å缺口还存在吗？

存在。按界面等权均值，Pedro原始为6.350%；匹配FoldBench后升到8.167%，实验为14.852%；匹配SAFES后为7.913%，实验为13.564%。匹配后仍分别差−6.685、−5.652个百分点，原差距描述性缩小21.4%、21.7%。Mgen全参考匹配后Pedro为8.642%，实验10.654%，仍差−2.013个百分点，差距缩小53.2%。这些是界面均值而非接触池化比例。

共同分层提供一致的直观支持：例如70–165层，Pedro总比例均值约7.773%，FoldBench/SAFES实验约19.179%/13.522%；≥166层约8.497%对13.405%/13.632%。宽层仍可能残留大小差，所以精确匹配比单靠这些箱更重要。最小样本层不作强解释。

### 3. BB-SC仍是最大的单一贡献吗？

**不能再不加限定地说“所有参考都是”。**匹配后界面均值贡献中，SAFES和Mgen及其两套敏感性仍是BB-SC最大；但FoldBench是BB-BB（−3.082个百分点），随后BB-SC（−2.122）和SC-SC（−1.481）。

这个排序并非单由大小匹配造成：匹配前若同样采用界面等权，FoldBench的BB-BB缺口已为−3.211个百分点，大于BB-SC的−2.919。上一轮“所有参考BB-SC最大”描述的是接触池化结果，不能直接搬到另一种权重下。FoldBench BB-BB实验贡献的均值3.450%高于中位数0.822%，显示异质性；按中位数看，BB-SC的实验/匹配Pedro为5.182%/2.778%，BB-BB为0.822%/0%。不同统计量回答不同问题，不能相加中位数来分摊总体缺口，也不能据此删除某类界面。

### 4. BB-SC是否仍然近距离少、中间距离多？

支持“3–4 Å较少、4–6 Å合计相对较多”，但必须细分。匹配后BB-SC类型内部3–4 Å比例，Pedro对FoldBench/SAFES/Mgen分别为5.455%对12.022%、5.294%对10.695%、5.548%对9.530%。两套Mgen敏感性差距更大。

类型内部4–6 Å合计的Pedro减实验差依次约+1.270、+2.405、+4.800个百分点，两套Mgen敏感性约+3.783、+4.004。**较稳健的富集来自5–6 Å；4–5 Å并非都增加**，例如FoldBench反而少1.920个百分点。因此不应泛称每个中间距离箱都增加，更不能说所有残基对都向外移动。

类型总量并未普遍减少：匹配Pedro的BB-SC平均份额约52.43%、52.24%、53.87%，相应全参考实验为49.61%、51.40%、53.68%。Mgen两套敏感性实验为57.64%/58.04%，对应Pedro约54.07%/53.93%，这两组确有一定总量差，但仍同时存在更大的类型内距离差。结论应是几何组成差异不能仅用BB-SC总量解释，而非断言每套参考都完全没有总量效应。

### 5. SC-SC差异还剩多少？

对FoldBench/SAFES，匹配后SC-SC近接触全分母贡献仍低1.481/1.424个百分点，类型内部3–4 Å比例仍低10.453/8.208个百分点。

对Mgen全参考，SC-SC全分母贡献从原先低0.938个百分点变成略高0.279个百分点，类型内部近接触比例只差−0.162个百分点，基本相近。两套Mgen敏感性中，全分母贡献仅低0.101/0.144个百分点，但类型内部仍低1.940/2.389个百分点。后者部分被Pedro的SC-SC类型份额较高所抵消。因此“SC-SC匹配后仍同样明显偏低”不适用于所有实验参考。

### 6. 大小与近接触比例的关系由谁共同体现？

Pedro总比例与log1p接触数Spearman为0.556；三类全分母贡献分别为BB-BB 0.458、BB-SC 0.611、SC-SC 0.497。BB-SC的全分母关联最强，但BB-SC类型内部近接触比例相关为0.485，SC-SC类型内部为0.564，后者更强。类型存在与否、类型数量及类型内几何共同参与，不能将相关系数分摊成贡献率。

匹配FoldBench大小分布时，Pedro均值提升约1.817个百分点，其中BB-SC提升0.796、SC-SC提升0.891、BB-BB仅0.130个百分点；SAFES也是两类共同贡献。因此大小组成的调整同时影响BB-SC和SC-SC，而非单一类别。

### 7. 基准AF3也一样吗？

有较小的差异，但幅度与Pedro不同。保留真实配对后，FoldBench AF3减实验的总比例平均差约−0.983个百分点，SAFES约−0.176；配对Δlog接触数线性描述调整后的截距分别−1.258和−0.032个百分点。

BB-SC配对调整截距FoldBench为−0.476个百分点、SAFES为+0.299；SC-SC为−0.530、−0.233，均明显小于Pedro对应的约−2.12/−2.61与−1.48/−1.42。配对模型只是简单描述，不能保证所有非线性大小效应已消除。独立把Pedro匹配到两套基准AF3自身大小分布，仍分别有约−5.493/−5.238个百分点总差，支持Pedro并非仅重现一个普遍且同样幅度的AF3偏差。明确映射子集另表保留。

### 8. Mgen既有敏感性支持吗？

支持总缺口和BB-SC差异保留。排除8PBZ/8PC1后，实验与匹配Pedro均值为12.583%对8.545%，差−4.038个百分点；再排除8PC0后为12.988%对8.571%，差−4.417。大小匹配使原差距缩小约35.2%/33.5%，未消除。SC-SC全分母差几乎消失则是需要保留的新信息。样本小、重复装配和黏附体系偏倚仍限制推广。

### 9. 最终更支持A、B还是C？

**更支持B：界面大小组成与一部分观察差距有关，但在接触数可比的界面中仍存在局部接触几何差异。**它不是纯大小组成解释（A），也不是大小几乎无关（C）。这个判断对总近接触比例及BB-SC较稳健，对SC-SC需区分一般基准与Mgen。

既有严重重叠QC剔除敏感性并未改变方向：匹配后总差FoldBench约−6.686、SAFES−5.647、Mgen全参考−1.952个百分点；k=20敏感性也相近。总界面残基数敏感性保留差异；界面比例敏感性亦有差异，但FoldBench的21个域外界面必须单独注明，不能当成全样本严格复制。

## 下一步建议（尚未执行）

最有依据的下一步是，在保持接触数匹配的同时，先检查结果对**复合物/蛋白对重复权重及界面可信度**是否稳健，并针对不同参考的差异分别解释。若进一步做化学几何，优先选择BB-SC最近原子涉及的残基与原子类别，沿此前疏水侧链线索作分层描述；不要直接把它解释成氢键不足。氢键需要供受体和角度证据，盐桥、芳香及能量分析当前没有足够理由自动启动。也不应假定Mgen还存在与一般基准同样的SC-SC缺口。

本轮未进行上述新分析，未改动上一轮接触类型分解及其他已有结果。
