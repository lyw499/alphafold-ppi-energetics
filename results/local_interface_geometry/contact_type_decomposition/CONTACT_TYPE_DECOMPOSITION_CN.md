# 最近原子接触类型分解：Pedro近接触缺口

## 为什么做这一步

此前知道Pedro的(3,4] Å比例偏低，但不知道是主链接触还是侧链接触贡献了差异。本轮将所有≤8 Å残基接触分类，而非只看仍留在3–4 Å的那部分。比较三层：类型在全部接触中的份额；类型内部的距离组成；每个界面的类型贡献。三者不同。

## 输入审计：无需重解析
{
  "status": "passed",
  "counts": {
    "FoldBench AF3": 40836,
    "FoldBench Experimental": 35539,
    "SAFES AF3": 9883,
    "SAFES Experimental": 8922,
    "Mgen Experimental": 10528,
    "Pedro AF3": 6671664
  },
  "n_files": 2159,
  "total_size_bytes": 558859383,
  "missing_required_cells": 0,
  "structural_reparse_required": false
}

逐文件路径、字节数、行数、全部列名见input_file_audit.csv；缺失值逐列核对见missing_field_audit.csv。所有数据来自已有接触CSV/Parquet和界面摘要，未读取任何CIF/ZIP或重新执行基准/参考结构分析。
## 定义、边界与验证

N、CA、C、O为主链BB，其余已保存蛋白残基重原子名称为SC；BB-SC不分链方向。甘氨酸CA仍为BB。标签仅描述实现d_min的那一对原子；同一残基对仍可能有稍远的其他类型相互作用。已有算法对相同最小距离使用确定性原子顺序选择，本轮不重新解析坐标或枚举平局替代标签。

OXT若存在，按本次明确给定的四原子主链规则记SC，但它是末端原子，不能在解释中硬说为生物学侧链；所有出现的原子名/残基组合及特殊标记列在atom_name_inventory.csv。非标准残基沿同一名称规则分类并保留审计，不凭空指定化学作用。

距离箱采用 <2、[2,3]、(3,4]、(4,5]、(5,6]、(6,8] Å，恰好2 Å进入第二箱，确保原[0,3]区间完整拆分。<2是严格阈值，8包含边界。标准名程序测试覆盖N/CA/C/O、CB/CG/CD/CE/NZ/OD1/OE1/SG及两侧方向，样本表保留原字段供检查。

fraction_all = 某类型某距离箱计数 / 全部接触数；fraction_within = 同一计数 / 该类型全部接触数。例如前者低，可能是类型少，也可能是类型内部近距离少。三种fraction_all在3–4 Å的和必须等于旧总近接触比例。界面零接触时所有比例NA；非零界面缺少某类型时贡献为0，但该类型内部比例和距离中位数为NA。

## 类型总量

```
                 dataset contact_type  n_contacts  total_contacts  fraction  unclassified
           FoldBench AF3        BB-BB        8250           40836  0.202028             0
           FoldBench AF3        BB-SC       20781           40836  0.508889             0
           FoldBench AF3        SC-SC       11805           40836  0.289083             0
  FoldBench Experimental        BB-BB        6554           35539  0.184417             0
  FoldBench Experimental        BB-SC       18468           35539  0.519654             0
  FoldBench Experimental        SC-SC       10517           35539  0.295928             0
               SAFES AF3        BB-BB        1652            9883  0.167156             0
               SAFES AF3        BB-SC        5311            9883  0.537387             0
               SAFES AF3        SC-SC        2920            9883  0.295457             0
      SAFES Experimental        BB-BB        1502            8922  0.168348             0
      SAFES Experimental        BB-SC        4714            8922  0.528357             0
      SAFES Experimental        SC-SC        2706            8922  0.303295             0
       Mgen Experimental        BB-BB        1087           10528  0.103248             0
       Mgen Experimental        BB-SC        5857           10528  0.556326             0
       Mgen Experimental        SC-SC        3584           10528  0.340426             0
     Mgen excl 8PBZ 8PC1        BB-BB         755            6491  0.116315             0
     Mgen excl 8PBZ 8PC1        BB-SC        3585            6491  0.552303             0
     Mgen excl 8PBZ 8PC1        SC-SC        2151            6491  0.331382             0
Mgen excl 8PBZ 8PC0 8PC1        BB-BB         723            6054  0.119425             0
Mgen excl 8PBZ 8PC0 8PC1        BB-SC        3356            6054  0.554344             0
Mgen excl 8PBZ 8PC0 8PC1        SC-SC        1975            6054  0.326231             0
               Pedro AF3        BB-BB      695115         6671664  0.104189             0
               Pedro AF3        BB-SC     3504759         6671664  0.525320             0
               Pedro AF3        SC-SC     2471790         6671664  0.370491             0
```

每组分类总数与已有总接触数一致；unclassified为0。绝对总数受数据集大小影响，跨数据集“减少”指比例，不能把百万与万的原始计数直接比较。

## 主结果：3–4 Å贡献

```
                 dataset contact_type distance_bin  count  fraction_all  fraction_within
           FoldBench AF3        BB-BB         3to4   1004      0.024586         0.121697
           FoldBench AF3        BB-SC         3to4   1939      0.047483         0.093306
           FoldBench AF3        SC-SC         3to4   2620      0.064159         0.221940
  FoldBench Experimental        BB-BB         3to4    784      0.022060         0.119622
  FoldBench Experimental        BB-SC         3to4   1843      0.051859         0.099794
  FoldBench Experimental        SC-SC         3to4   2477      0.069698         0.235523
               SAFES AF3        BB-BB         3to4    147      0.014874         0.088983
               SAFES AF3        BB-SC         3to4    517      0.052312         0.097345
               SAFES AF3        SC-SC         3to4    630      0.063746         0.215753
      SAFES Experimental        BB-BB         3to4    133      0.014907         0.088549
      SAFES Experimental        BB-SC         3to4    478      0.053575         0.101400
      SAFES Experimental        SC-SC         3to4    604      0.067698         0.223208
       Mgen Experimental        BB-BB         3to4     28      0.002660         0.025759
       Mgen Experimental        BB-SC         3to4    525      0.049867         0.089636
       Mgen Experimental        SC-SC         3to4    579      0.054996         0.161551
     Mgen excl 8PBZ 8PC1        BB-BB         3to4     15      0.002311         0.019868
     Mgen excl 8PBZ 8PC1        BB-SC         3to4    396      0.061008         0.110460
     Mgen excl 8PBZ 8PC1        SC-SC         3to4    375      0.057772         0.174338
Mgen excl 8PBZ 8PC0 8PC1        BB-BB         3to4     13      0.002147         0.017981
Mgen excl 8PBZ 8PC0 8PC1        BB-SC         3to4    386      0.063759         0.115018
Mgen excl 8PBZ 8PC0 8PC1        SC-SC         3to4    355      0.058639         0.179747
               Pedro AF3        BB-BB         3to4  26908      0.004033         0.038710
               Pedro AF3        BB-SC         3to4 193034      0.028933         0.055078
               Pedro AF3        SC-SC         3to4 332667      0.049863         0.134585
```

Pedro减实验的贡献差（负号为Pedro较少）：

```
               reference contact_type distance_bin  pedro_minus_reference_all_fraction  pedro_minus_reference_within_fraction  composition_component  within_type_component
  FoldBench Experimental        BB-BB         3to4                           -0.018027                              -0.080911              -0.006351              -0.011676
  FoldBench Experimental        BB-SC         3to4                           -0.022925                              -0.044717               0.000439              -0.023364
  FoldBench Experimental        SC-SC         3to4                           -0.019835                              -0.100938               0.013798              -0.033634
      SAFES Experimental        BB-BB         3to4                           -0.010874                              -0.049838              -0.004082              -0.006791
      SAFES Experimental        BB-SC         3to4                           -0.024642                              -0.046322              -0.000238              -0.024404
      SAFES Experimental        SC-SC         3to4                           -0.017835                              -0.088622               0.012021              -0.029856
       Mgen Experimental        BB-BB         3to4                            0.001374                               0.012951               0.000030               0.001343
       Mgen Experimental        BB-SC         3to4                           -0.020934                              -0.034559              -0.002243              -0.018690
       Mgen Experimental        SC-SC         3to4                           -0.005134                              -0.026966               0.004452              -0.009585
     Mgen excl 8PBZ 8PC1        BB-BB         3to4                            0.001722                               0.018843              -0.000355               0.002077
     Mgen excl 8PBZ 8PC1        BB-SC         3to4                           -0.032074                              -0.055383              -0.002233              -0.029841
     Mgen excl 8PBZ 8PC1        SC-SC         3to4                           -0.007910                              -0.039752               0.006041              -0.013950
Mgen excl 8PBZ 8PC0 8PC1        BB-BB         3to4                            0.001886                               0.020730              -0.000432               0.002318
Mgen excl 8PBZ 8PC0 8PC1        BB-SC         3to4                           -0.034826                              -0.059940              -0.002468              -0.032358
Mgen excl 8PBZ 8PC0 8PC1        SC-SC         3to4                           -0.008776                              -0.045161               0.006956              -0.015732
```

composition_component与within_type_component把贡献p×q的差作对称代数分解：(pP−pE)(qP+qE)/2 与 (qP−qE)(pP+pE)/2。p为该类型占全部接触的比例，q为类型内落入该箱的比例。这是构成分解，不是因果效应；对每个参考单独解释，不把多个参考当作独立复制。

## 类型内部是否偏向更长距离

```
                 dataset contact_type distance_bin   count  fraction_all  fraction_within
           FoldBench AF3        BB-BB          lt2       8      0.000196         0.000970
           FoldBench AF3        BB-BB         2to3     322      0.007885         0.039030
           FoldBench AF3        BB-BB         3to4    1004      0.024586         0.121697
           FoldBench AF3        BB-BB         4to5     808      0.019786         0.097939
           FoldBench AF3        BB-BB         5to6    1288      0.031541         0.156121
           FoldBench AF3        BB-BB         6to8    4820      0.118033         0.584242
           FoldBench AF3        BB-SC          lt2       6      0.000147         0.000289
           FoldBench AF3        BB-SC         2to3     373      0.009134         0.017949
           FoldBench AF3        BB-SC         3to4    1939      0.047483         0.093306
           FoldBench AF3        BB-SC         4to5    1779      0.043565         0.085607
           FoldBench AF3        BB-SC         5to6    2739      0.067073         0.131803
           FoldBench AF3        BB-SC         6to8   13945      0.341488         0.671046
           FoldBench AF3        SC-SC          lt2      10      0.000245         0.000847
           FoldBench AF3        SC-SC         2to3     535      0.013101         0.045320
           FoldBench AF3        SC-SC         3to4    2620      0.064159         0.221940
           FoldBench AF3        SC-SC         4to5    1716      0.042022         0.145362
           FoldBench AF3        SC-SC         5to6    1452      0.035557         0.122999
           FoldBench AF3        SC-SC         6to8    5472      0.133999         0.463532
  FoldBench Experimental        BB-BB          lt2       0      0.000000         0.000000
  FoldBench Experimental        BB-BB         2to3     255      0.007175         0.038908
  FoldBench Experimental        BB-BB         3to4     784      0.022060         0.119622
  FoldBench Experimental        BB-BB         4to5     569      0.016011         0.086817
  FoldBench Experimental        BB-BB         5to6     995      0.027997         0.151816
  FoldBench Experimental        BB-BB         6to8    3951      0.111174         0.602838
  FoldBench Experimental        BB-SC          lt2       0      0.000000         0.000000
  FoldBench Experimental        BB-SC         2to3     338      0.009511         0.018302
  FoldBench Experimental        BB-SC         3to4    1843      0.051859         0.099794
  FoldBench Experimental        BB-SC         4to5    1411      0.039703         0.076402
  FoldBench Experimental        BB-SC         5to6    2448      0.068882         0.132554
  FoldBench Experimental        BB-SC         6to8   12428      0.349700         0.672948
  FoldBench Experimental        SC-SC          lt2       2      0.000056         0.000190
  FoldBench Experimental        SC-SC         2to3     368      0.010355         0.034991
  FoldBench Experimental        SC-SC         3to4    2477      0.069698         0.235523
  FoldBench Experimental        SC-SC         4to5    1462      0.041138         0.139013
  FoldBench Experimental        SC-SC         5to6    1128      0.031740         0.107255
  FoldBench Experimental        SC-SC         6to8    5080      0.142942         0.483027
               SAFES AF3        BB-BB          lt2       0      0.000000         0.000000
               SAFES AF3        BB-BB         2to3      46      0.004654         0.027845
               SAFES AF3        BB-BB         3to4     147      0.014874         0.088983
               SAFES AF3        BB-BB         4to5     125      0.012648         0.075666
               SAFES AF3        BB-BB         5to6     275      0.027826         0.166465
               SAFES AF3        BB-BB         6to8    1059      0.107154         0.641041
               SAFES AF3        BB-SC          lt2       1      0.000101         0.000188
               SAFES AF3        BB-SC         2to3     100      0.010118         0.018829
               SAFES AF3        BB-SC         3to4     517      0.052312         0.097345
               SAFES AF3        BB-SC         4to5     440      0.044521         0.082847
               SAFES AF3        BB-SC         5to6     692      0.070019         0.130296
               SAFES AF3        BB-SC         6to8    3561      0.360316         0.670495
               SAFES AF3        SC-SC          lt2       4      0.000405         0.001370
               SAFES AF3        SC-SC         2to3     120      0.012142         0.041096
               SAFES AF3        SC-SC         3to4     630      0.063746         0.215753
               SAFES AF3        SC-SC         4to5     415      0.041991         0.142123
               SAFES AF3        SC-SC         5to6     311      0.031468         0.106507
               SAFES AF3        SC-SC         6to8    1440      0.145705         0.493151
      SAFES Experimental        BB-BB          lt2       0      0.000000         0.000000
      SAFES Experimental        BB-BB         2to3      46      0.005156         0.030626
      SAFES Experimental        BB-BB         3to4     133      0.014907         0.088549
      SAFES Experimental        BB-BB         4to5     123      0.013786         0.081891
      SAFES Experimental        BB-BB         5to6     232      0.026003         0.154461
      SAFES Experimental        BB-BB         6to8     968      0.108496         0.644474
      SAFES Experimental        BB-SC          lt2       0      0.000000         0.000000
      SAFES Experimental        BB-SC         2to3      87      0.009751         0.018456
      SAFES Experimental        BB-SC         3to4     478      0.053575         0.101400
      SAFES Experimental        BB-SC         4to5     378      0.042367         0.080187
      SAFES Experimental        BB-SC         5to6     565      0.063327         0.119856
      SAFES Experimental        BB-SC         6to8    3206      0.359336         0.680102
      SAFES Experimental        SC-SC          lt2       1      0.000112         0.000370
      SAFES Experimental        SC-SC         2to3      99      0.011096         0.036585
      SAFES Experimental        SC-SC         3to4     604      0.067698         0.223208
      SAFES Experimental        SC-SC         4to5     379      0.042479         0.140059
      SAFES Experimental        SC-SC         5to6     318      0.035642         0.117517
      SAFES Experimental        SC-SC         6to8    1305      0.146268         0.482262
       Mgen Experimental        BB-BB          lt2       0      0.000000         0.000000
       Mgen Experimental        BB-BB         2to3       5      0.000475         0.004600
       Mgen Experimental        BB-BB         3to4      28      0.002660         0.025759
       Mgen Experimental        BB-BB         4to5      68      0.006459         0.062557
       Mgen Experimental        BB-BB         5to6     113      0.010733         0.103956
       Mgen Experimental        BB-BB         6to8     873      0.082922         0.803128
       Mgen Experimental        BB-SC          lt2       8      0.000760         0.001366
       Mgen Experimental        BB-SC         2to3     123      0.011683         0.021001
       Mgen Experimental        BB-SC         3to4     525      0.049867         0.089636
       Mgen Experimental        BB-SC         4to5     470      0.044643         0.080246
       Mgen Experimental        BB-SC         5to6     775      0.073613         0.132320
       Mgen Experimental        BB-SC         6to8    3956      0.375760         0.675431
       Mgen Experimental        SC-SC          lt2      18      0.001710         0.005022
       Mgen Experimental        SC-SC         2to3      98      0.009309         0.027344
       Mgen Experimental        SC-SC         3to4     579      0.054996         0.161551
       Mgen Experimental        SC-SC         4to5     514      0.048822         0.143415
       Mgen Experimental        SC-SC         5to6     572      0.054331         0.159598
       Mgen Experimental        SC-SC         6to8    1803      0.171258         0.503069
     Mgen excl 8PBZ 8PC1        BB-BB          lt2       0      0.000000         0.000000
     Mgen excl 8PBZ 8PC1        BB-BB         2to3       1      0.000154         0.001325
     Mgen excl 8PBZ 8PC1        BB-BB         3to4      15      0.002311         0.019868
     Mgen excl 8PBZ 8PC1        BB-BB         4to5      43      0.006625         0.056954
     Mgen excl 8PBZ 8PC1        BB-BB         5to6      72      0.011092         0.095364
     Mgen excl 8PBZ 8PC1        BB-BB         6to8     624      0.096133         0.826490
     Mgen excl 8PBZ 8PC1        BB-SC          lt2       0      0.000000         0.000000
     Mgen excl 8PBZ 8PC1        BB-SC         2to3      76      0.011709         0.021199
     Mgen excl 8PBZ 8PC1        BB-SC         3to4     396      0.061008         0.110460
     Mgen excl 8PBZ 8PC1        BB-SC         4to5     276      0.042520         0.076987
     Mgen excl 8PBZ 8PC1        BB-SC         5to6     444      0.068402         0.123849
     Mgen excl 8PBZ 8PC1        BB-SC         6to8    2393      0.368664         0.667503
     Mgen excl 8PBZ 8PC1        SC-SC          lt2       1      0.000154         0.000465
     Mgen excl 8PBZ 8PC1        SC-SC         2to3      48      0.007395         0.022315
     Mgen excl 8PBZ 8PC1        SC-SC         3to4     375      0.057772         0.174338
     Mgen excl 8PBZ 8PC1        SC-SC         4to5     309      0.047604         0.143654
     Mgen excl 8PBZ 8PC1        SC-SC         5to6     330      0.050840         0.153417
     Mgen excl 8PBZ 8PC1        SC-SC         6to8    1088      0.167617         0.505811
Mgen excl 8PBZ 8PC0 8PC1        BB-BB          lt2       0      0.000000         0.000000
Mgen excl 8PBZ 8PC0 8PC1        BB-BB         2to3       0      0.000000         0.000000
Mgen excl 8PBZ 8PC0 8PC1        BB-BB         3to4      13      0.002147         0.017981
Mgen excl 8PBZ 8PC0 8PC1        BB-BB         4to5      40      0.006607         0.055325
Mgen excl 8PBZ 8PC0 8PC1        BB-BB         5to6      67      0.011067         0.092669
Mgen excl 8PBZ 8PC0 8PC1        BB-BB         6to8     603      0.099604         0.834025
Mgen excl 8PBZ 8PC0 8PC1        BB-SC          lt2       0      0.000000         0.000000
Mgen excl 8PBZ 8PC0 8PC1        BB-SC         2to3      69      0.011397         0.020560
Mgen excl 8PBZ 8PC0 8PC1        BB-SC         3to4     386      0.063759         0.115018
Mgen excl 8PBZ 8PC0 8PC1        BB-SC         4to5     263      0.043442         0.078367
Mgen excl 8PBZ 8PC0 8PC1        BB-SC         5to6     405      0.066898         0.120679
Mgen excl 8PBZ 8PC0 8PC1        BB-SC         6to8    2233      0.368847         0.665375
Mgen excl 8PBZ 8PC0 8PC1        SC-SC          lt2       1      0.000165         0.000506
Mgen excl 8PBZ 8PC0 8PC1        SC-SC         2to3      46      0.007598         0.023291
Mgen excl 8PBZ 8PC0 8PC1        SC-SC         3to4     355      0.058639         0.179747
Mgen excl 8PBZ 8PC0 8PC1        SC-SC         4to5     284      0.046911         0.143797
Mgen excl 8PBZ 8PC0 8PC1        SC-SC         5to6     295      0.048728         0.149367
Mgen excl 8PBZ 8PC0 8PC1        SC-SC         6to8     994      0.164189         0.503291
               Pedro AF3        BB-BB          lt2    2759      0.000414         0.003969
               Pedro AF3        BB-BB         2to3    9157      0.001373         0.013173
               Pedro AF3        BB-BB         3to4   26908      0.004033         0.038710
               Pedro AF3        BB-BB         4to5   49117      0.007362         0.070660
               Pedro AF3        BB-BB         5to6   96967      0.014534         0.139498
               Pedro AF3        BB-BB         6to8  510207      0.076474         0.733989
               Pedro AF3        BB-SC          lt2   15827      0.002372         0.004516
               Pedro AF3        BB-SC         2to3   71688      0.010745         0.020454
               Pedro AF3        BB-SC         3to4  193034      0.028933         0.055078
               Pedro AF3        BB-SC         4to5  315831      0.047339         0.090115
               Pedro AF3        BB-SC         5to6  539830      0.080914         0.154028
               Pedro AF3        BB-SC         6to8 2368549      0.355016         0.675809
               Pedro AF3        SC-SC          lt2   14998      0.002248         0.006068
               Pedro AF3        SC-SC         2to3  138275      0.020726         0.055941
               Pedro AF3        SC-SC         3to4  332667      0.049863         0.134585
               Pedro AF3        SC-SC         4to5  364931      0.054699         0.147638
               Pedro AF3        SC-SC         5to6  411949      0.061746         0.166660
               Pedro AF3        SC-SC         6to8 1208970      0.181210         0.489107
```

所有距离箱都保留。较长区间富集可支持“距离组成偏向更长的最近原子距离”，不能说同一批残基接触移动了；Pedro与这些实验集合没有直接残基对应。

## 界面层面与大小

```
                 dataset contact_type            metric  n_defined     mean   median      q25      q75
           FoldBench AF3        BB-BB      median_d_min        125 6.603890 6.734972 6.481733 7.000081
           FoldBench AF3        BB-BB fraction_all_3to4        131 0.031549 0.008772 0.000000 0.022623
           FoldBench AF3        BB-SC      median_d_min        131 6.440696 6.673821 6.535986 6.842167
           FoldBench AF3        BB-SC fraction_all_3to4        131 0.046196 0.045455 0.031924 0.057128
           FoldBench AF3        SC-SC      median_d_min        131 5.629156 5.654753 5.253739 6.122124
           FoldBench AF3        SC-SC fraction_all_3to4        131 0.060806 0.060000 0.044624 0.076681
  FoldBench Experimental        BB-BB      median_d_min        128 6.610155 6.809408 6.483458 7.018787
  FoldBench Experimental        BB-BB fraction_all_3to4        132 0.034499 0.008279 0.000000 0.022318
  FoldBench Experimental        BB-SC      median_d_min        132 6.399623 6.676721 6.494337 6.817285
  FoldBench Experimental        BB-SC fraction_all_3to4        132 0.050064 0.051990 0.034589 0.062918
  FoldBench Experimental        SC-SC      median_d_min        131 5.686929 5.901957 5.187461 6.316384
  FoldBench Experimental        SC-SC fraction_all_3to4        132 0.063953 0.061927 0.045369 0.082391
       Mgen Experimental        BB-BB      median_d_min         19 7.040040 6.995986 6.962396 7.142764
       Mgen Experimental        BB-BB fraction_all_3to4         21 0.001745 0.001321 0.000000 0.002594
       Mgen Experimental        BB-SC      median_d_min         21 6.737018 6.699026 6.644268 6.827981
       Mgen Experimental        BB-SC fraction_all_3to4         21 0.055181 0.052083 0.029499 0.068182
       Mgen Experimental        SC-SC      median_d_min         21 5.955019 6.000993 5.659744 6.151820
       Mgen Experimental        SC-SC fraction_all_3to4         21 0.049617 0.052023 0.041488 0.062500
Mgen excl 8PBZ 8PC0 8PC1        BB-BB      median_d_min         13 7.085089 7.060191 6.963695 7.163374
Mgen excl 8PBZ 8PC0 8PC1        BB-BB fraction_all_3to4         13 0.001266 0.001319 0.000000 0.002320
Mgen excl 8PBZ 8PC0 8PC1        BB-SC      median_d_min         13 6.606917 6.648284 6.594834 6.727345
Mgen excl 8PBZ 8PC0 8PC1        BB-SC fraction_all_3to4         13 0.074733 0.067371 0.056180 0.068742
Mgen excl 8PBZ 8PC0 8PC1        SC-SC      median_d_min         13 5.877551 6.000993 5.659744 6.135130
Mgen excl 8PBZ 8PC0 8PC1        SC-SC fraction_all_3to4         13 0.053879 0.062005 0.046404 0.064851
     Mgen excl 8PBZ 8PC1        BB-BB      median_d_min         14 7.052153 7.024335 6.961746 7.157941
     Mgen excl 8PBZ 8PC1        BB-BB fraction_all_3to4         14 0.001502 0.001320 0.000000 0.002526
     Mgen excl 8PBZ 8PC1        BB-SC      median_d_min         14 6.624062 6.664470 6.607192 6.800233
     Mgen excl 8PBZ 8PC1        BB-SC fraction_all_3to4         14 0.071030 0.066667 0.053107 0.068707
     Mgen excl 8PBZ 8PC1        SC-SC      median_d_min         14 5.893125 6.009250 5.682020 6.125246
     Mgen excl 8PBZ 8PC1        SC-SC fraction_all_3to4         14 0.053299 0.059412 0.045926 0.064469
               Pedro AF3        BB-BB      median_d_min      35992 7.009996 7.036250 6.728782 7.326386
               Pedro AF3        BB-BB fraction_all_3to4      54109 0.002386 0.000000 0.000000 0.000000
               Pedro AF3        BB-SC      median_d_min      48375 6.753544 6.738929 6.534815 6.964036
               Pedro AF3        BB-SC fraction_all_3to4      54109 0.020877 0.018182 0.000000 0.033898
               Pedro AF3        SC-SC      median_d_min      53597 6.151128 6.106895 5.730917 6.541478
               Pedro AF3        SC-SC fraction_all_3to4      54109 0.040234 0.038911 0.000000 0.058824
               SAFES AF3        BB-BB      median_d_min         42 6.718458 6.792482 6.503422 6.948148
               SAFES AF3        BB-BB fraction_all_3to4         43 0.018110 0.008000 0.000000 0.027652
               SAFES AF3        BB-SC      median_d_min         43 6.556469 6.676835 6.450532 6.767374
               SAFES AF3        BB-SC fraction_all_3to4         43 0.052765 0.051005 0.043011 0.060120
               SAFES AF3        SC-SC      median_d_min         43 5.958189 6.043469 5.614004 6.363962
               SAFES AF3        SC-SC fraction_all_3to4         43 0.059899 0.062500 0.050000 0.072034
      SAFES Experimental        BB-BB      median_d_min         42 6.678362 6.709276 6.464079 6.966749
      SAFES Experimental        BB-BB fraction_all_3to4         42 0.019540 0.004689 0.000000 0.023671
      SAFES Experimental        BB-SC      median_d_min         42 6.583642 6.669493 6.456871 6.792423
      SAFES Experimental        BB-SC fraction_all_3to4         42 0.054100 0.050688 0.041856 0.064099
      SAFES Experimental        SC-SC      median_d_min         42 5.829506 5.895874 5.579857 6.270000
      SAFES Experimental        SC-SC fraction_all_3to4         42 0.062004 0.064289 0.049797 0.076367
```

```
    version contact_type  n_interfaces  spearman  pearson
        raw        BB-BB         54109  0.457991 0.154374
        raw        BB-SC         54109  0.610739 0.227838
        raw        SC-SC         54109  0.496857 0.147781
existing_QC        BB-BB         53645  0.453497 0.153792
existing_QC        BB-SC         53645  0.611873 0.251244
existing_QC        SC-SC         53645  0.500405 0.167733
```

```
    version  size_quartile contact_type  n_interfaces  contact_min  contact_max   median     mean      q25      q75   pooled
        raw              1        BB-BB         13600            1           16 0.000000 0.000042 0.000000 0.000000 0.000101
        raw              1        BB-SC         13600            1           16 0.000000 0.003820 0.000000 0.000000 0.007365
        raw              1        SC-SC         13600            1           16 0.000000 0.022571 0.000000 0.000000 0.030727
        raw              2        BB-BB         13523           17           69 0.000000 0.001885 0.000000 0.000000 0.002207
        raw              2        BB-SC         13523           17           69 0.015625 0.021764 0.000000 0.035714 0.023043
        raw              2        SC-SC         13523           17           69 0.038462 0.041481 0.016949 0.060000 0.042096
        raw              3        BB-BB         13465           70          165 0.000000 0.003431 0.000000 0.000000 0.003465
        raw              3        BB-SC         13465           70          165 0.026087 0.028055 0.013889 0.039062 0.028242
        raw              3        SC-SC         13465           70          165 0.044118 0.046244 0.030303 0.060000 0.046538
        raw              4        BB-BB         13521          166        18018 0.000000 0.004206 0.000000 0.005814 0.004510
        raw              4        BB-SC         13521          166        18018 0.029155 0.030000 0.021622 0.037464 0.030245
        raw              4        SC-SC         13521          166        18018 0.050209 0.050768 0.039088 0.061674 0.052236
existing_QC              1        BB-BB         13598            1           16 0.000000 0.000042 0.000000 0.000000 0.000101
existing_QC              1        BB-SC         13598            1           16 0.000000 0.003815 0.000000 0.000000 0.007355
existing_QC              1        SC-SC         13598            1           16 0.000000 0.022569 0.000000 0.000000 0.030726
existing_QC              2        BB-BB         13503           17           69 0.000000 0.001873 0.000000 0.000000 0.002192
existing_QC              2        BB-SC         13503           17           69 0.015625 0.021733 0.000000 0.035714 0.023012
existing_QC              2        SC-SC         13503           17           69 0.038462 0.041504 0.016949 0.060000 0.042119
existing_QC              3        BB-BB         13358           70          165 0.000000 0.003396 0.000000 0.000000 0.003428
existing_QC              3        BB-SC         13358           70          165 0.025974 0.027991 0.013889 0.039062 0.028175
existing_QC              3        SC-SC         13358           70          165 0.044276 0.046330 0.030349 0.060000 0.046628
existing_QC              4        BB-BB         13186          166         4315 0.000000 0.004089 0.000000 0.005714 0.004192
existing_QC              4        BB-SC         13186          166         4315 0.029022 0.029867 0.021505 0.037324 0.030065
existing_QC              4        SC-SC         13186          166         4315 0.050459 0.051056 0.039409 0.061905 0.052699
```

这里相关性的分母是整个界面的接触数，因此三类型贡献相加为总近接触比例；各类型相关系数不能相加分摊0.556。尺寸分位层中的贡献均值可显示增量主要体现在哪一类，仍不能推出界面大小导致几何变化。

## 配对基准检查

```
  dataset               scope contact_type            metric  n_pairs  median_delta  q25_delta  q75_delta  spearman  pearson
FoldBench        full_matched        BB-BB fraction_all_3to4      131      0.000000  -0.004850   0.000689  0.775468 0.962814
FoldBench        full_matched        BB-BB      median_d_min      124      0.041048  -0.129936   0.198173  0.658209 0.852121
FoldBench unambiguous_mapping        BB-BB fraction_all_3to4      109      0.000000  -0.004907   0.000456  0.787944 0.963791
FoldBench unambiguous_mapping        BB-BB      median_d_min      103      0.072962  -0.125204   0.214052  0.693856 0.868760
FoldBench        full_matched        BB-SC fraction_all_3to4      131     -0.003704  -0.015375   0.007149  0.487125 0.541884
FoldBench        full_matched        BB-SC      median_d_min      131      0.023451  -0.081137   0.135146  0.630646 0.957344
FoldBench unambiguous_mapping        BB-SC fraction_all_3to4      109     -0.003493  -0.014875   0.007092  0.508400 0.555394
FoldBench unambiguous_mapping        BB-SC      median_d_min      109      0.033562  -0.092988   0.143981  0.598916 0.958405
FoldBench        full_matched        SC-SC fraction_all_3to4      131     -0.003471  -0.015259   0.007002  0.696135 0.800112
FoldBench        full_matched        SC-SC      median_d_min      131     -0.048511  -0.427276   0.256667  0.568884 0.592292
FoldBench unambiguous_mapping        SC-SC fraction_all_3to4      109     -0.003471  -0.015302   0.005782  0.744423 0.847941
FoldBench unambiguous_mapping        SC-SC      median_d_min      109     -0.047375  -0.381285   0.249540  0.591243 0.599255
    SAFES        full_matched        BB-BB fraction_all_3to4       42      0.000000  -0.001662   0.002758  0.816037 0.934512
    SAFES        full_matched        BB-BB      median_d_min       42      0.033133  -0.067860   0.202685  0.693218 0.784849
    SAFES unambiguous_mapping        BB-BB fraction_all_3to4       42      0.000000  -0.001662   0.002758  0.816037 0.934512
    SAFES unambiguous_mapping        BB-BB      median_d_min       42      0.033133  -0.067860   0.202685  0.693218 0.784849
    SAFES        full_matched        BB-SC fraction_all_3to4       42      0.001166  -0.009637   0.009626  0.576337 0.454367
    SAFES        full_matched        BB-SC      median_d_min       42     -0.022053  -0.151268   0.049464  0.504740 0.674830
    SAFES unambiguous_mapping        BB-SC fraction_all_3to4       42      0.001166  -0.009637   0.009626  0.576337 0.454367
    SAFES unambiguous_mapping        BB-SC      median_d_min       42     -0.022053  -0.151268   0.049464  0.504740 0.674830
    SAFES        full_matched        SC-SC fraction_all_3to4       42     -0.000097  -0.017467   0.009638  0.263038 0.424103
    SAFES        full_matched        SC-SC      median_d_min       42      0.078715  -0.210479   0.369400  0.581233 0.680366
    SAFES unambiguous_mapping        SC-SC fraction_all_3to4       42     -0.000097  -0.017467   0.009638  0.263038 0.424103
    SAFES unambiguous_mapping        SC-SC      median_d_min       42      0.078715  -0.210479   0.369400  0.581233 0.680366
```

配对贡献包含已形成界面中缺少某类型的零贡献；类型中位距离只比较双方确有该类型的界面。FoldBench另外保留有序链映射明确子集。这里是配对界面统计，不声称三类型中的每一个残基对都相互对应。

## 同种敏感性与解释边界

Mgen原始、排除8PBZ/8PC1、再排除8PC0各自保留，主原始结果未替换。样本集中于少数相关装配、低分辨率模型和黏附体系，不能推广为物种全局规律。本轮未删除原始接触；所有接触型比较为描述性，不做残基级p值。最近原子侧链几何不等于完整侧链堆积能、氢键网络、盐桥能、范德华能或结合自由能。

## 核心结果：先说清楚“缺的是什么”

本轮不是去找一个预设的侧链故事。数据表明：**BB-SC 是每个实验参照下 Pedro 3–4 Å 净缺口的最大单项贡献；SC-SC 是重要的第二项；相对一般基准 BB-BB 也减少，但它不解释同种参照下的缺口。**“最大单项”不表示对所有参照都独占过半，更不表示其他类无关。

### 1. 3–4 Å 分解表怎样读？

以下单位全部是占所有≤8 Å残基接触的百分比。每行三类相加等于原分析的近接触比例；所以它直接回答总体缺口由哪一类贡献，而非只看近接触子集内部谁最多。

```
                 dataset    BB-BB    BB-SC    SC-SC        合计
  FoldBench Experimental 2.206027 5.185852 6.969808 14.361687
           FoldBench AF3 2.458615 4.748261 6.415908 13.622784
      SAFES Experimental 1.490697 5.357543 6.769783 13.618023
               SAFES AF3 1.487403 5.231205 6.374583 13.093190
       Mgen Experimental 0.265957 4.986702 5.499620 10.752280
     Mgen excl 8PBZ 8PC1 0.231089 6.100755 5.777230 12.109074
Mgen excl 8PBZ 8PC0 8PC1 0.214734 6.375950 5.863892 12.454575
               Pedro AF3 0.403318 2.893341 4.986267  8.282926
```

例如Pedro为0.403% BB-BB + 2.893% BB-SC + 4.986% SC-SC = 8.283%。SC-SC仍是Pedro现存近接触中的大部分，但“现存最多”不等于“相对实验缺得最多”。比较差值才能回答后者：

```
               reference     BB-BB     BB-SC     SC-SC
  FoldBench Experimental -1.802710 -2.292511 -1.983541
       Mgen Experimental  0.137360 -2.093361 -0.513353
Mgen excl 8PBZ 8PC0 8PC1  0.188584 -3.482609 -0.877624
     Mgen excl 8PBZ 8PC1  0.172228 -3.207414 -0.790963
      SAFES Experimental -1.087379 -2.464202 -1.783515
```

负号表示Pedro少，单位为百分点。相对FoldBench实验，BB-SC少2.29个百分点，SC-SC少1.98，BB-BB少1.80；相对SAFES实验，分别少2.46、1.78、1.09。相对去掉8PBZ/8PC0/8PC1的同种实验，BB-SC少3.48、SC-SC少0.88，而BB-BB反而多0.19，部分抵消总缺口。因此BB-SC结论跨参照成立，但BB-BB结论具有明显数据集组成依赖。

### 2. 类型整体减少了吗，还是类型内的距离不同？

Pedro的所有接触中BB-SC占52.53%，FoldBench/SAFES实验分别51.97%/52.84%，非常接近。同种敏感性约55.43%。所以主要问题不是Pedro普遍没有主链—侧链接触，而是它们落在近距离箱的比例较低。

在BB-SC类型内部，(3,4] Å比例为Pedro **5.51%**，FoldBench实验 **9.98%**、SAFES实验 **10.14%**、同种扩展敏感性 **11.50%**。这里分母只含BB-SC，与上表的全接触贡献分母不同。Pedro的BB-SC在4–5 Å为9.01%，5–6 Å为15.40%，相应实验参照在这些区间较低；6–8 Å并非对所有参照都增加。主要是4–6 Å方向的相对富集，而不是每个更长区间一律增多。

SC-SC更不能简单说整体缺失：Pedro占所有接触的37.05%，反而高于一般实验约29.59%/30.33%。但类型内3–4 Å比例只有13.46%，一般实验为23.55%/22.32%。Pedro的SC-SC在5–6 Å占16.67%，高于一般实验10.73%/11.75%；同时<3 Å部分也偏多。因此它是距离组成改变，不是纯粹所有距离统一拉长。SC-SC整体占比更高，还会部分掩盖其近距离条件比例的下降。

代数分解进一步区分“类型占比”和“类型内距离组成”。例如对FoldBench，SC-SC类型增多给近接触贡献+1.38个百分点，但类型内距离变化贡献−3.36个百分点，净为−1.98。不能把这两个数学分量解释成真实过程或能量机制。

**“分布向较长最近原子距离偏移”是跨集合的描述，不代表同一批残基对真的从3–4 Å移动到了4–6 Å。**本轮没有建立Pedro与这些实验界面的残基对应。

### 3. 配对基准是否已有相同偏差？

三类型分布并非实验与AF3完全重合，但两套配对基准的近接触差异通常小于Pedro与实验之间的差异。例如BB-SC类型内3–4 Å，FoldBench实验/AF3为9.98%/9.33%，SAFES为10.14%/9.73%，都明显高于Pedro的5.51%。SC-SC实验/AF3则为23.55%/22.19%和22.32%/21.58%，方向上也有轻微下降，不能把“AF3偏低”说成Pedro独有；Pedro的幅度更大。

配对界面也不是完美重现。BB-SC近接触全分母贡献的配对差中位数（AF3−实验），FoldBench约−0.370个百分点，SAFES约+0.117个百分点。类型内中位距离配对差分别约+0.023/−0.022 Å，较小但不代表逐残基一致。全部类型的差值、IQR、相关性和有效分母在paired_benchmark_contact_type_summary.csv；没有将接触当独立样本检验。

### 4. 同种实验支持哪一种解释？

同种实验的BB-BB近接触本来就很少，不能用一般基准的较多BB-BB来直接定义Pedro错误。对主问题更有用的是：不论保留原始同种集合还是排除已有低分辨率疑点，BB-SC仍是Pedro净缺口的最大贡献项，且敏感性排除后这一缺口更明显。SC-SC差异也仍在，但幅度较小。实验参考仍是少数相关黏附复合物，不支持推广为整个物种的类型比例。

### 5. 界面层面是否也看得到？

每个非零界面先算该类型3–4 Å接触数/全部接触数，再报告分布；缺少某类型时该贡献是真正的0，不能只挑有该类型的界面来抬高典型值。Pedro的BB-SC贡献中位数为1.82%，IQR为0–3.39%；FoldBench实验为5.20%，IQR 3.46–6.29%；SAFES实验为5.07%，IQR 4.19–6.41%；同种扩展敏感性为6.74%，IQR 5.62–6.87%。因此主现象不只是几个大界面在池化中占优势。

SC-SC贡献的界面中位数Pedro约3.89%，一般实验约6.19%/6.43%。三类型中位数不能相加当成总体中位数；只有同一界面的贡献、或相同权重的均值可以相加。

### 6. 界面大小关联由哪一类驱动？

Pedro接触数与各类型3–4 Å全分母贡献的Spearman，BB-SC约0.611，SC-SC约0.497，BB-BB约0.453。按排序关联强度，BB-SC最突出，QC通过子集结论相近。但相关系数不是可加的贡献率，不能说BB-SC“解释了61.1%”或将三系数加起来还原旧0.556。

接触数从最低到最高四分位时，各类型贡献的**界面均值**如下，单位为百分点：

```
 size_quartile    BB-BB    BB-SC    SC-SC
             1 0.004165 0.381960 2.257056
             2 0.188491 2.176434 4.148105
             3 0.343056 2.805465 4.624375
             4 0.420550 2.999978 5.076781
```

BB-SC从约0.382%到3.000%，SC-SC从2.257%到5.077%，BB-BB从0.004%到0.421%。因此按均值增量看，BB-SC与SC-SC共同贡献主要变化，SC-SC增量还略大；不能把大小效应全部归给单一类型。大/小界面的离散分母与组成差异也可能参与，仍不作因果推断。

### 7. 化学补充为什么只做有限范围？

BB-SC是主要净缺口，SC-SC对一般基准缺口也有重要贡献，因此只对这两类作附属汇总，不为BB-BB展开额外化学矩阵。BB-SC按负责最近距离的侧链一侧所属既有残基化学组划分；SC-SC直接复用已有八类接触化学。没有新增芳香专属组：芳香残基沿用项目原分组，不能从这次汇总声称芳香堆积已被检测。类别只表示残基身份，不代表某个原子的供体/受体资格、氢键角度、盐桥或能量。

## 输入路径、大小及特殊原子核对

```
        input  n_files     bytes    rows                                                                                                                                                   path
FoldBench CSV        1  24467719   76375                                       ./results/local_interface_geometry/foldbench/residue_contacts.csv
    SAFES CSV        1   5330342   18805                                           ./results/local_interface_geometry/safes/residue_contacts.csv
 Mgen Parquet        1    183405   10528 ./results/local_interface_geometry/mgen_same_species_control/mgen_experimental_contact_geometry.parquet
Pedro Parquet     2156 528877917 6671664                                             ./results/local_interface_geometry/pedro/contacts/*.parquet
```

相关字段：interface_id、d_min、residue_a_type/residue_b_type、atom_a_name/atom_b_name、atom_a_element/atom_b_element；链、结构/作业、蛋白标识及既有化学分类亦已保存。没有缺失原子字段、未知残基或氢/氘最近原子；实际残基均为20种标准氨基酸，特殊原子名仅OXT。甘氨酸出现的非N/CA/C/O名称也仅OXT，没有把CA或虚构CB当侧链。严格四名称规则下的SC包含这个已明确标记的末端例外。

将OXT改归末端主链后，任一数据集、任一类型/距离箱的全接触比例最大变化为 0.0911 个百分点；完整结果见OXT_classification_sensitivity.csv。该量远小于主要BB-SC近接触缺口，具体近接触驱动排序另作核对，不能因此将OXT称为真正侧链。

## 聚焦化学结果及其解释边界

这里仍以所有≤8 Å接触作分母，比较每个化学类别对3–4 Å比例的贡献；另列类型/化学类别内部比例，避免把组成差异误认为几何改变。

**BB-SC：疏水侧链组是各参照中最大的缺口项。**相对FoldBench、SAFES、同种全参考、同种扩展敏感性，分别少1.135、1.862、0.987、1.536个百分点。这一组内部的3–4 Å比例也分别低4.66、7.35、4.63、8.01个百分点。因此不仅是类别数量不同，条件距离组成也有差别。其他组贡献随参照变化：同种扩展敏感性中，极性不带电组少0.887个百分点，正电组少0.625个百分点。此处“疏水”来自残基身份，最近原子不一定形成某种专属疏水相互作用。

**SC-SC：具体化学驱动具有参考集合依赖。**相对FoldBench，疏水–疏水贡献少0.933个百分点，为最大项；相对SAFES，最大项却是疏水–极性，少1.014个百分点。同种扩展敏感性中，特殊/其他、疏水–极性、极性–极性分别少0.629、0.497、0.490个百分点，而疏水–疏水全分母贡献略高0.048个百分点，即使该组内部3–4 Å比例低14.72个百分点。组成与条件距离可以相互抵消，不能只看一个分母。

因此，侧链化学值得继续作为解释线索，但当前最稳健的线索是BB-SC的疏水侧链接近方式，而不是所有参考都证明SC-SC疏水接触缺失。正负电残基类别不等于盐桥；本轮也没有氢键角度、供受体资格或芳香环方向信息。化学结果为池化探索描述，不能把数百万残基接触当作独立样本。

## 对七个问题的直接回答与下一步建议

1. **缺口主要来自哪类？** BB-SC在每个实验参照中均是最大的单项净缺口；SC-SC也重要。相对一般基准BB-BB减少，相对同种参考则不减少。
2. **是整体减少还是距离改变？** BB-SC整体占比与一般实验相近，主要是类型内3–4 Å较少、4–6 Å相对较多；不能概括为所有4–8 Å区间都增加。SC-SC整体更多，但近距离条件比例更低。跨集合比较不意味着同一残基对实际移动。
3. **配对基准AF3能否重现？** FoldBench/SAFES大体重现BB/SC几何，仍有轻微近接触偏低及界面差异；其幅度通常小于Pedro对实验的差异，因此不能把Pedro现象简单视为所有AF3共有偏差。
4. **同种实验是否支持？** 支持BB-SC缺口及其低分辨率敏感性稳健性，但不支持一般基准所暗示的BB-BB缺口。同种集合小且存在复合物重复和组成偏倚，不能代表全物种。
5. **大小关联由谁驱动？** BB-SC的排序相关最强；四分位均值增量则BB-SC与SC-SC共同贡献，SC-SC略大。不能据相关系数声称唯一驱动或因果解释比例。
6. **有必要进一步研究侧链化学吗？** 有，聚焦BB-SC疏水侧链组最有根据，同时保留SC-SC的参照依赖；本轮有限补充已给出，尚不足以推断能量机制。
7. **下一步优先什么？** 优先整体packing/interface-size effects：先在接触数可比的界面中比较，并控制结构/复合物重复与既有QC，判断BB-SC差异在大小匹配后是否仍存在。随后才考虑聚焦BB-SC疏水侧链接近/堆积及SC-SC；氢键几何需要供受体与角度等新证据，盐桥和芳香相互作用目前没有足够一致的优先依据。以上仅为建议，未自动开展。

## 最终核对与交付

所有原始接触均完成三类归属，未分类比例为0；每个界面的接触总数和原3–4 Å比例已对账。OXT改归主链后，逐一核对全部实验参照，最大缺口类别仍为BB-SC。零界面保留在界面表中，无法定义的比例和距离保持缺失值；任何QC/分辨率子集均不替代原始结果。

本轮只读取既有接触表与元数据，没有重新读取Pedro结构归档、下载结构或修改既有分析结果。最终检查记录见final_validation.json和OXT_driver_validation.json；输入清单见input_file_audit.csv。图01–09、界面汇总、配对基准、大小关联及化学表全部存于本目录，图08/09是附属化学结果。
