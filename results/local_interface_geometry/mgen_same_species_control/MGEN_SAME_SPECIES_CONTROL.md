# M. genitalium 同种实验对照

仅重新计算已下载并审计的 29 个参考链对；全部与旧接触数、链长和界面残基数精确一致。主结果只含 23 个 Mgen–Mgen 链对；其他物种单独保存。Pedro 和 FoldBench/SAFES 仅读取既有结果，输入文件哈希核对未变。

定义：重原子最小距离 ≤8 Å 的唯一残基对。复用原解析器、标签链 ID、首个坐标模型、默认替代构象；排除 H/D，标准蛋白残基处理不变。分箱为 [0,3]、(3,4]、(4,5]、(5,6]、(6,8] Å。零接触界面保留，条件距离与比例未定义记 NA。

## 主要池化比较

```
               dataset  n_interfaces  n_zero  n_contacts     mean   median      q25      q75  fraction_0_3  fraction_3_4  fraction_4_5  fraction_5_6  fraction_6_8  fraction_lt_0_5  fraction_lt_1_0  fraction_lt_1_5  fraction_lt_2_0
FoldBench Experimental           132       0       35539 6.022782 6.477215 4.798746 7.283099      0.027097      0.143617      0.096851      0.128619      0.603816         0.000000         0.000000         0.000028         0.000056
         FoldBench AF3           132       1       40836 5.998795 6.444984 4.779731 7.283448      0.030708      0.136228      0.105373      0.134171      0.593520         0.000000         0.000000         0.000049         0.000588
    SAFES Experimental            50       8        8922 6.063343 6.525052 4.894464 7.316783      0.026115      0.136180      0.098633      0.124972      0.614100         0.000000         0.000000         0.000000         0.000112
             SAFES AF3            50       7        9883 6.059536 6.525698 4.919326 7.313104      0.027421      0.130932      0.099160      0.129313      0.613174         0.000000         0.000000         0.000000         0.000506
             Pedro AF3        170644  116535     6671664 6.116708 6.507373 5.153185 7.347011      0.037877      0.082829      0.109400      0.157194      0.612700         0.000094         0.000627         0.001777         0.005034
     Mgen Experimental            23       2       10528 6.155029 6.574119 5.145430 7.362613      0.023936      0.107523      0.099924      0.138678      0.629939         0.000000         0.000475         0.001140         0.002470
```

## 敏感性分析

```
                  view      weighting  n_pairs  n_zero  n_nonzero  n_pdb  n_assemblies  n_contacts     mean   median      q25      q75  fraction_0_3  fraction_3_4  fraction_4_5  fraction_5_6  fraction_6_8  fraction_lt_0_5  fraction_lt_1_0  fraction_lt_1_5  fraction_lt_2_0
              all_Mgen         pooled       23       2         21      7             8       10528 6.155029 6.574119 5.145430 7.362613      0.023936      0.107523      0.099924      0.138678      0.629939              0.0         0.000475         0.001140         0.002470
              all_Mgen assembly_equal       23       2         21      7             8       10528 6.164225 6.576383 5.157130 7.367774      0.022142      0.107568      0.099254      0.142314      0.628723              0.0         0.000200         0.000480         0.001146
              all_Mgen      PDB_equal       23       2         21      7             8       10528 6.167600 6.571387 5.170346 7.368535      0.022585      0.103845      0.100225      0.146000      0.627345              0.0         0.000229         0.000549         0.001310
     exclude_8PBZ_8PC1         pooled       16       2         14      5             6        6491 6.153090 6.587886 5.108087 7.362618      0.019411      0.121091      0.096749      0.130334      0.632414              0.0         0.000000         0.000000         0.000154
     exclude_8PBZ_8PC1 assembly_equal       16       2         14      5             6        6491 6.166195 6.587932 5.150441 7.368338      0.019835      0.112977      0.097985      0.138947      0.630256              0.0         0.000000         0.000000         0.000193
     exclude_8PBZ_8PC1      PDB_equal       16       2         14      5             6        6491 6.171313 6.587886 5.170346 7.370796      0.019994      0.108847      0.099091      0.143434      0.628634              0.0         0.000000         0.000000         0.000232
exclude_8PBZ_8PC0_8PC1         pooled       15       2         13      4             5        6054 6.145547 6.584865 5.073512 7.359121      0.019161      0.124546      0.096961      0.126693      0.632640              0.0         0.000000         0.000000         0.000165
exclude_8PBZ_8PC0_8PC1 assembly_equal       15       2         13      4             5        6054 6.147919 6.581430 5.096199 7.361212      0.019225      0.120927      0.098818      0.130580      0.630449              0.0         0.000000         0.000000         0.000232
exclude_8PBZ_8PC0_8PC1      PDB_equal       15       2         13      4             5        6054 6.149747 6.576383 5.113592 7.363854      0.019271      0.117752      0.100409      0.134098      0.628470              0.0         0.000000         0.000000         0.000290
```

## 逐装配结果

```
  structure_id pdb_id  resolution_A  pairs  nonzero_pairs  n_contacts     mean   median      q25      q75  fraction_0_3  fraction_3_4  fraction_4_5  fraction_5_6  fraction_6_8  fraction_lt_0_5  fraction_lt_1_0  fraction_lt_1_5  fraction_lt_2_0
6rut-assembly1   6RUT          2.65      6              5        1793 6.142714 6.599538 5.066610 7.350869      0.020078      0.132181      0.090351      0.118238      0.639152              0.0         0.000000         0.000000         0.000000
6rut-assembly2   6RUT          2.65      6              5        1777 6.138514 6.582714 5.023557 7.352069      0.018008      0.135059      0.094541      0.114800      0.637591              0.0         0.000000         0.000000         0.000000
6yrk-assembly1   6YRK          4.10      1              1         757 6.157106 6.614237 5.026549 7.369454      0.019815      0.134742      0.093791      0.108322      0.643329              0.0         0.000000         0.000000         0.000000
8pbx-assembly1   8PBX          3.30      1              1         865 6.167744 6.547642 5.199561 7.379631      0.017341      0.104046      0.106358      0.150289      0.621965              0.0         0.000000         0.000000         0.000000
8pby-assembly1   8PBY          3.70      1              1         862 6.133515 6.522932 5.157391 7.333260      0.020882      0.098608      0.109049      0.161253      0.610209              0.0         0.000000         0.000000         0.001160
8pbz-assembly1   8PBZ         11.00      6              6        3122 6.158007 6.558447 5.212550 7.356702      0.032992      0.081038      0.106662      0.151826      0.627482              0.0         0.001602         0.003844         0.008008
8pc0-assembly1   8PC0         17.00      1              1         437 6.257579 6.599545 5.338417 7.443522      0.022883      0.073227      0.093822      0.180778      0.629291              0.0         0.000000         0.000000         0.000000
8pc1-assembly1   8PC1         18.00      1              1         915 6.158624 6.541344 5.143187 7.372523      0.025137      0.101639      0.099454      0.153005      0.620765              0.0         0.000000         0.000000         0.000000
```

## 逐界面分布

```
               dataset                   metric  n_defined       mean     median        q25        q75
FoldBench Experimental             median_d_min        132   6.320220   6.518340   6.395318   6.656518
FoldBench Experimental             fraction_0_3        132   0.030295   0.024691   0.016475   0.032728
FoldBench Experimental             fraction_3_4        132   0.148516   0.131010   0.118228   0.158104
FoldBench Experimental             fraction_4_5        132   0.105538   0.093966   0.079702   0.117343
FoldBench Experimental             fraction_5_6        132   0.127750   0.127292   0.111926   0.143486
FoldBench Experimental             fraction_6_8        132   0.587901   0.614065   0.589563   0.647900
FoldBench Experimental          fraction_lt_2_0        132   0.000218   0.000000   0.000000   0.000000
FoldBench Experimental               n_contacts        132 269.234848 231.500000 128.000000 299.000000
         FoldBench AF3             median_d_min        131   6.340375   6.509907   6.390940   6.611433
         FoldBench AF3             fraction_0_3        131   0.032139   0.029126   0.019324   0.039922
         FoldBench AF3             fraction_3_4        131   0.138551   0.126316   0.109028   0.155526
         FoldBench AF3             fraction_4_5        131   0.110389   0.099492   0.087583   0.116003
         FoldBench AF3             fraction_5_6        131   0.136370   0.134328   0.114062   0.152366
         FoldBench AF3             fraction_6_8        131   0.582550   0.612676   0.583767   0.636810
         FoldBench AF3          fraction_lt_2_0        131   0.000783   0.000000   0.000000   0.000000
         FoldBench AF3               n_contacts        132 309.363636 251.000000 178.000000 362.500000
    SAFES Experimental             median_d_min         42   6.451973   6.536874   6.314652   6.619101
    SAFES Experimental             fraction_0_3         42   0.029552   0.029196   0.020041   0.036704
    SAFES Experimental             fraction_3_4         42   0.135644   0.128346   0.115320   0.154315
    SAFES Experimental             fraction_4_5         42   0.098518   0.094883   0.079733   0.113679
    SAFES Experimental             fraction_5_6         42   0.131361   0.126532   0.113898   0.139026
    SAFES Experimental             fraction_6_8         42   0.604925   0.615528   0.577590   0.643272
    SAFES Experimental          fraction_lt_2_0         42   0.000031   0.000000   0.000000   0.000000
    SAFES Experimental               n_contacts         50 178.440000 129.500000  44.000000 245.250000
             SAFES AF3             median_d_min         43   6.491454   6.547957   6.386001   6.604333
             SAFES AF3             fraction_0_3         43   0.027663   0.027778   0.017189   0.036932
             SAFES AF3             fraction_3_4         43   0.130773   0.122807   0.114980   0.147581
             SAFES AF3             fraction_4_5         43   0.092647   0.093525   0.075553   0.108007
             SAFES AF3             fraction_5_6         43   0.136186   0.127986   0.113639   0.156188
             SAFES AF3             fraction_6_8         43   0.612731   0.620000   0.577219   0.638195
             SAFES AF3          fraction_lt_2_0         43   0.000608   0.000000   0.000000   0.000000
             SAFES AF3               n_contacts         50 197.660000 139.500000  51.750000 267.500000
             Pedro AF3             median_d_min      54109   6.617567   6.559880   6.411094   6.771893
             Pedro AF3             fraction_0_3      54109   0.027596   0.026415   0.000000   0.043321
             Pedro AF3             fraction_3_4      54109   0.063497   0.068966   0.026316   0.091463
             Pedro AF3             fraction_4_5      54109   0.095409   0.101266   0.064516   0.125000
             Pedro AF3             fraction_5_6      54109   0.151538   0.152941   0.117647   0.181818
             Pedro AF3             fraction_6_8      54109   0.661960   0.629560   0.592593   0.690476
             Pedro AF3          fraction_lt_2_0      54109   0.003272   0.000000   0.000000   0.002899
             Pedro AF3               n_contacts     170644  39.096974   0.000000   0.000000  11.000000
             Pedro AF3 total_interface_residues     170644  18.773353   0.000000   0.000000  12.000000
             Pedro AF3       interface_fraction     170644   0.025308   0.000000   0.000000   0.015819
     Mgen Experimental             median_d_min         21   6.608968   6.585218   6.541344   6.616879
     Mgen Experimental             fraction_0_3         21   0.021081   0.019815   0.017341   0.025137
     Mgen Experimental             fraction_3_4         21   0.106543   0.105882   0.087391   0.132895
     Mgen Experimental             fraction_4_5         21   0.101845   0.099454   0.092105   0.116519
     Mgen Experimental             fraction_5_6         21   0.131518   0.140118   0.112360   0.154639
     Mgen Experimental             fraction_6_8         21   0.639014   0.632959   0.620765   0.642480
     Mgen Experimental          fraction_lt_2_0         21   0.001722   0.000000   0.000000   0.000000
     Mgen Experimental               n_contacts         23 457.739130 678.000000  86.500000 765.500000
     Mgen Experimental total_interface_residues         23 194.652174 277.000000  46.000000 297.000000
     Mgen Experimental       interface_fraction         23   0.088860   0.124213   0.022222   0.144017
```

## 短距离诊断

```
               dataset  threshold  contacts  contact_fraction  affected_interfaces  fraction_all_interfaces  fraction_nonzero_interfaces
FoldBench Experimental        0.5         0          0.000000                    0                 0.000000                     0.000000
FoldBench Experimental        1.0         0          0.000000                    0                 0.000000                     0.000000
FoldBench Experimental        1.5         1          0.000028                    1                 0.007576                     0.007576
FoldBench Experimental        2.0         2          0.000056                    2                 0.015152                     0.015152
         FoldBench AF3        0.5         0          0.000000                    0                 0.000000                     0.000000
         FoldBench AF3        1.0         0          0.000000                    0                 0.000000                     0.000000
         FoldBench AF3        1.5         2          0.000049                    2                 0.015152                     0.015267
         FoldBench AF3        2.0        24          0.000588                   15                 0.113636                     0.114504
    SAFES Experimental        0.5         0          0.000000                    0                 0.000000                     0.000000
    SAFES Experimental        1.0         0          0.000000                    0                 0.000000                     0.000000
    SAFES Experimental        1.5         0          0.000000                    0                 0.000000                     0.000000
    SAFES Experimental        2.0         1          0.000112                    1                 0.020000                     0.023810
             SAFES AF3        0.5         0          0.000000                    0                 0.000000                     0.000000
             SAFES AF3        1.0         0          0.000000                    0                 0.000000                     0.000000
             SAFES AF3        1.5         0          0.000000                    0                 0.000000                     0.000000
             SAFES AF3        2.0         5          0.000506                    5                 0.100000                     0.116279
             Pedro AF3        0.5       629          0.000094                  464                 0.002719                     0.008575
             Pedro AF3        1.0      4184          0.000627                 2383                 0.013965                     0.044041
             Pedro AF3        1.5     11853          0.001777                 6126                 0.035899                     0.113216
             Pedro AF3        2.0     33584          0.005034                14665                 0.085939                     0.271027
     Mgen Experimental        0.5         0          0.000000                    0                 0.000000                     0.000000
     Mgen Experimental        1.0         5          0.000475                    2                 0.086957                     0.095238
     Mgen Experimental        1.5        12          0.001140                    4                 0.173913                     0.190476
     Mgen Experimental        2.0        26          0.002470                    5                 0.217391                     0.238095
```

## 化学组成（探索性）

```
               dataset           contact_class  pooled_fraction  mean_interface_fraction
FoldBench Experimental hydrophobic_hydrophobic         0.188159                 0.169949
FoldBench Experimental         opposite_charge         0.039759                 0.045689
FoldBench Experimental             same_charge         0.035623                 0.040982
FoldBench Experimental             polar_polar         0.046315                 0.049579
FoldBench Experimental       hydrophobic_polar         0.155013                 0.145797
FoldBench Experimental     hydrophobic_charged         0.182588                 0.176229
FoldBench Experimental           polar_charged         0.095641                 0.093954
FoldBench Experimental        special_or_other         0.256901                 0.277821
         FoldBench AF3 hydrophobic_hydrophobic         0.194412                 0.175411
         FoldBench AF3         opposite_charge         0.038838                 0.047861
         FoldBench AF3             same_charge         0.032643                 0.037871
         FoldBench AF3             polar_polar         0.047140                 0.047198
         FoldBench AF3       hydrophobic_polar         0.149280                 0.141574
         FoldBench AF3     hydrophobic_charged         0.183025                 0.179646
         FoldBench AF3           polar_charged         0.090778                 0.090018
         FoldBench AF3        special_or_other         0.263885                 0.280420
    SAFES Experimental hydrophobic_hydrophobic         0.158485                 0.143471
    SAFES Experimental         opposite_charge         0.022304                 0.021166
    SAFES Experimental             same_charge         0.014571                 0.014329
    SAFES Experimental             polar_polar         0.078682                 0.073295
    SAFES Experimental       hydrophobic_polar         0.211948                 0.211367
    SAFES Experimental     hydrophobic_charged         0.130240                 0.131266
    SAFES Experimental           polar_charged         0.097848                 0.098343
    SAFES Experimental        special_or_other         0.285922                 0.306764
             SAFES AF3 hydrophobic_hydrophobic         0.152788                 0.137905
             SAFES AF3         opposite_charge         0.025397                 0.024322
             SAFES AF3             same_charge         0.015683                 0.023461
             SAFES AF3             polar_polar         0.078822                 0.069511
             SAFES AF3       hydrophobic_polar         0.212790                 0.205791
             SAFES AF3     hydrophobic_charged         0.134676                 0.132611
             SAFES AF3           polar_charged         0.099767                 0.107236
             SAFES AF3        special_or_other         0.280077                 0.299165
             Pedro AF3 hydrophobic_hydrophobic         0.182764                 0.131890
             Pedro AF3         opposite_charge         0.048591                 0.071408
             Pedro AF3             same_charge         0.038495                 0.065949
             Pedro AF3             polar_polar         0.067736                 0.067129
             Pedro AF3       hydrophobic_polar         0.184538                 0.152189
             Pedro AF3     hydrophobic_charged         0.183381                 0.189003
             Pedro AF3           polar_charged         0.144937                 0.180176
             Pedro AF3        special_or_other         0.149558                 0.142256
     Mgen Experimental hydrophobic_hydrophobic         0.146562                 0.198842
     Mgen Experimental         opposite_charge         0.027356                 0.018215
     Mgen Experimental             same_charge         0.018332                 0.014045
     Mgen Experimental             polar_polar         0.120916                 0.113560
     Mgen Experimental       hydrophobic_polar         0.238032                 0.234742
     Mgen Experimental     hydrophobic_charged         0.127185                 0.107897
     Mgen Experimental           polar_charged         0.133359                 0.103424
     Mgen Experimental        special_or_other         0.188260                 0.209274
```

## 敏感性方法及限制

装配等权先在每个生物学装配内汇总接触，再使每个非零装配的总权重相同；PDB 等权以同一规则合并同一 PDB 的装配。没有删除或任意挑选链对。零接触装配没有条件距离分布，不赋予虚构距离。装配等权仍不能消除不同 PDB 来自同一蛋白复合物/研究的相关性。

低分辨率敏感性明确排除既有审计中的 8PBZ（11 Å，短距离几何疑点）和 8PC1（18 Å）；扩展视图另排除同一审计中的 8PC0（17 Å）。这不是新增分辨率阈值或更改原始 QC。原始结果全部保留。

组间不配对，不做残基级显著性检验。ECDF 显示至多 10,000 个经验秩点；直方图为 0.25 Å 同一分箱。化学分类复用已有宽泛类别。短距离不自动称作碰撞。

## 八个问题的直接回答

1. **同种实验是否类似 Pedro？** 两者中位数接近，但局部距离分箱不同。不能将相近中位数解释为相同几何。下表同时给出峰区和较远区比例。

比例单位为 %；中位数单位为 Å。

```
               dataset   median  fraction_3_4  fraction_4_5  fraction_5_6  fraction_6_8
FoldBench Experimental 6.477215     14.361687      9.685135     12.861926     60.381553
         FoldBench AF3 6.444984     13.622784     10.537271     13.417083     59.352042
    SAFES Experimental 6.525052     13.618023      9.863259     12.497198     61.409998
             SAFES AF3 6.525698     13.093190      9.916017     12.931296     61.317414
             Pedro AF3 6.507373      8.282926     10.939984     15.719407     61.269962
     Mgen Experimental 6.574119     10.752280      9.992401     13.867781     62.993921
```

2. **是否也有弱 3–4 Å 峰？** 原始 Mgen 实验 (3,4] Å 比例 10.75%，低于 FoldBench/SAFES 实验的 14.36%/13.62%，但高于 Pedro 的 8.28%。因此有部分相似性，程度并不相同。

3. **是否更像一般实验基准？** 针对 (3,4] Å：去掉预先标记的 8PBZ/8PC1 后为 12.11%，再去掉 8PC0 后为 12.45%，接近一般实验而非 Pedro。但不同分箱、完整曲线可能给出不同相似程度，不应据此把参考集定义为通用标准。

4. **重复观测的影响？** 以下同时列出装配等权与 PDB 等权。等权消除一个装配仅因链对/接触数多而主导，但仍保留同一复合物跨 PDB 的相关性；不等于获得独立生物学样本。

```
                  view      weighting  n_pairs  n_nonzero  n_pdb  n_assemblies  fraction_3_4  fraction_4_5  fraction_5_6  fraction_6_8
              all_Mgen         pooled       23         21      7             8     10.752280      9.992401     13.867781     62.993921
              all_Mgen assembly_equal       23         21      7             8     10.756750      9.925357     14.231391     62.872310
              all_Mgen      PDB_equal       23         21      7             8     10.384481     10.022469     14.600001     62.734521
     exclude_8PBZ_8PC1         pooled       16         14      5             6     12.109074      9.674935     13.033431     63.241411
     exclude_8PBZ_8PC1 assembly_equal       16         14      5             6     11.297715      9.798544     13.894668     63.025623
     exclude_8PBZ_8PC1      PDB_equal       16         14      5             6     10.884731      9.909137     14.343378     62.863381
exclude_8PBZ_8PC0_8PC1         pooled       15         13      4             5     12.454575      9.696069     12.669310     63.263958
exclude_8PBZ_8PC0_8PC1 assembly_equal       15         13      4             5     12.092727      9.881822     13.058041     63.044935
exclude_8PBZ_8PC0_8PC1      PDB_equal       15         13      4             5     11.775250     10.040884     13.409771     62.846960
```

5. **低分辨率是否影响答案？** 排除既有低分辨率/几何疑点结构后，3–4 Å 比例升高，原始“弱峰”部分由这些结构贡献。因此不能将弱峰直接解释为物种特性。主结果没有删除这些结构。

6. **界面层面如何？** 下表用等权界面统计；零接触界面不虚构距离或条件比例。界面大小、总界面残基数和覆盖比例另见 interface_level_distributions.csv。

```
          dataset          metric  n_defined     mean   median      q25      q75
        Pedro AF3    median_d_min      54109 6.617567 6.559880 6.411094 6.771893
        Pedro AF3    fraction_3_4      54109 0.063497 0.068966 0.026316 0.091463
        Pedro AF3    fraction_4_5      54109 0.095409 0.101266 0.064516 0.125000
        Pedro AF3    fraction_5_6      54109 0.151538 0.152941 0.117647 0.181818
        Pedro AF3    fraction_6_8      54109 0.661960 0.629560 0.592593 0.690476
        Pedro AF3 fraction_lt_2_0      54109 0.003272 0.000000 0.000000 0.002899
Mgen Experimental    median_d_min         21 6.608968 6.585218 6.541344 6.616879
Mgen Experimental    fraction_3_4         21 0.106543 0.105882 0.087391 0.132895
Mgen Experimental    fraction_4_5         21 0.101845 0.099454 0.092105 0.116519
Mgen Experimental    fraction_5_6         21 0.131518 0.140118 0.112360 0.154639
Mgen Experimental    fraction_6_8         21 0.639014 0.632959 0.620765 0.642480
Mgen Experimental fraction_lt_2_0         21 0.001722 0.000000 0.000000 0.000000
```

7. **极短接触是否主要来自已有问题结构？** 同种实验 10,528 个接触中有 26 个 <2 Å（0.247%）；25 个在 8PBZ，另 1 个在 8PBY。排除 8PBZ/8PC1 后只剩 1/6,491（0.0154%），再排除 8PC0 后为 1/6,054（0.0165%）。Pedro 原始为 0.5034%，已有 QC 后仍为 0.4337%；FoldBench/SAFES AF3 为 0.0588%/0.0506%。不能把所有 <2 Å 一律称作碰撞；接触层面的比例不是独立样本显著性证据。

8. **是否支持物种解释？** 证据不足。参考主集只有 7 个 PDB、8 个装配，23 个链对中 2 个零接触，而且基本集中于同一黏附复合物及相关研究；分辨率、建模方式、蛋白组成和装配重复混杂。弱峰随敏感性处理而减弱，不能归因为 M. genitalium 的普遍几何特征，也不能定义该物种真实全局 PPI 分布。

## 化学与短距离补充

```
          dataset           contact_class  pooled_fraction  mean_interface_fraction
        Pedro AF3 hydrophobic_hydrophobic         0.182764                 0.131890
        Pedro AF3         opposite_charge         0.048591                 0.071408
        Pedro AF3             same_charge         0.038495                 0.065949
        Pedro AF3             polar_polar         0.067736                 0.067129
        Pedro AF3       hydrophobic_polar         0.184538                 0.152189
        Pedro AF3     hydrophobic_charged         0.183381                 0.189003
        Pedro AF3           polar_charged         0.144937                 0.180176
        Pedro AF3        special_or_other         0.149558                 0.142256
Mgen Experimental hydrophobic_hydrophobic         0.146562                 0.198842
Mgen Experimental         opposite_charge         0.027356                 0.018215
Mgen Experimental             same_charge         0.018332                 0.014045
Mgen Experimental             polar_polar         0.120916                 0.113560
Mgen Experimental       hydrophobic_polar         0.238032                 0.234742
Mgen Experimental     hydrophobic_charged         0.127185                 0.107897
Mgen Experimental           polar_charged         0.133359                 0.103424
Mgen Experimental        special_or_other         0.188260                 0.209274
```
化学分类只是组成描述；参考集以少数黏附复合物为主，与全蛋白组预测不匹配，不应把比例差异直接归因于物种或 AF3 的化学偏好。

### 等权结论及分位数约定
原始同种实验 (3,4] Å 比例为 10.75%；装配等权为 10.76%，PDB 等权为 10.38%，原始视图的结论对这两种权重较稳健。排除 8PBZ/8PC1 后，池化/装配等权/PDB 等权分别为 12.11%/11.30%/10.88%；再排除 8PC0 为 12.45%/12.09%/11.78%。联合敏感性说明“更接近一般实验”的程度仍受权重影响，不应声称分布完全一致。
本次池化/加权距离分位数统一使用经验累积分布的逆函数（首次达到目标累计权重的观测值），偶数样本中位数可能与旧报告的线性插值差一个相邻距离间隔；这不改变距离定义。界面内统计仍直接复用原 summarize 函数。
