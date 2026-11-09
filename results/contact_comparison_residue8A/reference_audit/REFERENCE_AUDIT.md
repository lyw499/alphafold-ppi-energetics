# 实验参考集分布审计

本审计保留全部29条参考链对（包含2条零接触），不改变现有Pedro比较、物种范围或QC选择。当前图的25条非零参考 = 至少一条Mgen链的27条候选减去2条零接触；另外2条双方均非Mgen的非零链对仍在全部29条审计中。

## 1. 样本范围与物种组成

| scope | population | n_zero_in_parent_scope | pct_zero_in_parent_scope | n_pdb | n_assemblies | n | median | q1 | q3 | mean | min | max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_reference_pairs | all_including_zeros | 2 | 6.897 | 9 | 10 | 29 | 480 | 89 | 758 | 437.7 | 0 | 915 |
| all_reference_pairs | nonzero_only | 2 | 6.897 | 9 | 10 | 27 | 494 | 93.5 | 759 | 470.1 | 50 | 915 |
| Mgen_Mgen_only | all_including_zeros | 2 | 8.696 | 7 | 8 | 23 | 678 | 86.5 | 765.5 | 457.7 | 0 | 915 |
| Mgen_Mgen_only | nonzero_only | 2 | 8.696 | 7 | 8 | 21 | 699 | 89 | 771 | 501.3 | 50 | 915 |
| mixed_species_only | all_including_zeros | 0 | 0 | 2 | 2 | 4 | 285.5 | 127 | 444 | 285.5 | 91 | 480 |
| mixed_species_only | nonzero_only | 0 | 0 | 2 | 2 | 4 | 285.5 | 127 | 444 | 285.5 | 91 | 480 |
| non_Mgen_non_Mgen | all_including_zeros | 0 | 0 | 2 | 2 | 2 | 511 | 502.5 | 519.5 | 511 | 494 | 528 |
| non_Mgen_non_Mgen | nonzero_only | 0 | 0 | 2 | 2 | 2 | 511 | 502.5 | 519.5 | 511 | 494 | 528 |
| current_Mgen_containing_scope | all_including_zeros | 2 | 7.407 | 9 | 10 | 27 | 437 | 88.5 | 759 | 432.2 | 0 | 915 |
| current_Mgen_containing_scope | nonzero_only | 2 | 7.407 | 9 | 10 | 25 | 480 | 91 | 760 | 466.8 | 50 | 915 |

物种按蛋白基因来源判断，不把表达宿主当作来源。23条Mgen–Mgen、4条Mgen–human混合物种、2条human–human；后者不是混合物种，也不是Mgen–Mgen。最大的5条和10条全部是Mgen–Mgen，因此混合物种不是最大接口的主要来源。

## 2. 大接口集中程度与25条观测的实际来源

| scope | rank_count | n_unique_PDB | n_assemblies | n_mixed | sum_contacts | fraction_of_scope_contact_sum | PDBs |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all29 | 5 | 4 | 4 | 0 | 4240 | 0.3341 | 8PBX;8PBY;8PBZ;8PC1 |
| all29 | 10 | 6 | 7 | 0 | 8044 | 0.6338 | 6RUT;6YRK;8PBX;8PBY;8PBZ;8PC1 |
| current25 | 5 | 4 | 4 | 0 | 4240 | 0.3633 | 8PBX;8PBY;8PBZ;8PC1 |
| current25 | 10 | 6 | 7 | 0 | 8044 | 0.6893 | 6RUT;6YRK;8PBX;8PBY;8PBZ;8PC1 |

当前25条非零观测来自 **9 个PDB、10 个biological assembly**。贡献如下：

| pdb_id | structure_id | n_observations | sum_contacts | median_contacts | max_contacts | percent_of_25 |
| --- | --- | --- | --- | --- | --- | --- |
| 4NZR | 4nzr-assembly1 | 2 | 619 | 309.5 | 480 | 8 |
| 4NZT | 4nzt-assembly1 | 2 | 523 | 261.5 | 432 | 8 |
| 6RUT | 6rut-assembly1 | 5 | 1793 | 96 | 771 | 20 |
| 6RUT | 6rut-assembly2 | 5 | 1777 | 89 | 760 | 20 |
| 6YRK | 6yrk-assembly1 | 1 | 757 | 757 | 757 | 4 |
| 8PBX | 8pbx-assembly1 | 1 | 865 | 865 | 865 | 4 |
| 8PBY | 8pby-assembly1 | 1 | 862 | 862 | 862 | 4 |
| 8PBZ | 8pbz-assembly1 | 6 | 3122 | 688.5 | 801 | 24 |
| 8PC0 | 8pc0-assembly1 | 1 | 437 | 437 | 437 | 4 |
| 8PC1 | 8pc1-assembly1 | 1 | 915 | 915 | 915 | 4 |

6RUT的两个assembly分别贡献5条非零观测，8PBZ贡献6条，合计16/25（64%）。同一复合物的多个链对不能视作独立实验。前5大接口来自4个PDB/assembly，其中8PBZ贡献2条。23条Mgen–Mgen候选集中于同一复合物体系的多个结构/构象及链拷贝，而不是从Mgen蛋白组均匀抽样；4NZR/4NZT则是另一类混合物种复合物。蛋白名称和物种逐条列于主CSV。

## 3. 链长能解释多少

| group | metric | n | median | q1 | q3 | mean | min | max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Pedro_raw_nonzero | contact_count | 54109 | 69 | 16 | 165 | 123.3 | 1 | 1.802e+04 |
| Pedro_raw_nonzero | total_chain_length | 54109 | 745 | 535 | 1040 | 839.1 | 85 | 3588 |
| Pedro_raw_nonzero | total_interface_residues | 54109 | 41 | 15 | 81 | 59.21 | 2 | 1645 |
| Pedro_raw_nonzero | interface_fraction | 54109 | 0.05556 | 0.02012 | 0.1102 | 0.07981 | 0.0007249 | 0.9495 |
| reference_current25 | contact_count | 25 | 480 | 91 | 760 | 466.8 | 50 | 915 |
| reference_current25 | total_chain_length | 25 | 2049 | 2028 | 2383 | 1949 | 587 | 2760 |
| reference_current25 | total_interface_residues | 25 | 249 | 57 | 296 | 196.6 | 41 | 395 |
| reference_current25 | interface_fraction | 25 | 0.1251 | 0.02318 | 0.1491 | 0.1109 | 0.01998 | 0.2572 |
| reference_all27_nonzero | contact_count | 27 | 494 | 93.5 | 759 | 470.1 | 50 | 915 |
| reference_all27_nonzero | total_chain_length | 27 | 2048 | 1972 | 2285 | 1837 | 434 | 2760 |
| reference_all27_nonzero | total_interface_residues | 27 | 193 | 59.5 | 292 | 196.1 | 41 | 395 |
| reference_all27_nonzero | interface_fraction | 27 | 0.1363 | 0.02536 | 0.1583 | 0.1348 | 0.01998 | 0.4367 |

| scope | population | predictor | n | pearson | spearman | n_unique_PDB | note |
| --- | --- | --- | --- | --- | --- | --- | --- |
| all_reference_pairs | all_including_zeros | total_chain_length | 29 | 0.1619 | 0.2389 | 9 | Descriptive only; clustered observations, no independence-based significance test. |
| all_reference_pairs | all_including_zeros | total_interface_residues | 29 | 0.9741 | 0.9484 | 9 | Descriptive only; clustered observations, no independence-based significance test. |
| all_reference_pairs | nonzero_only | total_chain_length | 27 | 0.1355 | 0.146 | 9 | Descriptive only; clustered observations, no independence-based significance test. |
| all_reference_pairs | nonzero_only | total_interface_residues | 27 | 0.9702 | 0.936 | 9 | Descriptive only; clustered observations, no independence-based significance test. |
| Mgen_Mgen_only | all_including_zeros | total_chain_length | 23 | 0.1708 | 0.2957 | 7 | Descriptive only; clustered observations, no independence-based significance test. |
| Mgen_Mgen_only | all_including_zeros | total_interface_residues | 23 | 0.9772 | 0.9063 | 7 | Descriptive only; clustered observations, no independence-based significance test. |
| Mgen_Mgen_only | nonzero_only | total_chain_length | 21 | -0.1367 | 0.07376 | 7 | Descriptive only; clustered observations, no independence-based significance test. |
| Mgen_Mgen_only | nonzero_only | total_interface_residues | 21 | 0.9731 | 0.8769 | 7 | Descriptive only; clustered observations, no independence-based significance test. |
| mixed_species_only | all_including_zeros | total_chain_length | 4 | -0.3885 | 0 | 2 | Descriptive only; clustered observations, no independence-based significance test. |
| mixed_species_only | all_including_zeros | total_interface_residues | 4 | 0.9967 | 1 | 2 | Descriptive only; clustered observations, no independence-based significance test. |
| mixed_species_only | nonzero_only | total_chain_length | 4 | -0.3885 | 0 | 2 | Descriptive only; clustered observations, no independence-based significance test. |
| mixed_species_only | nonzero_only | total_interface_residues | 4 | 0.9967 | 1 | 2 | Descriptive only; clustered observations, no independence-based significance test. |
| current_Mgen_containing_scope | all_including_zeros | total_chain_length | 27 | 0.2322 | 0.2773 | 9 | Descriptive only; clustered observations, no independence-based significance test. |
| current_Mgen_containing_scope | all_including_zeros | total_interface_residues | 27 | 0.9751 | 0.9391 | 9 | Descriptive only; clustered observations, no independence-based significance test. |
| current_Mgen_containing_scope | nonzero_only | total_chain_length | 25 | 0.1871 | 0.1555 | 9 | Descriptive only; clustered observations, no independence-based significance test. |
| current_Mgen_containing_scope | nonzero_only | total_interface_residues | 25 | 0.9715 | 0.9232 | 9 | Descriptive only; clustered observations, no independence-based significance test. |

当前25条参考非零观测总链长中位数2049，Pedro非零为745（约2.75倍）；接触数中位数480 vs69（约6.96倍），接口比例中位数0.1251 vs0.0556。参考链更长有助于解释跨集合差异，但在参考内部，总链长与接触数相关性很弱：全部29条Pearson=0.162、Spearman=0.239；当前25条Pearson=0.187、Spearman=0.155。因此不能把大接触数主要归因于链长。总界面残基与接触数相关性强（29条Pearson=0.974、Spearman=0.948）。长度与接触数的相关性是描述性的，不是因果解释。总链长相近的同一四链assembly仍存在0接触、小接口和大接口，所以仅靠链长无法解释所有差异。总界面残基数与接触数来自同一个接触图，具有数学耦合；强相关不能作为独立生物学证据。样本小且重复/相关，未给出假设29条相互独立的显著性检验。

## 4. 重复链对、拷贝和对称性

所有原始链对均唯一，没有同一无序链对的重复行；但相同蛋白实体组成的不同链拷贝被当作不同观测。`chain_source_and_copy_audit.csv`保存实体/完整序列hash/原始链/应用operator；`within_assembly_pair_similarity.csv`对相同实体组合的链对进行共同CA的刚体拟合，记录RMSD和接触数差。此拟合用于识别相似拷贝，不等同于原位原子重叠。跨PDB的共同citation及序列pair key也保留在审计表，序列相同不等于结构或界面必然相同。

| structure_id | interface_1 | interface_2 | same_entity_pair | same_full_entity_sequence_pair | n_common_CA | rigid_fit_CA_RMSD_A | aligned_chain_order_2 | contact_count_1 | contact_count_2 | contact_difference | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6rut-assembly2 | 6rut-assembly2:E-2:F-2 | 6rut-assembly2:C:D | True | True | 2028 | 0.3412 | C:D | 760 | 758 | 2 | CA superposition similarity is not direct spatial overlap; neither record removed. |
| 6rut-assembly1 | 6rut-assembly1:G-2:H-2 | 6rut-assembly1:A:B | True | True | 2025 | 0.4265 | A:B | 771 | 758 | 13 | CA superposition similarity is not direct spatial overlap; neither record removed. |
| 6rut-assembly1 | 6rut-assembly1:A:H-2 | 6rut-assembly1:B:G-2 | True | True | 2025 | 0.6459 | G-2:B | 96 | 85 | 11 | CA superposition similarity is not direct spatial overlap; neither record removed. |
| 8pbz-assembly1 | 8pbz-assembly1:C:D | 8pbz-assembly1:A:B | True | True | 2383 | 0.6997 | A:B | 801 | 797 | 4 | CA superposition similarity is not direct spatial overlap; neither record removed. |
| 6rut-assembly2 | 6rut-assembly2:D:E-2 | 6rut-assembly2:C:F-2 | True | True | 2028 | 0.7143 | F-2:C | 89 | 88 | 1 | CA superposition similarity is not direct spatial overlap; neither record removed. |
| 8pbz-assembly1 | 8pbz-assembly1:B:C | 8pbz-assembly1:A:D | True | True | 2383 | 1.535 | D:A | 699 | 678 | 21 | CA superposition similarity is not direct spatial overlap; neither record removed. |
| 8pbz-assembly1 | 8pbz-assembly1:A:B | 8pbz-assembly1:B:C | True | True | 2383 | 27.85 | C:B | 797 | 699 | 98 | CA superposition similarity is not direct spatial overlap; neither record removed. |
| 8pbz-assembly1 | 8pbz-assembly1:C:D | 8pbz-assembly1:B:C | True | True | 2383 | 27.89 | C:B | 801 | 699 | 102 | CA superposition similarity is not direct spatial overlap; neither record removed. |
| 8pbz-assembly1 | 8pbz-assembly1:A:B | 8pbz-assembly1:A:D | True | True | 2383 | 27.91 | A:D | 797 | 678 | 119 | CA superposition similarity is not direct spatial overlap; neither record removed. |
| 8pbz-assembly1 | 8pbz-assembly1:C:D | 8pbz-assembly1:A:D | True | True | 2383 | 27.94 | A:D | 801 | 678 | 123 | CA superposition similarity is not direct spatial overlap; neither record removed. |

8PBZ的A:B与C:D都是同一两蛋白实体组合的拷贝；其两个大接口为797和801。6RUT的两套assembly也包含重复的两实体组合，既有大界面也有较小的拷贝间邻接界面；无证据支持把所有链对都称作独立、等价的预期生物学界面。

## 5. biological assembly 文件验证

复用下载manifest中明确的 RCSB `{pdb}-assembly{id}.cif.gz`，并将每个缓存文件的 `_pdbx_chain_remapping` 原始label链和applied_operations与entry header的 `_pdbx_struct_assembly_gen` 逐项核对。10/10文件的全部坐标链（包括非蛋白链）及operator集合精确吻合，蛋白链数也与声明的oligomeric_count一致。6RUT正确使用两个不同assembly，包含operator 2生成的后缀链，而非把整个非对称单元的所有蛋白链混入。

8PBZ的四蛋白链属于作者指定tetrameric assembly；8PC0/8PC1的两蛋白链分别属于声明的dimeric assembly。当biological assembly与非对称单元在链组成上相同时，不应仅因两者一致就认定为晶体接触。当前未发现误用非预期ASU文件的证据。

核验范围：下载来源、链/实体映射、声明的assembly生成操作与组成；没有重新下载原始ASU坐标并逐原子重建矩阵变换。正确使用数据库assembly不等于每个枚举链对的生物学意义已独立证明，尤其6RUT的assembly标为software_defined。

| structure_id | assembly_details | oligomeric_details | method | resolution_A | assembly_remapping_verified |
| --- | --- | --- | --- | --- | --- |
| 4nzr-assembly1 | author_and_software_defined_assembly | trimeric | X-RAY DIFFRACTION | 1.650 | True |
| 4nzt-assembly1 | author_and_software_defined_assembly | trimeric | X-RAY DIFFRACTION | 2.497 | True |
| 6rut-assembly1 | software_defined_assembly | tetrameric | X-RAY DIFFRACTION | 2.65 | True |
| 6rut-assembly2 | software_defined_assembly | tetrameric | X-RAY DIFFRACTION | 2.65 | True |
| 6yrk-assembly1 | author_and_software_defined_assembly | dimeric | ELECTRON MICROSCOPY | 4.1 | True |
| 8pbx-assembly1 | author_defined_assembly | dimeric | ELECTRON MICROSCOPY | 3.3 | True |
| 8pby-assembly1 | author_defined_assembly | dimeric | ELECTRON MICROSCOPY | 3.7 | True |
| 8pbz-assembly1 | author_defined_assembly | tetrameric | ELECTRON MICROSCOPY | 11 | True |
| 8pc0-assembly1 | author_and_software_defined_assembly | dimeric | ELECTRON MICROSCOPY | 17 | True |
| 8pc1-assembly1 | author_defined_assembly | dimeric | ELECTRON MICROSCOPY | 18 | True |

## 6. 前5大接口逐项检查

| interface_id | contact_rank | contact_count | min_distance_rechecked_A | n_crosschain_atom_pairs_distance_lt_0_5A | closest_atom_elements | n_crosschain_atom_pairs_lt_2A_diagnostic_only | n_declared_covalent_links_between_chains | local_geometry_note | malformed_duplicate_identity | parser_warnings | n_models | resolution_A | total_chain_length | total_interface_residues | interface_fraction | n_other_same_entity_pairs_in_assembly | closest_other_pair_CA_RMSD_A | assessment |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 8pc1-assembly1:A:B | 1 | 915 | 2.461 | 0 | O-O | 0 | 0 | No cross-chain heavy-atom pair below 2 A in diagnostic check | False |  | 1 | 18 | 2187 | 365 | 0.1669 | 0 | nan | No severe overlap detected by existing <0.5 A criterion; assembly membership verified. Large coordinate interface does not by itself establish atomic accuracy or independent biological evidence. |
| 8pbx-assembly1:A:B | 2 | 865 | 2.3 | 0 | O-O | 0 | 0 | No cross-chain heavy-atom pair below 2 A in diagnostic check | False |  | 1 | 3.3 | 2187 | 326 | 0.1491 | 0 | nan | No severe overlap detected by existing <0.5 A criterion; assembly membership verified. Large coordinate interface does not by itself establish atomic accuracy or independent biological evidence. |
| 8pby-assembly1:A:B | 3 | 862 | 1.973 | 0 | O-O | 1 | 0 | Short non-hydrogen atom separation; potential local clash, not an automatic exclusion | False |  | 1 | 3.7 | 2187 | 331 | 0.1513 | 0 | nan | No severe overlap detected by existing <0.5 A criterion; assembly membership verified. Large coordinate interface does not by itself establish atomic accuracy or independent biological evidence. |
| 8pbz-assembly1:C:D | 4 | 801 | 1.412 | 0 | C-C | 4 | 0 | Short non-hydrogen atom separation; potential local clash, not an automatic exclusion | False |  | 1 | 11 | 2383 | 395 | 0.1658 | 3 | 0.6997 | No severe overlap detected by existing <0.5 A criterion; assembly membership verified. Large coordinate interface does not by itself establish atomic accuracy or independent biological evidence. |
| 8pbz-assembly1:A:B | 5 | 797 | 1.122 | 0 | C-C | 4 | 0 | Short non-hydrogen atom separation; potential local clash, not an automatic exclusion | False |  | 1 | 11 | 2383 | 394 | 0.1653 | 3 | 0.6997 | No severe overlap detected by existing <0.5 A criterion; assembly membership verified. Large coordinate interface does not by itself establish atomic accuracy or independent biological evidence. |

前5条最小链间距离从缓存坐标独立复核，未发现已有 <0.5 Å 严重重叠或重复原子身份标记。前3条是声明二聚体中的唯一蛋白链对；第4/5条是同一8PBZ四链assembly的两组相似拷贝。它们在所用坐标模型中确有较大界面，不能简单归因于原子穿插；但不能因此证明近原子级接触数量准确。

局部几何诊断补充：8PBZ A:B和C:D的最短C–C距离分别为1.122和1.412 Å，每对均有4个重原子对小于2 Å；header没有这两条蛋白链之间的共价连接注释。这支持可能存在局部碰撞/建模几何问题，不能因未触发<0.5 Å就称其无碰撞。8PBY还有1个O–O对约1.973 Å，值得保留为局部短距离提示。2 Å仅用于透明展示短距离诊断数量，不是新QC排除规则；未删除任何案例。

尤其最高的8PC1分辨率18 Å，8PBZ为11 Å；另有8PC0为17 Å。这些低分辨率EM模型的8 Å逐原子距离判定高度依赖所建模型，不能等同于直接解析出的原子级实验接触。8PBX/8PBY为3.3/3.7 Å，也应按结构分辨率和建模背景解释。未把分辨率变为新的自动过滤阈值。

## 7. 结论

- **较大生物学复合物的界面：** 大接口存在于声明的biological assemblies中，坐标层面无严重重叠证据；但低分辨率模型限制了对精确接触数和每个链对生物学作用的确信。
- **数据集/选择偏差：** 参考集是少数已解析复合物的富集样本，Pedro是全蛋白组池化预测；即使都只看非零接触，二者也不是同一候选选择机制。参考链通常更长，界面残基也更多。
- **同一复合物的重复观测：** 当前25条中16条来自6RUT两套assembly及8PBZ；多个结构/构象还代表相同蛋白体系，不能把25当作25个独立复合物。
- **混合物种效应：** 混合物种仅4条，前10大接口中0条；不是参考大接口尾部的主要原因。另有2条human–human只属于29条全量审计，不在当前25条图中。
- **可能的结构/assembly问题：** 未发现文件选错ASU或已有<0.5 Å严重重叠。8PBZ的1.12/1.41 Å C–C跨链距离提示局部碰撞，8PBY有1.97 Å O–O短距离；加上低分辨率模型依赖，需要保留几何疑点。它们并不证明整个大界面都是伪接口；assembly生物学解释亦有边界。

本轮不删除、不重加权、不替换任何原结果；物种子集仅作为敏感性审计。

## 输出

- `reference_interface_audit.csv`：全部29条，按残基对接触数降序；链ID均为计算使用的label_asym_id，另列author/original链映射。
- `species_sensitivity_summary.csv` / `species_contact_distributions.png|pdf`：所有/非零两套物种范围分布。
- `contact_correlations.csv` / `contact_size_diagnostics.png|pdf`：长度及总界面残基的相关性和散点；图中1–5标记主表排名。
- `assembly_audit.csv`、`chain_source_and_copy_audit.csv`、`within_assembly_pair_similarity.csv`、`top5_structural_audit.csv`：assembly、拷贝和top5证据。
- `current25_assembly_contributions.csv`、`current25_pdb_contributions.csv`：当前图25条观测的来源计数。
- `top_rank_concentration.csv`、`pedro_reference_length_context.csv`：集中度及Pedro长度背景。
- `validation.json`：旧输入与图的SHA256未变，所有29条原指标保留。
- `run_audit.py`：复现脚本，输出仅限本审计目录。
