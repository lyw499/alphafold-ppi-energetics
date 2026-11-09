"""Report the completed table-only size-matching analysis without modifying inputs."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'results/local_interface_geometry/size_matched_packing'
REFS=['FoldBench Experimental','SAFES Experimental','Mgen Experimental','Mgen excl 8PBZ 8PC1','Mgen excl 8PBZ 8PC0 8PC1']
TYPES=['BB-BB','BB-SC','SC-SC']
def table(d):
    return '\n```\n'+d.to_string(index=False,float_format=lambda x:f'{x:.4f}')+'\n```\n'
def primary(d):
    return d[(d.size_variable=='n_contacts')&(d.pedro_version=='raw')&(d.nearest_k==1)&d.reference.isin(REFS)].copy()
def main():
    assert json.loads((P/'validation.json').read_text())['status']=='passed'
    audit=pd.read_csv(P/'data_audit.csv');balance=pd.read_csv(P/'matching_balance_summary.csv')
    diff=pd.read_csv(P/'size_matched_effect_differences.csv');out=pd.read_csv(P/'size_matched_close_contact_summary.csv')
    bins=pd.read_csv(P/'size_matched_distance_bin_summary.csv');strata=pd.read_csv(P/'size_stratified_contact_type_summary.csv')
    corr=pd.read_csv(P/'pedro_size_contact_type_correlations.csv');paired=pd.read_csv(P/'benchmark_paired_size_context.csv')
    manifest=pd.read_csv(P/'input_manifest.csv');outside=pd.read_csv(P/'outside_size_support.csv')
    main=primary(diff);total=main[main.metric=='total_close'].copy()
    for c in ['experimental_mean','pedro_before_mean','pedro_matched_mean','before_difference','matched_difference','change_after_matching']:total[c]*=100
    total['gap_attenuation_percent']=100*(1-total.matched_difference/total.before_difference)
    main_table=total[['reference','experimental_mean','pedro_before_mean','pedro_matched_mean','before_difference','matched_difference','gap_attenuation_percent']]
    contributions=main[main.metric.isin([t+'_close' for t in TYPES])].pivot(index='reference',columns='metric',values='matched_difference')*100
    drivers=contributions.idxmin(axis=1).str.replace('_close','').to_dict()
    conditional=main[main.metric.isin(['BB-SC_share','SC-SC_share']+[t+'_within_'+b for t in ['BB-SC','SC-SC'] for b in ['3to4','4to5','5to6','6to8']])].copy()
    conditional[['experimental_mean','pedro_matched_mean','matched_difference']]*=100
    lines=['# 控制界面大小后，Pedro 是否仍有局部接触几何差异？\n',
    '## 本轮回答什么\n\n此前池化残基接触中，Pedro的(3,4] Å比例较低。本轮把界面作为主要观察单位：每个实验界面先得到自身比例，再与接触数可比的Pedro界面比较。目标是判断大小组成能否充分解释差异，而不是建立因果匹配研究或评估能量函数。\n',
    '## 输入审计与零界面\n\n主输入是上一轮interface_contact_type_summary.parquet，包含每界面的三类计数、类型总占比和全部六个距离箱。Pedro/Mgen的两项次要大小指标直接读既有界面表；FoldBench/SAFES按既有接触CSV中的残基索引去重得到界面残基数，再除以既有两链长度之和。没有重新计算距离或读取结构。\n',
    table(audit[['dataset','n_interfaces','n_zero','n_nonzero','zero_percent']]),
    '所有必要计数、总体近接触比例及三项大小指标齐全。某非零界面没有某接触类型时，该类型的全接触贡献为0，类型内部比例为NA；这是定义域缺失而非坏数据。零界面保留在analysis_interface_table.parquet和审计表，比例分析只使用非零界面。未把FoldBench不可用AF3案例当作零。\n',
    '输入文件及哈希清单：\n',table(manifest[['path','bytes']]),
    '## 分母、权重和精确方法\n\n'
    '- 总近接触比例：界面内(3,4] Å残基对数/全部≤8 Å残基对数。\n'
    '- 类型贡献：界面内该类型(3,4] Å残基对数/全部≤8 Å残基对数；三类贡献逐界面相加为总比例。\n'
    '- 类型内部比例：界面内该类型某距离箱的计数/该类型所有接触数；只在类型存在时定义。\n'
    '- 主结果对界面等权；匹配Pedro使用下述权重。堆叠图采用均值，因此贡献可相加。中位数/IQR另存表，不能把三类中位数相加。\n'
    '- 全距离范围沿用<2、[2,3]、(3,4]、(4,5]、(5,6]、(6,8] Å六箱。图07/08显示全范围分箱概率，而非假定箱内连续密度或精细峰形。\n'
    '- BB名称仍为N、CA、C、O；OXT沿上一轮主定义作为明确的末端例外。没有改变类型或距离定义。\n\n'
    '每套实验参考分别匹配：对每个实验界面，在Pedro非零界面中寻找log1p(contact_count)绝对距离最小的候选；所有并列最近候选平均分配该目标的1份权重，可被不同实验界面重复使用。这样不根据几何结果挑选候选，也避免在相同接触数界面中随意挑一个。总体权重和等于支持范围内实验界面数。k=20敏感性扩大到第20近距离以内的所有并列候选。\n\n'
    '超出Pedro大小范围的实验界面单独记录，不作边界强行外推。每套参考重新匹配。权重有效样本数ESS=(Σw)²/Σw²仅描述权重集中度，不是独立生物学样本量；界面仍嵌套在结构/复合物内。分位数统一采用加权经验CDF反函数，故偶数样本的中位数可能不同于相邻两点平均的惯例。\n',
    '## 一、共同绝对大小分层及支持度\n\nPedro非零接触数四分位边界为16、69、165，统一应用为1–16、17–69、70–165、≥166。最大层较宽，所以分层仅用于直观检查，主要大小控制依据更精确的最近邻匹配。n<5只标为低支持，不删除观察，不作为统计显著性门槛。\n',
    table(strata[strata.metric.eq('total_close')][['dataset','stratum','n_interfaces','median_contact_count','median','q25','q75']]),
    '上表比例为0–1；图02–04显示百分比。FoldBench最小层仅3个实验界面，SAFES最小层没有实验界面，Mgen的17–69层仅1个。它们不能提供强分层证据。完整三类及类型内部指标均保存在size_stratified_contact_type_summary.csv。\n',
    '## 二、匹配是否真的改善了大小平衡\n\n先看大小，后看几何。下表n是有正权重的独特Pedro界面数或实验界面数；effective_n是上述权重ESS。标准化差用log1p接触数的均值差除以匹配前实验/Pedro合并标准差，前后使用同一分母。weighted_ecdf_max_gap为两条加权经验CDF的最大差，不是检验p值。\n',
    table(primary(balance)[['reference','phase','n','effective_n','median','q25','q75','standardized_mean_difference','weighted_ecdf_max_gap','n_experimental_outside_support']]),
    '## 三、匹配前后的总近接触差异\n\n下表均为界面均值的百分比/百分点；负差表示Pedro较低。gap_attenuation_percent=(1−匹配后差/匹配前差)×100，只是描述性差值衰减，不能解释为大小的因果解释率。\n',table(main_table),
    '这些数值不应直接与上一轮“全部残基接触池化比例”相减：池化隐含按接触数给大界面更大权重，本轮首先改成界面等权，再在同一估计口径内比较匹配前后。\n',
    '## 四、匹配后哪一种类型贡献最多\n\nPedro减实验的三类界面均值差，单位百分点：\n',table(contributions.reset_index()),
    '逐参考最大的单项净缺口为：'+ '；'.join(f'{r}：{drivers[r]}' for r in REFS)+'。\n\n'
    '必须允许主导类别随参考及统计权重改变，不能预设BB-SC一定最大。此处是匹配后的界面均值贡献，上一轮主结论来自池化接触；两者区别应通过同口径的匹配前后表判断，不能把所有排序变化归因于大小匹配。\n',
    '## 五、类型数量还是类型内部距离组成\n\nshare是类型占界面全部接触的比例；within是该类型内部距离比例。以下列出实验、匹配Pedro及两者差，全部换成百分比/百分点。\n',
    table(conditional[['reference','metric','experimental_mean','pedro_matched_mean','matched_difference']]),
    '图07/08及size_matched_distance_bin_summary.csv保留0–8 Å六箱，含短于3 Å部分。主列mean是等界面的类型内部比例均值；pooled_fraction另外给出加权计数池化比例，分母不同，不混用。类型不存在的界面不进入within均值，相关n及ESS逐指标保存。任何“较长区间富集”均是跨集合组成描述，不表示对应残基对真实移动。\n',
    '## 六、连续大小关联和既有QC\n\nSpearman在全部界面上计算，log1p为单调变换，故排序相关与原接触数相同；图09的分箱仅用于显示。\n',
    table(corr[corr.size_variable.eq('n_contacts')]),
    '相关系数不可相加，也不是解释方差。BB-SC/SC-SC各自的全分母贡献与条件比例都报告；不将一个较大的相关系数解释成唯一机制。\n',
    '## 七、敏感性：不同大小指标、最近邻数及原有QC\n\n各次敏感性重新匹配，原始结果始终保留。总界面残基数使用log1p，界面比例直接使用原比例距离。单独匹配一个次要指标不等于同时平衡三个变量；不能替代接触数主分析。QC仅沿用existing_severe_overlap，无新阈值。\n',
    table(diff[(diff.metric=='total_close')&diff.reference.isin(REFS)][['reference','size_variable','pedro_version','nearest_k','before_difference','matched_difference']]),
    '上表差值为0–1比例；乘100得百分点。次要指标的平衡、域外名单、样本数同样见matching_balance_summary.csv和outside_size_support.csv。界面比例是近连续变量，最近邻可能只有一个，需结合ESS看待，不能与接触数大量精确并列匹配具有同等稳定性。\n',
    '## 八、FoldBench和SAFES保留配对的控制\n\n使用相同interface_id配对，双方非零且该指标可定义才计算差值；零接触仍在配对ID审计。保留全映射和明确链映射子集。为了描述配对大小变化的影响，拟合Δfraction=截距+斜率×Δlog1p(contact_count)，截距代表该线性模型在大小不变时的预测差，不能称为实验控制或因果效应；不报告p值。\n',
    table(paired[(paired.metric.isin(['total_close','BB-SC_close','SC-SC_close']))][['dataset','scope','metric','n','mean','median','q25','q75','size_adjusted_intercept','spearman']]),
    '另外对FoldBench AF3/SAFES AF3各自的接触数分布匹配Pedro，结果保存在主结果表中，作为补充而非替代原始配对。所有控制差值均按界面定义，不把残基接触当独立样本。\n',
    '## 九、Mgen敏感性与限制\n\n沿用全23界面（21非零）、排除8PBZ/8PC1的16界面（14非零）、再排除8PC0的15界面（13非零）。各自独立匹配，没有新增过滤门槛。Mgen包含同一装配多个链对及相关黏附复合物，因此只作描述性支持，不能把13–21个界面理解成同样数量的独立实验复合物。\n\n'
    '只匹配接触数不平衡蛋白种类、构象状态、长度、实验方法、序列、复合物选择或可信度。接触数本身也是几何派生量；观察到的“匹配后差异”不是证明某种物理机制。最近原子BB-SC几何更不等同于氢键几何。\n',
    '## 文件与验证\n\n'
    '输入计数与旧比例逐界面对账、每目标匹配权重和、零界面NA、输入文件前后SHA-256均通过。没有读取CIF/ZIP、重新计算原子距离、下载数据或进行新的化学相互作用分析。结构文件读取数为0。\n\n'
    '主表：size_stratified_contact_type_summary.csv、matching_balance_summary.csv、size_matched_close_contact_summary.csv、size_matched_distance_bin_summary.csv、pedro_size_contact_type_correlations.csv。补充：size_matched_effect_differences.csv、matched_interface_ids.csv、outside_size_support.csv、benchmark_paired_size_context.csv、benchmark_pair_ids.csv、data_audit.csv、analysis_interface_table.parquet、input_manifest.csv、methods.json、validation.json。图01–09位于figures/。\n']
    target=P/'SIZE_MATCHED_PACKING_CN.md'
    if target.exists():raise RuntimeError('Refuse to overwrite report')
    target.write_text('\n'.join(lines))
    main_table.to_csv(P/'key_findings.csv',index=False)
    print(main_table.to_string(index=False));print(contributions.to_string())
if __name__=='__main__':main()
