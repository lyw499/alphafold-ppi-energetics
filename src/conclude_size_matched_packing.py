"""Add checked interpretation to the current new analysis only."""
from pathlib import Path
import pandas as pd
import json,hashlib
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'results/local_interface_geometry/size_matched_packing'
def main():
    f=P/'SIZE_MATCHED_PACKING_CN.md'
    marker='## 十、九个问题的明确回答'
    assert marker not in f.read_text()
    answers='''
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
'''
    f.write_text(f.read_text()+answers)
    provenance=dict(scripts={n:hashlib.sha256((ROOT/'src'/n).read_bytes()).hexdigest() for n in ['size_matched_packing.py','size_matched_packing.sbatch','report_size_matched_packing.py','finalize_size_matched_packing.py','conclude_size_matched_packing.py']},input_manifest='input_manifest.csv',batch_job=65825882,output_directory=str(P),structural_parsing=False)
    (P/'provenance.json').write_text(json.dumps(provenance,indent=2))
    files=[dict(path=str(q),bytes=q.stat().st_size) for q in sorted(P.rglob('*')) if q.is_file()]
    pd.DataFrame(files).to_csv(P/'OUTPUT_FILES.csv',index=False)
    figures=list((P/'figures').glob('*.png'));assert len(figures)==9
    (P/'COMPLETED.json').write_text(json.dumps(dict(status='completed',figures=len(figures),report=str(f),prior_results_preserved=True),indent=2))
    print('Output directory:',P);print(pd.read_csv(P/'key_findings.csv').to_string(index=False));print('Report and 9 figures completed.')
if __name__=='__main__':main()
