"""Three-panel percent histograms from existing paired residue metrics only."""
import os,json,hashlib
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/paired-residue-triptychs')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/contact_comparison_residue8A'
def main():
 for name,title in [('foldbench','FoldBench'),('safes_provisional_two_chain','SAFES — provisional two-chain subset')]:
  source=BASE/f'{name}_paired_metrics.csv';digest=hashlib.sha256(source.read_bytes()).hexdigest();d=pd.read_csv(source)
  assert d.interface_id.is_unique
  out=BASE/'figures'/name;out.mkdir(exist_ok=True,parents=True);bins_rows=[];meta=[]
  for include_zero in [False,True]:
   version='all_including_zeros' if include_zero else 'nonzero'
   fig,axes=plt.subplots(1,3,figsize=(18,5.4))
   groups=[];zero=[]
   for prefix,label,color in [('experimental','Experimental','#1f77b4'),('af3','AF3','#ff7f0e')]:
    count=d[prefix+'_contact_count'];assert count.notna().all()
    zero.append(f'{label}: {int(count.eq(0).sum())}/{len(d)} zero')
    groups.append((d if include_zero else d[count>0],prefix,label,color))
   for ax,metric,label,width in zip(axes,['contact_count','total_interface_residues','interface_fraction'],['Residue-pair contacts (8 Å)','Total interface residues','Interface fraction: (IA + IB) / (LA + LB)'],[50,20,.05]):
    maximum=max(g[p+'_'+metric].max() for g,p,_,_ in groups)
    upper=1 if metric=='interface_fraction' else max(width,np.ceil(maximum/width)*width)
    edges=np.linspace(0,upper,int(round(upper/width))+1)
    for g,p,label_group,color in groups:
     vals=g[p+'_'+metric];hist,_=np.histogram(vals,bins=edges);assert hist.sum()==len(vals)
     pct=hist/len(vals)*100;assert np.isclose(pct.sum(),100)
     ax.stairs(pct,edges,label=f'{label_group} (n={len(vals)})',color=color,linewidth=2)
     for left,right,n,percentage in zip(edges[:-1],edges[1:],hist,pct):bins_rows.append(dict(version=version,group=label_group,metric=metric,bin_left=left,bin_right=right,n=int(n),percent=percentage,denominator=len(vals)))
    ax.set(xlabel=label,ylabel='Observations (%)',xlim=(0,upper));ax.legend(frameon=False,fontsize=10);ax.spines[['right','top']].set_visible(False)
   heading='all valid paired cases, zeros retained' if include_zero else 'non-zero interface distributions'
   fig.suptitle(f'{title}: Experimental vs AF3; {len(d)} valid matched interfaces\n{heading}; raw results, QC flags retained',fontsize=13)
   fig.text(.5,.02,'; '.join(zero)+'. Full value ranges shown; each group normalized separately.\nFraction is one length-weighted value per interface, not two pooled chain-level observations.',ha='center',fontsize=9)
   fig.tight_layout(rect=(0,.085,1,.9));stem=f'{name}_exp_vs_af3_residue_triptych_{version}'
   for ext in ['png','pdf']:fig.savefig(out/f'{stem}.{ext}',dpi=180)
   plt.close(fig)
   meta.append(dict(version=version,matched_interfaces=len(d),experimental_n=len(groups[0][0]),af3_n=len(groups[1][0]),experimental_zeros=int(d.experimental_contact_count.eq(0).sum()),af3_zeros=int(d.af3_contact_count.eq(0).sum()),source_sha256=digest,source=str(source.relative_to(ROOT)),fraction_definition='(IA+IB)/(LA+LB)',full_range=True,provisional='safes' in name))
  pd.DataFrame(bins_rows).to_csv(out/f'{name}_residue_triptych_bins.csv',index=False)
  (out/f'{name}_residue_triptych_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
  assert hashlib.sha256(source.read_bytes()).hexdigest()==digest
  print(name,meta)
 (BASE/'PAIRED_TRIPTYCHS.md').write_text('''# 配对三联直方图\n\n直接读取已保存的配对残基指标，没有重算接触或修改底表。\n\n- FoldBench：132个有效配对；主非零图实验132、AF3 131；1个AF3零接触保留在底表及含零补充图。\n- SAFES：仅10个已核对的两蛋白链暂定案例，实验和AF3均10个非零；不能代表完整20案例的正式结合伙伴基准。\n- 横轴线性，纵轴每组观测百分比；两组共享分箱，展示全部范围，没有裁剪或QC删除。\n- 三项指标为8 Å残基对接触数、总界面残基、(IA+IB)/(LA+LB)接口级比例。第三项沿用本轮定义，不是示例旧图将两条链分别汇入的链级分布。\n- 新图位于各组figures目录，名称包含 `residue_triptych_nonzero`；`all_including_zeros` 为含零补充版本。PNG/PDF、分箱CSV和元数据JSON均保存。\n''')
if __name__=='__main__':main()
