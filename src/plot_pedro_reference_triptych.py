"""Replot existing residue tables as percent-weighted linear step histograms."""
import os,json
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/pedro-reference-triptych-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/contact_comparison_residue8A'
OUT=BASE/'figures/pedro_reference'
def main():
 p=pd.read_csv(BASE/'pedro_raw.csv');r=pd.read_csv(BASE/'reference_all_29_candidates.csv')
 r=r[r.reference_taxonomy_scope!='neither_Mgen']
 records=[];metadata=[]
 for version,raw in [('raw',p),('qc',p[p.qc_included])]:
  pp=raw[raw.contact_count>0];rr=r[r.contact_count>0]
  for view in ['focused','full_range']:
   fig,axes=plt.subplots(1,3,figsize=(18,5.4))
   for ax,metric,label,step,focus in zip(axes,['contact_count','total_interface_residues','interface_fraction'],['Residue-pair contacts (8 Å)','Total interface residues','Interface fraction: (IA + IB) / (LA + LB)'],[50,20,.05],[1200,500,1]):
    groups=[(rr,'Experimental reference','#1f77b4'),(pp,'Pedro AF3','#ff7f0e')]
    maximum=max(d[metric].max() for d,_,_ in groups)
    upper=max(step,np.ceil(maximum/step)*step) if metric!='interface_fraction' else 1
    bins=np.linspace(0,upper,int(round(upper/step))+1)
    for d,name,color in groups:
     counts,edges=np.histogram(d[metric],bins=bins)
     percent=counts/len(d)*100
     assert counts.sum()==len(d) and np.isclose(percent.sum(),100)
     ax.stairs(percent,edges,color=color,linewidth=2,label=name)
     beyond=int((d[metric]>focus).sum())
     if view=='focused' and metric!='interface_fraction':
      idx=0 if name=='Experimental reference' else 1
      ax.text(.98,.70-.075*idx,f'{name} > {focus}: {beyond:,} ({100*beyond/len(d):.2f}%)',transform=ax.transAxes,ha='right',fontsize=9,color=color)
     for left,right,n,v in zip(edges[:-1],edges[1:],counts,percent):records.append(dict(version=version,view=view,group=name,metric=metric,bin_left=left,bin_right=right,n=int(n),percent=v,n_denominator=len(d)))
    ax.set(xlabel=label,ylabel='Observations (%)',xlim=(0,focus if view=='focused' else upper));ax.legend(frameon=False,fontsize=10,loc='upper right');ax.spines[['top','right']].set_visible(False)
   qc_note='raw; QC-flagged cases retained' if version=='raw' else 'QC: severe overlaps excluded'
   fig.suptitle(f'Pedro AF3 vs experimental reference candidates; non-zero interfaces; {qc_note}\nPedro n={len(pp):,}; reference n={len(rr)}; distribution-level comparison (not paired)',fontsize=13)
   footer=f'Zeros retained in data, omitted here: Pedro {int((raw.contact_count==0).sum()):,}/{len(raw):,}; reference {int((r.contact_count==0).sum())}/{len(r)}. '
   footer+='Percentages use all non-zero interfaces; '+('linear display limits only, tails reported above.' if view=='focused' else 'full value ranges shown.')
   fig.text(.5,.015,footer,ha='center',fontsize=9)
   fig.tight_layout(rect=(0,.06,1,.88))
   name=f'pedro_vs_reference_residue_triptych_{version}_{view}'
   for ext in ['png','pdf']:fig.savefig(OUT/f'{name}.{ext}',dpi=180)
   plt.close(fig)
   metadata.append(dict(version=version,view=view,pedro_nonzero=len(pp),reference_nonzero=len(rr),fraction_definition='sum of interface residues divided by sum of coordinate chain lengths; one observation per interface',reference_scope='at least one Mgen chain; 27 computable, 25 nonzero',basename=name))
 pd.DataFrame(records).to_csv(OUT/'pedro_vs_reference_triptych_bins.csv',index=False)
 (OUT/'pedro_vs_reference_triptych_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
 (OUT/'TRIPTYCH_README.md').write_text('''# 三联线性百分比直方图\n\n直接复用已有残基表，不重算接触。每组每幅图的所有分箱百分比合计100%。三项指标为唯一残基对数、总界面残基数、接口级长度加权界面比例。第三项沿用本轮统一定义，每个接口贡献一个比例；不切换成示例旧图中两条链各贡献一次的链级比例。\n\nraw 与 QC 各有 focused 和 full_range PNG/PDF。focused 仅限制显示范围（接触数0–1200、界面残基0–500），越界数量/比例标在图中，分母仍为全部非零样本；full_range 展示全部尾部。界面比例为0–1。\n\n参考范围为至少一条链来自Mgen的27条候选，其中25条非零；不是逐一验证的独立界面基准。Pedro使用结构链对观测，包含重跑。零值保留在底表，raw不排除QC标记。图中不作配对推断。\n\n绘图分箱及百分比保存在 `pedro_vs_reference_triptych_bins.csv`。\n''')
 print(json.dumps(metadata,indent=2))
if __name__=='__main__':main()
