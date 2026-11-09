"""Existing audited reference coordinates vs saved Pedro outputs; never opens Pedro archives."""
import sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/local_interface_geometry';OUT=BASE/'mgen_same_species_control';sys.path.insert(0,str(BASE/'pedro/_deps'))
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from atom_contact_common import parse_structure
from local_geometry_common import pack_chain,min_contacts,summarize,BIN_NAMES,CHEM_CLASSES
EDGES=np.arange(0,8.00001,.25);TH=[.5,1.,1.5,2.];hashes={}
def record(p):hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
def save(fig,name):
 fig.tight_layout();fig.savefig(OUT/'figures'/name,dpi=180);plt.close(fig)
def main():
 path=ROOT/'results/contact_comparison_residue8A/reference_audit/reference_interface_audit.csv';record(path);audit=pd.read_csv(path)
 assert len(audit)==29 and audit.pair_species_class.eq('Mgen-Mgen').sum()==23
 contacts=[];sums=[];checks=[]
 for source,g in audit.groupby('source'):
  path=ROOT/source;record(path);model,info=parse_structure(path);cache={}
  for old in g.to_dict('records'):
   for c in [old['chain1'],old['chain2']]:
    if c not in cache:cache[c]=pack_chain(model[c])
   a,b=cache[old['chain1']],cache[old['chain2']];rows=min_contacts(a,b);stats=summarize(rows)
   ia=len({r['residue_index_a'] for r in rows});ib=len({r['residue_index_b'] for r in rows})
   assert len(rows)==old['contact_count'] and ia==old['interface_residues_chain1'] and ib==old['interface_residues_chain2']
   assert len(a['res'])==old['chain1_length'] and len(b['res'])==old['chain2_length']
   meta={k:old[k] for k in ['interface_id','structure_id','pdb_id','biological_assembly_id','chain1','chain2','pair_species_class','existing_qc_flags','resolution_A','protein_name_A','protein_name_B']};meta['parser_warnings']=info['parser_warnings']
   contacts.extend({**meta,**r} for r in rows)
   stats.update({f'fraction_lt_{str(t).replace(".","_")}':stats['n_lt_'+str(t).replace('.','_')]/len(rows) if rows else np.nan for t in TH})
   sums.append({**old,**meta,**stats});checks.append(dict(interface_id=old['interface_id'],old_count=old['contact_count'],new_count=len(rows),status='passed'))
 c=pd.DataFrame(contacts);s=pd.DataFrame(sums);pd.DataFrame(checks).to_csv(OUT/'reference_count_validation.csv',index=False)
 c.to_parquet(OUT/'all_reference_contact_geometry.parquet',index=False,compression='zstd');s.to_csv(OUT/'all_reference_interface_summary.csv',index=False)
 m=s[s.pair_species_class.eq('Mgen-Mgen')].copy();mc=c[c.pair_species_class.eq('Mgen-Mgen')].copy()
 mc.to_parquet(OUT/'mgen_experimental_contact_geometry.parquet',index=False,compression='zstd');m.to_csv(OUT/'mgen_experimental_interface_summary.csv',index=False)
 c[c.d_min<2].sort_values('d_min').to_csv(OUT/'reference_short_distance_audit.csv',index=False)
 # Existing benchmark outputs only; no mappings or benchmark coordinates are recomputed.
 groups={};summary={}
 for ds,label in [('foldbench','FoldBench'),('safes','SAFES')]:
  for file in ['residue_contacts.csv','interface_summaries.csv']:record(BASE/ds/file)
  bc=pd.read_csv(BASE/ds/'residue_contacts.csv');bs=pd.read_csv(BASE/ds/'interface_summaries.csv')
  for typ in ['Experimental','AF3']:
   name=label+' '+typ;groups[name]=(bc.loc[bc.structure_type==typ,'d_min'].to_numpy(),None);summary[name]=bs[bs.structure_type==typ].copy()
 path=BASE/'pedro/pedro_interface_geometry_summary.parquet';record(path);ped=pd.read_parquet(path)
 path=BASE/'pedro/qc/pooled_distances.npy';record(path);pv=np.load(path,mmap_mode='r');assert len(pv)==int(ped.n_contacts.sum())
 groups['Pedro AF3']=(pv,None);summary['Pedro AF3']=ped
 groups['Mgen Experimental']=(mc.d_min.to_numpy(),None);summary['Mgen Experimental']=m
 views={'all_Mgen':m,'exclude_8PBZ_8PC1':m[~m.pdb_id.isin(['8PBZ','8PC1'])],'exclude_8PBZ_8PC0_8PC1':m[~m.pdb_id.isin(['8PBZ','8PC0','8PC1'])]}
 # Each nonzero assembly (or PDB) receives equal total contact weight. Zero groups have no conditional contact spectrum.
 sensitivity=[];senscurves={};assembly=[]
 for view,frame in views.items():
  cc=mc[mc.interface_id.isin(frame.interface_id)]
  for unit in ['pooled','assembly_equal','PDB_equal']:
   w=None
   if unit!='pooled':
    key='structure_id' if unit=='assembly_equal' else 'pdb_id';size=cc.groupby(key).d_min.transform('size');w=1/size.to_numpy();w=w/w.sum()
   name=view+' / '+unit;senscurves[name]=(cc.d_min.to_numpy(),w)
   sensitivity.append(dict(view=view,weighting=unit,n_pairs=len(frame),n_zero=int(frame.n_contacts.eq(0).sum()),n_nonzero=int(frame.n_contacts.gt(0).sum()),n_pdb=frame.pdb_id.nunique(),n_assemblies=frame.structure_id.nunique(),**distribution(cc.d_min.to_numpy(),w)))
 for sid,g in m.groupby('structure_id'):
  cc=mc[mc.structure_id==sid];assembly.append(dict(structure_id=sid,pdb_id=g.pdb_id.iloc[0],resolution_A=g.resolution_A.iloc[0],pairs=len(g),nonzero_pairs=int(g.n_contacts.gt(0).sum()),**distribution(cc.d_min.to_numpy(),None)))
 pd.DataFrame(sensitivity).to_csv(OUT/'mgen_reference_sensitivity_summary.csv',index=False);pd.DataFrame(assembly).to_csv(OUT/'assembly_level_geometry.csv',index=False)
 pooled=[];ilevel=[];chem=[];qc=[]
 for name,(v,w) in groups.items():
  frame=summary[name];pooled.append(dict(dataset=name,n_interfaces=len(frame),n_zero=int(frame.n_contacts.eq(0).sum()),**distribution(v,w)))
  for metric in ['median_d_min']+BIN_NAMES+['fraction_lt_2_0','n_contacts','total_interface_residues','interface_fraction']:
   if metric=='fraction_lt_2_0':z=frame.n_lt_2_0/frame.n_contacts.replace(0,np.nan)
   elif metric in frame:z=frame[metric]
   else:continue
   z=z.dropna();ilevel.append(dict(dataset=name,metric=metric,n_defined=len(z),mean=z.mean(),median=z.median(),q25=z.quantile(.25),q75=z.quantile(.75)))
  for cls in CHEM_CLASSES:
   chem.append(dict(dataset=name,contact_class=cls,pooled_fraction=float((frame['chem_'+cls].fillna(0)*frame.n_contacts).sum()/frame.n_contacts.sum()),mean_interface_fraction=frame['chem_'+cls].mean()))
  for t in TH:
   n=frame['n_lt_'+str(t).replace('.','_')];qc.append(dict(dataset=name,threshold=t,contacts=int(n.sum()),contact_fraction=float(n.sum()/frame.n_contacts.sum()),affected_interfaces=int(n.gt(0).sum()),fraction_all_interfaces=n.gt(0).mean(),fraction_nonzero_interfaces=n.gt(0).sum()/frame.n_contacts.gt(0).sum()))
 pd.DataFrame(pooled).to_csv(OUT/'mgen_same_species_comparison_summary.csv',index=False);pd.DataFrame(ilevel).to_csv(OUT/'interface_level_distributions.csv',index=False);pd.DataFrame(chem).to_csv(OUT/'chemistry_comparison.csv',index=False);pd.DataFrame(qc).to_csv(OUT/'short_distance_comparison.csv',index=False)
 # Retain other species as clearly separated audit curves, never in primary Mgen curve.
 categories=[]
 for cat,g in s.groupby('pair_species_class'):
  cc=c[c.pair_species_class==cat];categories.append(dict(category=cat,n_pairs=len(g),n_zero=int(g.n_contacts.eq(0).sum()),**distribution(cc.d_min.to_numpy(),None)))
 pd.DataFrame(categories).to_csv(OUT/'species_category_sensitivity.csv',index=False)
 def plotdist(ax,items,ecdf=False):
  for label,(x,w) in items.items():
   if ecdf:
    order=np.argsort(x);xx=np.asarray(x)[order];yy=np.arange(1,len(x)+1)/len(x) if w is None else np.cumsum(w[order])/np.sum(w);ix=np.unique(np.linspace(0,len(x)-1,min(10000,len(x)),dtype=int));ax.plot(xx[ix],yy[ix],label=label)
   else:
    h,_=np.histogram(x,EDGES,weights=w);ax.stairs(h/h.sum()/.25,EDGES,label=label,lw=1.6)
  ax.set(xlim=(0,8),xlabel='Residue-pair minimum heavy-atom distance (Å)',ylabel='ECDF' if ecdf else 'Density (Å⁻¹)');ax.legend(fontsize=7)
 for ecdf,suffix in [(False,'distribution'),(True,'ecdf')]:
  fig,ax=plt.subplots(figsize=(8,4));plotdist(ax,{n:groups[n] for n in ['Mgen Experimental','Pedro AF3']},ecdf);ax.set_title('Descriptive pooled contacts; contacts are not independent');save(fig,f'mgen_exp_vs_pedro_distance_{suffix}.png')
  fig,axes=plt.subplots(1,3,figsize=(17,4),sharey=True)
  for ax,names in zip(axes,[['FoldBench Experimental','FoldBench AF3'],['SAFES Experimental','SAFES AF3'],['Mgen Experimental','Pedro AF3']]):plotdist(ax,{n:groups[n] for n in names},ecdf)
  fig.suptitle('Paired benchmarks | Same-species unpaired control — pooled contacts are descriptive')
  save(fig,'all_benchmarks_three_panel_distance_'+('ecdf' if ecdf else 'comparison')+'.png')
 for name,keys in [('redundancy',['all_Mgen / pooled','all_Mgen / assembly_equal','all_Mgen / PDB_equal']),('resolution',['all_Mgen / pooled','exclude_8PBZ_8PC1 / pooled','exclude_8PBZ_8PC0_8PC1 / pooled'])]:
  fig,ax=plt.subplots(figsize=(9,4));plotdist(ax,{k:senscurves[k] for k in keys});save(fig,f'mgen_reference_{name}_sensitivity.png')
 fig,axes=plt.subplots(1,4,figsize=(14,4))
 for ax,b in zip(axes,BIN_NAMES[1:]):
  ax.boxplot([m[b].dropna(),ped[b].dropna()],labels=['Mgen Exp','Pedro AF3'],showfliers=False);ax.scatter(np.ones(m[b].notna().sum()),m[b].dropna(),s=12,alpha=.6);ax.set(title=b.replace('fraction_',''),ylabel='Within-interface fraction',ylim=(0,1))
 save(fig,'mgen_exp_vs_pedro_interface_level_bins.png')
 fig,ax=plt.subplots(figsize=(8,4))
 for label,frame in [('Mgen Experimental',m),('Pedro AF3',ped)]:
  x=np.sort(frame.median_d_min.dropna());ii=np.unique(np.linspace(0,len(x)-1,min(10000,len(x)),dtype=int));ax.plot(x[ii],(ii+1)/len(x),label=label)
 ax.set(xlabel='Within-interface median distance (Å)',ylabel='ECDF',xlim=(0,8));ax.legend();save(fig,'interface_median_ecdf.png')
 for p,h in hashes.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
 (OUT/'provenance.json').write_text(json.dumps(dict(input_sha256=hashes,Pedro_coordinates_accessed=False,benchmark_recomputed=False,reference_pairs_validated=len(s),code_sha256={f:hashlib.sha256((ROOT/'src'/f).read_bytes()).hexdigest() for f in ['local_geometry_common.py','atom_contact_common.py','mgen_same_species_control.py']}),indent=2))
 print(pd.DataFrame(pooled).to_string(index=False));print(pd.DataFrame(sensitivity).to_string(index=False));print(pd.DataFrame(assembly).to_string(index=False))
 text='# M. genitalium 同种实验对照\n\n仅重新计算已下载并审计的 29 个参考链对；全部与旧接触数、链长和界面残基数精确一致。主结果只含 23 个 Mgen–Mgen 链对；其他物种单独保存。Pedro 和 FoldBench/SAFES 仅读取既有结果，输入文件哈希核对未变。\n\n'
 text+='定义：重原子最小距离 ≤8 Å 的唯一残基对。复用原解析器、标签链 ID、首个坐标模型、默认替代构象；排除 H/D，标准蛋白残基处理不变。分箱为 [0,3]、(3,4]、(4,5]、(5,6]、(6,8] Å。零接触界面保留，条件距离与比例未定义记 NA。\n\n'
 for title,table in [('主要池化比较',pooled),('敏感性分析',sensitivity),('逐装配结果',assembly),('逐界面分布',ilevel),('短距离诊断',qc),('化学组成（探索性）',chem)]:text+='## '+title+'\n\n```\n'+pd.DataFrame(table).to_string(index=False)+'\n```\n\n'
 text+='## 敏感性方法及限制\n\n装配等权先在每个生物学装配内汇总接触，再使每个非零装配的总权重相同；PDB 等权以同一规则合并同一 PDB 的装配。没有删除或任意挑选链对。零接触装配没有条件距离分布，不赋予虚构距离。装配等权仍不能消除不同 PDB 来自同一蛋白复合物/研究的相关性。\n\n低分辨率敏感性明确排除既有审计中的 8PBZ（11 Å，短距离几何疑点）和 8PC1（18 Å）；扩展视图另排除同一审计中的 8PC0（17 Å）。这不是新增分辨率阈值或更改原始 QC。原始结果全部保留。\n\n组间不配对，不做残基级显著性检验。ECDF 显示至多 10,000 个经验秩点；直方图为 0.25 Å 同一分箱。化学分类复用已有宽泛类别。短距离不自动称作碰撞。\n'
 (OUT/'MGEN_SAME_SPECIES_CONTROL.md').write_text(text)
 (OUT/'completion.json').write_text(json.dumps(dict(status='complete',reference_pairs=29,mgen_pairs=len(m),mgen_nonzero=int(m.n_contacts.gt(0).sum()),mgen_pdb=m.pdb_id.nunique(),mgen_assemblies=m.structure_id.nunique())))
def distribution(x,w):
 x=np.asarray(x);n=len(x)
 if not n:return dict(n_contacts=0)
 if w is None:w=np.full(n,1/n)
 w=w/w.sum();order=np.argsort(x);cum=np.cumsum(w[order]);weighted=lambda q:float(x[order][min(np.searchsorted(cum,q),n-1)])
 result=dict(n_contacts=n,mean=float(np.sum(w*x)),median=weighted(.5),q25=weighted(.25),q75=weighted(.75))
 idx=np.searchsorted([3,4,5,6,8],x,side='left')
 for i,name in enumerate(BIN_NAMES):result[name]=float(w[idx==i].sum())
 for t in TH:result['fraction_lt_'+str(t).replace('.','_')]=float(w[x<t].sum())
 return result
if __name__=='__main__':main()
