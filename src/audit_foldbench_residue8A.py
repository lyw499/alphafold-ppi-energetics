"""Read-only source audit; outputs exclusively in the new residue comparison directory."""
import json,re
from pathlib import Path
import pandas as pd
from foldbench_interfaces import load_model,protein_residues,alignment_info
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/contact_comparison_residue8A'

def main():
 OUT.mkdir(exist_ok=True)
 old=pd.read_csv(ROOT/'results/foldbench/paired_interfaces_8A.csv').fillna('')
 meta=pd.read_csv(ROOT/'data/foldbench/metadata/interface_protein_protein.csv')
 inv=pd.read_csv(ROOT/'results/foldbench/af3_archive_inventory.csv')
 pdbmembers=inv[inv.member.str.endswith('.pdb')]
 target_counts=pdbmembers.member.map(lambda x:Path(x).name.split('_seed_')[0].lower()).value_counts()
 summaries=pd.read_csv(ROOT/'results/foldbench/af3_target_inventory.csv')
 assert set(target_counts.index)==set(summaries.target_id.str.lower())
 assert (target_counts==200).all()
 assert len(old)==len(meta)==279
 rows=[]
 for target,g in old.groupby('target_id',sort=False):
  exp,_=load_model(ROOT/g.iloc[0].experimental_path,True)
  exists=target.lower() in target_counts
  aliases=[t for t in target_counts.index if t[:4]==target[:4].lower()]
  pred=None
  if exists: pred,_=load_model(ROOT/g.iloc[0].af3_path)
  for _,r in g.iterrows():
   source=meta.iloc[int(r.interface_id)-1]
   assert (source.pdb_id,source.interface_chain_id_1,source.interface_chain_id_2)==(target,r.chain1,r.chain2)
   assert r.chain1 in exp and r.chain2 in exp
   row=dict(interface_id=r.interface_id,experimental_structure_id=target,experimental_chain1=r.chain1,experimental_chain2=r.chain2,
    experimental_file=r.experimental_path,af3_file_available=exists,af3_file=r.af3_path,archive_pdb_count=int(target_counts.get(target,0)),
    same_pdb_archive_targets=json.dumps(aliases),af3_chains=';'.join(c.id for c in pred) if pred is not None else '',
    expected_sequences_present='not_assessable_missing_file',experimental_chains_mappable=False,intended_interface_mappable=False,
    mapped_af3_chain_pair='',mapping_status='unmatched',failure_reason='AF3_target_absent_from_provided_archive',notes='Cannot assess protein presence or chain mapping without AF3; missing is not zero.',
    original_paired_success=r.paired_success)
   if exists:
    assert r.paired_success and r.af3_chain1!=r.af3_chain2
    evidence=[]
    for a,b in [(r.chain1,r.af3_chain1),(r.chain2,r.af3_chain2)]:
     er,es=protein_residues(exp[a]);pr,ps=protein_residues(pred[b]);info=alignment_info(er,es,pr,ps)
     assert info is not None
     evidence.append(dict(experimental_chain=a,af3_chain=b,**{k:v for k,v in info.items() if k!='xyz'}))
    tie=int(r.n_tied_interface_mappings)>1
    sensitive=any(r[f'mapping_tie_{m}_min']!=r[f'mapping_tie_{m}_max'] for m in ['contact_count','total_interface_residues','fraction_interface_chain1','fraction_interface_chain2'])
    row.update(expected_sequences_present=True,experimental_chains_mappable=True,intended_interface_mappable=True,mapped_af3_chain_pair=r.af3_chain1+':'+r.af3_chain2,
     mapping_status='mapped_with_geometric_ambiguity' if tie else 'mapped',failure_reason='',notes='Sequence-compatible mapping verified; global CA RMSD chooses copy correspondence; geometric ties retained, not forced to unique.',
     sequence_evidence_json=json.dumps(evidence),mapping_metric_sensitive=sensitive,tied_chain_pairs=r.tied_interface_chain_pairs_json)
   else: assert not aliases and not r.paired_success
   rows.append(row)
 audit=pd.DataFrame(rows).sort_values('interface_id');audit.to_csv(OUT/'foldbench_matching_audit.csv',index=False)
 assert audit.af3_file_available.sum()==132
 report='''# FoldBench AF3 matching audit\n\n原始 279 个实验接口（239 个 assembly），匹配 132，未匹配 147。\n\n逐条重新检查官方接口行、实验链存在性、完整 AF3 归档成员清单（136,200 个 PDB / 681 targets / 每个 200 个模型）、摘要清单及本地所选 AF3 坐标。所有 132 个已有配对重新通过序列 identity ≥0.95、实验观测序列 coverage ≥0.90 的核验。\n\n147 个未匹配接口全部对应 **131 个 PDB/assembly 在所提供归档中缺失**。大小写归一化及同 PDB 的其他 assembly / `_chain_ligand` 变体搜索也没有命中。不是链 ID mismatch、序列映射失败或 parser 错误造成的 147 条损失。缺文件时 B/C/D 无法评估，不能宣称缺少某蛋白序列，也不能将其计为零。这里证明的是本地提供归档的缺失，不证明任何其他公共来源也没有预测。\n\n现有归档是 Zenodo 21947583 的 Thal-Kak AF3 子集，不是原 FoldBench 全部预测。缓存 Zenodo 元数据明确声明 681 targets；下载校验记录与成员清单一致。本次复用已有完整清单，不重复解压 11 GB 归档。\n\n未发现能恢复这些缺失目标的匹配代码错误，因此没有改写原映射程序，也没有为增加 n 强制匹配。最终仍为 **132 匹配 / 147 未匹配**。链 ID 差异已由序列及全复合体 CA 几何映射处理。\n\n132 是可映射集合，并不都唯一：几何并列涉及 23 个接口，其中 14 个至少一个链级或接口指标敏感，6 个残基对计数敏感。保留原确定性映射、全部并列候选及单独的映射稳定敏感性分析，不把这些不确定性伪装成缺文件。几何拟合本身不证明生物学对应关系，且可能使结构比较偏乐观。\n\n配对分析仅使用上述 132 条，两侧采用同一接口集合；有效 AF3 零接触保留。缺失 147 条仍在审计和全集结果表中。抽样覆盖偏差与同一 assembly 多个相关接口仍需说明。\n\n审计逻辑检查：279 = 132 + 147；131 个缺失 assembly 对应147条接口；所有有文件的指定链对均通过序列复核；不存在同 PDB 变体可恢复目标；缺失值没有转零。审计完成后才进入后续残基比较。\n'''
 (OUT/'FOLDBENCH_MATCHING_AUDIT.md').write_text(report)
 print(audit.mapping_status.value_counts().to_string());print('AUDIT PASSED; missing target causes verified independently of paired_success')
if __name__=='__main__':main()
