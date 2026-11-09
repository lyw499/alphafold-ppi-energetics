"""Checkpointed Pedro residue minima, reusing the frozen benchmark implementation."""
import sys,os,json,hashlib,zipfile,tempfile,shutil,argparse,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/local_interface_geometry/pedro'
sys.path.insert(0,str(OUT/'_deps'))
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from local_geometry_common import pack_chain,min_contacts,summarize,BIN_EDGES,BIN_NAMES,CHEM_CLASSES
from atom_contact_common import parse_structure

def fingerprint():
 return {n:hashlib.sha256((ROOT/'src'/n).read_bytes()).hexdigest() for n in ['local_geometry_common.py','atom_contact_common.py','parse_pedro_local_geometry.py']}
def atomic_json(path,obj):
 tmp=path.with_suffix('.partial');tmp.write_text(json.dumps(obj,indent=2,default=str)+'\n');tmp.replace(path)
def run(index,selected=None):
 protocol=json.loads((OUT/'protocol.json').read_text());assert protocol==fingerprint(),'Implementation changed since protocol freeze'
 if selected is None:
  assert json.loads((OUT/'validation/status.json').read_text())['status']=='passed','Validation gate not passed'
 meta=pd.read_parquet(OUT/f'metadata_parts/archive_{index:03d}.parquet')
 if selected is not None:meta=meta[meta.job_key.isin(selected)]
 groups=[(k,g) for k,g in meta.groupby('job_key',sort=False)]
 pending=[]
 for k,g in groups:
  marker=OUT/f'checkpoints/{k}.done.json'
  if marker.exists():
   done=json.loads(marker.read_text());assert done['protocol']==protocol
   for rel in done['files']:assert (OUT/rel).exists(),f'Missing committed output: {rel}'
  else:pending.append((k,g))
 if not pending:print(f'Archive {index}: all requested jobs already complete',flush=True);return
 member=meta.archive_member.iloc[0]
 scratch=Path(os.environ.get('SLURM_TMPDIR','/tmp'))
 with tempfile.TemporaryDirectory(prefix=f'pedro_geometry_{index}_',dir=scratch) as temp:
  innerpath=Path(temp)/'inner.zip'
  print(f'Archive {index}: copying {member}; pending jobs={len(pending)}',flush=True)
  with zipfile.ZipFile(ROOT/'data/AF3_completed_zips.zip') as outer,outer.open(member) as src,innerpath.open('wb') as dst:shutil.copyfileobj(src,dst,8*1024*1024)
  with zipfile.ZipFile(innerpath) as inner:
   names=inner.namelist()
   for key,g in pending:
    contactpath=OUT/f'contacts/contacts_{key}.parquet';partial=contactpath.with_suffix('.partial');writer=None
    try:
     first=g.iloc[0];modelnum=int(first.best_model)
     candidates=[n for n in names if n.startswith(first.job_id+'/') and n.endswith(f'_model_{modelnum}.cif')]
     assert len(candidates)==1,f'Expected one selected CIF, found {len(candidates)}'
     model,info=parse_structure(cif_text=inner.read(candidates[0]).decode(),author_chains=True)
     packed={};summaries=[];audit=[];nrows=0
     for old in g.to_dict('records'):
      ca,cb=old['chain_a'],old['chain_b']
      for c in [ca,cb]:
       if c not in packed:packed[c]=pack_chain(model[c])
      a,b=packed[ca],packed[cb];rows=min_contacts(a,b);stats=summarize(rows)
      ia=len({r['residue_index_a'] for r in rows});ib=len({r['residue_index_b'] for r in rows})
      globalmin=float(b['tree'].query(a['xyz'])[0].min())
      checks={'contact_count':len(rows)==int(old['contact_count']),'length_a':len(a['res'])==int(old['length_a']),'length_b':len(b['res'])==int(old['length_b']),'interface_residues_a':ia==int(old['interface_residues_a']),'interface_residues_b':ib==int(old['interface_residues_b']),'minimum_distance':abs(globalmin-float(old['min_interchain_distance']))<1e-4}
      status='passed' if all(checks.values()) else 'mismatch'
      context={k:old[k] for k in ['interface_id','job_key','source_file','job_id','protein_pair','protein_a','protein_b','chain_a','chain_b','mean_pair_iptm','best_model_pair_iptm','job_ranking_score','has_clash']}
      context.update(selected_model=modelnum,existing_severe_overlap=bool(old['min_interchain_distance']<.5),parser_warnings=str(info.get('parser_warnings','')),archive=member,cif_member=candidates[0])
      values=np.asarray([r['d_min'] for r in rows]);bins=np.bincount(np.searchsorted(BIN_EDGES,values,side='left'),minlength=5)
      s={**old,**context,**stats,'contact_count':len(rows),'old_contact_count':old['contact_count'],'interface_residues_a':ia,'interface_residues_b':ib,'total_interface_residues':ia+ib,'interface_fraction':(ia+ib)/(len(a['res'])+len(b['res'])),'minimum_interchain_heavy_atom_distance':globalmin,'total_chain_length':len(a['res'])+len(b['res']),'n_deuterium':a['n_deuterium']+b['n_deuterium'],'status':status}
      for name,count in zip(BIN_NAMES,bins):s[name.replace('fraction','count')]=int(count)
      for cls in CHEM_CLASSES:s['count_chem_'+cls]=sum(r['contact_class']==cls for r in rows)
      summaries.append(s)
      audit.append({**context,'number_of_chains':len(list(model)),'old_contact_count':old['contact_count'],'new_contact_count':len(rows),'difference':len(rows)-old['contact_count'],'status':status,'reason':';'.join(k for k,v in checks.items() if not v),'n_deuterium':s['n_deuterium']})
      if rows:
       table=pa.Table.from_pylist([{**context,**r} for r in rows])
       if writer is None:writer=pq.ParquetWriter(partial,table.schema,compression='zstd')
       writer.write_table(table,row_group_size=20000);nrows+=len(rows)
     pd.DataFrame(audit).to_csv(OUT/f'audit_parts/{key}.csv',index=False)
     assert all(s['status']=='passed' for s in summaries),'Legacy count/length/interface/minimum-distance discrepancy; see audit'
     if writer is not None:writer.close();writer=None;partial.replace(contactpath)
     summarypath=OUT/f'summary_parts/{key}.parquet';pd.DataFrame(summaries).to_parquet(summarypath.with_suffix('.partial'),index=False,compression='zstd');summarypath.with_suffix('.partial').replace(summarypath)
     files=[str(summarypath.relative_to(OUT)),f'audit_parts/{key}.csv']+([str(contactpath.relative_to(OUT))] if nrows else [])
     atomic_json(OUT/f'checkpoints/{key}.done.json',dict(job_key=key,protocol=protocol,interfaces=len(g),contacts=nrows,files=files))
     (OUT/f'checkpoints/{key}.failed.json').unlink(missing_ok=True)
     print(f'DONE {key}: interfaces={len(g)} contacts={nrows}',flush=True)
    except Exception as e:
     if writer is not None:writer.close()
     atomic_json(OUT/f'checkpoints/{key}.failed.json',dict(job_key=key,reason=str(e),traceback=traceback.format_exc()))
     if not (OUT/f'audit_parts/{key}.csv').exists():
      fail=g.copy();fail['status']='failed';fail['reason']=str(e);fail.to_csv(OUT/f'audit_parts/{key}.csv',index=False)
     raise
 if selected is None:atomic_json(OUT/f'checkpoints/archive_{index:03d}.done.json',dict(index=index,jobs=len(groups),protocol=protocol))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--index',type=int);p.add_argument('--validation',action='store_true');p.add_argument('--freeze',action='store_true');args=p.parse_args()
 if args.freeze:
  assert not (OUT/'protocol.json').exists();atomic_json(OUT/'protocol.json',fingerprint())
 elif args.validation:
  chosen=pd.read_csv(OUT/'validation/selected_jobs.csv')
  for idx,g in chosen.groupby('archive_index'):run(int(idx),set(g.job_key))
 else:run(args.index)
