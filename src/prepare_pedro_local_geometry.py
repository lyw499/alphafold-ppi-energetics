"""Inventory existing Pedro metadata; new files only under pedro local geometry."""
import sys,json,hashlib,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/local_interface_geometry/pedro';sys.path.insert(0,str(OUT/'_deps'))
import pandas as pd
import pyarrow

def main():
 dest=OUT/'metadata_parts';dest.mkdir(exist_ok=True)
 manifest=json.loads((ROOT/'results/contact_comparison_atom8A/pedro_manifest.json').read_text())
 conf=pd.read_csv(ROOT/'results/best_model_pair_iptm_raw.csv',usecols=['source_file','job_id','protein_pair','best_model','best_model_pair_iptm'])
 conf.source_file=conf.source_file.str.replace(r'\.zip$','.csv',regex=True)
 keys=['source_file','job_id','protein_pair','best_model'];assert not conf.duplicated(keys).any()
 archives=[];jobs=[];candidates=[];hashes={}
 with zipfile.ZipFile(ROOT/'data/AF3_completed_zips.zip') as z:
  for spec in manifest:
   path=ROOT/'results/pooled_parts_8A'/spec['source_file'];hashes[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
   d=pd.read_csv(path);d['source_file']=spec['source_file'];d['archive_index']=spec['index'];d['archive_member']=spec['member'];d['source_row']=range(len(d));d=d.merge(conf,on=keys,how='left',validate='one_to_one',sort=False)
   assert d.contact_cutoff_A.eq(8).all()
   for j,((job,model),g) in enumerate(d.groupby(['job_id','best_model'],sort=False)):
    key=f'a{spec["index"]:03d}_j{j:04d}'
    d.loc[g.index,'job_key']=key;d.loc[g.index,'interface_id']=[key+':r'+str(i) for i in g.source_row]
    jobs.append(dict(job_key=key,archive_index=spec['index'],archive_member=spec['member'],source_file=spec['source_file'],job_id=job,selected_model=int(model),n_interfaces=len(g),n_zero=int(g.contact_count.eq(0).sum()),n_severe=int(g.min_interchain_distance.lt(.5).sum()),max_contacts=int(g.contact_count.max()),min_positive_length_product=int((g[g.contact_count>0].length_a*g[g.contact_count>0].length_b).min()) if g.contact_count.gt(0).any() else None))
   d.to_parquet(dest/f'archive_{spec["index"]:03d}.parquet',index=False,compression='zstd')
   info=z.getinfo(spec['member']);archives.append(dict(**spec,compressed_bytes=info.compress_size,uncompressed_bytes=info.file_size,n_interfaces=len(d),n_jobs=d.job_key.nunique(),n_zero=int(d.contact_count.eq(0).sum()),n_severe=int(d.min_interchain_distance.lt(.5).sum()),sum_residue_contacts=int(d.contact_count.sum()),n_missing_best_model_iptm=int(d.best_model_pair_iptm.isna().sum())))
 a=pd.DataFrame(archives);j=pd.DataFrame(jobs);a.to_csv(OUT/'archive_manifest.csv',index=False);j.to_csv(OUT/'job_manifest.csv',index=False)
 assert a.n_interfaces.sum()==170644 and len(a)==70 and j.job_key.is_unique
 selection=j.sort_values('max_contacts').tail(1).job_key.tolist()+j.dropna(subset=['min_positive_length_product']).sort_values('min_positive_length_product').head(1).job_key.tolist()+j.head(1).job_key.tolist()+j[j.selected_model!=0].head(1).job_key.tolist()
 selected=j[j.job_key.isin(selection)];selected.to_csv(OUT/'validation/selected_jobs.csv',index=False)
 (OUT/'metadata_input_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
 print(a.drop(columns=['member','source_file']).sum(numeric_only=True).to_string());print('Selected validation jobs:');print(selected.to_string(index=False));print('Parquet pyarrow',pyarrow.__version__)
if __name__=='__main__':main()
