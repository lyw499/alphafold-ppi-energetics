"""Reconcile saved geometry with original confidence and QC tables; no coordinates."""
import sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/local_interface_geometry';OUT=BASE/'pedro_confidence_analysis';sys.path.insert(0,str(BASE/'pedro/_deps'))
import numpy as np
import pandas as pd
KEY=['source_file','job_id','protein_pair','best_model']
def main():
 path=BASE/'pedro/pedro_interface_geometry_summary.parquet';s=pd.read_parquet(path);hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()}
 path=ROOT/'results/best_model_pair_iptm_raw.csv';c=pd.read_csv(path);hashes[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest();c.source_file=c.source_file.str.replace(r'\.zip$','.csv',regex=True)
 dc=int(c.duplicated(KEY).sum());ds=int(s.duplicated(KEY).sum());assert dc==ds==0
 join=s.merge(c[KEY+['best_model_pair_iptm','best_model_ranking_score']],on=KEY,how='left',validate='one_to_one',indicator=True,suffixes=('','_source'))
 assert join._merge.eq('both').all() and len(join)==len(s)
 assert np.allclose(join.best_model_pair_iptm,join.best_model_pair_iptm_source,rtol=0,atol=1e-12)
 assert np.allclose(join.job_ranking_score,join.best_model_ranking_score,rtol=0,atol=1e-10)
 parts=[]
 for source in s.source_file.unique():
  path=ROOT/'results/pooled_parts_8A'/source;hashes[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
  d=pd.read_csv(path,usecols=['job_id','protein_pair','best_model','chain_a','chain_b','contact_count','min_interchain_distance','has_clash','mean_pair_iptm']);d['source_file']=source;parts.append(d)
 q=pd.concat(parts,ignore_index=True);assert not q.duplicated(KEY).any()
 j=s.merge(q,on=KEY,how='left',validate='one_to_one',indicator=True,suffixes=('','_source'))
 assert j._merge.eq('both').all()
 for col in ['chain_a','chain_b','contact_count','has_clash']:assert j[col].eq(j[col+'_source']).all(),col
 assert np.allclose(j.mean_pair_iptm,j.mean_pair_iptm_source,rtol=0,atol=1e-12)
 assert j.existing_severe_overlap.eq(j.min_interchain_distance_source.lt(.5)).all()
 s['fraction_lt_2']=np.where(s.n_contacts>0,s.n_lt_2_0/s.n_contacts,np.nan)
 s['geometry_defined']=s.n_contacts>0
 assert s.loc[~s.geometry_defined,['fraction_3_4','median_d_min','fraction_lt_2']].isna().all().all()
 s.to_csv(OUT/'pedro_confidence_interface_table.csv',index=False)
 audit=dict(status='passed',join_keys=KEY,interfaces_before=len(s),confidence_matched=len(join),confidence_unmatched=int(join._merge.ne('both').sum()),summary_duplicate_keys=ds,confidence_duplicate_keys=dc,qc_matched=len(j),qc_unmatched=int(j._merge.ne('both').sum()),interface_id_duplicates=int(s.interface_id.duplicated().sum()),nonzero=int(s.geometry_defined.sum()),zero=int((~s.geometry_defined).sum()),severe=int(s.existing_severe_overlap.sum()),confidence_values_exactly_reconciled=True,qc_flags_reconciled=True,primary_metric='best_model_pair_iptm',available_metrics=['best_model_pair_iptm','mean_pair_iptm','job_ranking_score'],ptm_available=False,selected_model_counts={str(k):int(v) for k,v in s.selected_model.value_counts().items()})
 (OUT/'join_audit.json').write_text(json.dumps(audit,indent=2));(OUT/'input_sha256.json').write_text(json.dumps(hashes,indent=2))
 fields=[]
 for col in audit['available_metrics']:
  v=s[col];fields.append(dict(metric=col,n_present=int(v.notna().sum()),n_missing=int(v.isna().sum()),n_unique=v.nunique(),min=v.min(),max=v.max(),median=v.median()))
 pd.DataFrame(fields).to_csv(OUT/'confidence_field_inventory.csv',index=False)
 print(json.dumps(audit,indent=2));print(pd.DataFrame(fields).to_string(index=False))
if __name__=='__main__':main()
