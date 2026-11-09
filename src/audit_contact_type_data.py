"""Read contact tables only: schema, missing-value and atom-name inventory."""
import sys,json,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/local_interface_geometry';OUT=BASE/'contact_type_decomposition';sys.path.insert(0,str(BASE/'pedro/_deps'))
import pandas as pd
import numpy as np
import pyarrow.parquet as pq
REQUIRED=['interface_id','d_min','residue_a_type','residue_b_type','atom_a_name','atom_b_name','atom_a_element','atom_b_element']
def sources():
 for ds in ['foldbench','safes']:yield BASE/ds/'residue_contacts.csv',ds
 yield BASE/'mgen_same_species_control/mgen_experimental_contact_geometry.parquet','Mgen Experimental'
 for p in sorted((BASE/'pedro/contacts').glob('*.parquet')):yield p,'Pedro AF3'
def frames(p,label,columns=None):
 if p.suffix=='.csv':
  for d in pd.read_csv(p,low_memory=False,chunksize=50000,usecols=columns):
   for typ,g in d.groupby('structure_type'):yield ('FoldBench' if label=='foldbench' else 'SAFES')+' '+typ,g
 else:
  for batch in pq.ParquetFile(p).iter_batches(batch_size=50000,columns=columns):yield label,batch.to_pandas()
def main():
 files=[];counts=collections.Counter();missing=collections.Counter();atoms=collections.Counter();samples=[]
 for num,(p,label) in enumerate(sources(),1):
  cols=pq.ParquetFile(p).schema_arrow.names if p.suffix=='.parquet' else pd.read_csv(p,nrows=0).columns.tolist()
  absent=sorted(set(REQUIRED)-set(cols));assert not absent,(str(p),absent)
  rows=0;blank=0;readcols=REQUIRED+(['structure_type'] if p.suffix=='.csv' else [])
  for dataset,g in frames(p,label,readcols):
   rows+=len(g);counts[dataset]+=len(g)
   assert np.isfinite(g.d_min).all() and g.d_min.between(0,8).all()
   for col in REQUIRED:
    n=int((g[col].isna() | g[col].astype(str).str.strip().eq('')).sum());missing[dataset,col]+=n;blank+=n
   for side in ['a','b']:
    vc=g.groupby([f'residue_{side}_type',f'atom_{side}_name',f'atom_{side}_element'],dropna=False).size()
    for (res,atom,element),n in vc.items():atoms[dataset,str(res),str(atom),str(element)]+=int(n)
   if len(samples)<60:samples.extend(g.head(2).assign(dataset=dataset).to_dict('records'))
  files.append(dict(file_path=str(p),size_bytes=p.stat().st_size,n_rows=rows,required_columns_present=True,missing_or_blank_required_cells=blank,columns=json.dumps(cols)))
  if num%200==0:print('Audited files',num,flush=True)
 pd.DataFrame(files).to_csv(OUT/'input_file_audit.csv',index=False)
 pd.DataFrame([dict(dataset=ds,column=col,n_missing=n) for (ds,col),n in missing.items()]).to_csv(OUT/'missing_field_audit.csv',index=False)
 pd.DataFrame([dict(dataset=ds,residue_type=r,atom_name=a,element=e,occurrences=n,rule='BB' if a in ['N','CA','C','O'] else 'SC',special_case=a=='OXT' or r not in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']) for (ds,r,a,e),n in atoms.items()]).to_csv(OUT/'atom_name_inventory.csv',index=False)
 pd.DataFrame(samples).to_csv(OUT/'classification_input_samples.csv',index=False)
 expected={'FoldBench Experimental':35539,'FoldBench AF3':40836,'SAFES Experimental':8922,'SAFES AF3':9883,'Mgen Experimental':10528,'Pedro AF3':6671664}
 assert dict(counts)==expected,(counts,expected)
 assert sum(missing.values())==0,'Missing essential data: stop; do not reparse structures'
 status=dict(status='passed',counts=dict(counts),n_files=len(files),total_size_bytes=sum(f['size_bytes'] for f in files),missing_required_cells=sum(missing.values()),structural_reparse_required=False)
 (OUT/'data_audit.json').write_text(json.dumps(status,indent=2));print(json.dumps(status,indent=2),flush=True)
if __name__=='__main__':main()
