"""Conditional secondary analysis, only SC-SC and existing chemistry categories."""
from compute_contact_type_decomposition import OUT,fullframes,classify,BINS,BB
import pandas as pd
import numpy as np
import collections,json

def main():
 assert json.loads((OUT/'classification_validation.json').read_text())['status']=='passed'
 # This task is run only after inspecting primary class contributions.
 counts=collections.Counter();sc_counts=collections.Counter();bbsc_counts=collections.Counter();bbsc_totals=collections.Counter()
 for num,(dataset,d) in enumerate(fullframes(),1):
  d=d.copy();d['contact_type']=classify(d.atom_a_name,d.atom_b_name);v=d.d_min.to_numpy();d['distance_bin']=np.select([v<2,v<=3,v<=4,v<=5,v<=6],BINS[:5],default='6to8')
  versions=[(dataset,d)]
  if dataset=='Mgen Experimental':versions += [('Mgen excl 8PBZ 8PC1',d[~d.pdb_id.isin(['8PBZ','8PC1'])]),('Mgen excl 8PBZ 8PC0 8PC1',d[~d.pdb_id.isin(['8PBZ','8PC0','8PC1'])])]
  for label,g in versions:
   ss=g[g.contact_type=='SC-SC'];sc_counts[label]+=len(ss)
   for (chem,b),n in ss.groupby(['contact_class','distance_bin']).size().items():counts[label,chem,b]+=int(n)
   bs=g[g.contact_type=='BB-SC'].copy();bbsc_totals[label]+=len(bs);bs['sidechain_chemistry']=np.where(bs.atom_a_name.isin(BB),bs.chemistry_b,bs.chemistry_a)
   for (chem,b),n in bs.groupby(['sidechain_chemistry','distance_bin']).size().items():bbsc_counts[label,chem,b]+=int(n)
  if num%300==0:print('SC-SC files',num,flush=True)
 totals=pd.read_csv(OUT/'contact_type_counts.csv').groupby('dataset').total_contacts.first().to_dict();classes=sorted({c for _,c,_ in counts});rows=[]
 for ds,total in totals.items():
  for cls in classes:
   nc=sum(counts[ds,cls,b] for b in BINS)
   for b in BINS:
    n=counts[ds,cls,b];rows.append(dict(dataset=ds,chemistry=cls,distance_bin=b,count=n,fraction_all_contacts=n/total,fraction_all_scsc=n/sc_counts[ds],fraction_within_chemistry=n/nc if nc else np.nan))
 pd.DataFrame(rows).to_csv(OUT/'scsc_chemistry_distance_summary.csv',index=False)
 bbrows=[]
 for ds,total in totals.items():
  for cls in sorted({c for _,c,_ in bbsc_counts}):
   nc=sum(bbsc_counts[ds,cls,b] for b in BINS)
   for b in BINS:
    n=bbsc_counts[ds,cls,b];bbrows.append(dict(dataset=ds,sidechain_chemistry=cls,distance_bin=b,count=n,fraction_all_contacts=n/total,fraction_all_bbsc=n/bbsc_totals[ds],fraction_within_chemistry=n/nc if nc else np.nan))
 pd.DataFrame(bbrows).to_csv(OUT/'bbsc_sidechain_chemistry_distance_summary.csv',index=False)
 print('Focused BB-SC and SC-SC chemistry complete',flush=True)
if __name__=='__main__':main()
