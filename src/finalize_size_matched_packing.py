"""Final numerical audit: stabilize weighted quantiles at exact CDF boundaries."""
from size_matched_packing import OUT, stat
import pandas as pd
import numpy as np
import json
def main():
    d=pd.read_parquet(OUT/'analysis_interface_table.parquet');d=d[d.n_contacts>0]
    links=pd.read_csv(OUT/'matched_interface_ids.csv',dtype={'experimental_interface_id':str,'pedro_interface_id':str})
    keys=['reference','size_variable','pedro_version','nearest_k']
    files=['matching_balance_summary.csv','size_matched_close_contact_summary.csv','size_matched_distance_bin_summary.csv']
    tables={f:pd.read_csv(OUT/f) for f in files};gap=[];changes=[]
    for vals,l in links.groupby(keys):
        ref,size,version,k=vals;p=d[d.dataset=='Pedro AF3'].copy()
        if version!='raw':p=p[~p.existing_severe_overlap.eq(True)]
        weights=l.groupby('pedro_interface_id').target_weight.sum().reindex(p.interface_id).fillna(0).to_numpy()
        gap.append(dict(zip(keys,vals),n_targets=l.experimental_interface_id.nunique(),n_exact_targets=l[l.absolute_transformed_gap<1e-12].experimental_interface_id.nunique(),max_absolute_size_gap=float((l.pedro_size-l.experimental_size).abs().max()),max_relative_size_gap=float(((l.pedro_size-l.experimental_size).abs()/l.experimental_size).max()),mean_absolute_transformed_gap=float(np.average(l.absolute_transformed_gap,weights=l.target_weight))))
        for f,t in tables.items():
            mask=np.ones(len(t),bool)
            for key,v in zip(keys,vals):mask &= t[key].eq(v)
            mask &= t.phase.eq('Pedro_after' if f.startswith('matching') else 'Pedro_matched')
            for i,r in t[mask].iterrows():
                metric=size if f.startswith('matching') else (r.metric if 'metric' in r.index else r.contact_type+'_within_'+r.distance_bin)
                z=stat(p[metric],weights)
                assert np.isclose(z['mean'],r['mean'],equal_nan=True)
                for col in ['median','q25','q75']:
                    if not np.isclose(z[col],r[col],equal_nan=True):changes.append(dict(file=f,row=int(i),column=col,old=r[col],new=z[col]))
                    t.loc[i,col]=z[col]
    for f,t in tables.items():t.to_csv(OUT/f,index=False)
    pd.DataFrame(gap).to_csv(OUT/'matching_pair_size_gaps.csv',index=False)
    (OUT/'weighted_quantile_numerical_audit.json').write_text(json.dumps(dict(tolerance=1e-12,reason='Floating point accumulation at exact 0.25/0.5/0.75 empirical CDF boundaries',mean_effects_unchanged=True,corrected_cells=changes),indent=2))
    print('Quantile boundary corrections:',len(changes))
    print(pd.DataFrame(gap).query("size_variable=='n_contacts' and pedro_version=='raw' and nearest_k==1").to_string(index=False))
if __name__=='__main__':main()
