"""Size control using saved interface/contact tables only; never open structures."""
import os, sys, json, hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/local_interface_geometry'
OLD=BASE/'contact_type_decomposition'
OUT=BASE/'size_matched_packing'
sys.path.insert(0,str(BASE/'pedro/_deps'))
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
TYPES=['BB-BB','BB-SC','SC-SC']
BINS=['lt2','2to3','3to4','4to5','5to6','6to8']
REFS=['FoldBench Experimental','SAFES Experimental','Mgen Experimental','Mgen excl 8PBZ 8PC1','Mgen excl 8PBZ 8PC0 8PC1']
MAJOR=['Pedro AF3','FoldBench Experimental','FoldBench AF3','SAFES Experimental','SAFES AF3','Mgen Experimental']
COLORS=['#4c78a8','#f58518','#54a24b']
METRICS=['total_close']+[t+'_close' for t in TYPES]+[t+'_share' for t in TYPES]+[t+'_within_'+b for t in TYPES for b in BINS]
INPUTS=[]
def record(p):
    INPUTS.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
def csv(p,**kw):
    record(p);return pd.read_csv(p,**kw)
def stat(v,w=None):
    v=np.asarray(v,float);w=np.ones(len(v)) if w is None else np.asarray(w,float)
    ok=np.isfinite(v)&np.isfinite(w)&(w>0);v=v[ok];w=w[ok]
    if not len(v):return dict(n=0,weight_sum=0.,effective_n=0.,mean=np.nan,median=np.nan,q25=np.nan,q75=np.nan)
    order=np.argsort(v);v=v[order];w=w[order]
    # Inverse weighted empirical CDF, including unweighted summaries for consistency.
    q=v[np.searchsorted(np.cumsum(w)/w.sum(),np.array([.25,.5,.75])-1e-12,side='left')]
    return dict(n=len(v),weight_sum=w.sum(),effective_n=w.sum()**2/(w@w),mean=np.average(v,weights=w),median=q[1],q25=q[0],q75=q[2])
def savefig(name):
    plt.tight_layout();plt.savefig(OUT/'figures'/name,dpi=170,bbox_inches='tight');plt.close()
def ecdf(ax,v,w=None,**kw):
    v=np.asarray(v);w=np.ones(len(v)) if w is None else np.asarray(w)
    ok=np.isfinite(v)&(w>0);v=v[ok];w=w[ok];ix=np.argsort(v)
    if len(v):ax.step(v[ix],np.cumsum(w[ix])/w.sum(),where='post',**kw)
def weighted_ks(a,b,w):
    grid=np.unique(np.r_[a,b]);a=np.sort(np.asarray(a));ix=np.argsort(b);b=np.asarray(b)[ix];w=np.asarray(w)[ix]
    fa=np.searchsorted(a,grid,side='right')/len(a);cb=np.r_[0,np.cumsum(w)/w.sum()]
    return np.max(np.abs(fa-cb[np.searchsorted(b,grid,side='right')]))
def match(e,p,size,k=1):
    """Equal target weights; use all tied kth-nearest candidates, with replacement."""
    transform=np.log1p if size!='interface_fraction' else lambda x:x
    pv=transform(p[size].to_numpy(float));weights=np.zeros(len(p));rows=[];unsupported=[]
    lo,hi=p[size].min(),p[size].max()
    for _,r in e.iterrows():
        if not np.isfinite(r[size]) or not lo<=r[size]<=hi:
            unsupported.append(dict(interface_id=r.interface_id,size_value=r[size],pedro_min=lo,pedro_max=hi));continue
        dist=np.abs(pv-transform(r[size]));cut=np.partition(dist,min(k,len(dist))-1)[min(k,len(dist))-1]
        ix=np.flatnonzero(dist<=cut+1e-12);weights[ix]+=1/len(ix)
        for j in ix:
            rows.append(dict(experimental_interface_id=r.interface_id,pedro_interface_id=p.iloc[j].interface_id,experimental_size=r[size],pedro_size=p.iloc[j][size],absolute_transformed_gap=dist[j],target_weight=1/len(ix),n_tied_candidates=len(ix)))
    return weights,pd.DataFrame(rows),unsupported
def main():
    if OUT.exists() and any(p.name!='logs' for p in OUT.iterdir()):raise RuntimeError('Output directory is nonempty; refusing to overwrite')
    (OUT/'figures').mkdir(parents=True,exist_ok=True)
    p=OLD/'interface_contact_type_summary.parquet';record(p);d=pd.read_parquet(p)
    d['interface_id']=d.interface_id.astype(str)
    assert not d.duplicated(['dataset','interface_id','contact_type']).any()
    assert (d.groupby(['dataset','interface_id']).size()==3).all()
    assert np.array_equal(d[[f'count_{b}' for b in BINS]].sum(axis=1).to_numpy(),d.n_type.to_numpy())
    key=['dataset','interface_id'];w=d.drop_duplicates(key).set_index(key)[['structure_id','n_contacts','fraction_3_4','existing_severe_overlap','chain_mapping_unambiguous','legacy_mapping_stable']].copy()
    for t in TYPES:
        z=d[d.contact_type==t].set_index(key).reindex(w.index)
        w[t+'_close']=z.fraction_all_3to4;w[t+'_share']=z.type_fraction
        w[t+'_n']=z.n_type
        for b in BINS:w[t+'_within_'+b]=z['fraction_within_'+b];w[t+'_count_'+b]=z['count_'+b]
    w['total_close']=w[[t+'_close' for t in TYPES]].sum(axis=1,min_count=3)
    assert np.allclose(w.total_close,w.fraction_3_4,equal_nan=True)
    assert np.array_equal(w[[t+'_n' for t in TYPES]].sum(axis=1).to_numpy(),w.n_contacts.to_numpy())
    w=w.reset_index();w['total_interface_residues']=np.nan;w['interface_fraction']=np.nan
    metas=[]
    for ds,label in [('foldbench','FoldBench'),('safes','SAFES')]:
        s=csv(BASE/ds/'interface_summaries.csv',dtype={'interface_id':str})
        c=csv(BASE/ds/'residue_contacts.csv',dtype={'interface_id':str},usecols=['interface_id','structure_type','residue_index_a','residue_index_b'])
        lengths=c.groupby(['interface_id','structure_type']).agg(ia=('residue_index_a','nunique'),ib=('residue_index_b','nunique')).reset_index()
        s=s.merge(lengths,on=['interface_id','structure_type'],how='left',validate='one_to_one')
        s['total_interface_residues']=s.ia.fillna(0)+s.ib.fillna(0)
        s['interface_fraction']=s.total_interface_residues/(s.chain_a_length+s.chain_b_length)
        s['dataset']=label+' '+s.structure_type
        metas.append(s[key+['n_contacts','total_interface_residues','interface_fraction']])
    p=BASE/'pedro/pedro_interface_geometry_summary.parquet';record(p)
    s=pd.read_parquet(p,columns=['interface_id','n_contacts','total_interface_residues','interface_fraction']);s['dataset']='Pedro AF3';metas.append(s)
    s=csv(BASE/'mgen_same_species_control/mgen_experimental_interface_summary.csv',dtype={'interface_id':str});s['dataset']='Mgen Experimental'
    metas.append(s[key+['n_contacts','total_interface_residues','interface_fraction']])
    for name,ex in [('Mgen excl 8PBZ 8PC1',['8PBZ','8PC1']),('Mgen excl 8PBZ 8PC0 8PC1',['8PBZ','8PC0','8PC1'])]:
        z=s[~s.pdb_id.isin(ex)].copy();z['dataset']=name;metas.append(z[key+['n_contacts','total_interface_residues','interface_fraction']])
    meta=pd.concat(metas,ignore_index=True);meta.interface_id=meta.interface_id.astype(str)
    w=w.drop(columns=['total_interface_residues','interface_fraction']).merge(meta,on=key,validate='one_to_one',suffixes=('','_check'))
    assert (w.n_contacts==w.n_contacts_check).all();w=w.drop(columns='n_contacts_check')
    audits=[]
    for ds,g in w.groupby('dataset'):
        row=dict(dataset=ds,n_interfaces=len(g),n_zero=int(g.n_contacts.eq(0).sum()),n_nonzero=int(g.n_contacts.gt(0).sum()),zero_percent=100*g.n_contacts.eq(0).mean())
        for col in ['n_contacts','total_interface_residues','interface_fraction']+METRICS:
            row['missing_nonzero_'+col]=int(g.loc[g.n_contacts>0,col].isna().sum())
        audits.append(row)
    pd.DataFrame(audits).to_csv(OUT/'data_audit.csv',index=False)
    assert not w[['n_contacts','total_interface_residues','interface_fraction']].isna().any().any()
    assert w.loc[w.n_contacts==0,'total_close'].isna().all()
    w.to_parquet(OUT/'analysis_interface_table.parquet',index=False)
    pd.DataFrame(INPUTS).to_csv(OUT/'input_manifest.csv',index=False)
    nz=w[w.n_contacts>0].copy();pedro=nz[nz.dataset=='Pedro AF3'].copy().reset_index(drop=True)
    edges=[0]+pedro.n_contacts.quantile([.25,.5,.75]).tolist()+[np.inf]
    labels=[f'({edges[i]:g}, {edges[i+1]:g}]' for i in range(4)]
    nz['size_stratum']=pd.cut(nz.n_contacts,edges,labels=labels)
    (OUT/'methods.json').write_text(json.dumps(dict(primary_size='n_contacts',strata_edges=[None if np.isinf(x) else x for x in edges],strata_labels=labels,matching='all tied nearest log1p count; each supported experimental interface contributes unit mass; replacement permitted',bins=BINS,low_support_label='n<5 is a display caution, not an exclusion',quantiles='inverse weighted empirical CDF, same convention for unweighted summaries',no_structural_reads=True),indent=2))
    strat=[]
    for ds in w.dataset.unique():
        for lab in labels:
            g=nz[(nz.dataset==ds)&(nz.size_stratum==lab)]
            for m in METRICS:
                strat.append(dict(dataset=ds,stratum=lab,n_interfaces=len(g),median_contact_count=g.n_contacts.median(),low_support=len(g)<5,metric=m,**stat(g[m])))
    pd.DataFrame(strat).to_csv(OUT/'size_stratified_contact_type_summary.csv',index=False)
    corr=[]
    for version,g in [('raw',pedro),('existing_QC',pedro[~pedro.existing_severe_overlap.fillna(False).astype(bool)])]:
        for size in ['n_contacts','total_interface_residues','interface_fraction']:
            x=np.log1p(g[size]) if size!='interface_fraction' else g[size]
            for m in ['total_close']+[t+'_close' for t in TYPES]+[t+'_within_3to4' for t in TYPES]:
                ok=x.notna()&g[m].notna();corr.append(dict(version=version,size_variable=size,metric=m,n=int(ok.sum()),spearman=spearmanr(x[ok],g.loc[ok,m]).statistic))
    pd.DataFrame(corr).to_csv(OUT/'pedro_size_contact_type_correlations.csv',index=False)
    balance=[];outcomes=[];distance=[];ids=[];outside=[];comparison_weights={};differences=[]
    specs=[(r,'n_contacts','raw',1) for r in REFS+['FoldBench AF3','SAFES AF3']]
    specs += [(r,size,'raw',1) for r in REFS for size in ['total_interface_residues','interface_fraction']]
    specs += [(r,'n_contacts','existing_QC',1) for r in REFS]
    specs += [(r,'n_contacts','raw',20) for r in REFS]
    for ref,size,version,k in specs:
        tag=dict(reference=ref,size_variable=size,pedro_version=version,nearest_k=k)
        p=pedro if version=='raw' else pedro[~pedro.existing_severe_overlap.fillna(False).astype(bool)].reset_index(drop=True)
        e=nz[nz.dataset==ref].copy();weights,links,unsupported=match(e,p,size,k)
        outside.extend([{**tag,**u} for u in unsupported]);supported=e[~e.interface_id.isin([u['interface_id'] for u in unsupported])]
        assert len(supported)>0
        assert np.isclose(weights.sum(),len(supported))
        assert np.allclose(links.groupby('experimental_interface_id').target_weight.sum(),1)
        links=links.assign(**tag);ids.append(links)
        if size=='n_contacts' and version=='raw' and k==1:comparison_weights[ref]=(supported,p,weights)
        transform=np.log1p if size!='interface_fraction' else lambda a:np.asarray(a)
        ev=transform(supported[size]);pv=transform(p[size]);den=np.sqrt((np.var(ev)+np.var(pv))/2)
        for phase,a,ww in [('Experimental_all',e,np.ones(len(e))),('Experimental_supported',supported,np.ones(len(supported))),('Pedro_before',p,np.ones(len(p))),('Pedro_after',p,weights)]:
            z=stat(a[size],ww)
            balance.append(dict(**tag,phase=phase,**z,standardized_mean_difference=(np.average(transform(a[size]),weights=ww)-ev.mean())/den if den else np.nan,weighted_ecdf_max_gap=weighted_ks(supported[size],a[size],ww),n_experimental_outside_support=len(unsupported)))
        for phase,a,ww in [('Experimental',supported,np.ones(len(supported))),('Pedro_before',p,np.ones(len(p))),('Pedro_matched',p,weights)]:
            for m in METRICS:outcomes.append(dict(**tag,phase=phase,metric=m,**stat(a[m],ww)))
            for t in TYPES:
                for b in BINS:
                    distance.append(dict(**tag,phase=phase,contact_type=t,distance_bin=b,**stat(a[t+'_within_'+b],ww),pooled_fraction=np.sum(a[t+'_count_'+b]*ww)/np.sum(a[t+'_n']*ww)))
        for m in METRICS:
            before=stat(p[m])['mean'];after=stat(p[m],weights)['mean'];exp=stat(supported[m])['mean']
            differences.append(dict(**tag,metric=m,experimental_mean=exp,pedro_before_mean=before,pedro_matched_mean=after,before_difference=before-exp,matched_difference=after-exp,change_after_matching=after-before))
        print('Matched',tag,'targets',len(supported),'Pedro unique',int((weights>0).sum()),flush=True)
    pd.DataFrame(balance).to_csv(OUT/'matching_balance_summary.csv',index=False)
    pd.DataFrame(outcomes).to_csv(OUT/'size_matched_close_contact_summary.csv',index=False)
    pd.DataFrame(distance).to_csv(OUT/'size_matched_distance_bin_summary.csv',index=False)
    pd.DataFrame(differences).to_csv(OUT/'size_matched_effect_differences.csv',index=False)
    pd.concat(ids).to_csv(OUT/'matched_interface_ids.csv',index=False)
    pd.DataFrame(outside,columns=['reference','size_variable','pedro_version','nearest_k','interface_id','size_value','pedro_min','pedro_max']).to_csv(OUT/'outside_size_support.csv',index=False)
    paired=[];pairids=[]
    for ds in ['FoldBench','SAFES']:
        e=w[w.dataset==ds+' Experimental'];a=w[w.dataset==ds+' AF3']
        z=e.merge(a,on='interface_id',suffixes=('_exp','_af3'),validate='one_to_one')
        for scope,zz in [('full',z),('unambiguous',z[z.chain_mapping_unambiguous_exp.eq(True)&z.chain_mapping_unambiguous_af3.eq(True)])]:
            for m in METRICS:
                ok=(zz.n_contacts_exp>0)&(zz.n_contacts_af3>0)&zz[m+'_exp'].notna()&zz[m+'_af3'].notna();g=zz[ok]
                dx=np.log1p(g.n_contacts_af3)-np.log1p(g.n_contacts_exp);dy=g[m+'_af3']-g[m+'_exp']
                fit=np.linalg.lstsq(np.c_[np.ones(len(g)),dx],dy,rcond=None)[0]
                paired.append(dict(dataset=ds,scope=scope,metric=m,**stat(dy),size_adjusted_intercept=fit[0],slope_delta_logcount=fit[1],spearman=spearmanr(g[m+'_exp'],g[m+'_af3']).statistic,median_delta_logcount=dx.median()))
        pairids.append(z[['interface_id','n_contacts_exp','n_contacts_af3','structure_id_exp','structure_id_af3']].assign(dataset=ds))
    pd.DataFrame(paired).to_csv(OUT/'benchmark_paired_size_context.csv',index=False)
    pd.concat(pairids).to_csv(OUT/'benchmark_pair_ids.csv',index=False)
    figures(nz,pedro,labels,pd.DataFrame(strat),comparison_weights,pd.DataFrame(outcomes),pd.DataFrame(distance))
    # Check prior inputs stayed byte-identical, in addition to refusing old-output writes.
    for r in INPUTS:assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
    (OUT/'validation.json').write_text(json.dumps(dict(status='passed',n_interface_rows=len(w),input_hashes_unchanged=True,counts_and_fractions_reconciled=True,weights_per_target_sum_to_one=True,zero_contact_fractions_undefined=True,structure_files_read=0),indent=2))
def figures(nz,pedro,labels,strat,weights,outcomes,distance):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(9,5))
    for ds in MAJOR:
        g=nz[nz.dataset==ds];ecdf(ax,g.n_contacts,label=f'{ds} (n={len(g)})')
    ax.set(xscale='log',xlabel='Residue-pair contact count (8 Å)',ylabel='Interface ECDF',title='Nonzero interfaces; raw data and existing QC flags retained');ax.legend(fontsize=9);savefig('01_contact_count_distributions_before_matching.png')
    for m,fn in [('total_close','02_size_stratified_total_close_fraction.png'),('BB-SC_close','03_size_stratified_bbsc_close_fraction.png'),('SC-SC_close','04_size_stratified_scsc_close_fraction.png')]:
        fig,axs=plt.subplots(1,3,figsize=(16,5),sharey=True)
        for ax,ref in zip(axs,REFS[:3]):
            series=['Pedro AF3',ref]+([ref.replace('Experimental','AF3')] if not ref.startswith('Mgen') else [])
            for j,ds in enumerate(series):
                g=strat[(strat.dataset==ds)&(strat.metric==m)].set_index('stratum').reindex(labels);x=np.arange(4)+(j-1)*.16
                ax.errorbar(x,100*g['median'],yerr=np.array([100*(g['median']-g.q25),100*(g.q75-g['median'])]),fmt='o-',capsize=3,label=ds,alpha=.85)
                for xx,(_,r) in zip(x,g.iterrows()):
                    if r.n:ax.annotate(f'n={int(r.n)}'+('*' if r.n<5 else ''),(xx,100*r['median']),xytext=(0,10+10*j),textcoords='offset points',fontsize=7,ha='center')
            ax.set(xticks=range(4),xticklabels=['1–16','17–69','70–165','≥166'],xlabel='Common contact-count strata',title=ref);ax.legend(fontsize=8)
        axs[0].set_ylabel(m+' (% of all contacts)');fig.suptitle('Interface median and IQR; * n<5: low support (not excluded)');savefig(fn)
    fig,axs=plt.subplots(1,3,figsize=(15,4))
    for ax,ref in zip(axs,REFS[:3]):
        e,p,ww=weights[ref];ecdf(ax,p.n_contacts,color='gray',label='Pedro before');ecdf(ax,e.n_contacts,label='Experimental');ecdf(ax,p.n_contacts,ww,ls='--',label='Pedro matched')
        ax.set(xscale='log',xlabel='Contact count',ylabel='Interface ECDF',title=ref);ax.legend(fontsize=8)
    savefig('05_matching_balance.png')
    raw=outcomes[(outcomes.size_variable=='n_contacts')&(outcomes.pedro_version=='raw')&(outcomes.nearest_k==1)]
    fig,axs=plt.subplots(1,5,figsize=(19,5),sharey=True)
    for ax,ref in zip(axs,REFS):
        bottom=np.zeros(3)
        for t,c in zip(TYPES,COLORS):
            y=[100*raw[(raw.reference==ref)&(raw.phase==ph)&(raw.metric==t+'_close')]['mean'].iloc[0] for ph in ['Experimental','Pedro_before','Pedro_matched']]
            ax.bar(range(3),y,bottom=bottom,label=t,color=c);bottom+=y
        ax.set(xticks=range(3),xticklabels=['Exp.','Pedro\nraw','Pedro\nmatched'],title=ref.replace(' Experimental','').replace(' excl ','\nexcl '));ax.legend(fontsize=8)
    axs[0].set_ylabel('Mean interface (3,4] Å contribution (%)');fig.suptitle('Equal experimental-interface weights; matched Pedro weights; means are additive');savefig('06_size_matched_3to4_decomposition.png')
    for t,fn in [('BB-SC','07_size_matched_bbsc_distance_bins.png'),('SC-SC','08_size_matched_scsc_distance_bins.png')]:
        fig,axs=plt.subplots(1,3,figsize=(16,4),sharey=True)
        z=distance[(distance.size_variable=='n_contacts')&(distance.pedro_version=='raw')&(distance.nearest_k==1)&(distance.contact_type==t)]
        for ax,ref in zip(axs,REFS[:3]):
            for ph,off in [('Experimental',-.17),('Pedro_matched',.17)]:
                g=z[(z.reference==ref)&(z.phase==ph)].set_index('distance_bin').reindex(BINS);ax.bar(np.arange(6)+off,100*g['mean'],width=.34,label=ph)
            ax.set(xticks=range(6),xticklabels=['<2','[2,3]','(3,4]','(4,5]','(5,6]','(6,8]'],xlabel='Minimum heavy-atom distance (Å)',title=ref);ax.legend(fontsize=8)
        axs[0].set_ylabel(f'Mean within-{t} bin fraction (%)');fig.suptitle('Full 0–8 Å support, existing six bins; equal-interface summaries (not a density)');savefig(fn)
    fig,ax=plt.subplots(figsize=(9,5));p=pedro.copy();p['bin']=pd.qcut(np.log1p(p.n_contacts),20,duplicates='drop')
    for m in ['total_close']+[t+'_close' for t in TYPES]:
        g=p.groupby('bin',observed=True).agg(x=('n_contacts','median'),y=(m,'mean'));ax.plot(g.x,100*g.y,'o-',label=m)
    ax.set(xscale='log',xlabel='Contact count (median in display bin)',ylabel='Mean interface close-contact contribution (%)',title='Pedro: descriptive size relationship; correlations use every interface');ax.legend();savefig('09_pedro_size_relationship_by_contact_type.png')
if __name__=='__main__':main()
