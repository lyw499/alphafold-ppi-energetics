"""Exact residue minima from <=8 A KD-tree neighbors; no global distance matrix."""
import json
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from Bio.Align import PairwiseAligner
from foldbench_interfaces import protein_residues
from parse_pooled_inner_zip_8A import get_residues
CUTOFF=8.0
BIN_EDGES=[3.,4.,5.,6.,8.]
BIN_NAMES=['fraction_0_3','fraction_3_4','fraction_4_5','fraction_5_6','fraction_6_8']
CHEM_GROUPS={'hydrophobic':set('AVILMFWY'),'polar_uncharged':set('STNQC'),'positive':set('KR'),'negative':set('DE'),'special':set('GPH')}
CHEM_CLASSES=['hydrophobic_hydrophobic','opposite_charge','same_charge','polar_polar','hydrophobic_polar','hydrophobic_charged','polar_charged','special_or_other']
AA3={'ALA':'A','VAL':'V','ILE':'I','LEU':'L','MET':'M','PHE':'F','TRP':'W','TYR':'Y','SER':'S','THR':'T','ASN':'N','GLN':'Q','CYS':'C','LYS':'K','ARG':'R','ASP':'D','GLU':'E','GLY':'G','PRO':'P','HIS':'H'}

def chemistry(resname):
 a=AA3.get(resname,'X')
 return next((k for k,v in CHEM_GROUPS.items() if a in v),'other')
def contact_class(a,b):
 if a in ('special','other') or b in ('special','other'):return 'special_or_other'
 if a==b=='hydrophobic':return 'hydrophobic_hydrophobic'
 if {a,b}=={'positive','negative'}:return 'opposite_charge'
 if a==b and a in ('positive','negative'):return 'same_charge'
 if a==b=='polar_uncharged':return 'polar_polar'
 if {a,b}=={'hydrophobic','polar_uncharged'}:return 'hydrophobic_polar'
 if 'hydrophobic' in (a,b):return 'hydrophobic_charged'
 return 'polar_charged'

def residue_id(r):return str(r.id[1])+r.id[2].strip()
def pack_chain(chain):
 res=get_residues(chain);atoms=[];rid=[];per=[]
 for i,r in enumerate(res):
  eligible=[a for a in r if a.element not in ('H','D')];per.append(eligible)
  atoms.extend(eligible);rid.extend([i]*len(eligible))
 xyz=np.asarray([a.coord for a in atoms],dtype=float).reshape((-1,3))
 if not len(atoms):raise ValueError('no_eligible_heavy_atoms')
 if not np.isfinite(xyz).all():raise ValueError('nonfinite_coordinates')
 return dict(chain=chain.id,res=res,atoms=atoms,rid=np.asarray(rid),xyz=xyz,tree=cKDTree(xyz),per=per,n_deuterium=sum(a.element=='D' for r in res for a in r))

def min_contacts(a,b):
 # Sparse atom-neighbor distances are temporary and all <=8 A.
 n=a['tree'].sparse_distance_matrix(b['tree'],CUTOFF,output_type='ndarray')
 if not len(n):return []
 ri=a['rid'][n['i']];rj=b['rid'][n['j']];keys=ri*len(b['res'])+rj
 order=np.lexsort((n['j'],n['i'],n['v'],keys));sorted_keys=keys[order]
 take=order[np.r_[True,sorted_keys[1:]!=sorted_keys[:-1]]]
 rows=[]
 for k in take:
  i,j=int(ri[k]),int(rj[k]);ra,rb=a['res'][i],b['res'][j];aa,ab=a['atoms'][n['i'][k]],b['atoms'][n['j'][k]]
  ga,gb=chemistry(ra.resname),chemistry(rb.resname)
  rows.append(dict(residue_index_a=i,residue_index_b=j,residue_a_id=residue_id(ra),residue_b_id=residue_id(rb),residue_a_type=ra.resname,residue_b_type=rb.resname,d_min=float(n['v'][k]),atom_a_name=aa.name,atom_b_name=ab.name,atom_a_element=aa.element,atom_b_element=ab.element,atom_a_altloc=aa.altloc.strip(),atom_b_altloc=ab.altloc.strip(),chemistry_a=ga,chemistry_b=gb,contact_class=contact_class(ga,gb)))
 return rows

def exact_targeted_distance(a,b,i,j):
 # Only one relevant residue pair (a few heavy atoms), never entire chains.
 x=np.asarray([v.coord for v in a['per'][i]],dtype=float);y=np.asarray([v.coord for v in b['per'][j]],dtype=float)
 return float(cdist(x,y).min()) if len(x) and len(y) else np.nan

def summarize(rows):
 x=np.asarray([r['d_min'] for r in rows]);n=len(x)
 result=dict(n_contacts=n,median_d_min=float(np.median(x)) if n else np.nan,mean_d_min=float(x.mean()) if n else np.nan,q25_d_min=float(np.quantile(x,.25)) if n else np.nan,q75_d_min=float(np.quantile(x,.75)) if n else np.nan,min_d_min=float(x.min()) if n else np.nan)
 for threshold in [4,5,6]:result[f'fraction_le_{threshold}']=float(np.mean(x<=threshold)) if n else np.nan
 counts=np.bincount(np.searchsorted(BIN_EDGES,x,side='left'),minlength=5)
 for key,count in zip(BIN_NAMES,counts):result[key]=count/n if n else np.nan
 for t in [.5,1.,1.5,2.]:result['n_lt_'+str(t).replace('.','_')]=int((x<t).sum())
 for c in CHEM_CLASSES:result['chem_'+c]=sum(r['contact_class']==c for r in rows)/n if n else np.nan
 return result

def residue_map(exp_chain,af3_chain,max_alignments=128):
 er,es=protein_residues(exp_chain);pr,ps=protein_residues(af3_chain)
 aligner=PairwiseAligner();aligner.mode='global';aligner.match_score=2;aligner.mismatch_score=-3;aligner.open_gap_score=-5;aligner.extend_gap_score=-.5;aligner.end_gap_score=0
 alns=aligner.align(es,ps)
 try:n=len(alns)
 except OverflowError:return {},dict(status='too_many_optimal_alignments',n_optimal_alignments='overflow')
 if n>max_alignments:return {},dict(status='too_many_optimal_alignments',n_optimal_alignments=n)
 stable=None;union=set()
 for aln in alns:
  pairs={(i,j) for ea,pa in zip(*aln.aligned) for i,j in zip(range(*ea),range(*pa)) if es[i]==ps[j] and es[i]!='X'}
  union|=pairs;stable=pairs if stable is None else stable & pairs
 stable=stable or set()
 assert len({j for _,j in stable})==len(stable)
 return dict(stable),dict(status='unique_optimal_alignment' if n==1 else 'only_positions_identical_across_all_optimal_alignments',n_optimal_alignments=n,n_exp_residues=len(er),n_af3_residues=len(pr),n_stable_identical_residues=len(stable),n_ambiguous_identical_positions=len(union-stable))
