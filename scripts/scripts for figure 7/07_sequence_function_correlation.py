#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
from utils import *
def main():
 p=argparse.ArgumentParser(); p.add_argument('--binary-matrix',required=True,type=Path); p.add_argument('--sequences-csv',type=Path); p.add_argument('--alignment-fasta',type=Path); p.add_argument('--outdir',default=Path('data/processed'),type=Path); p.add_argument('--permutations',default=5000,type=int); p.add_argument('--seed',default=42,type=int); p.add_argument('--reported-r',default=0.07,type=float); p.add_argument('--reported-p',default=0.020,type=float); a=p.parse_args(); out=ensure_dir(a.outdir)
 df,cats=load_binary_matrix(a.binary_matrix); names=df.miRNA.tolist(); seqs={}
 if a.sequences_csv: seqs.update(load_sequences_csv(a.sequences_csv))
 if a.alignment_fasta:
  for k,v in load_sequences_fasta(a.alignment_fasta).items(): seqs.setdefault(k,v)
 matched=match_sequences(names,seqs); missing=[n for n in names if n not in matched]
 if missing: raise ValueError('Missing sequences: '+str(missing))
 seqd=compute_sequence_distance_matrix(names,matched); fund=compute_functional_distance_matrix(df,cats); seqd.to_csv(out/'sequence_distance_matrix.csv',index_label='miRNA'); fund.to_csv(out/'functional_distance_matrix.csv',index_label='miRNA')
 x=upper_triangle_values(seqd); y=upper_triangle_values(fund); r,_=spearmanr(x,y); rng=np.random.default_rng(a.seed); vals=fund.values; iu=np.triu_indices(len(names),1); null=np.empty(a.permutations)
 for i in range(a.permutations):
  perm=rng.permutation(len(names)); null[i],_=spearmanr(x,vals[perm][:,perm][iu])
 pval=(np.sum(np.abs(null)>=abs(r))+1)/(a.permutations+1)
 pd.DataFrame([{'analysis':'Mantel-style sequence-function correlation','implementation':'Spearman correlation of flattened upper-triangular distance matrices with miRNA-label permutations','sequence_distance_metric':'normalized pairwise mature-sequence mismatch/gap distance','functional_distance_metric':'Jaccard distance from binary disease-category profiles','n_miRNAs':len(names),'n_pairwise_comparisons':len(x),'permutations':a.permutations,'random_seed':a.seed,'computed_spearman_r':r,'computed_permutation_p':pval,'manuscript_reported_spearman_r':a.reported_r,'manuscript_reported_permutation_p':a.reported_p}]).to_csv(out/'mantel_statistics.csv',index=False)
 pd.DataFrame({'permutation':np.arange(1,a.permutations+1),'null_spearman_r':null}).to_csv(out/'mantel_null_distribution.csv',index=False)
 print('Computed R=',r,'p=',pval,'; manuscript reported R=',a.reported_r,'p=',a.reported_p)
if __name__=='__main__': main()
