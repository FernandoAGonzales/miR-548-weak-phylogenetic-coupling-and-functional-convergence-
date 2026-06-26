#!/usr/bin/env python3
import argparse
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import pdist
from utils import ensure_dir, load_binary_matrix, compute_functional_distance_matrix
def main():
 p=argparse.ArgumentParser(); p.add_argument('--binary-matrix',required=True,type=Path); p.add_argument('--outdir',default=Path('data/processed'),type=Path); a=p.parse_args(); out=ensure_dir(a.outdir)
 df,cats=load_binary_matrix(a.binary_matrix); fun=compute_functional_distance_matrix(df,cats); fun.to_csv(out/'functional_distance_matrix.csv',index_label='miRNA')
 Z=linkage(pdist(df[cats].astype(int).values,metric='jaccard'),method='average'); names=df.miRNA.tolist(); order={names[i]:r+1 for r,i in enumerate(leaves_list(Z))}
 o=df[['miRNA']+cats].copy(); o['cluster_3']=fcluster(Z,3,criterion='maxclust'); o['cluster_4']=fcluster(Z,4,criterion='maxclust'); o['cluster_5']=fcluster(Z,5,criterion='maxclust'); o['dendrogram_order']=o.miRNA.map(order); o['number_of_categories']=o[cats].sum(axis=1); o.sort_values('dendrogram_order').to_csv(out/'cluster_assignments.csv',index=False)
 print('Wrote clustering outputs to',out)
if __name__=='__main__': main()
