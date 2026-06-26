#!/usr/bin/env python3
import argparse, pandas as pd
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
from utils import ensure_dir, load_binary_matrix
try:
 import umap; HAS_UMAP=True
except Exception: HAS_UMAP=False
def main():
 p=argparse.ArgumentParser(); p.add_argument('--binary-matrix',required=True,type=Path); p.add_argument('--outdir',default=Path('data/processed'),type=Path); p.add_argument('--seed',default=42,type=int); p.add_argument('--n-neighbors',default=15,type=int); p.add_argument('--min-dist',default=.25,type=float); a=p.parse_args(); out=ensure_dir(a.outdir)
 df,cats=load_binary_matrix(a.binary_matrix); X=df[cats].astype(int).values
 if HAS_UMAP:
  coords=umap.UMAP(metric='jaccard',random_state=a.seed,n_neighbors=min(a.n_neighbors,len(df)-1),min_dist=a.min_dist).fit_transform(X); cols=('UMAP1','UMAP2'); method='UMAP'
 else:
  coords=MDS(n_components=2,dissimilarity='precomputed',random_state=a.seed,n_init=4).fit_transform(squareform(pdist(X,metric='jaccard'))); cols=('Dim1','Dim2'); method='MDS_fallback'
 o=df[['miRNA']+cats].copy(); o[cols[0]]=coords[:,0]; o[cols[1]]=coords[:,1]; o['embedding_method']=method; o['random_state']=a.seed; o.to_csv(out/'umap_coordinates.csv',index=False); print('Wrote umap_coordinates.csv using',method)
if __name__=='__main__': main()
