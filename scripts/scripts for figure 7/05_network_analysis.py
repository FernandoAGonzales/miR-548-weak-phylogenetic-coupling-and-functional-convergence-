#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd, networkx as nx
from pathlib import Path
from utils import ensure_dir, load_binary_matrix
def main():
 p=argparse.ArgumentParser(); p.add_argument('--binary-matrix',required=True,type=Path); p.add_argument('--outdir',default=Path('data/processed'),type=Path); a=p.parse_args(); out=ensure_dir(a.outdir)
 df,cats=load_binary_matrix(a.binary_matrix); G=nx.Graph()
 for c in cats: G.add_node(c,node_type='disease_category')
 for _,r in df.iterrows():
  G.add_node(r.miRNA,node_type='miRNA')
  for c in cats:
   if int(r[c])==1: G.add_edge(r.miRNA,c)
 deg=dict(G.degree()); dc=nx.degree_centrality(G); bc=nx.betweenness_centrality(G); cc=nx.closeness_centrality(G)
 try: ec=nx.eigenvector_centrality_numpy(G)
 except Exception: ec={n:np.nan for n in G.nodes}
 rows=[]
 for n in G.nodes:
  rec={'node':n,'node_type':G.nodes[n]['node_type'],'degree':deg[n],'degree_centrality':dc[n],'betweenness_centrality':bc[n],'closeness_centrality':cc[n],'eigenvector_centrality':ec[n]}
  if G.nodes[n]['node_type']=='miRNA':
   r=df.loc[df.miRNA==n].iloc[0]
   for c in cats: rec[c]=int(r[c])
   rec['number_of_categories']=int(r[cats].sum())
  rows.append(rec)
 pd.DataFrame(rows).sort_values(['node_type','degree'],ascending=[True,False]).to_csv(out/'network_metrics.csv',index=False); print('Wrote network_metrics.csv')
if __name__=='__main__': main()
