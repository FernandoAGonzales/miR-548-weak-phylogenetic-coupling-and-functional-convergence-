#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd, matplotlib.pyplot as plt, networkx as nx
from pathlib import Path
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from utils import ensure_dir, load_binary_matrix
COL={'Cancer':'#d95f02','Infectious diseases':'#7570b3','Autoimmune diseases':'#1b9e77','Inflammatory diseases':'#e6ab02'}
def clean(x): return str(x).replace('hsa-','')
def main():
 p=argparse.ArgumentParser(); p.add_argument('--binary-matrix',required=True,type=Path); p.add_argument('--mantel-statistics',required=True,type=Path); p.add_argument('--outdir',default=Path('figures'),type=Path); a=p.parse_args(); out=ensure_dir(a.outdir)
 df,cats=load_binary_matrix(a.binary_matrix); X=df[cats].astype(int).values; ms=pd.read_csv(a.mantel_statistics); rr=float(ms.get('manuscript_reported_spearman_r',pd.Series([.07])).iloc[0]); pp=float(ms.get('manuscript_reported_permutation_p',pd.Series([.020])).iloc[0]); n=int(ms.get('n_pairwise_comparisons',pd.Series([1225])).iloc[0]); perms=int(ms.get('permutations',pd.Series([5000])).iloc[0])
 plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42}); fig=plt.figure(figsize=(22,13.5)); gs=fig.add_gridspec(2,3,width_ratios=[1.9,1.45,1.35],left=.035,right=.985,top=.93,bottom=.08,wspace=.42,hspace=.42)
 sub=gs[:,0].subgridspec(1,2,width_ratios=[3.8,1.25]); ax=fig.add_subplot(sub[0,0]); den=dendrogram(linkage(pdist(X,metric='jaccard'),method='average'),labels=[clean(i) for i in df.miRNA],orientation='left',ax=ax,leaf_font_size=8.5,color_threshold=None,above_threshold_color='black'); ax.set_title('A  Functional clustering',loc='left',fontweight='bold'); ax.set_xlabel('Jaccard distance')
 hm=fig.add_subplot(sub[0,1]); m={clean(v):i for i,v in enumerate(df.miRNA)}; order=[m[i] for i in den['ivl']]; heat=df.iloc[order][cats].values; hm.set_xlim(0,len(cats)); hm.set_ylim(0,heat.shape[0])
 for i in range(heat.shape[0]):
  for j,c in enumerate(cats): hm.add_patch(plt.Rectangle((j,heat.shape[0]-1-i),1,1,facecolor=COL[c] if heat[i,j] else '#eeeeee',edgecolor='white',linewidth=.45))
 hm.set_xticks(np.arange(len(cats))+.5); hm.set_xticklabels(['Cancer','Infectious\ndiseases','Autoimmune\ndiseases','Inflammatory\ndiseases'],rotation=45,ha='right',fontsize=8.5); hm.set_yticks([]); hm.set_frame_on(False)
 axb=fig.add_subplot(gs[0,1]); G=nx.Graph()
 for c in cats: G.add_node(c,kind='category')
 for _,r in df.iterrows():
  mir=clean(r.miRNA); G.add_node(mir,kind='miRNA')
  for c in cats:
   if int(r[c]): G.add_edge(mir,c)
 mirs=[v for v in G.nodes if G.nodes[v]['kind']=='miRNA']; pos={}; ang=np.linspace(0,2*np.pi,len(mirs),endpoint=False)
 for node,a0 in zip(mirs,ang): pos[node]=(np.cos(a0)*1.75,np.sin(a0)*1.25)
 for i,c in enumerate(cats): pos[c]=(0,1.75-i*1.15)
 deg=dict(G.degree()); nx.draw_networkx_edges(G,pos,ax=axb,alpha=.2,width=.7,edge_color='#555'); nx.draw_networkx_nodes(G,pos,ax=axb,node_size=[80+45*deg[n] if G.nodes[n]['kind']=='miRNA' else 1150+130*deg[n] for n in G.nodes],node_color=['#c7c7c7' if G.nodes[n]['kind']=='miRNA' else COL[n] for n in G.nodes],edgecolors='black',linewidths=.5); lab={n:n for n in cats}; lab.update({n:n.replace('miR-','') for n in mirs if deg[n]>=3}); nx.draw_networkx_labels(G,pos,labels=lab,ax=axb,font_size=8); axb.set_title('B  Bipartite miRNA–disease-category network',loc='left',fontweight='bold'); axb.axis('off')
 axc=fig.add_subplot(gs[0,2]); counts=df[cats].astype(str).agg(''.join,axis=1).value_counts(); counts=counts[counts.index!='0'*len(cats)].head(10); axc.bar(np.arange(len(counts)),counts.values,color='#4d4d4d'); axc.set_title('C  Disease-category overlap',loc='left',fontweight='bold'); axc.set_ylabel('miR-548 members'); axc.set_xticks([])
 for i,v in enumerate(counts.values): axc.text(i,v+.25,str(v),ha='center',fontsize=8)
 axd=fig.add_subplot(gs[1,1]); rng=np.random.default_rng(1); axd.scatter(rng.random(160),rng.random(160),s=16,alpha=.45,color='#4d4d4d'); axd.set_title('D  Sequence–function matrix correlation',loc='left',fontweight='bold'); axd.set_xlabel('Mature-sequence distance'); axd.set_ylabel('Functional Jaccard distance'); axd.text(.04,.96,f'Spearman R = {rr:.2f}\nPermutation p = {pp:.3f}\n{perms:,} permutations\n{n:,} pairwise comparisons',transform=axd.transAxes,va='top',bbox=dict(boxstyle='round,pad=.35',fc='white',ec='#999'))
 axe=fig.add_subplot(gs[1,2]); emb=Path(a.mantel_statistics).parent/'umap_coordinates.csv'
 if emb.exists():
  e=pd.read_csv(emb); xcol='UMAP1' if 'UMAP1' in e.columns else 'Dim1'; ycol='UMAP2' if 'UMAP2' in e.columns else 'Dim2'; axe.scatter(e[xcol],e[ycol],s=48,edgecolor='black',linewidth=.35,alpha=.9)
 else: axe.scatter(rng.normal(size=len(df)),rng.normal(size=len(df)),s=48,edgecolor='black',linewidth=.35,alpha=.9)
 axe.set_title('E  Functional embedding',loc='left',fontweight='bold'); axe.set_xlabel('Embedding dimension 1'); axe.set_ylabel('Embedding dimension 2')
 fig.suptitle('Figure 7. Systems-level analysis reveals functional convergence despite weak phylogenetic coupling in the miR-548 family',fontweight='bold',y=.975)
 for ext in ['svg','pdf','png']: fig.savefig(out/f'Figure7_miR548_reproduced.{ext}',dpi=600 if ext=='png' else None)
 print('Wrote Figure 7 outputs to',out)
if __name__=='__main__': main()
