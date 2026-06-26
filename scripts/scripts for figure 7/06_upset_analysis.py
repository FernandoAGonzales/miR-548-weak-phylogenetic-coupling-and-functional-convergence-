#!/usr/bin/env python3
import argparse, pandas as pd
from pathlib import Path
from utils import ensure_dir, load_binary_matrix
def main():
 p=argparse.ArgumentParser(); p.add_argument('--binary-matrix',required=True,type=Path); p.add_argument('--outdir',default=Path('data/processed'),type=Path); a=p.parse_args(); out=ensure_dir(a.outdir)
 df,cats=load_binary_matrix(a.binary_matrix); rows=[]; patterns=df[cats].astype(str).agg(''.join,axis=1)
 for pat,g in df.groupby(patterns):
  if pat=='0'*len(cats): continue
  active=[c for c,b in zip(cats,pat) if b=='1']; rows.append({'intersection_pattern':pat,'active_categories':'; '.join(active),'n_categories':len(active),'n_miRNAs':len(g),'miRNAs':'; '.join(g.miRNA)})
 pd.DataFrame(rows).sort_values(['n_miRNAs','n_categories','active_categories'],ascending=[False,False,True]).to_csv(out/'upset_intersections.csv',index=False); print('Wrote upset_intersections.csv')
if __name__=='__main__': main()
