#!/usr/bin/env python3
import argparse
from pathlib import Path
from utils import ensure_dir, load_binary_matrix
def main():
 p=argparse.ArgumentParser(); p.add_argument('--binary-matrix',required=True,type=Path); p.add_argument('--outdir',default=Path('data/processed'),type=Path); a=p.parse_args(); out=ensure_dir(a.outdir)
 df,cats=load_binary_matrix(a.binary_matrix); df.to_csv(out/'binary_disease_matrix.csv',index=False)
 for c in cats: df.loc[df[c]==1,['miRNA',c]].to_csv(out/(c.lower().replace(' ','_')+'_matrix.csv'),index=False)
 print('Wrote disease matrices to',out)
if __name__=='__main__': main()
