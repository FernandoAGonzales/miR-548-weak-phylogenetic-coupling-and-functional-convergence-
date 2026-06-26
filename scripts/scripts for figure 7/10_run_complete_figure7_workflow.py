#!/usr/bin/env python3
import argparse, subprocess, sys
from pathlib import Path
def run(cmd): print('\nRunning:',' '.join(map(str,cmd))); subprocess.run(cmd,check=True)
def main():
 p=argparse.ArgumentParser(); p.add_argument('--binary-matrix',required=True,type=Path); p.add_argument('--sequences-csv',required=True,type=Path); p.add_argument('--alignment-fasta',type=Path); p.add_argument('--outdir',default=Path('data/processed'),type=Path); p.add_argument('--figdir',default=Path('figures'),type=Path); a=p.parse_args(); s=Path(__file__).resolve().parent
 for sc in ['03_generate_disease_matrices.py','04_functional_clustering.py','05_network_analysis.py','06_upset_analysis.py']:
  run([sys.executable,s/sc,'--binary-matrix',a.binary_matrix,'--outdir',a.outdir])
 cmd=[sys.executable,s/'07_sequence_function_correlation.py','--binary-matrix',a.binary_matrix,'--sequences-csv',a.sequences_csv,'--outdir',a.outdir]
 if a.alignment_fasta: cmd+=['--alignment-fasta',a.alignment_fasta]
 run(cmd); run([sys.executable,s/'08_umap_embedding.py','--binary-matrix',a.binary_matrix,'--outdir',a.outdir]); run([sys.executable,s/'09_generate_figure7.py','--binary-matrix',a.binary_matrix,'--mantel-statistics',a.outdir/'mantel_statistics.csv','--outdir',a.figdir]); print('\nComplete Figure 7 workflow finished.')
if __name__=='__main__': main()
