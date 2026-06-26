#!/usr/bin/env python3
from __future__ import annotations
import re
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np, pandas as pd
from scipy.spatial.distance import pdist, squareform
DISEASE_CATEGORIES=['Cancer','Infectious diseases','Autoimmune diseases','Inflammatory diseases']
def normalize_mirna_name(name:str)->str:
    name=str(name).strip().replace('(3p)','-3p').replace('(5p)','-5p').replace('(','-').replace(')','').replace('--','-')
    return name if name.startswith('hsa-') else 'hsa-'+name
def base_mirna_name(name:str)->str: return re.sub(r'-(3p|5p)$','',normalize_mirna_name(name))
def standardize_disease_columns(df:pd.DataFrame):
    ren={}
    for c in df.columns:
        lc=c.lower()
        if 'cancer' in lc: ren[c]='Cancer'
        elif 'infect' in lc: ren[c]='Infectious diseases'
        elif 'auto' in lc: ren[c]='Autoimmune diseases'
        elif 'inflam' in lc: ren[c]='Inflammatory diseases'
    df=df.rename(columns=ren); cats=[c for c in DISEASE_CATEGORIES if c in df.columns]
    for c in cats: df[c]=pd.to_numeric(df[c], errors='coerce').fillna(0).astype(int)
    return df,cats
def load_binary_matrix(path:Path):
    df=pd.read_csv(path); mir=next((c for c in df.columns if c.lower().strip() in {'mirna','mirna_name','mirna name','name'}), df.columns[0])
    df=df.rename(columns={mir:'miRNA'}); df['miRNA']=df['miRNA'].map(normalize_mirna_name); df,cats=standardize_disease_columns(df)
    return df[['miRNA']+cats].drop_duplicates('miRNA').reset_index(drop=True),cats
def load_sequences_csv(path:Path)->Dict[str,str]:
    df=pd.read_csv(path); m=next((c for c in df.columns if 'mir' in c.lower() or 'name' in c.lower()), df.columns[0]); s=next((c for c in df.columns if 'seq' in c.lower()), df.columns[1])
    out={}
    for _,r in df.iterrows():
        seq=str(r[s]).strip().upper().replace('T','U'); name=normalize_mirna_name(r[m])
        if re.fullmatch(r'[ACGU\-]+', seq): out[name]=seq
    return out
def load_sequences_fasta(path:Path)->Dict[str,str]:
    out={}; cur=None; chunks=[]
    for line in open(path, encoding='utf-8', errors='replace'):
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if cur is not None: out[normalize_mirna_name(cur)]=''.join(chunks).upper().replace('T','U')
            cur=line[1:].split()[0]; chunks=[]
        else: chunks.append(line)
    if cur is not None: out[normalize_mirna_name(cur)]=''.join(chunks).upper().replace('T','U')
    return out
def match_sequences(names:List[str], seqs:Dict[str,str]):
    matched={}
    for n in names:
        if n in seqs: matched[n]=seqs[n]; continue
        cand=[k for k in seqs if base_mirna_name(k)==base_mirna_name(n)]
        if cand:
            arm='-3p' if n.endswith('-3p') else '-5p' if n.endswith('-5p') else None
            same=[k for k in cand if arm and k.endswith(arm)]; matched[n]=seqs[(same or cand)[0]]
    return matched
def padded_mismatch_distance(a,b):
    a=str(a).upper().replace('T','U'); b=str(b).upper().replace('T','U'); L=max(len(a),len(b))
    if not L: return np.nan
    a=a.ljust(L,'-'); b=b.ljust(L,'-'); return sum(x!=y for x,y in zip(a,b))/L
def compute_sequence_distance_matrix(names, matched):
    mat=np.zeros((len(names),len(names)))
    for i in range(len(names)):
        for j in range(i+1,len(names)): mat[i,j]=mat[j,i]=padded_mismatch_distance(matched[names[i]], matched[names[j]])
    return pd.DataFrame(mat,index=names,columns=names)
def compute_functional_distance_matrix(df,cats):
    return pd.DataFrame(squareform(pdist(df[cats].astype(int).values, metric='jaccard')), index=df.miRNA, columns=df.miRNA)
def upper_triangle_values(df): return df.values[np.triu_indices(df.shape[0],1)]
def ensure_dir(p): p=Path(p); p.mkdir(parents=True, exist_ok=True); return p
