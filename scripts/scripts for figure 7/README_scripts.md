# miR-548 Figure 7 reproducibility scripts

These scripts reproduce the final manuscript systems-level analyses from the finalized datasets.





\## Workflow overview



Figure 7 integrates multiple systems-level analyses generated from the finalized binary disease-category matrix and mature miR-548 sequence dataset.



The workflow performs:



1\. Disease-category matrix standardization

2\. Functional-distance matrix calculation using Jaccard distance

3\. Hierarchical clustering using average linkage

4\. Bipartite miRNA–disease-category network reconstruction

5\. Disease-category overlap analysis

6\. Sequence–function correlation analysis

7\. UMAP embedding of disease-association profiles

8\. Generation of publication-ready Figure 7 outputs





\## Scripts



\### 03\_generate\_disease\_matrices.py

Generates standardized binary disease-category matrices from the curated input matrix.



\### 04\_functional\_clustering.py

Calculates Jaccard functional distances and performs average-linkage hierarchical clustering.



\### 05\_network\_analysis.py

Reconstructs the bipartite miRNA–disease-category network and calculates node-level network metrics.



\### 06\_upset\_analysis.py

Quantifies overlap among Cancer, Infectious diseases, Autoimmune diseases, and Inflammatory diseases.



\### 07\_sequence\_function\_correlation.py

Performs Mantel-style sequence–function correlation analysis using Spearman correlation and 5,000 permutations.



\### 08\_umap\_embedding.py

Generates a two-dimensional functional embedding from the binary disease-category matrix.



\### 09\_generate\_figure7.py

Generates publication-ready Figure 7 outputs.



\### 10\_run\_complete\_figure7\_workflow.py

Runs the complete Figure 7 workflow.





## Run all analyses

```bash
python scripts/10\_run\_complete\_figure7\_workflow.py \\
  --binary-matrix data/raw/Supplementary\_Data\_S2\_Curated\_miR548\_binary\_disease\_matrix.csv \\
  --sequences-csv data/raw/Supplementary\_Table\_S4\_Mature\_sequences\_hsa\_miR548.csv \\
  --alignment-fasta data/raw/mature\_miR548\_alignment.fas \\
  --outdir data/processed \\
  --figdir figures
```

## Software requirements

Python 3.12

Required packages

- NumPy
- pandas
- SciPy
- NetworkX
- scikit-learn
- umap-learn
- matplotlib

The complete computational environment can be reproduced using the provided `environment.yml` or `requirements.txt` files.


The following parameters correspond to those used to generate the analyses and publication figures reported in the manuscript.

\## Key parameters



\### Functional clustering

\- Distance metric: Jaccard distance

\- Clustering method: average linkage



\### Sequence–function correlation

\- Sequence-distance metric: normalized pairwise mature-sequence mismatch/gap distance

\- Functional-distance metric: Jaccard distance

\- Correlation test: Spearman rank correlation

\- Permutations: 5,000

\- Random seed: 42

\- Pairwise comparisons: 1,225



\### UMAP embedding

\- Input: binary miRNA–disease category matrix

\- n\_components: 2

\- n\_neighbors: 5

\- min\_dist: 0.30

\- metric: jaccard

\- init: spectral

\- random\_state: 42





## Outputs



The complete workflow generates the following files in `data/processed/`:



\- `binary\_disease\_matrix.csv`

\- `functional\_distance\_matrix.csv`

\- `sequence\_distance\_matrix.csv`

\- `cluster\_assignments.csv`

\- `network\_metrics.csv`

\- `upset\_intersections.csv`

\- `umap\_coordinates.csv`

\- `mantel\_statistics.csv`

\- `mantel\_null\_distribution.csv`



Publication-ready figure outputs are generated in `figures/`:



\- `Figure7.svg`

\- `Figure7.pdf`

\- `Figure7.png`



\## Reproducibility note



All Figure 7 analyses are generated programmatically from the finalized binary disease-category matrix and mature miR-548 sequence dataset. The workflow uses fixed analysis parameters and a fixed random seed (random_state = 42) for stochastic procedures (UMAP and permutation testing), ensuring deterministic reproduction of the reported clustering, network, sequence–function correlation, and embedding results.

The scripts can be executed individually to reproduce specific analyses or collectively through `10_run_complete_figure7_workflow.py` to regenerate all Figure 7 results from the finalized datasets.



