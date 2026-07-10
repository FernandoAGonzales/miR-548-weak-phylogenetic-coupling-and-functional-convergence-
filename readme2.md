# Functional convergence of the miR-548 family reveals weak coupling between evolution and regulatory roles in human disease

## Overview

This repository contains the datasets, phylogenetic resources, processed matrices, statistical analyses, and figure-generation workflows associated with the manuscript:

> **Functional convergence of the miR-548 family reveals weak coupling between evolution and regulatory roles in human disease**

The study investigates the relationship between sequence diversification, phylogenetic structure, and disease-associated functional organization across the human miR-548 family through integrative evolutionary and systems-level analyses.

---

## Repository Contents

- Raw Data
- Phylogenetic Resources
- Processed Data
- Figure 6 Workflow
- Figure 7 Workflow
- Software Environment
- Reproducibility
- Data Availability
- Citation
- License

---

## Repository Structure

```text
.
├── data
│   ├── raw
│   └── processed
├── results
│   └── phylogeny
├── scripts
│   ├── scripts-for-figure-6
│   └── scripts-for-figure-7
├── CITATION
├── LICENSE
├── environment.yml
├── requirements.txt
└── README.md
```

## Raw Data

**Directory**

```text
data/raw/
```

| File | Description |
|------|-------------|
| `mature_miR548_alignment.fas` | Multiple sequence alignment of mature hsa-miR-548 family members used for phylogenetic reconstruction. |
| `master_curated_miR548_disease_associations.xlsx` | Unified manually curated disease-association dataset used as the primary source for downstream analyses. |
| `Supplementary_Table_S1_Autoimmune.xlsx` | Curated autoimmune disease associations. |
| `Supplementary_Table_S2_Infectious_Inflammatory.xlsx` | Curated infectious and inflammatory disease associations. |
| `Supplementary_Table_S3_Cancer.xlsx` | Curated cancer disease associations. |
| `Supplementary_Table_S4_Mature_Sequences.xlsx` | Curated mature hsa-miR-548 sequences. |

## Phylogenetic Resources

**Directory**

```text
results/phylogeny/
```

| File | Description |
|------|-------------|
| `phylogenetic_tree_ML.nwk` | Maximum Likelihood phylogenetic tree used in Figure 5. |
| `miR548_phylogeny.mtsx` | Original MEGA v12.1 project containing reconstruction settings and metadata. |
| `model_selection_results.csv` | Model-selection statistics used to identify the optimal substitution model. |

## Processed Data

These files were generated programmatically from the curated master association table.

**Directory**

- binary_disease_matrix.csv

- autoimmune_matrix.csv

- infectious_and_inflammatory_matrix.csv

- cancer_matrix.csv

- sequence_distance_matrix.csv

- functional_distance_matrix.csv

- cluster_assignments.csv

- network_metrics.csv

- upset_intersections.csv

- umap_coordinates.csv

- mantel_statistics.csv

- mantel_null_distribution.csv

### Disease-specific Matrices

| File | Description |
|------|-------------|
| `autoimmune_matrix.csv` | Autoimmune disease matrix. |
| `infectious_and_inflammatory_matrix.csv` | Infectious and inflammatory disease matrix. |
| `cancer_matrix.csv` | Cancer disease matrix. |

Disease categories:

- Cancer
- Infectious diseases
- Autoimmune diseases
- Inflammatory diseases

### Distance Matrices

| File | Description |
|------|-------------|
| `sequence_distance_matrix.csv` | Pairwise mature-sequence distances. |
| `functional_distance_matrix.csv` | Pairwise Jaccard functional distances. |

### Clustering and Embedding Outputs

| File | Description |
|------|-------------|
| `cluster_assignments.csv` | Cluster assignments generated from hierarchical clustering. |
| `umap_coordinates.csv` | Low-dimensional functional embedding coordinates used in Figure 7E. |

## Analysis Overview

### Network Analyses
network_metrics.csv
Node-level statistics calculated from the bipartite miRNA–disease category network, including:

- Degree
- Degree centrality
- Betweenness centrality
- Closeness centrality
- Eigenvector centrality

upset_intersections.csv
Disease-category overlap statistics used for UpSet analysis (Figure 7C).

### Sequence–Function Correlation Analyses

| File | Description |
|------|-------------|
| `mantel_statistics.csv` | Summary statistics for sequence–function coupling analyses. |
| `mantel_null_distribution.csv` | Null distribution generated through permutation testing. |

### Phylogenetic Resources and Reconstruction

The analyses performed in this study include:

- Maximum Likelihood phylogenetic reconstruction using MEGA v12.1.
- Model selection based on information-theoretic criteria.
- Disease-association matrix construction.
- Hierarchical clustering.
- Bipartite network analysis.
- UpSet overlap analysis.
- Sequence–function correlation analysis.
- Functional embedding using UMAP.

The corresponding datasets and outputs are available throughout this repository.

## Figure 6 Workflow

### Figure 6A – Entropy analysis of miR-548-3p strand alignments

- Shannon entropy (H) calculated for each nucleotide position.

### Figure 6B – Entropy analysis of miR-548-5p strand alignments

- Shannon entropy (H) calculated for each nucleotide position.
- Entropy distribution by functional region (seed and non-seed regions; α = 0.05).

## Software Environment

| Component | Version |
|-----------|---------|
| Python | 3.12 |
| MEGA | 12.1 |

### Primary Packages

- NumPy
- pandas
- SciPy
- NetworkX
- scikit-learn
- UMAP-learn
- matplotlib
- Biopython

### Reproducible Environments

- `environment.yml`
- `requirements.txt`

## Figure 7 Workflow

Figure 7 integrates complementary systems-level analyses to evaluate functional organization within the miR-548 family.

### Figure 7A – Functional Clustering

Hierarchical clustering of miR-548 members based on disease-association profiles using:

- Jaccard distance
- Average-linkage clustering

---

### Figure 7B – Bipartite Network Analysis

Bipartite miRNA–disease category network.

- Node size is proportional to degree centrality.

---

### Figure 7C – Disease Category Overlap

UpSet analysis quantifying overlap among:

- Cancer
- Infectious diseases
- Autoimmune diseases
- Inflammatory diseases

---

### Figure 7D – Sequence–Function Coupling Analysis

Sequence–function coupling analysis using a Mantel-style matrix correlation.

#### Workflow

1. Calculate pairwise mature-sequence distances.
2. Calculate Jaccard functional distances from binary disease profiles.
3. Convert matrices to upper-triangular vectors.
4. Calculate Spearman rank correlation.
5. Assess significance through permutation testing.

#### Results

| Metric | Value |
|---------|------:|
| Spearman correlation (R) | 0.07 |
| Permutation *p*-value | 0.020 |
| Pairwise comparisons | 1,225 |
| Permutations | 5,000 |

#### Interpretation

Phylogenetic relatedness explains only a limited fraction of the observed functional organization within the miR-548 family.

---

### Figure 7E – Functional Embedding

Functional embedding generated from disease-association profiles.

The low-dimensional representation demonstrates extensive mixing of phylogenetically distant miR-548 members within a shared functional space, consistent with widespread functional convergence.

## Data Availability

All raw datasets, sequence alignments, phylogenetic resources, processed matrices, analysis outputs, and figure source data are available in this repository and archived through Zenodo.

## Citation

If you use these datasets, analyses, or workflows, please cite:

> Gonzales-Zubiate FA, Ramos-Sanchez EM, Asanza-Sanmartin D, Solorzano-Toala EA, *et al.*
> *Functional convergence of the miR-548 family reveals weak coupling between evolution and regulatory roles in human disease.*

## License

This repository is distributed under the **MIT License**.
See the `LICENSE` file for details.
