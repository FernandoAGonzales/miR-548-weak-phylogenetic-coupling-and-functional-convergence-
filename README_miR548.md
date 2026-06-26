#### **Functional convergence of the miR-548 family reveals weak coupling between evolution and regulatory roles in human disease**

#### 

###### **Overview**



This repository contains the datasets, phylogenetic resources, processed matrices, statistical analyses, and figure-generation workflows associated with the manuscript:



**Functional convergence of the miR-548 family reveals weak coupling between evolution and regulatory roles in human disease**



The study investigates the relationship between sequence diversification, phylogenetic structure, and disease-associated functional organization across the human miR-548 family through integrative evolutionary and systems-level analyses.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

###### 

###### **Repository Contents**

###### 

###### **Raw Data**

###### 

###### **Directory structure**



data/raw/

•	mature\_miR548\_alignment.fas

•	master\_curated\_miR548\_disease\_associations.xlsx

•	Supplementary\_Table\_S1\_Autoimmune.xlsx

•	Supplementary\_Table\_S2\_Infectious\_Inflammatory.xlsx

•	Supplementary\_Table\_S3\_Cancer.xlsx

•	Supplementary\_Table\_S4\_Mature\_Sequences.xlsx



**Description**



**mature\_miR548\_alignment.fas**



Multiple sequence alignment of mature hsa-miR-548 family members used for phylogenetic reconstruction.



**master\_curated\_miR548\_disease\_associations.xlsx**



Unified manually curated literature-derived disease association dataset used as the primary source for all downstream analyses.



**Supplementary Tables S1–S3**



Disease-specific curated datasets.



**Supplementary Table S4**



Curated mature hsa-miR-548 sequences.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_





###### **Phylogenetic Resources**



###### **Directory structure**



results/phylogeny/

•	phylogenetic\_tree\_ML.nwk

•	miR548\_phylogeny.mtsx

•	model\_selection\_results.csv



**Description**



**phylogenetic\_tree\_ML.nwk**



Maximum Likelihood phylogenetic tree used in Figure 5.



**miR548\_phylogeny.mtsx**



Original MEGA v12.1 project/session file containing phylogenetic reconstruction settings and metadata.



**model\_selection\_results.csv**



Model-selection statistics used to identify the optimal substitution model.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_





###### **Processed Data**

###### 

###### **Directory structure**



data/processed/

•	binary\_disease\_matrix.csv

•	autoimmune\_matrix.csv

•	infectious\_and\_inflammatory\_matrix.csv

•	cancer\_matrix.csv

•	sequence\_distance\_matrix.csv

•	functional\_distance\_matrix.csv

•	cluster\_assignments.csv

•	network\_metrics.csv

•	upset\_intersections.csv

•	umap\_coordinates.csv

•	mantel\_statistics.csv

•	mantel\_null\_distribution.csv



**Description**

These files were generated programmatically from the curated master association table.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Disease Matrices**



**binary\_disease\_matrix.csv**



Binary disease-category matrix used throughout the systems-level analyses.



**Categories**

•	Cancer

•	Infectious diseases

•	Autoimmune diseases

•	Inflammatory diseases



**Disease-specific matrices**

•	autoimmune\_matrix.csv

•	infectious\_and\_inflammatory\_matrix.csv

•	cancer\_matrix.csv



These matrices were derived from the unified curated disease-association dataset.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Distance Matrices**



**sequence\_distance\_matrix.csv**



Pairwise mature-sequence distances among miR-548 family members.



**functional\_distance\_matrix.csv**



Pairwise functional distances calculated from disease-association profiles using Jaccard distance.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Clustering and Embedding Outputs**



**cluster\_assignments.csv**



Cluster assignments generated from hierarchical clustering.



**umap\_coordinates.csv**



Low-dimensional functional embedding coordinates used in Figure 7E.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Network Analyses**



**network\_metrics.csv**



Node-level statistics calculated from the bipartite miRNA–disease category network, including:

•	Degree

•	Degree centrality

•	Betweenness centrality

•	Closeness centrality

•	Eigenvector centrality



**upset\_intersections.csv**



Disease-category overlap statistics used for UpSet analysis (Figure 7C).



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Sequence–Function Correlation Analyses**



**mantel\_statistics.csv**



Summary statistics for sequence–function coupling analyses.



**mantel\_null\_distribution.csv**



Null distribution generated through permutation testing.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Phylogenetic Resources and Reconstruction**



Phylogenetic relationships among mature hsa-miR-548 family members were inferred using **Maximum Likelihood (ML)** analysis in **MEGA v12.1**.



Model testing was performed using information-theoretic criteria.



The **K2+G** substitution model was identified as the optimal evolutionary model and used for phylogenetic inference.



**Files provided**

•	mature\_miR548\_alignment.fas

•	model\_selection\_results.csv

•	phylogenetic\_tree\_ML.nwk

•	miR548\_phylogeny.mtsx



The deposited Newick file corresponds to the phylogenetic topology used in Figure 5 and all associated evolutionary analyses reported in the manuscript.



The MEGA session file preserves the original phylogenetic reconstruction and allows regeneration of the phylogenetic tree directly within MEGA.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### Functional and Systems-Level Analyses



Disease associations were manually curated from published studies and consolidated into a unified master dataset.



All downstream analyses were generated programmatically from this curated dataset, including:

•	Disease-category matrices

•	Hierarchical clustering

•	Network reconstruction

•	Overlap analyses

•	Sequence–function correlation analyses

•	Dimensionality-reduction embeddings



This workflow ensures consistency across all supplementary datasets and manuscript figures.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



##### Figure 6. Entropy-based variability analysis of microRNA alignments

Figure 6 analyzes sequence variability across aligned miR-548 family members.

**Figure 6A – Entropy analysis of miR-548-3p strand alignments**

**Figure 6B – Entropy analysis of miR-548-5p strand alignments**

• Shannon entropy (H) calculated for each nucleotide position.

• Entropy distribution by functional region (seed and non-seed regions), α = 0.05.





\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Figure 7 Workflow**



Figure 7 integrates multiple complementary systems-level analyses.



**Figure 7A – Functional Clustering**



Hierarchical clustering of miR-548 members based on disease-association profiles using:

•	Jaccard distance

•	Average-linkage clustering



**Figure 7B – Bipartite Network Analysis**



Bipartite miRNA–disease category network.



Node size is proportional to degree centrality.



**Figure 7C – Disease Category Overlap**



UpSet analysis quantifying overlap among:

•	Cancer

•	Infectious diseases

•	Autoimmune diseases

•	Inflammatory diseases



**Figure 7D – Sequence–Function Coupling Analysis**



Sequence–function coupling analysis using a Mantel-style matrix correlation.



**Workflow**

1\.	Calculate pairwise mature-sequence distances.

2\.	Calculate Jaccard functional distances from binary disease profiles.

3\.	Convert matrices to upper-triangular vectors.

4\.	Calculate Spearman rank correlation.

5\.	Assess significance through permutation testing.



**Final reported result**

•	Spearman R = 0.07

•	Permutation p = 0.020

•	Pairwise comparisons = 1,225

•	Permutations = 5,000



**Interpretation**



Phylogenetic relatedness explains only a limited fraction of the observed functional organization within the miR-548 family.



**Figure 7E – Functional Embedding**



Functional embedding generated from disease-association profiles.



Low-dimensional representation demonstrates extensive mixing of phylogenetically distant miR-548 members within shared functional space, consistent with widespread functional convergence.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Software Environment**



All analyses were implemented using:



**Python 3.12**



**Primary packages**

•	NumPy

•	pandas

•	SciPy

•	NetworkX

•	scikit-learn

•	UMAP-learn

•	matplotlib

•	Biopython



**Reproducible environments**

•	environment.yml

•	requirements.txt



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Reproducibility**



This repository contains all datasets and analysis outputs required to reproduce the principal findings reported in the manuscript.



**Included resources**

•	Mature miR-548 sequence datasets

•	Sequence alignments

•	Phylogenetic reconstruction files

•	Model-selection outputs

•	Curated disease-association datasets

•	Disease-category matrices

•	Clustering outputs

•	Network analyses

•	Overlap analyses

•	Sequence–function correlation analyses

•	Dimensionality-reduction embeddings

•	Source data underlying manuscript figures



Together, these resources enable independent verification of the evolutionary, functional, and systems-level analyses presented in the study.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Data Availability**



All sequence datasets, alignments, phylogenetic resources, processed matrices, and analysis outputs are available through this repository and archived through Zenodo.



The repository includes:

•	phylogenetic\_tree\_ML.nwk

•	miR548\_phylogeny.mtsx



which correspond to the original phylogenetic reconstruction used in the manuscript.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **Citation**



If you use these datasets, analyses, or workflows, please cite:

Gonzales-Zubiate FA, Ramos-Sanchez EM, Asanza-Sanmartin D, Solorzano-Toala EA, et al.



Functional convergence of the miR-548 family reveals weak coupling between evolution and regulatory roles in human disease.



\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_



###### **License**



This repository is distributed under the **MIT License**.



See **LICENSE** for details.

