
## Project Overview

This project focuses on the transcriptomic analysis of the human cerebral cortex and its comparison between two age groups: newborns and adults, using clustering and functional enrichment methods. The goal was to identify clusters of differentially expressed genes and characterize their functions, as well as explore their enrichment across different cell types and developmental stages.

## Project Steps

1. Data Preprocessing and Normalization  
   For each gene, the mean expression across all samples was subtracted to normalize the expression values. This allowed comparing the shape of gene expression profiles across two age groups.

2. Gene Clustering  
   Differentially expressed genes were hierarchically clustered based on normalized expression to estimate the optimal number of clusters. Then spectral clustering was applied using the selected number of clusters.

3. Average Expression Profiling  
   For each of the 8 gene clusters, average expression across cortical layers was calculated and visualized to explore spatial patterns.

4. Functional Enrichment Analysis  
   To understand the biological meaning of each cluster, overrepresentation analysis was performed using the GSEAPy package across multiple databases, including:
   - Gene Ontology (biological process, molecular function, cellular component)
   - KEGG
   - SynGO
   - Reactome
   - Azimuth cell types

5. Cell-Type Enrichment Analysis  
   Enrichment analysis was based on cell-type-specific marker genes using two independent human single-cell datasets:
   - Adult dataset: 17 cell types; differential gene expression was tested via logistic regression and used in GSEA.
   - Prenatal/neonatal dataset: 9 cell types; same approach was applied.

   The enrichment results were visualized as heatmaps of normalized enrichment scores for each cell type and gene cluster.

6. Visualization  
   Heatmaps and cluster expression profiles were created using Seaborn and Matplotlib. Figures are stored in the figures/ folder.

## Tools and Libraries

- R: edgeR
- Python: pandas, scikit-learn, scanpy, seaborn, matplotlib, GSEAPy
