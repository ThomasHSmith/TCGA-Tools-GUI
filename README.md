# TCGA-Tools-GUI
User-friendly, graphical interface for easily performing basic analysis and visualization of RNA-seq data publicly available from [The Cancer Genome Atlas (TCGA)](https://cancergenome.nih.gov/).  Currently only supports loading locally-available TCGA datasets that have been parsed and aggregated using [TCGA-Tools](https://github.com/ThomasHSmith/TCGA-Tools).

![image](screens/example.png)

## Features
- Easily view dataset metadata, including solid vs. metastatic samples and staging information when available.
- Extract gene expression data from datasets by entering Ensembl Gene ID.
	- Built-in search tool allows for easy retrieval of Ensembl Gene IDs from gene symbol or Uniprot ID.
- Easily log2-transform raw expression data (FPKM-UQ) data or represent as z-scores.
- Calculate Pearson and Spearman correlation values between genes of interest.
- Plot histograms of individual gene expression profiles.
- Generate scatter plots to visualize co-expression of two genes.
- Create heatmaps to visualize expression profiles of several genes of interest, and how they vary in Normal vs. Tumor samples.
- Export raw or transformed data to Excel
- Subset data based on tumor stage metadata to interrogate how gene expression changes with tumor progression.


## Authors

* **Thomas Smith** - [ThomasHSmith](https://github.com/ThomasHSmith)
