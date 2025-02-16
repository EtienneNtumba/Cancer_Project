# Computational Genomics Pipeline for Gene Expression Analysis

## Overview
This project is a **computational genomics pipeline** designed to analyze **gene expression** in **cancer and brain development**. It automates key RNA-Seq data analysis steps, including:
- **Alignment of sequencing reads** (using STAR) to map reads to a reference genome.
- **Read quantification** (using featureCounts) to count the number of reads mapped to each gene.
- **Differential expression analysis** (using DESeq2 in R) to identify significantly differentially expressed genes between conditions.
- **Functional enrichment analysis** (using gseapy for Gene Ontology analysis) to determine the biological pathways associated with differentially expressed genes.
- **Visualization of gene expression** (using Matplotlib & Seaborn) to interpret the results effectively.

## Features
- **Fully automated pipeline** for RNA-Seq analysis, reducing manual intervention.
- **Scalability** for large datasets, allowing high-throughput processing.
- **Modular implementation**, making it easy to modify or extend different components.
- **Integration with R** for robust statistical analysis of differentially expressed genes.
- **High-quality visualization outputs** for publication-ready figures.

## Installation
### Prerequisites
Ensure you have the following dependencies installed:
- Python (>=3.8)
- R (>=4.0)
- Required Python packages:
  ```bash
  pip install pandas seaborn matplotlib bioinfokit gseapy
  ```
- Required R packages:
  ```r
  install.packages("DESeq2")
  ```

## Usage
### Running the Pipeline
The pipeline is structured as follows:
```bash
python main.py
```

### Example Workflow
1. **Prepare raw FASTQ files and metadata** - Ensure sequencing data and metadata are formatted correctly.
2. **Set up genome index for STAR alignment** - Index the reference genome for efficient mapping.
3. **Run the alignment script** - Align raw reads to the reference genome using STAR.
4. **Perform read quantification** - Count gene expression levels using featureCounts.
5. **Run differential expression analysis** - Identify upregulated and downregulated genes with DESeq2.
6. **Conduct Gene Ontology enrichment analysis** - Understand biological significance using functional enrichment tools.
7. **Generate visualization of key genes** - Create heatmaps, PCA plots, and volcano plots to interpret results.

### Expected Output
- **Aligned BAM files** (stored in `aligned_data/`), used for read mapping.
- **Gene expression counts** (stored in `counts/`), representing raw expression levels.
- **Differential expression results** (`differential_expression/deseq2_results.csv`), showing significantly changed genes.
- **GO enrichment results** (`go_enrichment_results.csv`), providing functional insights into affected pathways.
- **Heatmap visualization of selected genes** [`View Heatmap`](figures/gene_expression_heatmap.png)
- **Volcano Plot** [`View Volcano Plot`](figures/volcano_plot.png)
- **PCA Plot** [`View PCA Plot`](figures/pca_plot.png)
- **MA Plot** [`View MA Plot`](figures/ma_plot.png)
- **Boxplot** [`View Boxplot`](figures/boxplot.png)

## Example Visualizations and Analysis
### **1. Heatmap of Gene Expression**
A **heatmap** represents gene expression levels across multiple samples, helping to identify patterns of **upregulated** or **downregulated** genes.
![Gene Expression Heatmap](figures/gene_expression_heatmap.png)  
[🔗 View Full Image](figures/gene_expression_heatmap.png)

### **2. Volcano Plot**
A **volcano plot** visualizes differentially expressed genes by plotting **log2 fold change** vs **-log10 p-value**. Genes with significant changes are highlighted.
![Volcano Plot](figures/volcano_plot.png)  
[🔗 View Full Image](figures/volcano_plot.png)

### **3. PCA Plot**
A **PCA (Principal Component Analysis) plot** shows clustering of samples based on gene expression. This helps in assessing batch effects or sample separability.
![PCA Plot](figures/pca_plot.png)  
[🔗 View Full Image](figures/pca_plot.png)

### **4. MA Plot**
An **MA plot** compares the mean expression of genes across conditions with **log2 fold change**, helping to visualize systematic shifts in gene expression.
![MA Plot](figures/ma_plot.png)  
[🔗 View Full Image](figures/ma_plot.png)

### **5. Boxplot**
A **boxplot** visualizes the distribution of gene expression levels across different samples, ensuring data normalization and variance comparison.
![Boxplot](figures/boxplot.png)  
[🔗 View Full Image](figures/boxplot.png)

## Interpretation of Results
- **Clusters in PCA plot** indicate strong differentiation between conditions.
- **Significantly differentially expressed genes** in the **volcano plot** highlight potential biomarkers.
- **GO enrichment analysis** provides insights into biological pathways affected by differentially expressed genes.
- **The heatmap** helps in identifying specific genes involved in disease mechanisms.

## Contributing
If you would like to contribute to this project, feel free to submit a **pull request** or open an **issue**.

## License
This project is licensed under the MIT License.

## Contact
For inquiries or collaborations, contact:
**Etienne Ntumba Kabongo**
Email: [etienne.ntumba.kabongo@umontreal.ca](mailto:etienne.ntumba.kabongo@umontreal.ca)

