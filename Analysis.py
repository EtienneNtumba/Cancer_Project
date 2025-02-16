import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from bioinfokit.analys import get_data, stat
from gseapy import enrichr

# Define paths
raw_data_dir = "./raw_data"
aligned_data_dir = "./aligned_data"
counts_dir = "./counts"
diff_expr_dir = "./differential_expression"

def run_alignment(sample_fastq, index_path, output_bam):
    """ Run RNA-Seq alignment using STAR """
    cmd = f"STAR --runThreadN 8 --genomeDir {index_path} --readFilesIn {sample_fastq} \
            --outFileNamePrefix {output_bam} --outSAMtype BAM SortedByCoordinate"
    subprocess.run(cmd, shell=True)

def count_reads(aligned_bam, gtf_file, output_counts):
    """ Count reads using featureCounts """
    cmd = f"featureCounts -T 8 -a {gtf_file} -o {output_counts} {aligned_bam}"
    subprocess.run(cmd, shell=True)

def differential_expression(counts_file, metadata_file, output_dir):
    """ Perform differential expression analysis with DESeq2 in R """
    r_script = f"""
    library(DESeq2)
    countData <- read.csv('{counts_file}', row.names=1)
    colData <- read.csv('{metadata_file}', row.names=1)
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
    dds <- DESeq(dds)
    res <- results(dds)
    write.csv(as.data.frame(res), file='{output_dir}/deseq2_results.csv')
    """
    subprocess.run(["Rscript", "-e", r_script], shell=True)

def enrichment_analysis(diff_expr_file):
    """ Perform GO enrichment analysis on differentially expressed genes """
    df = pd.read_csv(diff_expr_file)
    sig_genes = df[df['padj'] < 0.05]['gene'].tolist()
    enr = enrichr(gene_list=sig_genes, gene_sets='GO_Biological_Process_2021', organism='human')
    enr.results.to_csv("./go_enrichment_results.csv", index=False)

def plot_gene_expression(counts_file, genes_of_interest):
    """ Plot gene expression heatmap for selected genes """
    df = pd.read_csv(counts_file, index_col=0)
    df_selected = df.loc[genes_of_interest]
    plt.figure(figsize=(10, 8))
    sns.heatmap(df_selected, cmap="coolwarm", annot=True)
    plt.title("Gene Expression Heatmap")
    plt.show()

# Example usage
if __name__ == "__main__":
    # Example file paths
    sample_fastq = "sample_1.fastq"
    index_path = "./genome_index"
    gtf_file = "./genome_annotation.gtf"
    metadata_file = "./metadata.csv"
    output_bam = os.path.join(aligned_data_dir, "sample_1.bam")
    output_counts = os.path.join(counts_dir, "counts.txt")
    diff_expr_file = os.path.join(diff_expr_dir, "deseq2_results.csv")
    
    # Run the pipeline
    run_alignment(sample_fastq, index_path, output_bam)
    count_reads(output_bam, gtf_file, output_counts)
    differential_expression(output_counts, metadata_file, diff_expr_dir)
    enrichment_analysis(diff_expr_file)
    
    # Plot selected genes
    genes_of_interest = ["BRCA1", "TP53", "SOX2", "PAX6"]
    plot_gene_expression(output_counts, genes_of_interest)
