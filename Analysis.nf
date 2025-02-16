nextflow.enable.dsl=2

process ALIGNMENT {
    tag "Aligning reads using STAR"
    input:
        path sample_fastq
        path index_path
    output:
        path "sorted_output.bam"
    script:
    """
    STAR --runThreadN 8 --genomeDir ${index_path} --readFilesIn ${sample_fastq} \
         --outFileNamePrefix sorted_output --outSAMtype BAM SortedByCoordinate
    """
}

process COUNT_READS {
    tag "Counting reads using featureCounts"
    input:
        path aligned_bam
        path gtf_file
    output:
        path "counts.txt"
    script:
    """
    featureCounts -T 8 -a ${gtf_file} -o counts.txt ${aligned_bam}
    """
}

process DIFFERENTIAL_EXPRESSION {
    tag "Running DESeq2 in R"
    input:
        path counts_file
        path metadata_file
    output:
        path "deseq2_results.csv"
    script:
    """
    Rscript -e "\
    library(DESeq2); \
    countData <- read.csv('${counts_file}', row.names=1); \
    colData <- read.csv('${metadata_file}', row.names=1); \
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition); \
    dds <- DESeq(dds); \
    res <- results(dds); \
    write.csv(as.data.frame(res), file='deseq2_results.csv')\"
    """
}

process ENRICHMENT_ANALYSIS {
    tag "Performing GO enrichment analysis"
    input:
        path diff_expr_file
    output:
        path "go_enrichment_results.csv"
    script:
    """
    python -c "\
    import pandas as pd; \
    from gseapy import enrichr; \
    df = pd.read_csv('${diff_expr_file}'); \
    sig_genes = df[df['padj'] < 0.05]['gene'].tolist(); \
    enr = enrichr(gene_list=sig_genes, gene_sets='GO_Biological_Process_2021', organism='human'); \
    enr.results.to_csv('go_enrichment_results.csv', index=False)\"
    """
}

workflow {
    sample_fastq = file("sample_1.fastq")
    index_path = file("./genome_index")
    gtf_file = file("./genome_annotation.gtf")
    metadata_file = file("./metadata.csv")

    aligned_bam = ALIGNMENT(sample_fastq, index_path)
    counts_file = COUNT_READS(aligned_bam, gtf_file)
    diff_expr_file = DIFFERENTIAL_EXPRESSION(counts_file, metadata_file)
    ENRICHMENT_ANALYSIS(diff_expr_file)
}
