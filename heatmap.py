import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Generate synthetic gene expression data for cancer-related genes
np.random.seed(42)
genes = ["TP53", "BRCA1", "BRCA2", "MYC", "EGFR", "PTEN", "RB1", "KRAS", "BRAF", "CDKN2A"]
samples = [f"Sample_{i+1}" for i in range(10)]
expression_data = np.random.rand(len(genes), len(samples)) * 100

# Create DataFrame
df = pd.DataFrame(expression_data, index=genes, columns=samples)

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(df, cmap="coolwarm", annot=True, fmt=".1f", linewidths=0.5)

# Title and labels
plt.title("Gene Expression Heatmap for Cancer-related Genes")
plt.xlabel("Samples")
plt.ylabel("Genes")

# Save as PNG
heatmap_path = "gene_expression_heatmap.png"
plt.savefig(heatmap_path, dpi=300)
plt.show()

# Return the file path
heatmap_path
