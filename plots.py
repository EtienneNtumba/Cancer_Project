import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# Simulating gene expression data from the heatmap
genes = ["TP53", "BRCA1", "BRCA2", "MYC", "EGFR", "PTEN", "RB1", "KRAS", "BRAF", "CDKN2A"]
samples = [f"Sample_{i+1}" for i in range(10)]
expression_data = np.array([
    [37.5, 95.1, 73.2, 59.9, 15.6, 15.6, 5.8, 86.6, 60.1, 70.8],
    [2.1, 97.0, 83.2, 21.2, 18.2, 18.3, 30.4, 52.5, 43.2, 29.1],
    [61.2, 13.9, 29.2, 36.6, 45.6, 78.5, 20.0, 51.4, 59.2, 4.6],
    [60.8, 17.1, 6.5, 94.9, 96.6, 80.8, 30.5, 9.8, 68.4, 44.0],
    [12.2, 49.5, 3.4, 90.9, 25.9, 66.3, 31.2, 52.0, 54.7, 18.5],
    [97.0, 77.5, 93.9, 89.5, 59.8, 92.2, 8.8, 19.6, 4.5, 32.5],
    [38.9, 27.1, 82.9, 35.7, 28.1, 54.3, 14.1, 80.2, 7.5, 98.7],
    [77.2, 19.9, 0.6, 81.5, 70.7, 72.9, 77.1, 7.4, 35.8, 11.6],
    [86.3, 62.3, 33.1, 6.4, 31.1, 32.5, 73.0, 63.8, 88.7, 47.2],
    [12.0, 71.3, 76.1, 56.1, 77.1, 49.4, 52.3, 42.8, 2.5, 10.8]
])

df = pd.DataFrame(expression_data, index=genes, columns=samples)

# ----- 1. Volcano Plot -----
log2_fold_change = np.log2(df.mean(axis=1) / np.median(df.mean(axis=1)))  # Simulated log2 fold change
p_values = np.random.rand(len(genes))  # Simulated p-values
neg_log_p = -np.log10(p_values)

plt.figure(figsize=(8, 6))
plt.scatter(log2_fold_change, neg_log_p, color="gray", alpha=0.7)
plt.axhline(-np.log10(0.05), color='red', linestyle='dashed', label="p = 0.05")
plt.axvline(-1, color='blue', linestyle='dashed', label="Downregulated (log2 FC < -1)")
plt.axvline(1, color='green', linestyle='dashed', label="Upregulated (log2 FC > 1)")
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10(p-value)")
plt.title("Volcano Plot of Differential Gene Expression")
plt.legend()
volcano_plot_path = "volcano_plot.png"
plt.savefig(volcano_plot_path, dpi=300)
plt.show()

# ----- 2. PCA Plot -----
pca = PCA(n_components=2)
X_pca = pca.fit_transform(df.T)

plt.figure(figsize=(8, 6))
sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1], color="blue")
plt.xlabel("PCA Component 1")
plt.ylabel("PCA Component 2")
plt.title("PCA of Gene Expression")
pca_plot_path = "pca_plot.png"
plt.savefig(pca_plot_path, dpi=300)
plt.show()

# ----- 3. MA Plot -----
mean_expression = df.mean(axis=1)
plt.figure(figsize=(8, 6))
plt.scatter(mean_expression, log2_fold_change, alpha=0.7)
plt.axhline(0, color='red', linestyle='dashed')
plt.xlabel("Mean Expression Level")
plt.ylabel("Log2 Fold Change")
plt.title("MA Plot of Gene Expression")
ma_plot_path = "ma_plot.png"
plt.savefig(ma_plot_path, dpi=300)
plt.show()

# ----- 4. Boxplot -----
plt.figure(figsize=(10, 6))
sns.boxplot(data=df.T)
plt.xticks(rotation=90)
plt.title("Gene Expression Boxplot")
boxplot_path = "boxplot.png"
plt.savefig(boxplot_path, dpi=300)
plt.show()

# Return file paths for user download
volcano_plot_path, pca_plot_path, ma_plot_path, boxplot_path
