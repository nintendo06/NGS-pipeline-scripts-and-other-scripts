# DIFFERENTIAL EXPRESSION ANALYSIS IN PYTHON 
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import ols
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.packages import importr
deseq = importr('DESeq2')

# Step 1: Data Preparation
# Read count data into a pandas DataFrame
count_data = pd.read_csv('counts.txt', sep='\t', index_col=0)

# Create metadata table with conditions for each sample
metadata = pd.DataFrame({
    'Sample': count_data.columns,
    'Condition': ['tumor', 'tumor', 'normal', 'normal'] # Adjust according to your sample conditions
})

# Step 2: Normalization and Preprocessing
# Normalize the count data using DESeq2's variance stabilizing transformation (VST)
r_count_data = pandas2ri.py2ri(count_data)
dds = deseq.DESeqDataSetFromMatrix(countData=r_count_data, colData=metadata, design=~Condition)
dds = deseq.DESeq(dds)
vst = deseq.vst(dds)
vst_data = pandas2ri.ri2py(vst)

# Step 3: Statistical Testing
# Perform statistical test (e.g., t-test) to identify differentially expressed genes
design_matrix = deseq.modelMatrix(dds)
contrast_matrix = np.array([1, -1]).reshape(2,1)
contrasts = deseq.makeContrasts(condition_diff=contrast_matrix, levels=design_matrix)
contrast_fit = deseq.estimateContrast(dds, contrasts)
contrast_results = deseq.results(contrast_fit)
contrast_results = pandas2ri.ri2py(contrast_results)

# Step 4: Multiple Testing Correction
p_values = contrast_results['pvalue']
adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]
contrast_results['adjusted_pvalue'] = adjusted_p_values

# Step 5: Interpretation and Visualization
# Volcano Plot
plt.scatter(contrast_results['log2FoldChange'], -np.log10(contrast_results['adjusted_pvalue']), alpha=0.5)
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(Adjusted p-value)')
plt.title('Volcano Plot')
plt.show()

# Heatmap
sns.clustermap(vst_data, cmap='viridis')

# Gene Set Enrichment Analysis
gene_ids = count_data.index
de_genes = contrast_results.loc[contrast_results['adjusted_pvalue'] < 0.05, 'baseMean'].index.tolist()
de_genes = np.intersect1d(de_genes, gene_ids) # Ensure genes are present in original count data

# Perform gene set enrichment analysis using libraries like clusterProfiler or fgsea

# Principal Component Analysis (PCA)
scaler = StandardScaler()
scaled_data = scaler.fit_transform(vst_data.T)
pca = PCA(n_components=2)
pca_data = pca.fit_transform(scaled_data)
pca_df = pd.DataFrame(data=pca_data, columns=['PC1', 'PC2'])
pca_df['Condition'] = metadata['Condition']
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Condition')

# Step 6: Additional Analyses
# Additional downstream analysis such as pathway analysis, gene network analysis, etc.


VOLCANO PLOT 
import pandas as pd
import matplotlib.pyplot as plt

# Read the differential expression results file
df = pd.read_csv("differential_expression_results.txt", sep="\t")

# Calculate the negative logarithm of p-values
df['log10p'] = -1 * df['p-value'].apply(lambda x: math.log10(x))

# Set the significance threshold
significance_threshold = 0.05

# Create a boolean mask to identify significantly differentially expressed genes
sig_genes = df['p-value'] < significance_threshold
