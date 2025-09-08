adata = sc.read_h5ad('/seq/epiprod02/etorlait/p16/Published_Data/Gray_et_al_2022/Remap_VersionC_Gtf_Terra/Output_Cumulus/LowGeneMin.unknown-rna.h5ad')
df = adata[:,adata.var.index.isin(['CDKN2A_p14_ex1','CDKN2A_p16_ex1','CDKN2A_p14_p16'])].to_df()
# df_dialup has barcode,cellSubtype,dialup_p16_ratio_nCounts,dialup_p14_ratio_nCounts
# it's the number of umis mapping to p16 (resp p14) divided by the number of UMIs for that barcode
df_dialup = pd.read_csv("dialup_ratios.csv", index_col=0)
adata_original = sc.read_h5ad('/Volumes/seq_epiprod02/etorlait/p16/Published_Data/Gray_et_al_2022/Figure1B-2A_allCells_SingleCellExperiment_converted_to_Seurat_then_to_h5ad_withCoords.h5ad')
barcode_names = adata_original.obs.index.tolist()
# Extract sample IDs from df (AnnData-derived expression data)
df_sample_ids = df.index.str.split('-').str[0].unique()
# Extract sample IDs from dialup_df (Dial-up ratio data)
dialup_sample_ids = dfdialup.index.str.split('').str[0].unique()
sample_mapping = {
    "RM_A": "Human-WT-A", "PM_G": "Human-BRCA2-A", "PM_A": "Human-BRCA1-C", "PM_C": "Human-BRCA1-D",
    "PM_K": "Human-BRCA2-B", "PM_I": "Human-BRCA2-C", "PM_J": "Human-BRCA2-D", "RM_D": "Human-WT-C",
    "RM_C": "Human-WT-D", "PM_F": "Human-BRCA1-H", "PM_L": "Human-BRCA2-E"
}
# Apply the sample mapping to the 'new_barcode' column
def map_barcode(barcode):
    parts = barcode.split('-')
    if parts[0] in sample_mapping:
        return samplemapping[parts[0]] + '' + '_'.join(parts[1:])
    else:
        return barcode

# Create a new column 'new_barcode' in df
df['new_barcode'] = df.index
df['new_barcode'] = df['new_barcode'].apply(map_barcode)
# Merge df and df_dialup based on 'new_barcode' and index, respectively.
merged_df = pd.merge(df, df_dialup, left_on='new_barcode', right_index=True, how='inner')
# Rename columns in merged_df
merged_df = merged_df.rename(columns={
    'CDKN2A_p16_ex1': 'CDKN2A-p16',
    'CDKN2A_p14_ex1': 'CDKN2A-p14',
    'CDKN2A_p14_p16': 'CDKN2A-p14+p16',
    'dialup_p16_ratio_nCounts': 'dialup-p16',
    'dialup_p14_ratio_nCounts': 'dialup-p14'
})
# Create a new column 'sample' from the 'barcodekey' column
merged_df['sample'] = merged_df['newbarcode'].str.split('').str[0]
cdkn2a_expression = adata_original[:, 'CDKN2A'].X.toarray().flatten()
# Create a DataFrame from the expression data:
cdkn2a_df = pd.DataFrame({'CDKN2A': cdkn2a_expression}, index=adata_original.obs.index)
# Merge the DataFrames:
merged_df = merged_df.merge(cdkn2a_df, left_on='new_barcode', right_index=True, how='left')
genes = ['CDKN2A-p16', 'CDKN2A-p14', 'CDKN2A', 'dialup-p16', 'dialup-p14']
# Calculate the overall percentage of cells expressing each gene (i.e., where expression > 0)
total_cells = len(merged_df)
percentages = [np.sum(merged_df[gene] > 0) / total_cells * 100 for gene in genes]
# Calculate the standard error of the mean (SEM) per sample for each gene
sems = [merged_df.groupby('sample')[gene].apply(lambda x: (x > 0).mean()).sem() * 100 for gene in genes]
# Prepare the data for plotting
labels = genes
means = percentages
errors = sems
# Plot the overall bar chart with error bars
plt.figure(figsize=(8, 6))
bar_width = 0.5
plt.bar(labels, means, yerr=errors, capsize=5, color = ['blue', 'green', 'gray','blue', 'green'])
plt.xticks(x_positions, x_labels, rotation=45, ha='right')
plt.ylabel('Percentage of Cells (%)')
plt.tight_layout()
plt.show()
