library(dplyr)

outname <- "sceptre_guides_new_version_mixture_corrected"
# Read in the results file
discovery_set <- read.table(paste0(outname, ".tsv"), header = TRUE, sep = "\t")

print(paste("Total pairs in file:", nrow(discovery_set)))
print(paste("Pairs with p_val_adj <= 0.05:", sum(discovery_set$p_val_adj <= 0.05, na.rm = TRUE)))

# Filter for significant results (p_val_adj <= 0.05)
discovery_set_sig <- discovery_set %>%
  filter(!is.na(p_val_adj) & p_val_adj <= 0.05)

print(paste("Significant pairs to plot:", nrow(discovery_set_sig)))

# Create negative log10 adjusted values
discovery_set_sig$neg_log10 <- -log10(discovery_set_sig$p_val_adj)
discovery_set_sig$neg_log10_adj <- ifelse(discovery_set_sig$log_2_fold_change < 0, 
                                           discovery_set_sig$neg_log10 * -1, 
                                           discovery_set_sig$neg_log10)

# Replace target names with coordinates if needed
discovery_set_sig$grna_target <- gsub("p16_tss1", "sgRNA_21974644.21975049", discovery_set_sig$grna_target)
discovery_set_sig$grna_target <- gsub("p16_tss2", "sgRNA_21975050.21975381", discovery_set_sig$grna_target)
discovery_set_sig$grna_target <- gsub("p16_tss3", "sgRNA_21975473.21975815", discovery_set_sig$grna_target)
discovery_set_sig$grna_target <- gsub("p14_tss", "sgRNA_21994286.21994909", discovery_set_sig$grna_target)

# Extract coordinates for bedGraph (guide mode - using cut sites)
discovery_set_sig <- discovery_set_sig %>% 
  mutate(cut_site = sapply(strsplit(as.character(grna_target), "_"), `[`, 3),
         bedgraph_num1 = as.numeric(sapply(strsplit(cut_site, "\\."), `[`, 1)) - 5,
         bedgraph_num2 = bedgraph_num1 + 10,
         chr9 = "chr9")

# Remove rows with NA coordinates
discovery_set_sig <- discovery_set_sig %>%
  filter(!is.na(bedgraph_num1) & !is.na(bedgraph_num2))

print(paste("Pairs with valid coordinates:", nrow(discovery_set_sig)))

# Write sorted bedGraph files for each gene (no overlaps)
unique_gene_ids <- unique(discovery_set_sig$response_id)
print(paste("Writing bedGraph files for", length(unique_gene_ids), "genes"))

for (gene_id in unique_gene_ids) {
  file_name <- paste0(outname, "_", gene_id, ".bedGraph")
  
  # Filter and sort by position
  df_gene_id <- discovery_set_sig %>%
    filter(response_id == gene_id) %>%
    arrange(bedgraph_num1)
  
  # Check for and handle overlaps
  if (nrow(df_gene_id) > 1) {
    for (i in 2:nrow(df_gene_id)) {
      if (df_gene_id$bedgraph_num1[i] < df_gene_id$bedgraph_num2[i-1]) {
        # Overlap detected - adjust start position
        df_gene_id$bedgraph_num1[i] <- df_gene_id$bedgraph_num2[i-1]
        df_gene_id$bedgraph_num2[i] <- df_gene_id$bedgraph_num1[i] + 1  # minimum 1bp interval
      }
    }
  }
  
  # Write header
  header <- "track type=bedGraph"
  writeLines(header, file_name)
  
  # Write data
  write.table(df_gene_id[, c("chr9", "bedgraph_num1", "bedgraph_num2", "neg_log10_adj", "grna_target")], 
              file = file_name, 
              append = TRUE, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
  print(paste("  Written:", file_name, "with", nrow(df_gene_id), "intervals"))
}

print("Done! bedGraph files created for significant results only (p_val_adj <= 0.05)")
