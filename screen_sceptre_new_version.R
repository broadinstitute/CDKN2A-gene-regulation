library(dplyr)
library(ggplot2)
library(cowplot)
library(Matrix)
library(Seurat)
library(sceptre)
library(ggpubr)

outname <- "sceptre_guides_new_version_threshold"

# Load your Seurat object
screen.rna <- readRDS("data/screen.rna.rds")

# Extract matrices (updated for Seurat v5)
gene_matrix <- GetAssayData(screen.rna, layer = "counts")
gRNA_matrix <- t(screen.rna@meta.data[, 7:730])
gRNA_matrix <- gRNA_matrix[rowSums(gRNA_matrix) > 0, ]

# Filter cells with num_guides <= 50
binary_gRNA_matrix <- ifelse(gRNA_matrix >= 10, 1, 0)
num_guides <- colSums(binary_gRNA_matrix)
cells_to_keep <- names(num_guides)[num_guides <= 50]

# Subset matrices
gene_matrix <- gene_matrix[, cells_to_keep]
gRNA_matrix <- gRNA_matrix[, cells_to_keep]

# Create covariate matrix (sceptre will auto-compute some of these, but we add custom ones)
extra_covariates <- data.frame(
  lg_crop_numi = log(screen.rna@meta.data[cells_to_keep, "crop_nUMI"] + 0.0001),
  num_guides = num_guides[cells_to_keep],
  batch = screen.rna@meta.data[cells_to_keep, "orig.ident"]
)

# Load gene-gRNA pairs
gene_gRNA_pairs_raw <- read.table("data/genes_gRNAs_guides_MTAP_DMRTA1_IFN_corrected.tsv", 
                                   header = TRUE, sep = "\t")

# Filter out negative_control rows - sceptre will create these automatically
gene_gRNA_pairs_filtered <- gene_gRNA_pairs_raw %>%
  filter(pair_type != "negative_control")

# Create grna_target_data_frame for sceptre
# This maps each gRNA to its target (in guide mode, each guide is its own target)
# Map safe harbor gRNAs to "non-targeting"
grna_ids <- rownames(gRNA_matrix)
grna_targets <- ifelse(grepl("safe.harbor", grna_ids, ignore.case = TRUE), 
                       "non-targeting", 
                       grna_ids)

grna_target_data_frame <- data.frame(
  grna_id = grna_ids,
  grna_target = grna_targets
)

# Print summary
print(paste("Non-targeting gRNAs found:", sum(grna_target_data_frame$grna_target == "non-targeting")))
print(paste("Targeting gRNAs:", sum(grna_target_data_frame$grna_target != "non-targeting")))

# Import data into sceptre_object
sceptre_object <- import_data(
  response_matrix = gene_matrix,
  grna_matrix = gRNA_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "high",
  extra_covariates = extra_covariates
)

# Construct pairs to test (excluding negative controls)
# Your file has columns: gene_id, gRNA_group, pair_type
pos_control_pairs <- gene_gRNA_pairs_filtered %>%
  filter(pair_type == "positive_control") %>%
  select(grna_target = gRNA_group, response_id = gene_id)

discovery_pairs <- gene_gRNA_pairs_filtered %>%
  filter(pair_type == "candidate") %>%
  select(grna_target = gRNA_group, response_id = gene_id)

# Set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = pos_control_pairs,
  side = "both",  # You had this as "both" in your original script
  grna_integration_strategy = "singleton"  # Important: keeps gRNAs separate (guide mode)
)

# Assign gRNAs
# Option 1: Use mixture method (recommended for high-MOI, accounts for background)
#sceptre_object <- assign_grnas(
#  sceptre_object = sceptre_object,
#  method = "mixture",
#  parallel = TRUE
#)

# Option 2: Use thresholding method 
sceptre_object <- assign_grnas(
   sceptre_object = sceptre_object,
   method = "thresholding",
   threshold = 10,
   parallel = TRUE
 )

# Plot gRNA assignments
plot(sceptre_object)
ggsave(paste0(outname, "_grna_assignments.pdf"), width = 10, height = 8)

# Run calibration check
sceptre_object <- run_calibration_check(
  sceptre_object = sceptre_object,
  parallel = TRUE
)

# Plot calibration check
plot(sceptre_object)
ggsave(paste0(outname, "_calibration_check.pdf"), width = 10, height = 8)

# Run power check (on positive controls)
sceptre_object <- run_power_check(
  sceptre_object = sceptre_object,
  parallel = TRUE
)

# Plot power check
plot(sceptre_object)
ggsave(paste0(outname, "_power_check.pdf"), width = 10, height = 8)

# Run discovery analysis
sceptre_object <- run_discovery_analysis(
  sceptre_object = sceptre_object,
  parallel = TRUE
)

# Plot discovery results
plot(sceptre_object)
ggsave(paste0(outname, "_discovery.pdf"), width = 10, height = 8)

# Get results
result_full <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_discovery_analysis"
)

# Check for NA results
print("=== Discovery Analysis Summary ===")
print(paste("Total pairs:", nrow(result_full)))
print(paste("Pairs with NA p-values:", sum(is.na(result_full$p_value))))
print(paste("Pairs that passed QC:", sum(result_full$pass_qc)))
if(sum(!result_full$pass_qc) > 0) {
  print("Example pairs that failed QC:")
  print(head(result_full %>% filter(!pass_qc)))
}

# Add positive control results
pos_control_results <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_power_check"
)

# Combine and adjust p-values
# Add pair_type column to distinguish results
candidate_results <- result_full %>%
  mutate(pair_type = "candidate",
         p_val_adj = p.adjust(p_value, method = "BH"))

pos_control_results <- pos_control_results %>%
  mutate(pair_type = "positive_control",
         p_val_adj = p_value)

discovery_set <- bind_rows(candidate_results, pos_control_results)

# Save full results
write.table(discovery_set, paste0(outname, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

# Create negative log10 adjusted values
discovery_set$neg_log10 <- -log10(discovery_set$p_val_adj)
discovery_set$neg_log10_adj <- ifelse(discovery_set$log_2_fold_change < 0, 
                                       discovery_set$neg_log10 * -1, 
                                       discovery_set$neg_log10)

# Replace target names with coordinates if needed
discovery_set$grna_target <- gsub("p16_tss1", "sgRNA_21974644.21975049", discovery_set$grna_target)
discovery_set$grna_target <- gsub("p16_tss2", "sgRNA_21975050.21975381", discovery_set$grna_target)
discovery_set$grna_target <- gsub("p16_tss3", "sgRNA_21975473.21975815", discovery_set$grna_target)
discovery_set$grna_target <- gsub("p14_tss", "sgRNA_21994286.21994909", discovery_set$grna_target)

# Extract coordinates for bedGraph (guide mode - using cut sites)
discovery_set <- discovery_set %>% 
  mutate(cut_site = sapply(strsplit(as.character(grna_target), "_"), `[`, 3),
         bedgraph_num1 = as.numeric(sapply(strsplit(cut_site, "\\."), `[`, 1)) - 5,
         bedgraph_num2 = bedgraph_num1 + 10)

# Write bedGraph files for each gene
unique_gene_ids <- unique(discovery_set$response_id)

for (gene_id in unique_gene_ids) {
  file_name <- paste0(outname, "_", gene_id, ".bedGraph")
  
  header <- "track type=bedGraph"
  writeLines(header, file_name)
  
  df_gene_id <- discovery_set[discovery_set$response_id == gene_id, ]
  df_gene_id$chr9 <- "chr9"
  
  write.table(df_gene_id[, c("chr9", "bedgraph_num1", "bedgraph_num2", "neg_log10_adj")], 
              file = file_name, 
              append = TRUE, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
}

# Write all outputs to directory
write_outputs_to_directory(sceptre_object, paste0(outname, "_outputs"))

print("Analysis complete!")
print(paste("Significant discoveries:", sum(discovery_set$significant, na.rm = TRUE)))



