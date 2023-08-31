## Attempt to run High MOI experimental from Sceptre
## https://katsevich-lab.github.io/sceptre/reference/run_sceptre_highmoi_experimental.html

library(sceptre)
library(Seurat)
library(dplyr)

outname <- "sceptre_new_targets_nG50_no_dialuplib_new_ref"

screen.rna <- readRDS("screen.rna.rds")
# gene matrix
response_matrix <- GetAssayData(screen.rna, slot = "counts")
# guides matrix
gRNA_matrix <- t(screen.rna@meta.data[,7:730])
gRNA_matrix <- gRNA_matrix[rowSums(gRNA_matrix) > 0,]

# this maps the guide RNAs to their targets
grna_group_data_frame <- read.table("data/gRNAs_targets_new_version.tsv", head=TRUE, sep="\t")

#co-variates
lg_gene_lib_size <- log(colSums(response_matrix)+0.0001)
lg_gRNA_lib_size <- log(colSums(gRNA_matrix)+0.0001)
#lg_dialup_lib_size<-log(colSums(dialup_matrix)+0.0001)
lg_crop_numi <- log(screen.rna@meta.data$crop_nUMI+0.0001)

binary_gRNA_matrix <- ifelse(gRNA_matrix >= 10, 1, 0)

covariate_data_frame <- data.frame(lg_gene_lib_size = lg_gene_lib_size, 
                                   lg_gRNA_lib_size = lg_gRNA_lib_size, 
                                   #lg_dialup_lib_size = lg_dialup_lib_size, 
                                   lg_crop_numi = lg_crop_numi, 
                                   num_guides = colSums(binary_gRNA_matrix),
                                   p_mito = screen.rna@meta.data$percent.mt, 
                                   batch = screen.rna@meta.data$orig.ident)

# discovery set to test - genes versus which targets
gene_gRNA_group_pairs <- read.table("data/genes_gRNAs_targets_new_version.tsv", head=TRUE, sep="\t")

# creates covariate object
formula_object <- formula(~lg_gene_lib_size + 
                            lg_gRNA_lib_size +
                            lg_crop_numi +
                            num_guides +
                            p_mito +
                            batch)

calibration_result <- run_sceptre_highmoi_experimental(
  response_matrix = response_matrix,
  grna_matrix = gRNA_matrix,
  covariate_data_frame = covariate_data_frame,
  grna_group_data_frame = grna_group_data_frame,
  response_grna_group_pairs = gene_gRNA_group_pairs,
  formula_object = formula_object,
  side = "both",
  grna_assign_threshold = 10,
  calibration_check = TRUE
)


plot_calibration_result(calibration_result)

discovery_set <- run_sceptre_highmoi_experimental(
  response_matrix = response_matrix,
  grna_matrix = gRNA_matrix,
  covariate_data_frame = covariate_data_frame,
  grna_group_data_frame = grna_group_data_frame,
  response_grna_group_pairs = gene_gRNA_group_pairs,
  formula_object = formula_object,
  side = "both",
  grna_assign_threshold = 10,
  calibration_check = FALSE
)

# negative log 10 is nicer for viewing in tracks
discovery_set$neg_log10<--log10(discovery_set$p_value)
discovery_set <- discovery_set %>% mutate(num_str = sapply(strsplit(as.character(grna_group), "_"), `[`, 2),
                                          bedgraph_num1 = as.numeric(sapply(strsplit(num_str, "\\."), `[`, 1)),
                                          bedgraph_num2 = as.numeric(sapply(strsplit(num_str, "\\."), `[`, 2)))

# Adjust the value in neg_log10 to be positive or negative depending on fold change
discovery_set$neg_log10_adj <- ifelse(discovery_set$log_2_fold_change < 0, discovery_set$neg_log10 * -1, discovery_set$neg_log10)

# Getting unique gene_ids
unique_gene_ids <- unique(discovery_set$response_id)

# set is sorted by p-value when returned so need to resort by genome coordinates
discovery_set <- discovery_set %>% arrange(bedgraph_num1)

# Writing to different files
for (gene_id in unique_gene_ids) {
  file_name <- paste0(outname, gene_id, ".bedGraph")
  
  # Writing header first
  header <- "track type=bedGraph"
  writeLines(header, file_name)
  
  # Then writing data
  df_gene_id <- discovery_set[discovery_set$response_id == gene_id, ]
  df_gene_id$chr9<-"chr9"
  write.table(df_gene_id[, c("chr9","bedgraph_num1", "bedgraph_num2", "neg_log10_adj")], file = file_name, 
              append = TRUE, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
}


compare_calibration_and_discovery_results(
  calibration_result = calibration_result,
  discovery_result = discovery_set
)