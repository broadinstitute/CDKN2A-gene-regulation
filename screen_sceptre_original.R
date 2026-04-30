library(dplyr)
library(ggplot2)
library(cowplot)
library(Matrix)
library(Seurat)
library(sceptre)
library(ggpubr)

outname <- "sceptre_guides_nG50_new_ref_anril_guide_MTAP_DMRTA1"

screen.rna <- readRDS("data/screen.rna.rds")
gene_matrix <- GetAssayData(screen.rna, slot = "counts")
gRNA_matrix <- t(screen.rna@meta.data[,7:730])
gRNA_matrix <- gRNA_matrix[rowSums(gRNA_matrix) > 0,]

lg_gene_lib_size <- log(colSums(gene_matrix)+0.0001)
# removing this made QQ plot worse
lg_gRNA_lib_size <- log(colSums(gRNA_matrix)+0.0001)
# this covariate seems to induce a false inversion effect
# though the QQ plot looks very good
#lg_dialup_lib_size<-log(colSums(dialup_matrix)+0.0001)
lg_crop_numi <- log(screen.rna@meta.data$crop_nUMI+0.0001)

# binarize to get the number of guides
binary_gRNA_matrix <- ifelse(gRNA_matrix >= 10, 1, 0)

# might be good to add the number of nonzero entires per cell
# for both genes and gRNAs
covariate_data_frame <- data.frame(lg_gene_lib_size = lg_gene_lib_size, 
                                   lg_gRNA_lib_size = lg_gRNA_lib_size, 
                                   #lg_dialup_lib_size = lg_dialup_lib_size, 
                                   lg_crop_numi = lg_crop_numi, 
                                   num_guides = colSums(binary_gRNA_matrix),
                                   p_mito = screen.rna@meta.data$percent.mt, 
                                   batch = screen.rna@meta.data$orig.ident)

covariate_data_frame<-covariate_data_frame[which(covariate_data_frame$num_guides<=50),]

covariate_matrix <- covariate_data_frame
covariate_matrix <- covariate_matrix[covariate_matrix$lg_gRNA_lib_size > 0,]
covariate_matrix <- covariate_matrix[covariate_matrix$lg_gene_lib_size > 0,]

common_barcodes <- Reduce(intersect, list(colnames(gene_matrix), colnames(gRNA_matrix), rownames(covariate_matrix)))
gene_matrix <- gene_matrix[, common_barcodes]
gRNA_matrix <- gRNA_matrix[, common_barcodes]
covariate_matrix <- covariate_matrix[common_barcodes, ]

perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix,10)

# This is only used when running in target mode
gRNA_groups_table <- read.table("data/gRNAs_targets.tsv", head=TRUE, sep="\t")
combined_perturbation_matrix <- combine_perturbations(perturbation_matrix = perturbation_matrix,
                                                      gRNA_groups_table = gRNA_groups_table)
# If running in target mode, use "data/genes_gRNAs_targets_MTAP_DMRTA1.tsv" (below);
# else use "data/genes_gRNAs_guides.tsv" and pass perturbation_matrix instead of combined_perturbation_matrix.
gene_gRNA_group_pairs <- read.table("data/genes_gRNAs_targets_MTAP_DMRTA1.tsv", head=TRUE, sep="\t")
side <- "both"

# "full" because we ask for full_output, thus far this has not been useful
result_full <- run_sceptre_highmoi(gene_matrix = gene_matrix,
                                   # HERE change for target vs guide; if target use combined_perturbation_matrix
                                   combined_perturbation_matrix = combined_perturbation_matrix,
                                   covariate_matrix = covariate_matrix,
                                   gene_gRNA_group_pairs = gene_gRNA_group_pairs,
                                   side = side, 
                                   full_output = TRUE)

# save the result because in guide mode this takes forever and even longer if you add more genes
write.table(result_full, paste0(outname, ".tsv"), quote=F, sep="\t")

# negative control QQ plot should be all in gray, in practice has not been
neg_control_p_vals <- result_full |> dplyr::filter(pair_type == "negative_control") |> pull(p_value)
qq_plot <- make_qq_plot(neg_control_p_vals)
plot(qq_plot)

# sanity check that positive controls are significant (one p16 target was not)
pos_control_p_vals <- result_full |> filter(pair_type == "positive_control") |> pull(p_value)
pos_control_p_vals

# it can also be useful to look at the QQ plot with the candidate pairs
candidate_pair_p_vals <-  result_full |> filter(pair_type == "candidate") |> pull(p_value)
make_qq_plot(candidate_pair_p_vals)

candidate_pair_results <- result_full[,c(1:4,11)] |> filter(pair_type == "candidate")
candidate_pair_results_p_adj <- candidate_pair_results |>
  mutate(p_val_adj = p.adjust(p_value, method = "BH"))
# this is commented out because we like to look at all results
#discovery_set <- candidate_pair_results_p_adj |> filter(p_val_adj <= 0.1)
discovery_set <- candidate_pair_results_p_adj 
pos_control <- result_full |> filter(pair_type == "positive_control")
pos_control <- pos_control[,c(1:4,11)] 
pos_control$p_val_adj <- pos_control$p_value
# we add the positive controls into the discovery set for making graphs
discovery_set <- rbind(discovery_set, pos_control)

# negative log 10 is nicer for viewing in tracks
discovery_set$neg_log10<--log10(discovery_set$p_val_adj)
# if in target mode, need to get back the coordinates
discovery_set$gRNA_id <- gsub("p16_tss1", "sgRNA_21974644.21975049", discovery_set$gRNA_id)
discovery_set$gRNA_id <- gsub("p16_tss2", "sgRNA_21975050.21975381", discovery_set$gRNA_id)
discovery_set$gRNA_id <- gsub("p16_tss3", "sgRNA_21975473.21975815", discovery_set$gRNA_id)
discovery_set$gRNA_id <- gsub("p14_tss", "sgRNA_21994286.21994909", discovery_set$gRNA_id)

# Split the string in gRNA_id
# this is for target mode
discovery_set <- discovery_set %>% mutate(num_str = sapply(strsplit(as.character(gRNA_id), "_"), `[`, 2),
                    bedgraph_num1 = as.numeric(sapply(strsplit(num_str, "\\."), `[`, 1)),
                    bedgraph_num2 = as.numeric(sapply(strsplit(num_str, "\\."), `[`, 2)))

# this is for guide mode, where we want to use the cut site
#discovery_set <- discovery_set %>% mutate(cut_site = sapply(strsplit(as.character(gRNA_id), "_"), `[`, 3),
#                                          bedgraph_num1 = as.numeric(sapply(strsplit(cut_site, "\\."), `[`, 1))-5,
#                                          bedgraph_num2 = bedgraph_num1+10)



# Adjust the value in neg_log10 to be positive or negative depending on fold change
discovery_set$neg_log10_adj <- ifelse(discovery_set$log_fold_change < 0, discovery_set$neg_log10 * -1, discovery_set$neg_log10)

# Getting unique gene_ids
unique_gene_ids <- unique(discovery_set$gene_id)

# Writing to different files
for (gene_id in unique_gene_ids) {
  file_name <- paste0(outname, gene_id, ".bedGraph")
  
  # Writing header first
  header <- "track type=bedGraph"
  writeLines(header, file_name)
  
  # Then writing data
  df_gene_id <- discovery_set[discovery_set$gene_id == gene_id, ]
  df_gene_id$chr9<-"chr9"
  write.table(df_gene_id[, c("chr9","bedgraph_num1", "bedgraph_num2", "neg_log10_adj")], file = file_name, 
              append = TRUE, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
}

### The below is to make correlation plots for all the possible covariates
### Variable names have changed since first wrote it, did my best to correct but just be aware
# p15<-GetAssayData(object = screen.rna, slot = 'counts')["CDKN2B", ]
# p16<-GetAssayData(object = screen.rna, slot = 'counts')["p16.dialup", ]
# p14<-GetAssayData(object = screen.rna, slot = 'counts')["p14.dialup", ]
# screen.rna@meta.data$p15<-p15
# screen.rna@meta.data$p16<-p16
# screen.rna@meta.data$p14<-p14
# 
# covariate_data_frame$cells<-rownames(covariate_data_frame)
# covariate_data_frame <- merge(covariate_data_frame, new_seurat_obj@meta.data[c("cells", "p16", "p15", "p14")], by = "cells")
# 
# # Get the names of the columns in covariate_data_frame
# covariate_names <- setdiff(colnames(covariate_data_frame), c("cells", "p14", "p15", "p16"))
# 
# # Loop through each target column
# for (target in c("p16", "p15", "p14")) {
#   # Loop through each covariate column
#   for (covariate in covariate_names) {
#     # Create a data frame
#     df <- data.frame(x = covariate_data_frame[[covariate]], y = covariate_data_frame[[target]])
#     
#     # Create a ggplot
#     p <- ggplot(df, aes(x = x, y = y)) +
#       geom_point(alpha = 0.5, color = "red") +  # Scatter plot with points semi-transparent to show density
#       geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Add regression line
#       theme_minimal() +  # Minimal theme
#       labs(x = covariate, y = target)  # Add labels
#     
#     # Add correlation coefficient to the plot
#     p <- p + stat_cor(method = "pearson", label.x = Inf, label.y = Inf, hjust = 1, vjust = 1)
#     # Save the plot to a PDF or PNG file
#     filename <- paste(target, covariate, "plot50.pdf", sep = "_")  # For PDF
#     # filename <- paste(target, covariate, "plot.png", sep = "_")  # For PNG
#     ggsave(filename, plot = p, width = 7, height = 5) 
#     # Add the plot to the list
#     #plot_list[[paste(target, covariate, sep = "_")]] <- p
#   }
# }
# 
# 
# 
# 
