#######load required libraries######

library(dplyr)
library(Seurat)
library(ggplot2)
library(Matrix)

#' This function processes the mRNA array. It filters by min cells 3 and min
#' features 200, calculates percent mitochondria, then subsets on  
#' 500 <= nGenes <= 6000, nUMI >= 750, and percent mitochondra < 15 
#' 
#' @param filename The base directory, should contain 10X style barcodes, features, matrix
#' @returns Seurat object with filtered gene by cell data
processarray <- function(filename) {
  data <- Read10X(data.dir = filename)
  singlearray <- CreateSeuratObject(counts = data, project = basename(filename), min.cells = 3, min.features = 200)
  singlearray[["percent.mt"]] <- PercentageFeatureSet(singlearray, pattern = "^MT-")
  singlearray@meta.data <- singlearray@meta.data %>% dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
  
  singlearray <- subset(singlearray, subset = nGene > 499 & nGene < 6001 & nUMI > 749 & percent.mt < 15)
  # subset adds these back in for some reason
  singlearray$nCount_RNA <- NULL
  singlearray$nFeature_RNA <- NULL
  # Sceptre does not use normalized counts so we don't do this step, as it creates an enormous object
  #singlearray<- SCTransform(singlearray, verbose = TRUE)
  return(singlearray)
} 

#' This function processes the CROP array, returning matrices with 
#' guides per cell, specific genes per cell, and
#' the general crop number of umi per cell (including off target alignments) 
#' 
#' @param filename The base directory for CROP sequencing, should contain 10X style barcodes, features, matrix
#' @param prefix The prefix for this array, will need to match what was done for the mRNA
#' @param high.quality.cells The names of the cells kept from the mRNA
#' @returns Three matrices, one that is cells by guides, one that is cells by genes, and one that is cells by crop_numi
processcrop <- function(filename, prefix, high.quality.cells) {
  crop.data <- Read10X(data.dir = filename)
  n <- paste0(prefix,colnames(crop.data))
  inds <- which(n %in% high.quality.cells)
  mytable <- tail(crop.data, n=724)[,inds]
  n <- paste0(prefix,colnames(mytable))
  colnames(mytable) <- n
  mytable <- t(as.matrix(mytable))
  
  mytable1<-crop.data[which(rownames(crop.data)=="p16"):which(rownames(crop.data)=="ANRIL_long"),inds]
  n<-paste0(prefix,colnames(mytable1))
  colnames(mytable1)<-n
  mytable1 <- t(as.matrix(mytable1))
  
  mytable2<-crop.data[,inds]
  n<-paste0(prefix,colnames(mytable2))
  colnames(mytable2)<-n

  crop1 <- CreateSeuratObject(counts = mytable2)
  crop1@meta.data <- crop1@meta.data %>% dplyr::rename(crop_nUMI = nCount_RNA)
  crop1@meta.data$nFeature_RNA <- NULL
  crop1@meta.data$cells <- rownames(crop1@meta.data)
  
  return(list(guides = mytable, dialup_genes = mytable1,crop_meta = crop1@meta.data))
}

array01 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array01")
array02 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array02")
array03 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array03")
array04 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array04")
array05 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array05")
array06 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array06")
array07 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array07")
array08 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array08")
array09 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array09")
array10 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array10")
array11 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array11")
array12 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array12")
array13 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array13")
array14 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array14")
array15 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array15")
array16 <- processarray("/Users/neva/Documents/calico/mRNA_standard_reference/array16")

screen.rna <- merge(array01, y = c(array02, array03,array04,array05,array06,array07,array08,array09,array10,array11,array12,array13,array14,array15,array16), add.cell.ids = c("a01", "a02", "a03","a04", "a05", "a06", "a07","a08","a09","a10","a11","a12","a13","a14","a15","a16"), project = "screen")

# Add cell IDs to metadata
screen.rna@meta.data$cells <- rownames(screen.rna@meta.data)

# Dialup filenames
filenames <- c("/Users/neva/Documents/calico/newscreen/new-screen-CROP-orig-stitch-array1/", 
               "/Users/neva/Documents/calico/newscreen/new-screen-CROP-orig-stitch-array2/",
               "/Users/neva/Documents/calico/newscreen/new-screen-CROP-orig-stitch-array3/",
               "/Users/neva/Documents/calico/newscreen/new-screen-CROP-orig-stitch-array4/",
               "/Users/neva/Documents/calico/newscreen/new-screen-CROP-orig-stitch-array5/",
               "/Users/neva/Documents/calico/newscreen/new-screen-CROP-orig-stitch-array6/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array7/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array8/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array9/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array10/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array11/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array12/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array13/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array14/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array15/",
               "/Users/neva/Documents/calico/newscreen/Dialup-array16/")
prefixes <- c("a01_", "a02_", "a03_", "a04_", "a05_", "a06_", "a07_", "a08_", "a09_", "a10_", "a11_", "a12_","a13_", "a14_", "a15_", "a16_") 
hq_cells <- screen.rna@meta.data$cells

# Apply processcrop function on each dialup array
crop_outputs <- lapply(1:length(filenames), function(i) {
  processcrop(filename = filenames[i], prefix = prefixes[i], high.quality.cells = hq_cells)
})

# Concatenate tables for each of 'guides', 'dialup_genes', and 'crop_meta'
final_guides <- do.call(rbind, lapply(crop_outputs, function(x) x$guides))
final_dialup_genes <- do.call(rbind, lapply(crop_outputs, function(x) x$dialup_genes))
final_crop_meta <- do.call(rbind, lapply(crop_outputs, function(x) x$crop_meta))

colnames(final_dialup_genes) <- paste0(colnames(final_dialup_genes), ".dialup")

#write.table(final_guides, "guides_cells.tsv",col.names = NA, quote=F, sep="\t")
#write.table(final_dialup_genes, "p16_cells_evennew.tsv",col.names = NA, quote=F, sep="\t")

data <- GetAssayData(screen.rna, slot = "counts")
data2<-t(final_dialup_genes)
gene_matrix <- rbind(data, data2)
screen.rna <- CreateSeuratObject(counts = gene_matrix)
screen.rna[["percent.mt"]] <- PercentageFeatureSet(screen.rna, pattern = "^MT-")
screen.rna@meta.data$cells<-rownames(screen.rna@meta.data)

md <- left_join(screen.rna@meta.data, final_crop_meta, by = "cells")
screen.rna@meta.data$crop_nUMI<-md$crop_nUMI

## Can add dialup as columns, but we've already added it to the gene matrix
#final_dialup_genes <- as.data.frame(final_dialup_genes)
#final_dialup_genes$cells <- rownames(final_dialup_genes)
#md <- left_join(screen.rna@meta.data, final_dialup_genes, by = "cells")
#screen.rna@meta.data$p14.dialup <- md$p14
#screen.rna@meta.data$p16.dialup <- md$p16
#screen.rna@meta.data$ANRIL_long.dialup <- md$ANRIL_long
#screen.rna@meta.data$ANRIL_short.dialup <- md$ANRIL_short

final_guides <- as.data.frame(final_guides)
final_guides$cells <- rownames(final_guides)
metadata2 <- left_join(screen.rna@meta.data, final_guides, by = "cells")
rownames(metadata2)<-metadata2$cells
screen.rna <- AddMetaData(screen.rna, metadata = metadata2)

saveRDS(screen.rna, "screen.rna.rds")

################# QC Plots and Clustering
#########
#screen.rna@meta.data %>% ggplot(aes(x=nUMI)) + geom_density(alpha = 0.2,color="darkred", fill="red") + scale_x_log10() + theme_classic() + ylab("log10 cell density")
#screen.rna@meta.data %>% ggplot(aes(x=nGene)) + geom_density(alpha = 0.2,color="darkred", fill="red") + theme_classic() + scale_x_log10() + ylab("log10 cell density")
#screen.rna@meta.data %>% ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + geom_point() + scale_colour_gradient(low = "gray90", high = "black") + stat_smooth(method=lm) + scale_x_log10() + scale_y_log10() + theme_classic() + geom_vline(xintercept = 750, linetype="dashed") + geom_hline(yintercept = 500, linetype="dashed") + geom_hline(yintercept = 6000, linetype="dashed") 
#VlnPlot(screen.rna, features = "nGene", pt.size=0)
#VlnPlot(screen.rna, features = "nUMI", pt.size=0)
#VlnPlot(screen.rna, features ="percent.mt", pt.size=0)


## Should lognormalize before the following

## Run PCA, autochoose dimensions
#screen.rna <- RunPCA(screen.rna)

## Choose principal components to consider via two metrics:
## The point where the principal components only contribute 5% of standard deviation and
## the principal components cumulatively contribute 90% of the standard deviation.
#pct <- screen.rna[["pca"]]@stdev / sum(screen.rna[["pca"]]@stdev)*100
#cumu <- cumsum(pct)
#co1 <- which(cumu > 90 & pct < 5)[1]

## The point where the percent change in variation between the consecutive PCs is less than 0.1%.
#co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# Minimum of the two calculations
#pcs <- min(co1,co2)

######### Generate the Elbow plot to select number of dimensions #############
#pdf("screen.rna_ElbowPlot.pdf", height = 6, width = 6)
# Create a dataframe with values
#plot_df <- data.frame(pct = pct, 
#                      cumu = cumu, 
#                      rank = 1:length(pct))

# Elbow plot to visualize 
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw() +
#  ylab("Percent of std dev explained by PC") +
#  xlab("Cumulative std dev")
#dev.off()


######### Run UMAP and clustering algorithm, Select dimensions based on Elbow plot ############
#screen.rna <- RunUMAP(screen.rna, dims = 1:pcs)
#screen.rna <- FindNeighbors(screen.rna, dims = 1:pcs)
#screen.rna <- FindClusters(screen.rna, resolution=0.3)

#pdf("screen.rna_UMAP_Clusters.pdf", height = 6, width = 7)
#DimPlot(screen.rna, reduction="umap")
#dev.off()
#pdf("screen.rna_UMAP_Samples.pdf", height = 6, width = 7)
#DimPlot(screen.rna, reduction="umap", group.by = "Timepoints", cols=c("pink", "goldenrod3", "limegreen", "deepskyblue2"))
#dev.off()

#FeaturePlot(screen.rna, features = "MKI67")