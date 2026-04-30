#######load required libraries######

library(dplyr)
library(Seurat)
library(ggplot2)
library(Matrix)

#' This function processes the CROP array, returning matrices with
#' guides per cell, specific genes per cell, and
#' the general crop number of umi per cell (including off target alignments)
#'
#' @param geo_prefix The GEO file prefix for this array (e.g. "Dialup-array7")
#' @param prefix The cell-id prefix for this array, will need to match what was done for the mRNA
#' @param high.quality.cells The names of the cells kept from the mRNA
#' @returns Three matrices, one that is cells by guides, one that is cells by genes, and one that is cells by crop_numi
processcrop <- function(geo_prefix, prefix, high.quality.cells) {
  crop.data <- read10x_flat(geo_prefix)
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

## --------------------------------------------------------------------------
## Set this to the directory containing the GEO processed files
## (mRNA_arrayN_{barcodes,features,matrix} and the dial-up triplets).
## See README §"Reproducing from GEO" for the file layout.
## --------------------------------------------------------------------------
geo_dir <- "geo_screen"

## Read 10x-style triplets that GEO ships as flat files (prefix_{barcodes,features,matrix}).
read10x_flat <- function(geo_prefix) {
  d <- tempfile(); dir.create(d)
  file.copy(file.path(geo_dir, paste0(geo_prefix, "_barcodes.tsv.gz")), file.path(d, "barcodes.tsv.gz"))
  file.copy(file.path(geo_dir, paste0(geo_prefix, "_features.tsv.gz")), file.path(d, "features.tsv.gz"))
  file.copy(file.path(geo_dir, paste0(geo_prefix, "_matrix.mtx.gz")),   file.path(d, "matrix.mtx.gz"))
  Read10X(data.dir = d)
}

#' Process an mRNA array. Filter by min cells 3 and min features 200, calculate
#' percent mitochondria, subset on 500 <= nGenes <= 6000, nUMI >= 750, percent.mt < 15.
processarray <- function(geo_prefix) {
  data <- read10x_flat(geo_prefix)
  singlearray <- CreateSeuratObject(counts = data, project = geo_prefix, min.cells = 3, min.features = 200)
  singlearray[["percent.mt"]] <- PercentageFeatureSet(singlearray, pattern = "^MT-")
  singlearray@meta.data <- singlearray@meta.data %>% dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
  singlearray <- subset(singlearray, subset = nGene > 499 & nGene < 6001 & nUMI > 749 & percent.mt < 15)
  singlearray$nCount_RNA <- NULL
  singlearray$nFeature_RNA <- NULL
  return(singlearray)
}

## mRNA libraries: GEO names mRNA_array1 ... mRNA_array16
mrna_arrays <- lapply(1:16, function(i) processarray(paste0("mRNA_array", i)))
names(mrna_arrays) <- paste0("a", sprintf("%02d", 1:16))

screen.rna <- merge(mrna_arrays[[1]], y = mrna_arrays[-1],
                    add.cell.ids = names(mrna_arrays), project = "screen")

# Add cell IDs to metadata
screen.rna@meta.data$cells <- rownames(screen.rna@meta.data)

## Dial-up / CROP libraries: arrays 1-6 are named new-screen-CROP-orig-stitch-arrayN,
## arrays 7-16 are named Dialup-arrayN on GEO.
dialup_prefixes <- c(paste0("new-screen-CROP-orig-stitch-array", 1:6),
                     paste0("Dialup-array", 7:16))
cell_prefixes <- paste0("a", sprintf("%02d", 1:16), "_")
hq_cells <- screen.rna@meta.data$cells

# Apply processcrop function on each dialup array
crop_outputs <- lapply(seq_along(dialup_prefixes), function(i) {
  processcrop(geo_prefix = dialup_prefixes[i], prefix = cell_prefixes[i], high.quality.cells = hq_cells)
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

saveRDS(screen.rna, "data/screen.rna.rds")

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