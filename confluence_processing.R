### This starts out as how to create the Multiome object
### from the confluence experiment and then devolves into
### a bunch of scratch code, some of which is how we created
### the atac peaks correlated with p16 or p14.
###
#######load required libraries######

library(dplyr)
library(Seurat)
library(ggplot2)
library(Matrix)

########## Read 10X matrix and Create a Seurat object ###############

data <- Read10X(data.dir = "/Users/neva/Documents/calico/all-Broad-expression")

Multiome.rna <- CreateSeuratObject(counts = data, project = "Multiome_RNA", min.cells = 3, min.features = 200)

Multiome.rna[["percent.mt"]] <- PercentageFeatureSet(Multiome.rna, pattern = "^MT-")

# Add number of genes per UMI for each cell to metadata
Multiome.rna$log10GenesPerUMI <- log10(Multiome.rna$nFeature_RNA) / log10(Multiome.rna$nCount_RNA)
# Compute percent mito ratio
Multiome.rna$mitoRatio <- PercentageFeatureSet(object = Multiome.rna, pattern = "^MT-")
Multiome.rna$mitoRatio <- Multiome.rna@meta.data$mitoRatio / 100
# Add cell IDs to metadata
Multiome.rna@meta.data$cells <- rownames(Multiome.rna@meta.data)

# Rename columns
Multiome.rna@meta.data <- Multiome.rna@meta.data %>%
        dplyr::rename(nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

ttt<-Multiome.rna@meta.data
ttt<-left_join(ttt, all_atac, by="cells")
all_atac<-read.table("ppp", header=TRUE)
all_atac$cells <- all_atac$cell
all_atac$cell <- NULL
ttt<-left_join(ttt, all_atac, by="cells")
mytable2 <- read.table("all-p16-p14-Broad-dialout.tsv", head = TRUE)
ttt<-left_join(ttt, mytable2, by="cells")
Multiome.rna@meta.data$TSS<-ttt$TSS
Multiome.rna@meta.data$nFrags<-ttt$nFrags
Multiome.rna@meta.data$p16.dialout<-ttt$p16.dialout
Multiome.rna@meta.data$p14.dialout<-ttt$p14.dialout
Multiome.rna@meta.data$ANRIL_short.dialout<-ttt$ANRIL_short
Multiome.rna@meta.data$ANRIL_long.dialout<-ttt$ANRIL_long


Multiome.rna <- subset(Multiome.rna, subset = nUMI >= 1000 & TSS >= 4 & nFrags >= 300 & nGene >=500 & percent.mt <=30)

# Visualize the number UMIs/transcripts per cell
pdf("Multiome.rna_UMIs_per_cell.pdf", height = 6, width = 9)
Multiome.rna@meta.data %>%
     ggplot(aes(x=nUMI)) + 
     geom_density(alpha = 0.2,color="darkred", fill="red") + 
     scale_x_log10() + 
     theme_classic() +
     ylab("log10 cell density")
dev.off()

# Visualize the distribution of genes detected per cell via histogram
pdf("Multiome.rna_genes_per_cell.pdf", height = 6, width = 9)
Multiome.rna@meta.data %>% 
  	ggplot(aes(x=nGene)) + 
  	geom_density(alpha = 0.2,color="darkred", fill="red") + 
  	theme_classic() +
  	scale_x_log10()
dev.off()

# Visualize the correlation between genes detected and number of UMIs and determine whether strong
# presence of cells with low numbers of genes/UMIs
pdf("Multiome.rna_mitoratio.pdf", height = 6, width = 9)
Multiome.rna@meta.data %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) 
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
#Multiome.rna@meta.data %>% 
#  	ggplot(aes(x=mitoRatio)) + 
#  	geom_density(alpha = 0.2) + 
#  	scale_x_log10() + 
#  	theme_classic() +
#  	geom_vline(xintercept = 0.2)

########### Filter to remove low quality cells ############
pdf("Multiome.rna_BF.pdf", height = 6, width = 9)

VlnPlot(Multiome.rna, features = c("nGene", "nUMI", "percent.mt"), ncol = 3)
dev.off()

#Multiome.rna <- subset(Multiome.rna, subset = nGene > 499 & nGene < 6001 & nUMI > 999 & percent.mt < 15)

pdf("Multiome.rna_QC_AF.pdf", height = 6, width = 9)
VlnPlot(Multiome.rna, features = c("nGene", "nUMI", "percent.mt"), ncol = 3)
dev.off()

########### Normalize and Scale data using sctransform #################

library(sctransform)
#Multiome.rna <- SCTransform(Multiome.rna, vars.to.regress = "percent.mt", verbose = TRUE)
Multiome.rna <- SCTransform(Multiome.rna, verbose = TRUE)

## Run PCA, autochoose dimensions
Multiome.rna <- RunPCA(Multiome.rna)

## Choose principal components to consider via two metrics:
## The point where the principal components only contribute 5% of standard deviation and
## the principal components cumulatively contribute 90% of the standard deviation.
pct <- Multiome.rna[["pca"]]@stdev / sum(Multiome.rna[["pca"]]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]

## The point where the percent change in variation between the consecutive PCs is less than 0.1%.
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# Minimum of the two calculations
pcs <- min(co1,co2)

######### Generate the Elbow plot to select number of dimensions #############
pdf("Multiome.rna_ElbowPlot.pdf", height = 6, width = 6)
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))

# Elbow plot to visualize 
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw() +
  ylab("Percent of std dev explained by PC") +
  xlab("Cumulative std dev")
dev.off()

#write.table(Multiome.rna@meta.data, file="Multiome_RNA_MetaData.txt", sep="\t", quote=FALSE)

######### Run UMAP and clustering algorithm, Select dimensions based on Elbow plot ############
Multiome.rna <- RunUMAP(Multiome.rna, dims = 1:pcs)
Multiome.rna <- FindNeighbors(Multiome.rna, dims = 1:pcs)
Multiome.rna <- FindClusters(Multiome.rna, resolution=0.3)

pdf("Multiome.rna_UMAP_Clusters.pdf", height = 6, width = 7)
DimPlot(Multiome.rna, reduction="umap")
dev.off()
#pdf("Multiome.rna_UMAP_Samples.pdf", height = 6, width = 7)
#DimPlot(Multiome.rna, reduction="umap", group.by = "Timepoints", cols=c("pink", "goldenrod3", "limegreen", "deepskyblue2"))
#dev.off()

FeaturePlot(Multiome.rna, features = "MKI67")
s_genes = readLines("g1-s-cycling")
g2m_genes = readLines("g2-m-cycling")
FeaturePlot(Multiome.rna, features = s_genes)
FeaturePlot(Multiome.rna, features = g2m_genes)
FeaturePlot(Multiome.rna, features = "percent.mt")

#dialout <- Read10X(data.dir = "/Users/neva/Documents/calico/all-p16-p14-Broad-dialout")
#rn <- rownames(dialout)
#mytable <- dialout[which(rn == "p16"):which(rn=="p14"),]
#write.table(t(as.matrix(mytable)), "all-p16-p14-Broad-dialout.tsv")
### EDIT THE HEADER TO BE cells p16.dialout p14.dialout 
### or ANRIL.short.dialout ANRIL.long.dialout
mytable2 <- read.table("all-p16-p14-Broad-dialout.tsv", head=TRUE)
metadata <- left_join(Multiome.rna@meta.data, mytable2, by = "cells")
Multiome.rna@meta.data$p16.dialout <- metadata$p16.dialout
Multiome.rna@meta.data$p14.dialout <- metadata$p14.dialout

### senesenscene scoring
rn<-read.table("rn")
Multiome.rna<-AddModuleScore(Multiome.rna, features=rn,name = "p16NormHayTop")
FeaturePlot(Multiome.rna, features="p16NormHayTop1")+ scale_color_viridis(option = "magma")
Multiome.rna$senscore.class[which(Multiome.rna$p16NormHayLeast1>=0.10726883)]<-"least10"
Multiome.rna$senscore.class[which(Multiome.rna$p16NormHayTop1>=0.5295226)]<-"top10"
df<-data.frame(Multiome.rna@meta.data$senscore.class, Multiome.rna@meta.data$cnarchr)
DimPlot(Multiome.rna,group.by = "senscore.class")
write.table(df, "ordering",row.names=FALSE, col.names=FALSE)
p16.gene<-read.table("archr-col-sen2")
p16.gene<-unlist(p16.gene)
archr$p16.gene<-p16.gene
getGroupBW(archr, groupBy = "p16.gene")






cluster0.markers <- FindMarkers(Multiome.rna, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(Multiome.rna, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(Multiome.rna, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(Multiome.rna, ident.1 = 3, min.pct = 0.25)
cluster4.markers <- FindMarkers(Multiome.rna, ident.1 = 4, min.pct = 0.25)
cluster5.markers <- FindMarkers(Multiome.rna, ident.1 = 5, min.pct = 0.25)

VlnPlot(Multiome.rna, features = c("p16.dialout", "p14.dialout", "ANRIL.long.dialout","ANRIL.short.dialout"))
cl0 <- which(Multiome.rna@meta.data$seurat_clusters == 0)
cl1 <- which(Multiome.rna@meta.data$seurat_clusters == 1)
cl2 <- which(Multiome.rna@meta.data$seurat_clusters == 2)
cl3 <- which(Multiome.rna@meta.data$seurat_clusters == 3)
cl4 <- which(Multiome.rna@meta.data$seurat_clusters == 4)
cl5 <- which(Multiome.rna@meta.data$seurat_clusters == 5)
cl0p14avg<-mean(Multiome.rna@meta.data$p14.dialout[cl0])
cl1p14avg<-mean(Multiome.rna@meta.data$p14.dialout[cl1])
cl2p14avg<-mean(Multiome.rna@meta.data$p14.dialout[cl2])
cl3p14avg<-mean(Multiome.rna@meta.data$p14.dialout[cl3])
cl4p14avg<-mean(Multiome.rna@meta.data$p14.dialout[cl4])
cl5p14avg<-mean(Multiome.rna@meta.data$p14.dialout[cl5])
cl0p16avg<-mean(Multiome.rna@meta.data$p16.dialout[cl0])
cl1p16avg<-mean(Multiome.rna@meta.data$p16.dialout[cl1])
cl2p16avg<-mean(Multiome.rna@meta.data$p16.dialout[cl2])
cl3p16avg<-mean(Multiome.rna@meta.data$p16.dialout[cl3])
cl4p16avg<-mean(Multiome.rna@meta.data$p16.dialout[cl4])
cl5p16avg<-mean(Multiome.rna@meta.data$p16.dialout[cl5])
p14.avg.by.cluster <- c(cl0p14avg,cl1p14avg, cl2p14avg, cl3p14avg, cl4p14avg, cl5p14avg)
p16.avg.by.cluster <- c(cl0p16avg,cl1p16avg, cl2p16avg, cl3p16avg, cl4p16avg, cl5p16avg)
id <- c(0,1,2,3,4,5)
df <- data.frame(id, p14.avg.by.cluster, p16.avg.by.cluster)
ggplot(df, aes(x=p16.avg.by.cluster, y=p14.avg.by.cluster)) + geom_point() + geom_text(label=id,hjust=0, vjust=2)
umis.avg.by.cluster <- c(mean(Multiome.rna@meta.data$nUMI[cl0]),mean(Multiome.rna@meta.data$nUMI[cl1]), mean(Multiome.rna@meta.data$nUMI[cl2]), mean(Multiome.rna@meta.data$nUMI[cl3]), mean(Multiome.rna@meta.data$nUMI[cl4]), mean(Multiome.rna@meta.data$nUMI[cl5]))
df <- data.frame(id, p14.avg.by.cluster, p16.avg.by.cluster, umis.avg.by.cluster)
ggplot(df, aes(x=p16.avg.by.cluster, y=umis.avg.by.cluster)) + geom_point() + geom_text(label=id,hjust=0, vjust=2)
ggplot(df, aes(x=p14.avg.by.cluster, y=umis.avg.by.cluster)) + geom_point() + geom_text(label=id,hjust=0, vjust=2)

numis <- Multiome.rna@meta.data$nUMI
p16 <- Multiome.rna@meta.data$p16.dialout
p14 <- Multiome.rna@meta.data$p14.dialout
df3 <- data.frame(numis, p16,p14)
ggplot(df3,aes(x=p14, y=numis))+geom_jitter(alpha=0.5)
cor(p14,numis)
cor(p16,numis)
ggplot(df3,aes(x=p14, y=numis))+geom_jitter(alpha=0.3)
ggplot(df3,aes(x=p16, y=numis))+geom_jitter(alpha=0.3) + geom_smooth(method=lm)
ggplot(df3,aes(x=p16, y=numis))+geom_jitter(alpha=0.3)
ggplot(df3,aes(x=p16, y=numis))+geom_jitter(alpha=0.5) + geom_label(label="corr=0.312631", x=15,y=20000)
ggplot(df3,aes(x=p14, y=numis))+geom_jitter(alpha=0.5) + geom_label(label="corr=0.463566", x=20,y=20000)


dialout.pdl50<-Read10X("hayflick_paper_dialout_all/Calico_PDL_50Solo.out/Gene/raw/")
dialout.pdl50<-Read10X("hayflick_paper_dialout_all/Calico_PDL_50Solo.out/Gene/raw/")
dialout.pdl25<-Read10X("hayflick_paper_dialout_all/Calico_PDL_25Solo.out/Gene/raw/")
dialout.pdl29<-Read10X("hayflick_paper_dialout_all/Calico_PDL_29Solo.out/Gene/raw/")
dialout.pdl33<-Read10X("hayflick_paper_dialout_all/Calico_PDL_33Solo.out/Gene/raw/")
dialout.pdl37<-Read10X("hayflick_paper_dialout_all/Calico_PDL_37Solo.out/Gene/raw/")
dialout.pdl46<-Read10X("hayflick_paper_dialout_all/Calico_PDL_46Solo.out/Gene/raw/")
rn <- rownames(dialout.pdl50)

mytable <- dialout.pdl50[which(rn == "p16"):which(rn=="p14"),]
mytable<-t(mytable)
mytable<-as.data.frame(mytable)
mytable$cells<-rownames(mytable)
metadata <- left_join(Hayflick@meta.data, mytable, by = "cells")
Hayflick@meta.data$p16.dialout.pdl50 <- metadata$p16
Hayflick@meta.data$p14.dialout.pdl50 <- metadata$p14

mytable <- dialout.pdl25[which(rn == "p16"):which(rn=="p14"),]
mytable<-t(mytable)
mytable<-as.data.frame(mytable)
mytable$cells<-rownames(mytable)
metadata <- left_join(Hayflick@meta.data, mytable, by = "cells")
Hayflick@meta.data$p16.dialout.pdl25 <- metadata$p16
Hayflick@meta.data$p14.dialout.pdl25 <- metadata$p14

mytable <- dialout.pdl29[which(rn == "p16"):which(rn=="p14"),]
mytable<-t(mytable)
mytable<-as.data.frame(mytable)
mytable$cells<-rownames(mytable)
metadata <- left_join(Hayflick@meta.data, mytable, by = "cells")
Hayflick@meta.data$p16.dialout.pdl29 <- metadata$p16
Hayflick@meta.data$p14.dialout.pdl29 <- metadata$p14

mytable <- dialout.pdl33[which(rn == "p16"):which(rn=="p14"),]
mytable<-t(mytable)
mytable<-as.data.frame(mytable)
mytable$cells<-rownames(mytable)
metadata <- left_join(Hayflick@meta.data, mytable, by = "cells")
Hayflick@meta.data$p16.dialout.pdl33 <- metadata$p16
Hayflick@meta.data$p14.dialout.pdl33 <- metadata$p14

mytable <- dialout.pdl37[which(rn == "p16"):which(rn=="p14"),]
mytable<-t(mytable)
mytable<-as.data.frame(mytable)
mytable$cells<-rownames(mytable)
metadata <- left_join(Hayflick@meta.data, mytable, by = "cells")
Hayflick@meta.data$p16.dialout.pdl37 <- metadata$p16
Hayflick@meta.data$p14.dialout.pdl37 <- metadata$p14

mytable <- dialout.pdl46[which(rn == "p16"):which(rn=="p14"),]
mytable<-t(mytable)
mytable<-as.data.frame(mytable)
mytable$cells<-rownames(mytable)
metadata <- left_join(Hayflick@meta.data, mytable, by = "cells")
Hayflick@meta.data$p16.dialout.pdl46 <- metadata$p16
Hayflick@meta.data$p14.dialout.pdl46 <- metadata$p14

inds<-c(inds1,inds2,inds3,inds4,inds5,inds6)
inds1<-which(Hayflick@meta.data$PDL == "PDL_25")
inds2<-which(Hayflick@meta.data$PDL == "PDL_29")
inds3<-which(Hayflick@meta.data$PDL == "PDL_33")
inds4<-which(Hayflick@meta.data$PDL == "PDL_37")
inds5<-which(Hayflick@meta.data$PDL == "PDL_46")
inds6<-which(Hayflick@meta.data$PDL == "PDL_50")
Hayflick@meta.data$p14 <- 0
Hayflick@meta.data$p14[inds1] <- Hayflick@meta.data$p14.dialout.pdl25[inds1]
Hayflick@meta.data$p14[inds2] <- Hayflick@meta.data$p14.dialout.pdl29[inds2]
Hayflick@meta.data$p14[inds3] <- Hayflick@meta.data$p14.dialout.pdl33[inds3]
Hayflick@meta.data$p14[inds4] <- Hayflick@meta.data$p14.dialout.pdl37[inds4]
Hayflick@meta.data$p14[inds5] <- Hayflick@meta.data$p14.dialout.pdl46[inds5]
Hayflick@meta.data$p14[inds6] <- Hayflick@meta.data$p14.dialout.pdl50[inds6]
Hayflick@meta.data$p16 <- 0
Hayflick@meta.data$p16[inds1] <- Hayflick@meta.data$p16.dialout.pdl25[inds1]
Hayflick@meta.data$p16[inds2] <- Hayflick@meta.data$p16.dialout.pdl29[inds2]
Hayflick@meta.data$p16[inds3] <- Hayflick@meta.data$p16.dialout.pdl33[inds3]
Hayflick@meta.data$p16[inds4] <- Hayflick@meta.data$p16.dialout.pdl37[inds4]
Hayflick@meta.data$p16[inds5] <- Hayflick@meta.data$p16.dialout.pdl46[inds5]
Hayflick@meta.data$p16[inds6] <- Hayflick@meta.data$p16.dialout.pdl50[inds6]
Hayflick@meta.data$p14.norm <- 0

numis=Hayflick@meta.data$nCount_RNA[inds1]
Hayflick@meta.data$p14.norm[inds1] <- Hayflick@meta.data$p14.dialout.pdl25[inds1]/numis
numis=Hayflick@meta.data$nCount_RNA[inds2]
Hayflick@meta.data$p14.norm[inds2] <- Hayflick@meta.data$p14.dialout.pdl29[inds2]/numis
numis=Hayflick@meta.data$nCount_RNA[inds3]
Hayflick@meta.data$p14.norm[inds3] <- Hayflick@meta.data$p14.dialout.pdl33[inds3]/numis
numis=Hayflick@meta.data$nCount_RNA[inds4]
Hayflick@meta.data$p14.norm[inds4] <- Hayflick@meta.data$p14.dialout.pdl37[inds4]/numis
numis=Hayflick@meta.data$nCount_RNA[inds5]
Hayflick@meta.data$p14.norm[inds5] <- Hayflick@meta.data$p14.dialout.pdl46[inds5]/numis
numis=Hayflick@meta.data$nCount_RNA[inds6]
Hayflick@meta.data$p14.norm[inds6] <- Hayflick@meta.data$p14.dialout.pdl50[inds6]/numis


atac.multiome<-read.table("cellranger_atac/atac_fragments_p16.tsv")

tmptable<-Multiome.rna@meta.data
tmptable$V4<-Multiome.rna$cells
atac.multiome$cells <- atac.multiome$V4
tmptable<-left_join(tmptable, atac.multiome, by="cells")
my_tab_sort3 <- tmptable %>%                 # Order table with dplyr
as.data.frame() %>%
arrange(desc(p16.norm))
#my_tab_sort4$pos1 <- my_tab_sort4$V2.x-6400000
#my_tab_sort4$pos1 <- my_tab_sort4$pos1/1000
#my_tab_sort4$pos2 <- my_tab_sort4$V3-6400000
#my_tab_sort4$pos2 <- my_tab_sort4$pos2/1000
rotate <- function(x) t(apply(x, 2, rev))

my_tab_sort3$pos1 <- my_tab_sort3$V2-21956130
my_tab_sort3$pos1 <- my_tab_sort3$pos1/1000
my_tab_sort3$pos2 <- my_tab_sort3$V3-21956130
my_tab_sort3$pos2 <- my_tab_sort3$pos2/1000
ncols<-214
#ncols<-280
# chr9: 21.94-22.22
newdat <- data.frame(matrix(ncol = ncols, nrow = 0))
matrixind<-1
#vals<-rep(0,400)
#fvals<-rep(0,800)
while (nrow(newdat) <= 400 && matrixind <= nrow(my_tab_sort3)) {
  cellname<-my_tab_sort3$cell[matrixind]
  new <- rep(0,ncols)
  while (my_tab_sort3$cell[matrixind] == cellname && matrixind <= nrow(my_tab_sort3)) {
    if (my_tab_sort3$pos1[matrixind]>0 & my_tab_sort3$pos1[matrixind]<=ncols) {
      #new[my_tab_sort3$pos1[matrixind]] <-new[my_tab_sort3$pos1[matrixind]]+1
      new[my_tab_sort3$pos1[matrixind]] <-1
    }
    if (my_tab_sort3$pos2[matrixind]>0 & my_tab_sort3$pos2[matrixind]<=ncols) {
      #new[my_tab_sort3$pos2[matrixind]] <-new[my_tab_sort3$pos2[matrixind]]+1
      new[my_tab_sort3$pos2[matrixind]] <-1
    }
    matrixind<-matrixind+1
  }
  newdat[nrow(newdat)+1,] <- new;
}
new <- rep(0,ncols)
print(matrixind)

#x5 <- sample(7200:, 1000, replace=F)
#ind.to.ind<-1
matrixind <-1
while (nrow(newdat) <= 286) {
  cellname<-my_tab_sort3$cell[matrixind]
  new <- rep(0,ncols)
  while (my_tab_sort3$cell[matrixind] == cellname) {
    if (my_tab_sort3$pos1[matrixind]>0 & my_tab_sort3$pos1[matrixind]<=ncols) {
      #new[my_tab_sort3$pos1[matrixind]] <-new[my_tab_sort3$pos1[matrixind]]+1
      new[my_tab_sort3$pos1[matrixind]] <-1
    }
    if (my_tab_sort3$pos2[matrixind]>0 & my_tab_sort3$pos2[matrixind]<=ncols) {
      #new[my_tab_sort3$pos2[matrixind]] <-new[my_tab_sort3$pos2[matrixind]]+1
      new[my_tab_sort3$pos2[matrixind]] <- 1
    }
    matrixind<-matrixind+1
  }
  newdat[nrow(newdat)+1,] <- new;
}
image(rotate(as.matrix(newdat)), axes=FALSE, col=c("white", "black"))
axis(1, at=seq(0,1,0.5), labels=c("chr12:6.4","GADPH","chr12:6.6"))
y.lab<-c(0,0,round(vals[301],4),round(vals[201],4),round(vals[101],4), round(vals[1],4))
y.lab<-c("p14","p16")
axis(2, at=seq(0.5,1,0.5), labels=y.lab)

ggplot(df3, mapping = aes(x=Multiome.rna.p16.norm, y=Multiome.rna.p14.norm, color= I(ifelse(Multiome.rna.p14.norm<=0.0025 &Multiome.rna.p16.norm <= 0.0025, 'black', 'red'))))+geom_point()+xlim(0,0.025)+ylim(0,0.025)+geom_hline(yintercept=0.0025, linetype="dashed", color = "blue", size=0.5)+geom_vline(xintercept=0.0025, linetype="dashed", color = "blue", size=0.5)+xlab("p16 dialout umis / total umis") + ylab("p14 dialout umis / total umis") + ggtitle("Wi38 confluence p14 and p16 per single cell")
ggplot(df2, mapping = aes(x=Hayflick.p16.norm, y=Hayflick.p14.norm, color=I(ifelse(Hayflick.p16.norm<=0.0005 & Hayflick.p14.norm<=0.0005, 'black', 'red'))))+geom_point()+ylim(0,0.0035)+geom_hline(yintercept=0.0005, linetype="dashed", color = "blue", size=0.5)+geom_vline(xintercept=0.0005, linetype="dashed", color = "blue", size=0.5)+xlab("p16 dialout umis / total umis") + ylab("p14 dialout umis / total umis") + ggtitle("Wi38 senescence all PDLs per single cell")

# possibly didn't need to make this separate object
obj<-Hayflick
Idents(object=obj) <- "sample"

sub_obj <- subset(x = obj, idents = c(12,16), invert = TRUE)

p1 <- FeaturePlot(sub_obj, features = c("p16.norm", "p14.norm"), combine = FALSE )
fix.sc <- scale_color_gradientn( colours = c('white', 'red'),  limits = c(0.00001, 0.0015))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)


########## Perform DGE and ROC analyses #############
#DefaultAssay(Multiome.rna) <- "RNA"

#Reprogramming.rna <- NormalizeData(Reprogramming.rna)
#all.genes <- rownames(Reprogramming.rna)
#Reprogramming.rna <- ScaleData(Reprogramming.rna, vars.to.regress = "percent.mt", features = all.genes)
#Reprogramming.rna.markers <- FindAllMarkers(Reprogramming.rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc")
#write.table(Reprogramming.rna.markers, file="Reprogramming.rna_Markers_AUC.txt", sep="\t", quote=FALSE, col.names=NA)
#Reprogramming.rna.markers.DGE <- FindAllMarkers(Reprogramming.rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.table(Reprogramming.rna.markers.DGE, file="Reprogramming.rna_Markers_DGE.txt", sep="\t", quote=FALSE, col.names=NA)
#write.table(Reprogramming.rna@meta.data, file="Reprogramming.rna_MetaData.txt", sep="\t", col.names=NA, quote=FALSE)


# p16.umi.per.million.0 <-Multiome.rna@meta.data$p16.dialout[cl0]*1000000/Multiome.rna@meta.data$nUMI[cl0]
# p16.umi.per.million.1 <-Multiome.rna@meta.data$p16.dialout[cl1]*1000000/Multiome.rna@meta.data$nUMI[cl1]
# p16.umi.per.million.2 <-Multiome.rna@meta.data$p16.dialout[cl2]*1000000/Multiome.rna@meta.data$nUMI[cl2]
# p16.umi.per.million.3 <-Multiome.rna@meta.data$p16.dialout[cl3]*1000000/Multiome.rna@meta.data$nUMI[cl3]
# p16.umi.per.million.4 <-Multiome.rna@meta.data$p16.dialout[cl4]*1000000/Multiome.rna@meta.data$nUMI[cl4]
# p16.umi.per.million.5 <-Multiome.rna@meta.data$p16.dialout[cl5]*1000000/Multiome.rna@meta.data$nUMI[cl5]
# p16.umi.per.million <- c(mean(p16.umi.per.million.0),mean(p16.umi.per.million.1),mean(p16.umi.per.million.2),mean(p16.umi.per.million.3),mean(p16.umi.per.million.4),mean(p16.umi.per.million.5))
# p14.umi.per.million.0 <-Multiome.rna@meta.data$p14.dialout[cl0]*1000000/Multiome.rna@meta.data$nUMI[cl0]
# p14.umi.per.million.1 <-Multiome.rna@meta.data$p14.dialout[cl1]*1000000/Multiome.rna@meta.data$nUMI[cl1]
# p14.umi.per.million.2 <-Multiome.rna@meta.data$p14.dialout[cl2]*1000000/Multiome.rna@meta.data$nUMI[cl2]
# p14.umi.per.million.3 <-Multiome.rna@meta.data$p14.dialout[cl3]*1000000/Multiome.rna@meta.data$nUMI[cl3]
# p14.umi.per.million.4 <-Multiome.rna@meta.data$p14.dialout[cl4]*1000000/Multiome.rna@meta.data$nUMI[cl4]
# p14.umi.per.million.5 <-Multiome.rna@meta.data$p14.dialout[cl5]*1000000/Multiome.rna@meta.data$nUMI[cl5]
# p14.umi.per.million <- c(mean(p14.umi.per.million.0),mean(p14.umi.per.million.1),mean(p14.umi.per.million.2),mean(p14.umi.per.million.3),mean(p14.umi.per.million.4),mean(p14.umi.per.million.5))

# p14 0.00015 .0003  0.0006 0.00095
# p16 0.0001 0.0005
# i0 <- (Multiome.rna$p14.norm == 0 )
# i1 <- (Multiome.rna$p14.norm <=0.00015 & Multiome.rna$p14.norm > 0)
# i2 <- (Multiome.rna$p14.norm >0.00015 & Multiome.rna$p14.norm <0.0003)
# i3 <- (Multiome.rna$p14.norm >=0.0003 & Multiome.rna$p14.norm <.0006)
# i4 <- (Multiome.rna$p14.norm >=0.0006 & Multiome.rna$p14.norm <0.00095)
# i5 <- (Multiome.rna$p14.norm >=0.00095)
# 
# j0 <- (Multiome.rna$p16.norm == 0 )
# j1 <- (Multiome.rna$p16.norm <=0.0001 & Multiome.rna$p16.norm>0)
# j2 <- (Multiome.rna$p16.norm >0.0001 & Multiome.rna$p16.norm <0.0005)
# j3 <- (Multiome.rna$p16.norm >=0.0005 )
# 
# row1 <- c(sum(i0&j0),sum(i0&j1),sum(i0&j2), sum(i0&j3))
# row2 <- c(sum(i1&j0),sum(i1&j1),sum(i1&j2), sum(i1&j3))
# row3 <- c(sum(i2&j0),sum(i2&j1),sum(i2&j2), sum(i2&j3))
# row4 <- c(sum(i3&j0),sum(i3&j1),sum(i3&j2), sum(i3&j3))
# row5 <- c(sum(i4&j0),sum(i4&j1),sum(i4&j2), sum(i4&j3))
# row6 <- c(sum(i5&j0),sum(i5&j1),sum(i5&j2), sum(i5&j3))

pdls <- c(rep("PDL 25" , 2) , rep("PDL 29" , 2) , rep("PDL 33" , 2) , rep("PDL 37" , 2), rep("PDL 46", 2), rep("PDL 50", 2) )
expression <- rep(c("p14 dialout" , "p16 dialout") , 6)
values<-c(sum(Hayflick$p14.dialout.pdl25[inds1]/Hayflick$nCount_RNA[inds1]),sum(Hayflick$p16.dialout.pdl25[inds1]/Hayflick$nCount_RNA[inds1]), sum(Hayflick$p14.dialout.pdl29[inds2]/Hayflick$nCount_RNA[inds2]), sum(Hayflick$p16.dialout.pdl29[inds2]/Hayflick$nCount_RNA[inds2]), sum(Hayflick$p14.dialout.pdl33[inds3]/Hayflick$nCount_RNA[inds3]), sum(Hayflick$p16.dialout.pdl33[inds3]/Hayflick$nCount_RNA[inds3]), sum(Hayflick$p14.dialout.pdl37[inds4]/Hayflick$nCount_RNA[inds4]), sum(Hayflick$p16.dialout.pdl37[inds4]/Hayflick$nCount_RNA[inds4]), sum(Hayflick$p14.dialout.pdl46[inds5]/Hayflick$nCount_RNA[inds5]), sum(Hayflick$p16.dialout.pdl46[inds5]/Hayflick$nCount_RNA[inds5]), sum(Hayflick$p14.dialout.pdl50[inds6]/Hayflick$nCount_RNA[inds6]), sum(Hayflick$p16.dialout.pdl50[inds6]/Hayflick$nCount_RNA[inds6]))
databr <- data.frame(pdls,expression,values)
ggplot(databr, aes(fill=expression, y=values, x=pdls)) +  geom_bar(position="dodge", stat="identity")+ xlab("Time points")+ylab("Sum Normalized UMIs")+scale_fill_manual(values = c("#054C70","#05C3DE")) +theme(legend.position="top",legend.title =
element_blank())


quants<-quantile(Multiome.rna$p16NormHayTop1, c(.25, .50, .75))
Multiome.rna$senscore.16.class<-"none"
Multiome.rna$senscore.16.class[which(Multiome.rna$p16NormHayTop1<=quants[1])]<-"least25"
Multiome.rna$senscore.16.class[which(Multiome.rna$p16NormHayTop1>=quants[3])]<-"top25"
quants<-quantile(Multiome.rna$p14NormHayTop1, c(.25, .50, .75))
Multiome.rna$senscore.14.class<-"none"
Multiome.rna$senscore.14.class[which(Multiome.rna$p14NormHayTop1<=quants[1])]<-"least25"
Multiome.rna$senscore.14.class[which(Multiome.rna$p14NormHayTop1>=quants[3])]<-"top25"
DimPlot(Multiome.rna,group.by = "senscore.14.class")
df<-data.frame(Multiome.rna@meta.data$senscore.14.class, Multiome.rna@meta.data$cnarchr)
write.table(df, "ordering14", row.names=FALSE, col.names=FALSE)
# edit via awk code: 
# awk 'FNR==NR{a[$2]=$1}FNR!=NR && FNR>1{if ($2 in a){print a[$2]}else{print "none"}}' ordering14 archr-order-cellnames > archr-order-p14-class
p14.class<-read.table("archr-order-p14-class")
p14.class<-unlist(p14.class)
archr$p14.class<-p14.class
getGroupBW(archr, groupBy = "p14.class")
p16.class<-read.table("archr-order-p16-class")
p16.class<-unlist(p16.class)
archr$p16.class<-p16.class
getGroupBW(archr, groupBy = "p16.class")
df<-as.data.frame(t(assay(p16pks)))
mycorrs<-numeric(44)
mycol<-colData(p16pks)$p16.gene
for(i in 1:ncol(df)) {
c1<-df[ , i]
mycorrs[i]<-cor(c1,mycol, method="Spearman")
}
mycorrs<-numeric(44)
mycol<-colData(p16pks)$p16.gene
for(i in 1:ncol(df)) {
c1<-df[ , i]
mycorrs[i]<-cor(c1,mycol, method="spearman")
}
max(mycorrs)
mypvals<-numeric(44)
for(i in 1:ncol(df)) {
c1<-df[ , i]
res<-cor.test(c1,mycol, method="spearman")
pvalue <- res$p.value
mypvals[i]<-pvalue
}
#p16.gene<-read.table("archr-order-p16")
df<-data.frame(Multiome.rna@meta.data$p14NormHayTop1, Multiome.rna@meta.data$cnarchr)
write.table(df, "ordering2",row.names=FALSE, col.names=FALSE)
p14.gene<-read.table("archr-order-p14")
p14.gene<-unlist(p14.gene)
archr$p14.gene<-p14.gene
mycol<-colData(p16pks)$p14.gene
mypvals2<-numeric(44)
for(i in 1:ncol(df)) {
c1<-df[ , i]
mycorrs2[i]<-cor(c1,mycol)
}
df<-as.data.frame(t(assay(p16pks)))
for(i in 1:ncol(df)) {
c1<-df[ , i]
mycorrs2[i]<-cor(c1,mycol)
}
mtx <- getMatrixFromProject(archr, useMatrix = "PeakMatrix")
p16pks<-subsetByOverlaps(mtx, roi)
df<-as.data.frame(t(assay(p16pks)))
for(i in 1:ncol(df)) {
c1<-df[ , i]
mycorrs2[i]<-cor(c1,mycol)
}
View(df)
View(mycol)
mycol<-colData(p16pks)$p14.gene
View(mycol)
for(i in 1:ncol(df)) {
c1<-df[ , i]
mycorrs2[i]<-cor(c1,mycol)
}
for(i in 1:ncol(df)) {
c1<-df[ , i]
res<-cor.test(c1,mycol)
pvalue <- res$p.value
mypvals2[i]<-pvalue
}
plot(mycorrs2, mypvals2)
which(mypvals2<0.01)
mycorrs2[which(mypvals2<0.01)]
mycorrs[25]
mycorrs2[25]
mypvals[25]
mypvals2[25]
cor.test(mycol, df[,25])
mypvals[26]
mycorrs2[which(mypvals2<0.01)]
which(mypvals2<0.01)
plot(Multiome.rna$p14NormHayTop1, Multiome.rna$p16NormHayTop1)

