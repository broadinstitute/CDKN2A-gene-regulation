setwd('/Volumes/etorlait/p16/Plotgardener/')

library(plotgardener)
library(RColorBrewer)

hicDataChrom <- readHic(file = "Wi21.hic",
                        chrom = "chr9", assembly = "hg38",
                        resolution = 1000, res_scale = "BP", norm = "NONE")

hicDataChrom_latePDL <- readHic(file = "microC_Wi48_Rep2.hic",
                                chrom = "chr9", assembly = "hg38",
                                resolution = 1000, res_scale = "BP", norm = "NONE")


# Merge the data frames by chr9_A and chr9_B
merged_data <- merge(hicDataChrom, hicDataChrom_latePDL, by = c("chr9_A", "chr9_B"), all = TRUE)

# Replace NA values with 0 for subtraction
merged_data[is.na(merged_data)] <- 0

# Subtract the counts from the two data frames
merged_data$counts <- (merged_data$counts.y / sum(hicDataChrom_latePDL$counts) *1e6) - (merged_data$counts.x / sum(hicDataChrom$counts)*1e6)

# Select relevant columns and rename them
result <- merged_data[, c("chr9_A", "chr9_B", "counts")]
names(result) <- c("chr9_A", "chr9_B", "counts")


######


## pgParams
params <- pgParams(
  chrom = "chr9", chromstart = 21750000, chromend = 22500000,
  assembly = "hg38",
  x = 1, just = c("left", "top"),
  width = 8, length = 2, default.units = "inches"
)

params2 <- pgParams(
  chrom = "9", chromstart = 21750000, chromend = 22500000,
  assembly = "hg38",
  x = 1, just = c("left", "top"),
  width = 8, length = 2, default.units = "inches"
)

pdf('plotGardener_Wi38_late_min_early_chr9_21750000_22500000_with_isoforms.pdf', width =10 , height = 8.5)

## Create a plotgardener page
pageCreate(
  width = 10, height = 8.5, default.units = "inches",
  showGuides = F, xgrid = 0, ygrid = 0
)


## Plot Hi-C data in region
hicPlot<-plotHicTriangle(
  data = result,
  params= params,
  x = 1, y = 0, width = 8, height = 5,
  just = c("left", "top"), default.units = "inches", colorTrans='linear', zrange=c(-2,2), palette = colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))))


## Annotate Hi-C heatmap legend
annoHeatmapLegend(
  plot = hicPlot, fontsize = 7,
  orientation = "v",
  x = 8, y = 2.25,
  width = 0.07, height = 0.5, just = c("left", "top"),
  default.units = "inches"
)

## Add text labels
plotText(
  label = "MicroC PDL48 - PDL21, normalized by total reads in chr9, 1kb res", fontsize = 7, fontcolor = "darkred",
  x = 7, y = 1.75, just = c("left", "top"),
  default.units = "inches"
)


genePlot <- plotGenes(
  chrom = "chr9", chromstart = 21750000, chromend = 22500000,
  assembly = assembly(Genome = "hg38refGene", TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", OrgDb = "org.Hs.eg.db"), geneOrder = c('CDKN2A', 'CDKN2B', 'DMRTA1'),  
  x = 1, y = 5, width = 8, height = 1, 
  just = c("left", "top"), default.units = "inches", fill = c("black", "black"), fontcolor = c("black", "black")
)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

tbg <- transcriptsBy(txdb, by="gene")

hilite <- data.frame(transcript=c("ENST00000304494.10","ENST00000579755.2"), color=c("#1B76BB","#009545"))

plotTranscripts(
  params = params, x = 1, y = 6, width = 8, height = 1.1, fill = c("grey", "grey"),
  transcriptHighlight=hilite
)


## Plot genome label
plotGenomeLabel(
  params= params,
  x = 1, y = 7.5, length = 8, scale = "bp", at =seq(21750000, 22500000, by=50000), 
  just = c("left", "top"), default.units = "inches"
)

dev.off()


