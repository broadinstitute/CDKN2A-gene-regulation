setwd("/Volumes/etorlait/p16/Salmon_Remapping/")

samples <- read.table("samples.txt", header = TRUE)
samples

samples=samples[1:7,]

files <- file.path("/Volumes/etorlait/p16/Salmon_Remapping/",samples$dir, "/quant.sf")
all(file.exists(files))

names(files) <- paste0(samples$sample)

mart_gene_transc <- read.delim('/Volumes/etorlait/p16/Tyler_RS_Salmon_Quants_fromNeva/ENSEMBLT_to_Gene/mart_export_downloaded_20230303.txt', sep='\t')

gene_gs=subset(mart_gene_transc, select=c('Gene.stable.ID', 'Gene.name'))

#Ensembl library for hg38
library("EnsDb.Hsapiens.v86")

edb <- EnsDb.Hsapiens.v86

transcripts <- transcripts(edb, columns=c("tx_id", "gene_id"), return.type="DataFrame")

head(transcripts)
     
library(tximport)

txi <- tximport(files, type = "salmon", tx2gene = transcripts, txIn=TRUE, ignoreTxVersion = TRUE)

tpms=txi$abundance

colSums(tpms)
#note that the values now do not sum to 1M but at 963,700 - 966,784

#assign Gene symbols to each EnsemblGene ID
tpms2=unique(merge(tpms, gene_gs, by.x="row.names", by.y='Gene.stable.ID'))

cell_cycle=subset(tpms2, Gene.name %in% c('MKI67','TOP2A'))

write.table(cell_cycle, 'replSen_cellcycleGenes_fromTtoG_tximport.txt', sep = '\t', quote = F, row.names = F, col.names = T)

#SenMayo gene list from publication
senmayo = c('ACVR1B', 'ANG', 'ANGPT1', 'ANGPTL4', 'AREG', 'AXL', 'BEX3', 'BMP2', 'BMP6', 'C3', 'CCL1', 'CCL13', 'CCL16', 'CCL2', 'CCL20', 'CCL24', 'CCL26', 'CCL3', 'CCL3L1', 'CCL4', 'CCL5', 'CCL7', 'CCL8', 'CD55', 'CD9', 'CSF1', 'CSF2', 'CSF2RB', 'CST4', 'CTNNB1', 'CTSB', 'CXCL1', 'CXCL10', 'CXCL12', 'CXCL16', 'CXCL2', 'CXCL3', 'CXCL8', 'CXCR2', 'DKK1', 'EDN1', 'EGF', 'EGFR', 'EREG', 'ESM1', 'ETS2', 'FAS', 'FGF1', 'FGF2', 'FGF7', 'GDF15', 'GEM', 'GMFG', 'HGF', 'HMGB1', 'ICAM1', 'ICAM3', 'IGF1', 'IGFBP1', 'IGFBP2', 'IGFBP3', 'IGFBP4', 'IGFBP5', 'IGFBP6', 'IGFBP7', 'IL10', 'IL13', 'IL15', 'IL18', 'IL1A', 'IL1B', 'IL2', 'IL32', 'IL6', 'IL6ST', 'IL7', 'INHA', 'IQGAP2', 'ITGA2', 'ITPKA', 'JUN', 'KITLG', 'LCP1', 'MIF', 'MMP1', 'MMP10', 'MMP12', 'MMP13', 'MMP14', 'MMP2', 'MMP3', 'MMP9', 'NAP1L4', 'NRG1', 'PAPPA', 'PECAM1', 'PGF', 'PIGF', 'PLAT', 'PLAU', 'PLAUR', 'PTBP1', 'PTGER2', 'PTGES', 'RPS6KA5', 'SCAMP4', 'SELPLG', 'SEMA3F', 'SERPINB4', 'SERPINE1', 'SERPINE2', 'SPP1', 'SPX', 'TIMP2', 'TNF', 'TNFRSF10C', 'TNFRSF11B', 'TNFRSF1A', 'TNFRSF1B', 'TUBGCP2', 'VEGFA', 'VEGFC', 'VGF', 'WNT16', 'WNT2')

senmayo_genes=subset(tpms2, Gene.name %in% senmayo)

senmayo_scores=colSums(senmayo_genes[,2:8])

#PDL21     PDL33     PDL41     PDL48    PDL535     PDL54     PDL57 
#8817.053 11335.822 10164.158 13737.148 11607.014 12664.602 12380.462 

write.table(senmayo_scores, 'replSen_SenMayoScores_fromTtoG_tximport.txt', sep = '\t', quote = F, row.names = F, col.names = T)



