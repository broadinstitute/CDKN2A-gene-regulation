## This was an attempt to use DESeq to do differences between Hi-C files
## There are too many false positives with this method IMO
##
library(strawr)

read_microC <- function(name) {
  data <- straw(norm = "NONE", fname = paste0("https://storage.googleapis.com/broad-p16-calico/microC/microC_", name, ".hic"), chr1loc = "chr9", chr2loc = "chr9", unit = "BP", binsize = 5000)
  data <- data[abs(data$x - data$y) >= 15000,]
  data$counts <- round(data$counts)
  rownames(data) <- paste(data$x, data$y, sep = "-")
  data[[name]]<-data$counts
  data$x<-NULL
  data$y<-NULL
  data$counts<-NULL
  return(data)
}

esc_rep1<-read_microC("ESC_Rep1")
esc_rep2<-read_microC("ESC_Rep2")
huvec_rep1<-read_microC("HUVEC_Rep1")
huvec_rep2<-read_microC("HUVEC_Rep2")
wi21_rep1<-read_microC("Wi21_Rep1") 
wi21_rep2<-read_microC("Wi21_Rep2") 
wi33_rep1<-read_microC("Wi33_Rep1") 
wi33_rep2<-read_microC("Wi33_Rep2") 
wi48_rep2<-read_microC("Wi48_Rep2") 
wi53_rep1<-read_microC("Wi53_Rep1") 
wi57_rep1<-read_microC("Wi57_Rep1")
wi57_rep2<-read_microC("Wi57_Rep2") 


df_list <- list(wi21_rep1, wi21_rep2, wi33_rep1, wi33_rep2, wi48_rep2,wi53_rep1,wi57_rep1,wi57_rep2)
# Use Reduce to merge all data frames in the list by row names
merged_df <- Reduce(function(x, y) {
  df <- merge(x, y, by = "row.names", all = TRUE)
  rownames(df) <- df$Row.names
  df$Row.names <- NULL
  df
}, df_list)
merged_df[is.na(merged_df)] <- 0

coldata<-data.frame(condition=c("cycling", "cycling", "cycling", "cycling", "senescent", "senescent", "senescent", "senescent"))
rownames(coldata)<-colnames(merged_df)
dds <- DESeqDataSetFromMatrix(countData = merged_df,
colData = coldata,
design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)
plotDispEsts(dds)

coldata<-data.frame(condition=c("wi21", "wi21", "esc", "esc"))
rownames(coldata)<-colnames(merged_df2)
dds2 <- DESeqDataSetFromMatrix(countData = merged_df2, colData = coldata,design = ~ condition)
dds2 <- DESeq(dds2)
res <- results(dds2)
summary(res)
