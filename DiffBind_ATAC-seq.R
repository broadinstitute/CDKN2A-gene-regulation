## --------------------------------------------------------------------------
## DiffBind/DESeq2 differential accessibility analysis of the Chan et al. 2022
## (GSE175533) WI-38 ATAC-seq time-course (PDL20→PDL50, 7 timepoints, 3 reps).
##
## Outputs the differentially accessible regions (DARs) overlapping the p16
## locus that were used to seed the 173-region CRISPRa screen target list
## (data/CRISPRa_target_regions.bed). See README §"CRISPRa screen (Fig 6)".
##
## Inputs expected next to this script:
##   - Wi38_calico_ATAC-seq.csv             (DiffBind sample sheet, in repo)
##   - atac_bams/<sample>.shift.sorted.bam  (bowtie2 outputs from GSE175533)
##   - atac_peaks/<sample>_peaks.xls        (MACS2 outputs from GSE175533)
##   - hg38-blacklist.v2.bed                (ENCODE blacklist, Boyle-Lab/Blacklist)
##
## ATAC-seq read processing (bowtie2 + macs2) is described in Methods; that
## step is not included in this repo.
## --------------------------------------------------------------------------

library(DiffBind)

compDERs <- function(Wi38, contrast, comp, fc) {
    der <- dba.contrast(Wi38, design = "~ Condition", contrast)
    der <- dba.analyze(der, bParallel = TRUE,
                       bBlacklist = "hg38-blacklist.v2.bed")
    der <- dba.report(der, th = 0.05, fold = log2(fc), bNormalized = TRUE,
                      bCounts = TRUE)
    mcols(der, level = "within")$Comp <- comp
    return(der)
}

Wi38.counts <- dba(sampleSheet = "Wi38_calico_ATAC-seq.csv")
Wi38.counts <- dba.count(Wi38.counts, summits = 500, minOverlap = 2,
                         bParallel = TRUE)

info <- dba.show(Wi38.counts)
libsizes <- cbind(LibReads = info$Reads, FRiP = info$FRiP,
                  PeakReads = round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

Wi38.norm <- dba.normalize(Wi38.counts, method = DBA_ALL_METHODS,
                           normalize = DBA_NORM_NATIVE, background = TRUE)
Wi38.analyze <- dba.analyze(Wi38.norm, method = DBA_ALL_METHODS,
                            design = "~ Condition", bParallel = TRUE,
                            bBlacklist = "hg38-blacklist.v2.bed")
dba.show(Wi38.analyze, bContrasts = TRUE)

## consecutive time-point comparisons
TP2vTP1 <- compDERs(Wi38.norm, c("Condition", "TP2", "TP1"), "TP2vTP1", 1.5)
TP3vTP2 <- compDERs(Wi38.norm, c("Condition", "TP3", "TP2"), "TP3vTP2", 1.5)
TP4vTP3 <- compDERs(Wi38.norm, c("Condition", "TP4", "TP3"), "TP4vTP3", 1.5)
TP5vTP4 <- compDERs(Wi38.norm, c("Condition", "TP5", "TP4"), "TP5vTP4", 1.5)
TP6vTP5 <- compDERs(Wi38.norm, c("Condition", "TP6", "TP5"), "TP6vTP5", 1.5)
TP7vTP6 <- compDERs(Wi38.norm, c("Condition", "TP7", "TP6"), "TP7vTP6", 1.5)

## TP7 vs every earlier time point
TP7vTP1 <- compDERs(Wi38.norm, c("Condition", "TP7", "TP1"), "TP7vTP1", 1.5)
TP7vTP2 <- compDERs(Wi38.norm, c("Condition", "TP7", "TP2"), "TP7vTP2", 1.5)
TP7vTP3 <- compDERs(Wi38.norm, c("Condition", "TP7", "TP3"), "TP7vTP3", 1.5)
TP7vTP4 <- compDERs(Wi38.norm, c("Condition", "TP7", "TP4"), "TP7vTP4", 1.5)
TP7vTP5 <- compDERs(Wi38.norm, c("Condition", "TP7", "TP5"), "TP7vTP5", 1.5)
TP7vTP6 <- compDERs(Wi38.norm, c("Condition", "TP7", "TP6"), "TP7vTP6", 1.5)

## subset to the p16 locus
p16 <- GRanges(seqnames = "chr9", ranges = IRanges(start = 21889107,
                                                   end = 22285292))
TP2vTP1.p16 <- TP2vTP1[queryHits(findOverlaps(TP2vTP1, p16)), "Comp"]
TP3vTP2.p16 <- TP3vTP2[queryHits(findOverlaps(TP3vTP2, p16)), "Comp"]
TP4vTP3.p16 <- TP4vTP3[queryHits(findOverlaps(TP4vTP3, p16)), "Comp"]
TP5vTP4.p16 <- TP5vTP4[queryHits(findOverlaps(TP5vTP4, p16)), "Comp"]
TP6vTP5.p16 <- TP6vTP5[queryHits(findOverlaps(TP6vTP5, p16)), "Comp"]

TP7vTP1.p16 <- TP7vTP1[queryHits(findOverlaps(TP7vTP1, p16)), "Comp"]
TP7vTP2.p16 <- TP7vTP2[queryHits(findOverlaps(TP7vTP2, p16)), "Comp"]
TP7vTP3.p16 <- TP7vTP3[queryHits(findOverlaps(TP7vTP3, p16)), "Comp"]
TP7vTP4.p16 <- TP7vTP4[queryHits(findOverlaps(TP7vTP4, p16)), "Comp"]
TP7vTP5.p16 <- TP7vTP5[queryHits(findOverlaps(TP7vTP5, p16)), "Comp"]
TP7vTP6.p16 <- TP7vTP6[queryHits(findOverlaps(TP7vTP6, p16)), "Comp"]

ders <- c(TP2vTP1.p16, TP3vTP2.p16, TP4vTP3.p16, TP5vTP4.p16, TP6vTP5.p16,
          TP7vTP1.p16, TP7vTP2.p16, TP7vTP3.p16, TP7vTP4.p16, TP7vTP5.p16,
          TP7vTP6.p16)
names(ders) <- paste0(ders$Comp, "_", names(ders))

dir.create("ATAC-seq_DiffBind", showWarnings = FALSE)
rtracklayer::export.bed(ders, con = "ATAC-seq_DiffBind/ders.bed")
rtracklayer::export.bed(ders, con = "ATAC-seq_DiffBind/ders_uniq_TP7_IGV.bed")

dba.save(DBA = Wi38.counts,  file = "Wi38_ATAC-seq_counts",  dir = "ATAC-seq_DiffBind")
dba.save(DBA = Wi38.norm,    file = "Wi38_ATAC-seq_norm",    dir = "ATAC-seq_DiffBind")
dba.save(DBA = Wi38.analyze, file = "Wi38_ATAC-seq_analyze", dir = "ATAC-seq_DiffBind")
save(TP2vTP1, TP3vTP2, TP4vTP3, TP5vTP4, TP6vTP5,
     TP7vTP1, TP7vTP2, TP7vTP3, TP7vTP4, TP7vTP5, TP7vTP6,
     file = "ATAC-seq_DiffBind/time_point_comparison.RData")
