library(plotgardener)
library(RColorBrewer)
library("TxDb.Hsapiens.UCSC.hg38.refGene")
library("org.Hs.eg.db")

# --- Helper Function Definition ---

generate_hic_plots <- function(hic_file_path, sample_label) {

  message(paste("Processing:", hic_file_path, "| Label:", sample_label))
  
  # 1. Setup Common Parameters
  # ---------------------------
  bright_reds <- colorRampPalette(brewer.pal(9, "Reds")[1:7])
  
  # Define assembly once for reuse
  assembly_hg38 <- assembly(
    Genome = "hg38refGene", 
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", 
    OrgDb = "org.Hs.eg.db"
  )
  
  # Define the two plot configurations (Zoom In vs Zoom Out)
  # List structure: Resolution, Start, End, Z-Range Max
  configs <- list(
    zoom_in  = list(res = 1000, start = 21750000, end = 22500000, zmax = 5,  desc = "1kb resolution"),
    zoom_out = list(res = 5000, start = 20300000, end = 22600000, zmax = 50, desc = "5kb resolution")
  )
  
  # 2. Iterate through configurations to generate both PDFs
  # -----------------------------------------------------
  for (cfg in configs) {
    
    # Define Params for this specific plot
    chrom <- "chr9"
    
    # Read Hi-C Data with no normalization, then normalize by total chr9 reads
    hicData <- readHic(
      file = hic_file_path,
      chrom = chrom,
      assembly = "hg38",
      resolution = cfg$res,
      res_scale = "BP",
      norm = "NONE"
    )
    hicData$counts <- hicData$counts / sum(hicData$counts) * 1e6
    
    # Define Page Parameters
    params <- pgParams(
      chrom = chrom, 
      chromstart = cfg$start, 
      chromend = cfg$end,
      assembly = "hg38",
      x = 1, 
      just = c("left", "top"),
      width = 8,
      length = 2,
      default.units = "inches"
    )
    
    # Construct output filename dynamically
    # Uses basename() to strip directory paths and save locally
    out_name <- paste0(basename(hic_file_path), "_", cfg$res, "_", chrom, "_", cfg$start, "-", cfg$end, ".pdf")
    
    pdf(out_name, width = 12, height = 12)
    
    # Create Page
    pageCreate(
      width = 10, height = 10, default.units = "inches",
      showGuides = FALSE, xgrid = 0, ygrid = 0
    )
    
    # Plot Hi-C Triangle
    hicPlot <- plotHicTriangle(
      data = hicData,
      params = params,
      x = 1, y = -0.5, width = 8, height = 5,
      just = c("left", "top"), 
      default.units = "inches", 
      colorTrans = 'linear', 
      zrange = c(0, cfg$zmax), 
      palette = bright_reds
    )
    
    # Annotate Legend
    annoHeatmapLegend(
      plot = hicPlot, fontsize = 7,
      orientation = "v",
      x = 8, y = 2.75,
      width = 0.07, height = 0.5, 
      just = c("left", "top"),
      default.units = "inches"
    )
    
    # Add Text Label
    plotText(
      label = paste(sample_label, cfg$desc), 
      fontsize = 7, fontcolor = "darkred",
      x = 7, y = 1.75, 
      just = c("left", "top"),
      default.units = "inches"
    )
    
    # Plot Genes
    plotGenes(
      chrom = chrom, chromstart = cfg$start, chromend = cfg$end,
      assembly = assembly_hg38, 
      geneOrder = c('CDKN2A', 'CDKN2B', 'DMRTA1'),  
      x = 1, y = "0b", width = 8, height = 1, 
      just = c("left", "top"), default.units = "inches"
    )
    
    # Plot Genome Label
    plotGenomeLabel(
      params = params,
      x = 1, y = 5.65, length = 8, scale = "bp", 
      at = seq(cfg$start, cfg$end, by = 50000), 
      just = c("left", "top"), default.units = "inches"
    )
    
    dev.off()
    message(paste("  -> Generated PDF:", out_name))
  }
}


# --- Execution Section ---

# File 1: Wi21 combined (Rep1 + Rep2)
generate_hic_plots(
  hic_file_path = "microC_WI-38_PDL21.hic",
  sample_label = "PDL21"
)

# File 2: Wi48 Rep2
generate_hic_plots(
  hic_file_path = "microC_WI-38_PDL48_Rep2.hic",
  sample_label = "PDL48"
)
