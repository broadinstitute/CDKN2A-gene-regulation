# CDKN2A-gene-regulation

Code for: **Chromatin topology and distal elements underlie divergent cell type-specific regulation of 9p21 locus cell cycle genes**

Elena Torlai Triglia\*, Tyler E. Miller\*, Neva C. Durand, et al.

## Data availability

Raw and processed data are available at [GEO GSE309515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309515).

## Custom GTFs

To distinguish p14ARF and p16INK4A in single-cell expression data, we generated three custom GTF files based on [GENCODE v34 (GRCh38.p13)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz). See the Methods section of the paper for full details on how each GTF was constructed.

| GTF file | Purpose | Description |
|---|---|---|
| [GSE309515_gencode.v34.5prime.gtf](https://storage.googleapis.com/broad-p16-calico/public/GSE309515_gencode.v34.5prime.gtf) | 5’ scRNA-seq analysis (Fig 2) | CDKN2A split into two separate genes: p16INK4A (ENST00000304494.9) and p14ARF (ENST00000579755.1), each with its own gene_id |
| [GSE309515_gencode.v34.3primetest.gtf](https://storage.googleapis.com/broad-p16-calico/public/GSE309515_gencode.v34.3primetest.gtf) | Testing 3’ transcript specificity (Fig 3) | CDKN2A annotated with only exon 2 and exon 3 (shared between p14 and p16), to verify that 3’ RNA-seq CDKN2A counts are not transcript-specific |
| [GSE309515_gencode.v34.dialup.gtf](https://storage.googleapis.com/broad-p16-calico/public/GSE309515_gencode.v34.dialup.gtf) | Targeted enrichment (dial-up) from 3’ libraries (Fig 3) | Tailored to primer positioning for the p14/p16 enrichment protocol |

## Figure-to-code mapping

| Figure | Description | Scripts |
|---|---|---|
| Fig 1 | GTEx/FANTOM tissue expression | GTEx data preparation described in Methods |
| Fig 2 | Cell type expression from 5’ scRNA-seq | Notebooks in [`scRNA_5prime/`](#deconvolution-of-cdkn2a-in-5-scrna-seq-fig-2) |
| Fig 3 | p16/p14 enrichment from 3’ libraries | [`Figure_3C.py`](#enrichment-of-p16-and-p14-from-3-single-cell-libraries-fig-3), `Figure_3C_expanded.py`, `Brugge_comparison.ipynb`, `Brugge_umaps_dotplot.ipynb` |
| Fig 4 | WI-38 replicative senescence switch | Notebooks in [`scRNA_dialup/`](#p14-to-p16-switch-during-replicative-senescence-fig-4) |
| Fig 5 | RCMC and epigenomic landscape | [`microC_juicer_commands.sh`](#region-capture-micro-c-and-epigenomic-landscape-fig-5-7), `plotGardener_Wi38_*.R` |
| Fig 6 | CRISPRa screen and validation | [`screen_sceptre_original.R`](#crispra-screen-fig-6), `screen_processing.R`, `write_bedGraph.R` |
| Fig 7 | CRE characterization, virtual 4C, tissue accessibility | [`extract_virtual_4C.py`](#region-capture-micro-c-and-epigenomic-landscape-fig-5-7), `Plotgardener_PDL.R` |

## Analysis

### Deconvolution of CDKN2A in 5’ scRNA-seq (Fig 2)

Published 5’ scRNA-seq datasets were remapped with the [custom 5’ GTF](#custom-gtfs) using STAR solo (via Cumulus) to quantify p14ARF and p16INK4A per cell. Each notebook processes one dataset:

- [He_et_al_expressing_cells_per_celltype_v1.ipynb](scRNA_5prime/He_et_al_expressing_cells_per_celltype_v1.ipynb) — 15-organ human atlas (He et al. 2020, GSE159929)
- [He_et_al_heatmaps_from_table_v1.R](scRNA_5prime/He_et_al_heatmaps_from_table_v1.R) — heatmaps of expression percentages
- [Elmentaite_plot_CDKN2A_p14_p16_from_Cumulus_Object_versionB_on_publishedUMAP_plot_by_sample-onlyAdult_v1.ipynb](scRNA_5prime/Elmentaite_plot_CDKN2A_p14_p16_from_Cumulus_Object_versionB_on_publishedUMAP_plot_by_sample-onlyAdult_v1.ipynb) — human gut (Elmentaite et al. 2021, E-MTAB-9543)
- [Koenig_et_al_STARsolo_UMAPs_barplots_9p21-v1.ipynb](scRNA_5prime/Koenig_et_al_STARsolo_UMAPs_barplots_9p21-v1.ipynb) — human heart (Koenig et al. 2022, GSE183852)
- [Hashimoto_assign_cell_types_umap_plot_v1.ipynb](scRNA_5prime/Hashimoto_assign_cell_types_umap_plot_v1.ipynb) — CD4 T cells from supercentenarians (Hashimoto et al. 2019)
- [PBMC_10xdata_annotate_and_plot_v1.ipynb](scRNA_5prime/PBMC_10xdata_annotate_and_plot_v1.ipynb) — 10x Genomics PBMCs

### Enrichment of p16 and p14 from 3’ single cell libraries (Fig 3)

Analysis of the Gray et al. 2022 mammary tissue 3’ scRNA-seq data, comparing standard CDKN2A detection vs. targeted dial-up enrichment of p14ARF and p16INK4A:

- [Figure_3C.py](Figure_3C.py) — bar plot comparing p14/p16 detection in standard 3’ RNA-seq vs. dial-up enrichment
- [Figure_3C_expanded.py](Figure_3C_expanded.py) — expanded version showing that CDKN2A detection is identical between the standard GTF and the exon 2+3 only GTF (i.e. 3’ quantification does not rely on transcript-specific exons)
- [Brugge_comparison.ipynb](Brugge_comparison.ipynb) — comparison of remapped expression vs. dial-up enrichment
- [Brugge_umaps_dotplot.ipynb](Brugge_umaps_dotplot.ipynb) — UMAPs and dot plots of dial-up expression by cell type

### p14 to p16 switch during replicative senescence (Fig 4)

Analysis of the WI-38 Hayflick limit time-course (Chan et al. 2022, GSE175533):

- [Chan_et_al_Hayflick_PseudoTime_SASP_Senescence-p15_ANRIL-savePlotsPdf-newSignatures_v1.ipynb](scRNA_dialup/Chan_et_al_Hayflick_PseudoTime_SASP_Senescence-p15_ANRIL-savePlotsPdf-newSignatures_v1.ipynb) — pseudotime analysis showing the p14-to-p16 switch
- [Chan_et_al_CDKN2A_CDKN2B_correlated_genes_refine_signatures_v1.ipynb](scRNA_dialup/Chan_et_al_CDKN2A_CDKN2B_correlated_genes_refine_signatures_v1.ipynb) — gene correlation analysis and signature refinement
- [Explore_Chan_Hayflick_data-redoUMAP_extract_values_forHeatmap_newSign.ipynb](scRNA_dialup/Explore_Chan_Hayflick_data-redoUMAP_extract_values_forHeatmap_newSign.ipynb) — UMAP recomputation and heatmap data extraction
- [GSE175533_sceasy_hay_onlyWT_onlu_pstime_heatmap_v1.R](scRNA_dialup/GSE175533_sceasy_hay_onlyWT_onlu_pstime_heatmap_v1.R) — pseudotime heatmap
- [Correlated_Genes/](scRNA_dialup/Correlated_Genes/) — top 200 Pearson-correlated genes for p14, p16, and CDKN2B

### Region Capture Micro-C and epigenomic landscape (Fig 5, 7)

- [microC_juicer_commands.sh](microC_juicer_commands.sh) — commands for creating .hic files, calling loops with [HiCCUPS](https://github.com/aidenlab/juicer/wiki/HiCCUPS), and calling differential loops
- [microC.R](microC.R) — DESeq-based differential contact analysis
- [extract_virtual_4C.py](extract_virtual_4C.py) — generates virtual 4C profiles from merged .hic files using hicstraw
- [plotGardener_Wi38_late_min_early_chr9_21750000_22500000_with_isoforms.R](plotGardener_Wi38_late_min_early_chr9_21750000_22500000_with_isoforms.R) — 750kb view of the locus
- [plotGardener_Wi38_late_min_early_chr9_20300000_22600000_with_isoforms.R](plotGardener_Wi38_late_min_early_chr9_20300000_22600000_with_isoforms.R) — 2.3Mb view
- [Plotgardener_PDL.R](Plotgardener_PDL.R) — PDL21/PDL48 micro-C contact maps with CPM normalization

### CRISPRa screen (Fig 6)

The CRISPRa screen used [SCEPTRE](https://katsevich-lab.github.io/sceptre/index.html) for statistical analysis of single-cell CRISPR perturbations (High MOI).

- [screen_alignment.sh](screen_alignment.sh) — STAR solo alignment calls for screen mRNA and dial-up
- [screen_processing.R](screen_processing.R) — creates the Seurat object from mRNA and CROP-seq alignment output
- [screen_sceptre_original.R](screen_sceptre_original.R) — runs SCEPTRE in target or guide mode (switch marked in code). Uses [data/gRNAs_targets.tsv](data/gRNAs_targets.tsv) for target mode. Also outputs per-gene bedGraph tracks
- [screen_sceptre_new_version.R](screen_sceptre_new_version.R) — newer SCEPTRE version
- [write_bedGraph.R](write_bedGraph.R) — generates per-gene bedGraph tracks from SCEPTRE results for IGV visualization
- [screen_plotting_cells_guides.R](screen_plotting_cells_guides.R), [screen_plotting_promoter_dialup.R](screen_plotting_promoter_dialup.R) — exploratory screen plots
- [combine_max_p_value.py](combine_max_p_value.py) — combines p-values across targets

### Total RNA-seq

- [replSen_Salmon_mapping_merge_aggregateTPM_atGeneLevel_cellcycle_senmayo.R](replSen_Salmon_mapping_merge_aggregateTPM_atGeneLevel_cellcycle_senmayo.R) — Salmon quantification of total RNA-seq across WI-38 PDLs, aggregation to gene-level TPMs, cell cycle and SenMayo scoring

### Confluence multiome

Analysis of WI-38 multiome (ATAC + RNA) data:

- [confluence_processing.R](confluence_processing.R) — data processing
- [confluence_high_quality_cell_plots.R](confluence_high_quality_cell_plots.R), [confluence_p14_p16_histogram_plots.R](confluence_p14_p16_histogram_plots.R) — visualization

### Interactive genome browser

An [IGV](https://igv.org/) session file is provided for interactive exploration of the 9p21 locus, including RCMC loops, CTCF and H3K27ac ChIP-seq, ATAC-seq (PDL20 and PDL50), SCEPTRE guide significance tracks, and gRNA target regions (Fig 5, 6, 7):

- [CDKN2A_CDKN2B_9p21.xml](CDKN2A_CDKN2B_9p21.xml) — open in [IGV desktop](https://igv.org/) or load in [IGV web app](https://igv.org/app/)

<!-- TODO: add IGV web app direct link with session URL once XML is finalized -->

## Software

Key tools and versions used (see Methods for complete details):

- [**STAR**](https://github.com/alexdobin/STAR) 2.7.10b — single-cell alignment (via [Cumulus](https://cumulus.readthedocs.io/) wrapper)
- [**SCEPTRE**](https://katsevich-lab.github.io/sceptre/index.html) — CRISPRa screen analysis (High MOI)
- [**Juicer**](https://github.com/aidenlab/juicer) 2.0 — RCMC data processing and loop calling with [HiCCUPS](https://github.com/aidenlab/juicer/wiki/HiCCUPS)
- [**Salmon**](https://combine-lab.github.io/salmon/) 1.10.1 — total RNA-seq quantification
- [**ChromBPNet**](https://github.com/kundajelab/chrombpnet) 0.1.7 — ATAC-seq attribution scores
- [**PlotGardener**](https://phanstiellab.github.io/plotgardener/) — multi-panel genomic figure layout in R
- [**IGV**](https://igv.org/) — interactive genome browser for track visualization
- [**scanpy**](https://scanpy.readthedocs.io/) 1.9.1+ / [**Seurat**](https://satijalab.org/seurat/) — single-cell analysis
- R 4.4.2, Python 3.x
