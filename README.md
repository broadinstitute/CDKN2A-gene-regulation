# CDKN2A-gene-regulation

Code for: **Chromatin topology and distal elements underlie divergent cell type-specific regulation of 9p21 locus cell cycle genes**

Elena Torlai Triglia\*, Tyler E. Miller\*, Neva C. Durand, et al.

## Data availability

Raw and processed data are available at [GEO GSE309515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309515).

## Custom GTFs

To distinguish p14ARF and p16INK4A in single-cell expression data, we generated three custom GTF files based on [GENCODE v34 (GRCh38.p13)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz). See the Methods section of the paper for full details on how each GTF was constructed.

| GTF file | Purpose | Description |
|---|---|---|
| [GSE309515_gencode.v34.5prime.gtf](https://storage.googleapis.com/broad-p16-calico/public/GSE309515_gencode.v34.5prime.gtf) | 5’ scRNA-seq analysis (Fig 2) | CDKN2A split into two separate genes: p16INK4A (ENST00000304494.9) and p14ARF (ENST00000579755.1), each with its own gene_id; removed all other less-characterised isoforms |
| [GSE309515_gencode.v34.3primetest.gtf](https://storage.googleapis.com/broad-p16-calico/public/GSE309515_gencode.v34.3primetest.gtf) | Testing 3’ transcript specificity (Fig 3) | CDKN2A annotated with only exon 2 and exon 3 (shared between p14 and p16), exon 1s are annotated as different genes. This is the GTF to check if 3’ RNA-seq CDKN2A reads are mostly in the common exons. |
| [GSE309515_gencode.v34.dialup.gtf](https://storage.googleapis.com/broad-p16-calico/public/GSE309515_gencode.v34.dialup.gtf) | Targeted enrichment (dial-up) from 3’ libraries (Fig 3-4,6) | Tailored to primer positioning for the p14/p16 enrichment protocol, used on existing breast tissue and fibroblast 10x 3' libraries and CRISPRa Seq-well libraries |

## Figure-to-code mapping

| Figure | Description | Scripts |
|---|---|---|
| Fig 1 | GTEx/FANTOM tissue expression | [`GTEx_heatmaps_rev_order_addingMTAP.R`](#gtex-tissue-expression-fig-1) |
| Fig 2 | Cell type expression from 5’ scRNA-seq | Notebooks in [`scRNA_5prime/`](#deconvolution-of-cdkn2a-in-5-scrna-seq-fig-2) |
| Fig 3 | p16/p14 enrichment from 3’ libraries | [`Figure_3C_expanded.py`](#enrichment-of-p16-and-p14-from-3-single-cell-libraries-fig-3), `Brugge_comparison.ipynb`, `Brugge_umaps_dotplot.ipynb` |
| Fig 4 | WI-38 replicative senescence switch | Notebooks in [`scRNA_dialup/`](#p14-to-p16-switch-during-replicative-senescence-fig-4) |
| Fig 5 | RCMC and epigenomic landscape | [`microC_juicer_commands.sh`](#region-capture-micro-c-and-epigenomic-landscape-fig-5-7), `plotGardener_Wi38_*.R`, `Plotgardener_PDL.R` |
| Fig 6 | CRISPRa screen and validation | [`screen_sceptre_original.R`](#crispra-screen-fig-6), `screen_processing.R`, `write_bedGraph.R`, `senmayo_safeharbor_E1.R` |
| Fig 7 | CRE characterization, virtual 4C, tissue accessibility | [`extract_virtual_4C.py`](#region-capture-micro-c-and-epigenomic-landscape-fig-5-7), `Plotgardener_PDL.R` |

## Analysis

### GTEx tissue expression (Fig 1)

- [GTEx_heatmaps_rev_order_addingMTAP.R](GTEx_heatmaps_rev_order_addingMTAP.R) — ComplexHeatmap plots of CDKN2A (and p14/p16 separately), CDKN2B, CDKN2B-AS1, and MTAP expression across GTEx tissues, ordered by p16/p14 ratio. Reads `heatmap_data_Brugge.csv` (GTEx TPMs aggregated per tissue; preparation described in Methods).

### Quantification of p16 and p14 in 5’ scRNA-seq (Fig 2)

Published 5’ scRNA-seq datasets were remapped with the [custom 5’ GTF](#custom-gtfs) using STAR solo (via Cumulus) to quantify p14ARF and p16INK4A per cell. Each notebook processes one dataset:

- [He_et_al_expressing_cells_per_celltype_v1.ipynb](scRNA_5prime/He_et_al_expressing_cells_per_celltype_v1.ipynb) — 15-organ human atlas (He et al. 2020, GSE159929)
- [He_et_al_heatmaps_from_table_v1.R](scRNA_5prime/He_et_al_heatmaps_from_table_v1.R) — heatmaps of expression percentages
- [Elmentaite_plot_CDKN2A_p14_p16_from_Cumulus_Object_versionB_on_publishedUMAP_plot_by_sample-onlyAdult_v1.ipynb](scRNA_5prime/Elmentaite_plot_CDKN2A_p14_p16_from_Cumulus_Object_versionB_on_publishedUMAP_plot_by_sample-onlyAdult_v1.ipynb) — human gut (Elmentaite et al. 2021, E-MTAB-9543)
- [Koenig_et_al_STARsolo_UMAPs_barplots_9p21-v1.ipynb](scRNA_5prime/Koenig_et_al_STARsolo_UMAPs_barplots_9p21-v1.ipynb) — human heart (Koenig et al. 2022, GSE183852)
- [Hashimoto_assign_cell_types_umap_plot_v1.ipynb](scRNA_5prime/Hashimoto_assign_cell_types_umap_plot_v1.ipynb) — CD4 T cells from supercentenarians (Hashimoto et al. 2019)
- [PBMC_10xdata_annotate_and_plot_v1.ipynb](scRNA_5prime/PBMC_10xdata_annotate_and_plot_v1.ipynb) — 10x Genomics PBMCs

### Enrichment and quantification of p16 and p14 from 3’ single cell libraries (Fig 3)

Analysis of the Gray et al. 2022 mammary tissue 3’ scRNA-seq data, comparing standard CDKN2A detection vs. targeted dial-up enrichment of p14ARF and p16INK4A:

- [Figure_3C_expanded.py](Figure_3C_expanded.py) — bar plot comparing p14/p16 detection in standard 3’ RNA-seq vs. dial-up enrichment, also showing that CDKN2A detection is identical between the standard GTF and the exon 2+3 only GTF (i.e. 3’ quantification does not rely on transcript-specific exons). Reads `Brugge_5_prime_merged.csv`
- [Brugge_comparison.ipynb](Brugge_comparison.ipynb) — comparison of remapped expression vs. dial-up enrichment
- [Brugge_umaps_dotplot.ipynb](Brugge_umaps_dotplot.ipynb) — UMAPs and dot plots of dial-up expression by cell type

### Expression of p14 and p16 during replicative senescence (Fig 4)

Analysis of the WI-38 Hayflick limit time-course (Chan et al. 2022, GSE175533):

- [Chan_et_al_Hayflick_PseudoTime_SASP_Senescence-p15_ANRIL-savePlotsPdf-newSignatures_v1.ipynb](scRNA_dialup/Chan_et_al_Hayflick_PseudoTime_SASP_Senescence-p15_ANRIL-savePlotsPdf-newSignatures_v1.ipynb) — pseudotime analysis showing the p14-to-p16 switch
- [Chan_et_al_CDKN2A_CDKN2B_correlated_genes_refine_signatures_v1.ipynb](scRNA_dialup/Chan_et_al_CDKN2A_CDKN2B_correlated_genes_refine_signatures_v1.ipynb) — gene correlation analysis and signature refinement
- [Explore_Chan_Hayflick_data-redoUMAP_extract_values_forHeatmap_newSign.ipynb](scRNA_dialup/Explore_Chan_Hayflick_data-redoUMAP_extract_values_forHeatmap_newSign.ipynb) — UMAP recomputation and heatmap data extraction
- [GSE175533_sceasy_hay_onlyWT_onlu_pstime_heatmap_v1.R](scRNA_dialup/GSE175533_sceasy_hay_onlyWT_onlu_pstime_heatmap_v1.R) — pseudotime heatmap
- [Correlated_Genes/](scRNA_dialup/Correlated_Genes/) — top 200 Pearson-correlated genes for p14, p16, and CDKN2B

### Region Capture Micro-C and epigenomic landscape (Fig 5, 7)

- [microC_juicer_commands.sh](microC_juicer_commands.sh) — commands for creating .hic files, calling loops with [HiCCUPS](https://github.com/aidenlab/juicer/wiki/HiCCUPS), and calling differential loops
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
- [senmayo_safeharbor_E1.R](senmayo_safeharbor_E1.R) — compares SenMayo signature score in cells with E1 target guides vs safe-harbor guides, regressing out per-cell guide count. Uses [SenMayo_signature.txt](SenMayo_signature.txt)

### Total RNA-seq

- [replSen_Salmon_mapping_merge_aggregateTPM_atGeneLevel_cellcycle_senmayo.R](replSen_Salmon_mapping_merge_aggregateTPM_atGeneLevel_cellcycle_senmayo.R) — Salmon quantification of total RNA-seq across WI-38 PDLs, aggregation to gene-level TPMs, cell cycle and SenMayo scoring

### Interactive genome browser

[IGV](https://igv.org/) session files are provided per figure for interactive exploration of the 9p21 locus. Each session bundles the relevant tracks (CTCF and H3K27ac ChIP-seq, ATAC-seq, RCMC loops, SCEPTRE significance tracks, virtual 4C, gene annotations, etc.) and pulls them from public URLs.

| Figure | One-click in browser | Session XML |
|---|---|---|
| Fig 5B | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FFigure_5B.xml) | [`Figure_5B.xml`](igv_sessions/Figure_5B.xml) |
| Fig 6B-C | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FFigure_6B-C.xml) | [`Figure_6B-C.xml`](igv_sessions/Figure_6B-C.xml) |
| Fig 6F | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FFigure_6F.xml) | [`Figure_6F.xml`](igv_sessions/Figure_6F.xml) |
| Fig 7A | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FFigure_7A.xml) | [`Figure_7A.xml`](igv_sessions/Figure_7A.xml) |
| Suppl Fig 11 | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FSupplemental_Figure_11.xml) | [`Supplemental_Figure_11.xml`](igv_sessions/Supplemental_Figure_11.xml) |
| Suppl Fig 13 | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FSupplemental_Figure_13.xml) | [`Supplemental_Figure_13.xml`](igv_sessions/Supplemental_Figure_13.xml) |
| Suppl Fig 14 | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FSupplemental_Figure_14.xml) | [`Supplemental_Figure_14.xml`](igv_sessions/Supplemental_Figure_14.xml) |
| Suppl Fig 15A | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FSupplemental_Figure_15A.xml) | [`Supplemental_Figure_15A.xml`](igv_sessions/Supplemental_Figure_15A.xml) |
| Suppl Fig 15B | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FSupplemental_Figure_15B.xml) | [`Supplemental_Figure_15B.xml`](igv_sessions/Supplemental_Figure_15B.xml) |
| Suppl Fig 18 | [Open in IGV web app](https://igv.org/app/?sessionURL=https%3A%2F%2Fraw.githubusercontent.com%2Fbroadinstitute%2FCDKN2A-gene-regulation%2Fmain%2Figv_sessions%2FSupplemental_Figure_18.xml) | [`Supplemental_Figure_18.xml`](igv_sessions/Supplemental_Figure_18.xml) |

Working sessions used during analysis (not tied to a specific paper figure):

- [`Virtual_4C_session_2026-01-20.xml`](igv_sessions/Virtual_4C_session_2026-01-20.xml) — virtual 4C from CRE viewpoints
- [`Virtual_4C_session_2026-01-22.xml`](igv_sessions/Virtual_4C_session_2026-01-22.xml) — virtual 4C from promoter viewpoints

XMLs can also be opened directly in [IGV desktop](https://igv.org/). All track URLs in each session are publicly hosted at `gs://broad-p16-calico/`.

## Reproducing from GEO

All raw and processed data live at [GEO GSE309515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309515) as five sub-series (microC, dial-up scRNA-seq, CRISPRa screen, ChIP-seq, total RNA-seq). Scripts in this repo expect GEO file names as-is. Below is what each script reads.

### Region Capture Micro-C (Fig 5, 7)

Stage these processed files in the repo root (or symlink them):

- `merged.hic`, `microC_WI-38_PDL21.hic`, `microC_WI-38_PDL48_Rep2.hic`, `merged_hiccups_filtered.bedpe`

Then `extract_virtual_4C.py`, `Plotgardener_PDL.R`, and `plotGardener_Wi38_late_min_early_*.R` will run from the repo root.

### CRISPRa screen (Fig 6)

Download all dial-up scRNA-seq processed files (`mRNA_arrayN_{barcodes,features,matrix}.tsv.gz/.mtx.gz`, `Dialup-array{7..16}_{barcodes,features,matrix}.tsv.gz/.mtx.gz`, `new-screen-CROP-orig-stitch-array{1..6}_{barcodes,features,matrix}.tsv.gz/.mtx.gz`) into a single directory `geo_screen/` next to the scripts. Then:

```r
source("screen_processing.R")        # writes data/screen.rna.rds
source("screen_sceptre_original.R")  # or screen_sceptre_new_version.R
source("write_bedGraph.R")
```

### Total RNA-seq

Download `PDL{21,33,41,48,535,54,57}_quant.sf` into `geo_totalRNA/`. Also requires `mart_export_GRCh38.p13.txt` (download from Ensembl BioMart with columns "Gene stable ID" + "Gene name") next to the script. Then run `replSen_Salmon_mapping_merge_aggregateTPM_atGeneLevel_cellcycle_senmayo.R`.

### Dial-up scRNA-seq notebooks (Fig 3, 4)

The scRNA_5prime/ and scRNA_dialup/ notebooks read from public GEO/SRA accessions cited in the Methods. The Brugge breast tissue notebooks (`Brugge_*.ipynb`, `Figure_3C_expanded.py`) operate on per-cell summary tables derived from the dial-up matrices on GEO; the small derived inputs we use (`Brugge_5_prime_merged.csv`, `dialup_ratios.csv`, `expression_ratios.csv`, `heatmap_data_Brugge.csv`) are tracked in this repo. The cell-level p16/p14 enrichment counts (`*_p16_p14_per_cell.tsv`, expected by `Brugge_umaps_dotplot.ipynb` under `dialup_brugge/`) and the Gray et al. annotated h5ad must be regenerated from the GEO matrices and Synapse `syn26560310`.

### Analyses not in this repo

The following analyses are described in Methods but use standard pipelines or were performed by collaborators; no custom code is included here. Refer to the Methods section for parameters.

- ATAC-seq processing (bowtie2 + macs2 + deepTools CPM normalization)
- ChIP-seq processing (ENCODE uniform pipeline; HOMER peak calling)
- ChromBPNet 0.1.7 attribution scores + DeepLIFT + TOMTOM motif annotation
- DiffBind/DESeq2 selection of the 173 CRISPRa target regions
- qRT-PCR Z-score / fold-change combination across replicates
- GTEx violin plots (Fig 1F, Suppl Fig 3) and GTEx + FANTOM upset plots (Fig 1D-E)
- sci-ATAC-seq cCRE × cell type heatmap (Fig 7E)

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
