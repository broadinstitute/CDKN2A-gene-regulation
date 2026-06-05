# GWAS / eQTL / LD analysis at 9p21 (Fig 7C)

Code and condensed reference figure for the locus-summary panel showing fine-mapped GWAS variants, eQTLs, LD blocks, and conservation across the 9p21 region.

## Files

- [`p16_locus_summary_plot_simple.ipynb`](p16_locus_summary_plot_simple.ipynb) — generates the Fig 7C track summary. Pulls fine-mapped UKBB variants, GTEx V10 eQTLs, UKBB LD (337K British-ancestry individuals, lifted to hg38), phastCons100way conservation, and Ensembl common SNPs (MAF ≥ 1%); overlays them on the gene track for `chr9:21,750,000–22,500,000`.
- [`GWAS_catalog_cancer_locus_PIPs.condensed.pdf`](GWAS_catalog_cancer_locus_PIPs.condensed.pdf) — condensed reference plot of EBI GWAS catalog cancer-germline lead SNPs and posterior inclusion probabilities (PIPs) at the locus.

## Inputs

Sources for the source tracks (full details in Methods, available in /data):

- Fine-mapped UKBB variants: [Finucane lab catalog](https://www.finucanelab.org/data).
- Cancer germline GWAS lead SNPs: [EBI GWAS catalog](https://www.ebi.ac.uk/gwas).
- Fine-mapped eQTLs: GTEx V10, per-tissue, manually overlapped per transcript.

- LD: UKBB 337K British-ancestry individuals. -- available on https://storage.googleapis.com/broad-p16-calico/public/UKBB_LD.p16_expanded_locus.parquet

## Dependencies

Requires the latest version of [pyqtl](https://github.com/broadinstitute/pyqtl):

```
git clone git@github.com:broadinstitute/pyqtl.git
cd pyqtl
pip3 install -e .
```
