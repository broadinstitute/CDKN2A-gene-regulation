# CDKN2A-gene-regulation
This repository contains scripts for analysis of data in the CDKN2A gene
regulation project. There are three main branches of analysis:
- [CRISPRa screen](Screen)
- [Micro-C](Micro-C)
- [Confluence multiome single cell experiment](Confluence)

## Screen
The main scripts associated with the screen involve running [Sceptre](https://katsevich-lab.github.io/sceptre/index.html)
The initial processing is accomplished via [screen_processing.R](screen_processing.R), which
creates the rds object by reading and processing the mRNA alignment folders and the CROP alignment folders.
The version of Sceptre we have working, [screen_sceptre_original.R](screen_sceptre_original.R), 
dates from April 2023 (though the functions should still work
with the most recent versions, as the newer version has different function names). 

For this version of Sceptre, we are following the 
[High MOI tutorial](https://katsevich-lab.github.io/sceptre/articles/highmoi_tutorial.html).
There are two possible modes to run in - target mode or guide mode. Target mode is what is 
suggested and requires a "groups table" that assigns each guide to a target; our version
can be found [here](data/gRNAs_targets.tsv)

The other important input file is the group_pairs, essentially what we are testing. It is here
that one would modify to explore, e.g., the effect of the guides on MTAP. In this case we
test the effect of p16.dialup, p14.dialup, and CDKN2B on each of the targets or guides. To run
in target or guide mode, you will need to manually change these; it's been marked in the code 
where you would do so.

The newer version of Sceptre follows a [different tutorial](https://katsevich-lab.github.io/sceptre/reference/run_sceptre_highmoi_experimental.html); 
I have tried to run it but so far without success, see [screen_sceptre_new_version.R](screen_sceptre_new_version.R)
It is not clear how to run in "guide" mode, because negative controls must all belong to the same
"non-targeting" group (in the previous high MOI version, we set different non-targeting groups so that
there are the same number guides per target for negative controls and candidates).

The other screen scripts were used to produce plots that we looked at over the course of the
project and are included in case we want to use them in the future.

## Micro-C
There are two scripts associated with the Micro-C. [microC_juicer_commands.sh](microC_juicer_commands.sh)
gives some example commands to create the hic files, call loops, and call differential loops
using [hiccups](https://github.com/aidenlab/juicer/wiki/HiCCUPS). The R script 
[microC.R](microC.R) aimed to use DESeq to find differential areas between Hi-C contacts but this seems
to have a lot of false positives; my suggestion would be to use [HiC-DC+](https://pubmed.ncbi.nlm.nih.gov/34099725/)
to further pursue this line of inquiry. 

## Confluence
The confluence scripts are the oldest and least well-documented. 
There seem to be a lot of false starts. There is one useful plot, the Wi-38 inverse relationship plot.
We also colored the ATAC UMAP by low-middle-high association with the p14- and p16- signature scores. 
This is at the end of the file and the "ordering" files are in the Google bucket.
Here is the top level notes on how the peaks were correlated with senescence:

- Call ATAC peaks on all ATAC data (done using ArchR and MACS2; could be revisited, possibly too conservative based on how it’s done)
- Around 50 peaks found in multiome ATAC in our region of interest
- Correlate the counts in the peaks versus the senescence score on a per-peak basis (rows are the cell barcodes)
- Look peaks showing significant correlation or anti-correlation 
- Those bed files are in the Google bucket under 
p14_significant_color_atac_peaks.bed and p16_significant_color_atac_peaks.bed
p14 had both positive and negative correlations, p16 had only positive. p14 positive correlation is green and 
negative is blue; p16 positive correlation is blue.
