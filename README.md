# KRAS-RIT1-profiling

This repository contains the code used to process and analyze the data for the KRAS and RIT1 multi-omic profiling project (Lo et al., *in revision*).

Included:
- Pipeline for RNA-seq read alignment and transcript quantification: Snakefile, config.yaml
- QC analysis of RNA-sequencing, alignments, and quantification: analysis_QC.R
- Analyses of differential expression: do_DE.R, analysis_DE.R
- Analyses of transcription factor target enrichment: do_ChEA.R, analysis_ChEA.R
- Analyses of LC-MS/MS proteomic and phosphoproteomic data:
- Re-analysis of L1000 data of the same genes in a separate screen (Berger et al., *2016*)