# On the Impact of Batch Effect Correction in TCGA IsomiR Expression Data

Repository to reproduce analyses from "On the Impact of Batch Effect Correction in TCGA IsomiR Expression Data". Upon execution of the scripts, TCGA data is downloaded to the /data folder, and figures or analysis results are stored in a /results folder. 

## Author affiliations

Susanne Ibing<sup>1,a</sup>, Birgitta Michels<sup>2</sup>, Moritz Mosdzien<sup>2</sup>, Helen R. Meyer<sup>2</sup>,Lars Feuerbach<sup>1</sup>, Cindy Körner<sup>2,</sup>*

<sup>1</sup> Division of Applied Bioinformatics, German Cancer Research Center (DKFZ), Berliner Straße 41, 69120 Heidelberg, Germany

<sup>2</sup> Division of Molecular Genome Analysis, German Cancer Research Center (DKFZ), Im Neuenheimer Feld 580, 69120 Heidelberg, Germany

\* To whom correspondence should be addressed. Email: c.koerner@dkfz.de

<sup>a</sup> Current Address: Susanne Ibing, Digital Health Center, Hasso Plattner Institute, University of Potsdam, Prof.-Dr.-Helmert-Str. 2-3, 14482 Potsdam, Germany

## Requirements

R scripts were run using R versions 3.6.0 and made use of the R packages `TCGAbiolinks` (version 2.12.6), `stats` (version 3.6.0), `ggplot2` (version 3.3.2), `limma` (version 3.40.6), `sva` (version 3.32.1), `Rtsne` (version 0.15), and `heatmaps.2` (version 3.0.4). Further dependencies are stated at the header of each indivdual script. 

## Data (/data)

* Harmonized isomiR and miRNA Expression Quantification Data downloaded from the GDC portal (https://portal.gdc.cancer.gov/) (script `1_data_download.R`)
* Genomic locations of mature miRNAs were extracted from the miRBase (version 22.1, http://www.mirbase.org/) (data/hsa_isomirs_v22_R.csv)
* data/TCGA_confounder_table.txt: contains information about potential confounders (Information about Illumina platform used for sequencing of the samples was extracted from the supplementary information by Thorsson et al., 2018 (doi: 10.1016/j.immuni.2018.03.023, https://gdc.cancer.gov/about-data/publications/panimmune))
* data/isomiRs_sequences.csv: contains sequences of isomiRs obtained from miR base which are needed for Figure 4. 


## Execution (/scr)

* `functions.batch_cor.R`: Stores functions required for batch correction and subsequent analyses. 

### Data download

`1_data_download.R`: Download and preprocessing of isomiR and miRNA Expression Quantification data from 16 TCGA projects.

**Figures:**
* `1_figure_1c.R`: Script to reproduce figure 1c (hsa-let-7a-5p isomiR expression before and after annotation correction).
* `1_figure_1d.R`: Script to reproduce figure 1d (calculate distribution of canonical miRs among all isomiRs before and after annotation correction).
* `1_figure_2ghi.R`: Script to reproduce figure 2g-i (box plots of sequencing depth per plate colored by sequencing platform; scatter plot showing correlation between sequencing depth and detected isomiRs, box plots with relative isomiR detection per isomiR type colored by sequencing depth).

### Removal of batch effects

`2_remove_batch_effects.R`: Sequential removal of up to two batch variables using the ComBat function from the sva R package and the removeBatchEffect function from the limma R package. 

### Assessment of batch effects

`3_assess_batch_effects.R`: Assessment of batch effects in the isomiR expression data set by calculating association between potential confounding variables and principal components 1-10 of the log2 transformed expression data.

`3_assess_batch_effects_miRNA.R`: Assessment of batch effects in the miRNA expression data set by calculating association between potential confounding variables and principal components 1-10 of the log2 transformed expression data.

**Figures:**
* `3_assess_batch_effects.R`: Includes script for PCA plots (scatter plots (figure 2a-b) and boxplots (figure 2d-e)).
* `3_figures_comparison_batcheffects.R`: Script to reproduce comparative plots (showing -log(q) of statistical association between PCs and variables): Figure 2c/2f and Supplementary figures S2, S5, S7-S21.
* `3_figure_5.R`: Script to reproduce figure 5 and supplementary figure S22 (bubble plots providing pan-cancer overview of success of batch correction)
* `3_figure_6_supp_figures_23-32.R`: Script to reproduce figure 6 and Supplementary figures S23-S32. 

### Identification of outlier isomiRs 

* `4_define_outliers_platforms_LUSC.R`: Comparison between samples sequenced on the sequencing platforms Gene Analyzer (GA) or HiSeq before and after batch correction for the TCGA-LUSC dataset. 
* `4_define_outliers_platforms_all_entities.R`: Comparison between samples sequenced on the sequencing platforms using all datasets containing samples sequenced on GA to investigate patterns in the data before batch correction. 

**Figures:**
* `4_define_outliers_platforms_LUSC.R`: Includes script to reproduce Figure 3a,b,d and e and supplementary Figures 6a and b.
* `4_figure_3c.R`: Script to reproduce figure 3c.
* `4_figure4.R`: Script to reproduce figure 4.  


### Expression analysis tumor vs. normal 

`5_expression_analysis_tumor_vs_normal.R`: Analyse differences in expression between tumor and normal samples before and after batch correction. 

**Figures:**
* `5_figure_7.R`: Script to reproduce figure 7. 

## References 
The makeTransparent function was resumed from Nick Sabbe from a "Stack Overflow" post (https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color). 
The countCharOccurrences function was resumed from a techoverflow post (https://techoverflow.net/2012/11/10/r-count-occurrences-of-character-in-string/). 
