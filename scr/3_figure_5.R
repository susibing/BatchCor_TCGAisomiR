#!/usr/bin/env Rscript
library(optparse)

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'forcats_0.5.0', 'stringr_1.4.0', 'dplyr_1.0.2', 
'readr_1.3.1', 'tidyr_1.1.2', 'tibble_3.0.3', 'ggplot2_3.3.2', 'tidyverse_1.3.0', 
'doParallel_1.0.15', 'optparse_1.6.6', 'limma_3.40.6', 'sva_3.32.1', 'BiocParallel_1.18.1', 
'data.table_1.13.0', 'bsub_1.0.2', 'genefilter_1.66.0', 'mgcv_1.8-33', 'nlme_3.1-149', 
'wrapr_2.0.2', 'purrr_0.3.4')
Description = 'Script to create bubble plots of Figure 5 and supplementary figure S22.'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}


## LIBRARIES ----
library(RColorBrewer)
library(ggplot2)
library(scales)
library(data.table)
library(plyr)

## PARAMETER SETTINGS ----
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$output_dir = file.path(conf$base_dir, "results/bubbleplots/")

if(!dir.exists(conf$output_dir)){dir.create(conf$output_dir, recursive = T)}
source(file.path(conf$base_dir, "scr", "functions.batch_cor.R")) # load functions 

## Correction parameters: 
## Optimal combination of batch correction for the cohorts
conf.batch_cor = fread(file.path(conf$base_dir, "data", "config.batch_cor.txt"))

## TCGA projects including (sufficient) normal samples with isomiR Expression Quantification data 
reads_filter = T # Reads filter (threshold 1,000,000 mapped reads per sample) set to TRUE

## EXECUTION ----
pvals_tumor_med = get_pvals(include_all_samples = FALSE, conf.batch_cor = conf.batch_cor, 
                            filter.method = "median", filter.value = 15)
pvals_all_med = get_pvals(include_all_samples = TRUE, conf.batch_cor = conf.batch_cor, 
                          filter.method = "median", filter.value = 15)
pvals_tumor_sum = get_pvals(include_all_samples = FALSE, conf.batch_cor = conf.batch_cor, 
                            filter.method = "sum", filter.value = 0)
pvals_all_sum = get_pvals(include_all_samples = TRUE, conf.batch_cor = conf.batch_cor, 
                          filter.method = "sum", filter.value = 0)

variables = c("plate", "platform", "purity", "reads", "subtype")
PCs = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")


## MEDIAN FILTER
# only tumor: 
p_before1 = bubble_plot(df = pvals_tumor_med, legend = TRUE, ticks = TRUE, time_point = "before", 
                        PCs = PCs, variables = variables, filter = "median")
p_after1 = bubble_plot(df = pvals_tumor_med, legend = TRUE, ticks = TRUE, time_point = "after", 
                       PCs = PCs, variables = variables, filter = "median")
ggsave(plot=p_before1,height=20,width=6.7,dpi=300, 
       filename=paste0(conf$output_dir, "bubble_before_tumor_only_med.pdf"), useDingbats=FALSE)                      
ggsave(plot=p_after1,height=20,width=6.7,dpi=300, 
       filename=paste0(conf$output_dir, "bubble_after_tumor_only_med.pdf"), useDingbats=FALSE)                      

# including normal:
p_before2 = bubble_plot(df = pvals_all_med, legend = TRUE, ticks = TRUE, time_point = "before", PCs = PCs, 
                        variables = variables, filter = "median")
p_after2 = bubble_plot(df = pvals_all_med, legend = TRUE, ticks = TRUE, time_point = "after", PCs = PCs, 
                       variables = variables, filter = "median")
ggsave(plot=p_before2,height=20,width=6.7,dpi=300, 
       filename=paste0(conf$output_dir, "bubble_before_incl_normal_med.pdf"), useDingbats=FALSE)                      
ggsave(plot=p_after2,height=20,width=6.7,dpi=300, 
       filename=paste0(conf$output_dir, "bubble_after_incl_normal_med.pdf"), useDingbats=FALSE)                      

## SUM FILTER
## only tumor:
p_before3 = bubble_plot(df = pvals_tumor_sum, legend = TRUE, ticks = TRUE, time_point = "before", 
                        PCs = PCs, variables = variables, filter = "sum")
p_after3 = bubble_plot(df = pvals_tumor_sum, legend = TRUE, ticks = TRUE, time_point = "after", 
                       PCs = PCs, variables = variables, filter = "sum")
ggsave(plot=p_before3,height=20,width=6.7,dpi=300, 
       filename=paste0(conf$output_dir, "bubble_before_tumor_only_sum.pdf"), useDingbats=FALSE)
ggsave(plot=p_after3,height=20,width=6.7,dpi=300, 
       filename=paste0(conf$output_dir, "bubble_after_tumor_only_sum.pdf"), useDingbats=FALSE)

## including normal:
p_before4 = bubble_plot(df = pvals_all_sum, legend = TRUE, ticks = TRUE, time_point = "before", 
                        PCs = PCs, variables = variables, filter = "sum")
p_after4 = bubble_plot(df = pvals_all_sum, legend = TRUE, ticks = TRUE, time_point = "after", 
                       PCs = PCs, variables = variables, filter = "sum")
ggsave(plot=p_before4,height=20,width=6.7,dpi=300, 
       filename=paste0(conf$output_dir, "bubble_before_incl_normal_sum.pdf"), useDingbats=FALSE)                      
ggsave(plot=p_after4,height=20,width=6.7,dpi=300, 
       filename=paste0(conf$output_dir, "bubble_after_incl_normal_sum.pdf"), useDingbats=FALSE)                      