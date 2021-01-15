#!/usr/bin/env Rscript
library(optparse)
library(bsub)

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'plyr_1.8.6', 'optparse_1.6.6', 'data.table_1.13.0',
'genefilter_1.66.0', 'doParallel_1.0.15', 'iterators_1.0.12', 'foreach_1.5.0', 'wrapr_2.0.2',
'matrixStats_0.57.0', 'forcats_0.5.0', 'stringr_1.4.0', 'dplyr_1.0.2', 'purrr_0.3.4', 'readr_1.3.1',
'tidyr_1.1.2', 'tibble_3.0.3', 'ggplot2_3.3.2', 'tidyverse_1.3.0', 'bsub_1.0.2', 'optparse_1.6.6')
Description = 'This script is to compare tumor and normal samples before and after batch correction 
for the TCGA-LUSC dataset. Input: TCGA-LUSC count/rpm-normalized expression matrices before and 
after batch correction.'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## PARAMETER SETTINGS ----
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$output_dir = file.path(conf$base_dir, "results", "Figure7")
conf$normalization_type = "rpm"
conf$filter = "median" # median or sum filter
conf$filter_value = 15 # filter threshold (15: median, 0: sum)
conf$include_all_samples = TRUE # T in case normal samples should be included
conf$reads_filter = TRUE # T in case total reads filter should be applied
conf$tot.reads.filter = 1000000 # threshold of total reads filter

if(!dir.exists(conf$output_dir)){dir.create(conf$output_dir, recursive = T)}

data_type = "batch_cor" # "before_batch_cor" (before) or "batch_cor" (after)

## Submit job to LSF cluster
bsub_chunk(name = paste0("tumor_normal_comparison.", data_type),
           variables = c("conf",  "data_type"), 
           memory = 5, hour = 1, core = 2, packages = c("plyr", "tidyverse", "data.table"),{

  source(file.path(conf$base_dir, "scr", "functions.batch_cor.R")) # load functions 
  
  library(matrixStats,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(wrapr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')   
  # library(parallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(doParallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(genefilter,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  # library(matrixStats,lib='/home/michelsb/R/x86_64-pc-linux-gnu-library/4.0/')
  # library(wrapr,lib='/home/michelsb/R/x86_64-pc-linux-gnu-library/4.0/')   
  # library(parallel,lib='/home/michelsb/R/x86_64-pc-linux-gnu-library/4.0/')
  # library(doParallel,lib='/home/michelsb/R/x86_64-pc-linux-gnu-library/4.0/')
  # library(ggplot2, lib ='/home/michelsb/R/x86_64-pc-linux-gnu-library/4.0')
  # library(genefilter, lib ='/home/michelsb/R/x86_64-pc-linux-gnu-library/4.0')
  # library(data.table, lib = '/home/michelsb/R/x86_64-pc-linux-gnu-library/4.0')
  registerDoParallel(cores=8)
  
  ## Correction parameters: 
  ## Optimal combination of batch correction for the cohorts
  conf.batch_cor = fread(file.path(conf$base_dir, "data", "config.batch_cor.txt"))
  if(conf$include_all_samples){
   conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == TRUE,]
  } else {
   conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == FALSE,]
  }
  
  ## go though different parameter settings (in parallel)
  # only apply to LUSC dataset for this figure 
  llply(which(conf.batch_cor$Project == "TCGA-LUSC"), function(i){
   # Parameter settings
   conf$Project = conf.batch_cor$Project[i]
   conf$platform = conf.batch_cor$platform[i]
   conf$rm_batches = conf.batch_cor$rm_batches[i]
   conf$batch_1 = conf.batch_cor$batch_1[i]
   conf$batch_2 = conf.batch_cor$batch_2[i]
   conf$correction_method = conf.batch_cor$correction_method[i]
  
   # in case, all samples are used, always take sample as variable of interest
   if(conf$include_all_samples){
     conf$var_of_interest = "sample"
   } else {
     conf$var_of_interest = conf.batch_cor$var_of_interest[i]
   }
   
   if(data_type == "before_batch_cor"){ # no batch variables to put in folder name
     conf$batch_2 = "none"
     conf$var_of_interest = "none"
   } else if (data_type == "batch_cor"){
     # in case there are no or hardly any normal samples (concerns CESC, COAD, LGG and OV datasets), 
     # batch correction including normal samples was not performed.
     no_normal <- c("TCGA-CESC", "TCGA-COAD", "TCGA-LGG", "TCGA-OV")
     if(conf$Project %in% no_normal){
       conf$var_of_interest = "none"
       conf$include_all_samples = FALSE
     } 
     # correction for purity in the presence of normal samples is not possible 
     # as normal samples per definition should not contain tumor cells 
     # --> only the second confounder was used in this case. 
     if (conf$var_of_interest == "sample"){
       if (conf$batch_1 == "purity"){
         conf$batch_1 = conf$batch_2
         conf$batch_2 = "none"
       }
       if (conf$batch_2 == "purity"){conf$batch_2 = "none"}
     } 
   }
  
   ## read expression data and define output paths
   if(data_type == "before_batch_cor"){
     table_output = paste0(conf$output_dir, "/", conf$Project,"_", data_type, "_", 
                           "include.all.samples", conf$include_all_samples, "_", 
                           conf$normalization_type, "_", conf$filter, conf$filter_value)
     exprs_df = as.data.frame(fread(paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                           conf$Project, "_", conf$normalization_type, ".csv"), 
                                    sep=",", check.names = F, stringsAsFactors = F, header = T))
   } else {
     ## batch corrected data file
     table_output = paste0(conf$output_dir, "/", conf$Project, "_", data_type, "_", 
                           conf$correction_method, "_",  "include.all.samples", 
                           conf$include_all_samples, "_", conf$normalization_type, "_", conf$filter, 
                           conf$filter_value, "_",conf$var_of_interest, "_", conf$batch_1, "_", 
                           conf$batch_2)
     if(conf$reads_filter){
       exprs_df = as.data.frame(fread(paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                             data_type, "/", conf$correction_method, "/", 
                                             conf$Project,"_", conf$filter, conf$filter_value, "_",
                                             conf$var_of_interest, "_", conf$batch_1, "_", 
                                             conf$batch_2, "_", conf$platform, "_", 
                                             paste0(c(conf$rm_batches), collapse = "", sep="_"),
                                             "all_samples_",conf$include_all_samples, "_" , 
                                             conf$normalization_type, ".csv"), sep=",", 
                                      check.names = F, stringsAsFactors = F, header = T))
     }else{
       exprs_df = as.data.frame(fread(paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                             data_type, "_wo_reads_filter", "/", 
                                             conf$correction_method, "/", conf$Project,"_", 
                                             conf$filter, conf$filter_value, "_",
                                             conf$var_of_interest, "_", conf$batch_1, "_", 
                                             conf$batch_2, "_", conf$platform, "_", 
                                             paste0(c(conf$rm_batches), collapse = "", sep="_"),
                                             "all_samples_",conf$include_all_samples, "_" , 
                                             conf$normalization_type, ".csv"), sep=",", 
                                      check.names = F, stringsAsFactors = F, header = T))
     }
   }
   
   exprs_df <- makeRn(exprs_df)
   
   ## exclude normal samples if indicated in the config part
   ## sample types:
   ## 01 Primary Solid Tumor
   ## 02 Recurrent Solid Tumor
   ## 05 Additional - New Primary
   ## 06 Metastatic
   ## 10 Blood Derived Normal
   ## 11 Solid Tissue Normal
   if(!conf$include_all_samples){
     ## only include recurrent and primary solid tumor samples
     exprs_df = exprs_df[,which(grepl('[0-9A-Z]{2}.[0-9A-Z]{4}.02[0-9A-Z]{1}.',colnames(exprs_df)) | 
                                grepl('[0-9A-Z]{2}.[0-9A-Z]{4}.01[0-9A-Z]{1}.',colnames(exprs_df)))]
   }
   
   ## apply filter
   ## different filtering approaches are implemented and can be chosen in the config part 
   ## sum: row sum, median: row median
   
   if(conf$filter == "sum"){
     exprs_df$sum = rowSums(exprs_df)
     exprs_df = exprs_df[exprs_df$sum > conf$filter_value,]
     exprs_df$sum = NULL
   } else if(conf$filter == "median"){
     exprs_df = as.matrix(exprs_df)
     median.val = rowMedians(exprs_df)
     median_expr <- data.frame(row.names(exprs_df), rowMedians(as.matrix(exprs_df)))
     names(median_expr) <- c("ids", paste0("row_median_", conf$Project))
     counts = which(median.val > conf$filter_value)
     exprs_df = exprs_df[counts, ]
     exprs_df = as.data.frame(exprs_df)
   }
   
   ## LOG-transform data
   exprs_log = log2((as.matrix(exprs_df) + (min(exprs_df[exprs_df > 0]) / 10)))
   exprs_log = t(exprs_log)
   
   ## load confounders matrix
   confounders = read.table(paste0(conf$base_dir, "/data/TCGA_confounder_table.txt"), header=TRUE, 
                            stringsAsFactors = FALSE)
   confounder.df = prep_conf_df(exprs_log = exprs_log, confounders = confounders, conf = conf)
   exprs_df = exprs_df[,which(names(exprs_df) %in% confounder.df$ids)]
   exprs_log = exprs_log[which(rownames(exprs_log) %in% confounder.df$ids),]
   confounder.df = confounder.df[match_order(confounder.df$ids, rownames(exprs_log)),]
   
   if (data_type == "before_batch_cor"){
     write.csv2(exprs_df, paste0(table_output, "_exprs_df.csv"))
   }
   
   ## t-tests between tumour and normal samples -----
   exprs_df_samp_type <- t(exprs_df)
   exprs_df_samp_type <- make1col(exprs_df_samp_type)
   exprs_df_samp_type <-join(confounder.df[,c("ids", "sample")], exprs_df_samp_type, type = "right")
   
   # calculate ttest comparing tumour versus normal for each isomiR over all patients
   calc_ttest <- function(x){
     ttest <- t.test(x[exprs_df_samp_type$sample =="Tumor"],x[exprs_df_samp_type$sample =="Normal"])
     c(colnames(x), ttest$statistic, ttest$p.value,
       mean(log2(x[exprs_df_samp_type$sample =="Tumor"]+ 0.00001)) - 
         mean(log2(x[exprs_df_samp_type$sample=="Normal"]+ 0.0001)))
   }
   
   res_ttest <- apply(exprs_df_samp_type[3:ncol(exprs_df_samp_type)],2, calc_ttest) 
   res_ttest_df <- as.data.frame(t(res_ttest))
   names(res_ttest_df) <- c( "t", "p", "log2FC")
  
   # calculate bonferroni correction
   res_ttest_df$padj <- p.adjust(res_ttest_df$p, method = "BH") 
   names(res_ttest_df) <- c(paste0("t_",data_type), paste0("p_",data_type), 
                            paste0("log2TvsN_",data_type), paste0("padj_",data_type))
   
   # export results
   write.csv2(res_ttest_df, paste0(table_output, "_ttest.csv"))
  
  }, .parallel = T)  
})
