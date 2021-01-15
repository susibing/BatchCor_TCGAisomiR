#!/usr/bin/env Rscript
library(optparse)
library(data.table)
library(bsub)

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
Description = 'Remove isomiR batch effects from TCGA data. Usage of the limma R package and the 
ComBat function of the sva R package.'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## CONFIG PART ----
## Further parameter settings 
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$normalization_type = "rpm"
conf$filter = "sum" # parameter should be %in% c("median", "sum") 
conf$filter_value = 0 # 15 for median filter, 0 for sum filter
conf$include_all_samples = FALSE # if normal samples should be included in batch correction, set to T
conf$reads_filter = TRUE # set to FALSE in case no reads filter should be applied
conf$tot.reads.filter = 1000000 # threshold reads filter

data_type = "batch_cor"

## Submit job to LSF cluster
bsub_chunk(name = paste0("remove.batch.effects.", data_type, ".", conf$filter, conf$filter_value, 
                         conf$include_all_samples),
           variables = c("conf", "data_type"), 
           memory = 30, hour = 24, core = 8, 
           packages = c("tidyverse", "data.table", "stats", "RColorBrewer"),{
  
  ## enable more working memory 
  options("expressions"=500000)
  
  ## load libraries
  library(foreach,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(sva,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(iterators,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(matrixStats,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(wrapr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')   
  library(parallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(doParallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(plyr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(ggpubr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(matlab,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(limma,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  registerDoParallel(cores=8) 
  
  source(file.path(conf$base_dir, "scr", "functions.batch_cor.R")) # load functions
  
  ## Correction parameters: 
  ## Optimal combination of batch correction for the cohorts
  conf.batch_cor = fread(file.path(conf$base_dir, "data", "config.batch_cor.txt"))
  ## filter down to parameter combinations depending on whether normal data is included or not
  if(conf$include_all_samples == TRUE){
    conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == TRUE,]
  } else if(conf$include_all_samples == FALSE){
    conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == FALSE,]
  }
  
  ## go through different configurations (in parallel)
  llply(1:nrow(conf.batch_cor), function(i){
    conf$Project = conf.batch_cor$Project[i]
    conf$platform = conf.batch_cor$platform[i]
    conf$rm_batches = conf.batch_cor$rm_batches[i] # remove individual plates 
    conf$batch_1 = conf.batch_cor$batch_1[i]
    conf$batch_2 = conf.batch_cor$batch_2[i]
    if(conf$include_all_samples){
      conf$var_of_interest = "sample"
    } else {
      conf$var_of_interest = conf.batch_cor$var_of_interest[i]
    }
    
    ## create output folders
    if(conf$reads_filter){
      folder_name = paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", data_type)
    } else {
      folder_name = paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", data_type, 
                           "_wo_reads_filter")
    }

    ## LOAD DATA AND CONFOUNDER MATRIX ----
    exprs_df = load_data(path = paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                       conf$Project, "_", conf$normalization_type, ".csv"))

    ## remove individual plates 
    if(!is.na(conf$rm_batches)){
      for(j in 1:length(conf$rm_batches)){
        exprs_df = exprs_df[!(grepl(conf$rm_batches[j], rownames(exprs_df))),]
      }
    }
  
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
      exprs_df = exprs_df[,which(grepl('[0-9A-Z]{2}.[0-9A-Z]{4}.02[0-9A-Z]{1}.',colnames(exprs_df)) 
                             | grepl('[0-9A-Z]{2}.[0-9A-Z]{4}.01[0-9A-Z]{1}.',colnames(exprs_df)))]
    }
    
    ## apply filter
    ## different filtering approaches are implemented and can be chosen in the config part 
    ## sum: row sum, median: row median, percent: fraction of patients expressing miR
    if(conf$filter == "sum"){
      exprs_df$sum = rowSums(exprs_df)
      exprs_df = exprs_df[exprs_df$sum > conf$filter_value,]
      exprs_df$sum = NULL
    } else if(conf$filter == "median"){
      exprs_df = as.matrix(exprs_df)
      median.val = rowMedians(exprs_df)
      counts = which(median.val > conf$filter_value)
      exprs_df = exprs_df[counts, ]
      exprs_df = as.data.frame(exprs_df)
    } 
    
    ## LOG-transform data
    exprs_log = log2((as.matrix(exprs_df) + (min(exprs_df[exprs_df > 0]) / 10)))
    exprs_log = t(exprs_log)
    
    confounders = read.table(paste0(conf$base_dir, "/data/TCGA_confounder_table.txt"), header=T, 
                             stringsAsFactors = F)
    confounder.df = prep_conf_df(exprs_log = exprs_log, confounders = confounders, conf = conf)
    exprs_df = exprs_df[,which(names(exprs_df) %in% confounder.df$ids)]
    exprs_log = exprs_log[which(rownames(exprs_log) %in% confounder.df$ids),]
    confounder.df = confounder.df[match_order(confounder.df$ids, rownames(exprs_log)),]
    
    ## columns as factor variables (for later plotting of the data)
    for(j in 1:ncol(confounder.df)){
      if(names(confounder.df[j]) != "total.reads"){
        confounder.df[,j] = as.factor(confounder.df[,j])
        if(length(which(is.na(confounder.df[,j]))) > 0){
          confounder.df[,j] = fct_explicit_na(confounder.df[,j], na_level = "NA")
        }
      }
    }
    
    ## COMBAT ----
    ## In case normal samples are included in batch correction, second batch variable 
    ## cannot be purity => 100% confounding. 
    if(!(conf$batch_2 == "purity" & conf$var_of_interest == "sample")){ 
      batch = as.factor(confounder.df[,conf$batch_1])

      if(conf$var_of_interest == "none"){ # include variable of interest 
        modcombat = model.matrix(~1, data=confounder.df) 
      } else if(conf$var_of_interest == "Subtype_Selected"){
        modcombat = model.matrix(~as.factor(Subtype_Selected), data=confounder.df) 
      } else if (conf$var_of_interest == "sample"){
        modcombat = model.matrix(~as.factor(sample), data=confounder.df) 
      }
      
      mod0 = model.matrix(~1,data=confounder.df)
      
      combat_exprs = ComBat(dat=t(exprs_log), batch=batch, mod=modcombat, prior.plots = F) 

      if(conf$batch_2 != "none"){
        batch = as.factor(confounder.df[,conf$batch_2])
        combat_exprs = ComBat(dat=combat_exprs, batch=batch, mod=modcombat, prior.plots = F) # exprs.df: rownames= ids
      }
      combat_df = as.data.frame(combat_exprs)
      
      ## LIMMA ----
      if(conf$batch_2 != "none"){
        batch2 = as.factor(confounder.df[,conf$batch_2])
        limma_exprs = removeBatchEffect(t(exprs_log), batch=batch, design=modcombat) 
        limma_exprs = removeBatchEffect(limma_exprs, batch=batch2, design=modcombat) 
      } else{
        limma_exprs = removeBatchEffect(t(exprs_log), batch, design=modcombat)
      }
      
      limma_df = as.data.frame(limma_exprs)
      
      if(!dir.exists(paste0(folder_name, "/limma"))){dir.create(paste0(folder_name, "/limma"), 
                                                                recursive = T)}
      if(!dir.exists(paste0(folder_name, "/combat"))){dir.create(paste0(folder_name, "/combat"), 
                                                                 recursive = T)}
      
      limma_exp_df = 2^(limma_df)
      data.table::fwrite(limma_exp_df, paste0(folder_name,"/limma/", conf$Project,"_", conf$filter, 
                                              conf$filter_value, "_",conf$var_of_interest, "_", 
                                              conf$batch_1, "_", conf$batch_2, "_", conf$platform, 
                                              "_", paste0(c(conf$rm_batches), collapse = "", sep="_"), 
                                              "all_samples_",conf$include_all_samples, "_" , 
                                              conf$normalization_type, ".csv"), 
                         row.names = TRUE, col.names = TRUE)
      
      combat_exp_df = 2^(combat_df)
      data.table::fwrite(combat_exp_df, paste0(folder_name,"/combat/", conf$Project,"_", conf$filter, 
                                               conf$filter_value, "_",conf$var_of_interest, "_", 
                                               conf$batch_1, "_", conf$batch_2, "_", conf$platform, 
                                               "_", paste0(c(conf$rm_batches), collapse = "", sep="_"), 
                                               "all_samples_",conf$include_all_samples, "_" , 
                                               conf$normalization_type, ".csv"), 
                         row.names = TRUE, col.names = TRUE)
      
    }
  }, .parallel = T)    
})
