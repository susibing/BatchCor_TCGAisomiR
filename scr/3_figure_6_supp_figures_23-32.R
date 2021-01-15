#!/usr/bin/env Rscript
library(optparse)
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
Description = 'This script is to produce tSNE plots for figure 6 and supplementary figures 23-34.
Input: count/rpm-normalized expression matrices before and after batch correction using tumor and 
normal samples where available.'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## Script run with R 4.0.0

## PARAMETER SETTINGS ----
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$output_dir = file.path(conf$base_dir, "results", "tSNE_plots")
conf$normalization_type = "rpm"
conf$filter = "median" # median or sum filter
conf$filter_value = 15 # filter threshold (15: median, 0: sum)
conf$include_all_samples = T # TRUE in case normal samples should be included
conf$reads_filter = TRUE # TRUE in case total reads filter should be applied
conf$tot.reads.filter = 1000000 # threshold of total reads filter

if(!dir.exists(conf$output_dir)){dir.create(conf$output_dir, recursive = TRUE)}

data_type = "batch_cor" # "before_batch_cor" (before) or "batch_cor" (after)

## Submit job to LSF cluster
bsub_chunk(name = paste0("identify.batch.effects.", data_type, ".", conf$filter.method, 
                         conf$filter.value, conf$include_all_samples),
 variables = c("conf", "data_type"), 
 memory = 40, hour = 18, core = 8, 
 packages = c("data.table", "dplyr", "tidyverse",  "stats", "RColorBrewer"),{

   ## load functions and libraries
   source(file.path(conf$base_dir, "scr", "functions.batch_cor.R")) # load functions 
   
   library(matrixStats,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
   library(doParallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
   library(plyr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
   library(wrapr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')   
   library(Rtsne,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')   
   
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
   llply(1:nrow(conf.batch_cor), function(i){
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
       # in case there are no or hardly any normal samples 
       # (concerns CESC, COAD, LGG and OV datasets), 
       # batch correction including normal samples was not performed.
       no_normal <- c("TCGA-CESC", "TCGA-COAD", "TCGA-LGG", "TCGA-OV")
       if(conf$Project %in% no_normal){
         conf$var_of_interest = "none"
         conf$include_all_samples = F 
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
       figure_path = file.path(conf$output_dir, paste0(conf$Project, "_", data_type, "_", 
                                                       "include.all.samples", 
                                                       conf$include_all_samples, "_", 
                                                       conf$normalization_type, "_", 
                                                       conf$filter, conf$filter_value,"_", 
                                                       conf$platform, "_"))
       exprs_df = load_data(path = paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                          conf$Project, "_", conf$normalization_type, ".csv"))
                                          
     } else {
       ## batch corrected data file
       figure_path = file.path(conf$output_dir, paste0(conf$Project, "_", data_type, "_", 
                                                       conf$correction_method, "_", 
                                                       "include.all.samples", 
                                                       conf$include_all_samples, "_", 
                                                       conf$normalization_type, "_", conf$filter,
                                                       conf$filter_value, "_",
                                                       conf$var_of_interest, "_", conf$batch_1, "_", 
                                                       conf$batch_2, "_", conf$platform, "_", "rm",
                                                       paste0(c(conf$rm_batches), collapse = "", 
                                                              sep="_")))
                                                              
       if(conf$reads_filter){
         exprs_df = load_data(path = paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                            data_type, "/", conf$correction_method, "/", 
                                            conf$Project,"_", conf$filter, conf$filter_value, "_",
                                            conf$var_of_interest, "_", conf$batch_1, "_", 
                                            conf$batch_2, "_", conf$platform, "_", 
                                            paste0(c(conf$rm_batches), collapse = "", 
                                                   sep="_"),"all_samples_",conf$include_all_samples, 
                                            "_" , conf$normalization_type, ".csv"))
       }else{
         exprs_df = load_data(path = paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                            data_type, "_wo_reads_filter", "/", 
                                            conf$correction_method, "/", conf$Project,"_", 
                                            conf$filter, conf$filter_value, "_",
                                            conf$var_of_interest, "_", conf$batch_1, "_", 
                                            conf$batch_2, "_", conf$platform, "_", 
                                            paste0(c(conf$rm_batches), collapse = "", 
                                                   sep="_"),"all_samples_",conf$include_all_samples, 
                                            "_" , conf$normalization_type, ".csv"))
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
                                  | grepl('[0-9A-Z]{2}.[0-9A-Z]{4}.01[0-9A-Z]{1}.', 
                                          colnames(exprs_df)))]
     }
     
     ## apply filter
     ## different filtering approaches are implemented in the config part 
     ## sum: row sum, median: row median
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
     
     ## load confounders matrix; creation of confounder matrix described in 
     ## confounders_matrix_batch_cor.R
     confounders = read.table(paste0(conf$base_dir, "/data/TCGA_confounder_table.txt"), header=T, 
                              stringsAsFactors = F)
     confounder.df = prep_conf_df(exprs_log=exprs_log, confounders=confounders, conf=conf)
     exprs_df = exprs_df[,which(names(exprs_df) %in% confounder.df$ids)]
     exprs_log = exprs_log[which(rownames(exprs_log) %in% confounder.df$ids),]
     confounder.df = confounder.df[match_order(confounder.df$ids, rownames(exprs_log)),]
     
     ## columns as factor variables (for later plotting of the data)
     for(i in 1:ncol(confounder.df)){
       if(names(confounder.df[i]) != "total.reads"){
         confounder.df[,i] = as.factor(confounder.df[,i])
         if(length(which(is.na(confounder.df[,i]))) > 0){
           confounder.df[,i] = fct_explicit_na(confounder.df[,i], na_level = "NA")
         }
       }
     }
     
     # calculations for tSNEs
     set.seed(42)
     tsne_out <- Rtsne(exprs_log, perplexity = 25, check_duplicates= FALSE)
     symb <- c(17,16)
     if (length(unique(confounder.df$sample)) ==1){ symb <- 16 }         
     tsne_res <- tsne_out$Y
     
     ## Figure 6 ----
     # t-SNE plots for LUSC for sample type,  total reads and platform
     if (conf$Project == "TCGA-LUSC"){
       plot_sample <- function(tsne_res) {
         colors <- c( "darkorange", "darkred")
         plot(tsne_res, pch=16, xlab= "tSNE(1)", ylab= "tSNE(2)", 
              col=colors[confounder.df$sample], 
              main = c(conf$Project, " SampleType"),las=1, cex=1.2, cex.axis= 2,cex.lab=2) 
         plot(1,0, xaxt="n", yaxt="n", xlab = "", ylab= "", bty= "n", col = "white")
         legend("left",inset = c(-0.15,0), xpd = TRUE, ncol= 2, cex= 2,  pch= 16, bty="n", 
                col = c( "darkred", "darkorange"),
                legend= c("tumour", "normal"))
       }
       
       plot_totalReads <- function(tsne_res) {
         colors = c("#1f78b4", "#a6cee3", "#b2df8a", "#33a02c")
         names(colors) = sort(unique(confounder.df$quantile))
         par(mar=c(4.1, 4.1, 5.1, 5.1))
         plot(tsne_res, pch=16, xlab="tSNE(1)", ylab="tSNE(2)", 
              col=colors[confounder.df$quantile],
              main=c(conf$Project, " total_reads"), las=1, cex=1.2, cex.axis= 2,cex.lab=2)
         plot(1,0,       xaxt="n", yaxt="n", xlab = "", ylab= "", bty= "n", col = "white")
         legend("left", bty = "n", col = c("white", colors), pch= 16, 
                legend= c("quartile","1st", "2nd", "3rd", "4th"),
                inset = c(-0.15,0), xpd = TRUE, ncol= 1, cex= 2)
       }
       
       plot_platformLUSC <- function(tsne_res){
         colors <- c("#d73027","#4575b4", "#fc8d59")
         names(colors)<- sort(unique(confounder.df$Platform))
         plot(tsne_res, pch=16, xlab="tSNE(1)", ylab="tSNE(2)", 
              col=colors[confounder.df$Platform],
              main = c(conf$Project, " Platform"), las =1, cex=1.2, cex.axis= 2,cex.lab=2)
         plot(1,0, xaxt="n", yaxt="n", xlab = "", ylab= "", bty= "n", col = "white")
         legend("left",  pch= 16, bty = "n", col = colors, cex=2,inset = c(-0.15,0), 
                xpd = TRUE, ncol = 3, legend= names(colors))
       }
       
       pdf(paste0(figure_path, "_tSNEs_LUSC.pdf"))
       plot_sample(tsne_res)
       plot_totalReads(tsne_res)
       plot_platformLUSC(tsne_res)
       dev.off()
     }  
     
     
     ## Supplementary Figure 23-34 ----
     ## t-SNE plots for batch correction factors for all entities
     plot_plate <- function(tsne_res) {
       colors=colorRampPalette(c("#313695","#fee090","#a50026"))(length(unique(confounder.df$plate)))
       names(colors) = sort(unique(confounder.df$plate))
       par(mar=c(4.1, 4.1, 5.1, 5.1))
       plot(tsne_res, pch=symb[confounder.df$sample], xlab= "tSNE(1)", ylab= "tSNE(2)", 
            col=colors[confounder.df$plate],
            main = c(conf$Project, " Plate"), las =1, cex=1.2, cex.axis= 2,cex.lab=2)
       plot(1,0, xaxt="n", yaxt="n", xlab = "", ylab= "", bty= "n", col = "white")
       legend("left", bty = "n", col = c(colors,rep("black",length(symb))),
              pch= c(rep(16, length(unique(confounder.df$plate))), 16,17),
              legend= c(names(colors),c("tumor", "normal" )),
              inset = c(-0.15,0), xpd = TRUE, ncol= 3, cex= 2)
     }
     
     plot_purity <- function(tsne_res){
       colors <- colorRampPalette(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6",
           "#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"))(length(unique(confounder.df$purity)))
       names(colors) <- sort(unique(confounder.df$purity)) 
       par(mar=c(4.1, 4.1, 5.1, 5.1))
       plot(tsne_res, pch=symb[confounder.df$sample], xlab= "tSNE(1)", ylab= "tSNE(2)", 
            col=colors[confounder.df$purity],
            main = c(conf$Project, " Purity"), las =1, cex=1.2, cex.axis= 2,cex.lab=2)
       plot(1,0, xaxt="n", yaxt="n", xlab = "", ylab= "", bty= "n", col = "white")
       legend("left",  pch= c(rep(16, length(unique(confounder.df$purity))), 16,17),
              bty = "n", col=c(colors,rep("black",length(symb))), 
              legend=c(names(colors),c("tumor", "normal")),
              inset = c(-0.15,0), xpd = TRUE, ncol = 3, cex= 2)
     }
     
     plot_platform <- function(tsne_res){
       colors <- c("#d73027","#4575b4", "#fc8d59")
       names(colors)<- sort(unique(confounder.df$Platform))
       plot(tsne_res, pch=symb[confounder.df$sample], xlab= "tSNE(1)", ylab= "tSNE(2)", 
            col=colors[confounder.df$Platform],
            main = c(conf$Project, " Platform"), las =1, cex=1.2, cex.axis= 2,cex.lab=2)
       plot(1,0, xaxt="n", yaxt="n", xlab = "", ylab= "", bty= "n", col = "white")
       legend("left",pch=c(rep(16, length(unique(confounder.df$Platform))),16,17),bty="n",
              col = c(colors,rep("black",length(symb))), cex=2,inset = c(-0.15,0), 
              xpd = TRUE, ncol = 3, legend= c(names(colors),c("tumor", "normal")))
     }
     
     
     pdf(paste0(figure_path, "_tSNEs.pdf"))
     plot_plate(tsne_res)
     plot_purity(tsne_res)
     plot_platform(tsne_res)
     dev.off()

   }, .parallel = T)    
})
