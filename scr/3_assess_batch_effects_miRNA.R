#!/usr/bin/env Rscript
library(optparse)
library(bsub)

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'RColorBrewer_1.1-2', 'forcats_0.5.0', 'stringr_1.4.0', 
'dplyr_1.0.2', 'purrr_0.3.4', 'readr_1.3.1', 'tidyr_1.1.2', 'tibble_3.0.3', 'tidyverse_1.3.0',
'ggpubr_0.4.0', 'ggplot2_3.3.2', 'plyr_1.8.6', 'doParallel_1.0.15', 'iterators_1.0.12', 
'foreach_1.5.0', 'wrapr_2.0.2', 'matrixStats_0.57.0', 'bsub_1.0.2', 'data.table_1.13.0', 
'optparse_1.6.6')
Description = 'This script is to identify batch effects in the TCGA miRNA quantification data.
Input: count/rpm-normalized expression matrices before batch correction. (Only use of uncorrected 
data)'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## CONFIG PART ----
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$normalization_type = "rpm"
conf$correction_method = c("combat", "limma")
conf$filter = "median" # median or sum filter
conf$filter_value = 15 # filter threshold (15: median, 0: sum)
conf$include_all_samples = FALSE # T in case normal samples should be included
conf$reads_filter = FALSE # T in case total reads filter should be applied
conf$tot.reads.filter = 1000000 # threshold of total reads filter

data_type = "before_batch_cor" # only before_batch_cor


bsub_chunk(name = paste0("identify.batch.effects.miRNA.", conf$filter, conf$filter_value), 
           memory = 40, hour = 18, core = 8, 
           packages = c("tidyverse", "data.table", "RColorBrewer"),{
             
  # enable more working memory 
  options("expressions"=500000)
  source(file.path(conf$base_dir, "scr", "functions.batch_cor.R")) # load functions
  
  library(foreach,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(iterators,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(matrixStats,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(wrapr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')   
  library(doParallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(plyr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(ggpubr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(matlab,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(ggforce,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(ggsignif,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(forcats)
  
  registerDoParallel(cores=8)
  
  ## Correction parameters: 
  ## Optimal combination of batch correction for the cohorts
  conf.batch_cor = fread(file.path(conf$base_dir, "data", "config.batch_cor.txt"))
  ## filter down to parameter combinations depending on whether normal data is included or not
  if(conf$include_all_samples == TRUE){
    conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == TRUE,]
  } else if(conf$include_all_samples == FALSE){
    conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == FALSE,]
  }
  
  llply(1:nrow(conf.batch_cor), function(i){ # for all tumor entities
  #llply(11, function(i){ # for LUSC only
    conf$tot.reads.filter = 1000000
    conf$Project = conf.batch_cor$Project[i]
    conf$platform = conf.batch_cor$platform[i]
    conf$rm_batches = conf.batch_cor$rm_batches[i]
    conf$batch_1 = conf.batch_cor$batch_1[i]
    conf$batch_2 = conf.batch_cor$batch_2[i]
    # in case, all samples are used, always take sample as variable of interest
    if(conf$include_all_samples){
      conf$var_of_interest = "sample"
    } else {
      conf$var_of_interest = conf.batch_cor$var_of_interest[i]
    }
    if(data_type == "before_batch_cor"){  
      conf$batch_2 = "none"
      conf$var_of_interest = "none"
    }
    
    if(!(conf$batch_2 == "purity" & conf$var_of_interest == "sample")){
    llply(conf$correction_method, function(correction_method){

      if(data_type == "batch_cor"){
        figure.path = paste0(conf$base_dir,"/results/miRNA_only/", data_type, "/", 
                             correction_method, "/", conf$Project, "/", "include.all.samples", 
                             conf$include_all_samples, "/", conf$normalization, "_", conf$filter, 
                             conf$filter_value, "_",conf$var_of_interest, "_", conf$batch_1, "_", 
                             conf$batch_2, "_",conf$platform, "_", "rm",
                             paste0(c(conf$rm_batches), collapse = "", sep="_"))
        table_output = paste0(conf$Project, "_", data_type, "_", correction_method, "_",  
                              "include.all.samples", conf$include_all_samples, "_", 
                              conf$normalization, "_", conf$filter, conf$filter_value, "_",
                              conf$var_of_interest, "_", conf$batch_1, "_", conf$batch_2)
      } else { 
        figure.path = paste0(conf$base_dir,"/results/miRNA_only/", data_type, "/", conf$Project, 
                             "/", "include.all.samples", conf$include_all_samples, "/", 
                             conf$normalization, "_", conf$filter, conf$filter_value,"_", 
                             conf$platform)
        table_output = paste0(conf$Project,"_", data_type, "_", "include.all.samples", 
                              conf$include_all_samples, "_", conf$normalization, "_", conf$filter, 
                              conf$filter_value)
      }
      dir.create(figure.path, recursive = T)
      
      if(data_type == "before_batch_cor"){
        exprs_df = read.table(paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                     conf$Project, "_", conf$normalization, "_miRNA_only.csv"), 
                              sep=",", check.names = F, stringsAsFactors = F, header = T)
      } else {
        # batch corrected data file
        # correction method missing
        exprs_df = read.table(paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                     data_type, "/", correction_method, "/", conf$Project,"_", 
                                     conf$filter, conf$filter_value, "_",conf$var_of_interest,"_", 
                                     conf$batch_1, "_", conf$batch_2, "_", conf$platform, "_", 
                                     paste0(c(conf$rm_batches), collapse = "", sep="_"),
                                     "all_samples_",conf$include_all_samples, "_" , 
                                     conf$normalization, "_miRNA_only.csv"), sep=",", 
                              check.names = F, stringsAsFactors = F, header = T)
        exprs_df[,(ncol(exprs_df)-3):ncol(exprs_df)] = NULL
      }
      
      rownames(exprs_df) = exprs_df[,1]
      exprs_df[,1] = NULL
      
      # exclude normal samples if indicated in the config part
      # sample types:
      # 01 Primary Solid Tumor
      # 02 Recurrent Solid Tumor
      # 05 Additional - New Primary
      # 06 Metastatic
      # 10 Blood Derived Normal
      # 11 Solid Tissue Normal
      if(!conf$include_all_samples){
        # only include recurrent and primary solid tumor samples
        exprs_df = exprs_df[,which(grepl('[0-9A-Z]{2}.[0-9A-Z]{4}.02[0-9A-Z]{1}.',
                                         colnames(exprs_df)) | 
                              grepl('[0-9A-Z]{2}.[0-9A-Z]{4}.01[0-9A-Z]{1}.',colnames(exprs_df)))]
      }
      
      # apply filter
      # different filtering approaches are implemented and can be chosen in the config part 
      # sum: row sum, median: row median, percent: fraction of patients expressing miR
      
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
      } else if(conf$filter == "percent"){
        exprs_df$percentage = sapply(rownames(exprs_df), function(miR){
          length(which(exprs_df[rownames(exprs_df) == miR,] != 0))/ncol(exprs_df)
        })
        exprs_df = exprs_df[exprs_df$percentage > conf$filter_value,]
        exprs_df$percentage = NULL
      } else {
        # remove all rows with only 0 entries
        exprs_df$sum = rowSums(exprs_df)
        exprs_df = exprs_df[exprs_df$sum > 0,]
        exprs_df$sum = NULL
      }
      
      ## LOG-tranform data
      exprs_log = log2((as.matrix(exprs_df) + (min(exprs_df[exprs_df > 0]) / 10)))
      exprs_log = t(exprs_log)
      
      # load confounders matrix
      confounders = read.table(paste0(conf$base_dir, "/data/TCGA_confounder_table.txt"), header=T, 
                               stringsAsFactors = F)
      confounder.df = prep_conf_df(exprs_log = exprs_log, confounders = confounders, conf = conf)
      exprs_df = exprs_df[,which(names(exprs_df) %in% confounder.df$ids)]
      exprs_log = exprs_log[which(rownames(exprs_log) %in% confounder.df$ids),]
      confounder.df = confounder.df[match_order(confounder.df$ids, rownames(exprs_log)),]
      
      ## columns as factor variables (for later plotting of the data)
      for(x in 1:ncol(confounder.df)){
        if(names(confounder.df[x]) != "total.reads"){
          confounder.df[,x] = as.factor(confounder.df[,x])
          if(length(which(is.na(confounder.df[,x]))) > 0){
            confounder.df[,x] = fct_explicit_na(confounder.df[,x], na_level = "NA")
          }
        }
      }
      
      
      ## PCA confounder analysis (PC 1-10 in scope) ----
      pca = prcomp(exprs_log)
      pca_results = as.data.frame(pca$x)
      percentage = calculate_percentage(pca_df = pca, df = pca_results)
      pca_results$ids = rownames(pca_results)
      pca_results = merge(x= pca_results, y=confounder.df, by = "ids", all.x=T, all.y=F)
      
      confounders = c("plate", "Platform", "total.reads",  "Subtype_Selected",  "purity", 
                      "tumor_stage", "sample", "race", "ethnicity", "birth_year", "Protocol", 
                      "gender")
      # create PC plots for each of the defined confoundes (PC1 plotted against PC2-10)
      confounders.plots = lapply(confounders, function(variable){
        PCs.plots = lapply(c(2:10), function(i){
          pca_plot_batch(df = pca_results, colour_variable = variable, x_axis = "PC1", y_axis = i, 
                         title=paste0("PCA based on log2(exprs.) of ",ncol(exprs_log),
                                      " isomirs and ", nrow(exprs_log), " patients"), 
                         percentage=percentage, centroids = F)
        })
        return(PCs.plots)
      })
      
      # Scatterplots with centroids 
      confounders.centroids.plots = lapply(confounders, function(variable){
        PCs.plots = lapply(c(2:10), function(i){
          plot = pca_plot_batch(df = pca_results, colour_variable = variable, x_axis = "PC1", 
                                y_axis = i, title=paste0("PCA based on log2(exprs.) of ",
                                                         ncol(exprs_log)," isomirs and ", 
                                                         nrow(exprs_log), " patients"), 
                                percentage=percentage, centroids=T)
          return(plot)
        })
        return(PCs.plots)
      })

      #counfounders box plots 
      confounder.box.plots = lapply(confounders, function(variable){
        box.plots = lapply(c(1:10), function(i){
          y_name = paste0("PC",i)
          if(variable == "total.reads"){
            ggscatter(pca_results, x = "total.reads", y = y_name,
                      add = "reg.line",  # Add regressin line
                      add.params = list(color = "blue", fill = "lightgray"))+
              ylab(percentage[i])+
              stat_cor(method = "pearson", label.x = (max(pca_results$total.reads)-5000000), 
                       label.y = min(pca_results[,y_name]))
          }else {
            ggplot(data=pca_results, aes(x = reorder(pca_results[,variable], pca_results[,y_name], 
                                                     FUN = median), y = pca_results[,y_name]))+
              geom_boxplot()+
              theme_classic(base_size = 18)+
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
              labs(x = variable, y = percentage[i])
          }
        })
        return(box.plots)
      })
      
      pdf(paste0(figure.path, "/",conf$Project,"_PCA_boxplots_", conf$filter, conf$filter_value,
                 ".pdf"), width = 8, height= 6)
      for(k in 1:length(confounder.box.plots)){
        for(j in 1:length(confounder.box.plots[[k]])){
          print(confounder.box.plots[[k]][j])
        }
      }
      dev.off()
      
      # Save all PC plots (further analysis only with true confounders; rest to be ignored)
      pdf(paste0(figure.path, "/",conf$Project,"_PCA_scatterplots_", conf$filter, 
                 conf$filter_value,".pdf"), width = 8, height= 6)
      for(k in 1:length(confounders.plots)){
        for(j in 1:length(confounders.plots[[k]])){
          print(confounders.plots[[k]][j])
        }
      }
      dev.off()
      
      pdf(paste0(figure.path, "/",conf$Project,"_PCA_scatterplots_centroids_", conf$filter, 
                 conf$filter_value,".pdf"), width = 8, height= 6)
      for(k in 1:length(confounders.centroids.plots)){
        for(j in 1:length(confounders.centroids.plots[[k]])){
          print(confounders.centroids.plots[[k]][j])
        }
      }
      dev.off()
      
      
      # create table with p-values & each confounder
      # Plate batch effects if normal distributed: ANOVA, otherwise: Kruskal-Wallis
      # Platform batch effect if normal distributed: t-test, otherwise: Wilkoxon-rank-sum
      # Purity if normal distribution: logistic regression, otherwise Kendall tau rank correlation 
      # b (can handle ties contrary to Spearman)
      # total reads: spearman correlation
      
      # remove NA values for table
      for(header in names(pca_results)){
        pca_results[which(pca_results[,header] == "NA"), header] = NA
      }
      pca_results$purity = as.character(pca_results$purity)
      if(length(which(pca_results$purity == "0-10"))>0){
        pca_results$pur[pca_results$purity == "0-10"] = 0}
      if(length(which(pca_results$purity == "10-20"))>0){
        pca_results$pur[pca_results$purity == "10-20"] = 1}
      if(length(which(pca_results$purity == "20-30"))>0){
        pca_results$pur[pca_results$purity == "20-30"] = 2}
      if(length(which(pca_results$purity == "30-40"))>0){
        pca_results$pur[pca_results$purity == "30-40"] = 3}
      if(length(which(pca_results$purity == "40-50"))>0){
        pca_results$pur[pca_results$purity == "40-50"] = 4}
      if(length(which(pca_results$purity == "50-60"))>0){
        pca_results$pur[pca_results$purity == "50-60"] = 5}
      if(length(which(pca_results$purity == "60-70"))>0){
        pca_results$pur[pca_results$purity == "60-70"] = 6}
      if(length(which(pca_results$purity == "70-80"))>0){
        pca_results$pur[pca_results$purity == "70-80"] = 7}
      if(length(which(pca_results$purity == "80-90"))>0){
        pca_results$pur[pca_results$purity == "80-90"] = 8}
      if(length(which(pca_results$purity == "90-100"))>0){
        pca_results$pur[pca_results$purity == "90-100"] = 9}
      if(length(which(pca_results$purity == "100"))>0){
        pca_results$pur[pca_results$purity == "100"] = 10}
      pca_results$purity = pca_results$pur

      pval = lapply(1:10, function(pc){
        print(pc)
        PC.vector = pca_results[,paste0("PC", pc)]          
        plate.pval = kruskal.test(PC.vector ~ plate, data = pca_results)$p.value
        if(length(unique(na.omit(pca_results$Subtype_Selected)))>2){
          subtype.pval = kruskal.test(PC.vector ~ Subtype_Selected, data = pca_results)$p.value
        } else if(length(unique(na.omit(pca_results$Subtype_Selected)))==2){
          subtype.pval = wilcox.test(PC.vector ~ Subtype_Selected, data = pca_results)$p.value
        } else {
          subtype.pval = NA
        }
        if(length(unique(na.omit(pca_results$sample)))>2){
          sample.pval = kruskal.test(PC.vector ~ sample, data = pca_results)$p.value
        } else if(length(unique(na.omit(pca_results$sample)))==2){
          sample.pval = wilcox.test(PC.vector ~ sample, data = pca_results)$p.value
        } else {
          sample.pval = NA
        }
        if(length(unique(na.omit(pca_results$Platform)))>1){
          platform.pval = wilcox.test(PC.vector ~ Platform, data = pca_results)$p.value
        }else {
          platform.pval = NA
        }
        if(length(unique(na.omit(pca_results$gender)))>1){
          gender.pval = wilcox.test(PC.vector ~ gender, data = pca_results)$p.value
        } else {
          gender.pval = NA
        }
        if(length(unique(na.omit(pca_results$Protocol)))>1){
          protocol.pval = wilcox.test(PC.vector ~ Protocol, data = pca_results)$p.value
        } else {
          protocol.pval = NA
        }
        if(!unique(is.na(pca_results$purity))){
          purity2.pval = summary(glm(as.factor(pca_results$purity)~PC.vector, data = pca_results, 
                                     family = binomial(link = "logit")))$coefficients[2,4]
          purity.cor=cor.test(as.numeric(pca_results$purity),PC.vector,method="kendall")$estimate
          purity.pval=cor.test(as.numeric(pca_results$purity),PC.vector,method="kendall")$p.value
        } else {
          purity2.pval = NA
          purity.cor = NA
          purity.pval = NA
        }
        reads.cor = cor.test(pca_results$total.reads, PC.vector, method = "spearman")$estimate
        reads.pval = cor.test(pca_results$total.reads, PC.vector, method = "spearman")$p.value
        if(length(unique(na.omit(pca_results$tumor_stage)))>2){
          tumor_stage.pval = kruskal.test(PC.vector ~ tumor_stage, data = pca_results)$p.value
        } else if(length(unique(na.omit(pca_results$tumor_stage)))==2){
          tumor_stage.pval = wilcox.test(PC.vector ~ tumor_stage, data = pca_results)$p.value
        } else {
          tumor_stage.pval = NA
        }  
        if(length(unique(na.omit(pca_results$race)))>2){
          race.pval = kruskal.test(PC.vector ~ race, data = pca_results)$p.value
        } else if(length(unique(na.omit(pca_results$race)))==2){
          race.pval = wilcox.test(PC.vector ~ race, data = pca_results)$p.value
        } else {
          race.pval = NA
        } 
        
        
        
        pvals = data.frame(variable = c("plate", "platform", "purity", "purity.2", "reads", 
                                "subtype", "gender", "protocol", "sample", "tumor_stage", "race"),
                           p.value = c(plate.pval, platform.pval, purity.pval, purity2.pval, 
                                       reads.pval, subtype.pval, gender.pval, protocol.pval, 
                                       sample.pval, tumor_stage.pval, race.pval),
                           correlation = c(NA, NA, purity.cor, NA, reads.cor, NA , NA , NA, NA, 
                                           NA, NA),
                           component = rep(paste0("PC", pc), 11),
                           test = c("kruskal.wallis", "wilcox.test", "kendall.tau.cor",
                                    "logistic.regression", "spearman.cor", "kruskal.wilcox.test", 
                                    "wilcox.test", "wilcox.test", "kruskal.wilcox.test", 
                                    "kruskal.wilcox.test", "kruskal.wilcox.test"))
        return(pvals)
      })
      
      pval.df = do.call(rbind, pval)
      # adjust p values (for multiple testing)
      pval.df$fdr = p.adjust(pval.df$p.value, method = "fdr")
      # change p values and correlation formating: 2 digits, with scientific writing (e-X)
      
      pval.signif = pval.df[which(pval.df$fdr < 0.01),]
      pval.signif = pval.signif[-which(pval.signif$correlation < 0.25 & 
                                         pval.signif$correlation > -0.25),]
      pval.signif = pval.signif[order(pval.signif$fdr, decreasing = FALSE),]
      
      if(nrow(pval.signif) > 0){
        pval.signif$fdr.format = formatC(pval.signif$fdr, format = "e", digits = 2)
        pval.signif$p.value.format = formatC(pval.signif$p.value, format = "e", digits = 2)
        pval.signif$correlation.format = formatC(pval.signif$correlation, format = "e", digits =2)
        write.csv(pval.signif, paste0(figure.path, "/",conf$Project,"_signif.p.values.pcs_", 
                                      conf$filter, conf$filter_value,".csv"), row.names = FALSE)
        
      }
      
      pval.df$fdr.format = formatC(pval.df$fdr, format = "e", digits = 2)
      pval.df$p.value.format = formatC(pval.df$p.value, format = "e", digits = 2)
      pval.df$correlation.format = formatC(pval.df$correlation, format = "e", digits = 2)
      write.csv(pval.df, paste0(figure.path, "/",conf$Project,"_p.values.pcs_", conf$filter, 
                                conf$filter_value,".csv"), row.names = FALSE)
      
    }, .parallel = TRUE)  
    }
  }, .parallel = TRUE)    
})
