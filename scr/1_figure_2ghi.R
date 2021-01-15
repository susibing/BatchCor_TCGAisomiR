#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'tidyverse_1.3.0', 'matrixStats_0.57.0', 
'data.table_1.13.0', 'wrapr_2.0.2', 'ggpubr_0.4.0', 'doParallel_1.0.15')
Description = 'Visualization of sequencing depth per plate, colored by platform. (Boxplot figure 2g)
Annotation corrected isomiR expression quantification data serves as input.'
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
conf$include_all_samples = F # if normal samples should be included in batch correction, set to T
conf$reads_filter = FALSE # reads filter is set to FALSE as all samples should be included for this overview
conf$tot.reads.filter = 1000000 # threshold reads filter

data_type = "before_batch_cor" # define folder to put figure

## EXECUTION

## Submit job to LSF cluster
bsub_chunk(name = "figure_2g",
           variables = c("conf.batch_cor", "conf", "data_type"), 
           memory = 30, hour = 24, core = 8, 
           packages = c("tidyverse", "data.table", "matrixStats"),{
             
  library(ggpubr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(wrapr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/') 
  library(doParallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  
  ## load functions 
  source(file.path(conf$base_dir, "scr", "functions.batch_cor.R")) # load functions
  registerDoParallel(cores=8) 
  
  ## Correction parameters: 
  ## Only required to identify cohorts and folders in scope
  conf.batch_cor = fread(file.path(conf$base_dir, "data", "config.batch_cor.txt"))
  conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == FALSE,]        
       
  ## go through different configurations (in parallel)
  plyr::llply(1:nrow(conf.batch_cor), function(i){
    
    ## Load data
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
    folder_name = file.path(conf$base_dir, "results", data_type, conf$Project, 
                            paste0("include.all.samples", conf$include_all_samples), 
                            paste0(conf$normalization_type, "_", conf$filter, conf$filter_value,"_", 
                                   conf$platform))
    #folder_name = paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", data_type)
    if(!dir.exists(folder_name)){dir.create(folder_name, recursive = T)}
    
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
    
    ## figure 2g ----
    boxplot = ggplot(data=confounder.df, aes(x=reorder(confounder.df[,"plate"], 
                                                       confounder.df[,"total.reads"], FUN = median), 
                                             y=total.reads))+
      geom_boxplot(aes(fill=Platform))+
      theme_classic2(base_size = 18)+
      scale_fill_manual("Platform",
                        values = c("#4575b4", "#d73027", "#fc8d59"),
                        breaks = c("HiSeq", "GA", "NA"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
            legend.position = "bottom",
            axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
      labs(x = "Plate", y = "Number of total mapped reads")
    
    pdf(file.path(folder_name, "boxplot_seqdepth_plate.pdf"), width = 6.5, height = 7.15)
    plot(boxplot)
    dev.off()
    
    ## figure 2h ----
    quantile_df = confounder.df[,c("ids", "quantile", "total.reads")]
    quantile_df$nr.isomiRs = sapply(quantile_df$ids, function(pid){
      length(which(exprs_df[,pid, drop=FALSE] > 0))
    })
    
    median_exprs_df = as.data.frame(t(exprs_df))
    median_exprs_df[median_exprs_df == 0] <- NA
    median_exprs_df$median = rowMedians(as.matrix(median_exprs_df), na.rm = T)
    median_exprs_df$ids = rownames(median_exprs_df)
    
    quantile_df = merge(x = median_exprs_df[, c("ids", "median")], y = quantile_df, by = "ids")
    
    scatter.read.cor = ggscatter(quantile_df, x = "total.reads", y = "nr.isomiRs",
              color = "median", size = 3,
              add = "reg.line",  # Add regressin line
              add.params = list(color = "red", fill = "lightgray"),
              font.label = c(16, "italic"))+
      labs(x = "Number of total mapped reads", y = "Nr. of isomiRs with mapped reads", 
           col = "Median exprs. [RPM]")+
      stat_cor(method = "spearman")+
      theme_classic2(base_size = 18)+
      theme(legend.position = "bottom",
            axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
      scale_colour_gradientn(colors=rainbow(12)) 
    
    pdf(file.path(folder_name, "scatter_read.cor.pdf"), width = 6.5, height = 7.15)
    plot(scatter.read.cor)
    dev.off()
    
    ## figure 2i ----
    data = exprs_df
    data$isomiRs = rownames(data)
    
    data$ThreePrimeEnd = sapply(data$isomiRs, function(x){
      as.numeric(strsplit(x, split = "\\|")[[1]][3])})
    data$FivePrimeEnd = sapply(data$isomiRs, function(x){
      as.numeric(strsplit(x, split = "\\|")[[1]][2])})
    data$isomiR_status = sapply(1:nrow(data), function(x){
      print(x)
      if((data$ThreePrimeEnd[x] == 0) & (data$FivePrimeEnd[x] == 0)){
        isomiR_status = "Canonical"
      } else if((data$ThreePrimeEnd[x] != 0) & (data$FivePrimeEnd[x] == 0)){
        isomiR_status = "3' isomiR"
      } else if((data$ThreePrimeEnd[x] == 0) & (data$FivePrimeEnd[x] != 0)){
        isomiR_status = "5' isomiR"
      } else if((data$ThreePrimeEnd[x] != 0) & (data$FivePrimeEnd[x] != 0)){
        isomiR_status = "Combi isomiR"
      } 
      return(isomiR_status)
    })
    
    ## Summarize data per sample 
    pid_sum = data %>% gather(key = "pid", value = "exprs", 
                              -c(isomiRs, ThreePrimeEnd, FivePrimeEnd, isomiR_status))
    pid_summary = pid_sum[pid_sum$exprs != 0,] %>% group_by(pid, isomiR_status) %>% 
      tally() # pid_sum$exprs != 0 sum filter = 0
    
    pid_summary = pid_summary %>% spread(key = isomiR_status, value = n)
    pid_summary$Total = rowSums(pid_summary[,c("Canonical", "3' isomiR", "5' isomiR", 
                                               "Combi isomiR")])
    pid_summary = pid_summary %>% gather(key = isomiR_status, value = n, -pid)
    
    ## number of possible isomiRs per category
    all_canonical = length(data$isomiR_status[data$isomiR_status == "Canonical"])
    all_3prime = length(data$isomiR_status[data$isomiR_status == "3' isomiR"])
    all_5prime = length(data$isomiR_status[data$isomiR_status == "5' isomiR"])
    all_mixed = length(data$isomiR_status[data$isomiR_status == "Combi isomiR"])
    all_total = all_canonical + all_3prime + all_5prime + all_mixed
    
    ## number of isomiRs mapped per category with rowsum > 0
    data_expressed = data[which(rowSums(as.matrix(data[,1:(ncol(data)-4)]))> 0),]
    expr_canonical = length(data_expressed$isomiR_status[data_expressed$isomiR_status == "Canonical"])
    expr_3prime = length(data_expressed$isomiR_status[data_expressed$isomiR_status == "3' isomiR"])
    expr_5prime = length(data_expressed$isomiR_status[data_expressed$isomiR_status == "5' isomiR"])
    expr_mixed = length(data_expressed$isomiR_status[data_expressed$isomiR_status == "Combi isomiR"])
    expr_total = expr_canonical + expr_3prime + expr_5prime + expr_mixed
    
    pid_summary_2 = pid_summary
    pid_summary_2$isomiR_status = "Total"
    pid_summary = rbind(pid_summary, pid_summary_2)
    
    ## relative freq => per all possible isomiRs
    pid_summary$FreqRel = NA
    pid_summary$FreqRel[pid_summary$isomiR_status == "Canonical"] = 
      pid_summary$n[pid_summary$isomiR_status == "Canonical"] / all_canonical
    pid_summary$FreqRel[pid_summary$isomiR_status == "3' isomiR"] = 
      pid_summary$n[pid_summary$isomiR_status == "3' isomiR"] / all_3prime
    pid_summary$FreqRel[pid_summary$isomiR_status == "5' isomiR"] = 
      pid_summary$n[pid_summary$isomiR_status == "5' isomiR"] / all_5prime
    pid_summary$FreqRel[pid_summary$isomiR_status == "Combi isomiR"] = 
      pid_summary$n[pid_summary$isomiR_status == "Combi isomiR"] / all_mixed
    pid_summary$FreqRel[pid_summary$isomiR_status == "Total"] = 
      pid_summary$n[pid_summary$isomiR_status == "Total"] / all_total
    pid_summary$isomiR_status = factor(pid_summary$isomiR_status,
                                       levels=c("Canonical","3' isomiR","5' isomiR","Combi isomiR", 
                                                "Total"))
    
    ## rel freq per all mapped isomiRs 
    pid_summary$FreqRelExpr = NA
    pid_summary$FreqRelExpr[pid_summary$isomiR_status == "Canonical"] = 
      pid_summary$n[pid_summary$isomiR_status == "Canonical"] / expr_canonical
    pid_summary$FreqRelExpr[pid_summary$isomiR_status == "3' isomiR"] = 
      pid_summary$n[pid_summary$isomiR_status == "3' isomiR"] / expr_3prime
    pid_summary$FreqRelExpr[pid_summary$isomiR_status == "5' isomiR"] = 
      pid_summary$n[pid_summary$isomiR_status == "5' isomiR"] / expr_5prime
    pid_summary$FreqRelExpr[pid_summary$isomiR_status == "Combi isomiR"] = 
      pid_summary$n[pid_summary$isomiR_status == "Combi isomiR"] / expr_mixed
    pid_summary$FreqRelExpr[pid_summary$isomiR_status == "Total"] = 
      pid_summary$n[pid_summary$isomiR_status == "Total"] / expr_total

    ## merge with total reads info 
    pid_summary_merge = merge(x=pid_summary, 
                              y = confounder.df[,c("total.reads", "quantile", "ids")], 
                              by.x = "pid", by.y = "ids")

    ## Create plots
    plot_rel = ggplot(data = pid_summary_merge, 
                      aes(x = isomiR_status, y = FreqRelExpr, fill = quantile))+
      geom_boxplot()+
      scale_fill_manual(values = c("#1f78b4",	"#a6cee3",	"#b2df8a",	"#33a02c"))+
      theme_classic2(base_size = 18)+
      labs(x = "IsomiR Type", y = paste0("Relative isomiR detection per patient"), 
           fill = "Quartile")+
      theme(legend.position = "bottom",
            axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
    
    plot_abs = ggplot(data = pid_summary_merge[pid_summary_merge$isomiR_status != "Total",], 
                      aes(x = isomiR_status, y = n, fill = quantile))+
      geom_boxplot()+
      scale_fill_manual(values = c("#1f78b4",	"#a6cee3",	"#b2df8a",	"#33a02c"))+
      theme_classic2(base_size = 18)+
      labs(x = "IsomiR Type", y = paste0("Absolute isomiR detection per patient"), 
           fill = "Quartile")+
      theme(legend.position = "bottom",
            axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
    
    pdf(file.path(folder_name, "boxplot_isomiRcat.pdf"), width = 6.5, height = 6.8)
    plot(plot_rel)
    plot(plot_abs)
    dev.off()
    
  }, .parallel = TRUE)
})

