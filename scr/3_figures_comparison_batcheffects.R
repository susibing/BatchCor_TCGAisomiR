#!/usr/bin/env Rscript
library(optparse)

#### CHECK MIRNA CORRECTION!!!

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'ggplot2_3.3.2', 'optparse_1.6.6', 'data.table_1.13.0')
Description = 'This script is to produce Figure 2c / 2f and Supplementary figures S2, S5, S7-S21.
This script aims to visualize a comparison of quantified batch effects in:
before vs. after batch correction; mirna vs. isomiR data set; combat vs. limma correction method
Based on tables including statistical association between PCs and potential confounders, generated 
and stored in 3_assess_batch_effects.R.'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## LIBRARIES ---- 
library(ggplot2)
library(data.table)

## CONFIG PART ----
## Further parameter setting
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$output_dir = file.path(conf$base_dir, "results/pVal_forest_plots")

source(file.path(conf$base_dir, "scr", "functions.batch_cor.R")) # load functions 

## Correction parameters: 
## Optimal combination of batch correction for the cohorts
conf.batch_cor = fread(file.path(conf$base_dir, "data", "config.batch_cor.txt"))
## only tumor samples in scope
conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == FALSE,] 

## Function to load tables from 3_assess_batch_effects.R 
load_signif = function(path, project, data_type, include.all.samples, normalization_type, filter, 
                       conf.batch_cor, i, correction_method){
  if(data_type == "before_batch_cor"){
    file = file.path(path, data_type, project, paste0("include.all.samples", include.all.samples), 
                     paste0(normalization_type, "_", filter, "_", conf.batch_cor$platform[i]), 
                     paste0(project, "_p.values.pcs_", filter, ".csv"))
  } else if(data_type == "batch_cor"){
    if(correction_method == "conf.batch"){
      correction_method = conf.batch_cor$correction_method[i]
    }
    file = file.path(path, data_type, correction_method, project, 
                     paste0("include.all.samples", include.all.samples), 
                     paste0(normalization_type, "_", filter, "_", conf.batch_cor$var_of_interest[i], "_", 
                           conf.batch_cor$batch_1[i], "_", conf.batch_cor$batch_2[i], "_", 
                           conf.batch_cor$platform[i], "_", "rm", conf.batch_cor$rm_batches[i], "_"), 
                     paste0(project, "_p.values.pcs_", filter, ".csv"))
  } else if (data_type == "miRNA_only/before_batch_cor"){
    file = file.path(path, data_type, project, paste0("include.all.samples", include.all.samples), 
                     paste0(normalization_type, "_", filter, "_", conf.batch_cor$platform[i]), 
                     paste0(project, "_p.values.pcs_", filter, ".csv"))
  }
  print(file)
  pval.df = fread(file, data.table = F, stringsAsFactors = F)
  
  return(pval.df)
}

## Function for facet plot including all confounding variables (supplement)
forest.plot.all.conf = function(color.palette, confounders, df, y_val, col.values, line.val){
  df = df[which(df$variable %in% confounders),]
  df$correlation = 1
  
  # New facet label names for confounders variable
  if("subtype" %in% confounders){
    confounders.labs <- c("Gender", "Plate", "Platform", "Tumor Purity", "Number of Mapped Reads", "Tumor Subtype", "Tumor Stage")
    names(confounders.labs) <- c("gender", "plate", "platform", "purity", "reads", "subtype", "tumor_stage")
  } else {
    confounders.labs <- c("Gender", "Plate", "Platform", "Tumor Purity", "Number of Mapped Reads", "Tumor Stage")
    names(confounders.labs) <- c("gender", "plate", "platform", "purity", "reads", "tumor_stage")
  }
  
  plot = ggplot(data=df, aes(x = component,y = df[,y_val], ymin = df[,y_val], ymax = df[,y_val]))+
    geom_hline(yintercept = line.val, linetype=2)+
    geom_pointrange(aes(col=time.point), alpha = 0.8)+
    labs(x = "", y = "Association with variable (log10 q-value)", col = "Time point")+
    coord_flip()+
    scale_colour_manual(values= color.palette) +
    theme_classic2(base_size = 16)+
    facet_wrap(~variable, labeller = labeller(variable = confounders.labs)) +
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.x=element_text(face="bold", angle = 0, hjust = 0),
          axis.title=element_text(size=16,face="bold"),
          axis.text.y = element_text(hjust=0,vjust = 0,angle=0,face="bold"))
  
  return(plot)
}

## Function for individual plots for each potential confounding variable
forest.plot.individual = function(color.palette, df, variable, y_val, col.values, line.val){
  df = df[df$variable == variable,]
  df$correlation = 1
  plot = ggplot(data=df, aes(x = component,y = df[,y_val], ymin = df[,y_val], ymax = df[,y_val]))+
    geom_hline(yintercept = line.val, linetype=2)+
    geom_pointrange(aes(col=time.point), size = 1.5, alpha = 0.8, position = position_dodge(width = 0.25))+
    labs(x = "", y = paste0("Association with ", variable, " variable (log10 q-value)"), col = "Time point")+
    coord_flip()+
    scale_colour_manual(values= color.palette) +
    theme_classic2(base_size = 18)+    
    theme(legend.position = "bottom")
  return(plot)
}

## Function to save plots at defined path
forest.plots.save = function(color.palette, normalization_type, filter, conf.batch_cor, save.path, 
                             title, data_type1, data_type2, path, correction_method){
  for(i in 1:nrow(conf.batch_cor)){ 
    project = conf.batch_cor$Project[i]
    if(data_type1 == "before_batch_cor" & data_type2 == "batch_cor"){
      p.vals.before = load_signif(path = path, project = project, data_type = data_type1, 
                               include.all.samples = FALSE, normalization_type = normalization_type, 
                               filter = filter, conf.batch_cor = conf.batch_cor, i = i)
      p.vals.before$time.point = "before"
      p.vals.after = load_signif(path = path, project = project, data_type = data_type2, 
                               include.all.samples = FALSE, normalization_type = normalization_type, 
                               filter = filter, conf.batch_cor = conf.batch_cor, i = i, 
                               correction_method = correction_method)
      p.vals.after$time.point = "after"
    } else if(data_type1 == "before_batch_cor" & data_type2 == "miRNA_only/before_batch_cor"){
      p.vals.before = load_signif(path = path, project = project, data_type = data_type1, 
                               include.all.samples = FALSE, normalization_type = normalization_type, 
                               filter = filter, conf.batch_cor = conf.batch_cor, i = i)
      p.vals.before$time.point = "isomiR"
      p.vals.after = load_signif(path = path, project = project, data_type = data_type2, 
                               include.all.samples = FALSE, normalization_type = normalization_type, 
                               filter = filter, conf.batch_cor = conf.batch_cor, i = i)
      p.vals.after$time.point = "miRNA"
      p.vals.after$percentage = NA
    } else if(data_type1 == "batch_cor" & data_type2 == "batch_cor"){
      p.vals.before = load_signif(path = path, project = project, data_type = data_type2, 
                               include.all.samples = FALSE, normalization_type = normalization_type, 
                               filter = filter, conf.batch_cor = conf.batch_cor, i = i, 
                               correction_method = "combat")
      p.vals.before$time.point = "combat"
      p.vals.after = load_signif(path = path, project = project, data_type = data_type2, 
                               include.all.samples = FALSE, normalization_type = normalization_type, 
                               filter = filter, conf.batch_cor = conf.batch_cor, i = i, 
                               correction_method = "limma")
      p.vals.after$time.point = "limma"
    }
    
    forest.df = rbind(p.vals.before, p.vals.after)
    forest.df = forest.df[-which(is.na(forest.df$p.value)),]
    forest.df$log10.fdr = log10(forest.df$fdr + (min(forest.df$fdr[forest.df$fdr != 0]))/10)
    forest.df$component <- factor(forest.df$component,levels = c("PC1", "PC2", "PC3", "PC4", "PC5", 
                                                                 "PC6","PC7", "PC8", "PC9", "PC10"))
    forest.df$variable = as.character(forest.df$variable)

    variables = unique(forest.df$variable) 
    
    forest.plot_all.log = forest.plot.all.conf(color.palette = color.palette, df = forest.df, 
                                               confounders = variables, y_val = "log10.fdr", 
                                               col.values = col.values, line.val = log10(0.01))

    forest.plots.log = lapply(variables, function(variable){
      plot.v2 = forest.plot.individual(color.palette = color.palette, df = forest.df, 
                                       variable = variable, y_val = "log10.fdr", 
                                       col.values = col.values, line.val = log10(0.01))
      return(plot.v2)
    })
    
    ## define size of image depending on number of potential confounders examined and save facet plots
    if(length(forest.plots.log) %in% c(5,6)){
      pdf(paste0(save.path, "/",project, "_", title,"_facet.pdf"), width = 9, height= 6)
        print(forest.plot_all.log)
      dev.off()
    } else if(length(forest.plots.log) %in% c(7,8,9)){
      pdf(paste0(save.path, "/",project, "_", title,"_facet.pdf"), width = 9.5, height= 8.2)
        print(forest.plot_all.log)
      dev.off()
    } else if(length(forest.plots.log) %in% c(4)){
      pdf(paste0(save.path, "/",project, "_", title,"_facet.pdf"), width = 7, height= 6)
      print(forest.plot_all.log)
      dev.off()
    } else if(length(forest.plots.log) %in% c(3)){
      pdf(paste0(save.path, "/",project, "_", title,"_facet.pdf"), width = 9.5, height= 4)
      print(forest.plot_all.log)
      dev.off()
    }
  
    ## Save individual plots 
    pdf(paste0(save.path, "/",project, "_", title,"_singular.pdf"), width = 6.5, height= 6.5)
    for(k in 1:length(forest.plots.log)){print(forest.plots.log[[k]])}
    dev.off()
  }
  
}

## Forest plots before and after batch correction (chosen methods above)
output = file.path(conf$output_dir, "before_after_sum0")
if(!dir.exists(output)){dir.create(output, recursive = T)}
forest.plots.save(normalization_type = "rpm", 
                  path = file.path(conf$base_dir, "results"),
                  filter = "sum0",
                  data_type1 = "before_batch_cor",
                  data_type2 = "batch_cor",
                  conf.batch_cor = conf.batch_cor,
                  title = "pVal_forest_plts_before_after_cor",
                  save.path = output, 
                  correction_method = "conf.batch",
                  color.palette = c("#1b9e77", "#d95f02"))

output = file.path(conf$output_dir, "before_after_median15")
if(!dir.exists(output)){dir.create(output, recursive = T)}
forest.plots.save(normalization_type = "rpm", 
                  path = file.path(conf$base_dir, "results"),
                  filter = "median15",
                  data_type1 = "before_batch_cor",
                  data_type2 = "batch_cor",
                  conf.batch_cor = conf.batch_cor,
                  title = "pVal_forest_plts_before_after_cor",
                  save.path = output, 
                  correction_method = "conf.batch",
                  color.palette = c("#1b9e77", "#d95f02"))



## Forest plots miRNA and isomiR before batch correction
output = file.path(conf$output_dir, "miRNA_isomiR_sum0")
if(!dir.exists(output)){dir.create(output, recursive = T)}
forest.plots.save(normalization_type = "rpm", 
                  path = file.path(conf$base_dir, "results"),
                  filter = "sum0",
                  data_type1 = "before_batch_cor",
                  data_type2 = "miRNA_only/before_batch_cor",
                  conf.batch_cor = conf.batch_cor,
                  title = "pVal_forest_plts_miRNA_isomiR_cor",
                  save.path = output, 
                  correction_method = "conf.batch",
                  color.palette = c("#7b3294", "#008837"))

output = file.path(conf$output_dir, "miRNA_isomiR_median15")
if(!dir.exists(output)){dir.create(output, recursive = T)}
forest.plots.save(normalization_type = "rpm", 
                  path = file.path(conf$base_dir, "results"),
                  filter = "median15",
                  data_type1 = "before_batch_cor",
                  data_type2 = "miRNA_only/before_batch_cor",
                  conf.batch_cor = conf.batch_cor,
                  title = "pVal_forest_plts_miRNA_isomiR_cor",
                  save.path = output, 
                  correction_method = "conf.batch",
                  color.palette = c("#7b3294", "#008837"))


## Comparison combat - limma
output = file.path(conf$output_dir, "combat_limma_sum0")
if(!dir.exists(output)){dir.create(output, recursive = T)}
forest.plots.save(normalization_type = "rpm", 
                  path = file.path(conf$base_dir, "results"),
                  filter = "sum0",
                  data_type1 = "batch_cor",
                  data_type2 = "batch_cor",
                  conf.batch_cor = conf.batch_cor,
                  title = "pVal_forest_plts_combat_limma_cor",
                  save.path = output, 
                  correction_method = "test", 
                  color.palette = c("#d7191c", "#2c7bb6"))

output = file.path(conf$output_dir, "combat_limma_median15")
if(!dir.exists(output)){dir.create(output, recursive = T)}
forest.plots.save(normalization_type = "rpm", 
                  path = file.path(conf$base_dir, "results"),
                  filter = "median15",
                  data_type1 = "batch_cor",
                  data_type2 = "batch_cor",
                  conf.batch_cor = conf.batch_cor,
                  title = "pVal_forest_plts_combat_limma_cor",
                  save.path = output, 
                  correction_method = "test",
                  color.palette = c("#d7191c", "#2c7bb6"))
