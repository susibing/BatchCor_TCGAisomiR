#!/usr/bin/env Rscript
library(optparse)

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'optparse_1.6.6')
Description = 'This script includes all functions relevant for batch correction and visualization of 
batch effects in the TCGA IsomiR Quantificaion files. Structured depending on in which script they 
are called'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## general functions ----

#' negation of the %in% operator
#' \code{'%!in%'} returns boolean value (TRUE/FALSE) depending on whether input x is not element of 
#' input y
#' @param x: element (character, integer, ...)
#' @param y: vector or data frame 
#' @return boolean value (TRUE/FALSE) 
#' \dontrun{
#' x %!in% y
#' }
'%!in%' <- function(x,y){!('%in%'(x,y))}

#' make rownames out of first column 
#' \code{makeRn} returns a data frame where the first column is converted into rownames
#' @param x: input for this function is a data frame
#' @return If the first column does not contain duplicate values, it will be converted into rownames 
#'         The first column will then be deleted. Otherwise, an error message will be returned. 
#' \dontrun{
#' makeRn(df)
#' }
makeRn <- function(x) {
      rownames(x) <- x[, 1]
      x[, 1] <- NULL
      return(x)
    }    

#' make first column out of rownames
#' \code{make1col} returns a data frame where rownames are converted into the first column
#' @param x: input for this function is a data frame
#' @return Rownames will be converted into the first column of a dataframe and rownames are then 
#' deleted. The first column is named "ids".
#' \dontrun{
#' make1col(df)
#' }
make1col <- function(x) {
  x <- cbind(rownames(x), data.frame(x, row.names = NULL))
  colnames(x)[1] <- "ids"
  return(x)
}

#' make colors transparent
#' \code{makeTransparent} makes a colour transparent 
#' @param someColor: color name, e.g."red"
#' @return name of the transparent version of the colour
#' \dontrun{
#' makeTransparent("red")
#' }
makeTransparent <- function(someColor, alpha=100) #note: always pass alpha on the 0-255 scale
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}  

#' read csv file and have first column as rownames
#' \code{load_data} 
#' @param path: path to file to csv file to read
#' @return data frame with first column as row names
#' \dontrun{
#' load_data(path)
#' }
load_data = function(path){
  exprs = data.table::fread(path, sep=",", check.names = F, stringsAsFactors = F, header = T, 
                            data.table = F)
  rownames(exprs) = exprs[,1]
  exprs[,1] = NULL
  return(exprs)
}

#' Adaption of the theme_classic theme from ggplot2
#' \code{theme_classic2} 
#' @param base_size defines font size in ggplot; added to ggplot 
#' @return adapted classic theme in plot
#' @expample ggplot + theme_classic2(base_size = 18) 
#' \dontrun{
#' theme_classic2(base_size)
#' }
theme_classic2 = function(base_size){
  theme_bw(base_size = base_size, base_family = "", 
           base_line_size = 18/22, base_rect_size = 18/22) %+replace% 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", 
                                                                       size = rel(1)), legend.key = element_blank(), 
          strip.background = element_rect(fill = "white", colour = "black", 
                                          size = rel(2)), complete = TRUE)
}


## 1_data_download.R script ----

#' Date pre-processing of TCGA isomIiR/miRNA expression data downloaded with the TCGA Biolinks R 
#' package: Creation of the count matrix
#' \code{create_count_matrix} 
#' @param df: input data frame 
#' @return count matrix to save in /data/DataClean 
#' \dontrun{
#' create_count_matrix(df)
#' }
create_count_matrix = function(df){
  ## create count matrix
  count.matrix = df[,c("pid", "read_count", "Name", "chrom", "start", "end", "strand")] %>% 
    spread(key = pid, value = read_count)
  count.matrix[is.na(count.matrix)] = 0
  ## sum up same entries (isomiRs derived from different genomic locations) 
  ## => counts were devided by the number of loci
  count.matrix = data.table(count.matrix)[, lapply(.SD[, names(count.matrix)[6:ncol(count.matrix)], 
                                                       with = FALSE], sum), by = Name] 
  rownames(count.matrix) = count.matrix$Name
  count.matrix$Name = NULL
  return(count.matrix)
}

#' Date pre-processing of TCGA isomIiR/miRNA expression data downloaded with the TCGA Biolinks R 
#' package: Creation of the rpm matrix
#' \code{create_rpm_matrix} 
#' @param df: input data frame 
#' @return rpm matrix to save in /data/DataClean
#' \dontrun{
#' create_rpm_matrix(df)
#' }
create_rpm_matrix = function(df){
  ## create rpm matrix (provided by TCGA)
  rpm.matrix = df[,c("pid", "reads_per_million_miRNA_mapped", "Name", "chrom", "start", "end", "strand")] %>% spread(key = pid, value = reads_per_million_miRNA_mapped)
  rpm.matrix[is.na(rpm.matrix)] = 0
  rpm.matrix = data.table(rpm.matrix)[, lapply(.SD[, names(rpm.matrix)[6:ncol(rpm.matrix)], with = FALSE], sum), by = Name] # faster
  rownames(rpm.matrix) = rpm.matrix$Name
  rpm.matrix$Name = NULL
  return(rpm.matrix)
}

## 1_figure_1c.R

#' Read isomiR expression data and add columns (isomiR types, median expression, etc.)
#' \code{load_annotation_data} 
#' @param path: path to isomiR expression data to read
#' @return isomiR expression data with added column information
#' \dontrun{
#' load_annotation_data(path)
#' }
load_annotation_data = function(path){
  exprs = data.table::fread(path, sep=",", check.names = F, stringsAsFactors = F, header = T, 
                            data.table = F)
  rownames(exprs) = exprs[,1]
  exprs[,1] = NULL
  # add columns with expression and isomiR details
  exprs$median = rowMedians(as.matrix(exprs))
  exprs$miRNA.arm = sapply(rownames(exprs), function(x){strsplit(x, split = "\\|")[[1]][1]})
  exprs$isoform = sapply(rownames(exprs), function(x){paste0(strsplit(x, split = "\\|")[[1]][2],
                                                             strsplit(x, split = "\\|")[[1]][3])})
  exprs$three.prime= sapply(rownames(exprs), function(x){strsplit(x, split = "\\|")[[1]][3]})
  exprs$five.prime= sapply(rownames(exprs), function(x){strsplit(x, split = "\\|")[[1]][2]})
  return(exprs)
}

## 2_remove_batch_effects.R / 3_assess_batch_effects.R ----

#' Prepare confounder matrix (depending on expression data) 
#' \code{prep_conf_df} 
#' @param exprs_log: expression data frame of one TCGA project, used to extract pids
#' @param confounders: TCGA confounders file read from /data folder; includes information for all
#' samples from all TCGA projects under analysis 
#' @param conf: configuration list, includes information about further filtering of the dataset
#' @return filtered and pre-processed confounders file
#' \dontrun{
#' prep_conf_df(exprs_log, confounders, conf)
#' }
prep_conf_df = function(exprs_log, confounders, conf){
  ## load confounders matrix
  confounder.df = confounders[which(confounders$ids %in% rownames(exprs_log)),
                              c("ids", "purity", "plate", "race", "ethnicity", "Subtype_Selected", 
                                "Platform", "gender", "tumor_stage", "birth_year", "Protocol", 
                                "total.reads")]
  confounder.df$sample.id=sapply(confounder.df$ids, function(x){
    paste("TCGA",strsplit(x, split="-")[[1]][2],strsplit(x, split="-")[[1]][3],sep="-")
  })
  ## if by any chance some IDs are not present in the confounder matrix 
  ## => add those rows for coloring purposes 
  if(length(which(rownames(exprs_log) %!in% confounder.df$ids)) > 0){
    ids = rownames(exprs_log)[which(rownames(exprs_log) %!in% confounder.df$ids)]
    additional.pids = cbind.data.frame(ids, purity=NA, plate=NA, race=NA, ethnicity=NA, 
                                  Subtype_Selected=NA, Platform=NA, gender=NA, stringsAsFactors = F)
    additional.pids$sample.id=sapply(additional.pids$ids, function(x){
      paste("TCGA",strsplit(x, split="-")[[1]][2],strsplit(x, split="-")[[1]][3],sep="-")
    })
    additional.pids$plate = sapply(additional.pids$ids, function(x){strsplit(x,split="-")[[1]][6]})
    additional.pids$gender = sapply(additional.pids$sample.id, function(x){
      confounder.df$gender[which(confounder.df$sample.id == x)][[1]]
    })
    additional.pids$tumor_stage = sapply(additional.pids$sample.id, function(x){
      confounder.df$tumor_stage[which(confounder.df$sample.id == x)][[1]]
    })
    additional.pids$Protocol = sapply(additional.pids$sample.id, function(x){
      confounder.df$Protocol[which(confounder.df$sample.id == x)][[1]]
    })
    additional.pids$total.reads = sapply(additional.pids$sample.id, function(x){
      confounder.df$total.reads[which(confounder.df$sample.id == x)][[1]]
    })
    confounder.df = rbind(confounder.df, additional.pids)
  }
  confounder.df$sample = sapply(confounder.df$ids, function(x){strsplit(x, split="-")[[1]][4]})
  confounder.df$sample = substr(confounder.df$sample,1,2)
  ## exchange sample number with written explanation
  confounder.df$sample[confounder.df$sample == "01"] = "Tumor" #"PrimarySolidTumor"
  confounder.df$sample[confounder.df$sample == "02"] = "Tumor" #"RecurrentSolidTumor"
  confounder.df$sample[confounder.df$sample == "05"] = "Tumor" #"AdditionalNewPrimaryTumor"
  confounder.df$sample[confounder.df$sample == "06"] = "Tumor" #"MetastaticTumor"
  confounder.df$sample[confounder.df$sample == "10"] = "Normal" #"BloodDerivedNormal"
  confounder.df$sample[confounder.df$sample == "11"] = "Normal" #"SolidTissueNormal"
  confounder.df$plate <- sapply(confounder.df$ids, function(x){strsplit(x, split="-")[[1]][6]})
  ## in case, only patients sequenced with one specific Illumina platform are in scope
  if(conf$platform == "GA"){
    confounder.df = confounder.df[confounder.df$Platform == "GA",]
  } else if(conf$platform == "HiSeq"){
    confounder.df = confounder.df[confounder.df$Platform == "HiSeq",]
  }
  ## total reads filter
  if(conf$reads_filter){
    confounder.df = confounder.df[which(confounder.df$total.reads > conf$tot.reads.filter),]
  } 
  ## add total reads quantiles to confounder df
  confounder.df = confounder.df %>% mutate(quantile = ntile(total.reads, 4))
  return(confounder.df)
}

## 3_assess_batch_effects.R ----

#' calculate percentages explained variance of principle components for PCA plots
#' \code{calculate_percentage} 
#' @param pca_df: pca results
#' @param df: data frame with different PCs
#' @return return named vector with PC variance and PC names
#' \dontrun{
#' calculate_percentage(pca_df, df)
#' }
calculate_percentage <- function(pca_df, df){
  percentage <- round((pca_df$sdev*pca_df$sdev) / sum((pca_df$sdev*pca_df$sdev)) * 100, 2)
  percentage <- paste0( names(df), " (", as.character(percentage), "%", ")")
  return(percentage)
}

#' PCA scatter plots of batch effects
#' \code{pca_plot_batch} 
#' @param df:
#' @param colour_variable:
#' @param x_axis: principle component plotted on x axis
#' @param y_axis: principle component plotted on y axis
#' @param title: title of ggplot
#' @param percentage: names percentage vector (variance [%] explained by each PC)
#' @param centroids: boolean value (TRUE/FALSE) if all observations should be plotted or only 
#' centroids e.g. of plate variable
#' @return ggplot scatter plot of batch effects
#' \dontrun{
#' pca_plot_batch(df, colour_variable, x_axis, y_axis, title, percentage, centroids)
#' }
pca_plot_batch <- function(df, colour_variable, x_axis, y_axis, title, percentage, centroids){
  y = names(df)[which(names(df) == paste0("PC",y_axis))]
  if(centroids){
    df = aggregate(cbind(df[,"PC1"],df[,which(names(df) == y)])~df[,colour_variable],df,mean)
    names(df) = c(colour_variable, "PC1", y)
  }
  if(colour_variable == "total.reads"){
    plot <- ggplot(df,aes(x=df[,x_axis], y=df[,y])) + 
      geom_point(alpha=0.7, aes(colour=df[,colour_variable]), size = 3) + 
      labs(x = percentage[1], y = percentage[y_axis], title = title, col = colour_variable)+
      theme_classic(base_size = 18)
  }else{
    plot <- ggplot(df,aes(x=df[,x_axis], y=df[,y])) + 
      geom_point(alpha=0.7, aes(colour=df[,colour_variable]), size = 3) + 
      labs(x = percentage[1], y = percentage[y_axis], title = title, col = colour_variable)+
      guides(colour = guide_legend(override.aes = list(alpha=1), title.position="top"))+
      theme_classic(base_size = 18)
  }
  return(plot)
}

## 3_figure_5.R ----

#' function to get p values 
#' (association tested between variables and PCs in 3_assess_batch_effects.R)
#' \code{get_pvals} 
#' @param include_all_samples: boolean value (TRUE/FALSE); TRUE if all samples should be included, 
#' false if only tumor samples are in scope of analysis
#' @param conf.batch_cor: batch correction configuration data frame
#' @param filter.method: expression filter method ("sum" or "median")
#' @param filter.value: expression filter value (0 if filter.method = sum, 15 if filter.method = 
#' median)
#' @return data frame with p values of association between confounders and PCs before and after
#' batch correction
#' \dontrun{
#' get_pvals(include_all_samples, conf.batch_cor, filter.method, filter.value)
#' }
get_pvals= function(include_all_samples, conf.batch_cor, filter.method, filter.value){
  ## load p values & results 
  ## in case, all samples are used, always take sample as variable of interest
  if(include_all_samples){
    conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == TRUE]
  } else {
    conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == FALSE]
  }
  all_pvals = llply(1:nrow(conf.batch_cor), function(i){
    conf = list()
    conf$base_dir = "/abi/data/sibing/gitWorkspace/isomiRs/BatchCor_TCGAisomiR" # REPLACE PERSONAL PATH WITH ~ !!!
    conf$normalization = "rpm"
    conf$tot.reads.filter = 1000000
    conf$filter = filter.method
    conf$filter_value = filter.value
    conf$Project = conf.batch_cor$Project[i]
    conf$platform = conf.batch_cor$platform[i]
    conf$rm_batches = conf.batch_cor$rm_batches[i]
    conf$batch_1 = conf.batch_cor$batch_1[i]
    conf$batch_2 = conf.batch_cor$batch_2[i]
    conf$correction_method = conf.batch_cor$correction_method[i]
    conf$var_of_interest = conf.batch_cor$var_of_interest[i]
    
    if(reads_filter){
      intermediate.path = "/results/"
    } else {
      intermediate.path = "/results/wo_reads_filter/"
    }
    
    print(conf$Project)
    
    figure.path_before = paste0(conf$base_dir,intermediate.path, "before_batch_cor/",  conf$Project, 
                                "/", "include.all.samples", include_all_samples, "/", 
                                conf$normalization, "_", conf$filter, conf$filter_value,"_", 
                                conf$platform)
    figure.path_after = paste0(conf$base_dir,intermediate.path, "batch_cor", "/",
                               conf$correction_method, "/", conf$Project, "/", "include.all.samples", 
                               include_all_samples, "/", conf$normalization, "_", conf$filter, 
                               conf$filter_value, "_",conf$var_of_interest, "_", conf$batch_1, "_", 
                               conf$batch_2, "_", conf$platform, "_", "rm",
                               paste0(c(conf$rm_batches), collapse = "", sep="_"))
    
    p_val_before = fread(paste0(figure.path_before, "/",conf$Project,"_p.values.pcs_", 
                                   conf$filter, conf$filter_value,".csv"), stringsAsFactors = F,
                         data.table = F)
    p_val_before$TimePoint = "before"
    p_val_after = fread(paste0(figure.path_after, "/",conf$Project,"_p.values.pcs_", 
                                  conf$filter, conf$filter_value,".csv"), stringsAsFactors = F,
                           data.table = F)
    p_val_after$TimePoint = "after"
    
    p_vals = rbind(p_val_before, p_val_after)
    p_vals$Project = conf$Project
    p_vals$correction_method = conf$correction_method
    p_vals$batch1 = conf$batch_1
    p_vals$batch2 = conf$batch_2
    p_vals$all_samples = include_all_samples
    
    return(p_vals)  })
  
  all_pvals_df = do.call(rbind, all_pvals)
  
  return(all_pvals_df)
}

#' Function to create bubble plots (figure 5); pan cancer comparison of reduction of batch effects
#' \code{bubble_plot} 
#' @param df: p value data frame derived from get_pvals function
#' @param legend: boolean value => should legend be included?
#' @param ticks: boolean value => should axis tick be included?
#' @param time_point: which data set (before/after batch correction)
#' @param variables: vector of variables to plot, e.g. c("plate", "platform", "purity", "reads")
#' @param PCs: vector of PCs to plot, e.g. c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
#' @param filter: filter method ("median15" or "sum0")
#' @return ggplot (bubbleplot)
#' \dontrun{
#' bubble_plot(df, legend, ticks, time_point, variables, PCs, filter)
#' }
bubble_plot = function(df, legend, ticks, time_point, variables, PCs, filter){
  
  df = df[which(df$component %in% PCs & df$variable %in% variables & df$TimePoint == time_point),]
  ## size of points in bubble plot
  
  if(filter == "median"){
    df$size = 0
    df$size[df$percentage <= 5] = 10
    df$size[df$percentage > 5 & df$percentage <= 10] = 12.5
    df$size[df$percentage > 10 & df$percentage <= 15] = 15
    df$size[df$percentage > 15 & df$percentage <= 20] = 17.5
    df$size[df$percentage > 20] = 20
  } else {
    df$size = 0
    df$size[df$percentage <= 1] = 10
    df$size[df$percentage > 1 & df$percentage <= 2] = 12.5
    df$size[df$percentage > 2 & df$percentage <= 4] = 15
    df$size[df$percentage > 4 & df$percentage <= 6] = 17.5
    df$size[df$percentage > 6] = 20
  }
  
  
  df$size[is.na(df$p.value)] = NA
  df$Project = gsub("TCGA-", "", df$Project)
  
  df$log_fdr = -log10(df$fdr)
  
  
  df$log_fdr[df$log_fdr > 130] = 130
  ## significance threshold: fdr adjusted p-value below 0.01
  
  Farben <- colorRampPalette(c("darkblue", "white", "darkred"))
  bunt <- Farben(900)
  
  df$Project2 = as.factor(df$Project)
  
  p=ggplot(df, aes(x = component, y = Project2, size = size, fill = log_fdr)) +
    geom_point(shape = 21,na.rm = TRUE, show.legend= TRUE)+ 
    theme_bw()+ theme(aspect.ratio=1)+
    facet_wrap(variable~., ncol = 1)+
    scale_y_discrete(limits = rev(levels(df$Project2)))+
    scale_fill_gradientn(colors= bunt,
                         guide = "colorbar",
                         values = scales::rescale(c(0, 1.9999, 2, 2.0001, 130)),
                         breaks=c(2,25,50, 75, 100, 125),
                         labels=c("-log10(0.01)", 25,50, 75, 100, 125))+
    labs(x = "PC", y = "TCGA Project")+
    theme_classic2()
  
  if(!legend){
    p = p+
      theme(legend.position = "none") 
  } 
  
  if(!ticks){
    p = p+
      theme(axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.text.y=element_blank())
  }
  
  return(p)
}



## 4_Define_outliers_platforms[...].R ----
        
#' calculate nucleotide occurence in isomiR sequences
#' \code{countCharOccurrences} counts frequency of a given letter in a string 
#' @param a letter to look for , e.g. "A" and a string to find the letter in, e.g. "ATTGCATTGGT"
#' @return number of occurences of the vector in the string
#' \dontrun{
#' countCharOccurrences("A", "ATTGCATTGGT")
#' }
countCharOccurrences <- function(char, s){
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2)) 
} 