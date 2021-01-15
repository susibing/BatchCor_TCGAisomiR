#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'plyr_1.8.6', 'data.table_1.13.0', 'dplyr_1.0.2',
'bsub_1.0.2', 'optparse_1.6.6')
Description = 'Calculate distribution of canonical miRs among all isomiRs before and after 
annotation correction. Tumor samples from the TCGA-LUSC cohort used for 
figure'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## Load bsub library
library(bsub)
base_dir = "~/BatchCor_TCGAisomiR" 

## Submit job to LSF cluster
bsub_chunk(name = "figure_1d", memory = 5, hour = 1, core = 2, 
           packages = c("plyr", "dplyr", "data.table"),
           variables = "base_dir",{
  
  ## CONFIG PART ----
  conf = list()
  conf$Project = "TCGA-LUSC"
  conf$normalization_type = "rpm"
  conf$correction_method = "combat"
  conf$filter = "sum"
  conf$filter_value = 0
  conf$var_of_interest = "none"
  conf$include_all_samples = FALSE
  conf$reads_filter = TRUE # TRUE in case total reads filter should be applied
  conf$tot.reads.filter = 1000000 # threshold of total reads filter
  conf$annotation = c("before_annot_cor", "after_annot_cor")
  data_type = "before_batch_cor"
  conf$platform = "both"
  conf$output_dir = file.path(base_dir, "results", data_type, conf$Project, 
                              paste0("include.all.samples", conf$include_all_samples),
                              paste0(conf$normalization_type, "_", conf$filter, conf$filter_value,
                                     "_", conf$platform))
  
  if(!dir.exists(conf$output_dir)){dir.create(conf$output_dir, recursive = T)}
  
  ## load functions and libraries
  source(file.path(base_dir, "scr", "functions.batch_cor.R")) # load functions 
  
  ## LIBRARIES ----
  library(matrixStats,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(gplots,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(RColorBrewer,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(stringr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')

  llply(conf$annotation, function(annotation){
    ## read expression data and define output paths
    if(annotation == "before_annot_cor"){
      figure_path = paste0(conf$output_dir, "/Figure1d_", data_type, "_", conf$Project, "_", 
                           "include.all.samples", conf$include_all_samples, "_", annotation)
      exprs_df = as.data.frame(fread(paste0(base_dir, "/data/DataWrongAnnotation/", conf$Project, "/", 
                                            conf$Project, "_", conf$normalization_type, ".csv"), 
                                     sep=",", check.names = F, stringsAsFactors = F, header = T))
    } else {
      ## annotation corrected data file
      figure_path = paste0(conf$output_dir, "/Figure1d_", data_type, "_", conf$Project, "_", 
                           "include.all.samples", conf$include_all_samples, "_", annotation)
      exprs_df = as.data.frame(fread(paste0(base_dir, "/data/DataClean/", conf$Project, "/", 
                                            conf$Project, "_", conf$normalization_type, ".csv")))
      }

    rownames(exprs_df) = exprs_df[,1]
    exprs_df[,1] = NULL
    
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
    
    exprs_df <- as.data.frame(t(as.matrix(exprs_df)))
    exprs_df <- cbind(rownames(exprs_df), data.frame(exprs_df, row.names = NULL))
    colnames(exprs_df)[1] <- "ids"
    
    ## load confounders matrix; merge with exprs_df
    confounders = read.table(paste0(base_dir, "/data/TCGA_confounder_table.txt"), header=T, 
                             stringsAsFactors = F)
    exprs_conf_df <- join(confounders[,c("ids",  "total.reads")], exprs_df, type = "right")
    
    ## total reads filter
    if(conf$reads_filter){
      exprs_conf_df = exprs_conf_df[which(exprs_conf_df$total.reads > conf$tot.reads.filter),]
    } 
    print(exprs_conf_df[1:5,1:5])
    exprs_df_filt <- exprs_conf_df[,-2]
    print(exprs_df_filt[1:5,1:5])  
    # Identification canonical miRs, stems of miRs and isomiRs
    cans <- names(exprs_df_filt)[which(str_sub(names(exprs_df_filt), start= -5) == ".0.0.")]
    stems <- str_sub(cans, end = -6)

    ## get sum of isomiR expression and expression of the respective canonical isomiR
    getIsomiRsums <- function(x){
      isos <- grep(x, names(exprs_df_filt), fixed= TRUE, value = TRUE)
      isomiRs <- exprs_df_filt[,c("ids",  isos)]
      isomiRs.t <- as.data.frame(t(as.matrix(makeRn(isomiRs))))
      isomiRs.t$sum <- rowSums(as.matrix(isomiRs.t))
      ## coerce matrix containing sums of RPM expression values for all isomiRs and RPM expression 
      ## value for the respective canoncial isomiR
      c(colSums(as.matrix(isomiRs.t$sum)), subset(isomiRs.t$sum, isos==paste0(x, ".0.0.")))
    }
    
    sumsRPM <- sapply(stems, getIsomiRsums)
    sums_RPM_df <- as.data.frame(t(sumsRPM))
    names(sums_RPM_df) <- c("iso_sums", "can_sums")
    
    sums_RPM_df$quot <- (sums_RPM_df$can_sums/sums_RPM_df$iso_sums)*100

    #identify isomiRs where the canconcial reads account for less than 50% of the total reads:
    iso.higher.canon <- subset(sums_RPM_df, quot< 50)
    # calculation of the percentage of miRs where the canonical accounts for less that 50% of reads:
    perc.can <-round((nrow(iso.higher.canon)/nrow(sums_RPM_df))*100,2)
    perc.can
    
    ## Density plots 
    pdf(paste0(figure_path, ".pdf"))
    
    par(mai=c(1,1,0.5,0.5))
    
    c2 <- makeTransparent("darkred")
    c1 <- makeTransparent("grey75")
    
    d<- density(sums_RPM_df$quot )#from = 0, to = 120
    e<- density(sums_RPM_df$quot,from = 0, to = 50 )
    plot(d, main = paste("Canonical isomiR distribution", annotation), xlim= c(3.5,96.5),
         ylim= c(0,0.055),
         xlab= "%of isomiRs that are canoncial",ylab= "density", cex=1.5, cex.axis= 2,cex.lab=2 )
    polygon(d, col = c1, border= "black")
    polygon(c(0,e$x),c(0,e$y), col = c2, border= "black")
    legend("left", bty = "n", legend= paste(perc.can, "% of miRNAs"), cex=2)
    dev.off()

    })
    
})
