#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'data.table_1.13.0', 'tidyverse_1.3.0', 
'matrixStats_0.57.0', 'plotrix_3.7-8', 'bsub_1.0.2', 'optparse_1.6.6')
Description = 'Create colored expression matrix with differences before and after annotation 
correction. Tumor samples from the TCGA-LUSC cohort used for figure'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## Load bsub library
library(bsub)

base_dir = "~/BatchCor_TCGAisomiR" 

bsub_chunk(name = "figure_1c", memory = 40, hour = 18, core = 8, 
           packages = c("tidyverse", "data.table"),
           variables = "base_dir",{
             
  ## LIBRARIES ----
  library(plotrix,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  library(matrixStats,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
  
  ## CONFIG PART ----
  ## Correction parameters for TCGA-LUSC
  conf = list()
  conf$Project = "TCGA-LUSC"
  conf$platform = "both"
  conf$rm_batches = c()
  conf$filter = "sum"
  conf$filter_value = 0
  conf$include_all_samples = FALSE
  
  normalization_type = "rpm"
  data_type = "before_batch_cor"
  
  figure.path = file.path(base_dir, "results",data_type, conf$Project, paste0("include.all.samples", 
                          conf$include_all_samples), paste0(normalization_type, "_", conf$filter, 
                                                           conf$filter_value,"_", conf$platform))
  if(!dir.exists(figure.path)){dir.create(figure.path, recursive = T)}
  
  miRNA.arm = "hsa-let-7a-5p" # miRNA arm in scope of expression overview
  
  source(file.path(base_dir, "scr", "functions.batch_cor.R")) # load functions 
  
  ## LOAD DATA ----
  exprs.after.annot = load_annotation_data(path = file.path(base_dir, "data/DataClean", conf$Project, 
                                             paste0(conf$Project, "_", normalization_type, ".csv")))#, exprs.after = exprs.after)
  exprs.before.annot = load_annotation_data(path = file.path(base_dir, "data/DataWrongAnnotation", 
                               conf$Project, paste0(conf$Project, "_", normalization_type, ".csv")))#, exprs.after = exprs.after)
  
  # prepare matrix after correction 
  exprs_plot = exprs.after.annot[exprs.after.annot$miRNA.arm == miRNA.arm,
                                 c("median", "isoform", "three.prime", "five.prime")]
  exprs_plot = exprs_plot[,c("median", "three.prime", "five.prime")] %>% spread(key = five.prime, 
                                                                                value = median)
  rownames(exprs_plot) = exprs_plot$three.prime
  exprs_plot$three.prime = NULL
  # change order of the matrix to have it from -3 until +3 in numerical order
  matrix.after = as.matrix(exprs_plot[c(3,2,1,4:7),c(3,2,1,4:7)])
  
  # prepare matrix before correction
  exprs.plot.before = exprs.before.annot[exprs.before.annot$miRNA.arm == miRNA.arm,
                                         c("median", "isoform", "three.prime", "five.prime")]
  exprs.plot.before = exprs.plot.before[,c("median", "three.prime", "five.prime")] %>% 
    spread(key = five.prime, value = median)
  rownames(exprs.plot.before) = exprs.plot.before$three.prime
  exprs.plot.before$three.prime = NULL
  matrix.before = as.matrix(exprs.plot.before[c(3,2,1,4:7),c(3,2,1,4:7)])
  
  # pdf with median exprs. of all isomiRs from all genomic locations
  pdf(file.path(figure.path, paste0(conf$Project,"_", miRNA.arm, "_isoform_exprs.pdf")),
      width = 6, height= 5)
  # define color scheme for matrix => use log2 scale otherwise only canonical form is colors)
  ## (too large differences between high and low expressed isoforms)
  cellcol.before<-color.scale(extremes = c("white", "#0868ac"),
                          cbind(log2(matrix.before + (min(matrix.before[matrix.before > 0]) / 10)),
                                c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.before,
                  main="Before annotation correction",
                  cellcolors=cellcol.before,
                  show.values=2,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])),
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol,
                       cex = vcex))
  # add x and y axis with wright order
  axis(side=2,at=0.5:(ncol(matrix.before)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.before)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  
  cellcol.after<-color.scale(extremes = c("white", "#0868ac"),
                             cbind(log2(matrix.after + (min(matrix.after[matrix.after > 0]) / 10)),
                                   c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.after,
                  main="After annotation correction",
                  cellcolors=cellcol.after,
                  show.values=2,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])), 
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol, 
                       cex = vcex))
  axis(side=2,at=0.5:(ncol(matrix.after)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.after)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  dev.off()
  
  # for all without labels
  pdf(paste0(figure.path, "/", conf$Project,"_", miRNA.arm, "_isoform_exprs_color.pdf"), width = 3, 
      height= 3)
  cellcol.before<-color.scale(extremes = c("white", "#0868ac"),
                          cbind(log2(matrix.before + (min(matrix.before[matrix.before > 0]) / 10)),
                                c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.before,
                  main="Before annotation correction",
                  cellcolors=cellcol.before,
                  show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position")
  # add x and y axis with wright order
  axis(side=2,at=0.5:(ncol(matrix.before)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.before)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  cellcol.after<-color.scale(extremes = c("white", "#0868ac"),
                             cbind(log2(matrix.after + (min(matrix.after[matrix.after > 0]) / 10)),
                                   c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.after,
                  main="After annotation correction",
                  cellcolors=cellcol.after,
                  show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position")
  axis(side=2,at=0.5:(ncol(matrix.after)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.after)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  dev.off()
  
  ### exprs annot figures for hsa-let-7a-5p (different genetic location)
  exprs.after.annot = data.table::fread(paste0(base_dir, "/data/DataClean/", conf$Project, "/", 
                                              conf$Project, "_clean.csv"), sep=",", check.names = F, 
                                        stringsAsFactors = F, header = T, data.table = F)
  exprs.before.annot = data.table::fread(paste0(base_dir, "/data/DataWrongAnnotation/", conf$Project, "/", 
                                              conf$Project, "_clean.csv"), sep=",", check.names = F, 
                                         stringsAsFactors = F, header = T, data.table = F)
  
  expr_data_chr = function(exprs, miRNA.arm, chrom){
    ## function prepares expression matrices to plot for miRNAs from different genomic loci
    exprs = exprs[grep(miRNA.arm, exprs$Name),] # extract expression of miRNA arm of interest
    exprs[is.na(exprs)] = 0 # set NA values to zero
    exprs.chr = exprs[exprs$chrom == chrom,] # extract data from chrom of interest
    exprs.chr = exprs.chr[,c("reads_per_million_miRNA_mapped", "Name", "pid")] %>% 
      spread(key = pid, value = reads_per_million_miRNA_mapped)
    rownames(exprs.chr) = exprs.chr$Name
    exprs.chr$Name = NULL
    exprs.chr$median = rowMedians(as.matrix(exprs.chr)) # calculate median values
    exprs.chr$miRNA.arm = sapply(rownames(exprs.chr), function(x){
      strsplit(x, split = "\\|")[[1]][1]})
    exprs.chr$isoform = sapply(rownames(exprs.chr), function(x){
      paste0(strsplit(x, split = "\\|")[[1]][2],strsplit(x, split = "\\|")[[1]][3])})
    exprs.chr$three.prime= sapply(rownames(exprs.chr), function(x){
      strsplit(x, split = "\\|")[[1]][3]})
    exprs.chr$five.prime= sapply(rownames(exprs.chr), function(x){
      strsplit(x, split = "\\|")[[1]][2]})
    exprs.chr = exprs.chr[,c("median", "isoform", "three.prime", "five.prime")]
    exprs.chr = exprs.chr[,c("median", "three.prime", "five.prime")] %>% 
      spread(key = five.prime, value = median)
    rownames(exprs.chr) = exprs.chr$three.prime
    exprs.chr$three.prime = NULL
    # change order of the matrix to have it from -3 until +3 in numerical order
    matrix.chr = as.matrix(exprs.chr[c(3,2,1,4:7),c(3,2,1,4:7)])
    return(matrix.chr)
  }
  
  matrix.after.9 = expr_data_chr(exprs = exprs.after.annot, miRNA.arm = miRNA.arm, chrom = "chr9")
  matrix.after.11 = expr_data_chr(exprs = exprs.after.annot, miRNA.arm = miRNA.arm, chrom = "chr11")
  matrix.after.22 = expr_data_chr(exprs = exprs.after.annot, miRNA.arm = miRNA.arm, chrom = "chr22")
  matrix.before.9 = expr_data_chr(exprs = exprs.before.annot, miRNA.arm = miRNA.arm, chrom = "chr9")
  matrix.before.11 = expr_data_chr(exprs = exprs.before.annot, miRNA.arm = miRNA.arm, chrom = "chr11")
  matrix.before.22 = expr_data_chr(exprs = exprs.before.annot, miRNA.arm = miRNA.arm, chrom = "chr22")
  
  
  ## export colored matrices for isomiRs from different genomic loci
  pdf(paste0(figure.path, "/", conf$Project,"_hsa-let-7a-5p_individual_isoform_exprs.pdf"), 
      width = 6, height= 5)
  # define color scheme for matrix => use log2 scale otherwise only canonical form is colores )(too large differences between high and low expressed isoforms)
  cellcol.before.9<-color.scale(extremes = c("white", "#0868ac"),
                    cbind(log2(matrix.before.9 + (min(matrix.before.9[matrix.before.9 > 0]) / 10)),
                          c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.before.9,
                  main="Before annotation correction (chr 9, plus strand)",
                  cellcolors=cellcol.before.9,
                  #main="Blue to red correlations",
                  show.values=2,
                  #show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])),
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol,
                       cex = vcex))
  # add x and y axis with wright order
  axis(side=2,at=0.5:(ncol(matrix.before.9)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.before.9)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  
  cellcol.after.9<-color.scale(extremes = c("white", "#0868ac"),
                       cbind(log2(matrix.after.9 + (min(matrix.after.9[matrix.after.9 > 0]) / 10)),
                             c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.after.9,
                  main="After annotation correction (chr 9, plus strand)",
                  cellcolors=cellcol.after.9,
                  #main="Blue to red correlations",
                  show.values=2,
                  #show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])), 
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol, 
                       cex = vcex))
  axis(side=2,at=0.5:(ncol(matrix.after.9)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.after.9)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  
  cellcol.before.11<-color.scale(extremes = c("white", "#0868ac"),
                 cbind(log2(matrix.before.11 + (min(matrix.before.11[matrix.before.11 > 0]) / 10)),
                       c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.before.11,
                  main="Before annotation correction (chr 11, minus strand)",
                  cellcolors=cellcol.before.11,
                  #main="Blue to red correlations",
                  show.values=2,
                  #show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])),
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol,
                       cex = vcex))
  # add x and y axis with wright order
  axis(side=2,at=0.5:(ncol(matrix.before.11)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.before.11)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  
  cellcol.after.11<-color.scale(extremes = c("white", "#0868ac"),
                    cbind(log2(matrix.after.11 + (min(matrix.after.11[matrix.after.11 > 0]) / 10)),
                          c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.after.11,
                  main="After annotation correction (chr 11, minus strand)",
                  cellcolors=cellcol.after.11,
                  #main="Blue to red correlations",
                  show.values=2,
                  #show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])), 
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol, 
                       cex = vcex))
  axis(side=2,at=0.5:(ncol(matrix.after.11)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.after.11)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  cellcol.before.22<-color.scale(extremes = c("white", "#0868ac"),
                  cbind(log2(matrix.before.22 + (min(matrix.before.22[matrix.before.22 > 0]) / 10)),
                        c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.before.22,
                  main="Before annotation correction (chr 22, plus strand)",
                  cellcolors=cellcol.before.22,
                  #main="Blue to red correlations",
                  show.values=2,
                  #show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])),
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol,
                       cex = vcex))
  # add x and y axis with wright order
  axis(side=2,at=0.5:(ncol(matrix.before.22)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.before.22)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  
  cellcol.after.22<-color.scale(extremes = c("white", "#0868ac"),
                      cbind(log2(matrix.after.22 + (min(matrix.after.22[matrix.after.9 > 0]) / 10)),
                            c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.after.22,
                  main="After annotation correction (chr 22, plus strand)",
                  cellcolors=cellcol.after.22,
                  #main="Blue to red correlations",
                  show.values=2,
                  #show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])), 
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol, 
                       cex = vcex))
  axis(side=2,at=0.5:(ncol(matrix.after.22)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.after.22)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  dev.off()
  
  pdf(paste0(figure.path, "/", conf$Project,"_hsa-let-7a-5p_isoform_individual_exprs_colors.pdf"), 
      width = 3, height= 3)
  # define color scheme for matrix => use log2 scale otherwise only canonical form is colores )(too large differences between high and low expressed isoforms)
  cellcol.before.9<-color.scale(extremes = c("white", "#0868ac"),
                    cbind(log2(matrix.before.9 + (min(matrix.before.9[matrix.before.9 > 0]) / 10)),
                          c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.before.9,
                  main="Before annotation correction (chr 9, plus strand)",
                  cellcolors=cellcol.before.9,
                  show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])),
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol,
                       cex = vcex))
  # add x and y axis with wright order
  axis(side=2,at=0.5:(ncol(matrix.before.9)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.before.9)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  
  cellcol.after.9<-color.scale(extremes = c("white", "#0868ac"),
                        cbind(log2(matrix.after.9 + (min(matrix.after.9[matrix.after.9 > 0]) / 10)),
                              c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.after.9,
                  main="After annotation correction (chr 9, plus strand)",
                  cellcolors=cellcol.after.9,
                  show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])), 
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol, 
                       cex = vcex))
  axis(side=2,at=0.5:(ncol(matrix.after.9)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.after.9)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  
  cellcol.before.11<-color.scale(extremes = c("white", "#0868ac"),
                  cbind(log2(matrix.before.11 + (min(matrix.before.11[matrix.before.11 > 0]) / 10)),
                        c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.before.11,
                  main="Before annotation correction (chr 11, minus strand)",
                  cellcolors=cellcol.before.11,
                  show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])),
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol,
                       cex = vcex))
  # add x and y axis with wright order
  axis(side=2,at=0.5:(ncol(matrix.before.11)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.before.11)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  
  cellcol.after.11<-color.scale(extremes = c("white", "#0868ac"),
                    cbind(log2(matrix.after.11 + (min(matrix.after.11[matrix.after.11 > 0]) / 10)),
                          c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.after.11,
                  main="After annotation correction (chr 11, minus strand)",
                  cellcolors=cellcol.after.11,
                  show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])), 
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol, 
                       cex = vcex))
  axis(side=2,at=0.5:(ncol(matrix.after.11)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.after.11)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  cellcol.before.22<-color.scale(extremes = c("white", "#0868ac"),
                  cbind(log2(matrix.before.22 + (min(matrix.before.22[matrix.before.22 > 0]) / 10)),
                        c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.before.22,
                  main="Before annotation correction (chr 22, plus strand)",
                  cellcolors=cellcol.before.22,
                  show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])),
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol,
                       cex = vcex))
  # add x and y axis with wright order
  axis(side=2,at=0.5:(ncol(matrix.before.22)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.before.22)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  
  
  cellcol.after.22<-color.scale(extremes = c("white", "#0868ac"),
                      cbind(log2(matrix.after.22 + (min(matrix.after.22[matrix.after.9 > 0]) / 10)),
                            c(-1,rep(1,6))),c(0,1),0,c(1,0))[,1:7]
  color2D.matplot(matrix.after.22,
                  main="After annotation correction (chr 22, plus strand)",
                  cellcolors=cellcol.after.22,
                  show.values = F,
                  axes=F,
                  xlab = "5' position",
                  ylab = "3' position",
                  text(sort(rep((1:xdim[2]) - 0.5, xdim[1])), 
                       texty, sprintf("%.2f", round(x, show.values)), col = vcol, 
                       cex = vcex))
  axis(side=2,at=0.5:(ncol(matrix.after.22)-0.5),labels=c("3", "2", "1", "0", "-1", "-2", "-3"))
  axis(side=1,at=0.5:(ncol(matrix.after.22)-0.5),labels=c("-3", "-2", "-1", "0", "1", "2", "3"))
  dev.off()
})
