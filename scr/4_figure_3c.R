#!/usr/bin/env Rscript
library(optparse)

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'dplyr_1.0.2', 'optparse_1.6.6', 'plyr_1.8.6')
Description = 'This script is to produce Figure 3c. 
Input: Output files from 4_define_outliers.R'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## LIBRARIES ----
library(plyr)
library(dplyr)

## PARAMETER SETTINGS ----
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$output_dir = file.path(conf$base_dir, "results", "Figure3")
if(!dir.exists(conf$output_dir)){dir.create(conf$output_dir, recursive = T)}

## Figure 3c -----
files_ttest_LUSC <- dir(conf$output_dir, pattern= "^TCGA-LUSC.+ttest\\.csv$")
setwd(conf$output_dir)
data_ttest_LUSC <- lapply(files_ttest_LUSC, read.csv2)

t <- join_all(data_ttest_LUSC,  type = "inner")
t <- t[ with(t, order(p_batch_cor_TCGA.LUSC )), ]
t$num <- c(1:nrow(t))
x1 <- sort(t$padj_before_batch_cor )
x2 <- sort(t$padj_batch_cor )

pdf(file.path(conf$output_dir, "Figure3c.pdf"))
plot(c(1:length(x1)),x1, col = "#d95f02", pch= 16, xaxt='n', 
     ylab = "q-value", xlab = "isomiR", cex=1.5, cex.axis= 2,cex.lab=1.5) 
points(c(1:length(x2)),x2, pch = 20, col = "#1b9e77",  cex=1.2) #
legend("topleft",  pch= c(16,16), bty = "n", col = c( "#d95f02","#1b9e77"), 
       legend= c("before","after" ), cex= 2)
dev.off()    
