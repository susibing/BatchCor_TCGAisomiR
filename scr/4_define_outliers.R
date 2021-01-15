#!/usr/bin/env Rscript
library(optparse)
library(bsub)

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'ggplot2_3.3.2', 'optparse_1.6.6', 'data.table_1.13.0')
Description = 'This script is to define outliers comparing Genome Analyzer and HiSeq sequenced 
samples before and after batch correction for the TCGA-LUSC dataset and produces Figure 3a,b,d,e 
as well as Supplementary Figure 6a and b. Additionally, outliers between the two platforms before 
batch correction are determined for the 8 datasets comprising samples sequenced on Genome Analyzer
and HiSeq and results are exported for use in Figure 4. 
Input: count/rpm-normalized expression matrices before and after batch correction.'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## PARAMETER SETTINGS ----
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$output_dir = file.path(conf$base_dir, "results/Figure3")
conf$normalization_type = "rpm"
conf$filter = "median" # median or sum filter
conf$filter_value = 15 # filter threshold (15: median, 0: sum)
conf$include_all_samples = FALSE # T in case normal samples should be included
conf$reads_filter = TRUE # T in case total reads filter should be applied
conf$tot.reads.filter = 1000000 # threshold of total reads filter

if(!dir.exists(conf$output_dir)){dir.create(conf$output_dir, recursive = T)}

data_type = "before_batch_cor" # "before_batch_cor" (before) or "batch_cor" (after)

## Submit job to LSF cluster
bsub_chunk(name = paste0("figure.3.", data_type),
           variables = c("conf", "data_type"), 
           memory = 25, hour = 1, core = 8, packages = c("plyr", "dplyr", "tidyverse", 
                                                         "RColorBrewer"),{
                                                           
 source(file.path(conf$base_dir, "scr", "functions.batch_cor.R")) # load functions 
 
 library(wrapr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
 library(doParallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
 library(matrixStats, lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
 library(genefilter, lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
 library(data.table, lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
 library(gplots, lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
 
 registerDoParallel(cores=8)
 
 ## Define whether data before or after batch correction is in scope.
 ## Correction parameters: 
 ## Optimal combination of batch correction for the cohorts
 conf.batch_cor = fread(file.path(conf$base_dir,"data","config.batch_cor.txt"),stringsAsFactors=F)
 if(conf$include_all_samples){
   conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == TRUE,]
 } else {
   conf.batch_cor = conf.batch_cor[conf.batch_cor$include_all_samples == FALSE,]
 }
 
 ## filter down to projects sequenced both using the HiSeq and GA platform
 conf.batch_cor = conf.batch_cor[conf.batch_cor$Project %in% c("TCGA-BRCA", "TCGA-COAD", 
                                                               "TCGA-HNSC", "TCGA-KIRC", 
                                                               "TCGA-LUAD", "TCGA-LUSC", 
                                                               "TCGA-STAD", "TCGA-UCEC")]
 
 # samples are only further analysed for "after batch correction" for the LUSC dataset: 
 if (data_type == "batch_cor"){
   conf.batch_cor = conf.batch_cor[conf.batch_cor$Project %in% c("TCGA-LUSC")]
 }
 
 ## go though different parameter settings (in parallel)
 # only apply to LUSC dataset for this figure (number 11 in conf)
 ## ONLY 11 IF NORMAL SAMPLES ARE NOT INCLUDED ! 
 llply(1:nrow(conf.batch_cor), function(i){
   print(i)
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
       conf$include_all_samples = F 
     } 
     # correction for purity in the presence of normal samples is not possible 
     # as normal samples per definition should not contain tumor cells --> only the second confounder was used in this case. 
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
     figure_path = paste0(conf$output_dir, "/", conf$Project, "_", data_type,  "_", 
                          "include.all.samples", conf$include_all_samples, "_", 
                          conf$normalization_type, "_", conf$filter, conf$filter_value,"_", 
                          conf$platform)
     table_output = paste0(conf$output_dir, "/", conf$Project,"_", data_type, "_", 
                           "include.all.samples", conf$include_all_samples, "_", 
                           conf$normalization_type, "_", conf$filter, conf$filter_value)
     exprs_df = as.data.frame(fread(paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", 
                                           conf$Project, "_", conf$normalization_type, ".csv"), 
                                    sep=",", check.names = F, stringsAsFactors = F, header = T))
   } else {
     ## batch corrected data file
     figure_path = paste0(conf$output_dir,"/", conf$Project, "_",data_type, "_", 
                          conf$correction_method, "_", "include.all.samples", 
                          conf$include_all_samples, "_", conf$normalization_type, "_", conf$filter, 
                          conf$filter_value, "_",conf$var_of_interest, "_", conf$batch_1, "_", 
                          conf$batch_2, "_", conf$platform, "_", "rm",paste0(c(conf$rm_batches), 
                                                                            collapse = "", sep="_"))
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
                                             conf$normalization_type, ".csv"), sep=",", check.names = F, stringsAsFactors = F, header = T))
     }else{
       exprs_df = as.data.frame(fread(paste0(conf$base_dir, "/data/DataClean/", conf$Project, "/", data_type, "_wo_reads_filter", "/", conf$correction_method, "/", conf$Project,"_", conf$filter, conf$filter_value, "_",conf$var_of_interest, "_", conf$batch_1, "_", conf$batch_2, "_", conf$platform, "_", paste0(c(conf$rm_batches), collapse = "", sep="_"),"all_samples_",conf$include_all_samples, "_" , conf$normalization_type, ".csv"), sep=",", check.names = F, stringsAsFactors = F, header = T))
     }
   }
   
   exprs_df = makeRn(exprs_df)
   
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
     exprs_df = exprs_df[,which(grepl('[0-9A-Z]{2}.[0-9A-Z]{4}.02[0-9A-Z]{1}.',colnames(exprs_df)) | grepl('[0-9A-Z]{2}.[0-9A-Z]{4}.01[0-9A-Z]{1}.',colnames(exprs_df)))]
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
     counts = which(median.val > conf$filter_value)
     exprs_df = exprs_df[counts, ]
     exprs_df = as.data.frame(exprs_df)
     median_expr <- data.frame(row.names(exprs_df), rowMedians(as.matrix(exprs_df)))
     names(median_expr) <- c("ids", paste0("row_median_", conf$Project))
     median_expr$ids <-  gsub("-", '.', median_expr$ids )
     median_expr$ids <-  gsub("\\|", '.', median_expr$ids )
   }
   
   
   ## LOG-transform data
   exprs_log = log2((as.matrix(exprs_df) + (min(exprs_df[exprs_df > 0]) / 10)))
   exprs_log = t(exprs_log)
   
   ## load confounders matrix
   confounders = read.table(paste0(conf$base_dir, "/data/TCGA_confounder_table.txt"), header=T, 
                            stringsAsFactors = F)
   confounder.df = prep_conf_df(exprs_log = exprs_log, confounders = confounders, conf = conf)
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
   
   ## ttest comparisons between platforms
   exprs_log_platform <- as.data.frame(cbind(rownames(exprs_log),
                                             data.frame(exprs_log,row.names=NULL)))
   colnames(exprs_log_platform)[1] <- "ids"
   exprs_log_platform <- join(confounder.df[,c("ids","Platform")], exprs_log_platform, type="right")
   
   if (nrow(subset(exprs_log_platform, Platform == "GA"))>0){
     # calculate ttest and FDR
     t_exprs <- exprs_log_platform[with(exprs_log_platform, order(Platform)), ]
     t_exprs <- subset(t_exprs, Platform != "NA", ) 
     rownames(t_exprs) <- t_exprs[,1]
     numberGA <- nrow(subset(t_exprs, Platform =="GA"))
     numberHi <- nrow(subset(t_exprs, Platform == "HiSeq"))
     t_exprs<- t_exprs[,-c(1,2)]
     t_mat <- t(as.matrix(t_exprs))
     
     g <- factor(c(rep(0, numberGA),rep(1,numberHi)))
     ttest <- rowttests(t_mat, g)
     ttest <- cbind(rownames(ttest), data.frame(ttest, row.names=NULL))
     ttest$padj <- p.adjust(ttest$p.value,  method = "BH")
     ttest <- ttest[,c("rownames(ttest)", "p.value", "padj")]
     sign_ttest <- subset(ttest, padj < 0.05)
     names(sign_ttest)[1] <- "ids"
     names(ttest) <- c("ids", paste0("p_", data_type,"_", conf$Project),  
                       paste0("padj_", data_type,"_", conf$Project))
     
     
     # add mean expression of GA and HiSeq
     mean_expr <- sapply(c("GA","HiSeq"), function(x){
       exprs_log_platf <- subset(exprs_log_platform, Platform == x)
       rownames(exprs_log_platf) <- exprs_log_platf$ids
       exprs_log_platf <- exprs_log_platf[,-c(1:2)]
       exprs_log_platf <- as.data.frame(t(exprs_log_platf))
       rowMeans(exprs_log_platf)
     })
     colnames(mean_expr) <- c((paste(data_type,"GA_mean", conf$Project ,sep = "_")), 
                              (paste(data_type,"HiSeq_mean", conf$Project ,sep = "_")))
     statistics <- cbind(ttest, mean_expr)
     
     if (data_type=="before_batch_cor"){
       statistics <- join(statistics, median_expr)
       write.csv2(statistics, paste0(table_output, "_ttest.csv"), row.names= FALSE)
     } else {
       if (conf$Project == "TCGA-LUSC"){
         write.csv2(statistics, paste0(table_output, "_ttest.csv"), row.names= FALSE)
       }
     }
     
   }
   
   
   ## Figure 3 -----
   if (conf$Project == "TCGA-LUSC"){
     ## Figure 3a or b -----
     # a or b depending on data_type configuration
     if (data_type == "before_batch_cor"){ subfig = "a" } else { subfig = "b" }
     
     # calculate Z-scores for same number of GA as HiSeq samples
     # do random sampling for HiSeq data/ GA data (depending on where less samples are)
     platform.summary <- as.data.frame(table(exprs_log_platform$Platform), stringsAsFactors = F)
     platform.summary <- platform.summary[platform.summary$Var1 %in% c("GA", "HiSeq"),]# exclude NAs
     min.platform <- platform.summary$Var1[platform.summary$Freq == min(platform.summary$Freq)]
     number <- min(platform.summary$Freq)
     
     set.seed(42)
     if(min.platform == "GA"){
       HiSeq_rand <- sample_n(subset(exprs_log_platform, Platform == "HiSeq"), number)
       GA_HiSeq_rand <- rbind(subset(exprs_log_platform, Platform == "GA"), HiSeq_rand)
     } else {
       GA_rand <- sample_n(subset(exprs_log_platform, Platform == "GA"), number)
       GA_HiSeq_rand <- rbind(GA_rand, subset(exprs_log_platform, Platform == "HiSeq"))
     }
     
     rownames(GA_HiSeq_rand) <- GA_HiSeq_rand[,1]
     # Z-scaling
     z_exprs_log <- as.matrix(GA_HiSeq_rand[,5:ncol(GA_HiSeq_rand)])
     z_exprs_log <- scale(z_exprs_log, center= TRUE, scale = TRUE)
     z_exprs_log[1:5,1:5]
     z_exprs_log <- cbind(rownames(z_exprs_log), data.frame(z_exprs_log, row.names = NULL))
     colnames(z_exprs_log)[1] <- "ids"
     # annotate z-expressed data with platforms
     z_exprs_platform <- join(confounder.df[,c("ids","Platform")], as.data.frame(z_exprs_log), 
                              type = "right")
     
     # get mean zscore expression of GA and HiSeq
     mean_zscore <- sapply(c("GA","HiSeq"), function(x){
       z_exprs_platf <- subset(z_exprs_platform, Platform == x)
       rownames(z_exprs_platf) <- z_exprs_platf$ids
       z_exprs_platf <- z_exprs_platf[,-c(1:2)]
       z_exprs_platf <- as.data.frame(t(z_exprs_platf))
       rowMeans(z_exprs_platf)
     })
     
     # heatmap comparing Z-scores between the two different platforms
     png(filename = paste0(figure_path, "Figure3",subfig,".png"), width = 680, height = 680)
     my_palette <- colorRampPalette(c("#2668A8", "#ffffbf","#CB0017"))(n = 299)
     
     colors = c(seq(-1,-0.13,length=100),
                seq(-0.1299,0.1299,length= 100),
                seq(0.13,1,length=100))
     
     heatmap.2(mean_zscore,Colv = FALSE,
               dendrogram = "row",
               col= my_palette,
               breaks= colors,
               scale="none",
               trace="none",
               cexCol = 0.6,
               srtCol=45,labRow = FALSE,symkey = FALSE,key.title = NA)
     dev.off()
     
     ## Figure 3d or e , left panel-----
     # d or e depending on data_type configuration
     if (data_type == "before_batch_cor"){
       subfig = "d"
     } else {subfig ="e"}
     
     pdf(paste0(figure_path, "Figure3",subfig,".pdf"))
     mean_expr_df <- as.data.frame(mean_expr)
     mean_expr_df <- cbind(rownames(mean_expr_df), data.frame(mean_expr_df, row.names = NULL))
     names(mean_expr_df) <- c("ids", "GA_mean", "HiSeq_mean")
     
     # label example isomiRs
     labels_ex_iso <- filter(mean_expr_df,ids %in% c( "hsa.miR.143.3p.0.0.","hsa.miR.22.3p.2.0."))
     
     # define outlier
     outlier <- join(sign_ttest, mean_expr_df, type = "left")
     outlier$FC <- outlier$GA_mean - outlier$HiSeq_mean
     outlier_down <- subset(outlier, FC<=(-0.5))
     outlier_up <- subset(outlier, FC>=0.5)
     
     plot(mean_expr_df$HiSeq_mean, mean_expr_df$GA_mean,
          xlim = c(0,round(max(mean_expr_df$HiSeq_mean),2)),
          ylim = c(0,round(max(mean_expr_df$GA_mean),2)),
          pch= 16, col= "#fc8d59",cex=1.2, cex.axis= 2,cex.lab=1.5, 
          xlab= "HiSeq sequencing - isomiR expressions (log2)", 
          ylab= "GA sequencing - isomiR expression (log2)" )
     points(outlier_down$HiSeq_mean, outlier_down$GA_mean, pch= 16, col= "#4575b4", cex= 1.2)
     points(outlier_up$HiSeq_mean, outlier_up$GA_mean, pch= 16, col= "#d73027", cex= 1.2)
     legend("topleft", pch=c(16,16,16), bty = "n", col = c("#fc8d59", "#4575b4", "#d73027"), 
            legend = c("mean", "overrep. in HiSeq" ,"overrep. in GA"), cex=2)
     if(data_type== "before_batch_cor"){
       points(labels_ex_iso $HiSeq_mean,labels_ex_iso $GA_mean, pch= 1, col= "black")
       text(labels_ex_iso $HiSeq_mean,labels_ex_iso $GA_mean, labels = labels_ex_iso $ids , 
            pos = c(3,1), cex=2)
     } 
     
     # calculate spearman correlation of figure
     corri <- cor.test(mean_expr_df$HiSeq_mean, mean_expr_df$GA_mean,method = "spearman", 
                       exact = FALSE)
     corri2 <- capture.output(print(corri))
     write.table(corri2, paste0(table_output,"spearman_cor_", data_type, ".txt"), row.names = FALSE)
     
     ## Figure 3d or e, middle and right panel -----
     # d or e depending on data_type configuration
     # prepare vectors for density plots
     for (p in c("GA", "HiSeq")){
       expr_log_p <- subset(exprs_log_platform, Platform == p)
       if (p == "GA"){
         GA_hist1 <- expr_log_p$hsa.miR.143.3p.0.0.
         GA_hist2 <- expr_log_p$hsa.miR.22.3p.2.0.
       }else if (p =="HiSeq"){
         Hi_hist1 <- expr_log_p$hsa.miR.143.3p.0.0.
         Hi_hist2 <- expr_log_p$hsa.miR.22.3p.2.0.
       }}  
     
     c1 <- makeTransparent("#d73027")
     c2 <- makeTransparent("#4575b4")
     
     d<- density(GA_hist1, from = 2, to = 17.5)
     plot(d, main = "hsa-miR-143-3p|0|0|", xlab= "expression (log2)",ylim=c(0,0.7), cex=1.5, 
          cex.axis= 2,cex.lab=2 )
     polygon(d, col = c1, border= "black")
     polygon(density(Hi_hist1), col = c2, border= "black")
     legend("topleft",  pch= c(15,15), bty = "n", col = c(c1, c2), legend= c("GA","HiSeq" ), cex=2)
     
     d<- density(GA_hist2, from = 2, to = 17.5)
     plot(d, main = "hsa-miR-22-3p|2|0|", xlab= "expression (log2)", ylim= c(0,0.72), cex=1.5, 
          cex.axis= 2,cex.lab=2 )
     polygon(d, col = c1, border= "black")
     polygon(density(Hi_hist2), col = c2, border= "black")
     legend("topleft",  pch= c(15,15), bty = "n", col = c(c1, c2), legend= c("GA", "HiSeq"), cex=2)
     
     dev.off()
     
     ## Supplementary Figure 6 a and b ----
     ## Heatmaps comparing isomiR expression according to sequencing on different plates
     # a or b depending on data_type configuration
     if (data_type == "before_batch_cor"){
       subfig = "a"
     } else {subfig ="b"}
     
     # calculate Z-scores for same number of GA as HiSeq samples
     # do random sampling for HiSeq data
     # annotate z-expressed data with plates
     z_exprs_plate <- join(confounder.df[,c("ids","plate")], as.data.frame(z_exprs_log), 
                           type = "right")
     
     # get mean zscore expression for each plate
     plates <- levels(unique(z_exprs_plate$plate))
     mean_zscore_plate <- sapply(plates, function(x){
       z_exprs_plate <- subset(z_exprs_plate, plate == x)
       rownames(z_exprs_plate) <- z_exprs_plate$ids
       z_exprs_plate <- z_exprs_plate[,-c(1:2)]
       z_exprs_plate <- as.data.frame(t(z_exprs_plate))
       rowMeans(z_exprs_plate)
     })
     
     # heatmap comparing Z-scores between the two different platforms
     pdf(paste0(figure_path, "SupFigure6",subfig,".pdf"))
     my_palette <- colorRampPalette(c("#313695", "#fee090", "#a50026"))(n = 299) 
     
     colors = c(seq(-1,-0.13,length=100),
                seq(-0.1299,0.1299,length= 100),
                seq(0.13,1,length=100))
     
     heatmap.2(mean_zscore_plate,Colv = FALSE,
               dendrogram = "row",
               col= my_palette,
               breaks= colors,
               scale="none",
               trace="none",
               cexCol = 0.6,
               srtCol=45,labRow = FALSE,symkey = FALSE,key.title = NA)
     dev.off()
   }
   
   
 }, .parallel = FALSE)#TRUE)  
 
})
