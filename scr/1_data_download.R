#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018','forcats_0.5.0', 'stringr_1.4.0', 'dplyr_1.0.2', 
'readr_1.3.1', 'tidyr_1.1.2', 'tibble_3.0.3', 'ggplot2_3.3.2', 'tidyverse_1.3.0', 'plyr_1.8.6',
'doParallel_1.0.15', 'iterators_1.0.12', 'foreach_1.5.0', 'TCGAbiolinks_2.12.6', 'optparse_1.6.6'
'data.table_1.13.0', 'bsub_1.0.2')
Description = 'Data download and preparation of TCGA isomiR and miRNA expression quantification 
data.'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## LOAD BSUB LIBRARY ----
library(bsub)
library(data.table)
bsub_opt$R_version = "3.6.0"

base_dir = "~/BatchCor_TCGAisomiR" 
## load config data and extract projects
conf.batch_cor = fread(file.path(base_dir, "data", "config.batch_cor.txt"))
projects = unique(conf.batch_cor$Project)

for(Project in projects){
  
  ## EXECUTION ----
  ## bsub_chunk from bsub package (https://github.com/jokergoo/bsub) submits R code to LSF cluster 
  ## without leaving R
  bsub_chunk(name = paste0("prepare.matrices.", Project), 
             memory = 380, hour = 12, core = 8, 
             packages = c("tidyverse", "data.table"),
             variables = c("Project", "base_dir"),{
               
             options("expressions"=500000)            
             
             ## LIBRARIES  ----
             ## REMOVE LIB PATHS !!!
             library(TCGAbiolinks,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/') 
             library(parallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
             library(doParallel,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
             library(plyr,lib='/home/ibing/R/x86_64-pc-linux-gnu-library/3.6/')
             
             ## PARAMETER SETTINGS ----
             conf = list()
             ## Parameters for the download
             ## for all data.categories: TCGAbiolinks:::getProjectSummary("TCGA-BRCA")
             conf$DataCategory = "Transcriptome Profiling" 
             conf$Access = "open"
             conf$DataBase = "harmonized"
             conf$isomiR_dir = file.path(base_dir,"data/tcga_data_isomiR")
             conf$miRNA_dir = file.path(base_dir,"data/tcga_data_miRNA")
             conf$mapping_dir = file.path(base_dir,"data/mapping_tcga_data")
             conf$output_dir = file.path(base_dir, "data/DataClean")
             ## directory for data without annotation correction
             conf$output_wrong = file.path(base_dir, "data/DataWrongAnnotation")
             conf$function_dir = file.path(base_dir, "scr", "functions.batch_cor.R")
             
             source(conf$function_dir) # load functions
             
             if(!dir.exists(conf$isomiR_dir)){dir.create(conf$isomiR_dir, recursive = TRUE)}
             if(!dir.exists(conf$miRNA_dir)){dir.create(conf$miRNA_dir, recursive = TRUE)}
             
             data <- data.frame("data_category" = c("Transcriptome Profiling", 
                                                    "Transcriptome Profiling"),
                                "data_type" = c("Isoform Expression Quantification", 
                                                "miRNA Expression Quantification"),
                                "name" = c("isomiR", "miRNA"),
                                "workflow_type" = c("BCGSC miRNA Profiling", 
                                                    "BCGSC miRNA Profiling"),
                                stringsAsFactors = F)
             conf$DataType = list("Isoform Expression Quantification", 
                                  "miRNA Expression Quantification")
             
             ## go through projects and download data
             ## save all results files as mapping files (file names/sample ids/etc.)
             ## save queries as df in case of downloading them later on 

             if(!dir.exists(file.path(conf$output_wrong, Project))){
               dir.create(file.path(conf$output_wrong, Project), recursive = TRUE)
             }
             if(!dir.exists(file.path(conf$output_dir, Project))){
               dir.create(file.path(conf$output_dir, Project), recursive = TRUE)
             }
             if(!dir.exists(file.path(conf$mapping_dir,Project))){
               dir.create(file.path(conf$mapping_dir,Project), recursive = TRUE)
             }
             ## DOWNLOAD DATA ----
             
             ## download both isomiR and miRNA data
             ## sometimes data download is interrupted => in case so, please start data download for 
             ## cohort in scope separately
             for(i in 1:nrow(data)){
               ## get download query of the current project and data type
               query_down = GDCquery(project=Project,
                                     data.category=data$data_category[i],
                                     data.type=data$data_type[i],
                                     workflow.type=data$workflow_type[i],
                                     legacy = FALSE,
                                     access = conf$Access)
               results = query_down[[1]][[1]]
               ## download data and generate "results" df that is stored to enable later mapping of 
               ## file and case id:
               if(i == 1){ # isomiR data
                 GDCdownload(query = query_down, directory = conf$isomiR_dir, method = "api")
               } else if(i == 2){ # miRNA data
                 GDCdownload(query = query_down, directory = conf$miRNA_dir, method = "api")
               }
               assign(paste0("query_",data$name[i]),query_down)
               assign(paste0("results_",data$name[i]),results)
               ## results df: mapping file
               write.csv(results,file=file.path(conf$mapping_dir,paste0(Project,"/results_",
                                                                        data$name[i],".csv")))
             }
             
             ## ISOFORM EXPRESSION QUANTIFICATION FILES PREPROCESSING ----
             ## INTERSECT WITH ADAPTED ISOMIR FILE FROM miRBASE v22 
             ## Reference file with isomiR annotation is based on the hsa gff3 file from the miRBase 
             ## (version 22.1)
             ## IsomiR information (+- 3 nt at 3' and 5' end) was added to the miRNA annotation of 
             ## the gff3 file and saved as hsa_isomiR_v22_R.csv
             if(!file.exists(paste0(base_dir, "/data/hsa_isomirs_v22_R.csv"))){
               warning("Store first sheet of Supplementary_Table_1.xlsx in data folder")
             } else {
               ## isomiR file: first sheet in Supplementary_Table_1.xlsx
               hsa_isomiRs = read.csv(file.path(base_dir, "data/hsa_isomirs_v22_R.csv"), 
                                      header = TRUE, sep = " ", stringsAsFactors = F)
             }
             
             ## LOAD DATA FILES AS MATRIX; CORRECT ANNOTATION AND INTERSECT WITH ISOMIR POSITIONS --
             FilePath <- paste0("/", Project, "/", conf$DataBase, "/", conf$DataCategory,  "/", 
                                conf$DataType[[1]], "/")
             FilePath <- gsub(" ", "_", FilePath) # path to files
             dataFiles_list = lapply(1:dim(results_isomiR[1])[1], function(i){
               file_path = paste0(base_dir, "/data/tcga_data_isomiR", FilePath, 
                                  results_isomiR$file_id[i], "/",results_isomiR$file_name[i])
               pid = as.character(results_isomiR$cases[i])
               return(list(file_path, pid))
             })
             
             ## go through all samples, read data and build expression matrix
             exprs_ls = list()
             uncor.exprs_ls = list()
             for (i in 1:length(dataFiles_list)){
               pid = as.character(dataFiles_list[[i]][2])
               ## read data
               txt.file = read.table(file = as.character(dataFiles_list[[i]][1]), header = T)
               txt.file = txt.file %>% separate(isoform_coords, into = c("referenceGenome", 
                                                                         "chrom", "position", 
                                                                         "strand"), sep = ":") %>% 
                 separate(position, into = c("start", "end"), sep = "-")
               txt.file$start = as.numeric(txt.file$start)
               txt.file$end = as.numeric(txt.file$end)
               
               ## for this file, no annotation correction (comparison)
               uncor.file = txt.file
               uncor.intersect <- left_join(hsa_isomiRs, txt.file, by = c("chrom", "start", "end", 
                                                                          "strand"), copy = FALSE) 
               uncor.intersect$Name <- gsub("Name=", "", uncor.intersect$Name) # remove Name=
               uncor.intersect$pid = pid
               uncor.exprs_ls[[i]] = uncor.intersect
               
               ## adapt end positions to correct annotation
               txt.file$end <- txt.file$end - 1
               
               intersect <- left_join(hsa_isomiRs, txt.file, by = c("chrom", "start", "end", 
                                                                    "strand"), copy = FALSE) 
               intersect$Name <- gsub("Name=", "", intersect$Name) # remove Name=
               intersect$pid = pid
               exprs_ls[[i]] = intersect
             }
             exprs_df = do.call(rbind, exprs_ls)
             uncor_exprs_df = do.call(rbind, uncor.exprs_ls)
             
             ## CREATE MATRICES, NORMALIZE AND SAVE DATA (count file; normalized: RPM, TMM; UQ in 
             ## DataClean)
             ## filter out all non-mature miRNAs 
             exprs_df = exprs_df[exprs_df$miRNA == "miRNA",]
             uncor_exprs_df = uncor_exprs_df[uncor_exprs_df$miRNA == "miRNA",]
             
             ## create count matrix
             count.isomiRs = create_count_matrix(df = exprs_df)
             uncor.count.isomiRs = create_count_matrix(df = uncor_exprs_df)
             
             rpm.isomiRs = create_rpm_matrix(df = exprs_df)
             uncor.rpm.isomiRs = create_rpm_matrix(df = uncor_exprs_df)
             
             data.table::fwrite(exprs_df, 
                                file = file.path(conf$output_dir, Project, 
                                                 paste0(Project, "_clean.csv")), 
                                row.names = FALSE, col.names = TRUE)
             data.table::fwrite(count.isomiRs, 
                                file = file.path(conf$output_dir, Project, 
                                                 paste0(Project, "_count.csv")), 
                                row.names = TRUE, col.names = TRUE)
             data.table::fwrite(rpm.isomiRs, 
                                file = file.path(conf$output_dir, Project, 
                                                 paste0(Project, "_rpm.csv")), 
                                row.names = TRUE, col.names = TRUE)
             
             data.table::fwrite(uncor_exprs_df, 
                                file = file.path(conf$output_wrong, Project, 
                                                 paste0(Project, "_clean.csv")), 
                                row.names = FALSE, col.names = TRUE)
             data.table::fwrite(uncor.count.isomiRs, 
                                file = file.path(conf$output_wrong, Project, 
                                                 paste0(Project, "_count.csv")), 
                                row.names = TRUE, col.names = TRUE)
             data.table::fwrite(uncor.rpm.isomiRs, 
                                file = file.path(conf$output_wrong, Project, 
                                                 paste0(Project, "_rpm.csv")), 
                                row.names = TRUE, col.names = TRUE)
             
             ## MIRNA EXPRESSION QUANTIFICATION PREPROCESSING  ----
             ## LOAD DATA FILES AS MATRIX; CORRECT ANNOTATION AND INTERSECT WITH ISOMIR POSITIONS --
             FilePath <- paste0("/", Project, "/", conf$DataBase, "/", conf$DataCategory,  "/", 
                                conf$DataType[[2]], "/")
             FilePath <- gsub(" ", "_", FilePath) # path to files
             print(FilePath)
             dataFiles_list = lapply(1:dim(results_miRNA[1])[1], function(i){
               file_path = paste0(base_dir, "/data/tcga_data_miRNA", FilePath, 
                                  results_miRNA$file_id[i], "/",results_miRNA$file_name[i])
               pid = as.character(results_miRNA$cases[i])
               return(list(file_path, pid))
             })

             exprs_ls = list()
             for (i in 1:length(dataFiles_list)) {
               pid = as.character(dataFiles_list[[i]][2])
               txt.file = read.table(file = as.character(dataFiles_list[[i]][1]), header = TRUE)
               txt.file$pid = pid
               exprs_ls[[i]] = txt.file
             }
             exprs_df = do.call(rbind, exprs_ls)

             ## create count matrix
             count.matrix = exprs_df[,c("pid", "read_count", "miRNA_ID")] %>% 
               spread(key = pid, value = read_count)
             count.matrix[is.na(count.matrix)] = 0
             ## sum up same entries (isomiRs derived from different genomic locations) 
             ## => counts were devided by the number of loci
             count.matrix = 
               data.table(count.matrix)[, lapply(.SD[, names(count.matrix)[6:ncol(count.matrix)], 
                                                     with = FALSE], sum), by = miRNA_ID] # faster
             rownames(count.matrix) = count.matrix$miRNA_ID
             count.matrix$miRNA_ID = NULL

             ## create rpm matrix (provided by TCGA)
             rpm.matrix = exprs_df[,c("pid", "reads_per_million_miRNA_mapped", "miRNA_ID")] %>% 
               spread(key = pid, value = reads_per_million_miRNA_mapped)
             rpm.matrix[is.na(rpm.matrix)] = 0
             rpm.matrix = 
               data.table(rpm.matrix)[, lapply(.SD[, names(rpm.matrix)[6:ncol(rpm.matrix)], 
                                                   with = FALSE], sum), by = miRNA_ID]
             rownames(rpm.matrix) = rpm.matrix$miRNA_ID
             rpm.matrix$miRNA_ID = NULL

             data.table::fwrite(exprs_df, file = file.path(conf$output_dir, Project, 
                                                           paste0(Project, "_clean_miRNA_only.csv")), 
                                row.names = FALSE, col.names = TRUE)
             data.table::fwrite(count.matrix, file = file.path(conf$output_dir, Project, 
                                                               paste0(Project, "_count_miRNA_only.csv")), 
                                row.names = TRUE, col.names = TRUE)
             data.table::fwrite(rpm.matrix, file = file.path(conf$output_dir, Project, 
                                                             paste0(Project, "_rpm_miRNA_only.csv")), 
                                row.names = TRUE, col.names = TRUE)
             })
}
