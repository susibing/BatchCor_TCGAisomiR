#!/usr/bin/env Rscript
library(optparse)
library(plyr)
library(dplyr)

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'dplyr_1.0.2', 'optparse_1.6.6', 'plyr_1.8.6')
Description = 'This script is to reproduce Figure 4 graphs comparing samples sequenced on HiSeq or 
GA,respectively concering sequence features. Input: Results from running 
4_define_outliers_platforms_all_entities.R and the isomiRs_sequences.csv file 
(provided in the github repository /data). '
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## PARAMETER SETTINGS ----
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$data_dir = file.path(conf$base_dir, "results", "Figure3")
conf$output_dir = file.path(conf$base_dir, "results", "Figure4/")
if(!dir.exists(conf$output_dir)){dir.create(conf$output_dir, recursive = T)}
source(file.path(conf$base_dir, "scr","functions.batch_cor.R")) # load functions

## Figure 4 -----
# read files containing statistical results
files_wGA_seq <- dir(conf$data_dir, pattern= "^.+before_batch_cor.+ttest\\.csv$")
setwd(conf$data_dir)
data_wGA_seq <- lapply(files_wGA_seq, read.csv2)
# data_wGA_seq <- lapply(data_wGA_seq, function(x){x[,1] <- NULL; x })
data_wGA_seq <- lapply(data_wGA_seq, function(x){
  FC = x[,4] -x[,5]
  cbind(x, FC)
})

data_wGA_seq <- lapply(data_wGA_seq, function(x){
  colnames(x)[7] <- paste0("FC_", str_sub(colnames(x[6]),start=-4)) 
  x
})


# Compile list containing identifiers which are expressed with a median higher 15 in 6 out of 
# 8 datasets
df_wGA_seq <-  join_all(data_wGA_seq, by = "ids", type='full')
wGA_seq_row_medians <- df_wGA_seq[,which(str_sub(colnames(df_wGA_seq), 
                                                 start = 1, end = 10) == "row_median" | 
                                           colnames(df_wGA_seq)=="ids")]
wGA_seq_row_medians$na_count <- apply(is.na(wGA_seq_row_medians), 1, sum)
wGA_seq_row_medians <- subset(wGA_seq_row_medians, na_count<3)

write.csv2(wGA_seq_row_medians[,c(1:9)], paste0(conf$output_dir,"isogr15in6of8.csv"))

# filter t-tests
wGA_seq_ttests <- df_wGA_seq[,-which(str_sub(colnames(df_wGA_seq), 
                                             start = 1, end = 10) == "row_median")]
wGA_seq_ttests <- wGA_seq_ttests[,which(colnames(wGA_seq_ttests)=="ids"
                                    |str_sub(colnames(wGA_seq_ttests), start = 1, end = 2)=="FC"
                                    |str_sub(colnames(wGA_seq_ttests), start = 1, end = 4)=="padj")]
# set NAs to 0 for FC
wGA_seq_ttests[c("FC_BRCA", "FC_COAD", "FC_HNSC", "FC_KIRC", "FC_LUAD", "FC_LUSC", "FC_STAD", 
                 "FC_UCEC")][is.na(wGA_seq_ttests[c("FC_BRCA", "FC_COAD", "FC_HNSC", "FC_KIRC", 
                                                  "FC_LUAD", "FC_LUSC", "FC_STAD", "FC_UCEC")])]<-0
wGA_seq_ttests[c("padj_before_batch_cor_TCGA.BRCA", "padj_before_batch_cor_TCGA.COAD", 
                 "padj_before_batch_cor_TCGA.HNSC", "padj_before_batch_cor_TCGA.KIRC", 
                 "padj_before_batch_cor_TCGA.LUAD", "padj_before_batch_cor_TCGA.LUSC",
                 "padj_before_batch_cor_TCGA.STAD","padj_before_batch_cor_TCGA.UCEC")][is.na(wGA_seq_ttests[
                   c("padj_before_batch_cor_TCGA.BRCA", "padj_before_batch_cor_TCGA.COAD", 
                        "padj_before_batch_cor_TCGA.HNSC","padj_before_batch_cor_TCGA.KIRC", 
                        "padj_before_batch_cor_TCGA.LUAD", "padj_before_batch_cor_TCGA.LUSC",
                        "padj_before_batch_cor_TCGA.STAD", "padj_before_batch_cor_TCGA.UCEC")])]<-1


# Filtering conditions: padj < 0.05 AND FC > 0.5 in 6 out of 8 samples
#Check for which datasets filtering conditions are true/false:
difference = 0.5
wGA_seq_ttests$BRCA <- wGA_seq_ttests$padj_before_batch_cor_TCGA.BRCA < 0.05 & 
  abs(wGA_seq_ttests$FC_BRCA)>difference
wGA_seq_ttests$COAD <- wGA_seq_ttests$padj_before_batch_cor_TCGA.COAD<0.05 & 
  abs(wGA_seq_ttests$FC_COAD)> difference
wGA_seq_ttests$HNSC <- wGA_seq_ttests$padj_before_batch_cor_TCGA.HNSC<0.05 & 
  abs(wGA_seq_ttests$FC_HNSC)> difference
wGA_seq_ttests$KIRC <- wGA_seq_ttests$padj_before_batch_cor_TCGA.KIRC<0.05 & 
  abs(wGA_seq_ttests$FC_KIRC)> difference
wGA_seq_ttests$LUAD <- wGA_seq_ttests$padj_before_batch_cor_TCGA.LUAD<0.05 & 
  abs(wGA_seq_ttests$FC_LUAD)> difference
wGA_seq_ttests$LUSC <- wGA_seq_ttests$padj_before_batch_cor_TCGA.LUSC<0.05 & 
  abs(wGA_seq_ttests$FC_LUSC)> difference
wGA_seq_ttests$STAD <- wGA_seq_ttests$padj_before_batch_cor_TCGA.STAD<0.05 & 
  abs(wGA_seq_ttests$FC_STAD)> difference
wGA_seq_ttests$UCEC <- wGA_seq_ttests$padj_before_batch_cor_TCGA.UCEC<0.05 & 
  abs(wGA_seq_ttests$FC_UCEC)> difference

wGA_seq_ttests$true_count <- apply(wGA_seq_ttests[,c((ncol(wGA_seq_ttests)-7):ncol(wGA_seq_ttests))]
                                   =="TRUE", 1, sum)
wGA_seq_ttests <- subset(wGA_seq_ttests, true_count>5)
head(wGA_seq_ttests)
nrow(wGA_seq_ttests)
wGA_seq_ttests$mean_FC <- apply(wGA_seq_ttests[,c("FC_BRCA", "FC_COAD", "FC_HNSC", 
                                                  "FC_KIRC", "FC_LUAD", "FC_LUSC",
                                                  "FC_STAD", "FC_UCEC")],1,mean)

overrep <- subset(wGA_seq_ttests, mean_FC>0)
underrep <- subset(wGA_seq_ttests, mean_FC<0)

# read sequences for all isomiRs
seq_all <- fread(paste0(conf$base_dir,"/data/", "isomiRs_sequences.csv"), header = TRUE, 
                 stringsAsFactors = F, data.table = F)
head(seq_all)

overrep_seq <- join(seq_all, overrep, type = "right")
overrep_seq <- overrep_seq[,c(1,2)]

underrep_seq <- join(seq_all, underrep, type = "right")
underrep_seq <- underrep_seq[,c(1,2)]

non_aff_seq <- anti_join(seq_all, overrep_seq, by= "ids")
non_aff_seq <- anti_join(non_aff_seq, underrep_seq, by= "ids")

write.csv2(overrep_seq, paste0(conf$output_dir, "list_overrep_isos.csv"))
write.csv2(underrep_seq, paste0(conf$output_dir, "list_underrep_isos.csv"))
write.csv2(non_aff_seq, paste0(conf$output_dir, "list_non-aff_isos.csv"))

seq_list <- list(overrep_seq, underrep_seq, non_aff_seq)
names(seq_list) <- c("overrep_seq", "underrep_seq", "non_aff_seq")

seqLengthGC_list <- lapply(seq_list, function(x) {
  x$length <- nchar(x$sequence)
  x$G <- countCharOccurrences("G",x$sequence)
  x$C <- countCharOccurrences("C",x$sequence)
  x$A <- countCharOccurrences("A",x$sequence)
  x$T <- countCharOccurrences("T",x$sequence)
  x$GC <- x$G+x$C
  allnuc <- x$A+x$G+x$C+x$T
  x$percGC <- x$GC/allnuc*100
  x[,c(1,3,9)]
} )

seqLengthGC_list <-
  lapply(names(seqLengthGC_list), function(i){
    x <- seqLengthGC_list[[ i ]]
    names(x)[2] <- paste0("length_", i)
    names(x)[3] <- paste0("GC_", i)
    x
  })

seqLengthGC_df = Reduce(function(...) merge(..., all=T), seqLengthGC_list)

# calculate nucleotide distribution in first position and GC content in first position
firstPos_list <- llply(seq_list, function(x) {
  df <- as.data.frame(table(substr(x$sequence,1,1)))
  df$perc <- df$Freq/sum(df$Freq)*100
  rownames(df) <- df[,1]
  df[,1] <- NULL
  df <- as.data.frame(t(df))
  df$GC <- df$C+df$G
  df$AT<- df$A+df$T
  df <- as.data.frame(t(df))
  df <- cbind(rownames(df), data.frame(df, row.names = NULL))
  colnames(df)[1] <- "base"
  df
})
firstPos_list

firstPos_list <-
  lapply(names(firstPos_list), function(i){
    x <- firstPos_list[[ i ]]
    names(x)[2] <- paste0("abs_", i)
    names(x)[3] <- paste0("perc_", i)
    x
  })

firstPos_df = Reduce(function(...) merge(..., all=T), firstPos_list)
write.csv2(firstPos_df, paste0(conf$output_dir,"SuppTableFig4_FirstPosNucleotide.csv"))

pdf(paste0(conf$output_dir, "Figure4.pdf"))
## Figure 4a -----
## violin plot GC content
d2 <- seqLengthGC_df[,which(str_sub(names(seqLengthGC_df), start = 1, end = 2) == "GC" )]
names(d2) <- c(" overrep.", "underrep.", "non-aff.")
head(d2)

l1 = data.frame(x=c(1,1,2,2),y=c(77,78.5,78.5,77))
l2 = data.frame(x=c(2,2,3,3),y=c(82,83.5,83.5,82))

q<- d2 %>% 
  gather(key="isomiRs", value="GC") %>%
  ggplot( aes(x=isomiRs, y=GC, fill= isomiRs)) +
  geom_violin(trim= TRUE, na.rm = TRUE)+
  scale_fill_manual(values=c( "#d73027", "#fc8d59", "#4575b4"))+
  geom_boxplot(width=0.1, na.rm = TRUE) + 
  stat_summary(fun=median, geom="point", size=2, color="black", na.rm = TRUE) + 
  theme_classic(base_size = 30) +  theme(legend.position="none")+
  #scale_y_continuous(breaks = round(seq(range(na.omit(d2))[1], range(na.omit(d2))[2], by = 50),1))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 30),
        axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size= 30), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  labs( x= "isomiR selection",y = "% GC content") + 
  geom_path(data= l1, aes(x = x,y = y), inherit.aes = FALSE)+
  geom_path(data= l2, aes(x = x,y = y), inherit.aes = FALSE)+
  annotate("text",x=1.5,y=81,label="ns", size = 9)+
  annotate("text",x=2.5, y=85, label="**", size = 10)
q 

# calculate results ttests
(t_over_nonaff_GC <- t.test(d2$` overrep.`,d2$`non-aff.`))
(t_under_nonaff_GC <- t.test(d2$underrep.,d2$`non-aff.`))

## Figure 4b -----
# heatmap of nucleotide distribution in first position of isomiRs
pos1 <- firstPos_df[c(1,3,4,6),c("base","perc_overrep_seq","perc_underrep_seq", "perc_non_aff_seq")]
rownames(pos1) <- pos1[,1]
pos1<- pos1[,c(2,4,3)]
my_palette <- colorRampPalette(c("#2668A8", "#ffffbf","#CB0017"))(n = 299)

heatmap.2(as.matrix(pos1),Colv = FALSE,
          dendrogram = "row",
          col= my_palette,
          trace="none",
          cexCol = 0.6,
          srtCol=45,
          symkey = FALSE,key.title= NA,keysize = 2,
          main =  "")

# Fishers exact test for significance
pos1_abs <- firstPos_df[c(1,3,4,6),c("base","abs_overrep_seq","abs_non_aff_seq","abs_underrep_seq")]
rownames(pos1_abs) <- pos1_abs[,1]
pos1_abs[,1] <- NULL

(fish_A_pos1_overrep <- fisher.test(matrix(c(pos1_abs[1,1],pos1_abs[1,2], 
                                             sum(pos1_abs[2:4,1]),sum(pos1_abs[2:4,2])), 
                                           nrow=2, byrow=TRUE)))
(fish_C_pos1_overrep <- fisher.test(matrix(c(pos1_abs[2,1],pos1_abs[2,2], 
                                             sum(pos1_abs[c(1,3:4),1]),sum(pos1_abs[c(1,3:4),2])), 
                                           nrow=2, ncol=2, byrow=TRUE)))
(fish_G_pos1_overrep <- fisher.test(matrix(c(pos1_abs[3,1],pos1_abs[3,2], 
                                             sum(pos1_abs[c(1:2,4),1]),sum(pos1_abs[c(1:2,4),2])), 
                                           nrow=2, ncol=2, byrow=TRUE)))
(fish_T_pos1_overrep <- fisher.test(matrix(c(pos1_abs[4,1],pos1_abs[4,2], 
                                             sum(pos1_abs[c(1:3),1]),sum(pos1_abs[1:3,2])), 
                                           nrow=2, ncol=2, byrow=TRUE)))

(fish_A_pos1_underrep <- fisher.test(matrix(c(pos1_abs[1,2],pos1_abs[1,3], 
                                              sum(pos1_abs[2:4,2]),sum(pos1_abs[2:4,3])), 
                                            nrow=2, ncol=2, byrow=TRUE)))
(fish_C_pos1_underrep <- fisher.test(matrix(c(pos1_abs[2,2],pos1_abs[2,3], 
                                              sum(pos1_abs[c(1,3:4),2]),sum(pos1_abs[c(1,3:4),3])), 
                                            nrow=2, ncol=2, byrow=TRUE)))
(fish_G_pos1_underrep <- fisher.test(matrix(c(pos1_abs[3,2],pos1_abs[3,3], 
                                              sum(pos1_abs[c(1:2,4),2]),sum(pos1_abs[c(1:2,4),3])), 
                                            nrow=2, ncol=2, byrow=TRUE)))
(fish_T_pos1_underrep <- fisher.test(matrix(c(pos1_abs[4,2],pos1_abs[4,3], 
                                              sum(pos1_abs[c(1:3),2]),sum(pos1_abs[1:3,3])), 
                                            nrow=2, ncol=2, byrow=TRUE)))

## Figure 4c -----
# barplot showing GC content in first position of isomiRs
pos1gc <- firstPos_df[c(2,5),c("base","perc_overrep_seq" , "perc_underrep_seq", "perc_non_aff_seq")]
pos1gc_abs <- firstPos_df[c(2,5),c("abs_overrep_seq" , "abs_non_aff_seq", "abs_underrep_seq")]
rownames(pos1gc) <- pos1gc[,1]
names(pos1gc) <- c("base", "overrep.", "underrep.", "non-aff.")
pos1gc <- pos1gc[c(2,1),c(2,4,3)]

barplot(as.matrix(pos1gc),
        #xlab=pos1$ids, 
        density = c(5,40),
        col=c( rep("#d73027",2),rep("#fc8d59",2), rep("#4575b4",2)),
        beside=TRUE, ylim= c(0,100), cex.axis = 2.5, cex.names = 2.5)
legend("topleft", legend= paste0("%",rownames(pos1gc)),  density=c(5,40), cex=2.5,  bty = "n" )

(fish_GC_pos1_overrep <- fisher.test(pos1gc_abs[,1:2]))
(fish_GC_pos1_underrep <- fisher.test(pos1gc_abs[,2:3]))

## Figure 4d ----
## violin plot sequence length isomiRs
d1 <- seqLengthGC_df[,which(str_sub(names(seqLengthGC_df), start = 1, end = 6) == "length" )]
names(d1) <- c(" overrep.", "underrep.","non-aff.")
head(d1)

l1 = data.frame(x=c(1,1,2,2),y=c(27.5,27.7,27.7,27.5))
l2 = data.frame(x=c(2,2,3,3),y=c(28.0,28.2,28.2,28.0))

p<- d1 %>% 
  gather(key="isomiRs", value="length") %>%
  ggplot( aes(x=isomiRs, y=length, fill= isomiRs)) +
  geom_violin( trim= TRUE,na.rm = TRUE, adjust=2 ) + 
  scale_fill_manual(values=c("#d73027", "#fc8d59", "#4575b4"))+
  geom_boxplot(width=0.1, na.rm = TRUE) + 
  stat_summary(fun=median, geom="point", size=2, color="black", na.rm = TRUE) +
  theme_classic(base_size = 30)+theme(legend.position="none")+
  scale_y_continuous(breaks = round(seq(range(na.omit(d1$`non-aff.`))[1], 
                                        range(na.omit(d1$`non-aff.`))[2], by = 2),1))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 30),
        axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size= 30),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  labs( x= "isomiR selection",y = "length [nt]")+
  geom_path(data= l1, aes(x = x,y = y), inherit.aes = FALSE)+
  geom_path(data= l2, aes(x = x,y = y), inherit.aes = FALSE)+
  annotate("text",x=1.5,y=27.8,label="**", size = 10)+
  annotate("text",x=2.5,y=28.3,label="****", size = 10)
p
dev.off()

# t-tests
(t_over_nonaff_length <- t.test(d1$` overrep.`,d1$`non-aff.`))
(t_under_nonaff_length <- t.test(d1$underrep.,d1$`non-aff.`))

significances <- c(t_over_nonaff_GC$p.value, t_under_nonaff_GC$p.value, 
                   fish_A_pos1_overrep$p.value, fish_C_pos1_overrep$p.value,
                   fish_G_pos1_overrep$p.value, fish_T_pos1_overrep$p.value,
                   fish_A_pos1_underrep$p.value, fish_C_pos1_underrep$p.value,
                   fish_G_pos1_underrep$p.value, fish_T_pos1_underrep$p.value,
                   fish_GC_pos1_overrep$p.value, fish_GC_pos1_underrep$p.value,
                   t_over_nonaff_length$p.value, t_under_nonaff_length$p.value)
names(significances) <- c("t_over_nonaff_GC", "t_under_nonaff_GC", 
                          "fish_A_pos1_overrep", "fish_C_pos1_overrep",
                          "fish_G_pos1_overrep", "fish_T_pos1_overrep",
                          "fish_A_pos1_underrep", "fish_C_pos1_underrep",
                          "fish_G_pos1_underrep", "fish_T_pos1_underrep",
                          "fish_GC_pos1_overrep" , "fish_GC_pos1_underrep",
                          "t_over_nonaff_length", "t_under_nonaff_length")
significances <- as.data.frame(significances)
write.csv2(significances, paste0(conf$output_dir, "Figure4_calc_significances.csv"))

