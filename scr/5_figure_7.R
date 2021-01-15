#!/usr/bin/env Rscript
library(optparse)
library(gplots)
library(ggridges)
library(scales)
library(plyr)
library(dplyr)
library(ggplot2)

opt_parser = OptionParser(description = "
##_Title_###########################################################################################
Version = '0.0.1'
Date = '2020-10-22'
Dependencies = c('R version 3.6.0 (2019-04-26)',
'RStudio Version 1.1.463 – © 2009-2018', 'optparse_1.6.6', 'ggplot2_3.3.2', 'dplyr_1.0.2', 
'plyr_1.8.6', 'scales_1.1.1', 'ggridges_0.5.2', 'gplots_3.1.0')
Description = 'This script is to produce Figure 7 graphs comparing tumor and normal samples before 
and after batch correction for the TCGA-LUSC dataset. Input: Results from running 
5_expression_analysis_tumor_vs_normal.R once using data_type = before_batch_cor and once using 
data_type = batch_cor'
####################################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

## PARAMETER SETTINGS ----
conf = list()
conf$base_dir = "~/BatchCor_TCGAisomiR" 
conf$output_dir = file.path(conf$base_dir, "results", "Figure7")
if(!dir.exists(conf$output_dir)){dir.create(conf$output_dir, recursive = T)}

source(file.path(conf$base_dir, "scr", "functions.batch_cor.R")) # load functions

## Figure 7 ----
# read files containing statistical results
files_ttest <- dir(conf$output_dir, pattern= "ttest\\.csv$")
setwd(conf$output_dir)
data_ttest <- lapply(files_ttest, read.csv2)

t_test <-  join_all(data_ttest, by = "X", type='full')
names(t_test)[1] <- "id"

# identify significant samples 
log2FCboth <- subset(t_test, abs(log2TvsN_batch_cor ) >0.4 | abs(log2TvsN_before_batch_cor) > 0.4)

sig_after <- subset(t_test, padj_before_batch_cor > 0.05 & padj_batch_cor  < 0.05& abs(log2TvsN_batch_cor )>0.4) # significant after correction
sig_before <- subset(t_test, padj_before_batch_cor < 0.05 & padj_batch_cor  > 0.05& abs(log2TvsN_before_batch_cor)>0.4) # significant before correction 
sig_both <- subset(log2FCboth, padj_before_batch_cor < 0.05 & padj_batch_cor  < 0.05)

sig_p_after <- subset(t_test, p_before_batch_cor  > 0.05 & p_batch_cor  < 0.05& abs(log2TvsN_batch_cor )>0.4) # significant after correction
sig_p_before <- subset(t_test, p_before_batch_cor  < 0.05 & p_batch_cor  > 0.05& abs(log2TvsN_before_batch_cor)>0.4) # significant before correction 

sig_never_a <- t_test %>% anti_join(sig_before,by="id")
sig_never_b <- sig_never_a %>% anti_join(sig_after,by="id")
sig_never <- sig_never_b %>% anti_join(sig_both,by="id")

## Figure 7a -----
## scatter plot log2 FC before and after batch correction
pdf(file.path(conf$output_dir, "Figure7a.pdf"))
par(mai=c(1,1,0.5,0.5))

plot(t_test$log2TvsN_before_batch_cor,
     t_test$log2TvsN_batch_cor , 
     xlab= "log2 FC before correction",
     ylab="log2 FC after correction",
     pch=16,
     cex.lab=2,
     cex.axis =2,
     cex=1.2,
     abline(h=0, v=0),
     las= 1,
     xlim=c(-4,8),
     ylim=c(-4,8), col= "gray65",
     main = "log2 before/after correction")

fit<-(lm(t_test$log2TvsN_before_batch_cor~t_test$log2TvsN_batch_cor))
legend("topleft", bty="n", cex=2, legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=2)))
points(sig_both$log2TvsN_before_batch_cor,sig_both$log2TvsN_batch_cor,col= "#2668A8",pch=16, cex= 1.2)
points(sig_before$log2TvsN_before_batch_cor,sig_before$log2TvsN_batch_cor,col= "#d95f02",pch=16, cex= 1.2)
points(sig_after$log2TvsN_before_batch_cor,sig_after$log2TvsN_batch_cor,col= "#1b9e77",pch=16, cex= 1.2)
legend ("bottomright", pch=c(16,16), c("significance", paste0( "before (", nrow(sig_before), ")"), 
                                       paste0( "after (", nrow(sig_after), ")"), 
                                       paste0( "always (", nrow(sig_both), ")"),
                                       paste0( "never (", nrow(sig_never), ")")), 
        bty = "n",cex=1.5,
        col=c("white", "#d95f02","#1b9e77","#2668A8","gray65")) 

dev.off()

## Figure 7b -----
## ridgeline plot mean readcounts before and after selection 
# read expression values before correction
exprs_df <- read.csv2(file.path(conf$output_dir, "TCGA-LUSC_before_batch_cor_include.all.samplesTRUE_rpm_median15_exprs_df.csv"))

exprs_df <- makeRn(exprs_df)
exprs_df$meanExp <- rowMeans(exprs_df)
exprs_df <- make1col(exprs_df)
names(exprs_df)[1] <- "id"
exprs_df$id <- gsub('-|\\|', '.', exprs_df$id)

# Coerce a dataframe containing values (doesn't matte which ones) when condition is significant
expr_sign <- join(exprs_df[c("id", "meanExp")], sig_after[c("id", "t_before_batch_cor")], type= "full")
expr_sign <- join(expr_sign, sig_before[c("id", "p_before_batch_cor")], type= "full")
expr_sign <- join(expr_sign, sig_both[c("id", "t_batch_cor")], type= "full")
expr_sign <- join(expr_sign, sig_never[c("id", "p_batch_cor")], type= "full")
names(expr_sign)[3:6] <- c("sign_after","sign_before", "sign_both", "sig_never" )

# get expression values for the isomiRs when significant
for (row in 1:nrow(expr_sign)){
  if (!is.na(expr_sign[row, "sign_after"])){
    expr_sign[row, "sign_after"]<- expr_sign$meanExp[row]
  } else if (!is.na(expr_sign[row, "sign_before"])){
    expr_sign[row, "sign_before"]<- expr_sign$meanExp[row]
  }else if (!is.na(expr_sign[row, "sign_both"])){
    expr_sign[row, "sign_both"]<- expr_sign$meanExp[row]
  }else if (!is.na(expr_sign[row, "sig_never"])){
    expr_sign[row, "sig_never"]<- expr_sign$meanExp[row]
  }
}

expr_sign <- expr_sign[order(expr_sign$meanExp),]
expr_sign<-(t(makeRn(expr_sign)))

exp_df <- data.frame(rep(c("1never", "3after", "4before", "2both"), # names are plotted according to alphabet
                         c(ncol(expr_sign),ncol(expr_sign),ncol(expr_sign),ncol(expr_sign))), 
                     c(expr_sign[5,],expr_sign[2,],expr_sign[3,],expr_sign[4,]), 
                     rep(c("gray65","#2668A8", "#d95f02", "#1b9e77"), 
                         c(ncol(expr_sign),ncol(expr_sign),ncol(expr_sign),ncol(expr_sign))))
names(exp_df) <- c("type", "meanExpr", "colour")
exp_df <- na.omit(exp_df)
exp_df[,2]<- log2(exp_df[,2])

my_scale <- scale_fill_manual( values = exp_df$colour) 

pdf(file.path(conf$output_dir, "Figure7b.pdf"))
ggplot(exp_df, aes(x = meanExpr , y = type, fill = type)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")+
  scale_fill_manual( values = c("gray65","#2668A8", "#1b9e77" , "#d95f02" ))+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 30),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size= 25))

dev.off()

## Figure 7c -----
## Bubbleheatmap showing log2 fold changes and p-values of significant samples after batch correction
isos <- sig_p_after
isos <- isos[,c(1,3,4,7,8)]

# Add known functions of the isomiRs (derived from literature search)
names_isos <- c("hsa.let.7f.5p.0.0.", "hsa.let.7f.5p.0.1.", "hsa.miR.10a.5p.0.0.", "hsa.miR.10a.5p.0.1.", 
                "hsa.miR.145.5p.0.1.", "hsa.miR.146b.5p.0.1.", "hsa.miR.15a.5p.0.2.", "hsa.miR.199a.5p.0..3.", 
                "hsa.miR.26b.5p.0.0.", "hsa.miR.27a.3p.0..3.", "hsa.miR.27b.3p.0..3.", "hsa.miR.29a.3p..1..2.", 
                "hsa.miR.30c.5p.0.0.", "hsa.miR.30e.5p.0..3.", "hsa.miR.342.3p.0.0.", "hsa.miR.92a.3p.0..1.", 
                "hsa.miR.92b.3p.0..2.","hsa.miR.99b.5p.1..1.")
function_isos <- c("TSG in various entities", "TSG in various entities", "not classified", "not classified", 
                   "TSG in lung cancer", "TSG in lung cancer", "not classified", "TSG in lung cancer", 
                   "TSG in various entities", "OG in lung cancer", "TSG in various entities", "not classified", 
                   "TSG in lung cancer",  "TSG in lung cancer", "TSG in lung cancer", "OG in lung cancer", 
                   "OG in lung cancer", "not classified")



function_isos <- data.frame(names_isos, function_isos)
names(function_isos)[1] <-"id"

isos <- join(isos, function_isos, type = "left")
isos <- isos[order(isos$function_isos),]


# introduction of fake row to control bubble size
fake <- c("fake", 1.2, 0.123, 1.2, 0.123, "not available")
isos <- rbind(isos, c("fake", 1.2, 0.123, 1.2, 0.123, "not available"))
names(isos) <- c("ids", "p", "log2", "p", "log2", "fct")

(bubble_struct <- rbind(isos[,c(1:3,6)], isos[,c(1,4:6)]))
bubble_struct$x <- c(rep(0.5, nrow(isos)), rep(1.5, nrow(isos)))
bubble_struct$y <- c(rep(0.5:(nrow(isos)-0.5),2))

bubble_struct$p <- as.numeric(as.character(bubble_struct$p))
bubble_struct$log2 <- as.numeric(as.character(bubble_struct$log2))

bubble_struct$Bubblesize <- NA


for (row in 1:nrow(bubble_struct)){
  if (bubble_struct$p[row]> 1 ){ bubble_struct$Bubblesize[row] = 1}
  else if (bubble_struct$p[row]> 0.05 & bubble_struct$p[row]<=1){ bubble_struct$Bubblesize[row] = 2}
  else if (bubble_struct$p[row]< 0.05 & bubble_struct$p[row]>=0.001){ bubble_struct$Bubblesize[row] = 3}
  else if (bubble_struct$p[row]< 0.001 ){ bubble_struct$Bubblesize[row] = 4}
}

col <- colorRampPalette(c("blue", "white", "red"))
colorful <- col(900)

p <- ggplot(bubble_struct, aes(x = x, y = y, size = Bubblesize, fill = log2)) +
  geom_point(shape = 21,na.rm = TRUE, show.legend= TRUE)+ 
  theme_bw()+ theme(aspect.ratio=1)+
  scale_x_continuous(breaks = seq(0, 30, 60)) +
  scale_y_continuous(breaks = seq(0, 30, 60)) +
  scale_fill_gradientn(colors= colorful,values= rescale(c(-1,0,1)),guide = "colorbar", limits=c(-1,1))
p

ggsave(plot=p,height=6,width=6,dpi=200, filename=file.path(conf$output_dir,"Figure7c.pdf"), useDingbats=FALSE)
