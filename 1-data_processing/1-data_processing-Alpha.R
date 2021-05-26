#!/netscratch/dep_psl/grp_psl/ThomasN/tools/bin/bin/Rscripts
# R script to plot all RNAseq data on my hands
# 29 Dec 2020
# Ryohei Thomas Nakano; nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(edgeR,    quietly=T, warn.conflicts=F)
library(stringr,  quietly=T, warn.conflicts=F)
library(dplyr,    quietly=T, warn.conflicts=F)

source("/biodata/dep_psl/grp_psl/ThomasN/scripts/ggplot-themes_RTN.R")

# paths
dat_dir        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/201210all/original_data/"

processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/processed_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/fig/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/statistics/"
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/scripts/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# ========= paramters ========= #

target   <- "Alphaproteobacteria"
seq_type <- "SE"

# ============================= #

# import 
dat    <- read.table(paste(dat_dir, "featureCounts_table-", seq_type, ".txt", sep=""), header=T, row.names=1, sep="\t", comment.char="#", stringsAsFactors=F)
design <- read.table(paste(dat_dir, "design.txt", sep=""), header=T, stringsAsFactors=F, sep="\t")

colnames(dat) <- str_replace(colnames(dat), "_hisat2_sorted.bam", "") %>% str_replace(., ".*hisat2_map\\.", "")

# filter to target samples 
idx <- design$project == target & design$included
design <- design[idx,]

idx <- match(design$Library, colnames(dat))
dat <- as.matrix(dat[, idx])





# ========= paramters ========= #

# model matrix
Rep   <- factor(design$rep)
group <- factor(design$treatment_1)

model <- model.matrix( ~ 0 + group + Rep)
colnames(model) <- str_replace(colnames(model), "group", "")

# model <- model[, -ncol(model)]
contrasts <- makeContrasts(

	Alpha_R129_E   =  (   R129_E - axenic  ),
	Alpha_R13_A    =  (    R13_A - axenic  ),
	Alpha_Root142  =  (  Root142 - axenic  ),
	Alpha_Root491  =  (  Root491 - axenic  ),
	Alpha_Root1497 =  ( Root1497 - axenic  ),
	Alpha_Root700  =  (  Root700 - axenic  ),

	levels=model)

# ============================= #


source(paste(scripts, "1-GLM_run.R", sep=""))







