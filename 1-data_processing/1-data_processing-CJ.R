#!/netscratch/dep_psl/grp_psl/ThomasN/tools/bin/bin/Rscripts

# R script to process each dataset to calculate logFC by microbial treatments compared to axenic control
# Originally by Ryohei Thomas Nakano; nakano@mpipz.mpg.de
# 02 Feb 20201

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

target   <- "CJ"
seq_type <- "PE"

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

# rename names
design$genotype <- str_replace(design$genotype, "Col-0", "Col")
design$treatment_1 <- str_replace(design$treatment_1, "B-SynCom", "B_SynCom")

# model matrix
Rep   <- factor(design$rep)
group <- factor(paste(design$nutrient, design$genotype, design$treatment_1, sep="_"))

model <- model.matrix( ~ 0 + group + Rep)
colnames(model) <- str_replace(colnames(model), "group", "")

# model <- model[, -ncol(model)]
contrasts <- makeContrasts(

	CJ_normal_Col_BSynCom  = (   Normal_Col_B_SynCom  - Normal_Col_axenic   ),
	CJ_normal_f6h1_BSynCom = (  Normal_f6h1_B_SynCom  - Normal_f6h1_axenic  ),
	CJ_lowFe_Col_BSynCom   = (   low_Fe_Col_B_SynCom  - low_Fe_Col_axenic   ),
	CJ_lowFe_f6h1_BSynCom  = (  low_Fe_f6h1_B_SynCom  - low_Fe_f6h1_axenic  ),

	levels=model)

# ============================= #


source(paste(scripts, "1-GLM_run.R", sep=""))





