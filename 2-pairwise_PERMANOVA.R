#!/netscratch/dep_psl/grp_psl/ThomasN/tools/R403/bin/Rscript
# R script for comparing differenrt RNAseq datasets
# Originally by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de
# 02 Feb 2020

# logFC tables were created by scripts/compare_dataset/Fig2-*_data_processing.R

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(lsa,           quietly=T, warn.conflicts=F)
library(stringr,       quietly=T, warn.conflicts=F)
library(reshape2,      quietly=T, warn.conflicts=F)
library(RVAideMemoire, quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))

system(paste("mkdir -p", fig))

# =========== import and organize data =========== #
# import
design       <- read.table(paste(processed_data, "design.txt", sep=""), sep="\t", header=T, row.names=NULL, stringsAsFactors=F)

# filter out mutants
idx <- design$genotype %in% c("Col-0", "MYC2-FLAG")
design <- design[idx,]

logFC_P_melt <- lapply(list.files(processed_data, pattern="*logFC.P.txt"), function(x){
	temp <- as.matrix(read.table(paste(processed_data, x, sep=""), sep="\t", header=T, row.names=1, stringsAsFactors=F))
	temp_melt <- melt(temp)
	return(temp_melt)
}) %>% do.call(rbind, .)

# logFCs
idx <- str_detect(logFC_P_melt$Var2, "_logFC")
logFC_melt <- logFC_P_melt[idx,]
logFC_melt$Var2 <- str_replace(logFC_melt$Var2, "_logFC", "")

idx <- logFC_melt$Var2 %in% design$name
logFC_melt <- logFC_melt[idx,]

# sort
idx <- order(logFC_melt$Var2, logFC_melt$Var1)
logFC_melt <- logFC_melt[idx,]

# dcast
logFC_dcast <- dcast(logFC_melt, Var1 ~ Var2)

# replace NA with 0
idx <- is.na(logFC_dcast)
logFC_dcast[idx] <- 0

# convert to simple dcast matrix (mat)
mat <- as.matrix(logFC_dcast[, -1])
rownames(mat) <- logFC_dcast$Var1

# sort, just in case
idx <- match(design$name, colnames(mat))
mat <- mat[, idx]

# cor and pcoa
cor <- cosine(mat)

# pairwise permanova between lifestyles
d <- as.dist(1-cor)
set.seed(0)

pair_permanova <- pairwise.perm.manova(d, factor(design$lifestyle, levels=c("pathogenic", "SynCom", "beneficial", "commensal")), nperm=999, p.method="fdr", R2=T)
results <- cbind(melt(pair_permanova$R2.value), melt(pair_permanova$p.value))
results <- na.omit(results)[, c(2, 1, 3, 6)]
colnames(results) <- c("lifestyle_1", "lifestyle_2", "R2", "FDR")
write.table(results, file=paste(stat, "pairwise_PERMANOVA.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")











