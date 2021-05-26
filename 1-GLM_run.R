# A common script for GLM fitting and DEG picking
# Originally by Ryohei Thomas Nakano, nakano@mpipz.mpg.de
# 02 Feb 2021

# sort data table accorindg to AGI code
idx <- order(rownames(dat))
dat <- dat[idx,]

# count list for edgeR
y <- DGEList(counts=dat)

# normalized read counts
log2cpm <- cpm(y, prior.count=.25, log=TRUE)

# log2cpm
write.table(log2cpm, file=paste(processed_data, target, "-log2cpm.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)

# extract genes with a least 5 cpm over all samples AND expressed more than 5 counts at least in two samples 
idx <- (rowSums(2^log2cpm) > 10)  & (rowSums(2^log2cpm > 1) >= 2)
y <- y[idx, , keep.lib.sizes=F]

# TMM normalization
y <- calcNormFactors(y)

# estimate dispersion
y <- estimateGLMCommonDisp(y, model)
y <- estimateGLMTrendedDisp(y, model)
y <- estimateGLMTagwiseDisp(y, model)

# Generalized linear model fitting
fit <- glmFit(y, model)

# LRT contrasts
contrast_names <- colnames(contrasts)
n <- length(contrast_names)

# LRT for each contrasts
LRT.list <- lapply(1:n, function(x) glmLRT(fit, contrast=contrasts[,x]))
names(LRT.list) <- contrast_names

# logFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
	table <- LRT.list[[x]]$table[,c(1,4)]
	table$PValue <- p.adjust(table$PValue, method=p.adj.method)
	colnames(table) <- paste(contrast_names[x], colnames(table), sep="_")
	return(table)
	})
logFC_P <- do.call(cbind, logFC_P.list)
write.table(logFC_P, file=paste(processed_data, target, "-logFC.P.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)



# ##################### DEG extraction ############################

# # Significance picking for each tested model
DE.list <- lapply(1:n, function(x) decideTestsDGE(LRT.list[[x]], adjust.method=p.adj.method, p.value=alpha, lfc=log2(FC_threshold)))
names(DE.list) <- contrast_names

# # significance table
DE <- sapply(1:n, function(x) DE.list[[x]][,1])
colnames(DE) <- contrast_names
write.table(DE, file=paste(stat, target, "-significance_table.txt", sep=""), sep="\t", quote=T, row.names=T, col.names=NA)

# export log2cpm, logFC, pvalue of DEGs
idx <- rowSums(DE != 0) > 0
DEG <- rownames(DE)[idx]
write.table(DEG, file=paste(stat, target, "-DEG_list.txt", sep=""), sep="\n", quote=F, col.names=F, row.names=F)

