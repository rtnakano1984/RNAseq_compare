#!/netscratch/dep_psl/grp_psl/ThomasN/tools/bin/bin/Rscripts

# R script for comparing differenrt RNAseq datasets
# Originally by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de
# 02 Feb 2020

# logFC tables were created by scripts/compare_dataset/Fig2-*_data_processing.R

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(lsa,           quietly=T, warn.conflicts=F)
library(reshape2,      quietly=T, warn.conflicts=F)
library(patchwork,     quietly=T, warn.conflicts=F)
library(ggplot2,       quietly=T, warn.conflicts=F)
library(stringr,       quietly=T, warn.conflicts=F)
library(ggrepel,       quietly=T, warn.conflicts=F)
library(psych,         quietly=T, warn.conflicts=F)
library(dplyr,         quietly=T, warn.conflicts=F)
library(vegan,         quietly=T, warn.conflicts=F)
library(ape,           quietly=T, warn.conflicts=F)
library(ggtree,        quietly=T, warn.conflicts=F)
library(tidytree,      quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))

system(paste("mkdir -p", fig))

# functions
variability_table <- function(cca){

    chi <- c(cca$tot.chi,
                   cca$CCA$tot.chi, cca$CA$tot.chi)
    variability_table <- cbind(chi, chi/chi[1])
    colnames(variability_table) <- c("inertia", "proportion")
    rownames(variability_table) <- c("total", "constrained", "unconstrained")
    return(variability_table)

}

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

# === create DEG_list === #
# P values
idx <- str_detect(logFC_P_melt$Var2, "_PValue")
P_melt <- logFC_P_melt[idx, ]
P_melt$Var2 <- str_replace(P_melt$Var2, "_PValue", "")

idx <- P_melt$Var2 %in% design$name
P_melt <- P_melt[idx,]

# sort
idx <- order(P_melt$Var2, P_melt$Var1)
P_melt <- P_melt[idx,]

# logFCs
idx <- str_detect(logFC_P_melt$Var2, "_logFC")
logFC_melt <- logFC_P_melt[idx,]
logFC_melt$Var2 <- str_replace(logFC_melt$Var2, "_logFC", "")

idx <- logFC_melt$Var2 %in% design$name
logFC_melt <- logFC_melt[idx,]

# sort
idx <- order(logFC_melt$Var2, logFC_melt$Var1)
logFC_melt <- logFC_melt[idx,]



# significance threshold
idx <- (abs(logFC_melt$value) > log2(FC_threshold)) & (P_melt$value < alpha) 
# idx <- (P_melt$value < alpha) 
logFC_melt$sig <- as.numeric(idx)
P_melt$sig     <- as.numeric(idx)
# list of DEGs
DEG_list <- split(as.character(P_melt$Var1[idx]), P_melt$Var2[idx])

# DEGs
DEGs <- sort(unique(unlist(DEG_list)))

# DEGs without ShijiHou datasets
idx <- str_detect(names(DEG_list), "^Shiji")
DEG_Shiji <- sort(unique(unlist(DEG_list[!idx])))

write.table(DEGs, file=paste(stat, "DEGs.txt", sep=""), col.names=F, row.names=F, sep="\n", quote=F)

count <- as.data.frame(sapply(DEG_list, length))
write.table(count, file=paste(stat, "DEG_count.txt", sep=""), col.names=F, row.names=T, sep="\t", quote=F)

idx <- match(design$name, rownames(count))
design$count <- count[idx, 1]
design$label <- paste(design$short, " (", design$count, ")", sep="")

# === create logFC matrix, replacing missing values with 0 === #
# Missing genes are likely due to the low expression levels #

# dcast
logFC_dcast <- dcast(logFC_melt, Var1 ~ Var2)

# replace NA with 0
idx <- is.na(logFC_dcast)
logFC_dcast[idx] <- 0

# convert to simple dcast matrix (mat) and melted matrix (melt)
mat <- as.matrix(logFC_dcast[, -1])
rownames(mat) <- logFC_dcast$Var1

melt <- melt(mat)

# export logFC table
write.table(mat, paste(processed_data, "logFC_matrix.txt", sep=""), quote=F, sep="\t", col.names=NA, row.names=T)


# === filter to DEGs === #
idx <- rownames(mat) %in% DEGs
mat_DEG <- mat[idx,]

idx <- melt$Var1 %in% DEGs
melt_DEG <- melt[idx,]


# === filter to DEGs without Shiji === #
idx <- str_detect(colnames(mat), "^Shiji")
mat_Shiji <- mat[, !idx]

idx <- rownames(mat_Shiji) %in% DEG_Shiji
mat_Shiji <- mat_Shiji[idx,]

idx <- melt$Var1 %in% DEG_Shiji & melt$Var2 %in% colnames(mat_Shiji)
melt_Shiji <- melt[idx,]







# pcoa ###################################
# sort, just in case
idx <- match(design$name, colnames(mat))
mat <- mat[, idx]

# cor and pcoa
cor <- cosine(mat)

source(paste(scripts, "2-PCoA_run.R", sep=""))

ggsave(p3, file=paste(fig, "pcoa/pcoa-colour_lifestyle.pdf", sep=""), width=10, height=4.5, bg="transparent")
ggsave(p6, file=paste(fig, "pcoa/pcoa-colour_project.pdf", sep=""), width=10, height=4.5, bg="transparent")
ggsave(p9, file=paste(fig, "pcoa/pcoa-colour_treatment_1.pdf", sep=""), width=10, height=4.5, bg="transparent")

# statistics
d <- as.dist(1-cor)
set.seed(0)
# adonis2 <- adonis2(d ~ project + time + nutrient + lifestyle, design, sqrt.dist=T, by="term")
# adonis2 <- adonis2(d ~ system + time + nutrient + lifestyle, design, sqrt.dist=T, by="term")
adonis2 <- adonis2(d ~ system + time + nutrient + kingdom + lifestyle + project, design, sqrt.dist=T, by="term")

sink(paste(stat, "adonis2_var-all.txt", sep=""))
adonis2
sink()

# plot adnois2 results
df <- as.data.frame(adonis2)
residual <- df$R2[nrow(df)-1]

df <- df[1:(nrow(df)-2), ]

df$term <- factor(rownames(df), levels=rev(rownames(df)))
df$sig  <- as.numeric(df$'Pr(>F)' < alpha)

p <- ggplot(df, aes(y=term, x=R2*100, color=factor(sig))) +
	geom_bar(stat="identity", size=1) +
	scale_colour_manual(values=c("0"=NA, "1"="black"), guide=F) +
	labs(y="", x="Explained variance (%)", title=paste("Residual R^2: ", format(residual*100, digits=4), "%", sep="")) +
	theme_RTN +
	theme(
		title=element_text(size=7),
		axis.text.x=element_text(size=10, color="black"),
		axis.text.y=element_text(size=10, color="black"),
		axis.title.x=element_text(size=12),
		axis.title.y=element_text(size=12)
		)
ggsave(p, file=paste(fig, "adonis2_var-all-bar.pdf", sep=""), width=3, height=3, bg="transparent")

# statistics
d <- as.dist(1-cor)
set.seed(0)
# adonis2 <- adonis2(d ~ project + time + nutrient + lifestyle, design, sqrt.dist=T, by="term")
# adonis2 <- adonis2(d ~ system + time + nutrient + lifestyle, design, sqrt.dist=T, by="term")
adonis2 <- adonis2(d ~ system + time + nutrient + kingdom + patho + project, design, sqrt.dist=T, by="term")

sink(paste(stat, "adonis2_var_patho-all.txt", sep=""))
adonis2
sink()

# plot adnois2 results
df <- as.data.frame(adonis2)
residual <- df$R2[nrow(df)-1]

df <- df[1:(nrow(df)-2), ]

df$term <- factor(rownames(df), levels=rev(rownames(df)))
df$sig  <- as.numeric(df$'Pr(>F)' < alpha)

p <- ggplot(df, aes(y=term, x=R2*100, color=factor(sig))) +
	geom_bar(stat="identity", size=1) +
	scale_colour_manual(values=c("0"=NA, "1"="black"), guide=F) +
	labs(y="", x="Explained variance (%)", title=paste("Residual R^2: ", format(residual*100, digits=4), "%", sep="")) +
	theme_RTN +
	theme(
		title=element_text(size=7),
		axis.text.x=element_text(size=10, color="black"),
		axis.text.y=element_text(size=10, color="black"),
		axis.title.x=element_text(size=12),
		axis.title.y=element_text(size=12)
		)
ggsave(p, file=paste(fig, "adonis2_var_patho-all-bar.pdf", sep=""), width=3, height=3, bg="transparent")



# cpcoa
capscale <- capscale(d ~ treatment_1 + nutrient + time + light + Condition(system + project), data=design, add=FALSE, sqrt.dist=TRUE)


# ANOVA-like permutation analysis
set.seed(0)
perm_anova <- anova.cca(capscale)


				    
# generate variability tables and calculate confidence intervals for the variance
var_tbl <- variability_table(capscale)
eig <- capscale$CCA$eig
variance <- var_tbl["constrained", "proportion"]
p.val <- perm_anova[1, 4]

# extract the weighted average (sample) scores
points <- capscale$CCA$wa[, 1:3]
points <- as.data.frame(points)
colnames(points) <- c("x", "y", "z")

points <- cbind(points, design[match(rownames(points), design$name), ])
points$lifestyle   <- factor(points$lifestyle,   levels=lifestyles$names)
# points$project     <- factor(points$project,     levels=project$names)
points$treatment_1 <- factor(points$treatment_1, levels=treatment_1$names)

# by lifestyle
p1 <- ggplot(points, aes(x=x, y=y, colour=lifestyle, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=lifestyles$colours[-1]) +
	labs(x=paste("CAP1: ", format(eig[1]/sum(eig)*100, digits=4), "%", sep=""),
		 y=paste("CAP2: ", format(eig[2]/sum(eig)*100, digits=4), "%", sep=""),
		 colour="", label="") +
	ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
		      format(p.val, digits=2),
		      sep="")) +
	theme_RTN +
	theme(legend.position="none")

p2 <- ggplot(points, aes(x=x, y=z, colour=lifestyle, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=lifestyles$colours[-1]) +
	labs(x=paste("CAP1: ", format(eig[1]/sum(eig)*100, digits=4), "%", sep=""),
		 y=paste("CAP3: ", format(eig[3]/sum(eig)*100, digits=4), "%", sep=""),
		 colour="", label="") +
	ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
		      format(p.val, digits=2),
		      sep="")) +
	theme_RTN +
	theme(legend.position="right")

p3 <- p1 + p2 & theme(plot.background=element_blank())
ggsave(p3, file=paste(fig, "pcoa/cpcoa_treatment_1_nut_time_light-colour_lifestyle.pdf", sep=""), width=10, height=4.5, bg="transparent")



# by microbial inoculum
p7 <- ggplot(points, aes(x=x, y=y, colour=treatment_1, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=treatment_1$colours) +
	labs(x=paste("CAP1: ", format(eig[1]/sum(eig)*100, digits=4), "%", sep=""),
		 y=paste("CAP2: ", format(eig[2]/sum(eig)*100, digits=4), "%", sep=""),
		 colour="", label="") +
	ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
		      format(p.val, digits=2),
		      sep="")) +
	theme_RTN +
	theme(legend.position="none")

p8 <- ggplot(points, aes(x=x, y=z, colour=treatment_1, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=treatment_1$colours) +
	labs(x=paste("CAP1: ", format(eig[1]/sum(eig)*100, digits=4), "%", sep=""),
		 y=paste("CAP3: ", format(eig[3]/sum(eig)*100, digits=4), "%", sep=""),
		 colour="", label="") +
	ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
		      format(p.val, digits=2),
		      sep="")) +
	theme_RTN +
	theme(legend.position="right")

p9 <- p7 + p8 & theme(plot.background=element_blank())
ggsave(p9, file=paste(fig, "pcoa/cpcoa_treatment_1_nut_time_light-colour_treatment_1.pdf", sep=""), width=10, height=4.5, bg="transparent")








# pcoa - DEGs ##################################
cor <- cosine(mat_DEG)

source(paste(scripts, "2-PCoA_run.R", sep=""))

ggsave(p3, file=paste(fig, "pcoa/pcoa-colour_lifestyle.DEG.pdf", sep=""), width=10, height=4.5, bg="transparent")
ggsave(p6, file=paste(fig, "pcoa/pcoa-colour_project.DEG.pdf", sep=""), width=10, height=4.5, bg="transparent")
ggsave(p9, file=paste(fig, "pcoa/pcoa-colour_treatment_1.DEG.pdf", sep=""), width=10, height=4.5, bg="transparent")



# pcoa - DEGs - without Shiji ##################################
cor <- cosine(mat_Shiji)

source(paste(scripts, "2-PCoA_run.R", sep=""))

ggsave(p3, file=paste(fig, "pcoa/pcoa-colour_lifestyle.DEG.agar.pdf", sep=""), width=10, height=4.5, bg="transparent")
ggsave(p6, file=paste(fig, "pcoa/pcoa-colour_project.DEG.agar.pdf", sep=""), width=10, height=4.5, bg="transparent")
ggsave(p9, file=paste(fig, "pcoa/pcoa-colour_treatment_1.DEG.agar.pdf", sep=""), width=10, height=4.5, bg="transparent")




# ================== copmpare by pairwise PCC - DEG ================== #
n <- nrow(design)

# cor tests
cor_test <- lapply(1:n, function(i){
	lapply(1:n, function(j){

		data_1 <- design$name[i]
		data_2 <- design$name[j]

		# DEGs based on data_1
		genes <- DEG_list[[which(names(DEG_list) == data_1)]]

		# extract logFC
		logFC_1 <- mat[rownames(mat) %in% genes, data_1]
		logFC_2 <- mat[rownames(mat) %in% genes, data_2]

		# cor test
		cor <- cosine(logFC_1, logFC_2)

		out <- data.frame(data_1=data_1, data_2=data_2, PCC=cor)	
	}) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)


# sort for plotting by hclust
# dcast
dcast <- dcast(cor_test, data_1 ~ data_2, value.var="PCC")
rownames(dcast) <- dcast[, 1]
mat <- as.matrix(dcast[, -1])
diag(mat) <- 1

# mean of data_1 -> data_2 and data_2 -> data_1 to make it a symmetric dissimilarity matrix
for(i in 1:n){
	for(j in 1:n){
		mean <- (mat[i,j] + mat[j,i])/2
		mat[i,j] <- mean
		mat[j,i] <- mean
	}
}

# hclust and sort
hclust <- hclust(as.dist(1-mat), "average")
sorted <- rownames(mat)[hclust$order]

cor_test$data_1 <- factor(cor_test$data_1, levels=sorted)
cor_test$data_2 <- factor(cor_test$data_2, levels=rev(sorted))





# sort for plotting
idx <- match(cor_test$data_1, design$name)
cor_test$data_1_lifestyles <- factor(design$lifestyle[idx], levels=lifestyles$names[-1])

idx <- match(cor_test$data_2, design$name)
cor_test$data_2_lifestyles <- factor(design$lifestyle[idx], levels=lifestyles$names[-1])








# ============== PLOTTING EXCERCISE ============== #
# # simple heatmap
# p_heat <- ggplot(cor_test, aes(x=data_2, y=data_1_lab, fill=PCC)) +
# 	geom_tile(colour="black", alpha=1) +
# 	scale_fill_gradient2(low="magenta4", mid="white", high="green4", midpoint=0, guide="legend", na.value=c_grey, breaks=seq(-1, 1, .25)) +
# 	scale_size_manual(values=c("0"=1, "1"=2)) +
# 	theme_RTN +
# 	theme(
# 		axis.text.x=element_text(size=7, angle=90, hjust=1,  vjust=.5),
# 		axis.text.y=element_text(size=7, angle=0,  hjust=1, vjust=.5),
# 		axis.line.x=element_blank(),
# 		axis.line.y=element_blank(),
# 		axis.ticks.x=element_blank(),
# 		axis.ticks.y=element_blank(),
# 		legend.position="top") +
# 	labs(y="Query dataset (DEGs)", x="Target dataset", filL="PCC")
# ggsave(p_heat, file=paste(fig, "PCC_heatmap.pdf", sep=""), width=12, height=13.5, bg="transparent")



# hclust dendrogram, color by lifestyle
group <- split(design$name, design$lifestyle)

tree <- as.phylo(hclust)
write.tree(tree, file=paste(processed_data, "dendrogram.tree", sep=""))

tree <- groupOTU(tree, group, 'Lifestyle')

color_label <- lifestyles$colours
names(color_label) <- lifestyles$names

idx <- match(tree$tip.label, design$name)
tree$tip.label <- design$label

# tree plot #####
p_tree <- ggtree(tree, aes(color=Lifestyle), size=.8) +
	scale_color_manual(values=color_label) +
	geom_tiplab(size=5, offset=.1, align=T, linetype="blank") +
	theme_tree() +
	theme(legend.position="none",
		plot.margin=margin(t=5.5, b=5.5, r=8, l=0, unit="pt"))



# for y axis alignment by ylim2()
xy <- filter(p_tree, isTip)[, c("label", "x", "y")]

x_sorted <- xy$label[order(xy$y, decreasing=T)]

idx <- match(x_sorted, design$label)
x_sorted       <- design$name[idx]
x_sorted_short <- design$short[idx]

cor_test$data_2 <- factor(cor_test$data_2, levels=x_sorted)
levels(cor_test$data_2) <- x_sorted_short

idx <- match(cor_test$data_1, design$name)
cor_test$data_1_lab <- design$label[idx]

idx <- match(cor_test$data_1_lab, xy$label)
cor_test$y <- xy$y[idx]

x_max <- max(xy$x)
p_tree <- p_tree + xlim(NA, x_max + 1.2)




# heatmap plot #####
p_heat <- ggplot(cor_test, aes(x=data_2, y=y, fill=PCC)) +
	geom_tile(colour="black", alpha=1) +
	scale_fill_gradient2(low="magenta4", mid="white", high="green4", midpoint=0, guide="legend", na.value=c_grey, breaks=seq(-1, 1, .25)) +
	scale_size_manual(values=c("0"=1, "1"=2)) +
	theme_RTN +
	theme(
		axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
		axis.text.y=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		legend.position="top",
		plot.margin=margin(t=5.5, b=5.5, r=0, l=0, unit="pt")) +
	labs(y="", x="", fill="Cosine Similarity") +
	ylim2(p_tree)

# point plots
p_point <- ggplot(cor_test, aes(y=y, x=PCC, colour=data_2_lifestyles)) +
	geom_point(size=1) +
	geom_vline(xintercept=0, colour=c_grey, linetype="solid") +
	scale_colour_manual(values=lifestyles$colours[-1], guide=F) +
	scale_x_continuous(position="top") +
	labs(x="Cosine Similarity", y="") +
	theme_RTN +
	theme(
		axis.text.x=element_text(size=10),
		axis.title.x=element_text(size=12),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.line.y=element_blank(),
		plot.margin=margin(t=5.5, b=5.5, r=2, l=0, unit="pt")) +
	ylim2(p_tree) 

# patchwork
p <- p_tree + p_heat + p_point + plot_layout(ncol=3, width=c(1.5, 3, .5)) & theme(plot.background=element_blank())
ggsave(p, file=paste(fig, "PCC.dendro_heat_point.pdf", sep=""), width=18, height=12, bg="transparent")




# boxplot
# group by comparison
cor_test$group <- NA

idx <- cor_test$data_1_lifestyles != "pathogenic" & cor_test$data_2_lifestyles != "pathogenic"
cor_test$group[idx] <- "Within non-pathogenic"

idx <- cor_test$data_1_lifestyles == "pathogenic" & cor_test$data_2_lifestyles == "pathogenic"
cor_test$group[idx] <- "Within pathogenic"

idx <- is.na(cor_test$group)
cor_test$group[idx] <- "Between pathogenic and non-pathogenic"

levels(cor_test$group) <- c("Within non-pathogenic", "Within pathogenic", "Between pathogenic and non-pathogenic")

# color by shared lifestyle
cor_test$lifestyle <- "0"

idx <- cor_test$data_1_lifestyles == cor_test$data_2_lifestyles
cor_test$lifestyle[idx] <- as.character(cor_test$data_1_lifestyles[idx])

cor_test$lifestyle <- factor(cor_test$lifestyle, levels=lifestyles$names)

# plot
p5 <- ggplot(cor_test, aes(x=group, y=PCC, colour=lifestyle)) +
	geom_point(position=position_jitter(width=.2)) +
	geom_boxplot(outlier.shape=NA, fill=c_white, colour=c_black) +
	scale_colour_manual(values=lifestyles$colours, labels=c("", lifestyles$names[-1])) +
	theme_RTN +
	theme(axis.text.x=element_text(angle=75, hjust=1, vjust=1),
		legend.position="right") +
	labs(x="", y="Cosine Similarity", colour="Within: ")
ggsave(p5, file=paste(fig, "PCC.boxplot.pdf", sep=""), width=4.5, height=6, bg="transparent")













