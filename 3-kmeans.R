#!/netscratch/dep_psl/grp_psl/ThomasN/tools/bin/bin/Rscripts

# R script for comparing differenrt RNAseq datasets
# Originally by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de
# 02 Feb 2020

# logFC tables were created by scripts/compare_dataset/Fig2-*_data_processing.R

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(lsa,         quietly=T, warn.conflicts=F)
library(reshape2,    quietly=T, warn.conflicts=F)
library(patchwork,   quietly=T, warn.conflicts=F)
library(ggplot2,     quietly=T, warn.conflicts=F)
library(stringr,     quietly=T, warn.conflicts=F)
# library(ggrepel,     quietly=T, warn.conflicts=F)
# library(psych,       quietly=T, warn.conflicts=F)
library(dplyr,       quietly=T, warn.conflicts=F)
library(ape,         quietly=T, warn.conflicts=F)
library(ggtree,      quietly=T, warn.conflicts=F)
library(tidytree,    quietly=T, warn.conflicts=F)
library(sva,       quietly=T, warn.conflicts=F) # the devel version of sva package has to be installed from github ("zhangyuqing/sva-devel"), after normal installation from BiocManager

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/201210all/original_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))

# ==================== data preparation ==================== #
# load
design    <- read.table(file=paste(original_data, "design.txt", sep=""),        sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
DEGs      <- read.table(file=paste(stat, "DEGs.txt", sep=""),                   sep="\t", header=F, row.names=NULL, check.names=F, stringsAsFactors=F)$V1
logFC     <- read.table(file=paste(processed_data, "logFC_matrix.txt", sep=""), sep="\t", header=T, row.names=1,    check.names=F, stringsAsFactors=F) %>% as.matrix
tree      <-  read.tree(file=paste(processed_data, "dendrogram.tree", sep=""))
contrasts <- read.table(file=paste(processed_data, "design.txt", sep=""),       sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)

Alpha_cluster <- read.table(file="/biodata/dep_psl/grp_psl/ThomasN/RNAseq/project3598/statistics/k_means-zero_centered-9.sorted_hclust.txt", header=T, row.names=NULL, stringsAsFactors=F, sep="\t")

# refine design
design <- design[design$included,]

idx <- design$genotype %in% c("Col-0", "MYC2-FLAG")
design <- design[idx,]

idx <- contrasts$genotype %in% c("Col-0", "MYC2-FLAG")
contrasts <- contrasts[idx,]



# ==================== Prepare logFC melt table ==================== #
idx <- rownames(logFC) %in% DEGs
logFC <- logFC[idx,]

logFC_sat <- saturate(logFC, .05)

logFC_melt <- melt(logFC)

idx <- logFC_melt$Var2 %in% contrasts$name
logFC_melt <- logFC_melt[idx,]

logFC_melt$fill <- saturate(logFC_melt$value, .01)

idx <- match(logFC_melt$Var2, contrasts$name)
logFC_melt <- cbind(logFC_melt, contrasts[idx,])




# ==================== kemans ==================== #
# # best k
# min <- 1
# max <- 100
# range <- c(min:max)

# AIC <- t(sapply(range, function(x){

# 	set.seed(0)
# 	fit <- kmeans(logFC_sat, centers=x)

# 	m <- ncol(fit$centers)
# 	n <- length(fit$cluster)
# 	k <- nrow(fit$centers)
# 	D <- fit$tot.withinss
# 	return(c(
# 		AIC=(D + 2*m*k),
# 		BIC=(D + log(n)*m*k)))
# }))

# # plot AIC/BIC
# AIC <- data.frame(n=range, AIC)
# melt <- melt(AIC, id.vars="n")
# p <- ggplot(melt, aes(x=n, y=value, colour=variable)) +
# 	geom_line() +
# 	geom_point() +
# 	labs(colour="", x="Number of clusters", y="") +
# 	theme_RTN +
# 	theme(legend.position="top")
# ggsave(p, file=paste(fig, "logFC-AIC_BIC.pdf", sep=""), width=4.5, height=3.5)

# # best k
# k <- range[which.min(AIC$BIC)]
k <- 21
set.seed(0)
fit <- kmeans(logFC_sat, centers=k)


cluster <- data.frame(ID=rownames(logFC), Cluster=fit$cluster)


# merge kmeans info with the data
idx <- match(logFC_melt$Var1, cluster$ID)
logFC_melt$cluster <- cluster$Cluster[idx]


# ==================== Sort clusters for plotting ==================== #
# sort clusters
cluster_summary <- logFC_melt %>% group_by(Var2, cluster) %>% summarise(value=mean(value)) %>% data.frame(stringsAsFactors=F)

cor <- cosine(as.matrix(dcast(cluster_summary, Var2 ~ cluster)[,-1]))
d <- as.dist(1-cor)
hclust <- hclust(d, "ward.D2")
sorted_cl <- colnames(cor)[hclust$order]

# sort genes
cluster$Cluster <- factor(cluster$Cluster, levels=sorted_cl)

idx <- order(cluster$Cluster)
sorted_gene <- cluster$ID[idx]

# implement sorts
logFC_melt$Var1    <- factor(logFC_melt$Var1,    levels=sorted_gene)
logFC_melt$cluster <- factor(logFC_melt$cluster, levels=sorted_cl)

cluster_summary$cluster <- factor(cluster_summary$cluster, levels=sorted_cl)



# ==================== export clusters ==================== #
# sort
idx <- order(cluster$Cluster)
cluster <- cluster[idx,]

# export
write.table(cluster, file=paste(stat, "compare_dataset-logFC-kmeans_", k, "clusters.seed0.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")

# export csv for metascape
max_N <- max(table(cluster$Cluster))
metascape <- lapply(1:k, function(x){
	idx <- cluster$Cluster == x
	target <- as.character(cluster$ID[idx])
	out <- rep(NA, max_N)
	out[1:length(target)] <- target
	return(out)
}) %>% do.call(cbind, .)
colnames(metascape) <- paste("Cluster_", 1:k, sep="")
write.table(metascape, file=paste(stat, "compare_dataset-logFC-kmeans_", k, "clusters.seed0.for_metascape.csv", sep=""), sep=",", quote=F, col.names=T, row.names=F)





# ==================== Prepare plot for plotting k-means cluster from Alphaproteobacteria dataset ==================== #
# merge with kmeans info from Alphaproteobacteria dataset
idx <- match(logFC_melt$Var1, Alpha_cluster$ID)
logFC_melt$Alpha_cluster <- Alpha_cluster$Cluster[idx]

# profile
Alpha_cluster_profile <- dcast(logFC_melt, Var1 + cluster ~ Alpha_cluster, value.var="Alpha_cluster", length) %>% melt(id.var=c("Var1", "cluster"))

idx <- Alpha_cluster_profile$variable == "NA"
Alpha_cluster_profile <- Alpha_cluster_profile[!idx,]

# count
idx <- Alpha_cluster_profile$value > 0
Alpha_cluster_profile$value[idx] <- 1

count <- as.data.frame(table(Alpha_cluster_profile$cluster[idx]))
count$prop <- (count$Freq / sum(count$Freq))*100

idx <- order(count$prop, decreasing=T)
count <- count[idx,]

count$cum_prop <- cumsum(count$prop)
count$Var1 <- factor(count$Var1, levels=count$Var1)


p1 <- ggplot(count, aes(x=Var1, y=prop)) +
	geom_bar(stat="identity") +
	theme_RTN + 
	labs(x="Cluster", y="Proportion of Rhizobiales DEGs (%)")

p2 <- ggplot(count, aes(x=Var1, y=cum_prop, group=1)) +
	geom_line(colour=c_dark_red, size=.75) +
	geom_hline(yintercept=50, colour=c_dark_red, size=.5, linetype="dashed") +
	scale_y_continuous(position = "right") +
	theme_RTN + 
	theme(axis.text.y=element_text(colour=c_dark_red)) + 
	labs(x="Cluster", y="Cumulative sum of proportion")

p <- p1 + p2 + plot_layout(ncol=1, heights=c(1,1)) & theme(plot.background=element_blank())
ggsave(p, file=paste(fig ,"compare-kmeans-all21_vs_alpha9.pdf", sep=""), width=6, height=6, bg="transparent")



# factor
Alpha_cluster_profile$value    <- factor(Alpha_cluster_profile$value,    levels=c(0, 1))
Alpha_cluster_profile$variable <- factor(Alpha_cluster_profile$variable, levels=unique(Alpha_cluster$Cluster))
Alpha_cluster_profile$Var1     <- factor(Alpha_cluster_profile$Var1,     levels=sorted_gene)
Alpha_cluster_profile$cluster  <- factor(Alpha_cluster_profile$cluster,  levels=sorted_cl)

levels(Alpha_cluster_profile$variable) <- str_replace(levels(Alpha_cluster_profile$variable), "Cl_0", "")

# plot
p_alpha <- ggplot(Alpha_cluster_profile, aes(x=Var1, y=variable, fill=value)) +
	geom_tile() +
	scale_fill_manual(values=c("0"=NA, "1"=c_black), guide=F) +
	facet_grid(variable ~ cluster, drop=T, scales="free", space="free", switch="y") +
	theme_RTN +
	theme(strip.text.x=element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		panel.spacing=unit(.1, "lines"),
		plot.margin=margin(t=5.5, b=0, r=0, l=0, unit="pt"),
		panel.border=element_rect(fill=NA, size=.5, colour=c_black)) +
	labs(x="", y="k-means cluster in Figure 1")




# ==================== Prepare tree plott ==================== #
# dendrogram plot
group <- split(contrasts$name, contrasts$lifestyle)
tree <- groupOTU(tree, group, 'Lifestyle')

color_label <- lifestyles$colours
names(color_label) <- lifestyles$names

# tree plot #####
p_tree <- ggtree(tree, aes(color=Lifestyle), size=.8) +
	scale_color_manual(values=color_label) +
	geom_tiplab(size=3, offset=.1, align=T, linetype="blank") +
	theme_tree() +
	theme(legend.position="none",
		plot.margin=margin(t=5.5, b=5.5, r=5, l=0, unit="pt"))

# for y axis alignment by ylim2()
xy <- filter(p_tree, isTip)[, c("label", "x", "y")]

idx <- match(logFC_melt$Var2, xy$label)
logFC_melt$y <- xy$y[idx]

idx <- match(cluster_summary$Var2, xy$label)
cluster_summary$y <- xy$y[idx]

x_max <- max(xy$x)
p_tree <- p_tree + xlim(NA, x_max + 1.8)



# ==================== Plotting logFC heatmap ==================== #
# heatmap
p_heat <- ggplot(logFC_melt, aes(y=y, x=Var1, fill=fill)) +
	geom_tile(alpha=1, size=0) +
	scale_fill_gradient2(low=c_cudo_magenta, mid=c_white, high=c_very_dark_green, midpoint=0, breaks=seq(ceiling(min(logFC_melt$value)), floor(max(logFC_melt$value)), by=1), guide="legend") +
	# scale_size_manual(values=c("0"=.5, "1"=1)) +
	facet_grid(. ~ cluster, drop=T, scales="free", space="free") +
	theme_RTN +
	theme(legend.position="top",
		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=12),
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		panel.spacing=unit(.1, "lines"),
		plot.margin=margin(t=5.5, b=0, r=0, l=0, unit="pt"),
		panel.border=element_rect(fill=NA, size=.5, colour=c_black)) +
	ylim2(p_tree) +
	labs(x="", y="", fill="log2-FC") +
	guides(fill=guide_legend(override.aes=list(size=.5, colour=c_black)))



# heatmap, summarized by clusters
p_heat2 <- ggplot(cluster_summary, aes(y=y, x=cluster, fill=value)) +
	geom_tile(alpha=1, size=0) +
	scale_fill_gradient2(low=c_cudo_magenta, mid=c_white, high=c_very_dark_green, midpoint=0, breaks=seq(ceiling(min(logFC_melt$value)), floor(max(logFC_melt$value)), by=1), guide="legend") +
	# scale_size_manual(values=c("0"=.5, "1"=1)) +
	theme_RTN +
	theme(legend.position="top",
		strip.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
		axis.text.y=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		panel.spacing=unit(.1, "lines"),
		plot.margin=margin(t=5.5, b=0, r=0, l=0, unit="pt"),
		panel.border=element_rect(fill=NA, size=.5, colour=c_black)) +
	ylim2(p_tree) +
	labs(x="", y="", fill="log2-FC") +
	guides(fill=guide_legend(override.aes=list(size=.5, colour=c_black)))







# ==================== Plotting ==================== #

p <- p_tree + p_heat + plot_layout(widths=c(1,5)) & theme(plot.background=element_blank())
ggsave(p, file=paste(fig, "compare_dataset-logFC-kmeans_", k, "clusters.heatmap.png", sep=""), width=16, height=8, bg="transparent")

p <- plot_spacer() + p_alpha + p_tree + p_heat + plot_layout(ncol=2, nrow=2, heights=c(1,4), widths=c(1,5)) & theme(plot.background=element_blank())
ggsave(p, file=paste(fig, "compare_dataset-logFC-kmeans_", k, "clusters.heatmap_with_AlphaK.png", sep=""), width=16, height=10, bg="transparent")

p <- p_tree + p_heat2 + plot_layout(widths=c(1,2)) & theme(plot.background=element_blank())
ggsave(p, file=paste(fig, "compare_dataset-logFC-kmeans_", k, "clusters-summarized_heatmap.png", sep=""), width=9, height=8, bg="transparent")



# # ==================== After GO enrichment analysis by Metascape ==================== #

# # GO enrichment resulta manually curated from fig/metascape-logFC/Enrichment_heatmap/HeatmapSelectedGOTop100.pdf
# GO_df <- do.call(rbind, list(
# 	data.frame(GO="Cellular response to hypoxia",                   cluster=sorted_cl, idx=as.numeric(sorted_cl %in% c(2, 4, 9, 13, 5, 12, 11))),
# 	data.frame(GO="immune response",                                cluster=sorted_cl, idx=as.numeric(sorted_cl %in% c(4, 1, 14, 5))),
# 	data.frame(GO="response to toxic substances",                   cluster=sorted_cl, idx=as.numeric(sorted_cl %in% c(14, 16, 13, 11))),
# 	data.frame(GO="cellular detoxification",                        cluster=sorted_cl, idx=as.numeric(sorted_cl %in% c(16, 2))),
# 	data.frame(GO="negative regulation of cell cylce",              cluster=sorted_cl, idx=as.numeric(sorted_cl  ==  10)),
# 	data.frame(GO="meristem growth",                                cluster=sorted_cl, idx=as.numeric(sorted_cl  ==  14)),
# 	data.frame(GO="suberin biosynthetic process",                   cluster=sorted_cl, idx=as.numeric(sorted_cl  ==  16)),
# 	data.frame(GO="Carbon metabolism",                              cluster=sorted_cl, idx=as.numeric(sorted_cl %in% c(14, 16, 2))),
# 	data.frame(GO="Thioglycoside biosynthetic process",             cluster=sorted_cl, idx=as.numeric(sorted_cl  ==  20)),
# 	data.frame(GO="cellular response to oxgen-containing compound", cluster=sorted_cl, idx=as.numeric(sorted_cl %in% c(1, 13, 5))),
# 	data.frame(GO="reactive oxygen species metabolic process",      cluster=sorted_cl, idx=as.numeric(sorted_cl  ==  8))
# 	))
# GO_df$cluster <- factor(GO_df$cluster, levels=sorted_cl)
# GO_df$GO      <- factor(GO_df$GO,      levels=unique(GO_df$GO))

# p_GO <- ggplot(GO_df, aes(y=GO, x=cluster, fill=factor(idx))) +
# 	geom_tile(colour="black", size=0.6, alpha=1) +
# 	scale_fill_manual(values=c("1"="black", "0"="white"), labels=c("1"="+", "0"="-"), guide="legend") +
# 	theme_RTN +
# 	theme(
# 		legend.position="bottom",
# 		axis.text.y=element_text(colour="black", size=9),
# 		axis.text.x=element_blank(),
# 		legend.text=element_text(size=8),
# 		legend.title=element_text(size=10),
# 		axis.line.y=element_blank(),
# 		axis.line.x=element_blank(),
# 		axis.ticks=element_blank()) +
# 	labs(y="Enriched GOs", x="", fill="")


# p <- p_tree + p_heat2 + plot_spacer() + p_GO + plot_layout(ncol=2, widths=c(1,3), height=c(2, 1)) & theme(plot.background=element_blank())
# ggsave(p, file=paste(fig, "compare_dataset-logFC-kmeans_", k, "clusters-summarized_heatmap.GO_enrichment.png", sep=""), width=12, height=10, bg="transparent")




