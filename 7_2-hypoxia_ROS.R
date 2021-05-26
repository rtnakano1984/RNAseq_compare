#!/netscratch/dep_psl/grp_psl/ThomasN/tools/bin/bin/Rscripts

# R script for comparing differenrt RNAseq datasets
# Originally by Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de
# 02 Feb 2020

# logFC tables were created by scripts/compare_dataset/Fig2-*_data_processing.R

options(warn=-1)

# clean up
rm(list=ls())

# packages
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
# library(org.At.tair.db,  quietly=T, warn.conflicts=F)

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/201210all/original_data/"
ref_data       <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/ref_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/fig/ind_gene_set/"
source(paste(scripts, "plotting_parameters.R", sep=""))

# ==================== data preparation ==================== #
# load
# design    <- read.table(file=paste(original_data, "design.txt", sep=""),        sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
# DEGs      <- read.table(file=paste(stat, "DEGs.txt", sep=""),                   sep="\t", header=F, row.names=NULL, check.names=F, stringsAsFactors=F)$V1
logFC     <- read.table(file=paste(processed_data, "logFC_matrix.txt", sep=""),   sep="\t", header=T, row.names=1,    check.names=F, stringsAsFactors=F) %>% as.matrix
contrasts <- read.table(file=paste(processed_data, "design.txt", sep=""),         sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
# goi       <- read.table(file=paste(processed_data, "goi.txt", sep=""),            sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
tree      <-  read.tree(file=paste(processed_data, "dendrogram.tree", sep=""))
cluster   <- read.table(file=paste(stat, "compare_dataset-logFC-kmeans_21clusters.seed0.txt", sep=""), sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F) 
go_list   <- read.table(file=paste(ref_data, "tair.gaf", sep=""),                 sep="\t", header=F, row.names=NULL, check.names=F, stringsAsFactors=F, comment.char="!", quote="")

# filter out mutants
idx <- contrasts$genotype %in% c("Col-0", "MYC2-FLAG")
contrasts <- contrasts[idx,]




# ==================== make a gene list ==================== #

target_id      <- c("GO:0001666", "GO:0000302", "GO:0071456", "GO:0034614", "GO:0042542", "GO:0098754", "GO:1990748", "GO:0009407")
target_names   <- c("hypoxia", "ROS", "ROS", "H2O2", "H2O2", "detoxification", "detoxification", "detoxification")
target_names_u <- unique(target_names)

idx <- str_detect(go_list$V10, "^AT")
go_list <- go_list[idx,]

goi_list <- lapply(target_id, function(x) {
	idx <- go_list$V5 == x
	genes <- go_list$V10[idx]
	Name  <- go_list$V3[idx]
	if(Name != genes) Name <- paste(genes, Name, sep="_")
	out_df <- data.frame(ID=genes, Name=Name, group=target_names[which(target_id==x)], Category="Hypoxia_ROS", stringsAsFactors=F)
	return(out_df)
}) %>% do.call(rbind, .) %>% unique

query <- sort(unique(goi_list$ID))

goi <- sapply(target_names_u, function(x){
	idx <- goi_list$group == x
	genes <- goi_list$ID[idx]

	idx <- query %in% genes
	out <- as.numeric(idx)
}) %>% data.frame(ID=query, ., stringsAsFactors=F)

idx <- goi$ID %in% cluster$ID
goi <- goi[idx,]

idx <- cluster$ID %in% goi$ID
cluster <- cluster[idx,]

sorted_genes <- goi$ID[hclust(dist(as.matrix(goi[, -1])), "ward.D2")$order]

idx <- match(goi$ID, goi_list$ID)
goi$Name <- goi_list$Name[idx]

idx <- goi$Name == ""
goi$Name[idx] <- goi$ID[idx]

idx <- match(goi$ID, cluster$ID)
goi$cluster <- cluster$Cluster[idx]

#
rownames(logFC) <- toupper(rownames(logFC))
goi$ID <- toupper(goi$ID)

n <- nrow(goi)



write.table(goi, file=paste(stat, "ROS_hypoxia_detox-gene_list.txt", sep=""), col.names=T, row.names=F, sep="\t", quote=F)


# ==================== Prepare logFC melt table ==================== #

idx <- match(goi$ID, rownames(logFC))
goi <- data.frame(goi, logFC[idx,], stringsAsFactors=F)

logFC_melt <- melt(goi, id.vars=c("ID", "Name", "cluster", target_names_u)) %>% melt(id.vars=c("ID", "Name", "cluster", "variable", "value"))
names(logFC_melt) <- c("ID", "Name", "cluster", "variable", "value", "GO", "GO_hit")

idx <- logFC_melt$GO_hit == 1
logFC_melt <- logFC_melt[idx,]

logFC_melt$fill <- saturate(logFC_melt$value, .01)

# idx <- match(logFC_melt$variable, contrasts$name)
# logFC_melt <- cbind(logFC_melt, contrasts[idx,])

logFC_melt$cluster <- factor(logFC_melt$cluster, levels=unique(cluster$Cluster))


logFC_melt$ID <- factor(logFC_melt$ID, levels=sorted_genes)

idx <- match(sorted_genes, goi$ID)
logFC_melt$Name <- factor(logFC_melt$Name, levels=goi$Name[idx])



# ==================== Prepare tree for plotting ==================== #
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

idx <- match(logFC_melt$variable, xy$label)
logFC_melt$y <- xy$y[idx]

x_max <- max(xy$x)
p_tree <- p_tree + xlim(NA, x_max + 1.8)



# # ==================== Plotting ==================== #
# go_color <- c(c_black, c_blue, c_red)
# names(go_color) <- target_names_u

# p_point <- ggplot(logFC_melt, aes(x=value, y=y, colour=fill)) +
# 	geom_point(position=position_jitter(width=.2)) +
# 	geom_vline(xintercept=0) +
# 	scale_colour_gradient2(low=c_cudo_magenta, mid="black", high=c_very_dark_green, midpoint=0, guide=F) +
# 	facet_grid(. ~ GO, scales="free", drop=T) +
# 	theme_RTN +
# 	ylim2(p_tree)
# p <- p_tree + p_point + plot_layout(widths=c(1,4)) & theme(plot.background=element_blank())
# ggsave(p, file=paste(fig, "hypoxia_ROS-point.png", sep=""), width=8, height=12, bg="transparent", limitsize=F)




# ==================== Plotting ==================== #
# GO strip
goi_melt <- melt(goi[, c("ID", "Name", "cluster", target_names_u)], id.vars=c("ID", "Name", "cluster"))
goi_melt$Name    <- factor(goi_melt$Name,    levels(logFC_melt$Name))
goi_melt$cluster <- factor(goi_melt$cluster, levels(logFC_melt$cluster))

p_go <- ggplot(goi_melt, aes(x=Name, y=variable, fill=factor(value))) +
	geom_tile(colour=NA, size=2) +
	facet_grid(. ~ cluster, space="free", scales="free", drop=T) +
	scale_fill_manual(values=c("0"="white", "1"="black"), guide=F) +
	theme_RTN +
	labs(x="", y="") +
	theme(axis.text.x=element_blank(),
		axis.text.y=element_text(colour="black"),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		panel.spacing=unit(.1, "lines"),
		plot.margin=margin(t=5.5, b=0, r=0, l=0, unit="pt"),
		panel.border=element_rect(fill=NA, size=.5, colour=c_black))
# ggsave(p_go, file=paste(fig, "hypoxia_ROS-strip.png", sep=""), width=sqrt(n)+4, height=8, bg="transparent", limitsize=F)

p_heat <- ggplot(unique(logFC_melt[, c("Name", "y", "fill", "cluster")]), aes(y=y, x=Name, fill=fill)) +
	geom_tile(alpha=1, size=0) +
	scale_fill_gradient2(low=c_cudo_magenta, mid=c_white, high=c_very_dark_green, midpoint=0, na.value=c_black, breaks=seq(ceiling(min(logFC_melt$value)), floor(max(logFC_melt$value)), by=1), guide="legend") +
	facet_grid(. ~ cluster, space="free", scales="free", drop=T) +
	# scale_size_manual(values=c("0"=.5, "1"=1)) +
	theme_RTN +
	theme(legend.position="right",
		axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=4),
		axis.text.y=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		strip.text.x=element_blank(),
		strip.background.x=element_blank(),
		panel.spacing=unit(.1, "lines"),
		plot.margin=margin(t=5.5, b=0, r=0, l=0, unit="pt"),
		panel.border=element_rect(fill=NA, size=.5, colour=c_black)) +
	ylim2(p_tree) +
	# xlim2(p_go) +
	labs(x="", y="", fill="log2-FC") +
	guides(fill=guide_legend(override.aes=list(size=.5, colour=c_black)))

p <- plot_spacer() + p_go + p_tree + p_heat + plot_layout(nrow=2, heights=c(.1, 1), widths=c(1,5)) & theme(plot.background=element_blank())
ggsave(p, file=paste(fig, "hypoxia_ROS-heatmap.png", sep=""), width=sqrt(n)+4, height=8, bg="transparent", limitsize=F)


















