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

#diretoreis
scripts        <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/scripts/"
stat           <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/statistics/"
processed_data <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/processed_data/"
original_data  <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/201210all/original_data/"
ref_data       <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/ref_data/"
fig            <- "/biodata/dep_psl/grp_psl/ThomasN/RNAseq2/210202all_logFC/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))

# ==================== data preparation ==================== #
# load
# design    <- read.table(file=paste(original_data, "design.txt", sep=""),        sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
# DEGs      <- read.table(file=paste(stat, "DEGs.txt", sep=""),                   sep="\t", header=F, row.names=NULL, check.names=F, stringsAsFactors=F)$V1
logFC     <- read.table(file=paste(processed_data, "logFC_matrix.txt", sep=""),   sep="\t", header=T, row.names=1,    check.names=F, stringsAsFactors=F) %>% as.matrix
contrasts <- read.table(file=paste(processed_data, "design.txt", sep=""),         sep="\t", header=T, row.names=NULL, check.names=F, stringsAsFactors=F)
tree      <-  read.tree(file=paste(processed_data, "dendrogram.tree", sep=""))

# list of GST genes, based on Pislewska-Bednarek et al., 2018, Plant Physiol.
gst       <- read.table(file=paste(ref_data, "TAIR10_GST.txt", sep=""),           sep="\t",  header=T, row.names=NULL, check.names=F, stringsAsFactors=F, comment.char="")
gst_tree  <-  read.tree(file=paste(ref_data, "RAxML_bipartitions.size_refined_GST-2.txt", sep=""))

# filter tree to Ath GSTs
idx <- str_detect(gst_tree$tip.label, "Athaliana")
gst_tree <- drop.tip(gst_tree, gst_tree$tip.label[!idx])

# rename tip label (from PacID to gene ID)
gst_tree$tip.label <- str_replace(gst_tree$tip.label, "Athaliana\\|", "")

idx <- match(gst_tree$tip.label, gst$PAC_ID)
gst_tree$tip.label <- gst$ID[idx]

# filter gst list to Ath GSTs, those on the tree
idx <- gst$ID %in% gst_tree$tip.label
gst <- gst[idx, ]

idx <- gst_tree$tip.label %in% gst$ID
gst_tree <- drop.tip(gst_tree, gst_tree$tip.label[!idx])

# replace tip label
gst$label <- paste(gst$ID, gst$Name, sep="_")
idx <- gst$ID == gst$Name
gst$label[idx] <- gst$ID[idx]

idx <- match(gst_tree$tip.label, gst$ID)
gst_tree$tip.label <- gst$label[idx]

# filter out mutants
idx <- contrasts$genotype %in% c("Col-0", "MYC2-FLAG")
contrasts <- contrasts[idx,]






# ==================== Prepare logFC melt table ==================== #
idx <- gst$ID %in% rownames(logFC)
gst <- gst[idx,]

idx <- match(gst$ID, rownames(logFC))
gst <- data.frame(gst, logFC[idx,], stringsAsFactors=F)

logFC_melt <- melt(gst, id.vars=c("Organism", "ID", "PAC_ID", "Name", "Description", "label"))
logFC_melt$fill <- saturate(logFC_melt$value, .01)

idx <- match(logFC_melt$variable, contrasts$name)
logFC_melt <- cbind(logFC_melt, contrasts[idx,])





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
	theme_tree(bgcolor=NA) +
	theme(legend.position="none",
		plot.margin=margin(t=5.5, b=5.5, r=5, l=0, unit="pt"))

# for y axis alignment by ylim2()
xy <- filter(p_tree, isTip)[, c("label", "x", "y")]

idx <- match(logFC_melt$variable, xy$label)
logFC_melt$y <- xy$y[idx]

x_max <- max(xy$x)
p_tree <- p_tree + xlim(NA, x_max + 1.8)






# GST phylogenetic plot
# tree plot #####
p_gst_tree <- ggtree(gst_tree, size=.8, branch.length="none") +
	geom_tiplab(size=3, hjust=1, angle=90, offset=-3, align=T, linetype="blank") +
	layout_dendrogram()	+
	theme_tree(bgcolor=NA) +
	theme(plot.margin=margin(t=5.5, b=50, r=5, l=0, unit="pt"))


# for y axis alignment by ylim2()
xy <- filter(p_gst_tree, isTip)[, c("label", "x", "y")]

idx <- match(logFC_melt$label, xy$label)
logFC_melt$x <- xy$y[idx]










# ==================== Plotting ==================== #
# heatmap
p_heat <- ggplot(logFC_melt, aes(y=y, x=x, fill=fill)) +
	geom_tile(alpha=1, size=0) +
	scale_fill_gradient2(low=c_cudo_magenta, mid=c_white, high=c_very_dark_green, midpoint=0, na.value=c_black, breaks=seq(ceiling(min(logFC_melt$value)), floor(max(logFC_melt$value)), by=1), guide="legend") +
	# scale_size_manual(values=c("0"=.5, "1"=1)) +
	theme_RTN +
	theme(legend.position="bottom",
		axis.text.x=element_blank(),
		axis.text.y=element_blank(),
		axis.line.x=element_blank(),
		axis.line.y=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y=element_blank(),
		panel.spacing=unit(.1, "lines"),
		panel.background=element_rect(fill="grey25"),
		plot.margin=margin(t=50, b=0, r=0, l=0, unit="pt"),
		panel.border=element_rect(fill=NA, size=.5, colour=c_black)) +
	ylim2(p_tree) +
	xlim2(p_gst_tree) +
	labs(x="", y="", fill="log2-FC") +
	guides(fill=guide_legend(override.aes=list(size=.5, colour=c_black)))


p <- plot_spacer() + p_gst_tree + p_tree + p_heat + plot_layout(ncol=2, nrow=2, widths=c(2,5), height=c(1, 8)) & theme(plot.background=element_blank())
ggsave(p, file=paste(fig, "GSTs-heatmap.png", sep=""), width=12, height=8, bg="transparent")


