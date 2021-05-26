# partial R script for PCoA plotting
# a part of 2-compare_dataset.R


pcoa_1 <- cmdscale(as.dist(1-cor), k=3, eig=T)

point_1 <- as.data.frame(pcoa_1$points)
colnames(point_1) <- c("x", "y", "z")

idx <- match(rownames(point_1), design$name)
point_1 <- data.frame(point_1, design[idx,], row.names=NULL, stringsAsFactors=F)

point_1$lifestyle   <- factor(point_1$lifestyle,   levels=lifestyles$names)
point_1$project     <- factor(point_1$project,     levels=project$names)
point_1$treatment_1 <- factor(point_1$treatment_1, levels=treatment_1$names)

# by lifestyle
p1 <- ggplot(point_1, aes(x=x, y=y, colour=lifestyle, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=lifestyles$colours[-1]) +
	labs(x=paste("PC1: ", format(pcoa_1$eig[1]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 y=paste("PC2: ", format(pcoa_1$eig[2]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 title="PCoA (logFC)") +
	theme_RTN +
	theme(legend.position="none")

p2 <- ggplot(point_1, aes(x=x, y=z, colour=lifestyle, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=lifestyles$colours[-1]) +
	labs(x=paste("PC1: ", format(pcoa_1$eig[1]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 y=paste("PC3: ", format(pcoa_1$eig[3]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 title="PCoA (logFC)") +
	theme_RTN +
	theme(legend.position="right")

p3 <- p1 + p2 & theme(plot.background=element_blank())


# by project
p4 <- ggplot(point_1, aes(x=x, y=y, colour=project, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=project$colours) +
	labs(x=paste("PC1: ", format(pcoa_1$eig[1]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 y=paste("PC2: ", format(pcoa_1$eig[2]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 title="PCoA (logFC)") +
	theme_RTN +
	theme(legend.position="none")

p5 <- ggplot(point_1, aes(x=x, y=z, colour=project, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=project$colours) +
	labs(x=paste("PC1: ", format(pcoa_1$eig[1]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 y=paste("PC3: ", format(pcoa_1$eig[3]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 title="PCoA (logFC)") +
	theme_RTN +
	theme(legend.position="right")

p6 <- p4 + p5 & theme(plot.background=element_blank())



# by microbial inoculum
p7 <- ggplot(point_1, aes(x=x, y=y, colour=treatment_1, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=treatment_1$colours) +
	labs(x=paste("PC1: ", format(pcoa_1$eig[1]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 y=paste("PC2: ", format(pcoa_1$eig[2]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 title="PCoA (logFC)") +
	theme_RTN +
	theme(legend.position="none")

p8 <- ggplot(point_1, aes(x=x, y=z, colour=treatment_1, label=short)) +
	geom_point(size=2) +
	geom_text_repel(size=2, show.legend=F) +
	scale_colour_manual(values=treatment_1$colours) +
	labs(x=paste("PC1: ", format(pcoa_1$eig[1]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 y=paste("PC3: ", format(pcoa_1$eig[3]/sum(pcoa_1$eig)*100, digits=4), "%", sep=""),
		 title="PCoA (logFC)") +
	theme_RTN +
	theme(legend.position="right")

p9 <- p7 + p8 & theme(plot.background=element_blank())
