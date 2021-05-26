
require(ggplot2)
source("/biodata/dep_psl/grp_psl/ThomasN/scripts/ggplot-themes_RTN.R")

# parameter
p.adj.method <- "fdr"
FC_threshold <- 1.5
alpha <- 0.05

# saturate function
saturate <- function(value, cut){
	if(cut >= .5){
		message("cut should be lower quantile, please")
		return(NULL)
	} else{
		max <- quantile(value, (1-cut))
		idx <- value > max
		value[idx] <- max

		min <- quantile(value, cut)
		idx <- value < min
		value[idx] <- min

		return(value)
	}
}



# colours #################################
lifestyles <- data.frame(
	names=c("0", "commensal", "beneficial", "pathogenic", "SynCom"),
	colours=c(c_black, c_dark_green, c_blue, c_cudo_magenta, c_dark_red),
	stringsAsFactors=F)

treatment_1 <- data.frame(
	names=c("R129_E", "R13_A", "Root491", "Root142", "Root1497", "Root700", "Bglumae", "Bphyt", "Ct", "Ci", "B_SynCom", "BFO_SynCom"),
	lifestyles=c("commensal", "commensal", "commensal", "commensal", "commensal", "commensal", "pathogenic", "beneficial", "beneficial", "pathogenic", "SynCom", "SynCom"),
	colours=c("#ff4500b3", "#ff4500b3", "#c80000b3", "#0000ffb3", "#9A0079b3", "#808000b3", "#32dc32b3", "#32dc32b3", "#53DEFFb3", "#53DEFFb3", "#654321b3", "#654321b3"),
	stringsAsFactors=F)

project <- data.frame(
	names=c("Alphaproteobacteria", "R129_E-flg22", "R129_E-time", "Kopriva", "Colletotrichum", "CJ", "ShijiHou", "ShijiHou_myc"),
	colours=c(c_very_dark_green, c_dark_green, c_green, c_dark_red, c_red, c_blue, c_cudo_skyblue, c_cudo_skyblue),
	stringsAsFactors=F)




