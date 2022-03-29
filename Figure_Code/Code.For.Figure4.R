###########################################################################
#####    Code for Figure.4                                            #####
#####    Hongbo Liu (hongbo919@gmail.com)                             #####
###########################################################################


###########################################################################
#####    Code for Figure.4b                                           #####
###########################################################################
##### Note: Estimated proportion of heritability (h2med/h2g) mediated by DNA methylation and gene expression in 414 kidneys across 34 GWAS traits.
library(dplyr)
library(ggplot2)
library(Rmisc)
library(ggrepel)

read.table("h2med.Estimate_over_h2.Combined.meQTL.eQTL.txt", header=T, sep="\t")-> h2med_over_h2
read.table("h2med.Estimate_over_h2_SE.Combined.meQTL.eQTL.txt", header=T, row.names = 1, sep="\t")-> h2med_over_h2_SE

pdf(file='Figure.4b.pdf', height = 4, width = 4)
ggplot(h2med_over_h2, aes(y=meQTL.h2med_over_h2, x=eQTL.h2med_over_h2)) + theme(legend.position="none") +
geom_pointrange(aes(xmax = eQTL.h2med_over_h2 + eQTL.h2med_over_h2.SE, xmin = eQTL.h2med_over_h2 - eQTL.h2med_over_h2.SE, color = TraitType)) +
geom_pointrange(aes(ymax = meQTL.h2med_over_h2 + meQTL.h2med_over_h2.SE, ymin = meQTL.h2med_over_h2 - meQTL.h2med_over_h2.SE, color = TraitType)) +
geom_point(aes(color = TraitType, size = 3)) + scale_color_manual(values=c("gray57", "darkgreen", "gold3", "dodgerblue", "magenta")) + geom_abline(intercept = 0, slope = 1, linetype="dashed", color = "gray") +
scale_x_continuous(limits = c(-0.08,0.7), breaks=seq(0, 0.6, by = 0.2)) +
scale_y_continuous(limits = c(-0.08,0.7),breaks=seq(0, 0.6, by = 0.2)) +
theme_bw(base_size = 12) +
theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = c(2, 2))
dev.off()


###########################################################################
#####    Code for Figure.4e                                           #####
###########################################################################
##### Note: Enrichment of kidney methylation-mediated heritability for kidney function GWAS traits (x-axis) in different ChromHMM regulatory elements (y-axis).
##### White to red indicates h2med enrichment (nominal two-sided p < 0.05 calculated by MESC).
##### Asterisk indicates h2med enrichment passing FDR q < 0.05 (accounting for 374 tests for 11 chromatin state CpG sets and 34 GWAS traits).
library(ggplot2)
library(gplots)
h2med_enrichment <- read.table("Kidney.meQTL.Heritability.HMM.h2med_enrichment.txt",header=T, sep="\t")
Kidney.CpGset.h2med.HMM_State <- read.table("Kidney.CpGset.h2med.HMM_State.txt",header=T, sep="\t")

# A 4 x 4 matrix of characters (stars)
labs <- h2med_enrichment
for (i in 1:nrow(labs)) {
	for (j in 1:ncol(labs)) {
		labs[i,j] <- ""
		if (subset(Kidney.CpGset.h2med.HMM_State, Gene_category == rownames(labs)[i] & GWAS == colnames(labs)[j])$h2med_enrichment_qvalue < 0.05 & h2med_enrichment[i,j] > 0){
		labs[i,j] <- "*"
		}
	}
}

my_palette <- colorRampPalette(c("gray", "white", "red"))(n = 74)
pdf("Figure.4e.pdf", width = 7, height = 8)
x<- as.matrix(h2med_enrichment); 
heatmap <- heatmap.2(x, scale = "none", col=my_palette, trace = "none", density.info = "none", Rowv=FALSE, Colv=FALSE, cellnote=labs, notecol="black", notecex=2, margins = c(20, 20))
dev.off()





