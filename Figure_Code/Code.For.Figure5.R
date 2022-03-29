###########################################################################
#####    Code for Figure.5                                            #####
#####    Hongbo Liu (hongbo919@gmail.com)                             #####
###########################################################################


###########################################################################
#####    Code for Figure.5b                                           #####
###########################################################################
##### Note: Enrichment of kidney methylation-mediated heritability for kidney function GWAS traits (x-axis) in kidney cell type-specific accessible regions (y-axis).
##### White to red indicates h2med enrichment (nominal p < 0.05 calculated by MESC).
##### Asterisk indicates h2med enrichment passing FDR q < 0.05 (accounting for 408 tests for 12 cell type CpG sets and 34 GWAS traits).
library(ggplot2)
library(gplots)
h2med_enrichment <- read.table("Kidney.meQTL.Heritability.snATAC_Cluster.h2med_enrichment.txt",header=T, sep="\t")
Kidney.CpGset.h2med.snATAC_Cluster <- read.table("Kidney.CpGset.h2med.snATAC_Cluster.txt",header=T, sep="\t")

# A 4 x 4 matrix of characters (stars)
labs <- h2med_enrichment
for (i in 1:nrow(labs)) {
	for (j in 1:ncol(labs)) {
		labs[i,j] <- ""
		if (subset(Kidney.CpGset.h2med.snATAC_Cluster, Gene_category == rownames(labs)[i] & GWAS == colnames(labs)[j])$h2med_enrichment_qvalue < 0.05 & h2med_enrichment[i,j] > 0){
		labs[i,j] <- "*"
		}
	}
}

my_palette <- colorRampPalette(c("gray", "white", "red"))(n = 74)
pdf("Figure.5e.pdf", width = 7, height = 8)
x<- as.matrix(h2med_enrichment); 
heatmap <- heatmap.2(x, scale = "none", col=my_palette, trace = "none", density.info = "none", Rowv=FALSE, Colv=FALSE, cellnote=labs, notecol="black", notecex=2, margins = c(20, 20))
dev.off()



###########################################################################
#####    Code for Figure.5c                                           #####
###########################################################################
##### Note: Single cell GWAS trait enrichment in human kidney cells by gchromVAR.
library(ggplot2)
library(gplots)
UKBB2KidATAC.snATAC_Cluster.Mean <- read.table("UKBB2KidATAC.snATAC_Cluster.Mean.tsv", header=T, sep=",")

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 256)
library("gplots")
pdf("Figure.5c.pdf", width = 7, height = 7)
x<- as.matrix(UKBB2KidATAC.snATAC_Cluster.Mean); 
heatmap <- heatmap.2(x, scale = "none", Rowv = FALSE, Colv = FALSE, col=my_palette, trace = "none", density.info = "none", margins = c(7, 7))
dev.off()








