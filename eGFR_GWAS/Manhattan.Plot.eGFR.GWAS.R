###########################################################################
#####    Manhattan plot of Kidney eQTL                                #####
#####    Software: ggman                                              #####
###########################################################################

###########################################################################
##### Manhattan plot of eGFR GWAS (meta-analysis with highlighted novel) ##
###########################################################################
library(ggman)
library(gg.gap)
GWAS <- read.table("Meta.eGFR.GWAS.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.forManhattanPlot.txt", header=T, sep="\t")
colnames(GWAS) <- c("CHR","POS","SNP","PVAL")
GWAS.Lead <- read.table("GWAS.LeadSNP.Infor.Status.txt", header=T, sep="\t")
GWAS.Lead$Group <- "Known"
GWAS.Lead[which(GWAS.Lead$Status == "Novel"),"Group"] <- "Novel"
GWAS.Lead.Novel <- subset(GWAS.Lead, Status == "Novel")

GWAS.plot <- GWAS

### Set 1e-150 as the smallest p value
GWAS.plot[which(GWAS.plot$PVAL < 1e-150),4] = 1e-150

### GWAS.plot
p1 <- ggman(GWAS.plot, sigLine = 7.30103, snp = "SNP", bp = "POS", chrom = "CHR", pvalue = "PVAL", relative.positions = TRUE) + scale_color_manual(values=c("deepskyblue4", "goldenrod3")) + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2 <- ggmanHighlight(p1, highlight = as.character(GWAS.Lead$MarkerName), colour = "darkcyan", size = 2, stroke = 0)
p3 <- ggmanHighlight(p2, highlight = as.character(GWAS.Lead.Novel$MarkerName), colour = "red", size = 2, stroke = 0)

### Add closest gene of lead SNP (Novel loci only)
jpeg('Manhattan plot of p values of meta-analysis eGFR GWAS.jpg', width = 12, height = 4.5, units = 'in', res = 500, quality = 100)
ggmanLabel(p3, labelDfm = GWAS.Lead.Novel, snp = "MarkerName", label = "Cloest.Gene", type = "text", colour = "black", fontface = 'italic', size = 3.5,
	nudge_y = 75,
	max.iter = 5000,
    vjust=1,
    direction='y',
    nudge_x=0.1,
    segment.color = "grey50",
    segment.alpha = 0.5,
    segment.size  = 0.2)
dev.off()


###########################################################################
##### Manhattan plot of eGFR GWAS of each datasets separately         #####
###########################################################################
library(dplyr)
library(tidyr)
library(ggman)

#### Read data for Manhattan plot
Meta <- read.table("Meta.eGFR.GWAS.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.forManhattanPlot.txt", header=T, sep="\t")
CKDGen <- read.table("CKDGen.GWAS.txt", header=T, sep="\t")
UKB <- read.table("UKB.GWAS.txt", header=T, sep="\t")
MVP <- read.table("MVP.GWAS.txt", header=T, sep="\t")
PAGE <- read.table("PAGE.GWAS.txt", header=T, sep="\t")
SUMMIT <- read.table("SUMMIT.GWAS.txt", header=T, sep="\t")

#### Manhattan plot
jpeg('manhattan of p values of Meta-analysis eGFR GWAS.jpg', width = 20, height = 4, units = 'in', res = 500, quality = 100)
ggman(Meta, sigLine = 7.30103, snp = "SNP", bp = "POS", chrom = "CHR", pvalue = "PVAL", relative.positions = TRUE) + scale_color_manual(values=c("deepskyblue4", "goldenrod3")) + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

jpeg('manhattan of p values of CKDGen eGFR GWAS.jpg', width = 20, height = 4, units = 'in', res = 500, quality = 100)
ggman(CKDGen, sigLine = 7.30103, snp = "SNP", bp = "POS", chrom = "CHR", pvalue = "PVAL", relative.positions = TRUE) + scale_color_manual(values=c("deepskyblue4", "goldenrod3")) + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

jpeg('manhattan of p values of UKB eGFR GWAS.jpg', width = 20, height = 4, units = 'in', res = 500, quality = 100)
ggman(UKB, sigLine = 7.30103, snp = "SNP", bp = "POS", chrom = "CHR", pvalue = "PVAL", relative.positions = TRUE) + scale_color_manual(values=c("deepskyblue4", "goldenrod3")) + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

jpeg('manhattan of p values of MVP eGFR GWAS.jpg', width = 20, height = 4, units = 'in', res = 500, quality = 100)
ggman(MVP, sigLine = 7.30103, snp = "SNP", bp = "POS", chrom = "CHROM", pvalue = "PVAL", relative.positions = TRUE) + scale_color_manual(values=c("deepskyblue4", "goldenrod3")) + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

jpeg('manhattan of p values of PAGE eGFR GWAS.jpg', width = 20, height = 4, units = 'in', res = 500, quality = 100)
ggman(PAGE, sigLine = 7.30103, snp = "SNP", bp = "POS", chrom = "CHR", pvalue = "PVAL", relative.positions = TRUE) + scale_color_manual(values=c("deepskyblue4", "goldenrod3")) + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

jpeg('manhattan of p values of SUMMIT eGFR GWAS.jpg', width = 20, height = 4, units = 'in', res = 500, quality = 100)
ggman(SUMMIT, sigLine = 7.30103, snp = "SNP", bp = "POS", chrom = "CHR", pvalue = "PVAL", relative.positions = TRUE) + scale_color_manual(values=c("deepskyblue4", "goldenrod3")) + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


