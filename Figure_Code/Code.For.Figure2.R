###########################################################################
#####    Code for Figure.2                                            #####
#####    Hongbo Liu (hongbo919@gmail.com)                             #####
###########################################################################


###########################################################################
#####    Code for Figure.2a                                           #####
###########################################################################
##### Note: Manhattan plot of eGFRcrea GWAS of 1,508,659 individuals. X-axis is chromosomal location of SNP
library(ggman)
GWAS <- read.table(gzfile("GWAS.forManhattanPlot.txt.gz"), header=T, sep="\t")
colnames(GWAS) <- c("CHR","POS","SNP","PVAL")
GWAS[which(GWAS$PVAL < 1e-150),4] = 1e-150
IndependentLoci <- read.table(gzfile("Independent.Loci.txt.gz"), header=T, sep="\t")

p1 <- ggman(GWAS, sigLine = 7.30103, snp = "SNP", bp = "POS", chrom = "CHR", pvalue = "PVAL", relative.positions = TRUE) + scale_color_manual(values=c("deepskyblue4", "goldenrod3")) + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2 <- ggmanHighlight(p1, highlight = as.character(GWAS.Lead$MarkerName), colour = "darkcyan", size = 2, stroke = 0)
p3 <- ggmanHighlight(p2, highlight = as.character(GWAS.Lead.Close2Known$MarkerName), colour = "darkcyan", size = 2, stroke = 0)
p3 <- ggmanHighlight(p3, highlight = as.character(GWAS.Lead.Novel$MarkerName), colour = "red", size = 3, stroke = 0)

GWAS.Lead.Novel <- subset(GWAS.Lead, Group == "2_Novel")
GWAS.Lead.Novel.TopPrioritizedGene <- subset(GWAS.Lead.Novel, TopPrioritizedGene != "")
GWAS.Lead.Novel.NoTopPrioritizedGene <- subset(GWAS.Lead.Novel, TopPrioritizedGene == "")
GWAS.Lead.Novel.NoTopPrioritizedGene$TopPrioritizedGene <- GWAS.Lead.Novel.NoTopPrioritizedGene$Cloest.Gene
GWAS.Lead.Novel <- rbind(GWAS.Lead.Novel.TopPrioritizedGene, GWAS.Lead.Novel.NoTopPrioritizedGene)

p4 <- ggmanLabel(p3, labelDfm = GWAS.Lead.Novel, snp = "MarkerName", label = "TopPrioritizedGene", type = "text", colour = "white", fontface = 'italic', size = 3.5,
	nudge_y = 75,
	max.iter = 5000,
    vjust=1,
    direction='y',
    nudge_x=0.1,
    segment.color = "grey50",
    segment.alpha = 0.5,
    segment.size  = 0.2)

jpeg('Figure.2a.jpg', width = 12, height = 4.5, units = 'in', res = 500, quality = 100)
p4
dev.off()



###########################################################################
#####    Code for Figure.2b                                           #####
###########################################################################
##### Note: Scatter plot of minor allele frequency and effect size (absolute) of lead variants from 878 eGFRcrea GWAS loci.
##### Cyan represents loci overlapping or being physically nearby (within 500kb or LD R2 > 0.001) previously reported sentinel variants, and red novel loci. 
library(ggplot2)
IndependentLoci <- read.table("Independent.Loci.txt", header=T, sep="\t")
dim(IndependentLoci)
IndependentLoci$absBETA <- abs(IndependentLoci$BETA)

pdf(file='Figure.2b.pdf', width=7,height=6)
ggplot(IndependentLoci, aes(x = MAF, y = absBETA)) +
  geom_point(aes(color = Status),size = 1.5) + 
  scale_y_continuous(breaks=seq(0, 0.16, by = 0.02), limits = c(0, 0.16)) +
  scale_color_manual(values = c("darkcyan", "red")) +
  theme_bw(base_size = 16) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = c(0.9, 0.90))
dev.off()



###########################################################################
#####    Code for Figure.2c                                           #####
###########################################################################
##### Note: Human kidney eQTL Manhattan plot (n= 686 kidney samples)
library(qqman)

eQTL <- read.table("eQTL.forManhattanPlot.txt", header=T, sep="\t"); colnames(eQTL) <- c("CHR","POS","SNP","PVAL")
eQTL[which(eQTL$PVAL < 4.27e-308),4] = 4.27e-308

library(qqman)
jpeg('Figure.2c.jpg', width = 20, height = 5, units = 'in', res = 500, quality = 100)
manhattan(eQTL, chr="CHR", bp="POS", snp="SNP", p="PVAL", col = c("deepskyblue4", "goldenrod3"), suggestiveline = F, genomewideline = F, cex = 0.1, main = "Manhattan Plot")
dev.off()


###########################################################################
#####    Code for Figure.2d                                           #####
###########################################################################
##### Note: Enrichment of kidney specific eQTL SNPs to GWAS traits
library(ggplot2)
library(ggrepel)

Kidney.Specific.eQTL.2GWAS <- read.table("Kidney.Specific.eQTL.2GWAS.txt", header=T, sep="\t")
Kidney.Specific.eQTL.2GWAS$Minus_log10P <- -log10(Kidney.Specific.eQTL.2GWAS$chisq.test.p)

p1 <- ggplot(Kidney.Specific.eQTL.2GWAS, aes(x=fisher.test.oddsRatio, y=Minus_log10P, label = GWAS_Traits)) + geom_point() + ggtitle("Enrichment of Kidney Specific eQTL SNPs in GWAS traits") + theme(legend.position="none") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2 <- p1 + xlab("Odds Ratio") + ylab("-log10(p)")
p3 <- p2 + geom_point(aes(size = SNPNum, color = Trait_category)) + scale_size_continuous(name = "SNP Num", range = c(1, 5)) + geom_hline(yintercept= -log10(0.05), linetype="dashed", color = "gray") + geom_vline(xintercept= 1, linetype="dashed", color = "gray")

pdf(file='Figure.2d.pdf',width=4.6,height=3.7)
p3 + geom_text_repel(data = subset(Kidney.Specific.eQTL.2GWAS,Minus_log10P > 100 & fisher.test.oddsRatio>1.5 & Trait_category == "Renal"), size = 4) + ylim(0,305)
dev.off()

