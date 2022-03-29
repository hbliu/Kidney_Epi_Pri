###########################################################################
#####    Code for Figure.7                                            #####
#####    Hongbo Liu (hongbo919@gmail.com)                             #####
###########################################################################


###########################################################################
#####    Code for Figure.7a                                           #####
###########################################################################
##### http://locuszoom.org/genform.php?type=yourdata


###########################################################################
#####    Code for Figure.7b                                           #####
###########################################################################
##### Note: Genotype (rs2252281, x-axis) and normalized CpG methylation (cg15971010, y-axis) association in human kidneys (n=443)
library(ggplot2)
library(ggsci)
library(ggpubr)
S443.Methyl.Genotype <- read.table("Chr17.rs2252281.cg15971010.S443.Methyl.Genotype.infor.txt", header=T, sep="\t")
meQTL <- read.table("Chr17.rs2252281.cg15971010.S443.meQTL.infor.txt", header=T, sep="\t")
S443.Methyl.Genotype <- na.omit(S443.Methyl.Genotype)

pdf(file='Figure.7b.pdf',width=3,height=3.8)
ggboxplot(S443.Methyl.Genotype, x = "rs2252281", y = "cg15971010", fill = "rs2252281", title = paste ("meQTL\nbeta = ", round(meQTL$beta,2), "\np = ", signif(meQTL$pvalue,2), sep = ""), add = "jitter", add.params = list(size = 0.1, jitter = 0.1, color = "black"), legend = "none") + scale_color_npg() + scale_fill_npg() + theme(plot.title = element_text(size=12))
dev.off()


###########################################################################
#####    Code for Figure.7c                                           #####
###########################################################################
##### Note: Genotype (rs2252281 x-axis) and normalized gene expression (SLC47A1, y-axis) association in human kidney tubule samples (n=356).
library(ggplot2)
library(ggsci)
library(ggpubr)
S356.Gex.Genotype <- read.table("Gex.Genotype.infor.txt", header=T, sep="\t")
eQTL <- read.table("eQTL.infor.txt", header=T, sep="\t")
S356.Gex.Genotype <- na.omit(S356.Gex.Genotype)

pdf(file='Figure.7b.pdf',width=3,height=3.8)
ggboxplot(S356.Gex.Genotype, x = "rs2252281", y = "SLC47A1", fill = "rs2252281", title = paste ("eQTL\nbeta = ", round(eQTL$beta,2), "\np = ", signif(eQTL$pvalue,2), sep = ""), add = "jitter", add.params = list(size = 0.1, jitter = 0.1, color = "black"), legend = "none") + scale_color_npg() + scale_fill_npg() + theme(plot.title = element_text(size=12))
dev.off()


###########################################################################
#####    Code for Figure.7d                                           #####
###########################################################################
##### Note: Manhattan plot of phenome-wide association study of predicted loss-of-function (pLOF) variants in SLC47A1 in UKBB.
library(dplyr)
library(PheWAS)
load("SLC47A1_pLOF_burden_PheWAS.RData")
pdf(file='Figure.7d.pdf',width=7,height=3.2)
phewasManhattan(SLC47A1_pLOF_burden_PheWAS, Phenotypes, OR.direction = , annotate.angle= 35, y.axis.interval = 1, sort.by.category.value=T, annotate.list = annoteLable.pLOF, base.labels =T, annotate.size = 0, title="pLOF burden PheWAS manhattan plot")
dev.off()








