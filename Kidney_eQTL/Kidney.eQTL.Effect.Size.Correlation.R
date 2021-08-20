###########################################################################
#####    Effect size betwwen meta-anlayis kidney eQTL and each eQTL study #
#####    Only assocations with p < 0.00001 was extracted for comparion  ###
###########################################################################
R
library(dplyr)
library(tidyr)
library(ggplot2)
library(Rmisc)
library(ggpubr)

Meta.eQTL <- read.table("Kidney.eQTL.Meta.Analsyis.MAF0.05.txt", header=F, sep="\t")
colnames(Meta.eQTL) <- c("RSID_GeneSymbol","Meta_REF","Meta_ALT","Meta_Beta","Meta_PVAL")
Meta.eQTL.p00001 <- subset(Meta.eQTL, Meta_PVAL < 1e-5)

Sheng.eQTL <- read.table("Sheng.eQTL.txt", header=F, sep="\t")
colnames(Sheng.eQTL) <- c("RSID_GeneSymbol","Sheng_REF","Sheng_ALT","Sheng_BETA_Raw","Sheng_SE_Raw","Sheng_PVAL","Sheng_N")
Sheng.eQTL.p00001 <- subset(Sheng.eQTL, Sheng_PVAL < 1e-5)

Ko.eQTL <- read.table("Ko.eQTL.txt", header=F, sep="\t")
colnames(Ko.eQTL) <- c("RSID_GeneSymbol","Ko_REF","Ko_ALT","Ko_BETA_Raw","Ko_SE","Ko_PVAL","Ko_N")
Ko.eQTL.p00001 <- subset(Ko.eQTL, Ko_PVAL < 1e-5)

GTEx.eQTL <- read.table("GTEx.eQTL.txt", header=F, sep="\t")
colnames(GTEx.eQTL) <- c("RSID_GeneSymbol","GTEx_REF","GTEx_ALT","GTEx_BETA_Raw","GTEx_SE","GTEx_PVAL","GTEx_N")
GTEx.eQTL.p00001 <- subset(GTEx.eQTL, GTEx_PVAL < 1e-5)

NephQTL.eQTL <- read.table("NephQTL.eQTL.txt", header=F, sep="\t")
colnames(NephQTL.eQTL) <- c("RSID_GeneSymbol","NephQTL_REF","NephQTL_ALT","NephQTL_BETA_Raw","NephQTL_SE","NephQTL_PVAL","NephQTL_N")
NephQTL.eQTL.p00001 <- subset(NephQTL.eQTL, NephQTL_PVAL < 1e-5)

Sheng.Meta <- inner_join(Sheng.eQTL.p00001, Meta.eQTL.p00001, by = c("RSID_GeneSymbol" = "RSID_GeneSymbol"))
Ko.Meta <- inner_join(Ko.eQTL.p00001, Meta.eQTL.p00001, by = c("RSID_GeneSymbol" = "RSID_GeneSymbol"))
GTEx.Meta <- inner_join(GTEx.eQTL.p00001, Meta.eQTL.p00001, by = c("RSID_GeneSymbol" = "RSID_GeneSymbol"))
NephQTL.Meta <- inner_join(NephQTL.eQTL.p00001, Meta.eQTL.p00001, by = c("RSID_GeneSymbol" = "RSID_GeneSymbol"))

Sheng.Meta <- subset(Sheng.Meta, as.character(Sheng_REF) == as.character(Meta_REF) & as.character(Sheng_ALT) == as.character(Meta_ALT))
Ko.Meta <- subset(Ko.Meta, as.character(Ko_REF) == as.character(Meta_REF) & as.character(Ko_ALT) == as.character(Meta_ALT))
GTEx.Meta <- subset(GTEx.Meta, as.character(GTEx_REF) == as.character(Meta_REF) & as.character(GTEx_ALT) == as.character(Meta_ALT))
NephQTL.Meta <- subset(NephQTL.Meta, as.character(NephQTL_REF) == as.character(Meta_REF) & as.character(NephQTL_ALT) == as.character(Meta_ALT))

###### PLot
pdf(file='BETA_Distribution_between_Metal_BETA_Beta.AfterReverse.Used.v5.pdf',width=12,height=3)
theme_set(theme_bw())
p2 <- ggplot(Sheng.Meta, aes(x=Sheng_BETA_Raw, y=Meta_Beta)) + xlab("Effect size of eQTL (Sheng et al.)") + ylab("Effect size of Meta-analysis eQTL")+ geom_point(size = 0.001) + stat_cor(method = "spearman", r.digits = 3, label.y = 2) + xlim(-2,2) + ylim(-2,2) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + gradient_fill(c("yellow", "red")) + theme(legend.position="none") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none") 
p3 <- ggplot(Ko.Meta, aes(x=Ko_BETA_Raw, y=Meta_Beta)) + xlab("Effect size of eQTL (Ko et al.)") + ylab("Effect size of Meta-analysis eQTL")+ geom_point(size = 0.001) + stat_cor(method = "spearman", r.digits = 3, label.y = 2) + xlim(-2,2) + ylim(-2,2) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + gradient_fill(c("yellow", "red")) + theme(legend.position="none") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none") 
p4 <- ggplot(GTEx.Meta, aes(x=GTEx_BETA_Raw, y=Meta_Beta)) + xlab("Effect size of eQTL (GTEx v8)") + ylab("Effect size of Meta-analysis eQTL")+ geom_point(size = 0.001) + stat_cor(method = "spearman", r.digits = 3, label.y = 2) + xlim(-2,2) + ylim(-2,2) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + gradient_fill(c("yellow", "red")) + theme(legend.position="none") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none") 
p5 <- ggplot(NephQTL.Meta, aes(x=NephQTL_BETA_Raw, y=Meta_Beta)) + xlab("Effect size of eQTL (NephQTL)") + ylab("Effect size of Meta-analysis eQTL")+ geom_point(size = 0.001) + stat_cor(method = "spearman", r.digits = 3, label.y = 2) + xlim(-2,2) + ylim(-2,2) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + gradient_fill(c("yellow", "red")) + theme(legend.position="none") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none") 
p2 <- p2 + geom_hline(yintercept=0, linetype="dashed", color = "gray") + geom_vline(Shengtercept=0, linetype="dashed", color = "gray")
p3 <- p3 + geom_hline(yintercept=0, linetype="dashed", color = "gray") + geom_vline(Shengtercept=0, linetype="dashed", color = "gray")
p4 <- p4 + geom_hline(yintercept=0, linetype="dashed", color = "gray") + geom_vline(Shengtercept=0, linetype="dashed", color = "gray")
p5 <- p5 + geom_hline(yintercept=0, linetype="dashed", color = "gray") + geom_vline(Shengtercept=0, linetype="dashed", color = "gray")
multiplot(p2, p3, p4, p5, layout = matrix(c(1,2,3,4), nrow= 1, byrow=TRUE))
dev.off()

