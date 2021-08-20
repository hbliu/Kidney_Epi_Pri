###########################################################################
#####    Effect size betwwen meta-anlayis eGFR GWAS and each GWAS study ###
#####    Only assocations with p < 0.00001 was extracted for comparion  ###
#####    Only assocations with p < 0.001 was used for PAGE and SUMMIT   ###
###########################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(Rmisc)
library(ggpubr)

CKDGen <- read.table("CKDGen.GWAS.formated.p00001.txt", header=T, sep="\t")
UKB <- read.table("UKB.GWAS.formated.p00001.txt", header=T, sep="\t")
MVP <- read.table("MVP.GWAS.formated.p00001.txt", header=T, sep="\t")
PAGE <- read.table("PAGE.GWAS.formated.p001.txt", header=T, sep="\t")
SUMMIT <- read.table("SUMMIT.GWAS.formated.p001.txt", header=T, sep="\t")
Metal <- read.table("Meta.eGFR.GWAS.CKDGen.UKB.MVP.PAGE.SUMMIT.p00001.txt", header=T, sep="\t")

colnames(CKDGen) <- c("RSID","CKDGen_REF","CKDGen_ALT","CKDGen_REF_Freq","CKDGen_BETA_Raw","CKDGen_PVAL","CKDGen_SE","CKDGen_N")
colnames(UKB) <- c("RSID","UKB_REF","UKB_ALT","UKB_REF_Freq","UKB_BETA_Raw","UKB_PVAL","UKB_SE","UKB_N")
colnames(MVP) <- c("RSID","MVP_REF","MVP_ALT","MVP_REF_Freq","MVP_BETA_Raw","MVP_SE","MVP_PVAL","MVP_N")
colnames(PAGE) <- c("RSID","PAGE_REF","PAGE_ALT","PAGE_REF_Freq","PAGE_BETA_Raw","PAGE_PVAL","PAGE_SE","PAGE_N")
colnames(SUMMIT) <- c("RSID","SUMMIT_REF","SUMMIT_ALT","SUMMIT_REF_Freq","SUMMIT_BETA_Raw","SUMMIT_PVAL","SUMMIT_SE","SUMMIT_N")
colnames(Metal) <- c("RSID","Metal_REF","Metal_ALT","Metal_ALT_Freq","Metal_BETA_Raw","Metal_Zscore","Metal_PVAL","Metal_N")

### Combine
Metal.CKDGen <- inner_join(Metal[,c(1,2,3,5,6)], CKDGen[,c(1,2,3,5)], by = c("RSID" = "RSID"))
Metal.UKB <- inner_join(Metal[,c(1,2,3,5,6)], UKB[,c(1,2,3,5)], by = c("RSID" = "RSID"))
Metal.MVP <- inner_join(Metal[,c(1,2,3,5,6)], MVP[,c(1,2,3,5)], by = c("RSID" = "RSID"))
Metal.PAGE <- inner_join(Metal[,c(1,2,3,5,6)], PAGE[,c(1,2,3,5)], by = c("RSID" = "RSID"))
Metal.SUMMIT <- inner_join(Metal[,c(1,2,3,5,6)], SUMMIT[,c(1,2,3,5)], by = c("RSID" = "RSID"))

### Change direction if reversed alleles
Metal.CKDGen.SameREF <- subset(Metal.CKDGen, as.character(Metal_REF) == as.character(CKDGen_REF) & as.character(Metal_ALT) == as.character(CKDGen_ALT)); Metal.CKDGen.SameREF$CKDGen_BETA <- Metal.CKDGen.SameREF$CKDGen_BETA_Raw * (-1)
Metal.CKDGen.Reverse <- subset(Metal.CKDGen, as.character(Metal_REF) == as.character(CKDGen_ALT) & as.character(Metal_ALT) == as.character(CKDGen_REF)); Metal.CKDGen.Reverse$CKDGen_BETA <- Metal.CKDGen.Reverse$CKDGen_BETA_Raw
Metal.CKDGen <- rbind(Metal.CKDGen.SameREF, Metal.CKDGen.Reverse)

Metal.UKB.SameREF <- subset(Metal.UKB, as.character(Metal_REF) == as.character(UKB_REF) & as.character(Metal_ALT) == as.character(UKB_ALT)); Metal.UKB.SameREF$UKB_BETA <- Metal.UKB.SameREF$UKB_BETA_Raw * (-1)
Metal.UKB.Reverse <- subset(Metal.UKB, as.character(Metal_REF) == as.character(UKB_ALT) & as.character(Metal_ALT) == as.character(UKB_REF)); Metal.UKB.Reverse$UKB_BETA <- Metal.UKB.Reverse$UKB_BETA_Raw
Metal.UKB <- rbind(Metal.UKB.SameREF, Metal.UKB.Reverse)

Metal.MVP.SameREF <- subset(Metal.MVP, as.character(Metal_REF) == as.character(MVP_REF) & as.character(Metal_ALT) == as.character(MVP_ALT)); Metal.MVP.SameREF$MVP_BETA <- Metal.MVP.SameREF$MVP_BETA_Raw * (-1)
Metal.MVP.Reverse <- subset(Metal.MVP, as.character(Metal_REF) == as.character(MVP_ALT) & as.character(Metal_ALT) == as.character(MVP_REF)); Metal.MVP.Reverse$MVP_BETA <- Metal.MVP.Reverse$MVP_BETA_Raw
Metal.MVP <- rbind(Metal.MVP.SameREF, Metal.MVP.Reverse)

Metal.PAGE.SameREF <- subset(Metal.PAGE, as.character(Metal_REF) == as.character(PAGE_REF) & as.character(Metal_ALT) == as.character(PAGE_ALT)); Metal.PAGE.SameREF$PAGE_BETA <- Metal.PAGE.SameREF$PAGE_BETA_Raw * (-1)
Metal.PAGE.Reverse <- subset(Metal.PAGE, as.character(Metal_REF) == as.character(PAGE_ALT) & as.character(Metal_ALT) == as.character(PAGE_REF)); Metal.PAGE.Reverse$PAGE_BETA <- Metal.PAGE.Reverse$PAGE_BETA_Raw
Metal.PAGE <- rbind(Metal.PAGE.SameREF, Metal.PAGE.Reverse)

Metal.SUMMIT.SameREF <- subset(Metal.SUMMIT, as.character(Metal_REF) == as.character(SUMMIT_REF) & as.character(Metal_ALT) == as.character(SUMMIT_ALT)); Metal.SUMMIT.SameREF$SUMMIT_BETA <- Metal.SUMMIT.SameREF$SUMMIT_BETA_Raw * (-1)
Metal.SUMMIT.Reverse <- subset(Metal.SUMMIT, as.character(Metal_REF) == as.character(SUMMIT_ALT) & as.character(Metal_ALT) == as.character(SUMMIT_REF)); Metal.SUMMIT.Reverse$SUMMIT_BETA <- Metal.SUMMIT.Reverse$SUMMIT_BETA_Raw
Metal.SUMMIT <- rbind(Metal.SUMMIT.SameREF, Metal.SUMMIT.Reverse)

###### PLot
pdf(file='Beta_values_between_studies.pdf',width=13.6,height=2.7)
theme_set(theme_bw())
# p1 <- ggplot(Metal.Metal.Reversed.Zscore, aes(x=Metal_BETA, y=Metal_BETA_Raw)) + xlab("Metal BETA") + ylab("Metal BETA after reversing") + geom_point(size = 0.1) + geom_smooth(method = "lm", se = FALSE) + stat_cor(method = "spearman")
p2 <- ggplot(Metal.CKDGen, aes(y=CKDGen_BETA, x=Metal_BETA_Raw)) + ylab("Effect size of GWAS (CKDGen)") + xlab("Effect size\nof Meta-analysis eGFR GWAS")+ geom_point(size = 0.05) + geom_smooth(method = "lm", se = FALSE) + stat_cor(method = "spearman", label.y = 0.015*0.9) + xlim(-0.1,0.1) + ylim(-0.015,0.015) + theme(legend.position="none") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none") + geom_hline(yintercept=0, linetype="dashed", color = "gray") + geom_vline(xintercept=0, linetype="dashed", color = "gray")
p3 <- ggplot(Metal.UKB, aes(y=UKB_BETA, x=Metal_BETA_Raw)) + ylab("Effect size of GWAS (UKB)") + xlab("Effect size\nof Meta-analysis eGFR GWAS")+ geom_point(size = 0.05) + geom_smooth(method = "lm", se = FALSE) + stat_cor(method = "spearman", label.y = 0.12*0.9) + xlim(-0.1,0.1) + ylim(-0.12,0.12) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none") + geom_hline(yintercept=0, linetype="dashed", color = "gray") + geom_vline(xintercept=0, linetype="dashed", color = "gray")
p6 <- ggplot(Metal.MVP, aes(y=MVP_BETA, x=Metal_BETA_Raw)) + ylab("Effect size of GWAS (MVP)") + xlab("Effect size\nof Meta-analysis eGFR GWAS")+ geom_point(size = 0.05) + geom_smooth(method = "lm", se = FALSE) + stat_cor(method = "spearman", label.y = 2*0.9) + xlim(-0.1,0.1) + ylim(-2,2) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none") + geom_hline(yintercept=0, linetype="dashed", color = "gray") + geom_vline(xintercept=0, linetype="dashed", color = "gray") 
p4 <- ggplot(Metal.PAGE, aes(y=PAGE_BETA, x=Metal_BETA_Raw)) + ylab("Effect size of GWAS (PAGE)") + xlab("Effect size\nof Meta-analysis eGFR GWAS")+ geom_point(size = 0.05) + geom_smooth(method = "lm", se = FALSE) + stat_cor(method = "spearman", label.y = 2.2*0.9) + xlim(-0.05,0.05) + ylim(-2.2,2.2) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none") + geom_hline(yintercept=0, linetype="dashed", color = "gray") + geom_vline(xintercept=0, linetype="dashed", color = "gray") 
p5 <- ggplot(Metal.SUMMIT, aes(y=SUMMIT_BETA, x=Metal_BETA_Raw)) + ylab("Effect size of GWAS (SUMMIT)") + xlab("Effect size\nof Meta-analysis eGFR GWAS")+ geom_point(size = 0.05) + geom_smooth(method = "lm", se = FALSE) + stat_cor(method = "spearman", label.y = 5*0.9) + xlim(-0.05,0.05) + ylim(-5,5) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none") + geom_hline(yintercept=0, linetype="dashed", color = "gray") + geom_vline(xintercept=0, linetype="dashed", color = "gray")
#multiplot(p1, p3, p5, p7, p2, p4, p6, p8, layout = matrix(c(1,2,3,4,5,6,7,8), nrow= 2, byrow=TRUE))
multiplot(p2, p3, p6, p4, p5, layout = matrix(c(1,2,3,4,5), nrow= 1, byrow=TRUE))
dev.off()

