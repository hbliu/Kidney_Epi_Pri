###########################################################################
#####    Manhattan plot of Kidney eQTL                                #####
#####    Software: qqman                                              #####
###########################################################################
library(dplyr)
library(tidyr)
library(qqman)

######### 
Meta.eQTL <- read.table("Meta.eQTL.Flip_Reversed.forManhattanPlot.txt", header=T, sep="\t"); colnames(Meta.eQTL) <- c("CHR","POS","SNP","PVAL")
jpeg('manhattan of p values of Meta.eQTL.jpg', width = 20, height = 5, units = 'in', res = 500, quality = 100)
manhattan(Meta.eQTL, chr="CHR", bp="POS", snp="SNP", p="PVAL", col = c("deepskyblue4", "goldenrod3"), suggestiveline = F, genomewideline = F, cex = 0.1, main = "Manhattan")
dev.off()
