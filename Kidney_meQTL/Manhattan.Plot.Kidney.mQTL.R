###########################################################################
#####    Manhattan plot of kidney mQTL                                #####
#####    Sample size: 443 human kidney samples                        #####
#####    Due to limitation of R, only mQTL with p < 0.01 was used     #####
###########################################################################

library(dplyr)
library(qqman)
cis_mQTL_forManhattan <- read.table("cis_mQTL_forManhattan.txt", header=F, sep="\t")
colnames(cis_mQTL_forManhattan) <- c("CHR","SNP_loci","snps","pvalue")
jpeg('Manhattan of p values of kidney cis_mQTL.jpg', width = 20, height = 6, units = 'in', res = 300)
manhattan(cis_mQTL_forManhattan, chr="CHR", bp="SNP_loci", snp="snps", p="pvalue", col = c("deepskyblue4", "goldenrod3"), suggestiveline = F, genomewideline = F, cex = 0.05, main = "cis_mQTL_forManhattan")
dev.off()

