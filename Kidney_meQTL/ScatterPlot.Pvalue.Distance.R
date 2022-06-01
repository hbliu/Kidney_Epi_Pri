###########################################################################
#####    ScatterPlot of mQTL Pvalue and Distance                      #####
#####    Sample size: 443 human kidney samples                        #####
#####    For each mCpG, the top variant was used                      #####
###########################################################################
library(dplyr)
library(Rmisc)
library(ggplot2)
library(viridis)
library(ggpointdensity)

Best_mQTL_SNP.Sig <- read.table("cis_mQTL_Sig.mCpG.topSNPs.txt", header=F, sep="\t")

#### Plot
p1 <- ggplot(Best_mQTL_SNP.Sig, aes(x=Dis_SNP2CpG/1000000, y= -log10(pvalue))) + geom_point(size = 0.3) + xlab("Distance from SNP to CpG (Mb)") + ylab("-log10(p value)") + theme_bw() + theme_classic(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2 <- ggplot(Best_mQTL_SNP.Sig, aes(x=Dis_SNP2TSS/1000000, y= -log10(pvalue))) + geom_point(size = 0.3) + xlab("Distance from SNP to TSS (Mb)") + ylab("-log10(p value)") + theme_bw() + theme_classic(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p3 <- ggplot(Best_mQTL_SNP.Sig, aes(x=Dis_CpG2TSS/1000000, y= -log10(pvalue))) + geom_point(size = 0.3) + xlab("Distance from CpG to TSS (Mb)") + ylab("-log10(p value)") + theme_bw() + theme_classic(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#### Add density
p1 <- p1 + geom_rug(alpha = 1/5, size=0.01, colour = "brown")
p2 <- p2 + geom_rug(alpha = 1/5, size=0.01, colour = "brown")
p3 <- p3 + geom_rug(alpha = 1/5, size=0.01, colour = "brown")

jpeg('ScatterPlot of mQTL Pvalue and Distance.jpg', width = 12, height = 4, units = 'in', res = 300)
theme_set(theme_bw())
multiplot(p1, p2, p3, layout = matrix(c(1,2,3), nrow=1, byrow=TRUE))
dev.off()

