###########################################################################
#####    Code for Figure.3                                            #####
#####    Hongbo Liu (hongbo919@gmail.com)                             #####
###########################################################################


###########################################################################
#####    Code for Figure.3b                                           #####
###########################################################################
##### Note: Manhattan plot of human kidney meQTL data (n=443).
library(qqman)

meQTL <- read.table("meQTL.forManhattanPlot.txt", header=T, sep="\t"); colnames(meQTL) <- c("CHR","POS","SNP","PVAL")
jpeg('Figure.3c.jpg', width = 20, height = 5, units = 'in', res = 500, quality = 100)
manhattan(meQTL, chr="CHR", bp="POS", snp="SNP", p="PVAL", col = c("deepskyblue4", "goldenrod3"), suggestiveline = F, genomewideline = F, cex = 0.1, main = "Manhattan Plot")
dev.off()


###########################################################################
#####    Code for Figure.3e                                           #####
###########################################################################
##### Note: Enrichment of tissue-specific meQTL CpGs in tissue enhancers
library(ggplot2)
Tissue.Specific.meQTL.CpG.2HMM <- read.table("Tissue.Specific.meQTL.CpG.2HMM.txt", header=T, sep="\t"); 

pdf(Figure.3d.pdf",width=5,height=3)
ggplot(data = Tissue.Specific.meQTL.CpG.2HMM, aes(x=Enhancer, y=meQTL, fill=OddsRatio)) + 
 geom_tile(color = "white", size = 1)+
 scale_fill_gradient2(low = "#3C5488FF", high = "#DC0000FF", mid = "white", 
   midpoint = 1, limit = c(0.5,1.52), space = "Lab", 
   name="") +
  theme_minimal()
dev.off()
