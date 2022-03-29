###########################################################################
#####    Code for Figure.6                                            #####
#####    Hongbo Liu (hongbo919@gmail.com)                             #####
###########################################################################


###########################################################################
#####    Code for Figure.6c                                           #####
###########################################################################
##### Note: Manhattan plot highlighting 330 genes with evidence of multiple traits colocalization.
library(ggman)
GWAS <- read.table("GWAS.forManhattanPlot.txt", header=T, sep="\t")
HighlightSNPs <- read.table("HighlightSNPs.txt", header=T, sep="\t")
HighlightGMESNPs <- read.table("HighlightGMESNPs.txt", header=T, sep="\t")

GWAS[which(GWAS$PVAL < 1e-300),4] = 1e-300

p1 <- ggman(GWAS, sigLine = 7.30103, snp = "RSID", bp = "POS", chrom = "CHR", pvalue = "PVAL", relative.positions = TRUE) + theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2 <- ggmanHighlightGroup(p1, highlightDfm = HighlightSNPs, snp = "RSID", group = "Group", size = 2, alpha = 1, stroke = 0, legend.title = "Coloc Type")
p3 <- ggmanHighlight(p2, highlight = as.character(HighlightGMESNPs$RSID), size = 2, stroke = 0)

jpeg('Figure.6c.jpg', width = 15, height = 8, units = 'in', res = 500, quality = 100)
ggmanLabel(p3, labelDfm = Gene_Highlit, snp = "RSID", label = "Gene", type = "text", colour = "black", fontface = 'italic', size = 3.5,
	nudge_y = 50,
	max.iter = 5000,
    vjust=1,
    direction='y',
    nudge_x=0.1,
    segment.color = "grey50",
    segment.alpha = 0.5,
    segment.size  = 0.2)
dev.off()



