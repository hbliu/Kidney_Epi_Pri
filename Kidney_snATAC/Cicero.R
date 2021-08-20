###########################################################################
#####    Single cell GWAS enrichment by gchromVAR and snATAC-seq      #####
#####    GWAS: UKBB GWAS fine mapping (overlap mQTL SNPs)             #####
#####    Software: Cicero (v1.0.15)                                   #####
###########################################################################
library(cicero)
### Loading data from a simple sparse matrix format converted from SnapATAC
input_cds <- make_atac_cds("Kidney.snATAC.pmat.triples.txt", binarize = TRUE)

### Create a Cicero CDS
set.seed(2017)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6, reduction_method = 'tSNE', norm_method = "none", check_duplicates = FALSE)
tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

### Run Cicero
data("human.hg19.genome")
conns <- run_cicero(cicero_cds, human.hg19.genome)
Kidney.snATAC.Cicero.conns <- subset(conns, coaccess > 0.2)
write.table(Kidney.snATAC.Cicero.conns,"Kidney.snATAC.Cicero.conns.txt", quote=FALSE, row.names=FALSE, col.names=T, sep="\t")
