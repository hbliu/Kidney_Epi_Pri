###########################################################################
#####    Determation of cutoff for significant varaints identification#####
#####    Sample size: 443 human kidney samples                        #####
###########################################################################
opt_fdr <- 0.05
### according to the code by cis_permutation_pass.cpp, opt_n is the number of degrees of freedom used to compute the P-values
sample_count = 443
opt_n = sample_count - 2;

D = read.table("FastQTL.permutations.txt.gz", hea=F, stringsAsFactors=F)
colnames(D) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope","ppval", "bpval")
MASK=!is.na(D[,11])
Dnas=D[!MASK,]
D = D[MASK,]
library(qvalue)
D$st = qvalue(D$bpval)$qvalues
dim(D[which(D$st <= 0.05), ])
write.table(D[which(D$st <= 0.05), ], "FastQTL.permutations.storey.txt", quote=F, row.names=F, col.names=T,sep="\t")
Q = qvalue(D[,11]);
D$qval = NA;
D$qval= Q$qvalue;
set0 = D[which(D$qval <= opt_fdr),]
set1 = D[which(D$qval > opt_fdr),]
pthreshold = (sort(set1[,11])[1] - sort(-1.0 * set0[,11])[1]) / 2
pval0 = qbeta(pthreshold, D[,3], D[,4], ncp = 0, lower.tail = TRUE, log.p = FALSE) 
test0 = qf(pval0, 1, D[,5], ncp = 0, lower.tail = FALSE, log.p = FALSE)            
corr0 = sqrt(test0 / (D[,5] + test0))                                          
test1 = opt_n * corr0 * corr0 / (1 - corr0 * corr0)                                
pval1 = pf(test1, 1, opt_n, ncp = 0, lower.tail = FALSE, log.p = FALSE)            
D$nthresholds = pval1
D1=D[, c(1, 14)]
D2=Dnas[, c(1,10)]
names(D2)=names(D1)
D3=rbind(D1, D2)  
write.table(D3, "FastQTL.permutations.Cutoff_eVariant.txt", quote=FALSE, row.names=FALSE, col.names=FALSE,sep="\t")


##### Identification of significant mQTLs using the cutoff in FastQTL.permutations.Cutoff_eVariant.txt for each CpG site
