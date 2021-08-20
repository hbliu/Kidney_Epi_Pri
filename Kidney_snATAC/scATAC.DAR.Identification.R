###########################################################################
#####    Identification of differentially accessible regions (DARs)   #####
#####    by one-sided Fisher’s exact test for each cluster            #####
#####    Following is the code for cluster 1 (PT_S1)                  #####
#####    Software: SnapATAC (v2.0)                                    #####
###########################################################################
library(SnapATAC)
load("SnapATAC.object.RData")
load("pmat.RData")
load("peak.gr.RData")

### PT_S1.cluster.1
TestCluster = "1"
idy.PT_S1.cluster.1 = lapply(levels(x.after.sp@cluster), function(cluster_i){
	DARs = findDAR(
		obj=x.after.sp,
		input.mat="pmat",
		cluster.pos= TestCluster,
		cluster.neg=cluster_i,
		bcv=0.1,
		test.method="exactTest",
		seed.use=10
		);
	DARs$FDR = p.adjust(DARs$PValue, method="BH");
	idy = which(DARs$FDR < 0.05 & DARs$logFC > 0);
	if((x=length(idy)) < 2000L){
			PValues = DARs$PValue;
			PValues[DARs$logFC < 0] = 1;
			idy = order(PValues, decreasing=FALSE)[1:2000];
			rm(PValues); # free memory
	}
	idy
  })
### 
names(idy.PT_S1.cluster.1) = levels(x.after.sp@cluster);

idy.PT_S1.cluster.1.used <- idy.PT_S1.cluster.1[setdiff(c("1","2","3","4","5","6","8","9","10","11","12","13"),c(TestCluster))]
idy.PT_S1.cluster.1.Count <- as.data.frame(table(unlist(idy.PT_S1.cluster.1.used, recursive = FALSE)));
colnames(idy.PT_S1.cluster.1.Count) <- c("PeakNum","DAR_Count")
PT_S1.cluster.1.DAR.Count <- cbind(data.frame(peak.gr[as.numeric(as.character(idy.PT_S1.cluster.1.Count$PeakNum)),]), idy.PT_S1.cluster.1.Count)
PT_S1.cluster.1.DAR.Count$seqnames <- paste("chr",PT_S1.cluster.1.DAR.Count$seqnames, sep = "")

PT_S1.cluster.1.DAR.7 <- subset(PT_S1.cluster.1.DAR.Count, DAR_Count >= 7)

### output
write.table(PT_S1.cluster.1.DAR.7[,c(1,2,3,7)],"PT_S1.cluster.1.DAR.7.bed",col.names=F,row.names=F,quote=F,sep="\t")


