###########################################################################
#####    Colollization of each two traits among GWAS, mQTL and eQTL   #####
#####    Software: coloc (v5.1.0)                                     #####
###########################################################################
options(scipen = 1, digits = 2)
library(dplyr)
library(coloc)
library(gCMAP)
IndexSNP <- read.table(gzfile("Meta.GWAS.IndexSNP.txt.gz"), header=F, sep="\t"); colnames(IndexSNP) <- c("CHR","IndexSNP","IndexRSID")
GWAS <- read.table(gzfile("Meta.GWAS.txt.gz"), header=F, sep="\t");   colnames(GWAS) <- c("IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N")
mQTL <- read.table(gzfile("Kidney.mQTL.txt.gz"), header=F, sep="\t"); colnames(mQTL) <- c("RSID","CpG","mQTL_BETA","mQTL_SE","mQTL_PVAL","mQTL_N")
eQTL <- read.table(gzfile("Kidney.eQTL.txt.gz"), header=F, sep="\t"); colnames(eQTL) <- c("RSID","Gene","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N")
Sig.mQTL <- read.table(gzfile("Kidney.mQTL.Sig.txt.gz"), header=T, sep="\t")
Sig.eQTL <- read.table(gzfile("Kidney.eQTL.Sig.txt.gz"), header=T, sep="\t")

### Significant mQTLs based on FASTQTL permutation
mQTL_forPair <- mQTL
mQTL_forPair$RSID_CpG <- paste(mQTL_forPair$RSID, mQTL_forPair$CpG, sep="_")
mQTL_forPair <- subset(mQTL_forPair, RSID_CpG %in% Sig.mQTL$RSID_CpG)[,1:6]

### Significant eQTLs based on gene level FDR (q value)
eQTL_forPair <- eQTL
eQTL_forPair$RSID_GeneID <- paste(eQTL_forPair$RSID, eQTL_forPair$Gene, sep="_")
eQTL_forPair <- subset(eQTL_forPair, RSID_GeneID %in% Sig.eQTL$RSID_GeneID)[,1:6]


####################################################################################################
##### Colollization between eGFR GWAS and kidney eQTL 
GWAS.eQTL.coloc <- data.frame(matrix(nrow = 0, ncol = 25))  ### for final results
colnames(GWAS.eQTL.coloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","Gene","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")

### Run each indexSNP
for (Index in 1:nrow(IndexSNP)){
	### Combine all data based on POS because all in the same chromsome
	eQTL_forCurrent <- eQTL_forPair
	
	CurrentIndexSNP <- as.character(IndexSNP[Index,2])
    CurrentIndex <- IndexSNP[Index,]
    CurrentIndex <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));

    CurrentIndex_GWAS_eQTL <- inner_join(CurrentIndex, eQTL_forPair, by = c("RSID" = "RSID"));
    
    if(nrow(CurrentIndex_GWAS_eQTL) == 0) {
    	CurrentIndex_GWAS_eQTL <- inner_join(CurrentIndex, eQTL, by = c("RSID" = "RSID"));
    	Min_PVAL <- min(CurrentIndex_GWAS_eQTL$eQTL_PVAL)
    	Min_PVAL_Gene <- subset(CurrentIndex_GWAS_eQTL, eQTL_PVAL == Min_PVAL)
    	Min_PVAL_Gene <- as.character(Min_PVAL_Gene[1,]$Gene)
    	eQTL_forCurrent <- subset(eQTL, Gene == Min_PVAL_Gene)
    }
    	
	### Combine
    CurrentIndex <- inner_join(CurrentIndex, eQTL_forCurrent, by = c("RSID" = "RSID"));
    CurrentIndex <- na.omit(CurrentIndex);
    
    Pair <- unique(CurrentIndex[,c("IndexRSID","Gene")]); nrow(Pair)   ### Sig eQTL only
     
    Number <- nrow(Pair)
    currentIndex.GWAS.eQTL.coloc <- data.frame(matrix(nrow = 0, ncol = 25))
    colnames(currentIndex.GWAS.eQTL.coloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","Gene","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
	
	### Run each pair
    for (i in 1:Number){
		Current.eQTL <- subset(eQTL, Gene == Pair[i,2])
		
		CurrentIndex <- IndexSNP[Index,]
    	Current.pair <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));
   		Current.pair <- inner_join(Current.pair, Current.eQTL, by = c("RSID" = "RSID"));
    	Current.pair <- na.omit(Current.pair);
    	Current.pair$GWAS_Var <- (Current.pair$GWAS_SE)^2
    	Current.pair$eQTL_Var <- (Current.pair$eQTL_SE)^2
    			
		### Only do coloc if more than 50 SNP available 
  		if (nrow(Current.pair) >= 50) {
  			### Conbine SNP_Gene as ID
  			Current.pair$CHR <- paste(Current.pair$IndexRSID, Current.pair$Gene, sep="|"); colnames(Current.pair)[1] = "Pair"
  			
  			### IndexingSNP
			Current.pair.IndexSNP <- subset(Current.pair, as.character(IndexRSID) == as.character(RSID))
			if (nrow(Current.pair.IndexSNP) >=1){
				Current.GWAS.eQTL.coloc <- Current.pair.IndexSNP[,1:19]
			}else{ ## Indexing SNP does not have at least one kind of three data (maybe eQTL), use the SNP with lowest p value instead
				Current.GWAS.eQTL.coloc <- Current.pair[which.min(Current.pair$GWAS_PVAL),1:19]
			}
			Current.GWAS.eQTL.coloc[,20:25] <- NA
			colnames(Current.GWAS.eQTL.coloc) <- colnames(currentIndex.GWAS.eQTL.coloc)
						
			### coloc calculation
  			Current.gwas <- Current.pair[,c("RSID","POS","GWAS_BETA","GWAS_Var","GWAS_PVAL","GWAS_N","MAF")]; colnames(Current.gwas) <- c("snp","position","beta","varbeta","pvalues","N","MAF"); Current.gwas$snp <- as.character(Current.gwas$snp)
  			Current.eqtl <- Current.pair[,c("RSID","POS","eQTL_BETA","eQTL_Var","eQTL_PVAL","eQTL_N","MAF")]; colnames(Current.eqtl) <- c("snp","position","beta","varbeta","pvalues","N","MAF"); Current.eqtl$snp <- as.character(Current.eqtl$snp)
						
			Current.gwas <- as.list(as.data.frame(Current.gwas)); Current.gwas[["type"]] <- "quant"; str(Current.gwas)
			Current.eqtl <- as.list(as.data.frame(Current.eqtl)); Current.eqtl[["type"]] <- "quant"; str(Current.eqtl)
  			coloc <- coloc.abf(dataset1 = Current.gwas, dataset2 = Current.eqtl)     ### used for further analysis (default parameters)
  			
  			### coloc.Result
  			Current.GWAS.eQTL.coloc[1,20:25] <- coloc$summary
			currentIndex.GWAS.eQTL.coloc <- rbind(currentIndex.GWAS.eQTL.coloc,Current.GWAS.eQTL.coloc)
		} ## if
	} ### Run each pair
	GWAS.eQTL.coloc <- rbind(GWAS.eQTL.coloc,currentIndex.GWAS.eQTL.coloc)  ## Add current indexSNP result to finial result
} ### Run each indexSNP

##### Output final result
write.table(GWAS.eQTL.coloc, "GWAS.eQTL.coloc.txt", quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")



####################################################################################################
##### Colollization between eGFR GWAS and kidney mQTL 
GWAS.mQTL.coloc <- data.frame(matrix(nrow = 0, ncol = 25))  ### for final results
colnames(GWAS.mQTL.coloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","CpG","mQTL_BETA","mQTL_SE","mQTL_PVAL","mQTL_N","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")

### Run each indexSNP
for (Index in 1:nrow(IndexSNP)){
	### Combine all data based on POS because all in the same chromsome
	mQTL_forCurrent <- mQTL_forPair
	
	CurrentIndexSNP <- as.character(IndexSNP[Index,2])
    CurrentIndex <- IndexSNP[Index,]
    CurrentIndex <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));

    CurrentIndex_GWAS_mQTL <- inner_join(CurrentIndex, mQTL_forPair, by = c("RSID" = "RSID"));
    
    if(nrow(CurrentIndex_GWAS_mQTL) == 0) {
    	CurrentIndex_GWAS_mQTL <- inner_join(CurrentIndex, mQTL, by = c("RSID" = "RSID"));
    	Min_PVAL <- min(CurrentIndex_GWAS_mQTL$mQTL_PVAL)
    	Min_PVAL_CpG <- subset(CurrentIndex_GWAS_mQTL, mQTL_PVAL == Min_PVAL)
    	Min_PVAL_CpG <- as.character(Min_PVAL_CpG[1,]$CpG)
    	mQTL_forCurrent <- subset(mQTL, CpG == Min_PVAL_CpG)
    }
    	
	### Combine
    CurrentIndex <- inner_join(CurrentIndex, mQTL_forCurrent, by = c("RSID" = "RSID"));
    CurrentIndex <- na.omit(CurrentIndex);
    
    Pair <- unique(CurrentIndex[,c("IndexRSID","CpG")]); nrow(Pair)   ### Sig mQTL only
     
    Number <- nrow(Pair)
    currentIndex.GWAS.mQTL.coloc <- data.frame(matrix(nrow = 0, ncol = 25))
    colnames(currentIndex.GWAS.mQTL.coloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","CpG","mQTL_BETA","mQTL_SE","mQTL_PVAL","mQTL_N","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
	
	### Run each pair
    for (i in 1:Number){
		Current.mQTL <- subset(mQTL, CpG == Pair[i,2])
		
		CurrentIndex <- IndexSNP[Index,]
    	Current.pair <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));
   		Current.pair <- inner_join(Current.pair, Current.mQTL, by = c("RSID" = "RSID"));
    	Current.pair <- na.omit(Current.pair);
    	Current.pair$GWAS_Var <- (Current.pair$GWAS_SE)^2
    	Current.pair$mQTL_Var <- (Current.pair$mQTL_SE)^2
    			
		### Only do coloc if more than 50 SNP available 
  		if (nrow(Current.pair) >= 50) {
  			### Conbine SNP_CpG as ID
  			Current.pair$CHR <- paste(Current.pair$IndexRSID, Current.pair$CpG, sep="|"); colnames(Current.pair)[1] = "Pair"
  			
  			### IndexingSNP
			Current.pair.IndexSNP <- subset(Current.pair, as.character(IndexRSID) == as.character(RSID))
			if (nrow(Current.pair.IndexSNP) >=1){
				Current.GWAS.mQTL.coloc <- Current.pair.IndexSNP[,1:19]
			}else{ ## Indexing SNP does not have at least one kind of three data (maybe mQTL), use the SNP with lowest p value instead
				Current.GWAS.mQTL.coloc <- Current.pair[which.min(Current.pair$GWAS_PVAL),1:19]
			}
			Current.GWAS.mQTL.coloc[,20:25] <- NA
			colnames(Current.GWAS.mQTL.coloc) <- colnames(currentIndex.GWAS.mQTL.coloc)
						
			### coloc calculation
  			Current.gwas <- Current.pair[,c("RSID","POS","GWAS_BETA","GWAS_Var","GWAS_PVAL","GWAS_N","MAF")]; colnames(Current.gwas) <- c("snp","position","beta","varbeta","pvalues","N","MAF"); Current.gwas$snp <- as.character(Current.gwas$snp)
  			Current.mqtl <- Current.pair[,c("RSID","POS","mQTL_BETA","mQTL_Var","mQTL_PVAL","mQTL_N","MAF")]; colnames(Current.mqtl) <- c("snp","position","beta","varbeta","pvalues","N","MAF"); Current.mqtl$snp <- as.character(Current.mqtl$snp)
						
			Current.gwas <- as.list(as.data.frame(Current.gwas)); Current.gwas[["type"]] <- "quant"; str(Current.gwas)
			Current.mqtl <- as.list(as.data.frame(Current.mqtl)); Current.mqtl[["type"]] <- "quant"; str(Current.mqtl)
  			coloc <- coloc.abf(dataset1 = Current.gwas, dataset2 = Current.mqtl)     ### used for further analysis (default parameters)
  			
  			### coloc.Result
  			Current.GWAS.mQTL.coloc[1,20:25] <- coloc$summary
			currentIndex.GWAS.mQTL.coloc <- rbind(currentIndex.GWAS.mQTL.coloc,Current.GWAS.mQTL.coloc)
		} ## if
	} ### Run each pair
	GWAS.mQTL.coloc <- rbind(GWAS.mQTL.coloc,currentIndex.GWAS.mQTL.coloc)  ## Add current indexSNP result to finial result
} ### Run each indexSNP

##### Output final result
write.table(GWAS.mQTL.coloc, "GWAS.mQTL.coloc.txt", quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")



####################################################################################################
##### Colollization between kidney mQTL and kidney eQTL 
mQTL.eQTL.coloc <- data.frame(matrix(nrow = 0, ncol = 30))  ### for final results
colnames(mQTL.eQTL.coloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","CpG","mQTL_BETA","mQTL_SE","mQTL_PVAL","mQTL_N","Gene","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")

### Run each indexSNP
for (Index in 1:nrow(IndexSNP)){
	### Combine all data based on POS because all in the same chromsome
	mQTL_forCurrent <- mQTL_forPair
	eQTL_forCurrent <- eQTL_forPair
	
	CurrentIndexSNP <- as.character(IndexSNP[Index,2])
    CurrentIndex <- IndexSNP[Index,]
    CurrentIndex <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));
    
    CurrentIndex_GWAS_mQTL <- inner_join(CurrentIndex, mQTL_forPair, by = c("RSID" = "RSID"));
    CurrentIndex_GWAS_eQTL <- inner_join(CurrentIndex, eQTL_forPair, by = c("RSID" = "RSID"));
    
    if(nrow(CurrentIndex_GWAS_eQTL) == 0) {
    	CurrentIndex_GWAS_eQTL <- inner_join(CurrentIndex, eQTL, by = c("RSID" = "RSID"));
    	Min_PVAL <- min(CurrentIndex_GWAS_eQTL$eQTL_PVAL)
    	Min_PVAL_Gene <- subset(CurrentIndex_GWAS_eQTL, eQTL_PVAL == Min_PVAL)
    	Min_PVAL_Gene <- as.character(Min_PVAL_Gene[1,]$Gene)
    	eQTL_forCurrent <- subset(eQTL, Gene == Min_PVAL_Gene)
    }
    
    if(nrow(CurrentIndex_GWAS_mQTL) == 0) {
    	CurrentIndex_GWAS_mQTL <- inner_join(CurrentIndex, mQTL, by = c("RSID" = "RSID"));
    	Min_PVAL <- min(CurrentIndex_GWAS_mQTL$mQTL_PVAL)
    	Min_PVAL_CpG <- subset(CurrentIndex_GWAS_mQTL, mQTL_PVAL == Min_PVAL)
    	Min_PVAL_CpG <- as.character(Min_PVAL_CpG[1,]$CpG)
    	mQTL_forCurrent <- subset(mQTL, CpG == Min_PVAL_CpG)
    }
	
	### Combine
    CurrentIndex <- inner_join(CurrentIndex, mQTL_forCurrent, by = c("RSID" = "RSID"));
    CurrentIndex <- inner_join(CurrentIndex, eQTL_forCurrent, by = c("RSID" = "RSID"));
    CurrentIndex <- na.omit(CurrentIndex);
    
    Pair <- unique(CurrentIndex[,c("IndexRSID","CpG","Gene")]); nrow(Pair)   ### Sig mQTL and eQTL only
     
    Number <- nrow(Pair)
    currentIndex.mQTL.eQTL.coloc <- data.frame(matrix(nrow = 0, ncol = 30))
    colnames(currentIndex.mQTL.eQTL.coloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","CpG","mQTL_BETA","mQTL_SE","mQTL_PVAL","mQTL_N","Gene","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
	
	### Run each pair
    for (i in 1:Number){
		Current.mQTL <- subset(mQTL, CpG == Pair[i,2])
		Current.eQTL <- subset(eQTL, Gene == Pair[i,3])
		
		CurrentIndex <- IndexSNP[Index,]
    	Current.GWAS.mQTL.eQTL <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));
    	Current.GWAS.mQTL.eQTL <- inner_join(Current.GWAS.mQTL.eQTL, Current.mQTL, by = c("RSID" = "RSID"));
   		Current.GWAS.mQTL.eQTL <- inner_join(Current.GWAS.mQTL.eQTL, Current.eQTL, by = c("RSID" = "RSID"));
    	Current.GWAS.mQTL.eQTL <- na.omit(Current.GWAS.mQTL.eQTL);
    	Current.GWAS.mQTL.eQTL$mQTL_Var <- (Current.GWAS.mQTL.eQTL$mQTL_SE)^2
    	Current.GWAS.mQTL.eQTL$eQTL_Var <- (Current.GWAS.mQTL.eQTL$eQTL_SE)^2
    			
		### Only do coloc if more than 50 SNP available 
  		if (nrow(Current.GWAS.mQTL.eQTL) >= 50) {
  			### Conbine SNP_CpG as ID
  			Current.GWAS.mQTL.eQTL$CHR <- paste(Current.GWAS.mQTL.eQTL$IndexRSID, Current.GWAS.mQTL.eQTL$CpG, Current.GWAS.mQTL.eQTL$Gene, sep="|"); colnames(Current.GWAS.mQTL.eQTL)[1] = "Pair"
  			
  			### IndexingSNP
			Current.GWAS.mQTL.eQTL.IndexSNP <- subset(Current.GWAS.mQTL.eQTL, as.character(IndexRSID) == as.character(RSID))
			if (nrow(Current.GWAS.mQTL.eQTL.IndexSNP) >=1){
				Current.mQTL.eQTL.coloc <- Current.GWAS.mQTL.eQTL.IndexSNP[,1:24]
			}else{ ## Indexing SNP does not have at least one kind of three data (maybe eQTL), use the SNP with lowest p value instead
				Current.mQTL.eQTL.coloc <- Current.GWAS.mQTL.eQTL[which.min(Current.GWAS.mQTL.eQTL$GWAS_PVAL),1:24]
			}
			Current.mQTL.eQTL.coloc[,25:30] <- NA
			colnames(Current.mQTL.eQTL.coloc) <- colnames(currentIndex.mQTL.eQTL.coloc)
						
			### coloc calculation
  			Current.mqtl <- Current.GWAS.mQTL.eQTL[,c("RSID","POS","mQTL_BETA","mQTL_Var","mQTL_PVAL","mQTL_N","MAF")]; colnames(Current.mqtl) <- c("snp","position","beta","varbeta","pvalues","N","MAF"); Current.mqtl$snp <- as.character(Current.mqtl$snp)
  			Current.eqtl <- Current.GWAS.mQTL.eQTL[,c("RSID","POS","eQTL_BETA","eQTL_Var","eQTL_PVAL","eQTL_N","MAF")]; colnames(Current.eqtl) <- c("snp","position","beta","varbeta","pvalues","N","MAF"); Current.eqtl$snp <- as.character(Current.eqtl$snp)
						
			Current.mqtl <- as.list(as.data.frame(Current.mqtl)); Current.mqtl[["type"]] <- "quant"; str(Current.mqtl)
			Current.eqtl <- as.list(as.data.frame(Current.eqtl)); Current.eqtl[["type"]] <- "quant"; str(Current.eqtl)
  			coloc <- coloc.abf(dataset1 = Current.mqtl, dataset2 = Current.eqtl)     ### used for further analysis (default parameters)
  			
  			### coloc.Result
  			Current.mQTL.eQTL.coloc[1,25:30] <- coloc$summary
			currentIndex.mQTL.eQTL.coloc <- rbind(currentIndex.mQTL.eQTL.coloc,Current.mQTL.eQTL.coloc)
		} ## if
	} ### Run each pair
	mQTL.eQTL.coloc <- rbind(mQTL.eQTL.coloc,currentIndex.mQTL.eQTL.coloc)  ## Add current indexSNP result to finial result
} ### Run each indexSNP

##### Output final result
write.table(mQTL.eQTL.coloc, "mQTL.eQTL.coloc.txt", quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")


