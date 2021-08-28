###########################################################################
#####    Multiple trait colollization of GWAS, mQTL and eQTL          #####
#####    Software: moloc (v0.1.0)                                     #####
###########################################################################
options(scipen = 1, digits = 2)
library(dplyr)
library(moloc)
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

### Save result
res.moloc <- data.frame(matrix(nrow = 0, ncol = 54))  ### for final results

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
    currentIndex.res.moloc <- data.frame(matrix(nrow = 0, ncol = 54))
    colnames(currentIndex.res.moloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","CpG","mQTL_BETA","mQTL_SE","mQTL_PVAL","mQTL_N","Gene","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N","nsnps","PPA_a","PPA_a_b","PPA_a_c","PPA_a_bc","PPA_a_b_c","PPA_b","PPA_b_c","PPA_ac_b","PPA_c","PPA_ab_c","PPA_ab","PPA_ac","PPA_bc","PPA_abc","PPA_zero","Coloc_ppas_a","Coloc_ppas_b","Coloc_ppas_ab","Coloc_ppas_c","Coloc_ppas_ac","Coloc_ppas_bc","Coloc_ppas_abc","Best_snp_coloc_a","Best_snp_coloc_b","Best_snp_coloc_ab","Best_snp_coloc_c","Best_snp_coloc_ac","Best_snp_coloc_bc","Best_snp_coloc_abc")
	
	### Run each pair
    for (i in 1:Number){
		Current.mQTL <- subset(mQTL, CpG == Pair[i,2])
		Current.eQTL <- subset(eQTL, Gene == Pair[i,3])
		
		CurrentIndex <- IndexSNP[Index,]
    	Current.GWAS.mQTL.eQTL <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));
    	Current.GWAS.mQTL.eQTL <- inner_join(Current.GWAS.mQTL.eQTL, Current.mQTL, by = c("RSID" = "RSID"));
   		Current.GWAS.mQTL.eQTL <- inner_join(Current.GWAS.mQTL.eQTL, Current.eQTL, by = c("RSID" = "RSID"));
    	Current.GWAS.mQTL.eQTL <- na.omit(Current.GWAS.mQTL.eQTL);
		
		### Only do moloc if more than 10 SNP available 
  		if (nrow(Current.GWAS.mQTL.eQTL) >= 50) {
  			### Conbine SNP_CpG_Gene as ID
  			Current.GWAS.mQTL.eQTL$CHR <- paste(Current.GWAS.mQTL.eQTL$IndexRSID, Current.GWAS.mQTL.eQTL$CpG, Current.GWAS.mQTL.eQTL$Gene, sep="|"); colnames(Current.GWAS.mQTL.eQTL)[1] = "Pair"
  			
  			### IndexingSNP
			Current.GWAS.mQTL.eQTL.IndexSNP <- subset(Current.GWAS.mQTL.eQTL, as.character(IndexRSID) == as.character(RSID))
			if (nrow(Current.GWAS.mQTL.eQTL.IndexSNP) >=1){
				Current.res.moloc <- Current.GWAS.mQTL.eQTL.IndexSNP
			}else{ ## Indexing SNP does not have at least one kind of three data (maybe eQTL), use the SNP with lowest p value instead
				Current.res.moloc <- Current.GWAS.mQTL.eQTL[which.min(Current.GWAS.mQTL.eQTL$GWAS_PVAL),]
			}
			Current.res.moloc[,25:54] <- NA
			colnames(Current.res.moloc) <- colnames(currentIndex.res.moloc)
			
			### coloc calculation
  			gwas <- Current.GWAS.mQTL.eQTL[,c("RSID","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","MAF")]; gwas[,7] <- "quant"; colnames(gwas) <- c("SNP","BETA","SE","PVAL","N","MAF","type");
  			mqtl <- Current.GWAS.mQTL.eQTL[,c("RSID","mQTL_BETA","mQTL_SE","mQTL_PVAL","mQTL_N","MAF")]; mqtl[,7] <- "quant"; colnames(mqtl) <- c("SNP","BETA","SE","PVAL","N","MAF","type");
  			eqtl <- Current.GWAS.mQTL.eQTL[,c("RSID","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N","MAF")]; eqtl[,7] <- "quant"; colnames(eqtl) <- c("SNP","BETA","SE","PVAL","N","MAF","type");

  			Current.GWAS.mQTL.eQTL.moloc  <- list(gwas,mqtl,eqtl)
  			moloc <- moloc_test(Current.GWAS.mQTL.eQTL.moloc)       ### used for further analysis (default parameters)
  			
  			# moloc
  			Current.res.moloc[1,25] <- moloc$nsnps
  			Current.res.moloc[1,26:40] <- moloc$priors_lkl_ppa[,4]
  			Current.res.moloc[1,41:47] <- moloc$best_snp[,1]
			Current.res.moloc[1,48:54] <- as.vector(moloc$best_snp[,2])
			currentIndex.res.moloc <- rbind(currentIndex.res.moloc,Current.res.moloc)
		} ## if
	} ### Run each pair
	res.moloc <- rbind(res.moloc,currentIndex.res.moloc)  ## Add current indexSNP result to finial result
} ### Run each indexSNP

### To match previous code with Zscore
### https://rdrr.io/bioc/gCMAP/src/R/CMAPResults-accessors.R
res.moloc$mQTL_Zscore <- zScores(res.moloc$mQTL_PVAL, direction=res.moloc$mQTL_BETA); 
res.moloc$eQTL_Zscore <- zScores(res.moloc$eQTL_PVAL, direction=res.moloc$eQTL_BETA); 
res.moloc <- res.moloc[,c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","CpG","mQTL_Zscore","mQTL_PVAL","mQTL_N","Gene","eQTL_Zscore","eQTL_PVAL","eQTL_N","mQTL_BETA","mQTL_SE","eQTL_BETA","eQTL_SE","nsnps","PPA_a","PPA_a_b","PPA_a_c","PPA_a_bc","PPA_a_b_c","PPA_b","PPA_b_c","PPA_ac_b","PPA_c","PPA_ab_c","PPA_ab","PPA_ac","PPA_bc","PPA_abc","PPA_zero","Coloc_ppas_a","Coloc_ppas_b","Coloc_ppas_ab","Coloc_ppas_c","Coloc_ppas_ac","Coloc_ppas_bc","Coloc_ppas_abc","Best_snp_coloc_a","Best_snp_coloc_b","Best_snp_coloc_ab","Best_snp_coloc_c","Best_snp_coloc_ac","Best_snp_coloc_bc","Best_snp_coloc_abc")]

##### Output final result
write.table(res.moloc, "res.moloc.txt", quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")

