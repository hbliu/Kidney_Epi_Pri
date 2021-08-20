###########################################################################
##### Target genes prioritization for eGFR GWAS SNPs using multiple data ##
###########################################################################
##### Step1: Collect 8 Datasets for Gene Prioritization for eGFR GWAS SNPs
# (1) Significant SNP~gene associations by kidney eQTL (FDR < 0.05); 
# (2) Significant SNP~CpG~gene associations by kidney mQTL (FDR < 0.05) and eQTM (CpG level FDR <0.05); 
# (3) SNP~gene pairs by coloc analysis between eGFR GWAS and eQTL (H4 > 0.8); 
# (4) SNP~gene pairs by moloc analysis among eGFRGWAS, eQTL and mQTL (PPA.abc > 0.8); 
# (5) Significant SNP~gene pairs by mendelian randomization analysis between eGFR GWAS and eQTL (PSMR < 1.38×10-4); 
# (6) SNP~gene pairs passing HEIDI test between eGFR GWAS and eQTL (PHEIDI > 0.01); 
# (7) Cicero connections identified using 57,262 snATAC-seq cells (co-accessibility score > 0.2); 
# (8) Enhancer-promoter contacts identified by Activity-by-Contact (ABC) Model which predicts enhancers regulating genes based on estimating enhancer activity and enhancer-promoter contact frequency from epigenomic datasets (ABC scores >= 0.015).


###########################################################################
##### Step2: Extract all potential target genes within 1Mb for each significant eGFR GWAS SNP
cp Meta.eGFR.GWAS.CKDGen.UKB.MVP.PAGE.SUMMIT.2Study.autosome.MAF0.001.Sig5e8.NoMHC.txt GWAS.Meta.Sig5e8.txt
awk 'NR>1{print "chr"$2"\t"$3-1"\t"$3"\t"$1}' GWAS.Meta.Sig5e8.txt | sort -t $'\t' -k1,1 -k2,2n -V > GWAS.Meta.Sig5e8.bed

awk 'NR>1{if($3-1000000 >= 0) print "chr"$2"\t"$3-1000000"\t"$3+1000000"\tchr"$2"\t"$3-1"\t"$3"\t"$1; else if($3-1000000 < 0) print "chr"$2"\t"0"\t"$3+1000000"\tchr"$2"\t"$3-1"\t"$3"\t"$1}' GWAS.Meta.Sig5e8.txt |\
sort -t $'\t' -k1,1 -k2,2n -V > GWAS.Meta.Sig5e8.1M.bed

bedtools intersect -a GWAS.Meta.Sig5e8.1M.bed -b /home/hongbol/Genome_Annotation/Reference/gtf/hg19/CodingGene_Transcript_gencode_v35lift37_hg19_sorted.TSS.autosome.bed -wa -wb |\
awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$6-$10}' > GWAS.Meta.Sig5e8.1M.AllPotential.bed

head -n 1 GWAS.Meta.Sig5e8.1M.AllPotential.bed | awk '{print "chr\tSNP.Loc\tRSID\tGeneID\tGeneSymbol\tPair_ID\tPair_ID1\tPair_ID2"}' > GWAS.Meta.Sig5e8.1M.AllPotential.txt
awk '{print $1"\t"$3"\t"$4"\t"$9"\t"$10"\t"$4"_"$9"_"$10"\t"$4"_"$9"\t"$4"_"$10}' GWAS.Meta.Sig5e8.1M.AllPotential.bed | sort | uniq >> GWAS.Meta.Sig5e8.1M.AllPotential.txt


###########################################################################
##### Step3: Combine SNP~Gene pairs
R
library(dplyr)
library(tidyr)
AllPotential <- read.table("GWAS.Meta.Sig5e8.1M.AllPotential.txt", header=T, sep="\t")
AllPotential.Infor <- AllPotential

#### eQTL gene
eQTL <- read.table("Kidney.eQTL.MAF0.05.Sig.Gene.q0.05.Pair.txt", header=T, sep="\t")

#### mQTL gene
mQTL <- read.table("mQTL.Sig.eQTM.top1.Pair.txt", header=T, sep="\t")

#### coloc gene
coloc <- read.table("res.coloc.GWAS.eQTL.PP4.0.8.Pair.txt", header=T, sep="\t")

#### moloc gene
moloc <- read.table("res.moloc.GWAS.mQTL.eQTL.ppas_abc.0.8.Pair.txt", header=T, sep="\t")

#### SMR gene
SMR <- read.table("eGFR.GWAS.Kidney.eQTL.smr.coloc.SMR.txt", header=T, sep="\t")

#### SMR_HEIDI gene
SMR_HEIDI <- read.table("eGFR.GWAS.Kidney.eQTL.smr.coloc.SMR_HEIDI.txt", header=T, sep="\t")

#### Cicero gene
Cicero <- read.table(gzfile("GWAS.Meta.Sig5e8.1M.AllPotential.TSS2kbEach.Cicero.conns.coaccess0.2.Uniq.pair.txt.gz"), header=T, sep="\t")

#### ABC gene
ABC <- read.table("ABC.Enhancer.eGFR.Meta5GWAS.SNP2Gene.bed", header=T, sep="\t")

########################################################################
AllPotential.Infor <- left_join(AllPotential.Infor, eQTL, by = c("Pair_ID" = "RSID_GeneID_GeneSymbol"))
AllPotential.Infor <- left_join(AllPotential.Infor, mQTL, by = c("Pair_ID1" = "RSID_Gene"))
AllPotential.Infor <- left_join(AllPotential.Infor, coloc, by = c("Pair_ID" = "Pair_GeneSymbol"))
AllPotential.Infor <- left_join(AllPotential.Infor, moloc, by = c("Pair_ID" = "IndexRSID_Gene_GeneSymbol"))
AllPotential.Infor <- left_join(AllPotential.Infor, SMR, by = c("Pair_ID1" = "GE.Pair"))
AllPotential.Infor <- left_join(AllPotential.Infor, SMR_HEIDI, by = c("Pair_ID1" = "GE.Pair"))
AllPotential.Infor <- left_join(AllPotential.Infor, Cicero, by = c("Pair_ID" = "Pair_ID"))
AllPotential.Infor <- left_join(AllPotential.Infor, ABC, by = c("Pair_ID2" = "RSID_GeneSymbol"))

#### Summary
AllPotential.Infor.AtLeast1 <- subset(AllPotential.Infor, !is.na(eQTL_Pvalue) | !is.na(mQTL_pvalue) | !is.na(PP.H4.abf) | !is.na(Coloc_ppas_abc) | !is.na(GE.p_SMR) | !is.na(GE.p_HEIDI) | !is.na(Cicero_Coaccess) | !is.na(ABC.Score))

#### Score
AllPotential.Infor.AtLeast1$eQTL <- 0
AllPotential.Infor.AtLeast1$mQTL_eQTM <- 0
AllPotential.Infor.AtLeast1$coloc <- 0
AllPotential.Infor.AtLeast1$moloc <- 0
AllPotential.Infor.AtLeast1$SMR <- 0
AllPotential.Infor.AtLeast1$SMR_HEIDI <- 0
AllPotential.Infor.AtLeast1$Coaccess <- 0
AllPotential.Infor.AtLeast1$ABC <- 0

AllPotential.Infor.AtLeast1[which(!is.na(AllPotential.Infor.AtLeast1$eQTL_Pvalue)), "eQTL"] <- 1 
AllPotential.Infor.AtLeast1[which(!is.na(AllPotential.Infor.AtLeast1$mQTL_pvalue)), "mQTL_eQTM"] <- 1 
AllPotential.Infor.AtLeast1[which(!is.na(AllPotential.Infor.AtLeast1$PP.H4.abf)), "coloc"] <- 1 
AllPotential.Infor.AtLeast1[which(!is.na(AllPotential.Infor.AtLeast1$Coloc_ppas_abc)), "moloc"] <- 1 
AllPotential.Infor.AtLeast1[which(!is.na(AllPotential.Infor.AtLeast1$GE.p_SMR)), "SMR"] <- 1 
AllPotential.Infor.AtLeast1[which(!is.na(AllPotential.Infor.AtLeast1$GE.p_HEIDI)), "SMR_HEIDI"] <- 1 
AllPotential.Infor.AtLeast1[which(!is.na(AllPotential.Infor.AtLeast1$Cicero_Coaccess)), "Coaccess"] <- 1 
AllPotential.Infor.AtLeast1[which(!is.na(AllPotential.Infor.AtLeast1$ABC.Score)), "ABC"] <- 1 

AllPotential.Infor.AtLeast1$Priority_Score <- AllPotential.Infor.AtLeast1$eQTL + AllPotential.Infor.AtLeast1$mQTL_eQTM + AllPotential.Infor.AtLeast1$coloc + AllPotential.Infor.AtLeast1$moloc + AllPotential.Infor.AtLeast1$SMR + AllPotential.Infor.AtLeast1$SMR_HEIDI + AllPotential.Infor.AtLeast1$Coaccess + AllPotential.Infor.AtLeast1$ABC
write.table(AllPotential.Infor.AtLeast1,"AllPotential.Infor.AtLeast1.txt", col.names=T, row.names = F, quote=F,sep="\t")


################ Select target genes with priority score >= 3
AllPotential.Infor.AtLeast3 <- subset(AllPotential.Infor.AtLeast1, Priority_Score >= 3)
write.table(AllPotential.Infor.AtLeast3,"AllPotential.Infor.AtLeast3.txt", col.names=T, row.names = F, quote=F,sep="\t")


