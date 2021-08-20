###########################################################################
#####    Meta-analysis of five eGFR GWAS studies                      #####
#####    Kidney eQTLs by Sheng et al. (n=356)                         #####
#####    Kidney eQTLs by Ko et al. (n=91)                             #####
#####    Kidney eQTLs by GTEx v8 (n=73)                               #####
#####    Kidney eQTLs by NephQTL (n=166)                              #####
#####    Software: METAL (version 2011-03-25)                         #####
###########################################################################
##### Step1: Run meta-analysis using metal
metal Metal.Parameters.txt

###########################################################################
##### Step2: Association summary statistics were canonicalized to make sure effect size was always reported with respect to the alternate allele as above. 
R
library(dplyr)
library(tidyr)
library(stringr)
G1000SNP <- read.table(gzfile("./PLINK2_1000GP_Phase3/EUR_phase3.checkRef.rmdup.vcf.RSID.txt.gz"), header=F, sep="\t");
colnames(G1000SNP) <- c("SNP","RSID","G1000_REF","G1000_ALT")
G1000SNP.Freq <- read.table(gzfile("./PLINK2_1000GP_Phase3/EUR_phase3.SNPTEST.formatted.gz"), header=T, sep="\t");
G1000SNP.Freq <- G1000SNP.Freq[,c("RSID","MAF")]
G1000SNP <- left_join(G1000SNP, G1000SNP.Freq, by = c("RSID" = "RSID"))

### Add MAF infor
Meta.Kidney.eQTL <- read.table("METAANALYSIS1.TBL.txt", header=T, sep="\t")
Meta.Kidney.eQTL <- left_join(Meta.Kidney.eQTL, G1000SNP, by = c("RSID" = "RSID"))

### Reverse flipping
Meta.Kidney.eQTL.NoNA.REF.Effect <- subset(Meta.Kidney.eQTL.NoNA, as.character(Allele1) == as.character(G1000_REF))  ### need to reverse
Meta.Kidney.eQTL.NoNA.ALT.Effect <- subset(Meta.Kidney.eQTL.NoNA, as.character(Allele2) == as.character(G1000_REF))  ### Not need to reverse
Meta.Kidney.eQTL.NoNA.REF.Effect$Effect_Harmonized <- Meta.Kidney.eQTL.NoNA.REF.Effect$Effect * -1
Meta.Kidney.eQTL.NoNA.REF.Effect$Direction_Harmonized <- str_replace_all(Meta.Kidney.eQTL.NoNA.REF.Effect$Direction, "[+]", "*")
Meta.Kidney.eQTL.NoNA.REF.Effect$Direction_Harmonized <- str_replace_all(Meta.Kidney.eQTL.NoNA.REF.Effect$Direction_Harmonized, "[-]", "+")
Meta.Kidney.eQTL.NoNA.REF.Effect$Direction_Harmonized <- str_replace_all(Meta.Kidney.eQTL.NoNA.REF.Effect$Direction_Harmonized, "[*]", "-")
Meta.Kidney.eQTL.NoNA.ALT.Effect$Effect_Harmonized <- Meta.Kidney.eQTL.NoNA.ALT.Effect$Effect
Meta.Kidney.eQTL.NoNA.ALT.Effect$Direction_Harmonized <- Meta.Kidney.eQTL.NoNA.ALT.Effect$Direction

### Combine
Meta.Kidney.eQTL.Flip_Reversed <- rbind(Meta.Kidney.eQTL.NoNA.REF.Effect, Meta.Kidney.eQTL.NoNA.ALT.Effect)

### Add Gene Symbol
Gene <- read.table(gzfile("./Gene_gencode_v35lift37_hg19_ID.txt.gz"), header=F, sep="\t")
Gene <- Gene %>% separate(V1, c("GeneID","version"), sep = "([.])")
Gene.Simple <- Gene[,c(1,6,7)]
colnames(Gene.Simple) <- c("GeneID","GeneSymbol","GeneType")
Meta.Kidney.eQTL.Flip_Reversed.infor <- left_join(Meta.Kidney.eQTL.Flip_Reversed, Gene.Simple, by = c("GeneID" = "GeneID"))
Meta.Kidney.eQTL.Flip_Reversed.infor.WithSymbol <- subset(Meta.Kidney.eQTL.Flip_Reversed.infor, !is.na(GeneSymbol))
Meta.Kidney.eQTL.Flip_Reversed.infor.WithoutSymbol <- subset(Meta.Kidney.eQTL.Flip_Reversed.infor, is.na(GeneSymbol))
Meta.Kidney.eQTL.Flip_Reversed.infor.WithoutSymbol$GeneSymbol <- Meta.Kidney.eQTL.Flip_Reversed.infor.WithoutSymbol$GeneID
Meta.Kidney.eQTL.Flip_Reversed.infor <- rbind(Meta.Kidney.eQTL.Flip_Reversed.infor.WithSymbol, Meta.Kidney.eQTL.Flip_Reversed.infor.WithoutSymbol)

### Save
write.table(Meta.Kidney.eQTL.Flip_Reversed.infor, "METAANALYSIS1.TBL.ALT.Effect.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


###########################################################################
##### Step3: Only keep associations with MAF > 0.05
awk '{if($7 == "MAF" || $7>0.05) print $0}' METAANALYSIS1.TBL.ALT.Effect.txt > Kidney.eQTL.Meta.Analsyis.MAF0.05.txt


###########################################################################
##### Step4: Identification of significant eQTLs by q values calculated by the Storey approach
R
library(dplyr)
library(qvalue)
### Read eQTLs
Meta.Kidney.eQTL <- read.table("Kidney.eQTL.Meta.Analsyis.MAF0.05.txt",header=T, sep="\t")
### Calculate qvalue
Meta.Kidney.eQTL <- Meta.Kidney.eQTL %>% group_by(GeneSymbol) %>% mutate(Gene.qvalue=qvalue(PVAL, pi0.method="bootstrap")$qvalues) %>% as.data.frame  
### Significant eQTLs
Meta.Kidney.eQTL.Sig.Gene.q0.05 <- subset(Meta.Kidney.eQTL, Gene.qvalue < 0.05)
### Save
write.table(Meta.Kidney.eQTL.Sig.Gene.q0.05, "Meta.Kidney.eQTL.Sig.Gene.q0.05.txt", quote=FALSE, row.names=F, col.names=TRUE, sep="\t")


