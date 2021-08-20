###########################################################################
#####    Meta-analysis of five eGFR GWAS studies                      #####
#####    CKDGen eGFR GWAS by Wuttke et al. (n=765,348)                #####
#####    UKBB eGFR GWAS by Sinnott-Armstrong et al. (n=355,731)       #####
#####    MVP eGFR GWAS by Hellwege et al. (n=292,258)                 #####
#####    PAGE eGFR GWAS by Wojcik et al. (n=27,900)                   #####
#####    SUMMIT eGFR GWAS by van Zuydam et al. (n=13,158)             #####
#####    Software: METAL (version 2011-03-25)                         #####
###########################################################################
##### Step1: Run meta-analysis using metal
metal Metal.Parameters.txt


###########################################################################
##### Step2: Extract associations tested in at least two sutdies  
awk '{if($10 == "P-value" || ($14 > 0)) print $0}' METAANALYSIS1.TBL > Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.2Study.txt


###########################################################################
##### Step3: Only keep associations from autosomes and MAF >= 0.001
cat CHROM 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 > Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.2Study.autosome.txt
awk '{if($20 == "MAF" || $20 >= 0.001) print $0}' Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.2Study.autosome.txt > Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.2Study.autosome.MAF0.001.txt


###########################################################################
##### Step4: Association summary statistics were canonicalized to make sure effect size was always reported with respect to the alternate allele as above. 
R
library(dplyr)
library(tidyr)
library(stringr)
######### Check Metal.eGFR (Metal.CKDGen.UKB.MVP.PAGE.SUMMIT)
Metal.CKDGen.UKB.MVP.PAGE.SUMMIT <- read.table("Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.2Study.autosome.MAF0.001.txt", header=T, sep="\t")
SNP_b151_All <- read.table("dbSNP_build_151.txt", header=T, sep="\t");

#### Add SNP information from dbSNP_build_151
Metal.CKDGen.UKB.MVP.PAGE.SUMMIT <- left_join(Metal.CKDGen.UKB.MVP.PAGE.SUMMIT, SNP_b151_All, by = c("MarkerName" = "ID"))

### Check REF allele
EffectSize <- Metal.CKDGen.UKB.MVP.PAGE.SUMMIT
EffectSize$Allele1 <- toupper(EffectSize$Allele1)
EffectSize$Allele2 <- toupper(EffectSize$Allele2)

### Allele1 is the effect ellele
### Reverse effect size of GWAS if Allele1 is REF
EffectSize.Canonicalized <- subset(EffectSize, as.character(Allele1) == as.character(REF))
EffectSize.Canonicalized <- as.data.frame(EffectSize.Canonicalized)
EffectSize.Canonicalized$RSID_REF <- EffectSize.Canonicalized$Allele1; EffectSize.Canonicalized$RSID_ALT <- EffectSize.Canonicalized$Allele2; EffectSize.Canonicalized$Fre_ALT <- (1 - EffectSize.Canonicalized$Freq1)
EffectSize.Canonicalized$Zscore_Canonicalized <- EffectSize.Canonicalized$Zscore * -1
EffectSize.Canonicalized$Direction_Canonicalized <- str_replace_all(EffectSize.Canonicalized$Direction, "[+]", "*")
EffectSize.Canonicalized$Direction_Canonicalized <- str_replace_all(EffectSize.Canonicalized$Direction_Canonicalized, "[-]", "+")
EffectSize.Canonicalized$Direction_Canonicalized <- str_replace_all(EffectSize.Canonicalized$Direction_Canonicalized, "[*]", "-")

EffectSize.Keep <- subset(EffectSize, as.character(Allele2) == as.character(REF))
EffectSize.Keep$RSID_REF <- EffectSize.Keep$Allele2; EffectSize.Keep$RSID_ALT <- EffectSize.Keep$Allele1; EffectSize.Keep$Fre_ALT <- EffectSize.Keep$Freq1
EffectSize.Keep$Zscore_Canonicalized <- EffectSize.Keep$Zscore
EffectSize.Keep$Direction_Canonicalized <- EffectSize.Keep$Direction

Meta.eGFR.GWAS.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized <- rbind(EffectSize.Canonicalized, EffectSize.Keep)

### Save
write.table(Meta.eGFR.GWAS.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized, "Meta.eGFR.GWAS.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.txt",col.names=T,row.names=F,quote=F,sep="\t")
