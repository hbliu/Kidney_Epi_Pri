###########################################################################
#####    Overlap between Kidney cell type-specific DARs and eGFR GWAS SNPs#
#####    Software: bedtools (v2.29.2)                                 #####
###########################################################################
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b PT_S1.cluster.1.DAR.7.bed -u > eGFR_SNP.to.PT_S1.cluster.1.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b PT_S2.cluster.4.DAR.7.bed -u > eGFR_SNP.to.PT_S2.cluster.4.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b PT_S3.cluster.3.DAR.7.bed -u > eGFR_SNP.to.PT_S3.cluster.3.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b LOH.cluster.2.DAR.7.bed -u > eGFR_SNP.to.LOH.cluster.2.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b DCT.cluster.8.DAR.7.bed -u > eGFR_SNP.to.DCT.cluster.8.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b PC.cluster.6.DAR.7.bed -u > eGFR_SNP.to.PC.cluster.6.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b IC.cluster.10.DAR.7.bed -u > eGFR_SNP.to.IC.cluster.10.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b Endo.cluster.5.DAR.7.bed -u > eGFR_SNP.to.Endo.cluster.5.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b Podo.cluster.11.DAR.7.bed -u > eGFR_SNP.to.Podo.cluster.11.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b Stroma.cluster.9.DAR.7.bed -u > eGFR_SNP.to.Stroma.cluster.9.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b Immune.cluster.12.DAR.7.bed -u > eGFR_SNP.to.Immune.cluster.12.DAR.7.bed
bedtools intersect -a Metal.CKDGen.UKB.MVP.PAGE.SUMMIT.Canonicalized.SNP.bed -b Lymph.cluster.13.DAR.7.bed -u > eGFR_SNP.to.Lymph.cluster.13.DAR.7.bed
