###########################################################################
#####    Summary-data-based Mendelian Randomization (SMR)             #####
#####    Heterogeneity in dependent instruments analysis (HEIDI)      #####
#####    Software: SMR (v1.03)                                        #####
###########################################################################

##### SMR and HEIDI test for eGFR GWAS and kidney eQTL
##### Target SNP~eGene was those identified in coloc analysis between eGFR GWAS and kidney eQTL
for SNP~eGene in ${SNP~eGene.List[@]};
./Software/SMR/smr_Linux --bfile ./PLINK2_1000GP_Phase3/EUR_phase3_MAF05 \
--gwas-summary eGFR.GWAS.ma \
--beqtl-summary Kidney.eQTL \
--diff-freq-prop 0.1 \
--peqtl-heidi 0.05 \
--extract-target-snp-probe ${SNP~eGene} \
--out eGFR.GWAS.Kidney.eQTL.${SNP~eGene}.smr
done


##### SMR and HEIDI test for eGFR GWAS and kidney mQTL
##### Target SNP~mCpG was those identified in coloc analysis between eGFR GWAS and kidney mQTL
for SNP~mCpG in ${SNP~mCpG.List[@]};
do
./Software/SMR/smr_Linux --bfile ./PLINK2_1000GP_Phase3/EUR_phase3_MAF05 \
--gwas-summary eGFR.GWAS.ma \
--beqtl-summary Kidney.mQTL \
--diff-freq-prop 0.1 \
--peqtl-heidi 0.05 \
--extract-target-snp-probe ${SNP~mCpG} \
--out eGFR.GWAS.Kidney.mQTL.${SNP~mCpG}.smr
done


##### SMR and HEIDI test for kidney mQTL and kidney eQTL
##### Test all pairs of eGene~mCpG identified in coloc analysis between kidney mQTL and kidney eQTL
./Software/SMR/smr_Linux --bfile ./PLINK2_1000GP_Phase3/EUR_phase3_MAF05 \
--peqtl-smr 1e-4 \
--peqtl-heidi 0.05 \
--beqtl-summary Kidney.mQTL \
--beqtl-summary Kidney.eQTL \
--out Kidney.mQTL.eQTL.smr
