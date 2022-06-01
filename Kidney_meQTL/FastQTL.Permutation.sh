###########################################################################
#####    FastQTL permutation                                          #####
#####    Sample size: 443 human kidney samples                        #####
#####    Software: fastQTL (v2.0)                                     #####
###########################################################################
fastQTL \
--vcf ../Genotype.dosage.vcf.gz \
--bed ../Methylation.matrix.bed.gz \
--cov ../Covariates.txt.gz \
--permute 1000 \
--out FastQTL.permutations.txt


# In the file of FastQTL.permutations.txt, the 10 columns correspond to:
# 1. pid, ID of the tested molecular phenotype (in this particular case, the gene ID)
# 2. nvar, Number of variants tested in cis for this phenotype
# 3. shape1, MLE of the shape1 parameter of the Beta distribution
# 4. shape2, MLE of the shape2 parameter of the Beta distribution
# 5. dummy, Dummy [To be described later]
# 6. sid, ID of the best variant found for this molecular phenotypes (i.e. with the smallest p-value)
# 7. dist, Distance between the molecular phenotype - variant pair
# 8. npval, The nominal p-value of association that quantifies how significant from 0, the regression coefficient is
# 9. slope, The slope associated with the nominal p-value of association [only in version > v2-184]
# 10. ppval, A first permutation p-value directly obtained from the permutations with the direct method. This is basically a corrected version of the nominal p-value that accounts for the fact that multiple variants are tested per molecular phenotype.
# 11. bpval, A second permutation p-value obtained via beta approximation. We advice to use this one in any downstream analysis.
