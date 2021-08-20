###########################################################################
#####    Heritability mediated by assayed gene expression             #####
#####    GWAS datasets: 43 traits collected with sample size > 200000 #####
#####    mQTL datasets: Kidney, Whole blood and Skeletal muscle       #####
#####    eQTL datasets: Kidney and 48 GTEx tissues                    #####
#####    Software: MESC (commit: 3312dda)                             #####
###########################################################################
##### Estimate overall h2med for each QTL dataset and each GWAS trait
##### Step1: Download GTEx v8 expression scores
# https://github.com/douglasyao/mesc/wiki/Download-expression-scores

##### Step2: For other QTL datasets, the expression/methylation scores was calculated by following command, in which QTL.Summary.Statistics.txt is summary statistics formatted as "GENE	SNP_Loc	CHR	SNP_COORD	N	Z"
./Software/mesc/run_mesc.py --compute-expscore-sumstat --eqtl-sumstat QTL.Summary.Statistics.txt --out QTL.mesc

##### Step3: Estimate overall h2med for each QTL dataset and each GWAS trait
./Software/mesc/run_mesc.py --h2med GWAS.sumstats.gz --exp-chr QTL.mesc.@ --out GWAS.QTL.mesc.h2med

