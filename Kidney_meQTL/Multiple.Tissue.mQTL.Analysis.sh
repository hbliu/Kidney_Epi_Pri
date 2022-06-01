###########################################################################
#####    Multiple-tissue mQTL analysis                                #####
#####    Three mQTL datasets based on EPIC-array methylation          #####
#####    Kidney mQTL (n=443, this study)                              #####
#####    Whole blood mQTL (n=473, Sheng et al. PNAS 2020)             #####
#####    Skeletal muscle mQTL (n=265, Taylor et al. PNAS 2019)        #####
#####    Software: METASOFT (v2.0.1)                                  #####
#####    Input: mQTLs identified in at least one dataset              #####
###########################################################################

#### First run to get the lambda_hetero and lambda_mean
java -jar Metasoft.jar \
-input Kidney_Blood_SkeletalMus.for.metasoft.meta \
-log Kidney_Blood_SkeletalMus.meta.metasoft.log \
-pvalue_table HanEskinPvalueTable.txt \
-output Kidney_Blood_SkeletalMus.meta.metasoft.txt

# Newly calculated inflation factor lambda for   mean effect part: 62.232925
# Newly calculated inflation factor lambda for heterogeneity part: 43.903213

#### Second run to calculate m value for each mQTL in each dataset 
java -jar Metasoft.jar \
-input Kidney_Blood_SkeletalMus.for.metasoft.meta \
-log Kidney_Blood_SkeletalMus.meta.metasoft.log \
-lambda_mean 62.232925 -lambda_hetero 43.903213 \
-mvalue -mvalue_method mcmc -mvalue_p_thres 1 \
-pvalue_table HanEskinPvalueTable.txt \
-output Kidney_Blood_SkeletalMus.meta.metasoft.txt

