###########################################################################
#####    GWAS traits collection and formatting for heritability analysis ##
#####    Only traits with sample size > 200,000 were collected        #####
#####    Software: LDSC (v1.0.1)                                      #####
###########################################################################

########### Make summary data for meta-analysis, MVP and UKB eGFR GWAS datasets
./Software/ldsc/munge_sumstats.py --sumstats Meta_eGFR_GWAS_summary.ma --snp MarkerName --out Meta_eGFR_TransAncestry.sumstats
./Software/ldsc/munge_sumstats.py --sumstats MVP_eGFR_GWAS_summary.ma --snp SNP --out MVP_eGFR_TransAncestry.sumstats
./Software/ldsc/munge_sumstats.py --sumstats UKB_eGFR_GWAS_summary.ma --snp SNP --out UKB_eGFR_TransAncestry.sumstats

########### Make summary data for GWAS datasets obtained from CKDGen
./Software/ldsc/munge_sumstats.py --sumstats ../Wuttke/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.txt --snp RSID --out CKDgen_Wuttke_eGFR_overall_ALL.sumstats
./Software/ldsc/munge_sumstats.py --sumstats ../Wuttke/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt --snp RSID --out CKDgen_Wuttke_eGFR_overall_EA.sumstats
./Software/ldsc/munge_sumstats.py --sumstats ../Teumer/formatted_20170711-UACR_overall-ALL-nstud_27-sumMac_400.tbl.rsid --snp RSID --out CKDgen_Teumer_UACR_overall_All.sumstats
./Software/ldsc/munge_sumstats.py --sumstats ../Teumer/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid --snp RSID --out CKDgen_Teumer_UACR_overall_EA.sumstats
./Software/ldsc/munge_sumstats.py --sumstats ../Teumer/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid --snp RSID --out CKDgen_Teumer_MA_overall_ALL.sumstats
./Software/ldsc/munge_sumstats.py --sumstats ../Tin/urate_chr1_22_LQ_IQ06_mac10_all_741_nstud37_summac400_rsid.New.txt --snp RSID --out CKDgen_Tin_Urate_All.sumstats
./Software/ldsc/munge_sumstats.py --sumstats ../Tin/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.New.txt --snp RSID --out CKDgen_Tin_Urate_EA.sumstats

########### Download public available GWAS sumstats files from Alkes group
### https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.blood_EOSINOPHIL_COUNT.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.blood_PLATELET_COUNT.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.blood_RBC_DISTRIB_WIDTH.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.blood_RED_COUNT.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.blood_WHITE_COUNT.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.bmd_HEEL_TSCOREz.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.body_BALDING1.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.body_BMIz.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.body_HEIGHTz.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.body_WHRadjBMIz.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.bp_SYSTOLICadjMEDz.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.cov_EDU_YEARS.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.cov_SMOKING_STATUS.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.disease_HI_CHOL_SELF_REP.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.disease_T2D.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.lung_FEV1FVCzSMOKE.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.lung_FVCzSMOKE.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.mental_NEUROTICISM.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.other_MORNINGPERSON.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.pigment_HAIR.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.pigment_SKIN.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.pigment_SUNBURN.sumstats.gz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats/UKB_460K.pigment_TANNING.sumstats.gz

########### Download public available GWAS sumstats files from Neale lab
## https://nealelab.github.io/UKBB_ldsc/downloads.html
wget https://www.dropbox.com/s/erscb0kgs9c10v3/30510_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=0 -O UKB_Nealelab.Creatinine_enzymatic_urine.sumstats.bgz
wget https://www.dropbox.com/s/rlyuzamxja4j6lg/30700_irnt.imputed_v3.ldsc.both_sexes.tsv.gz?dl=0 -O UKB_Nealelab.Creatinine_umolperL.sumstats.gz
wget https://www.dropbox.com/s/svaojchyg01aikw/2443.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=0 -O UKB_Nealelab.Diabetes.sumstats.bgz
wget https://www.dropbox.com/s/qee9ruazsxa8hj5/4079_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=0 -O UKB_Nealelab.Diastolic_blood_pressure.sumstats.bgz
wget https://www.dropbox.com/s/lsjxn08zxneuxcc/4080_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=0 -O UKB_Nealelab.Systolic_blood_pressure.sumstats.bgz
wget https://www.dropbox.com/s/jtmcu2qq8gx4wpq/6150_4.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=0 -O UKB_Nealelab.High_blood_pressure.sumstats.bgz
wget https://www.dropbox.com/s/f1dmxqpjjo89z3l/30520_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=0 -O UKB_Nealelab.Potassium_urine.sumstats.bgz
wget https://www.dropbox.com/s/b4947ktkiza104h/30530_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=0 -O UKB_Nealelab.Sodium_urine.sumstats.bgz


