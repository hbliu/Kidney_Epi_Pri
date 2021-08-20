###########################################################################
#####    Multiple-tissue eQTL analysis                                #####
#####    eQTLs from Kidney and 48 GTEx (v7) tissues                   #####
#####    Software: METASOFT (v2.0.1)                                  #####
#####    Input: mQTLs identified in at least one dataset              #####
###########################################################################

###########################################################################
######### Step1: Download eQTLs from GTEx (v7)
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_all_associations.tar.gz


###########################################################################
######### Step2: Combine Kidney eQTL and GTEx eQTLs
SampleList=(Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra Breast_Mammary_Tissue Cells_EBV-transformed_lymphocytes Cells_Transformed_fibroblasts Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Liver Lung Minor_Salivary_Gland Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood)
for SamName in ${SampleList[@]};
do
awk 'NR==FNR{a[$0]++;next}a[$1]' Kidney.Sig.eQTL.ID.txt SamName.allpairs.Simple.txt > SamName.allpairs.Simple.KidID.txt
R
library(dplyr)
library(tidyr)
Kidney <- read.table("Kidney.Sig.eQTL.Uniq.txt", header=T, sep="\t"); colnames(Kidney) <- c("Pair","Kidney_Beta","Kidney_SE")

SamName <- read.table("SamName.allpairs.Simple.KidID.txt", header=F, sep="\t"); colnames(SamName) <- c("Pair","SamName_Beta","SamName_SE")
SamName <- SamName[!duplicated(SamName$Pair), ]
SamName <- left_join(Kidney, SamName, by = c("Pair" = "Pair"))
write.table(SamName[,4:5], "SamName.2Kidney.txt", quote=FALSE, row.names=FALSE, col.names=F, sep="\t")
done

####### Combine all
awk 'NR>1{print $0}' Kidney.Sig.eQTL.Uniq.txt > Kidney.Sig.eQTL.Uniq.Noheader.txt
paste Kidney.Sig.eQTL.Uniq.Noheader.txt Adipose_Subcutaneous.2Kidney.txt Adipose_Visceral_Omentum.2Kidney.txt Adrenal_Gland.2Kidney.txt Artery_Aorta.2Kidney.txt Artery_Coronary.2Kidney.txt Artery_Tibial.2Kidney.txt Brain_Amygdala.2Kidney.txt Brain_Anterior_cingulate_cortex_BA24.2Kidney.txt Brain_Caudate_basal_ganglia.2Kidney.txt Brain_Cerebellar_Hemisphere.2Kidney.txt Brain_Cerebellum.2Kidney.txt Brain_Cortex.2Kidney.txt Brain_Frontal_Cortex_BA9.2Kidney.txt Brain_Hippocampus.2Kidney.txt Brain_Hypothalamus.2Kidney.txt Brain_Nucleus_accumbens_basal_ganglia.2Kidney.txt Brain_Putamen_basal_ganglia.2Kidney.txt Brain_Spinal_cord_cervical_c_1.2Kidney.txt Brain_Substantia_nigra.2Kidney.txt Breast_Mammary_Tissue.2Kidney.txt Cells_EBV_transformed_lymphocytes.2Kidney.txt Cells_Transformed_fibroblasts.2Kidney.txt Colon_Sigmoid.2Kidney.txt Colon_Transverse.2Kidney.txt Esophagus_Gastroesophageal_Junction.2Kidney.txt Esophagus_Mucosa.2Kidney.txt Esophagus_Muscularis.2Kidney.txt Heart_Atrial_Appendage.2Kidney.txt Heart_Left_Ventricle.2Kidney.txt Liver.2Kidney.txt Lung.2Kidney.txt Minor_Salivary_Gland.2Kidney.txt Muscle_Skeletal.2Kidney.txt Nerve_Tibial.2Kidney.txt Ovary.2Kidney.txt Pancreas.2Kidney.txt Pituitary.2Kidney.txt Prostate.2Kidney.txt Skin_Not_Sun_Exposed_Suprapubic.2Kidney.txt Skin_Sun_Exposed_Lower_leg.2Kidney.txt Small_Intestine_Terminal_Ileum.2Kidney.txt Spleen.2Kidney.txt Stomach.2Kidney.txt Testis.2Kidney.txt Thyroid.2Kidney.txt Uterus.2Kidney.txt Vagina.2Kidney.txt Whole_Blood.2Kidney.txt > Kidney.GTEx48.meta



###########################################################################
######### Step3: Multiple-tissue eQTL analysis using Metasoft
#### First run to get the lambda_hetero and lambda_mean
java -jar Metasoft.jar \
-input Kidney.GTEx48.meta \
-log Kidney.GTEx48.NoHeader.meta.metasoft.log \
-pvalue_table HanEskinPvalueTable.txt \
-output Kidney.GTEx48.meta.metasoft.txt

# Newly calculated inflation factor lambda for   mean effect part: 367.445612
# Newly calculated inflation factor lambda for heterogeneity part: 76.601144

#### Second run to calculate m value for each eQTL in each dataset 
java -jar /home/hongbol/Software/Metasoft/Metasoft/Metasoft.jar \
-input Kidney.GTEx48.meta \
-log Kidney.GTEx48.NoHeader.meta.metasoft.log \
-lambda_mean 367.445612 -lambda_hetero 76.601144 \
-mvalue -mvalue_method mcmc -mvalue_p_thres 1 \
-pvalue_table HanEskinPvalueTable.txt \
-output Kidney.GTEx48.meta.metasoft.txt



###########################################################################
######### Step4: Identification Kidney-specific eQTL
R
library(dplyr)
library(tidyr)
Kidney.GTEx48.meta.Mvalue <- read.table("Kidney.GTEx48.metasoft.Mvalue.txt", header=T, sep="\t")
Kidney.Sig.eQTL <- read.table("Meta.Kidney.eQTL.Sig.Gene.q0.05.txt", header=T, sep="\t")

### Calculate number and fraction of tissues with m value > 0.9
Kidney.GTEx48.meta.Mvalue <- cbind(Kidney.GTEx48.meta.Mvalue, Mvalue.0.9.Num = rowSums(Kidney.GTEx48.meta.Mvalue[,3:51] > 0.9, na.rm = TRUE))
Kidney.GTEx48.meta.Mvalue$Mvalue.0.9.Frac <- Kidney.GTEx48.meta.Mvalue$Mvalue.0.9.Num / Kidney.GTEx48.meta.Mvalue$STUDY

### Identify kidney-specific eQTL
Kidney.GTEx48.meta.Mvalue.Kid.Specific.S49.N5 <- subset(Kidney.GTEx48.meta.Mvalue.Kid.Specific, Kidney_Mvalue > 0.9 & STUDY >= 49 & Mvalue.0.9.Num < 5); Kidney.GTEx48.meta.Mvalue.Kid.Specific.S49.N5$Specificity <- "Kidney_Specific"
write.table(Kidney.GTEx48.meta.Mvalue.Kid.Specific.S49.N5, "Kidney.GTEx48.meta.Mvalue.Kid.Specific.S49.N5.txt", quote=FALSE, row.names=FALSE, col.names=T, sep="\t")

### Add eQTL information
Kidney.GTEx48.meta.Mvalue.Kid.Specific.S49.N5.infor <- left_join(Kidney.GTEx48.meta.Mvalue.Kid.Specific.S49.N5, Kidney.Sig.eQTL, by = c("GeneID" = "GeneID", "SNP" = "SNP_POS"))
write.table(Kidney.GTEx48.meta.Mvalue.Kid.Specific.S49.N5.infor, "Kidney.GTEx48.meta.Mvalue.Kid.Specific.S49.N5.infor.txt", quote=FALSE, row.names=FALSE, col.names=T, sep="\t")



