###########################################################################
#####    Single cell GWAS enrichment by gchromVAR and snATAC-seq      #####
#####    GWAS: UKBB GWAS fine mapping (overlap mQTL SNPs)             #####
#####    Software: gchromVAR (v0.3.2)                                 #####
###########################################################################

###########################################################################
#### Step1: Download 94 fine mapped UKBB GWAS traits (UKBB_94traits_release1.bed.gz) from https://www.finucanelab.org/data

###########################################################################
#### Step2: Extract fine mapped SNPs included in the 95% credible sets by SusieR 
zcat UKBB_94traits_release1.bed.gz |\
awk '{if($11 == "SUSIE" && $19 != -1) print $5"\t"$1"\t"$2"\t"$3"\t"$13"\t"$18"\t"$12}' > UKBB_94traits_release1.txt

###########################################################################
#### Step3: Extract fine mapped SNPs identified as significant mQTL
awk 'NR==FNR{a[$0]++;next}a[$1]' Kidney.mQTL_S443_n35.1M.Sig.mQTL.SNP.txt UKBB_94traits_release1.txt | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'> UKBB_94traits_release1.mQTL.txt

###########################################################################
#### Step4: Separate file according to the IDs in one column
awk '{print > $6}' UKBB_94traits_release1.mQTL.txt

###########################################################################
#### Step5: Make bed format file for each trait
GWASList=(LipoA Smoking_CPD Glaucoma_Combined Inguinal_Hernia Depression_GP Hypothyroidism Fibroblastic_Disorders BrC Loneliness Cholelithiasis AFib CAD T2D_BMI Asthma T2D Worry_Too_Long AID_Combined Suffer_from_Nerves Insomnia MCP Risk_Taking Mood_Swings Tense MCHC Baso Glucose VitD Testosterone LOY Sensitivity Smoking_Ever_Never Guilty_Feelings Balding_Type4 Age_at_Menopause DVT FedUp_Feelings Miserableness Nervous_Feelings Irritability TC Morning_Person Worrier Neuroticism Urea ApoB TBil Ca ALT ApoA CRP Alb Age_at_Menarche SBP DBP TG AG WHRadjBMI AST UA TP LDLC GGT MAP SHBG FEV1FVC BFP PP Neutro Mono HDLC MCH HbA1c BMI Lym WBC MCV eGFRcys eGFR Hb Ht ALP Eosino RBC eBMD Plt BW IGF1 Height)
for GWAS in ${GWASList[@]};
do
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' ${GWAS} > ./bed/${GWAS}_PP001.bed
done


###########################################################################
#### Step6: Perform single cell GWAS enrichment analysis using gchromVAR
R
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
set.seed(123)
register(MulticoreParam(2))


### Building a RangedSummarizedExperiment
countsFile <- "Kidney.snATAC.Counts.tsv.gz"
peaksFile<- "Kidney.snATAC.Peaks.bed.gz"

peaksdf <- data.frame(fread(paste0("zcat < ", peaksFile)))
peaks <- makeGRangesFromDataFrame(peaksdf, seqnames = "V1", start.field = "V2", end.field = "V3")
counts <-  data.matrix(fread(paste0("zcat < ", countsFile)))
counts <- counts + 0.001

counts.filtered <- counts[which(Matrix::rowSums(counts)>0),]
peaksdf.filtered <- peaksdf[which(Matrix::rowSums(counts)>0),]
peaks.filtered <- makeGRangesFromDataFrame(peaksdf.filtered, seqnames = "V1", start.field = "V2", end.field = "V3")

SE <- SummarizedExperiment(assays = list(counts = counts.filtered),
                               rowData = peaks.filtered, 
                               colData = DataFrame(names = colnames(counts.filtered)))

SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg19)

### Importing GWAS summary statistics
files <- list.files("UKBB_88traits_bed", full.names = TRUE, pattern = "*.bed$")
head(read.table(files[1]))

ukbb <- importBedScore(rowRanges(SE), files, colidx = 5)
str(ukbb)

### Computing weighted deviations
ukbb_wDEV <- computeWeightedDeviations(SE, ukbb)
zdf <- reshape2::melt(t(assays(ukbb_wDEV)[["z"]]))
zdf[,2] <- gsub("_PP001", "", zdf[,2])
colnames(zdf) <- c("ct", "tr", "Zscore")
head(zdf)
write.table(zdf,"UKBB.GWAS.Kidney.snATAC.zdf.tsv", col.names=T, row.names = F, quote=F,sep=",")




###########################################################################
#### Step7: Heatmap and Kruskal test of single cell GWAS enrichment score by gchromVAR
R
library(dplyr)
library(tidyr)
library(ggplot2)
library(matrixStats)
library(diffloop)
library(viridis)

### Single-cell GWAS Enrichment Zscore (mQTL)
UKBB2KidATAC.zdf <- read.table("UKBB.GWAS.Kidney.snATAC.zdf.tsv", header=T, sep=",")
UKBB2KidATAC.zdf$ct <- gsub("\\.", "-", UKBB2KidATAC.zdf$ct)

### Median for each pair
UKBB2KidATAC.zdf$ZscoreLn <- log(abs(UKBB2KidATAC.zdf$Zscore)+1) * sign(UKBB2KidATAC.zdf$Zscore)
UKBB2KidATAC.zdf.Summary <- UKBB2KidATAC.zdf %>% dplyr::group_by(tr,cluster)%>% dplyr::summarise(Mean=mean(ZscoreLn), Max=max(ZscoreLn), Min=min(ZscoreLn), Median=median(ZscoreLn), Std=sd(ZscoreLn))
library(reshape2)
UKBB2KidATAC.zdf.Summary.Mean <- acast(UKBB2KidATAC.zdf.Summary, cluster~tr, value.var="Mean")
UKBB2KidATAC.zdf.Summary.Mean[which(UKBB2KidATAC.zdf.Summary.Mean > 0.1)] <- 0.1
UKBB2KidATAC.zdf.Summary.Mean[which(UKBB2KidATAC.zdf.Summary.Mean < -0.1)] <- -0.1
UKBB2KidATAC.zdf.Summary.Mean <- UKBB2KidATAC.zdf.Summary.Mean[c("PT-S1","PT-S2","PT-S3","LOH","DCT","PC","IC","Podo","Endo","Stroma","Immune","Lymph"),c("eGFR","eGFRcys","MCHC","GGT","IGF1","TC","RBC","Smoking_Ever_Never","ApoA","HDLC","TG","BW","Hb","T2D_BMI","T2D","TBil","Neutro","Mono","Plt","eBMD","ALP","Eosino","SHBG","Age_at_Menopause","AG","FEV1FVC","ALT","CRP","HbA1c","AST","Worrier","Alb","WHRadjBMI","LipoA","Glucose","BrC","Urea","Age_at_Menarche","Sensitivity","WBC","Lym","TP","BMI","Height","MCH","Baso","Ca","Ht","BFP","Worry_Too_Long","Miserableness","Depression_GP","Morning_Person","Nervous_Feelings","MAP","Loneliness","SBP","DBP","Neuroticism","Hypothyroidism","Tense","LDLC","Testosterone","ApoB","DVT","MCP","AID_Combined","Risk_Taking","Smoking_CPD","CAD","PP","MCV","UA","LOY","VitD","Insomnia","Irritability","Mood_Swings","FedUp_Feelings","Guilty_Feelings","Inguinal_Hernia","Suffer_from_Nerves","Fibroblastic_Disorders","Balding_Type4","Asthma","Glaucoma_Combined","Cholelithiasis","AFib")]
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 256)
library("gplots")
x<- as.matrix(UKBB2KidATAC.zdf.Summary.Mean);
pdf("Single Cell GWAS enrichment Mean heatmap.pdf", width = 20, height = 7)
heatmap <- heatmap.2(x, scale = "none", Rowv = FALSE, Colv = FALSE, col=my_palette, trace = "none", density.info = "none", margins = c(20, 20))
dev.off()
write.table(heatmap$carpet,"Single Cell GWAS enrichment Mean heatmap.txt",col.names=T,row.names=T,quote=F,sep="\t")


####### Statistic for significance of difference among cell types for each trait
UKBB2KidATAC.Enrich.Statistic <- UKBB2KidATAC.zdf.Summary.Mean[c("PT-S1","PT-S2","PT-S3","LOH","DCT","PC","IC","Podo","Endo","Stroma","Immune","Lymph"),c("eGFR","eGFRcys","MCHC","GGT","IGF1","TC","RBC","Smoking_Ever_Never","ApoA","HDLC","TG","BW","Hb","T2D_BMI","T2D","TBil","Neutro","Mono","Plt","eBMD","ALP","Eosino","SHBG","Age_at_Menopause","AG","FEV1FVC","ALT","CRP","HbA1c","AST","Worrier","Alb","WHRadjBMI","LipoA","Glucose","BrC","Urea","Age_at_Menarche","Sensitivity","WBC","Lym","TP","BMI","Height","MCH","Baso","Ca","Ht","BFP","Worry_Too_Long","Miserableness","Depression_GP","Morning_Person","Nervous_Feelings","MAP","Loneliness","SBP","DBP","Neuroticism","Hypothyroidism","Tense","LDLC","Testosterone","ApoB","DVT","MCP","AID_Combined","Risk_Taking","Smoking_CPD","CAD","PP","MCV","UA","LOY","VitD","Insomnia","Irritability","Mood_Swings","FedUp_Feelings","Guilty_Feelings","Inguinal_Hernia","Suffer_from_Nerves","Fibroblastic_Disorders","Balding_Type4","Asthma","Glaucoma_Combined","Cholelithiasis","AFib")]
UKBB2KidATAC.Enrich.Statistic <- t(UKBB2KidATAC.Enrich.Statistic)
UKBB2KidATAC.Enrich.Statistic <- as.data.frame(UKBB2KidATAC.Enrich.Statistic)
UKBB2KidATAC.Enrich.Statistic$kruskal.test.p = 1
for (i in 1:88) {
UKBB2KidATAC.zdf.tmp <- subset(UKBB2KidATAC.zdf, tr == rownames(UKBB2KidATAC.Enrich.Statistic)[i])
UKBB2KidATAC.Enrich.Statistic[i,]$kruskal.test.p <- kruskal.test(ZscoreLn ~ cluster, data = UKBB2KidATAC.zdf.tmp)$p.value
}
write.table(UKBB2KidATAC.Enrich.Statistic,"UKBB2KidATAC.Enrich.Statistic.txt",col.names=T,row.names=T,quote=F,sep="\t")


###### Selected traits for main figure
UKBB2KidATAC.zdf.Summary.Mean.MainFigure <- UKBB2KidATAC.zdf.Summary.Mean[c("PT-S1","PT-S2","PT-S3","LOH","DCT","PC","IC","Podo","Endo","Stroma","Immune","Lymph"),c("eGFR","eGFRcys","IGF1","MCHC","GGT","DBP","SBP","Ca","Urea")]
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 256)
library("gplots")
pdf("Single Cell GWAS enrichment Mean heatmap MainFigure.pdf", width = 7, height = 7)
x<- as.matrix(UKBB2KidATAC.zdf.Summary.Mean.MainFigure); 
heatmap <- heatmap.2(x, scale = "none", Rowv = FALSE, Colv = FALSE, col=my_palette, trace = "none", density.info = "none", margins = c(7, 7))
dev.off()
write.table(heatmap$carpet,"Single Cell GWAS enrichment Mean heatmap MainFigure.txt",col.names=T,row.names=T,quote=F,sep="\t")








