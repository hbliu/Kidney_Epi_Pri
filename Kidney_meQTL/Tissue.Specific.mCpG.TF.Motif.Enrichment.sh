###########################################################################
#####    TF motif Enrichment Analysis of Tissue-specific mQTL         #####
#####    Kidney mQTL (n=443, this study)                              #####
#####    Whole blood mQTL (n=473, Sheng et al. PNAS 2020)             #####
#####    Skeletal muscle mQTL (n=265, Taylor et al. PNAS 2019)        #####
#####    Software: HOMER (v4.10.3)                                    #####
#####    Input: Tissue-specific mQTL CpGs                             #####
#####    Background: All CpGs tested in Multiple Tissue mQTL Analysis #####
###########################################################################

##### Test: Kidney specific mQTL CpGs, Background: All CpGs tested in Multiple Tissue mQTL Analysis
findMotifsGenome.pl Kidney.Specific.mQTL.CpG.bed hg19 Kidney.mQTL.CpGs.Background_CpG.MotifOutput -bg Background.CpG.bed

##### Test: Blood specific mQTL CpGs, Background: genome wide
findMotifsGenome.pl Blood.Specific.mQTL.CpG.bed hg19 Blood.mQTL.CpGs.Background_CpG.MotifOutput -bg Background.CpG.bed

##### Test: SkeMus specific mQTL CpGs, Background: All CpGs tested in Multiple Tissue mQTL Analysis
findMotifsGenome.pl SkeMus.Specific.mQTL.CpG.bed hg19 SkeMus.mQTL.CpGs.Background_CpG.MotifOutput -bg Background.CpG.bed

