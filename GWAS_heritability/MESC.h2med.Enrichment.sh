###########################################################################
#####    Enrichment of heritability mediated by assayed gene expression ###
#####    GWAS datasets: 43 traits collected with sample size > 200000 #####
#####    mQTL datasets: Kidney, Whole blood                           #####
#####    CpG set: Kidney HMM states, Enhancers, cell type-specific DARs ###
#####    Software: MESC (commit: 3312dda)                             #####
###########################################################################

###########################################################################
##### Estimation of CpG set (Kidney HMM state) methylation scores and h2med for Kidney mQTL
##### This procedure was also applied to Kidney mQTL and enhancer regions of 127 Roadmap tissues 
##### This procedure was also applied to Blood mQTL and enhancer regions of 127 Roadmap tissues
##### Step1: Make CpG Set annotation file based on kidney HMM state
bedtools intersect -a mQTL.CpG.bed -b ./HMM/1_TssA -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_1_TssA\t"$0}' > CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/2_TssAFlnk -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_2_TssAFlnk\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/4_Tx -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_4_Tx\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/5_TxWk -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_5_TxWk\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/6_EnhG -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_6_EnhG\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/7_Enh -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_7_Enh\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/8_ZNF_or_Rpts -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_8_ZNF_or_Rpts\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/9_Het -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_9_Het\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/10_TssBiv -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_10_TssBiv\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/13_ReprPC -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_13_ReprPC\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/14_ReprPCWk -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_14_ReprPCWk\t"$0}' >> CpG_Set.HMM.txt
bedtools intersect -a mQTL.CpG.bed -b ./HMM/15_Quies -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "HMM_15_Quies\t"$0}' >> CpG_Set.HMM.txt

##### Step2: Estimation of CpG set (Kidney HMM state) methylation scores for Kidney mQTL
./Software/mesc/run_mesc.py --compute-expscore-sumstat --eqtl-sumstat Kidney.mQTL.Summary.Statistics.txt --gene-sets CpG_Set.HMM.txt --out Kidney.mQTL.mesc.CpG_Set.HMM

##### Step3: Estimate h2med enrichment of Kidney mQTL for each GWAS trait and each HMM state
./Software/mesc/run_mesc.py --h2med GWAS.sumstats.gz --exp-chr Kidney.mQTL.mesc.CpG_Set.HMM.@ --out Kidney.mQTL.mesc.CpG_Set.HMM.h2med


###########################################################################
##### Estimation of CpG set (Kidney and Roadmap tissue enhancers) methylation scores and h2med for Kidney mQTL and Blood mQTL, respectively
##### Step1: Make CpG Set annotation file based on Kidney and Roadmap tissue enhancers
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E000_Adult_Kidney_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E000_Adult_Kidney_HMM_Enhancer\t"$0}' > CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E001.ESC.I3_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E001.ESC.I3_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E002.ESC.WA7_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E002.ESC.WA7_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E003.ESC.H1_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E003.ESC.H1_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E004.ESDR.H1.BMP4.MESO_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E004.ESDR.H1.BMP4.MESO_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E005.ESDR.H1.BMP4.TROP_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E005.ESDR.H1.BMP4.TROP_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E006.ESDR.H1.MSC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E006.ESDR.H1.MSC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E007.ESDR.H1.NEUR.PROG_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E007.ESDR.H1.NEUR.PROG_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E008.ESC.H9_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E008.ESC.H9_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E009.ESDR.H9.NEUR.PROG_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E009.ESDR.H9.NEUR.PROG_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E010.ESDR.H9.NEUR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E010.ESDR.H9.NEUR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E011.ESDR.CD184.ENDO_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E011.ESDR.CD184.ENDO_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E012.ESDR.CD56.ECTO_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E012.ESDR.CD56.ECTO_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E013.ESDR.CD56.MESO_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E013.ESDR.CD56.MESO_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E014.ESC.HUES48_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E014.ESC.HUES48_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E015.ESC.HUES6_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E015.ESC.HUES6_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E016.ESC.HUES64_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E016.ESC.HUES64_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E017.LNG.IMR90_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E017.LNG.IMR90_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E018.IPSC.15b_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E018.IPSC.15b_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E019.IPSC.18_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E019.IPSC.18_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E020.IPSC.20B_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E020.IPSC.20B_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E021.IPSC.DF.6.9_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E021.IPSC.DF.6.9_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E022.IPSC.DF.19.11_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E022.IPSC.DF.19.11_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E023.FAT.MSC.DR.ADIP_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E023.FAT.MSC.DR.ADIP_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E024.ESC.4STAR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E024.ESC.4STAR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E025.FAT.ADIP.DR.MSC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E025.FAT.ADIP.DR.MSC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E026.STRM.MRW.MSC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E026.STRM.MRW.MSC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E027.BRST.MYO_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E027.BRST.MYO_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E028.BRST.HMEC.35_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E028.BRST.HMEC.35_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E029.BLD.CD14.PC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E029.BLD.CD14.PC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E030.BLD.CD15.PC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E030.BLD.CD15.PC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E031.BLD.CD19.CPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E031.BLD.CD19.CPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E032.BLD.CD19.PPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E032.BLD.CD19.PPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E033.BLD.CD3.CPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E033.BLD.CD3.CPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E034.BLD.CD3.PPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E034.BLD.CD3.PPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E035.BLD.CD34.PC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E035.BLD.CD34.PC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E036.BLD.CD34.CC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E036.BLD.CD34.CC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E037.BLD.CD4.MPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E037.BLD.CD4.MPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E038.BLD.CD4.NPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E038.BLD.CD4.NPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E039.BLD.CD4.CD25M.CD45RA.NPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E039.BLD.CD4.CD25M.CD45RA.NPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E040.BLD.CD4.CD25M.CD45RO.MPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E040.BLD.CD4.CD25M.CD45RO.MPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E041.BLD.CD4.CD25M.IL17M.PL.TPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E041.BLD.CD4.CD25M.IL17M.PL.TPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E042.BLD.CD4.CD25M.IL17P.PL.TPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E042.BLD.CD4.CD25M.IL17P.PL.TPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E043.BLD.CD4.CD25M.TPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E043.BLD.CD4.CD25M.TPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E044.BLD.CD4.CD25.CD127M.TREGPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E044.BLD.CD4.CD25.CD127M.TREGPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E045.BLD.CD4.CD25I.CD127.TMEMPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E045.BLD.CD4.CD25I.CD127.TMEMPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E046.BLD.CD56.PC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E046.BLD.CD56.PC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E047.BLD.CD8.NPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E047.BLD.CD8.NPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E048.BLD.CD8.MPC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E048.BLD.CD8.MPC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E049.STRM.CHON.MRW.DR.MSC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E049.STRM.CHON.MRW.DR.MSC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E050.BLD.MOB.CD34.PC.F_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E050.BLD.MOB.CD34.PC.F_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E051.BLD.MOB.CD34.PC.M_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E051.BLD.MOB.CD34.PC.M_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E052.MUS.SAT_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E052.MUS.SAT_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E053.BRN.CRTX.DR.NRSPHR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E053.BRN.CRTX.DR.NRSPHR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E054.BRN.GANGEM.DR.NRSPHR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E054.BRN.GANGEM.DR.NRSPHR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E055.SKIN.PEN.FRSK.FIB.01_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E055.SKIN.PEN.FRSK.FIB.01_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E056.SKIN.PEN.FRSK.FIB.02_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E056.SKIN.PEN.FRSK.FIB.02_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E057.SKIN.PEN.FRSK.KER.02_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E057.SKIN.PEN.FRSK.KER.02_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E058.SKIN.PEN.FRSK.KER.03_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E058.SKIN.PEN.FRSK.KER.03_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E059.SKIN.PEN.FRSK.MEL.01_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E059.SKIN.PEN.FRSK.MEL.01_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E061.SKIN.PEN.FRSK.MEL.03_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E061.SKIN.PEN.FRSK.MEL.03_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E062.BLD.PER.MONUC.PC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E062.BLD.PER.MONUC.PC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E063.FAT.ADIP.NUC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E063.FAT.ADIP.NUC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E065.VAS.AOR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E065.VAS.AOR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E066.LIV.ADLT_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E066.LIV.ADLT_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E067.BRN.ANG.GYR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E067.BRN.ANG.GYR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E068.BRN.ANT.CAUD_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E068.BRN.ANT.CAUD_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E069.BRN.CING.GYR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E069.BRN.CING.GYR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E070.BRN.GRM.MTRX_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E070.BRN.GRM.MTRX_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E071.BRN.HIPP.MID_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E071.BRN.HIPP.MID_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E072.BRN.INF.TMP_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E072.BRN.INF.TMP_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E073.BRN.DL.PRFRNTL.CRTX_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E073.BRN.DL.PRFRNTL.CRTX_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E074.BRN.SUB.NIG_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E074.BRN.SUB.NIG_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E075.GI.CLN.MUC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E075.GI.CLN.MUC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E076.GI.CLN.SM.MUS_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E076.GI.CLN.SM.MUS_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E077.GI.DUO.MUC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E077.GI.DUO.MUC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E078.GI.DUO.SM.MUS_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E078.GI.DUO.SM.MUS_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E079.GI.ESO_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E079.GI.ESO_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E080.ADRL.GLND.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E080.ADRL.GLND.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E081.BRN.FET.M_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E081.BRN.FET.M_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E082.BRN.FET.F_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E082.BRN.FET.F_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E083.HRT.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E083.HRT.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E084.GI.L.INT.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E084.GI.L.INT.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E085.GI.S.INT.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E085.GI.S.INT.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E086.KID.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E086.KID.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E087.PANC.ISLT_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E087.PANC.ISLT_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E088.LNG.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E088.LNG.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E089.MUS.TRNK.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E089.MUS.TRNK.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E090.MUS.LEG.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E090.MUS.LEG.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E091.PLCNT.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E091.PLCNT.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E092.GI.STMC.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E092.GI.STMC.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E093.THYM.FET_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E093.THYM.FET_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E094.GI.STMC.GAST_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E094.GI.STMC.GAST_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E095.HRT.VENT.L_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E095.HRT.VENT.L_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E096.LNG_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E096.LNG_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E097.OVRY_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E097.OVRY_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E098.PANC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E098.PANC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E099.PLCNT.AMN_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E099.PLCNT.AMN_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E100.MUS.PSOAS_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E100.MUS.PSOAS_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E101.GI.RECT.MUC.29_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E101.GI.RECT.MUC.29_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E102.GI.RECT.MUC.31_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E102.GI.RECT.MUC.31_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E103.GI.RECT.SM.MUS_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E103.GI.RECT.SM.MUS_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E104.HRT.ATR.R_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E104.HRT.ATR.R_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E105.HRT.VNT.R_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E105.HRT.VNT.R_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E106.GI.CLN.SIG_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E106.GI.CLN.SIG_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E107.MUS.SKLT.M_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E107.MUS.SKLT.M_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E108.MUS.SKLT.F_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E108.MUS.SKLT.F_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E109.GI.S.INT_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E109.GI.S.INT_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E110.GI.STMC.MUC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E110.GI.STMC.MUC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E111.GI.STMC.MUS_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E111.GI.STMC.MUS_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E112.THYM_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E112.THYM_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E113.SPLN_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E113.SPLN_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E114.LNG.A549.ETOH002.CNCR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E114.LNG.A549.ETOH002.CNCR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E115.BLD.DND41.CNCR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E115.BLD.DND41.CNCR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E116.BLD.GM12878_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E116.BLD.GM12878_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E117.CRVX.HELAS3.CNCR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E117.CRVX.HELAS3.CNCR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E118.LIV.HEPG2.CNCR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E118.LIV.HEPG2.CNCR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E119.BRST.HMEC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E119.BRST.HMEC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E120.MUS.HSMM_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E120.MUS.HSMM_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E121.MUS.HSMMT_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E121.MUS.HSMMT_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E122.VAS.HUVEC_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E122.VAS.HUVEC_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E123.BLD.K562.CNCR_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E123.BLD.K562.CNCR_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E124.BLD.CD14.MONO_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E124.BLD.CD14.MONO_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E125.BRN.NHA_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E125.BRN.NHA_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E126.SKIN.NHDFAD_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E126.SKIN.NHDFAD_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E127.SKIN.NHEK_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E127.SKIN.NHEK_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E128.LNG.NHLF_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E128.LNG.NHLF_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt
bedtools intersect -a mQTL.CpG.bed -b ./Roadmap_HMM_Enhancer/E129.BONE.OSTEO_HMM_Enhancer.bed.gz -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "E129.BONE.OSTEO_HMM_Enhancer\t"$0}' >> CpG_Set.Enhancer.txt

##### Step2: Estimation of CpG set (Kidney and Roadmap tissue enhancers) methylation scores for Kidney mQTL and Blood mQTL, respectively
./Software/mesc/run_mesc.py --compute-expscore-sumstat --eqtl-sumstat Kidney.mQTL.Summary.Statistics.txt --gene-sets CpG_Set.Enhancer.txt --out Kidney.mQTL.mesc.CpG_Set.Enhancer
./Software/mesc/run_mesc.py --compute-expscore-sumstat --eqtl-sumstat Blood.mQTL.Summary.Statistics.txt --gene-sets CpG_Set.Enhancer.txt --out Kidney.mQTL.mesc.CpG_Set.Enhancer

##### Step3: Estimate h2med enrichment of Kidney/Blood mQTL for each GWAS trait and each HMM state
./Software/mesc/run_mesc.py --h2med GWAS.sumstats.gz --exp-chr Kidney.mQTL.mesc.CpG_Set.Enhancer.@ --out Kidney.mQTL.mesc.CpG_Set.Enhancer.h2med
./Software/mesc/run_mesc.py --h2med GWAS.sumstats.gz --exp-chr Blood.mQTL.mesc.CpG_Set.Enhancer.@ --out Blood.mQTL.mesc.CpG_Set.Enhancer.h2med



###########################################################################
##### Estimation of CpG set (Kidney cell type-specific DARs) methylation scores and h2med for Kidney mQTL
##### Step1: Make CpG Set annotation file based on kidney cell type-specific DARs
bedtools intersect -a mQTL.CpG.bed -b ./DAR/PT_S1.cluster.1.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "PT_S1.cluster.1.DAR.7\t"$0}' > CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/PT_S2.cluster.4.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "PT_S2.cluster.4.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/PT_S3.cluster.3.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "PT_S3.cluster.3.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/LOH.cluster.2.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "LOH.cluster.2.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/DCT.cluster.8.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "DCT.cluster.8.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/PC.cluster.6.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "PC.cluster.6.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/IC.cluster.10.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "IC.cluster.10.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/Podo.cluster.11.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "Podo.cluster.11.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/Endo.cluster.5.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "Endo.cluster.5.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/Stroma.cluster.9.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "Stroma.cluster.9.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/Lymph.cluster.13.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "Lymph.cluster.13.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt
bedtools intersect -a mQTL.CpG.bed -b ./DAR/Immune.cluster.12.DAR.7.bed -u | sort | uniq | awk 'BEGIN { ORS = "\t" } { print $4}' | awk '{print "Immune.cluster.12.DAR.7\t"$0}' >> CpG_Set.snATAC.DARs.txt

##### Step2: Estimation of CpG set (Kidney HMM state) methylation scores for Kidney mQTL
./Software/mesc/run_mesc.py --compute-expscore-sumstat --eqtl-sumstat Kidney.mQTL.Summary.Statistics.txt --gene-sets CpG_Set.snATAC.DARs.txt --out Kidney.mQTL.mesc.CpG_Set.DARs

##### Step3: Estimate h2med of Kidney mQTL for each GWAS trait and each HMM state
./Software/mesc/run_mesc.py --h2med GWAS.sumstats.gz --exp-chr Kidney.mQTL.mesc.CpG_Set.DARs.@ --out Kidney.mQTL.mesc.CpG_Set.DARs.h2med


