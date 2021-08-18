###########################################################################
#####    mQTL analysis based on imputed genotyping and methylation    #####
#####    Sample size: 443 human kidney samples                        #####
###########################################################################
library(ggplot2)
library("MatrixEQTL")

### Parameters for QTL
SNP_file_name = ("SNP_matrix.txt")
expression_file_name = ("Methylation_matrix.txt")
covariates_file_name = ("Covariates.txt")
file.create("mQTL.txt")
output_file_name = ("mQTL.txt")

useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
pvOutputThreshold = 1e-2;
errorCovariance = numeric();

### SNP
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the space character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10000;     # read file in pieces of 10000 rows
snps$LoadFile( SNP_file_name );

### Gene
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 10000;     # read file in pieces of 10000 rows
gene$LoadFile( expression_file_name );

### Covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 10000;     # read file in pieces of 10000 rows
cvrt$LoadFile( covariates_file_name );


### Parameters for Cis-QTL
#http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html#cis

snps_location_file_name = ("snpsloc.txt")
gene_location_file_name = ("Probeloci.txt");
# Output file name
file.create("cis_mQTL.txt")
file.create("trans_mQTL.txt")
output_file_name_cis =("cis_mQTL.txt")
output_file_name_tra = ("trans_mQTL.txt")

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;     ### Output all Cis-QTL
pvOutputThreshold_tra = 0;     ### Do No output Trans-QTL 

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

### Cis-QTL and Trans-QTL
me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
nrow(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
nrow(me$trans$eqtls)

## Output cis-mQTL results
me$cis$eqtls[,3] <- signif(me$cis$eqtls[,3], 4)
me$cis$eqtls[,4] <- signif(me$cis$eqtls[,4], 4)
me$cis$eqtls[,5] <- signif(me$cis$eqtls[,5], 4)
me$cis$eqtls[,6] <- signif(me$cis$eqtls[,6], 4)
write.table(me$cis$eqtls,output_file_name_cis, col.names=F, row.names = F, quote=F,sep="\t")

## Plot the Q-Q plot of local and distant p-values

pdf(file='QQ plot local and distant p-values.pdf')
plot(me)
dev.off()
