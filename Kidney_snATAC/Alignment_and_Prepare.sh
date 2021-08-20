###########################################################################
#####    Alignment of snATAC-seq reads to human genome                #####
#####    Software: Cell Ranger ATAC (v1.1.0) and snaptools (v1.4.8)   #####
###########################################################################
#Reads Align
cellranger-atac count --id=$SAMPLEID --reference=$REF_GENOME --fastqs=$OUTPUT_PATH --sample=$SAMPLE_ID
#Format prepare
snaptools snap-add-bmat --snap-file=$SNAP_FILENAME.snap --bin-size-list 1000 5000 10000 --verbose=True

