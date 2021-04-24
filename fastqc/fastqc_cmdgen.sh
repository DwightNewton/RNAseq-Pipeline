#!/bin/sh
#SRC_DIR=/genome/scratch/Sibille/dnewton/RapidRun/


#hisat_cmd="hisat2 -x ${idx_dir} --known-splicesite-infile ${splice_dir}"
SRC_DIR="/external/mgmt1/scratch01/tmpdata02/dnewton/PITT_tetrad_SCT_RNAseq/fastq/pool_1/Samples/"
out_dir="/external/mgmt1/scratch01/tmpdata02/dnewton/PITT_tetrad_SCT_RNAseq/bam_pool1_5-15bp_3-15-bp/"

cd $SRC_DIR
for file in *.fastq.gz;
do
Base_Name=${file%.fastq.gz}
echo "fastqc –o ${SRC_DIR}/${Base_Name}.fastq.gz"

done




