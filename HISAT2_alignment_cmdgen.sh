#!/bin/sh

#idx_dir is the reference genome
idx_dir="/external/mgmt1/scratch01/tmpdata02/dnewton/PITT_tetrad_SCT_RNAseq/Reference/hisat2//Hs_release_98"
splice_dir="/external/mgmt1/scratch01/tmpdata02/dnewton/PITT_tetrad_SCT_RNAseq/Reference/splicesites.tsv"

#Location of .fastq files and desired output directory
SRC_DIR="/external/rprshnas01/eslab/dnewton/pool_5/"
out_dir="/external/rprshnas01/eslab/dnewton/PITT_tetrad_bam_5_3__3_15/pool5/"

cd $SRC_DIR
for file in *L001_R1_001.fastq.gz;
do
Base_Name=${file%L001_R1_001.fastq.gz}
out_file=${file%_S*_L001_R1_001.fastq.gz}
out_name=${out_dir}/${out_file}.bam

L1R1_File=${Base_Name}L001_R1_001.fastq.gz
L2R1_File=${Base_Name}L002_R1_001.fastq.gz
L1R2_File=${Base_Name}L001_R2_001.fastq.gz
L2R2_File=${Base_Name}L002_R2_001.fastq.gz

#Alignment step with default parameters, for paired-end reads. 6 bases on 5' end need to be trimmed as automatic adapter trimming was turned off in basespace.
#Manual review of fastqc results showed variable quality of first 5 bases on 5' end, and last 15 bases of 3' ends - thus these were trimmed
echo "hisat2 -x ${idx_dir} --known-splicesite-infile ${splice_dir} -p 8 -1 ${SRC_DIR}/${L1R1_File},${SRC_DIR}/${L2R1_File} -2 ${SRC_DIR}/${L1R2_File},${SRC_DIR}/${L2R2_File} --trim5 11 --trim3 15 | samtools view -bS - > ${out_name}"

done