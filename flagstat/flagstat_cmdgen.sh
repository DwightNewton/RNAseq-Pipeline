#!/bin/sh
SRC_DIR="/external/rprshnas01/eslab/dnewton/PITT_tetrad_bam_5_3__3_15/all_bam_files/"
out_dir="/external/rprshnas01/eslab/dnewton/PITT_tetrad_bam_5_3__3_15/all_bam_files/flagstat_output/"

cd $SRC_DIR
for file in *.bam;
do
out_file=${file%}
out_name=${out_dir}/${out_file}.txt
Stat_File=${Base_Name}

echo "samtools flagstat -@ 8 ${file} > ${out_name}"

done
   