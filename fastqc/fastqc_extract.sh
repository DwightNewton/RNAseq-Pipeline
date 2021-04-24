#First, cd to fastqc output location and run # unzip "*.zip"
SRC_DIR="/external/mgmt1/scratch01/tmpdata02/dnewton/PITT_tetrad_SCT_RNAseq/fastq/pool_1/Samples/fastqc/"
out_dir="/external/mgmt1/scratch01/tmpdata02/dnewton/PITT_tetrad_SCT_RNAseq/fastq/pool_1/Samples/fastqc/fqc_res/"
cd $SRC_DIR
for file in *.zip;
do
Base_Name=${file%.zip}
echo "cp $SRC_DIR/${Base_Name}/summary.txt $out_dir/$Base_Name.txt"
done