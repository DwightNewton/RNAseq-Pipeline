# RNAseq-Pipeline

Collection of general-purpose scripts for RNAseq data analysis. Nearly all these scripts are intended for use in a high-performance computing environment (CAMH SCC) with the SLURM workload manager, for both feasibility and speed.

1. Prepare reference genome, annotations, and splicesite files.



2. Copy .fastq files from Basespace to SCC working directory.
   *

3. QC .fastq files for any technical issues using FastQC.
   *

4. Generate .bam files, aligning reads to reference genome.
   * `HISAT2_alignment_cmdgen.sh > hisat_human.cmdlist` to generate parallel processing list.
   * Submit batch jobs via SLURM: `sbatch HISAT2_parallel_script.sh`.
   * Modify files accordingly based on filenames, library #, trimming parameters, and single-end v.s. paired-end reads.

5. Generate alignment statistics using flagstat command.

6.

