# RNAseq-Pipeline

Collection of general-purpose scripts for RNAseq data analysis. Nearly all these scripts are intended for use in a high-performance computing environment (CAMH SCC) with the SLURM workload manager, for both feasibility and speed.

1. Prepare reference genome, annotations, and splicesite files.
   * Refer to `reference_genome.sh`.

2. Copy .fastq files from Basespace to SCC working directory.
   * Mount Basespace in command line: `basemount basespace`
   * Navigate to Projects/, then appropriate run then /Samples/ folder.
   * Copy .fastq files to working directory: `find ./ -name "*.fastq.gz" -exec cp \{\} "path/to/output/directory" \;`
   * If permissions error occurs, `sudo` the `find` and `cp` commands.

3. QC .fastq files for any technical issues using FastQC.
   * `cd` to fastqc output location and run `unzip *.zip`
   * Run fastqc in parallel with `fastqc_cmdgen.sh` and `fastqc_parallelization.sh`
   * Run `fastqc_extract.sh`

4. Generate .bam files, aligning reads to reference genome.
   * `HISAT2_alignment_cmdgen.sh > hisat_human.cmdlist` to generate parallel processing list.
   * Submit batch jobs via SLURM: `sbatch HISAT2_parallel_script.sh`.
   * Modify files accordingly based on filenames, library #, trimming parameters, and single-end v.s. paired-end reads.

5. Generate alignment statistics using flagstat command.

6.

