# PerturbSci_Kinetics
Reads processing scripts for PerturbSci-Kinetics.

### Background SNP calling (/bg_SNP_calling/SNP_calling_processing.sh)
##### Key parameter
Fastq input: Paired-end full-coverage bulk RNA-seq
Sample ID: a text file containing prefixes of each sample. R1 and R2 files from the sample should share the prefixes.
Reference fasta: the fasta file of the reference genome.
Index: the STAR index of the reference genome.
Output folder: the directory for all output files
script folder: the folder for all sub scripts

##### Steps
1. Trim adapter sequences by automatic detection.
2. STAR alignment.
3. Filter aligned reads.
4. Merge bams and sort the merged bam.
5. Summarize the base identities of reads mapped to each genomic location. 
6. Inherent SNP calling.

### Single cell whole/nascent transcriptomes reprocessing steps (/whole_tx_processing/Main_processing.sh)
Fastq input: Paired-end PerturbSci-Kinetics demultiplexed whole transcriptome fastq files. 
1. Change file names of fastq to make them callable in the following steps.
2. Attach UMI sequences on R1 to headers of R2.
3. Trim potential polyA sequences from the 3'end of R2.
4. STAR alignment.
5. Filter aligned reads.
6. PCR duplicates removal based on both mapped genomic coordinates and UMI.
7. Single-cell sam files generation.
8. Transform the alignment information in single-cell sams to tables at the single-base level.
9. Identify T>C mutations on each single read and extract read names of nascent reads.
10. Extract nascent reads from single-cell sams.
11. Gene-level feature counting and re-format the single-cell gene expression matrix.

### Single cell sgRNA reprocessing steps (/sgRNA_reads_processing/sgRNA_processing.sh)
Fastq input: Paired-end PerturbSci-Kinetics demultiplexed sgRNA fastq files. 
1. Change file names of fastq to make them callable in the following steps.
2. One-step sgRNA identification, de-duplication, and counting.
3. Re-format the single-cell sgRNA expression matrix.
