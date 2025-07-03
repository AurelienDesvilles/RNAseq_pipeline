# RNAseq_pipeline

#### Author: Romain Piucco
#### Co-author: Aurélien Desvilles
#### Supervisor: Sébastien Hergalant

## Small RNAseq pipeline 
This repository contains a series of Bash and R scripts designed for processing small RNA-seq data. The pipeline includes quality control (QC), trimming, alignment, read counting, analysis using miRDeep2 and miRge3.0, count and normalization.
 
##  Project Structure
.
├── 1-QC_trimming.sh         # QC and trimming of raw FASTQ files
├── 2-Alignment.sh           # Alignment of trimmed reads using HISAT2
├── 3-CountStats.sh          # Read counting and alignment statistics
├── 4-miRDeep_analysis.sh    # Novel and known miRNA discovery using miRDeep2
├── 5-miRge_analysis.sh      # isomiR and A-to-I editing analysis using miRge3.0
├── 6-Count_miRDeep.r        # Aggregates miRDeep2 results into count matrices
├── 7-DistribQuality.r       # Quality control and expression distribution summaries
├── raw_data/                # Input FASTQ files
└── outputs/                 # Generated count tables,  distribution and quality files
 
⚙️ Requirements
•	conda
•	Conda environments:
o	mirge: miRDeep2
o	mirge3: miRge3.0
•	Bioinformatics tools : 
o	Fastqc (version 0.11.9)
o	 trim_galore (version 0.6.10)
o	hisat2 (version 2.2.1)
o	 samtools (version 1.13)
o	 htseq-count (version 1.99.2)
o	Bowtie (version 1.1.2)
•	R packages:
o	ggplot2
o	stringr
o	dplyr
o	tidyverse
o	ComplexHeatmap
o	reshape2

•	Reference genome to index :
Her ewe use Homo_sapiens.GRCh38.dna.primary_assembly.fa ( ref: https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/ )
•	miRNA library miRge3_Lib (ref: https://sourceforge.net/projects/mirge3/files/miRge3_Lib/ )
•	mature and hairpin miRNA sequences (ref: https://www.mirbase.org/download/ )
•	Genome coordinates file .gff3  (ref: https://www.mirbase.org/download/ )
 
🚀 Usage
Preliminary step: 
1: configure reference genome index
$ sed '/^>/ s/ /_/g'  <reference_genome.fa>  > tmp && mv tmp <reference_genome.fa>
    bowtie-build <reference_genome.fa> hgBowtie
    mkdir hgBowtie
    mv hgBowtie*ewbt hgBowtie

•	Replace spaces by “_” in headers so that script recognize file as fasta 
•	Create several index files with “hgBowtie” prefix
•	Create hgBowtie directory
•	Move all index files in this directory


2: convert genome coordinates .gff3 file to .gtf
If you don’t have a .gtf file, you can convert an other genome coordinate file with AGAT.
$ agat_convert_sp_gff2gtf.pl  --gff  <input_genome_coordinates.gff3> -o  <output_genome_coordinates.gtf>

 
1. QC and Trimming
$ bash 1-QC_trimming.sh <raw_data_directory>  <trimgalore>  <”config”>

•	Performs FastQC on raw FASTQ files (in raw data directory)
•	Trims adapters using Trim Galore
•	Post-trimming QC
•	Organizes output into:
o	0-Data_raw/
o	1-Data_QC/
o	2-Data_Trimming/


 
2. Alignment
$ bash 2-Alignment.sh <Data_Trimming_directory > <miRNA_files_index >  <config>
•	Aligns *_trimmed.fq reads using HISAT2
•	Moves SAM files to 3-Data_Alignment/
•	Converts SAM to sorted BAM and indexes them


 
3. Count and Stats
$ ./3-CountStats.sh <data_alignment_directory>  <genome_coordinates.gtf> <config>
•	Runs:
o	samtools idxstats
o	htseq-count (with provided GTF)
o	samtools flagstat
•	Organizes results into:
4-Stats_Count/
├── A-Stats_Counts/
├── B-HTseq_Counts/
└── C-Flagstats_Counts/

 
4. miRDeep2 Analysis 
$ ./4-miRDeep_analysis.sh  <data_trimming_directory> <config>
•	Automatically finds *_trimmed.fq files
•	Performs mapping with mapper.pl
•	Runs miRDeep2.pl on collapsed reads using:
o	Human genome
o	Mature miRNA reference (human and other species)
o	Hairpin sequences
•	Outputs one folder per input file in miRDeep_Analysis/
Requires:
•	refBowtie.fa: reference indexed genome
•	miRNAs_mature_refHg38.fa: Homo sapiens mature miRNAs sequences
•	mature_other.fa: other species mature miRNAs sequences
•	miRNAs_hairpin_refHg38.fa: Homo sapiens hairpin miRNAs sequences
Warning: This step may take several hours. Also, make sure you have enough memory space.


 
5. miRge3.0 Analysis to find isomiRNA
$ ./5-miRge_analysis.sh <data_trimming_directory> <config>
Finds *_trimmed.fq, renames to *_mirge.fastq
•	Runs miRge3.0 analysis with:
o	A-to-I editing (-ai), isomiRs (-ie), BAM output, and other detailed reports
o	Library: ~<your path> /miRge3_Lib
o	To find isomiRNA
•	Output in miRge/
Outputs include: .bam, .gff, .tsv, and quality reports.
Warning: Make sure you have enough memory space.


 
6. Count Matrix Aggregation (miRDeep2 Results)
$ Rscript 6-Count_miRDeep.r <dir> <split_index> <file_pattern>
•	Aggregates all miRDeep2 .csv results into combined count matrices.
•	Parameters:
o	dir: File is usually miRDeep_Analysis 
o	split_index: Integer to extract sample names from file paths
o	file_pattern: File pattern to match result files (e.g. result_*.csv)
•	Outputs:
o	All_miRDeep2_results_raw_read_counts.csv: raw counts per miRNA
o	All_miRDeep2_results_normalizedRPM_read_counts.csv: RPM-normalized counts


 
7. Distribution and CQ Summary (Count Matrix)
$ Rscript 7-DistribQuality.r <count_matrix.tsv> <output_dir>
•	Performs extensive quality control on miRNA count matrices (raw or RPM)
•	Generates:
o	Correlation plots of samples to the global median
o	Expression intensity summaries (percentiles, quantiles)
o	RLE plots (miRNAs_CQ_RLE.png)
o	Density plots (miRNAs_CQ_densities.png)
o	Quantile distribution histograms per sample
o	Barplots of 0-value proportions per sample and variable
o	Binary heatmap of 0 vs. non-zero values (0Values_Heatmap.png)
Input matrix must have sample columns starting from column 3.

 
✅ Output Overview
Each step generates specific subfolders:
•	0-Data_raw/ : Original FASTQ files
•	1-Data_QC/ : FastQC reports pre-trimming
•	2-Data_Trimming/ : Trimmed FASTQ + post-trimming QC
•	3-Data_Alignment/ : Aligned and sorted BAM files
•	4-Stats_Count/ : Count and QC metrics in TSV format
•	miRDeep_Analysis/ : Per-sample miRDeep2 folders
•	miRge/ : miRge3.0 output reports and files
•	7-Quentiles_plots/:  QC files
 
📬 Contact
For questions or suggestions, please contact Romain Piucco at romain.piucco@univ-lorraine.fr, Aurélien Desvilles at desvillesaurelien@gmail.com.


