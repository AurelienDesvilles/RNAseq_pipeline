#!/bin/bash 

## ./Pipeline.sh raw_seq_directory trimgalore "--small_rna --quality 28 --max_length 25 --phred33" ##
dir=$1
trim_tool=$2 #trimgalore
config=$3 # "--small_rna --quality 28 --max_length 25 --phred33"
mature_ref=$4
gtf_file=$5


echo "Trimming launch..."
bash 1-QC_trimming.sh $dir $trim_tool $config
echo "Alignment of trimmed sequences..."
bash 2-Alignement.sh 2-Data_Trimming $mature_ref #default config
echo "Stat count, it could take a litle time..."
bash 3-CountStats.sh 3-Data_Alignment $gtf_file #default config
echo "miRDeep2 analysis, it could take a lot of time..."
bash 4-miRDeep_analysis.sh 2-Data_Trimiming #default config
echo "Counting miRDeep2 results..."
Rscript --vanilla ~/Transcriptomics/6-Count_miRDeep.r miRDeep_Analysis 2 miRNA_expressed_all_samples
Rscript --vanilla ~/Transcriptomics/7-DistribQuality.r All_miRDeep2_results_normalizedRPM_read_counts.csv miRDeep_Quality
echo "miRge3.0 analysis, it could take a litle time..."
bash ~/Transcriptomics/5-miRge_analysis.sh 2-Data_Trimiming #default config
echo "Jobs completed. Thank you for your patience."

