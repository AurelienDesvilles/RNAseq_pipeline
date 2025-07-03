#!/bin/bash


## In the dir raw_data ##
## Script de QC / Trimming (+ post trimming QC) ##
## Choisir la config pour trim galore / $2=trimgalore or -a / $3=config (without output dir) ##
## Range les fichiers .fastq dans un dossier 0-Raw_Data, le QC dans 1-Data_QC et le trimming dans 2-Data_Trimming ##

## /home/romain/Transcriptomics/1-QC_trimming.sh ./ trimgalore "--small_rna --quality 28 --max_length 25 --phred33" ##

## Améliorations : Choix du nom de l'output dir,, faire appel à cutadapt ##

dir=$1
trim_tool=$2
config=$3

if [ -z "$1" ];then
    echo "No argument 1 supplied"
fi
if [ -z "$2" ];then
    echo "No argument 2 supplied"
fi
if [ -z "$3" ];then
    echo "No argument 3 supplied"
fi

function QC(){
	echo 'Phase I - QC starting'
	~/Transcriptomics/Tools/FastQC/fastqc *.fastq
	mkdir 1-Data_QC 
	mv *fastqc* -t 1-Data_QC/ 
}

QC $dir

function trimming(){
	echo 'Phase II - trimming phase - using trimgalore'
	
	for i in *.fastq 
		do  
			echo $i
			NEWNAME=`echo $i | sed "s/.fastq/_trimmed/g"`
			echo $NEWNAME
			
			if [ "$2" = "trimgalore" -o "$2" = "-a" ]
			then
				 echo $config
				 ~/Transcriptomics/Tools/TrimGalore-0.6.10/trim_galore $config --output_dir ./$NEWNAME $i
			else 
				echo 'non c pas bon'
			fi
		
		done
		
	mkdir 2-Data_Trimming/ ### Création dir & déplacement des fichiers ###
	mv *_trimmed -t 2-Data_Trimming/
	cd 2-Data_Trimming/ 
	
	trim_dirs=$(ls -d SRR[0-9]*) ### fastQC for trimming data ###
	for j in $trim_dirs
		do
			#echo $j
			cd $j 
			~/Transcriptomics/Tools/FastQC/fastqc *_trimmed.fq
			cd ../
		done  	   
	cd ../	
}

eval "$(conda shell.bash hook)"
conda activate mirge

trimming $dir $trim_tool $config

mkdir 0-Data_raw
mv *.fastq ./0-Data_raw

conda deactivate  
