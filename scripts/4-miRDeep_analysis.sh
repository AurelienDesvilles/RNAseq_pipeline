#!/bin/bash

## Sur données trim
## Va directement chercher les fichiers _trimmed.fq issus du script #1
## Fichiers miRNAs de référence dans le dossier Transcriptomics
## Crée un dossier pour chaque fichier avec les résultats de mapping et des comptes
## Usage:
## ~/Transcriptomics/4-miRDeep_analysis.sh <data_dir> [optional_mapper_params]
## Exemple:
## ~/Transcriptomics/4-miRDeep_analysis.sh ./ "-h -i -m"

dir=$1
optional_params=$2

## Vérification des arguments
if [ -z "$dir" ]; then
    echo "No argument 1 supplied (data directory)"
    exit 1
fi

## Paramètres par défaut si $2 vide
if [ -z "$optional_params" ]; then
    optional_params="-e -h -i -m -v -o 8"
    echo "No optional mapper.pl parameters supplied, using defaults: $optional_params"
else
    echo "Using user-provided mapper.pl parameters: $optional_params"
fi

## Activer l'environnement mirge
eval "$(conda shell.bash hook)"
conda activate mirge

function mapping(){
    mkdir -p miRDeep_Analysis
    trim_data=$(find "$dir" -name "*_trimmed.fq")

    cd miRDeep_Analysis/

    for i in $trim_data; do
        NEWNAME=${i##*/}
        DIRNAME=$(echo "$NEWNAME" | sed "s/_trimmed.fq//g")
        mkdir -p "$DIRNAME"
        cd "$DIRNAME"

        echo "Processing: $i"
        echo "Directory: $DIRNAME"

        mapper.pl ../../"$i" $optional_params -p ~/Transcriptomics/genome/hgBowtie/refBowtie.fa \
            -s reads_collapsed.fa -t reads_vs_refdb.arf

        cd ../
    done
}

mapping

function miRDeepAnalyze(){
    DIRS=$(ls -d */)
    for f in $DIRS; do
        echo "Analyzing: $f"
        cd "$f"
        miRDeep2.pl reads_collapsed.fa \
            ~/Transcriptomics/genome/new_Homo_sapiens.GRCh38.dna.primary_assembly.fa \
            reads_vs_refdb.arf \
            ~/Transcriptomics/genome/miRNAs/miRNAs_mature_refHg38.fa \
            ~/Transcriptomics/genome/miRNAs/mature_other.fa \
            ~/Transcriptomics/genome/miRNAs/miRNAs_hairpin_refHg38.fa \
            -t hsa 2> report.log
        cd ../
    done
}

miRDeepAnalyze

conda deactivate

exit 0

