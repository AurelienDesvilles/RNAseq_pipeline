#!/bin/bash

### Récupère données trimmées
### Pour analyse isomiRs + A-to-I
### Usage:
### ./5-miRge_analysis.sh <data_dir> [optional_mirge3_params]
### Exemple:
### ./5-miRge_analysis.sh ./ "-q 20 -db miRBase -on mouse"

dir=$1
optional_params=$2

## Vérification des arguments
if [ -z "$dir" ]; then
    echo "No argument 1 supplied (data directory)"
    exit 1
fi

## Paramètres par défaut si $2 vide
if [ -z "$optional_params" ]; then
    optional_params="-q 28 -lib ~/miRge3_Lib -on human -db miRBase -gff -bam -tcf -nmir -trf -ai -ie"
    echo "No optional miRge3 parameters supplied, using defaults: $optional_params"
else
    echo "Using user-provided miRge3 parameters: $optional_params"
fi

function FilesPreparation() {
    mkdir -p miRge_Analysis
    trim_data=$(find "$dir" -name "*_trimmed.fq")

    cp $trim_data miRge_Analysis/
    cd miRge_Analysis/

    for a in *.fq; do
        NEWNAME=$(echo "$a" | sed "s/_trimmed.fq/_mirge.fastq/g")
        mv "$a" "$NEWNAME"
    done
}

FilesPreparation

find -name "*_mirge.fastq" > mirge_files.txt

## Activation environnement conda
eval "$(conda shell.bash hook)"
conda activate mirge3

## Exécution miRge3 avec paramètres dynamiques
miRge3.0 $optional_params -s mirge_files.txt -o miRge

conda deactivate

exit 0

