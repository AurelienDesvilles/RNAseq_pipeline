#!/bin/bash

## Usage:
## ~/Transcriptomics/2-Alignment.sh <data_dir> <hisat2_index> [optional_hisat2_params]
## Example:
## ~/Transcriptomics/2-Alignment.sh ./ mature_hg38hisat2_mirnas "-k 5 --no-spliced-alignment"

dir=$1
mirna_files=$2
optional_params=$3

# Vérifications des arguments obligatoires
if [ -z "$dir" ]; then
    echo "No argument 1 supplied (data directory)"
    exit 1
fi
if [ -z "$mirna_files" ]; then
    echo "No argument 2 supplied (hisat2 index)"
    exit 1
fi

# Définir les paramètres par défaut si $3 vide
if [ -z "$optional_params" ]; then
    optional_params="--dta --rna-strandness F"
    echo "No optional hisat2 parameters supplied, using defaults: $optional_params"
else
    echo "Using user-provided hisat2 parameters: $optional_params"
fi

alignment(){
    trim_data=$(find "$dir" -name "*_trimmed.fq")
    echo "You chose this index to align your miRNAs: $mirna_files"

    for z in $trim_data; do
        echo "Aligning: $z"
        OPNAME=$(echo "$z" | sed "s/_trimmed.fq/_align/")
        hisat2 -p 8 -q -U "$z" -x "$mirna_files" $optional_params -S "$OPNAME.sam"
    done
}

alignment

sam_files=$(find "$dir" -name "*.sam")
mkdir -p 3-Data_Alignment
mv $sam_files 3-Data_Alignment/

ConvertIndex(){
    align_data=$(find 3-Data_Alignment -name "*.sam")

    for y in $align_data; do
        samtools sort -@ 8 -o "$y.bam" "$y"
        GOODNAME=$(echo "$y.bam" | sed "s/_align.sam.bam/_align.bam/")
        mv "$y.bam" "$GOODNAME"
    done

    find 3-Data_Alignment/*.bam -exec samtools index {} \;
}

ConvertIndex

exit 0

