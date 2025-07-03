#!/bin/bash

## Compte la sortie idxstats + htseq + flagstats
## Range le tout dans des dossiers et renomme les fichiers
## Usage:
## ~/Transcriptomics/3-CountStats.sh <data_dir> <gtf_file> [optional_htseq_params]
## Example:
## ~/Transcriptomics/3-CountStats.sh ./ hsa_mir.gtf "--stranded yes --type exon"

dir=$1
mirna_gtf=$2
optional_params=$3

# Vérifications des arguments obligatoires
if [ -z "$dir" ]; then
    echo "No argument 1 supplied (data directory)"
    exit 1
fi

if [ -z "$mirna_gtf" ]; then
    echo "No argument 2 supplied (GTF file)"
    exit 1
fi

# Définir les paramètres par défaut si $3 vide
if [ -z "$optional_params" ]; then
    optional_params="--quiet --format bam --order pos --mode union --stranded no --type exon --idattr gene_id"
    echo "No optional htseq-count parameters supplied, using defaults: $optional_params"
else
    echo "Using user-provided htseq-count parameters: $optional_params"
fi

function CountStats(){
    align_data=$(find "$dir" -name "*.bam")
    echo "Files to process: $align_data"

    for o in $align_data; do
        echo "Processing: $o"

        samtools idxstats "$o" > "${o}.stats.count.tsv"
        htseq-count $optional_params "$o" ~/Bureau/Genomes/"$mirna_gtf" > "${o}.count.whole.tsv"
        samtools flagstat "$o" > "${o}.flagstats.tsv"
    done
}

CountStats

function Arrange(){
    mkdir -p 4-Stats_Count/A-Stats_Counts 4-Stats_Count/B-HTseq_Counts 4-Stats_Count/C-Flagstats_Counts

    A_data=$(find "$dir" -name "*.stats.count.tsv")
    B_data=$(find "$dir" -name "*.count.whole.tsv")
    C_data=$(find "$dir" -name "*.flagstats.tsv")

    mv $A_data 4-Stats_Count/A-Stats_Counts/
    mv $B_data 4-Stats_Count/B-HTseq_Counts/
    mv $C_data 4-Stats_Count/C-Flagstats_Counts/
}

Arrange

function renameCounts(){
    tsv_files=$(find 4-Stats_Count -name "*.tsv")

    for u in $tsv_files; do
        NEWNAME=$(echo "$u" | sed "s/.bam//g")
        mv "$u" "$NEWNAME"
    done
}

renameCounts

exit 0

