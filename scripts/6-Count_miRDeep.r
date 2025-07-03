#!/bin/Rscript

## Prend tous les DIR en sortie du script miRDeep 2 ##
## Créé les matrices de comptes pour les échantillons ## 

suppressMessages(library(stringr))

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){print('No args provided - Need dir')}
if (length(args)==1){print('I need to know how to split path')}
if (length(args)==2){print('I need to know what is the pattern for find files')}
if (length(args)==3){dir <- args[1];split=args[2];pattern=args[3]}

dir=dir
split=as.numeric(split)
pattern=pattern

Trscp_files <- list.files(list.dirs(dir), pattern=pattern, full.names=T)
data_files <- lapply(Trscp_files, function(x) read.csv(x, sep='\t',dec='.',header=T, stringsAsFactors=F))
vec <- vector();for(i in Trscp_files){vec <- c(vec, str_split(i, '/')[[1]][[split]])} ##7 pour le '/', peut changer en fonction de la racine et du dir
names(data_files) <- vec
print(vec)

for (i in 1:length(data_files)){
	data <- data_files[[i]]
	data <- data[,-c(4:5)]
	data <- data[,c(1,3,2,4)]
	data_files[[i]] <- data
}

results <- data.frame(lapply(data_files, function(x) cbind(x)))
results <- results[,c(grep('X.miRNA',colnames(results))[1], grep('precursor',colnames(results))[1], grep('read_count', colnames(results)), grep('seq.norm', colnames(results)))]

colnames(results)[grep('X.miRNA', colnames(results))] <- 'miRNA'
colnames(results)[grep('precursor', colnames(results))] <- 'precursor'
colnames(results) <- str_replace(colnames(results), '.read_count', '_Read_Count')
colnames(results) <- str_replace(colnames(results), '.seq.norm', '_Norm_Count_RPM')
colnames(results) <- str_replace(colnames(results), '\\.', '')

raw_results <- results[,grep('miRNA|precursor|Read_Count', colnames(results))]
norm_results <- results[,grep('miRNA|precursor|Norm_Count', colnames(results))]
colnames(raw_results) <- str_replace(colnames(raw_results), '_Read_Count', '')
colnames(norm_results) <- str_replace(colnames(norm_results), '_Norm_Count_RPM', '')

write.table(raw_results, './All_miRDeep2_results_raw_read_counts.csv', col.names=NA, row.names=T, sep='\t',dec=',')
write.table(norm_results, './All_miRDeep2_results_normalizedRPM_read_counts.csv', col.names=NA, row.names=T, sep='\t',dec=',')
