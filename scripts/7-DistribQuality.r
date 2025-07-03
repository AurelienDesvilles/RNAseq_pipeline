#!/bin/Rscript 

## Script pour Contrôle Qualité en sortie des comptes ## 
## Prend une matrice de comptes (normalisée ou non) ## 
## 2 colonnes en ID + samples ## 

suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(reshape2))

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){print('No args provided - Need matrix of counts')}
if (length(args)==1){print('No Output dir supply')}
if (length(args)==2){data <- args[1];Output <- args[2];print('Checking data distribution.................................................................................')}

results <- read.csv(data, sep='\t', dec=',', header=T,stringsAsFactors=F)[,-1]

dir.create(Output)
setwd(Output)

if(length(results[,3:max(ncol(results))])!=1){
     mediane <- apply(results[,3:max(ncol(results))],1,median,na.rm=TRUE)
     nb.categories <- trunc(log10(dim(results[,3:max(ncol(results))])[1]))
     stats.matrix <- matrix(nrow=1+2*(nb.categories+1),ncol=dim(results[,3:max(ncol(results))])[2])
     colnames(stats.matrix) <- colnames(results[,3:max(ncol(results))])
     intervals <- 0:nb.categories
     petites.valeurs <- paste("pv-",as.character(10^intervals),sep="")
     grandes.valeurs <- paste("gv-",as.character(10^rev(intervals)),sep="")
     rownames(stats.matrix) <- c("coeff.corr",petites.valeurs,grandes.valeurs)
     mat_cor <- results
     mat_dis <- results
     mat_dis[,3:ncol(mat_dis)] <- apply(mat_dis[,3:ncol(mat_dis)],1, function(x) log2(x+1))
     for(i in (1:dim(results[,3:max(ncol(results))])[2])){
         stats.matrix[1,i] <- cor(mat_cor[,3:max(ncol(mat_cor))][,i],mediane,use="complete.obs")
         mat <- sort(mat_dis[,3:max(ncol(mat_dis))][,i],na.last=NA)
         for(j in intervals){
             stats.matrix[j+2,i] <- mat[10^j]
             stats.matrix[dim(stats.matrix)[1]-j,i] <- mat[length(mat)-10^j]
         }
     }
}
write.table(stats.matrix, 'QualityControl_Data_Distribution.txt')

labels <- colnames(results[,3:ncol(results)])
labels <- unlist(lapply(colnames(results[,3:ncol(results)]), function(x) str_split(x, '_')[[1]][[1]]))

pdf('QC_Correlation.pdf')
print(plot(stats.matrix[1,],ylab="correlation coeff.",ylim = c(0,1), xlab="",lwd=1,col="blue",type="l", xaxt="n"))
axis(side=1, labels = labels, at = 1:length(results[,3:ncol(results)]), cex.axis = 0.75, las=2)
dev.off()

#pdf('Data_Distribution.pdf')
#print(plot(stats.matrix[2,],ylab="Intensité d'expression.",ylim = c(0,30), xlab="",lwd=1,col="blue",type="l", xaxt="n"))
#lines(stats.matrix[3,], col = "Red" )
#lines(stats.matrix[4,], col = "Darkgreen" )
#lines(stats.matrix[5,], col = "Green" )
#lines(stats.matrix[6,], col = "Orange" )
#lines(stats.matrix[7,], col = "Grey" )
#lines(stats.matrix[8,], col = "Black" )
#lines(stats.matrix[9,], col = "Purple" )
#legend(1,30,legend=c("PV", "PV10", "PV100", "PV1000", "PG1000", "PG100", "PG10", "PG"), fill = c("Blue", "Red", "Darkgreen", "Green", "Orange", "Grey","Black","Purple"))
#axis(side=1, labels = labels, at = 1:length(results[,3:ncol(results)]), cex.axis = 0.75, las=2)
#dev.off()

## RLE ##
data_rle <- results 
data_rle[,3:dim(data_rle)[2]] <-  t(apply(data_rle[,3:dim(data_rle)[2]],1, function(x) log2(x+1)-log2(median(x+1))))
mel <- melt(data_rle)
png(filename="miRNAs_CQ_RLE.png", width=1280, height=800)
print(ggplot(mel, aes(x=variable, y=value)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90)))
dev.off()

## Density ## 
data_density <- results 
data_density[,3:dim(data_density)[2]] <- t(apply(data_density[,3:dim(data_density)[2]],1, function(x) log2(x+1)))
png(filename="miRNAs_CQ_densities.png", width=1200, height=800)
melta <- melt(data_density)
print(ggplot(melta, aes(x=value, color=variable)) + geom_density())
dev.off()

## Quantiles ##
matrix_qtle <- data.frame(matrix(ncol=ncol(results[,3:dim(results)[2]]), nrow=100))
colnames(matrix_qtle) <- colnames(results[,3:dim(results)[2]])

dir.create('Quantiles_plots')
setwd('./Quantiles_plots/')
for (i in colnames(results[,3:dim(results)[2]])){
	NAME=str_split(i, '_')[[1]][1]
	png(filename=paste0(NAME, "_miRNAs_CQ_quantiles.png"), width=1200, height=800)
	hist(quantile(log2(results[,i]+1), probs=seq(0.01, by=0.01)), breaks=25, xlim=c(0,20), xlab=paste0(NAME, '_QUANTILES'), main=paste0('Histogram of quantiles count for ', NAME))
	dev.off()
	qtle <- matrix(quantile(log2(results[,i]+1), probs=seq(0.01, by=0.01)))
	matrix_qtle[,i] <- matrix(quantile(log2(results[,i]+1), probs=seq(0.01, by=0.01)))
}
matrix_qtle <- cbind('Frequency'=paste0(rownames(matrix_qtle), '%'), matrix_qtle)
melt_qtle <- melt(matrix_qtle)

setwd('../')

## Contrôle des 0 / plots / HM ## 
is_empty=function(x){
	new_vec <- vector()
	for (i in x){
		if(i==0){
			new_vec <- c(new_vec, 'TRUE')
		}else{	new_vec <- c(new_vec, 'FALSE')
		}
	}
return(new_vec)
}

missing.values <- as.data.frame(results[,3:dim(results)[2]]) %>% 
  gather(key = "key", value = "val") %>%
  mutate(iszero = is_empty(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, iszero) %>%
  summarise(num.iszero = n()) %>%
  mutate(pct = num.iszero / total * 100)
levels <- (missing.values  %>% filter(iszero == T) %>% arrange(desc(pct)))$key
percentage.plot <- missing.values %>%
      ggplot() +
        geom_bar(aes(x = reorder(key, desc(pct)), 
                     y = pct, fill=iszero), 
                 stat = 'identity', alpha=0.8) +
      scale_x_discrete(limits = levels) +
      scale_fill_manual(name = "", 
                        values = c('steelblue', 'tomato3'), labels = c("Present", "Missing")) +
      coord_flip() +
      labs(title = "Percentage of 0 values", x =
             'Variable', y = "% of 0 values")
pdf('Samples_Repart_Valeurs_0.pdf')
print(percentage.plot)
dev.off()

missing.values <- as.data.frame(t(results[,3:dim(results)[2]])) %>% 
  gather(key = "key", value = "val") %>%
  mutate(iszero = is_empty(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, iszero) %>%
  summarise(num.iszero = n()) %>%
  mutate(pct = num.iszero / total * 100)
levels <- (missing.values  %>% filter(iszero == T) %>% arrange(desc(pct)))$key
percentage.plot <- missing.values %>%
      ggplot() +
        geom_bar(aes(x = reorder(key, desc(pct)), 
                     y = pct, fill=iszero), 
                 stat = 'identity', alpha=0.8) +
      scale_x_discrete(limits = levels) +
      scale_fill_manual(name = "", 
                        values = c('steelblue', 'tomato3'), labels = c("Present", "Missing")) +
      coord_flip() +
      labs(title = "Percentage of 0 values", x =
             'Variable', y = "% of 0 values")
pdf('Var_Repart_Valeurs_0.pdf')
print(percentage.plot)
dev.off()

# HM # 

hm <- results
hm_map <- hm[,3:dim(hm)[2]]
hm_map[hm_map>0] <- 1
hm_map[hm_map==0] <- 0
hm[,3:dim(hm)[2]] <- hm_map 
#print(dim(hm))
#print(head(hm))
png('0Values_Heatmap.png', width=1280, height=880)
print(Heatmap(hm[,3:dim(hm)[2]], name='0 values Repartition', cluster_columns=F, col=colorRampPalette(c("white", "grey", "black"))(50), row_dend_width=unit(30.5,'mm'),show_row_names=F, row_names_gp=gpar(fontsize=8.5)))
dev.off()

# Répartition de l'expression des miRNAs par ech et dans toute la cohorte #

#for (i in colnames(results[,3:dim(results)[2]])){
#	print(i)	
#	mat <- melt(data.frame('miRNA'=rownames(results), 'Value'=results[,i]))
#	mat <- mat[order(mat[,3], decreasing=T),]
#	perc <- vector()
#	val_sum <- sum(mat[,3])
#	mat <- cbind(mat, 'perc'=(mat$value/val_sum)*100)
#	mat <- mat[mat[,3]>0,]	
#	group <- vector()
#	for (j in 1:nrow(mat)){if(mat[j,4]>1){group <- c(group, mat[j,1])}else{group <- c(group, "Below1")}}
#	mat <- cbind(mat, 'Label'=group)
