#poojasingh
#coadaptree 2020
#r2 LD calculation per scaffold

library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
options(stringsAsFactors = FALSE)
setwd(".")



### read in SNPs

data1 <- read.table("/data/projects/pool_seq/DF_datasets/DF_pooled_GEA/DF_pooled/snpsANDindels/coastal_variety/03_maf-p05_RD-recalculated_FDC/DF_pooled-varscan_all_bedfiles_SNP_FDC_maf_RD-recalculated.txt", header=T, sep="\t")



### select columns with FREQ and order the populations and conver FREQ % to decimal

data2 <- data1  %>% select(contains(".FREQ"))
colhead <- sub(".FREQ", "", colnames(data2))
colnames(data2) <- colhead
data3 <- data2[ , order(colnames(data2))]
snps <- data3 %>% mutate_each(funs(as.numeric(gsub("%", "", ., fixed = TRUE))/100))
snps$CHROM <- data1$CHROM
snps$POS <- data1$POS


### get unique set of contigs/scaffolds/LGs to loop through

chroms <- unique(data1$CHROM)

### correlation of pairwise snps per scaffold ###
args <- commandArgs(trailingOnly = TRUE) 


#list of results per scaffold

Res<-list()

#set system time to calculate how long the loop takes

system.time(

## loop over each scaffold

for(j in 1:length(chroms)){
		out <- snps[snps$CHROM == chroms[j],]
		cormat <- cor(t(out[,1:37]), method="spearman", use = "pairwise.complete.obs") #1:37 picks the columns with the population FREQ data and ignores the chrom and pos columns this should be changed based on your input number of populations
		cormat1 <- as.data.frame(cormat)
		colnames(cormat1) <- out$POS
		cormat1$snp1 <- out$POS
		result <- gather(data = cormat1, key = "snp2", value = "correlation", -snp1)
		result$chrom <- chroms[j]
		Res[[j]]<- result
	
		}

) 
	
Res1 <- do.call(rbind, Res)


write.table(Res1, "snp_pairwise_maf_correlation.txt", sep="\t", col.names = T, row.names = F, quote = F)
