#poojasingh
#nov2020
#prepping files for FST calculation using poolfstat and running the calc
#this if for JP #Coadaptree


options(stringsAsFactors = FALSE)

# libaries
library(dplyr)
library(poolfstat)
library(reshape2)

# read in your snps

a <- read.table("/data/projects/pool_seq/pangenome/JP_pangenome/JP_pooled/snpsANDindels/03_maf-p05_RD-recalculated/JP_pooled-varscan_all_bedfiles_SNP_maf_RD-recalculated.txt", header=T, sep="\t")


# number of snps 
nsnp <- nrow(a)

# number of pools
b <- a  %>% select(contains(".FREQ"))
npools <- ncol(b)

# ref allele read count

refallele.readcount <- a %>% select(contains(".RD"))
colnames(refallele.readcount) <- gsub(".RD", "", colnames(refallele.readcount))
refallele.readcount <- data.matrix(refallele.readcount)

# total depth

readcoverage <- a %>% select(contains(".DP"))
colnames(readcoverage) <- gsub(".DP", "", colnames(readcoverage))
readcoverage <- data.matrix(readcoverage)

#snp chrom, pos, alt, ref info

snp.info <- a[,1:4]
snp.info <- data.matrix(snp.info)

#haploidy
c <- read.table("/data/projects/pool_seq/pangenome/JP_pangenome/JP_pooled/datatable.txt", header=T, sep="\t")
d <- c[,c(1,4)]
d <- unique(d)
poolsize <- d$ploidy/2


#names of pools

poolnames <- colnames(readcoverage)
poolnames <- as.character(poolnames)



## create S4 class vector for poolfstat: https://www.rdocumentation.org/packages/poolfstat/versions/1.2.0/topics/pooldata-class

setClass("pooldata", slots=list(nsnp="numeric", npools="numeric",refallele.readcount= "matrix", readcoverage="matrix", snp.info="matrix", poolsizes="vector", poolnames="vector"))

pooldata <- new("pooldata",nsnp=nsnp, npools=npools, refallele.readcount=refallele.readcount, readcoverage=readcoverage, snp.info=snp.info, poolsizes=poolsize, poolnames=poolnames)

#calcualte pairiwse FST

out <- computePairwiseFSTmatrix(pooldata, method = "Anova", min.maf =0.05, output.snp.values = FALSE)

out1 <- melt(out$PairwiseFSTmatrix)
colnames(out1) <- c("pop1", "pop2", "FST")
write.table(out1, "JP_pairwiseFST.txt", row.names=F, quote=F, sep="\t")

out2 <- melt(out$NbOfSNPs)
colnames(out1) <- c("pop1", "pop2", "SNPs")
write.table(out1, "JP_pairwiseFST_snps_used_for_calc.txt", row.names=F, quote=F, sep="\t")

#calcualte pairwise FST per SNP, apparenlty this takes FOREVER, so I recommend slicing your input data for the two pops that you are particularly interested in based on your hypothesis

out_snp <- computePairwiseFSTmatrix(pooldata, method = "Anova", min.maf =0.05, output.snp.values = TRUE)
