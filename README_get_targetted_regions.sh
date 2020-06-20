## Get targetted exons: Neutral_redesigned_2016-0415.fasta to 30.neutral.bed 
# The shortest target is 180bp, so we'll require the match be at least 175 bp while using blat, and then filter for unique & best matching later. 
module load blat
blat ./hg19/hg19.fa /project2/lbarreiro/users/tauras/pestis_aDNA/Neutral/Neutral_redesigned_2016-0415.fasta neutral.psl -minScore=175 -minIdentity=98 -fastMap
blat ./hg19/hg19.fa /project2/lbarreiro/users/tauras/pestis_aDNA/GWAS/Immune_seqs_2016-0912.fasta immune.psl -minScore=175 -minIdentity=98 -fastMap
blat ./hg19/hg19.fa /project2/lbarreiro/users/tauras/pestis_aDNA/Exon/ExonFinal_2018-0307.fasta exon.psl -minScore=175 -minIdentity=98 

## Neutral: 250 loci all of which map uniquely and completely
## Immune: 494 loci all of which map uniquely and completely
## Exon: 3556 loci, we keep 3433. CCL3L3CCL3L1CCL4L1CCL4L2-1-chr17-31648783-31649098 is duplicated with full mapping. there are 10 where there's an incomplete mapping. Some of these are for chrY, which I expect to be losing. 
# write.table(d4[,c(14,16,17)], "exon.bed", row.names=F, col.names=F, sep="\t", quote=F); write.table(d[,c(14,16,17)], "neutral.bed", row.names=F, col.names=F, sep="\t", quote=F); write.table(d[,c(14,16,17)], "immune.bed", row.names=F, col.names=F, sep="\t", quote=F)


## source file locations
/project2/lbarreiro/users/tauras/pestis_aDNA/Exon/ExonFinal_2018-0307.fasta
/project2/lbarreiro/users/tauras/pestis_aDNA/GWAS/Immune_seqs_2016-0912.fasta
/project2/lbarreiro/users/tauras/pestis_aDNA/Neutral/Neutral_redesigned_2016-0415.fasta





## If I wanted to post-process more in R and get the gene name by the exon
library(data.table); library(tidyverse); fread("hg19_exons.psl") -> data
colnames(data)[1] <- "V1"
data$rat <- data$V1/data$V11
for (i in 1:nrow(data)) {data$chrom_match[i] <- str_detect(pattern=data$V14[i], string=data$V10[i]); print(i)}
subset(data, data$chrom_match == T & data$rat > 0.9) -> d2

dim(d2); length(unique(d2$V10)) #There's some number of multiply mapped sections, but we'll ignore those for now. There are only 6... 
write.table(d2[,c(14, 16, 17)], "./30.regions_hg19.bed", sep="\t", quote=F, col.names=F, row.names=F) # output bed file 
library(plyr); strsplit(data$V10, "-") -> l1; ldply(l1[][1:4]) -> l1
l1$V3 -> data$init_chrom
new <- NULL; for (i in 1:length(l1)) { rbind(new, data.frame(matrix(nrow=1,ncol=5,l1[[i]]))[,1:5]) -> new; print(i)}
