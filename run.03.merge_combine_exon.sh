#!/bin/bash

index=${SLURM_ARRAY_TASK_ID}

n=`head -$index 01_chroms | tail -1`
#head -$n genome.35Mb.bed | tail -1 > tmp.$n.bed


module load java; module load samtools; module load python
module load htslib


/project2/lbarreiro/users/tauras/Programs/gatk-4.1.4.1/gatk CombineGVCFs -R hg19/hg19.fa -O $n.exon.g.vcf.gz -V 01_exon_gvcf.list -L $n

/project2/lbarreiro/users/tauras/Programs/gatk-4.1.4.1/gatk GenotypeGVCFs -V $n.exon.g.vcf.gz -R hg19/hg19.fa -O $n.exon.vcf.gz --tmp-dir /scratch/midway2/rdd6247/

