#!/bin/bash

module load java; module load samtools


g=${SLURM_ARRAY_TASK_ID}

f=`head -$g 00_bam.list | tail -1`;


#samtools addreplacerg -r ID:$f -r PL:Illumina -r SM:$f bams/$f.sort.bam -o bams/rg.$f.bam

#samtools index bams/rg.$f.bam


/project2/lbarreiro/users/tauras/Programs/gatk-4.1.4.1/gatk HaplotypeCaller -ERC GVCF -I bams/$f.sort.bam  -R hg19/hg19.fa  -O gVCF/$f.g.vcf.gz -mbq 20

echo $f $g

