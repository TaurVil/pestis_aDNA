# pestis_aDNA
Identifying selection using changes in allele frequency during the black plague


# Starting point
Analyses in this repository begin with bams mapped by Jennifer Klunk. 


# Get uniquely mapping reads
Files were originally mapped to only the regions of the genome targetted for enrichment. This produced weird patterns of heterozygosity due to calling as genetic variants regions which would have actually mapped to another part of the genome. To correct this, we will remap to the entire genome, and then retain reads which mapped uniquely. 

```console
ls */*bam > 00_bam.list; sed -i 's/.min24.MQ30.merged.RG.bam//g' 00_bam.list; sed -i 's/_20.*//g' 00_bam.list; sed -i 's/Final//g' 00_bam.list; sed -i 's/_redesigned//g' 00_bam.list; sed -i 's/bams\///g' 00_bam.list; sed -i 's/_seqs//g' 00_bam.list

# sort original bam by queryname, and convert to fastq files
samtools sort -n bams/$f*min24*.bam -o $f.sortname.bam

# remove PE reads from SE dataset 
# convert fq to bam; merge the two bams
java -jar /project2/lbarreiro/users/tauras/Programs/picard.jar FastqToSam FASTQ=$f.R1.fq FASTQ2=$f.R2.fq OUTPUT=tmp5.umap_fq.bam READ_GROUP_NAME=LD05 SAMPLE_NAME=LD05 LIBRARY_NAME=ILL
java -jar /project2/lbarreiro/users/tauras/Programs/picard.jar FastqToSam FASTQ=$f.R.fq OUTPUT=tmp5.umap_fq.bam READ_GROUP_NAME=LD05 SAMPLE_NAME=LD05 LIBRARY_NAME=ILL
samtools merge
# map 
/project2/lbarreiro/users/tauras/Programs/bwa_bam2bam/network-aware-bwa-master/bwa bam2bam -g hg19/hg19.fa -l 16500 -f bams/$f.mapped.bam tmp.$f.umap_merged.bam; samtools sort -o bams/$f.sort.bam bams/$f.mapped.bam; samtools index bams/$f.sort.bam; rm bams/$f.mapped.bam; samtools flagstat bams/$f.sort.bam


sbatch --array=1-1003%225 --mem=8G --account=pi-lbarreiro --partition=lbarreiro run.01.bwa_map.sh
  for index in `seq 1 10`; do f=`head -$index 00_bam.list | tail -1`; samtools sort -n bams/$f*min24*.bam -o $f.sortname.bam; /project2/lbarreiro/users/tauras/Programs/bwa_bam2bam/network-aware-bwa-master/bwa bam2bam -g hg19/hg19.fa -l 16500 -f bams/$f.mapped.bam $f.sortname.bam; rm $f.sortname.bam; samtools sort -o bams/$f.sort.bam bams/$f.mapped.bam; samtools index bams/$f.sort.bam; rm bams/$f.mapped.bam; echo $index; done


bedtools bamtobed -i exons/LM16_ExonFinal_2018-0307.min24.MQ30.merged.RG.bam  | cut -f 4 > $name.
```


# 
