# pestis_aDNA
Identifying selection using changes in allele frequency during the black plague


# Starting point
Analyses in this repository begin with bams mapped by Jennifer Klunk. 


# Get uniquely mapping reads
Files were originally mapped to only the regions of the genome targetted for enrichment. This produced weird patterns of heterozygosity due to calling as genetic variants regions which would have actually mapped to another part of the genome. To correct this, we will remap to the entire genome, and then retain reads which mapped uniquely within the target region. 

```console
# Get list of individuals and genomes
ls */*bam > 00_bam.list; sed -i 's/.min24.MQ30.merged.RG.bam//g' 00_bam.list; sed -i 's/_20.*//g' 00_bam.list; sed -i 's/Final//g' 00_bam.list; sed -i 's/_redesigned//g' 00_bam.list; sed -i 's/bams\///g' 00_bam.list; sed -i 's/_seqs//g' 00_bam.list

module load samtools; module load bedtools; module load java
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/project2/lbarreiro/users/tauras/Programs/zeromq/libzmq-master/build
export LIBRARY_PATH=$LIBRARY_PATH:/project2/lbarreiro/users/tauras/Programs/zeromq/libzmq-master/build/lib:/project2/lbarreiro/users/tauras/Programs/luarocks/lib/lua/5.3/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project2/lbarreiro/users/tauras/Programs/zeromq/libzmq-master/build/lib:/project2/lbarreiro/users/tauras/Programs/luarocks/lib/lua/5.3/
export CPATH=$CPATH:/project2/lbarreiro/users/tauras/Programs/zeromq/libzmq-master/include/
module load gcc; module load cmake/3.6.2; module load lua; module load bwa


for index in `seq 951 1003`; do f=`head -$index 00_bam.list | tail -1`;
# sort original bam by queryname, and convert to fastq files
samtools sort -n bams/$f*min24*.bam -o $f.sortname.bam
bedtools bamtofastq -i $f.sortname.bam -fq $f.R.fq
bedtools bamtofastq -i $f.sortname.bam -fq $f.R1.fq -fq2 $f.R2.fq
# remove PE reads from SE dataset 
awk 'NR%4==1 {print substr($1,2)}' $f.R1.fq | sed 's/\/1//g'> $f.names #extract read name from fastq
/project2/lbarreiro/users/tauras/Programs/bbmap/filterbyname.sh in=$f.R.fq names=$f.names out=$f.remaining.fq exclude qin=33 overwrite=T
# convert fq to bam; merge the two bams
java -jar /project2/lbarreiro/users/tauras/Programs/picard.jar FastqToSam FASTQ=$f.R1.fq FASTQ2=$f.R2.fq OUTPUT=tmp5.$f.umap_fq.bam READ_GROUP_NAME=$f SAMPLE_NAME=$f LIBRARY_NAME=ILL
java -jar /project2/lbarreiro/users/tauras/Programs/picard.jar FastqToSam FASTQ=$f.remaining.fq OUTPUT=tmp5.$f.umap_fq2.bam READ_GROUP_NAME=$f SAMPLE_NAME=$f LIBRARY_NAME=ILL
samtools merge tmp6.$f.tomap.bam tmp5.$f.umap_fq2.bam tmp5.$f.umap_fq.bam
# map 
/project2/lbarreiro/users/tauras/Programs/bwa_bam2bam/network-aware-bwa-master/bwa bam2bam -g hg19/hg19.fa -l 16500 -f bams/$f.mapped.bam tmp6.$f.tomap.bam; samtools sort -o bams/$f.sort.bam bams/$f.mapped.bam; samtools index bams/$f.sort.bam; rm bams/$f.mapped.bam; samtools flagstat bams/$f.sort.bam
rm tmp*$f.*bam; rm $f.*fq; rm $f.names; rm $f.sortname.bam; done




sbatch --array=1-1003%225 --mem=8G --account=pi-lbarreiro --partition=lbarreiro run.01.bwa_map.sh
  for index in `seq 1 10`; do f=`head -$index 00_bam.list | tail -1`; samtools sort -n bams/$f*min24*.bam -o $f.sortname.bam; /project2/lbarreiro/users/tauras/Programs/bwa_bam2bam/network-aware-bwa-master/bwa bam2bam -g hg19/hg19.fa -l 16500 -f bams/$f.mapped.bam $f.sortname.bam; rm $f.sortname.bam; samtools sort -o bams/$f.sort.bam bams/$f.mapped.bam; samtools index bams/$f.sort.bam; rm bams/$f.mapped.bam; echo $index; done


bedtools bamtobed -i exons/LM16_ExonFinal_2018-0307.min24.MQ30.merged.RG.bam  | cut -f 4 > $name.
```

Now we need to call genotypes for each dataset, starting with the gVCF files then assembling the combined call sets. 
```console
mkdir gVCF
	sbatch --array=1-250%80 --mem=16G --account=pi-lbarreiro --partition=lbarreiro run.02.gVCF.sh
	sbatch --array=501-750%100 --mem=16G run.02.gVCF.sh
  

ls gVCF/*Neutral*gz > 105_neutral_cohortmap; awk 'BEGIN { FS="/t"; OFS="/t" } { print $1 $1}'  105_neutral_cohortmap > tmp4 ; sed 's/gzgVCF/gz \t gVCF/g' tmp4 > tmp2; sed 's/ gVCF/gVCF/g' tmp2 > 105_neutral_cohortmap; sed -i 's/^gVCF\///g' 105_neutral_cohortmap; sed -i 's/.g.vcf.gz //g' 105_neutral_cohortmap; rm tmp2; rm tmp4
ls gVCF/*Immune*gz > 105_immune_cohortmap; awk 'BEGIN { FS="/t"; OFS="/t" } { print $1 $1}'  105_immune_cohortmap > tmp4 ; sed 's/gzgVCF/gz \t gVCF/g' tmp4 > tmp2; sed 's/ gVCF/gVCF/g' tmp2 > 105_immune_cohortmap; sed -i 's/^gVCF\///g' 105_immune_cohortmap; sed -i 's/.g.vcf.gz //g' 105_immune_cohortmap; rm tmp2; rm tmp4
ls gVCF/*Exon*gz > 105_exon_cohortmap; awk 'BEGIN { FS="/t"; OFS="/t" } { print $1 $1}'  105_exon_cohortmap > tmp4 ; sed 's/gzgVCF/gz \t gVCF/g' tmp4 > tmp2; sed 's/ gVCF/gVCF/g' tmp2 > 105_exon_cohortmap; sed -i 's/^gVCF\///g' 105_exon_cohortmap; sed -i 's/.g.vcf.gz //g' 105_exon_cohortmap; rm tmp2; rm tmp4

ls gVCF/*Exon*gz > 01_exon_gvcf.list
sbatch --array=1-22 --mem=16G --account=pi-lbarreiro --partition=lbarreiro run.03.merge_combine_exon.sh
ls gVCF/*Neutral*gz > 01_neutral_gvcf.list
sbatch --array=1-22 --mem=16G run.03.merge_combine_neutral.sh
ls gVCF/*Immune*gz > 01_immune_gvcf.list
sbatch --array=1-22 --mem=16G run.03.merge_combine_immune.sh

```

# 
