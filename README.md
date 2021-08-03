# pestis_aDNA
Identifying selection using changes in allele frequency during the black plague


# Starting point
Analyses in this repository begin with bams mapped by Jennifer Klunk. 


# Get uniquely mapping reads
Files were originally mapped to only the regions of the genome targetted for enrichment. This produced weird patterns of heterozygosity due to calling as genetic variants regions which would have actually mapped to another part of the genome. To correct this, we remap reads to the entire genome, and then retain reads which mapped uniquely within the target region. 

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
/project2/lbarreiro/users/tauras/Programs/bwa_bam2bam/network-aware-bwa-master/bwa bam2bam -g hg19/hg19.fa -l 16500 -f bams/$f.mapped.bam tmp6.$f.tomap.bam; samtools sort -o bams/$f.sort.bam bams/$f.mapped.bam; samtools index bams/$f.sort.bam; rm bams/$f.mapped.bam; samtools flagstat bams/$f.sort.bam; rm tmp*$f.*bam; rm $f.*fq; rm $f.names; rm $f.sortname.bam; done

sbatch --array=1-1003%225 --mem=8G --account=pi-lbarreiro --partition=lbarreiro run.01.bwa_map.sh

# time /project2/lbarreiro/users/tauras/Programs/bwa_bam2bam/network-aware-bwa-master/bwa bam2bam -g hg19/hg19.fa -n 0.01 -o 2 -l 16500 -f test.out bams/cut.LP372a_S252_L001.leehom.bam
## Let's remove "-n" and "-o" options which cuts down the time to ~1/5
				
```

# Get genotype calls
Now we need to call genotypes for each dataset, starting with the gVCF files then assembling the combined call sets. 
```console
mkdir gVCF
## make gVCF file for each individual
	sbatch --array=1-250%80 --mem=16G --account=pi-lbarreiro --partition=lbarreiro run.02.gVCF.sh
	sbatch --array=501-750%100 --mem=16G run.02.gVCF.sh
  
ls gVCF/*Neutral*gz > 105_neutral_cohortmap; awk 'BEGIN { FS="/t"; OFS="/t" } { print $1 $1}'  105_neutral_cohortmap > tmp4 ; sed 's/gzgVCF/gz \t gVCF/g' tmp4 > tmp2; sed 's/ gVCF/gVCF/g' tmp2 > 105_neutral_cohortmap; sed -i 's/^gVCF\///g' 105_neutral_cohortmap; sed -i 's/.g.vcf.gz //g' 105_neutral_cohortmap; rm tmp2; rm tmp4
ls gVCF/*Immune*gz > 105_immune_cohortmap; awk 'BEGIN { FS="/t"; OFS="/t" } { print $1 $1}'  105_immune_cohortmap > tmp4 ; sed 's/gzgVCF/gz \t gVCF/g' tmp4 > tmp2; sed 's/ gVCF/gVCF/g' tmp2 > 105_immune_cohortmap; sed -i 's/^gVCF\///g' 105_immune_cohortmap; sed -i 's/.g.vcf.gz //g' 105_immune_cohortmap; rm tmp2; rm tmp4
ls gVCF/*Exon*gz > 105_exon_cohortmap; awk 'BEGIN { FS="/t"; OFS="/t" } { print $1 $1}'  105_exon_cohortmap > tmp4 ; sed 's/gzgVCF/gz \t gVCF/g' tmp4 > tmp2; sed 's/ gVCF/gVCF/g' tmp2 > 105_exon_cohortmap; sed -i 's/^gVCF\///g' 105_exon_cohortmap; sed -i 's/.g.vcf.gz //g' 105_exon_cohortmap; rm tmp2; rm tmp4

ls gVCF/*Exon*gz > 01_exon_gvcf.list; sbatch --array=1-22 --mem=16G --account=pi-lbarreiro --partition=lbarreiro run.03.merge_combine_exon.sh
ls gVCF/*Neutral*gz > 01_neutral_gvcf.list; sbatch --array=1-22 --mem=16G run.03.merge_combine_neutral.sh
ls gVCF/*Immune*gz > 01_immune_gvcf.list; sbatch --array=1-22 --mem=16G run.03.merge_combine_immune.sh

```

# Filter genotype calls for target regions
```console
module load bcftools
bcftools concat chr*neutral.vcf.gz -a > neutral.vcf.gz
bcftools concat chr*exon.vcf.gz -a > exon.vcf.gz
bcftools concat chr*immune.vcf.gz -a > immune.vcf.gz
vcftools --gzvcf neutral.vcf.gz --mac 3 --max-alleles 2 --bed ./neutral.bed --minQ 30 --recode --out neutral_mac3_biallelic
vcftools --gzvcf exon.vcf.gz --mac 3 --max-alleles 2 --bed ./exon.bed --minQ 30 --recode --out exon_mac3_biallelic
vcftools --gzvcf immune.vcf.gz --mac 3 --max-alleles 2 --bed ./immune.bed --minQ 30 --recode --out immune_mac3_biallelic

sed -i 's/_Immune//g' immune_mac3_biallelic.recode.vcf
sed -i 's/_Exon//g' exon_mac3_biallelic.recode.vcf
sed -i 's/_Neutral//g' neutral_mac3_biallelic.recode.vcf

bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL immune_mac3_biallelic.recode.vcf | bcftools +setGT -- -ta -nu > joint.immune.forgenolik.vcf
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL exon_mac3_biallelic.recode.vcf | bcftools +setGT -- -ta -nu > joint.exon.forgenolik.vcf
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL neutral_mac3_biallelic.recode.vcf | bcftools +setGT -- -ta -nu > joint.neutral.forgenolik.vcf
```

# Get genotype likelihoods for each subset
```console
module load vcftools
## Exons
vcftools --vcf joint.exon.forgenolik.vcf --keep 05.pre_exon_london.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_pre_exons
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 38 > genolik.exons_london_pre.genolik
vcftools --vcf joint.exon.forgenolik.vcf --keep 05.post_exon_london.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_post_exons
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 63 > genolik.exons_london_post.genolik
vcftools --vcf joint.exon.forgenolik.vcf --keep 05.BD_exon_london.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_during_exons
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 41 > genolik.exons_london_during.genolik

	vcftools --vcf joint.exon.forgenolik.vcf --keep 05.pre_exon_denmark.txt --recode --out temp2; vcftools --vcf temp2.recode.vcf --012 --out denmark_pre_exons; sed '/^#/d' temp2.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../../Programs/LCLAE/filtbaboon1b 42 > genolik.exons_denmark_pre.genolik
	vcftools --vcf joint.exon.forgenolik.vcf --keep 05.post_exon_denmark.txt --recode --out temp2; vcftools --vcf temp2.recode.vcf --012 --out denmark_post_exons
	sed '/^#/d' temp2.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../../Programs/LCLAE/filtbaboon1b 58 > genolik.exons_denmark_post.genolik
	vcftools --vcf joint.exon.forgenolik.vcf --keep 05.BD_exon_denmark.txt --recode --out temp2; vcftools --vcf temp2.recode.vcf --012 --out  denmark_during_exons
	sed '/^#/d' temp2.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../../Programs/LCLAE/filtbaboon1b 24 > genolik.exons_denmark_during.genolik


## Immune
vcftools --vcf joint.immune.forgenolik.vcf --keep 00_London_pre_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_pre_immune; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 65 > genolik.gwas_london_pre.genolik
vcftools --vcf joint.immune.forgenolik.vcf --keep 00_London_post_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_post_immune; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 98 > genolik.gwas_london_post.genolik
vcftools --vcf joint.immune.forgenolik.vcf --keep 00_London_during_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_during_immune; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 64 > genolik.gwas_london_during.genolik

vcftools --vcf joint.immune.forgenolik.vcf --keep 00_Denmark_pre_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_pre_immune
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 42 > genolik.gwas_denmark_pre.genolik

vcftools --vcf joint.immune.forgenolik.vcf --keep 00_Denmark_post_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_post_immune
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 57 > genolik.gwas_denmark_post.genolik

vcftools --vcf joint.immune.forgenolik.vcf --keep 00_Denmark_during_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_during_immune
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 24 > genolik.gwas_denmark_during.genolik

## Immune
vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_London_pre_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_pre_neutral; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 65 > genolik.neutral_london_pre.genolik
vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_London_post_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_post_neutral; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 100 > genolik.neutral_london_post.genolik
vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_London_during_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_during_neutral; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 63 > genolik.neutral_london_during.genolik

vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_Denmark_pre_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_pre_neutral
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 42 > genolik.neutral_denmark_pre.genolik

vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_Denmark_post_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_post_neutral
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 57 > genolik.neutral_denmark_post.genolik

vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_Denmark_during_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_during_neutral
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | /project2/lbarreiro/users/tauras/Programs/LCLAE/filtbaboon1b 24 > genolik.neutral_denmark_during.genolik

```
Exons, London: pre=38, post=63, BD=41. Denmark: pre=42, post=58, BD=24 (n=264)

Immune & Neutral, London: pre=65, post=98 & 100, BD=64 & 63. Denmark: pre=42, post=57, BD=24 (n=350 & 352)

These results feed into aDNA_enrichment_analysis1.Rmd (which is stored locally)


```console
# Get list of sites that are human variants
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz
vcftools --gzvcf All_20180418.vcf.gz --remove-indels --kept-sites --out dbsnp.all
cat dbsnp.all.kept.sites | awk -v OFS='\t' '{print $1, $2, $2}' > dbsnp.all.bed

# Get list of sites known to be common human variants
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/common_all_20180418.vcf.gz
vcftools --gzvcf common_all_20180418.vcf.gz --remove-indels --kept-sites --out dbsnp.common
cat dbsnp.common.kept.sites | awk -v OFS='\t' '{print $1, $2, $2}' > dbsnp.common.bed


## pull out the variants I called which overlap common or all variants
module load vcftools
vcftools --gzvcf neutral.vcf.gz --mac 3 --max-alleles 2 --bed ./neutral.bed --minQ 30 --kept-sites --out int_dbsnp.neutral 
vcftools --gzvcf exon.vcf.gz --mac 3 --max-alleles 2 --bed ./exon.bed --minQ 30 --kept-sites --out int_dbsnp.exon 
vcftools --gzvcf immune.vcf.gz --mac 3 --max-alleles 2 --bed ./immune.bed --minQ 30 --kept-sites --out int_dbsnp.immune

cat int_dbsnp.immune.kept.sites | awk -v OFS='\t' '{print $1, $2, $2}' > int_dbsnp.immune.bed
cat int_dbsnp.exon.kept.sites | awk -v OFS='\t' '{print $1, $2, $2}' > int_dbsnp.exon.bed
cat int_dbsnp.neutral.kept.sites | awk -v OFS='\t' '{print $1, $2, $2}' > int_dbsnp.neutral.bed

module load bedtools
bedtools intersect -u -a int_dbsnp.immune.bed -b dbsnp.common.bed > common.immune.bed
bedtools intersect -u -a int_dbsnp.exon.bed -b dbsnp.common.bed > common.exon.bed
bedtools intersect -u -a int_dbsnp.neutral.bed -b dbsnp.common.bed > common.neutral.bed

for chr in `cat 01_chroms`; do grep -P "${chr}\t" dbsnp.all.bed > tmp.bed
	bedtools intersect -u -a int_dbsnp.immune.bed -b tmp.bed >> all.immune.bed; 
	bedtools intersect -u -a int_dbsnp.exon.bed -b tmp.bed >> all.exon.bed; 
	bedtools intersect -u -a int_dbsnp.neutral.bed -b tmp.bed >> all.neutral.bed; done

# VQSR?
# Hard filters?
/project2/lbarreiro/users/tauras/Programs/gatk-4.1.4.1/gatk VariantFiltration -V refilt.immune_mac3_biallelic.recode.vcf -R hg19/hg19.fa -O refilt.immune_gatk.vcf.gz -filter-ame "FS" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

```
