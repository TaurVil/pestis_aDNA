## 1k genomes data for finland and gbr to compare aDNA with modern individuals

## Get thousand genomes data 
for chr in `seq 1 22`; do wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr$chr'.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'; done

## Combine individuals
module load bcftools; bcftools concat -O z -o 1kgenomes_phased.vcf.gz ALL.chr*.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

## List of samples from https://www.internationalgenome.org/data-portal/population/FIN and https://www.internationalgenome.org/data-portal/population/GBR
## Use cut -f 1 to get just the names

## Need vcf to be in verison 4.2 to use vcftools installed on midway2
module load htslib; zcat 1kgenomes_phased.vcf.gz | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' | bgzip > 1kgenomes_v4.2.vcf.gz

## Get target regions 
cat ../remapped_data/exon.bed ../remapped_data/neutral.bed ../remapped_data/immune.bed | sort > targets.bed; sed -i s/chr//g targets.bed

## Calculate allele frequencies for each population
module load vcftools; vcftools --gzvcf 1kgenomes_v4.2.vcf.gz --keep gbr_names.list --freq --out gbr --bed targets.bed
module load vcftools; vcftools --gzvcf 1kgenomes_v4.2.vcf.gz --keep fin_names.list --freq --out fin --bed targets.bed

for f in `ls *frq`; do sed -i '1d' $f; sed -i 's/:/\t/g' $f; sed -i 's/^/chr/g' $f; done 

## Transport fin.frq & gbr.frq to local directory & integrate into analyses
