### Convert fastq to bams 
module load gcc; module load cmake/3.6.2; module load lua; module load bwa  
	for f in `cat 01_singleend_names.txt`; do /project2/lbarreiro/users/tauras/Programs/BCL2BAM2FASTQ/BCL2BAM2FASTQ-master/fastq2bam/fastq2bam -fq1 $f.fastq.gz -o $f.bam; done 
	for f in `head -1000 01_pairedend_names.txt`; do ls $f* | grep -v 'bam' > tmp; g=`head -1 tmp`; h=`tail -1 tmp`; /project2/lbarreiro/users/tauras/Programs/BCL2BAM2FASTQ/BCL2BAM2FASTQ-master/fastq2bam/fastq2bam -fq1 $g -fq2 $h -o $f.bam; done 
	for f in `tail -1002 01_pairedend_names.txt`; do ls $f* | grep -v 'bam' > tmp2; g=`head -1 tmp2`; h=`tail -1 tmp2`; /project2/lbarreiro/users/tauras/Programs/BCL2BAM2FASTQ/BCL2BAM2FASTQ-master/fastq2bam/fastq2bam -fq1 $g -fq2 $h -o $f.bam; done 


### Trim with leehom  
module load gcc; module load cmake/3.6.2
	for f in `head -200 00_all_names.txt`; do echo $f; /project2/lbarreiro/users/tauras/Programs/LeeHom/leeHom-master/src/leeHom -o $f.leehom.bam  --ancientdna $f.bam; done 
	for f in `head -2066 00_all_names.txt | tail -54`; do echo $f; /project2/lbarreiro/users/tauras/Programs/LeeHom/leeHom-master/src/leeHom -o $f.leehom.bam  --ancientdna $f.bam; done


### map with BWA 
# /project2/lbarreiro/users/tauras/Programs/bwa_bam2bam/network-aware-bwa-master/bwa index hg19/hg19.fa ## make a new index 
ls bams > 02_bams.txt; sed -i 's/.leehom.bam//g' 02_bams.txt 

sbatch --array=1-1003%225 --mem=8G --account=pi-lbarreiro --partition=lbarreiro run.01.bwa_map.sh
				

### Merge and dedup files from same library
	for index in `seq 101 200`; do f=`head -$index 102_uniq_libs | tail -1`; echo $index $f; samtools merge -f bams/merge1.$f.bam bams/*$f*.sort.bam; samtools index bams/merge1.$f.bam; samtools collate bams/merge1.$f.bam -o t2.sort.$f.bam; samtools fixmate -m t2.sort.$f.bam t2.fix.$f.bam; samtools sort t2.fix.$f.bam -o t2.sort.$f.bam; samtools markdup -r t2.sort.$f.bam bams/nodup.$f.bam; samtools index bams/nodup.$f.bam; samtools flagstat bams/nodup.$f.bam; done
	
	sbatch --array=700-716%80 --mem=6G --account=pi-lbarreiro --partition=lbarreiro run.02.merge_dedup.sh
	
	
	for index in  `seq 671 690`; do   f=`head -$index 13_unique_names.txt | tail -1`; samtools merge -f bams/merge1.$f.bam bams/*$f*.sort.bq.bam; samtools index bams/merge1.$f.bam; echo $index; done

### Merge and addRG from each individual 
	cp 102_uniq_libs 103_uniq_indivs; sed -i 's/cut.//g' 103_uniq_indivs; sed -i 's/[a-b,R,E]*$//g' 103_uniq_indivs; sed -i 's/[a-c]M3$//g' 103_uniq_indivs; uniq 103_uniq_indivs > tmp2; mv tmp2 103_uniq_indivs
	
	# 293 individuals 
	
	for g in `seq 1 293`; do f=`head -$g 103_uniq_indivs | tail -1`;  echo $g $f ; samtools merge -f bams/merge2.$f.bam bams/nodup.$f[a-c,E,R]*bam; samtools addreplacerg -r ID:$f -r PL:Illumina -r SM:$f bams/merge2.$f.bam -o bams/rg.$f.bam; samtools index bams/rg.$f.bam; echo $g; done 

### make gVCF 
mkdir gVCF ; module load java; module load samtools
sbatch --array=231-250 --mem=10G run.04.gVCF.sh
	sbatch --array=251-260%80 --mem=16G --account=pi-lbarreiro --partition=lbarreiro run.04.gVCF.sh
	# up to 291 

