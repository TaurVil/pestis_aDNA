#!/bin/bash
#SBATCH --get-user-env

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/project2/lbarreiro/users/tauras/Programs/zeromq/libzmq-master/build
export LIBRARY_PATH=$LIBRARY_PATH:/project2/lbarreiro/users/tauras/Programs/zeromq/libzmq-master/build/lib:/project2/lbarreiro/users/tauras/Programs/luarocks/lib/lua/5.3/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project2/lbarreiro/users/tauras/Programs/zeromq/libzmq-master/build/lib:/project2/lbarreiro/users/tauras/Programs/luarocks/lib/lua/5.3/
export CPATH=$CPATH:/project2/lbarreiro/users/tauras/Programs/zeromq/libzmq-master/include/

module load gcc; module load cmake/3.6.2; module load lua; module load bwa
module load samtools


index=${SLURM_ARRAY_TASK_ID}
f=`head -$index 00_bam.list | tail -1`

samtools sort -n bams/$f*min24*.bam -o $f.sortname.bam

/project2/lbarreiro/users/tauras/Programs/bwa_bam2bam/network-aware-bwa-master/bwa bam2bam -g hg19/hg19.fa -l 16500 -f bams/$f.mapped.bam $f.sortname.bam

rm $f.sortname.bam

samtools sort -o bams/$f.sort.bam bams/$f.mapped.bam; samtools index bams/$f.sort.bam

rm bams/$f.mapped.bam

samtools flagstat bams/$f.sort.bam

