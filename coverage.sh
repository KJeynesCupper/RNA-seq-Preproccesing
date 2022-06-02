#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --ntasks 6 # request 8 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly

module load Python/3.8.2-GCCcore-9.3.0
module load HTSeq/0.13.5-foss-2020a-Python-3.8.2


FILES=data/3_sorted/*

for f in $FILES
do
f=${f##*/}
f=${f%_sorted.bam}


### Create a Bedgraph summaring mapped reads at each bp (there is not longer the need of the chromosome file)
echo create the Bedgraph file for $f...
bedtools genomecov -split -bg -ibam data/3_sorted/${f}_sorted.bam \
-g reference/chomosome_size/chromosome_size_reference.txt > data/coverage/${f}.bedgraph

### Estimate average reads across the genome
echo estimate average reads across genome for $f...
t=data/3_sorted/${f}_sorted.bam;
c=$(samtools depth $t | awk -v var=$g '{sum+=$3;cnt++}END{print sum/var}'); #-> this is an average read per bp
echo -e ${f}\t${c} >> data/coverage/genome_coverage.txt

### Use a simple Perl command to normalise read counts by genome-wide average counts (simple coverage and log2)
echo normalise read counts by genome-wide average counts for $f...
perl -ne 'chomp($_); @a=split(/\t/,$_);print $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]/'$c'."\t"."\n";' data/coverage/${f}.bedgraph > data/coverage/${f}_norm.bedgraph

done
