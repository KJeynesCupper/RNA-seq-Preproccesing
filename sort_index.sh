#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --ntasks 6 # request 6 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly

module purge;
module load bluebear
module load Python/2.7.16-GCCcore-8.3.0
module load SAMtools/1.10-GCC-8.3.0
# module load SAMtools/1.12-GCC-10.3.0


FILES="data/2_mapped/"*".bam"

for f in $FILES
do
f=${f##*/}
f=${f%.bam}

#sort
echo sorting .bam file for  $f...
samtools sort --threads 4 -o "data/3_sorted/$f""_sorted.bam" "data/2_mapped/$f"".bam"

# indexing
echo indexing $f...
samtools index "data/3_sorted/$f""_sorted.bam"

done
