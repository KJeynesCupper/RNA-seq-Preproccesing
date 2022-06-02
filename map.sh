#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --ntasks 6 # request 6 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly

module purge;
module load bluebear
module load SAMtools/1.12-GCC-10.3.0
module load HISAT2/2.2.1-foss-2019b

REF_INDEX=reference/hisat_index/*
FILES=data/2_trimmed/*1_trimmo.fq.gz
MAPPED_FILES=data/3_mapped/


for f in $FILES
do

f=${f##*/}
f=${f%_1_trimmo*}

hisat2 -p 6 \
-x $REF_INDEX \
--summary-file $MAPPED_FILES/${f}_summary.txt \
-1 $FILES/${f}_1.fastq.gz \
-2 $FILES/${f}_2.fastq.gz | samtools view -bS -o $MAPPED_FILES/${f}.bam
done;
