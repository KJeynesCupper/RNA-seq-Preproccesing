#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --ntasks 6 # request 6 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly



FILES=1_raw/*_1.fq
for f in $FILES
do
f=${f##*/}
f=${f%_1.fq}

gunzip $f # First uncompress your forward and reverse reads:

sortmerna --ref $SORTMERNA_DB --reads $f --aligned read-rRNA-hits --other read-sortmerna --log -a 16 -v --fastx


done
