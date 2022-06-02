#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --ntasks 6 # request 8 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly

module load Python/3.8.2-GCCcore-9.3.0
module load HTSeq/0.13.5-foss-2020a-Python-3.8.2
#module load Python/3.8.6-GCCcore-10.2.0
#module load HTSeq/0.12.4-foss-2019b-Python-3.7.4

FILES=data/3_sorted/*_sorted.bam

for f in $FILES
do
f=${f##*/}
f=${f%_sorted*}

### raw count with HTSeq
#Genes tomato
python -m HTSeq.scripts.count \
--order=pos --stranded=no --mode=union --nonunique=none --type=gene --idattr=Name \
data/3_sorted/${f}_sorted.bam \
reference/1_raw_ref/*.gff \
> data/4_counts/${f}_counts.txt

done
