#!/bin/bash
#SBATCH --qos bbdefault
#SBATCH --ntasks 6 # request 8 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly

module purge;
module load bluebear
module load fastp/0.20.1-GCC-8.3.0

FILES=raw/*_1.fq
for f in $FILES
do
f=${f##*/}
f=${f%_1.fq}


fastp --thread 4 \
-i 1_raw/${f}_1.fq \
-I 1_raw/${f}_2.fq \
-o 2_trimmed/${f}_1_trimmo.fq.gz \
-O 2_trimmed/${f}_2_trimmo.fq.gz

# check overall quality with fastqc after trimming
fastqc 2_trimmed/${f}_1_trimmo.fq.gz -o qc/2_trim_stat/
fastqc 2_trimmed/${f}_2_trimmo.fq.gz -o qc/2_trim_stat/

multiqc 2_trimmed/${f}_1_trimmo.fq.gz qc/2_trim_stat/
multiqc 2_trimmed/${f}_2_trimmo.fq.gz qc/2_trim_stat/

done
