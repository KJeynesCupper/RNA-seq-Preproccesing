#!/bin/bash

module purge;
module load bluebear

################################################################################
#  Preprocessing for RNAseq raw reads (paired)
################################################################################

# create a workplace ie. your working directory.
###### in workplace
mkdir -p data
mkdir -p data/1_raw
mkdir -p data/2_trimmed
mkdir -p data/2_mapped
mkdir -p data/3_sorted
mkdir -p data/4_counts

mkdir -p reference
mkdir -p reference/1_raw_ref
mkdir -p reference/2_index
mkdir -p reference/4_splice
mkdir -p reference/hista_index
mkdir -p reference/chromosome_size


mkdir -p qc
mkdir -p qc/1_raw_read_stats
mkdir -p qc/2_trim_stat

mkdir -p data/coverage

set -e  # this ensures that the script fails if any of the commands result in an error


################################################################
# change file names
################################################################

# mv name_fastq.gz  new_name_1.fastq.gz
# mv name_fastq.gz  new_name_2.fastq.gz

################################################################
# download reference genome and annoatation gff file
################################################################
module load HISAT2/2.2.1-foss-2019b
module load Python/3.8.6-GCCcore-10.2.0

cd reference

#Tomato - check for updated versions.
wget ftp://ftp.solgenomics.net/tomato_genome/Heinz1706/annotation/ITAG4.1_release/ITAG4.1_gene_models.gff
wget ftp://ftp.solgenomics.net/tomato_genome/Heinz1706/annotation/REPET/SL2.50_REPET.gff3

wget ftp://ftp.solgenomics.net/tomato_genome/assembly/current_build/S_lycopersicum_chromosomes.4.00.fa



#generate the chromosome size file - change names for respective
samtools faidx Sl4_Pl_TYLCSV.fa #converts .fa to .fa.fai
cut -f1,2 Sl4_Pl_TYLCSV.fa.fai > chromosome_size/chromosome_size_reference.txt # produces a text



################################################################################
#  1.Check quality of raw reads with FastQC
################################################################################
module load FastQC/0.11.9-Java-11

fastqc --outdir qc/1_raw_read_stats/ 1_raw/*.fastq.gz

################################################################################
# 2. Removing rRNA READS using Sortmerna
################################################################################

#source remove_rRNA_reads.sh

################################################################################
# 3. Trim raw reads using Fastp & condense report with multiqc
################################################################################

source trim.sh

################################################################################
#  Prepare reference genome
################################################################################

module load HISAT2/2.2.1-foss-2019b
module load Python/3.8.6-GCCcore-10.2.0
module load R/4.0.3-foss-2020b
module load R-bundle-Bioconductor/3.12-foss-2020b-R-4.0.3

# index genome

gunzip reference/1_raw_ref/*.fa.gz

hisat2-build reference/1_raw_ref/*.fa reference/1_raw_ref/hisat_index/index

#  generate splice junction for HISAT2 (SKIPPED)

wget https://github.com/DaehwanKimLab/hisat2/blob/d4499b468564a85430f96f260050d53541f798e5/hisat2_extract_splice_sites.py


# #generate gff file in R
# R
	#library(rtracklayer)
	# gff <- import("smb://its-rds.bham.ac.uk/rdsprojects/c/catonim-babatom/Grafting/signal-mol-project/workplace/reference/1_raw_ref/ITAG4.0_gene_models.gff")
	# export(gff, "./ITAG4.1_gene_models.gtf")
# quit(save = "no")


 #generate splice junction for HISAT2
./hisat2_extract_splice_sites.py reference/1_raw_ref/*.gtf > reference/4_splice/splice_junction.txt

awk '{if ($3=="intron") {print $1"\t"$4-2"\t"$5}}' reference/1_raw_ref/*.gtf > reference/4_splice/ssFile.table

awk '{if ($3=="intron") {print $1"\t"$4-2"\t"$5}}' reference/1_raw_ref/*.gtf > ssFile.table


################################################################################
# 4. Align/Map with Histat2
################################################################################
source map.sh

################################################################################
# 5.  Sorting, indexing
################################################################################

source sort_index.sh

################################################################################
# 6. Generating count data
################################################################################

source counts.sh

################################################################
# 7. coverage
################################################################
# setting g - this needs to be set each time. This is the size of the genome.
# for hthe Solanum Lycopersicum genome ITAG3.0 http://plants.ensembl.org/Solanum_lycopersicum/Info/Annotation/#assembly
#  == 827,747,456

g=827747456

source coverage.sh
