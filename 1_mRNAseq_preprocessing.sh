#!/bin/bash
#SBATCH --qos castles
#SBATCH --ntasks 72
#SBATCH --nodes 2
#SBATCH --time 6000

module purge; module load bluebear
module load HISAT2/2.2.1-foss-2019b
module load SAMtools/1.10-GCC-8.3.0
module load HTSeq/0.12.4-foss-2019b-Python-3.7.4
module load Trimmomatic/0.39-Java-11
module load Java/11
module load BioPerl/1.7.2-GCCcore-8.3.0-Perl-5.30.0
module load HTSeq/0.12.4-foss-2019b-Python-3.7.4


# index genome
hisat2-build ./1_reference/genome_reference.fa ./1_reference/genome_reference


# Define paths to genome index, annotation file, and raw data directory:
genome_index="./1_reference/genome_reference"
gff="./1_annotation/genome_annotation.gff"
data_dir="./1_raw"

# Create the output directories if they don't exist
mkdir -p 2_trimmed
mkdir -p 3_mapped
mkdir -p 4_counts

# generate fastqc reports for raw files
mkdir -p $data_dir/qc/
for i in $data_dir/*.fq.gz; do fastqc $i -o $data_dir/qc/ ; done

# Loop over each file in the data directory
for replicate1 in $data_dir/*_L1_1.fq.gz; do
    # Extract the sample name from the filename
    sample_name=$(basename $replicate1 _L1_1.fq.gz)

    # Define the corresponding R2 read filename based on the sample name
    replicate2="$data_dir/${sample_name}_L1_2.fq.gz"
    replicate3="$data_dir/${sample_name}_L2_1.fq.gz"
    replicate4="$data_dir/${sample_name}_L2_2.fq.gz"

    # Step 2: Alignment with HISAT2
    trimmed_dir="./2_trimmed"
    mapped_dir="./3_mapped"

    trimmed_replicate1="$trimmed_dir/$(basename $replicate1 .fq.gz)_trimmed.fq.gz"
    trimmed_replicate2="$trimmed_dir/$(basename $replicate2 .fq.gz)_trimmed.fq.gz"
    trimmed_replicate3="$trimmed_dir/$(basename $replicate3 .fq.gz)_trimmed.fq.gz"
    trimmed_replicate4="$trimmed_dir/$(basename $replicate4 .fq.gz)_trimmed.fq.gz"

    if [ -e "$replicate3" ] && [ -e "$replicate4" ]; then
        # Both R1 and R2 are present, trim:
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 $replicate1 $replicate2 $trimmed_replicate1 /dev/null $trimmed_replicate2 /dev/null ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 $replicate3 $replicate4 $trimmed_replicate3 /dev/null $trimmed_replicate4 /dev/null ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        #perform paired-end mapping:
        echo "Performing paired-end mapping for sample $sample_name"
        hisat2 -p 6 -x $genome_index -1 $trimmed_replicate1,$trimmed_replicate3 -2 $trimmed_replicate2,$trimmed_replicate4 -S $mapped_dir/"${sample_name}_aligned.sam"
    else
        # Only R1 is present, perform single-end mapping
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 $replicate1 $replicate2 $trimmed_replicate1 /dev/null $trimmed_replicate2 /dev/null ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        echo "Performing single-end mapping for sample $sample_name"
        hisat2 -p 6 -x $genome_index -U $trimmed_replicate1,$trimmed_replicate2 -S $mapped_dir/"${sample_name}_aligned.sam"
    fi
    # Continue with the common steps for both cases
    output_sam="$mapped_dir/${sample_name}_aligned.sam"
    output_bam="$mapped_dir/${sample_name}_aligned.bam"
    sorted_bam="$mapped_dir/${sample_name}_sorted.bam"

    samtools view -b -o $output_bam $output_sam
    samtools sort -o $sorted_bam $output_bam
    samtools index $sorted_bam

    # Step 3: Read Counting with HTSeq
    counts_dir="./4_counts"
    count_file="$counts_dir/${sample_name}_gene_counts.txt"
    python -m HTSeq.scripts.count --format=bam --order=pos --stranded=no --mode=union --nonunique=none --type=gene --idattr=Name $sorted_bam $gff > $count_file

    # Clean up intermediate files if needed
    rm $output_sam $output_bam

    echo "Sample $sample_name completed. Gene counts are saved in $count_file."
    echo "....  .... .... .... .... .... .... .... .... .... .... .... .... ...."
done

# fastqc for trimmed files
mkdir -p $trimmed_dir/qc/
for i in $trimmed_dir/*.fq.gz; do fastqc $i -o $trimmed_dir/qc/ ; done

echo "All samples processed."
