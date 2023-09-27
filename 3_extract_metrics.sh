#!/bin/bash
#SBATCH --qos castles
#SBATCH --ntasks 72
#SBATCH --nodes 1
#SBATCH --time 6000
module purge; module load bluebear
module load SAMtools/1.10-GCC-8.3.0

# Generate txt file containing number of  reads in each file:
for file in *.fq; do
    read_count=$(wc -l < "$file")
    # Calculate the number of reads by dividing the total lines by 4 (each read has 4 lines)
    read_count=$((read_count / 4))
    echo "$file $read_count" >> read_counts.txt
done



# Output file to store the results
output_file="3_mapped_reads_info.txt"

# Clear the output file if it already exists
> "$output_file"

# Loop through sample folders
for sample in ./3_mapped/*_sorted.bam; do
        # BAM file path
         bam_file="$sample""_sorted.bam"

        # Count unmapped reads
        unmapped=$(samtools view -c -f 4 "$bam_file")

        # count multimappers
        multi_mapped=$(samtools view -c -q 1 "$bam_file")

        # unique
         uniquely_mapped=$(samtools view -c -F 0x4 "$bam_file")

        # Total reads in the BAM file
        total_reads=$(samtools view -c "$bam_file")



        # Print the results to the output file
        echo "Sample: $sample_name" >> "$output_file"
        echo "Uniquely Mapped Reads: $uniquely_mapped" >> "$output_file"
        echo "Multi Mapped Reads: $multi_mapped" >> "$output_file"
        echo "Unmapped Reads: $unmapped" >> "$output_file"
        echo "Total Reads: $total_reads" >> "$output_file"

        echo "" >> "$output_file"  # Add a blank line for separation
done

echo "All sample information has been extracted and saved to $output_file"
