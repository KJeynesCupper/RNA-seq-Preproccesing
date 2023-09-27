# RNA-seq-Preproccesing
This contains RNAseq preprocessing scripts. Keep in mind that directories and parameters may need amending. The scripts contain job pars for BEAR. 

## Files 
* `1_mRNAseq_preprocessing.sh`
* `2_coverage.sh`
* `3_extract_metrics.sh`
* sequential:
* * `RNAseq_preprocessing.sh`
* * `remove_rRNA_reads.sh`
* * `trim.sh`
* * `map.sh`
* * `sort_index.sh`
* * `count.sh`
* * `coverage.sh`

## Running Option 1
Download, amend and run the `mRNAseq_preprocessing.sh` script which in total will do QC with fastQC, trim with trimmomatic, map with HISAT2 and count raw read with HTSeq.  
Then Download, amend and run the `coverage.sh` script to generate a Bedgraph summarising the mapped reads at each bp, estimate the average reads across the genome and use the Perl command to normalise read counts by genome-wide average counts. 


## Running Option 2: sequential 
The sequential folder stores scripts to run the processing slightly differently. The `RNAseq_preprocessing.sh` file is the primary router and needs to be placed in a new "workplace" directory where all proccess will be executed. This file runs through all the processes and will stop is an errors occurs. 

## Output
This will produce the raw count files (non-normalised) for each replicate. As well as summary statistics for different steps. 

## Pull metrics 
Download, amend and run the `3_extract_metrics.sh` script to generate txt files containing information on number of reads in raw/trimmed FASTQ files, and get mapping stats for each sample. 
