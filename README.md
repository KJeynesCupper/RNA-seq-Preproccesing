# RNA-seq-Preproccesing
Bash scripts to automatically pre-process RNA-seq data in fastq format into raw count data. Keep in mind that some parameters may need amending, and you need to manually access the quality from the QC output.

This has been set up to run primarily on the BEAR server - meeting its requirements. 
## Files 

* RNAseq_preprocessing.sh
* remove_rRNA_reads.sh
* trim.sh
* map.sh
* sort_index.sh
* count.sh
* coverage.sh

## Running 

The "RNAseq_preprocessing.sh" file is the primary router and needs to be placed in a new "workplace" directory where all proccess will be executed. 
This file runs through all the processes and will stop is an errors occurs. 

For quality control, it is important that you check the outputs and amend the parameters where required. 

Each function within the individual shells have parameters that may need adjusting, please refer the the specific packages' vignette. 

## Output

This will produce the raw count files (non-normalised) for each replicate. As well as summary statistics for different steps. 
