## Quick script to generate the input table for FIRE or other snakemake pipelines.
## Table includes 2 columns: Sample_name + Bam file full path. 
## Rscript /project/spott/dveracruz/fiberseq_QC/workflow/generate_samples_table.r /path/to/bam_files/ type suffix > output_table.txt   
## type argument can be either 'tab' or 'smk'

args <- commandArgs(trailingOnly=TRUE)

bam_dir <- args[1]
output_type <- args[2]
## Check if suffix is provided, if not, default to bam. 
if(length(args) > 2) suffix <- args[3] else suffix = 'bam'

## Check if output_type is provided, if not, default to tab.
if( is.na(output_type) ) output_type = 'tab'
if(!output_type %in% c('tab', 'smk')){
  stop("Invalid type argument. Must be either 'tab' for tabular table or 'smk' for snakemake config yaml")
}

## Check bam_dir exists. 
if (!file.exists(bam_dir)) {
  stop("Input directory does not exist.")
}

## Get all the bam files in the directory
bam_files <- list.files(bam_dir, pattern = paste0(suffix,"$"), full.names = TRUE)

## Check if bam_files is empty.
if (length(bam_files) == 0) {
  stop("No bam files found in the input directory.")
}

## Get sample names. 
sample_names <- gsub(".bam|.aligned|.sort|.5mC|.6mA|.nuc|.phased", "", basename(bam_files))

## Write table or print to console instead. 
if(output_type == 'tab'){
    write.table(data.frame(sample = sample_names, bam = bam_files), file = "", sep = "\t", row.names = FALSE, quote = FALSE)
}
if(output_type == 'smk'){
    cat("samples:\n")
    for(i in 1:length(sample_names)){
        cat(paste0("\t", sample_names[i], ": ", bam_files[i], "\n"))
    }
}

