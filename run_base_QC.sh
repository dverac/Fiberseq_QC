#!/bin/bash

#SBATCH --job-name=QC-base
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=4:00:00
#SBATCH --account=pi-spott
#SBATCH --partition=spott
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

###############
##  Parameters 
## Use the 
## Rscript /project/spott/dveracruz/fiberseq_QC/workflow/generate_config_run.r to generate the config file.
###############

## Uncomment the following lines to load conda if you do not have it in your path.
#module load gcc/13.2.0 
#module load cmake/3.26
#module load python/anaconda-2024.10

## Set configuration file.
config_file=workflow/config_qc.yaml 

## Start output: Show start time & configuration. 
echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

## Check the config file exists. 
if [ ! -f $config_file ]; then
    echo "Configuration file not found!, expected $config_file"
    exit 1
fi

## Pipeline path.
smk=/project/spott/dveracruz/fiberseq_QC/workflow
## Pipeline to run: 
pp=hifi_bam_qc.smk

## Run snakemake pipeline. 
source activate /project/spott/dveracruz/bin/miniconda3/envs/fiber_sq  

snakemake -s $smk/$pp --unlock

snakemake -s $smk/$pp --configfile $config_file \
    --cluster "sbatch --account=pi-spott --partition=caslake --nodes=1 --mem-per-cpu=60000" \
    --jobs 50 --keep-target-files --keep-going --rerun-incomplete


echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins

