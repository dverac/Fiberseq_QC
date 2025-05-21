## Snakemake file to HiFi Bam QC
## Author: Diana Vera Cruz
## Produces QC tables & plots for BAM files after m6A calling. (HiFi BAM file)
## This script will run mosdepth, samtools stats, and a custom python script to get QC stats by read.
## Running commands:
##   tmux new -s this
##   conda activate fiber_sq
##   snakemake --cluster "sbatch --account=pi-spott --partition=caslake --nodes=1 --mem-per-cpu=60000" --jobs 10 --restart-times 3 --keep-going
## Requirements are low when performed at the level of chromosome.
## snakemake -s ./workflow/hifi_bam_qc.smk  --cluster "sbatch --account=pi-spott --partition=spott --nodes=1 --mem-per-cpu=50000" --jobs 50 --keep-target-files --keep-going --rerun-incomplete

## Dry run: snakemake -s ./workflow/hifi_bam_qc.smk  -np --cluster "sbatch --account=pi-spott --partition=caslake --nodes=1 --mem-per-cpu=30000"
## DAG: snakemake -s ./workflow/hifi_bam_qc.smk  --dag | dot -Tpdf > dag_base.pdf

## Configuration files 
configfile: "/project/spott/dveracruz/fiberseq_QC/workflow/config_qc.yaml"

samples = config["samples"]

## Global variables: 
# 1. bam_file: HiFi BAM file
# sam_bam = "test.bam"
# sam_name = 'test'

# 2. out_dir: output directory
#w_dir = '/project/spott/dveracruz/mini_projects/QC_plots/results'
w_dir = config["w_dir"]

## If conda environment is not activated, use the following line to activate it:
    #module load gcc/13.2.0 
    #module load cmake/3.26
    #module load python/anaconda-2024.10
    #source "$(conda info --base)/etc/profile.d/conda.sh"


## Rules

rule all: 
    input:
        expand("{dir}/{sample}/qc_report.html", sample=samples.keys(), dir=w_dir),
        expand("{dir}/{sample}/read_summary.tsv", sample=samples.keys(), dir=w_dir),


rule mosdepth_cov:
    input:
        bam=lambda wildcards: samples[wildcards.sample],
    output:
        table = "{dir}/{sample}/QC.mosdepth.summary.txt",
    params:
        dir = "{dir}/{sample}",
    shell:
        """
        source activate /project/spott/dveracruz/bin/miniconda3/envs/fiber_sq

        cd {params.dir}
        mosdepth -n -x QC {input.bam} 
        ## This might not work if the bam is not indexed!!
        """
        
rule samtools_stats:
    input:
        bam=lambda wildcards: samples[wildcards.sample],
    output:
        table = "{dir}/{sample}/samtools_stats.txt",
    params:
        dir = "{dir}/{sample}",
    shell:
        """
        source activate /project/spott/dveracruz/bin/miniconda3/envs/fiber_sq
        
        samtools stats {input.bam} > {output.table}
        """

rule modkit_metrics:
    input:
        bam=lambda wildcards: samples[wildcards.sample],
    output:
        probs_table = "{dir}/{sample}/thresholds.tsv",
        summary_table = "{dir}/{sample}/modkit_summary.txt",
    params:
        dir = "{dir}/{sample}",
    shell:
        """
        source activate /project/spott/dveracruz/bin/miniconda3/envs/dimelo_modkit_parsing
        
        modkit sample-probs --out-dir {params.dir}  --log-filepath {params.dir}/modkit.log --suppress-progress {input.bam} 
        ## Get the thresholds for the modkit summary.
        threshold=$(grep A {output.probs_table} | head -n 1 | awk '{{print $NF}}')
        modkit summary --filter-threshold $threshold --log-filepath {params.dir}/modkit.log --suppress-progress {input.bam}  > {output.summary_table}
        """

rule qc_read_pysam:
    input:
        bam=lambda wildcards: samples[wildcards.sample],
    output:
        table = "{dir}/{sample}/QC_stats_byRead.tsv",
    shell:
        """
        source activate /project/spott/dveracruz/bin/miniconda3/envs/fiber_sq

        ## Python script to get QC stats by read
        python /project/spott/dveracruz/fiberseq_QC/workflow/scripts/qc_by_read.py --bam_file {input.bam} --output_file {output.table} --adenine_threshold 160 --cytosine_threshold 180
        """

    
rule r_wrapper:
    input:
        mosdepth = "{dir}/{sample}/QC.mosdepth.summary.txt",
        modkit_probs = "{dir}/{sample}/thresholds.tsv",
        modkit_summary = "{dir}/{sample}/modkit_summary.txt",
        byread = "{dir}/{sample}/QC_stats_byRead.tsv",
        samtools = "{dir}/{sample}/samtools_stats.txt",
        bam =lambda wildcards: samples[wildcards.sample],
    output:
        html = "{dir}/{sample}/qc_report.html",
        table = "{dir}/{sample}/read_summary.tsv",
    params:
        dataset = "{dir}/{sample}",
        sam_name = "{sample}",
        platform = lambda wildcards: config['platform'],
    shell:
        """
        echo "Running R script to summarize and plot QC stats."
        echo "HTML report will be saved at {output.html}"
        echo "Table will be saved at {output.table}"

        source activate /project/spott/dveracruz/bin/miniconda3/envs/fiber_sq

        ## R markdown script to summarize and plot QC stats. 
        Rscript -e "rmarkdown::render('/project/spott/dveracruz/fiberseq_QC/workflow/scripts/QC_report.Rmd', output_file = '{params.dataset}/qc_report.html', params = list(sample = '{params.sam_name}', w_dir = '{params.dataset}', bam = '{input.bam}', platform = '{params.platform}'))"

        ## Print time/date at the end of the script
        echo "Done"
        """
