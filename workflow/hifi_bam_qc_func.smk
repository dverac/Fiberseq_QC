## Snakemake file to HiFi Bam QC
## Author: Diana Vera Cruz
## Produces QC tables & plots for BAM files after m6A calling. (HiFi BAM file)
## This script will run mosdepth, samtools stats, and a custom python script to get QC stats by read.
## Running commands:
##   tmux new -s this
##   conda activate fiber_sq
## snakemake -s ./workflow/hifi_bam_qc_func.smk  --cluster "sbatch --account=pi-spott --partition=spott --nodes=1 --mem-per-cpu=50000" --jobs 50 --keep-target-files --keep-going --rerun-incomplete
## Dry run: snakemake -s ./workflow/hifi_bam_qc_func.smk  -np --cluster "sbatch --account=pi-spott --partition=caslake --nodes=1 --mem-per-cpu=20000"
## DAG:  snakemake -s ./workflow/hifi_bam_qc_func.smk  --dag | dot -Tpdf > dag.pdf

## Configuration files 
configfile: "/project/spott/dveracruz/fiberseq_QC/workflow/config_qc_func.yaml"

## Load variables: Samples, fire_dirs and TSS, platform. 
samples = config["samples"]
fire_dirs = config["fire_dirs"]
tss = config['TSS']
platform = config['platform'] ## Used only to include in the HTML Report. 

w_dir = config['w_dir']

# 2. out_dir: output directory
#w_dir = '/project/spott/dveracruz/mini_projects/QC_plots/results'

## Rules

rule all: 
    input:
        expand("{dir}/{sample}/qc_report_func.html", sample=samples.keys(), dir=w_dir),
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
        module load gcc/13.2.0 
        module load cmake/3.26
        module load python/anaconda-2024.10

        source activate /project/spott/dveracruz/bin/miniconda3/envs/fiber_sq

        cd {params.dir}
        mosdepth -n -x QC {input.bam} 
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

## Stats per read: Length, methylation, haplotypes and also nucleosomes. 
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

rule promoter_stats:
    input:
        bam=lambda wildcards: samples[wildcards.sample],
        tss = {tss},
    output:
        table = "{dir}/{sample}/promoter_stats.tsv",
    shell:
        """
        source activate /project/spott/dveracruz/bin/miniconda3/envs/fiber_sq

        python /project/spott/dveracruz/fiberseq_QC/workflow/scripts/func_qc.py --bed {input.tss} --bam {input.bam} --output {output.table}
        """

rule fire_bed: 
    input:
        fire = lambda wildcards: fire_dirs[wildcards.sample],
    output: 
        "{dir}/{sample}/fire.bed",
    shell:
        """
        ## Convert bigbed to bed, temporary file.
        /project/spott/dveracruz/bin/UCSC_kent/bigBedToBed {input.fire}/trackHub-v0.1.1/bb/fire-fiber-decorators.bb {output}
        """

rule fire_stats: 
    input:
        bed = "{dir}/{sample}/fire.bed",
        #fire = lambda wildcards: fire_dirs[wildcards.sample],
    output: 
        fire = "{dir}/{sample}/fire_stats.tsv",
    shell:
        """
        ## Filter FIRE regions: Get class, length and methylation. 
        # zcat {input.bed}/fiber-calls/FIRE.bed.gz |  awk -F'\t' '{{print $9 "\t" $3 - $2 "\t" $10}}' > {output.fire}

        ## FIRE regions: Get class and length. 
        grep 'FIRE' {input.bed} | awk -F'\t' '{{print $15 "\t" $3 - $2 }}' > {output.fire}
        echo 'FIRE done'
        """

rule r_wrapper:
    input:
        mosdepth = "{dir}/{sample}/QC.mosdepth.summary.txt",
        byread = "{dir}/{sample}/QC_stats_byRead.tsv",
        samtools = "{dir}/{sample}/samtools_stats.txt",
        promoter = "{dir}/{sample}/promoter_stats.tsv",
        fire = "{dir}/{sample}/fire_stats.tsv",
        modkit_probs = "{dir}/{sample}/thresholds.tsv",
        modkit_summary = "{dir}/{sample}/modkit_summary.txt",
        bam =lambda wildcards: samples[wildcards.sample],
    output:
        html = "{dir}/{sample}/qc_report_func.html",
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
        Rscript -e "rmarkdown::render('/project/spott/dveracruz/fiberseq_QC/workflow/scripts/QC_report_func.Rmd', output_file = '{params.dataset}/qc_report_func.html', params = list(sample = '{params.sam_name}', w_dir = '{params.dataset}', bam = '{input.bam}', platform = '{params.platform}'))"

        ## Print time/date at the end of the script
        echo "Done"
        """
