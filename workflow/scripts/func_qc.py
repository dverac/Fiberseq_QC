## Python script to generate pileup: Coverage + methylation fraction per promoter. 
## Author: Diana Vera Cruz
## Last updated: 2025-01-30
## Example run: python func_qc.py --bed path/to/promoter.bed --bam path/to/input.bam --output path/to/output.tsv
## Try to see if the new version for methylation proportion works. 

import os
import pysam
import pandas as pd
import numpy as np
from collections import Counter
import statistics


def methylation_per_region(bed_file, bam_file, output_file):
    """
    Process a BED-like file and BAM file to compute coverage and adenine methylation per position.
    
    Parameters:
    - bed_file: Path to the BED-like input file (tab-separated).
    - bam_file: Path to the BAM file (indexed).
    - output_file: Path to save the processed results.
    """
    
    # Open the BAM file using pysam
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Read BED-like file
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, comment="#")
    bed_df = bed_df.iloc[:, :7]
    
    # Ensure correct columns are used
    bed_df.columns = ['chrom', 'stf', 'end', 'name', 'score', 'strand', 'tss']  

    ## Start output file with column names: 
    with open(output_file, "w") as f:
        f.write("chrom\tstart\tend\tstrand\tTSS\tMean_methylation_prop\tCoverage\n")

    results = []
    bed_n = 0

    for _, row in bed_df.iterrows():
        chrom = row["chrom"]
        coord = int(row["tss"])
        strand = row["strand"]

        print(f"Processing {chrom}:{coord} ({strand})")

        # Define the region of interest based on strand
        if strand == "+":
            start = coord - 200
            end = coord + 50
        else:
            start = coord - 50
            end = coord + 200
        
        # Ensure valid genomic coordinates
        if start < 0:
            start = 0

        bed_n += 1
        #if bed_n > 50:
        #    break
    
        ## Start counter for coverage in the region. 
        cov_counter = 0
        methylation_prop = []
        methylation_matrix = {}

        ##############################
        ##  Fetch reads in the region of interest. 
        ##############################
        for read in bam.fetch(chrom, start, end):
            try:
                mm_tag = read.get_tag("MM")  # Base modifications
                ml_tag = read.get_tag("ML")  # Modification likelihoods
            except KeyError:
                continue  # Skip reads without MM/ML tags

            ##############################
            ## 1. Calculate AT content in subregion of read within range of interest.
            ##############################
            read_seq = read.query_sequence # Get read sequence. 
            aligned_pairs = read.get_aligned_pairs(matches_only=True) ## Get aligned pairs. (query_pos, reference_pos)
            # aligned_pairs = read.get_aligned_pairs(with_seq=True) ## Get aligned pairs. (query_pos, reference_pos) Reads do not have MD tag

            ## Get the read sequence in region of interest.  start = aligned_pairs[0], end = aligned_pairs[-1]
            read_start_pos = None
            read_end_pos = None
            ## Iterate over aligned pairs to find the start and end positions of the read within the region
            for query_pos, ref_pos in aligned_pairs:
                if ref_pos is None:
                    continue   ## Pass unmapped bases. 
                if ref_pos >= start and read_start_pos is None:
                    #print("start", query_pos, ref_pos)
                    read_start_pos = query_pos  # First match inside the region
                if ref_pos >= end - 1:
                    #print("end", query_pos, ref_pos)
                    read_end_pos = query_pos + 1  # Last match before exceeding the region
                    break
            # If end was not set, take full read length
            if read_start_pos is not None and read_end_pos is None:
                read_end_pos = len(read_seq)

            ## If region within read cannot be detected in at least 10 bases, skip to next read.
            if read_start_pos is None or read_end_pos is None or (read_end_pos - read_start_pos) <= 10:
                continue

            ## Calculate AT content in subregion of read within range of interest, and increment coverage counter.
            sub_read_seq = read_seq[read_start_pos:read_end_pos]
            at_sites = sub_read_seq.count('A') #+ sub_read_seq.count('T')
            t_sites = sub_read_seq.count('T') ## Only A for now, adding T if PacBio T-a label is used. 
            cov_counter += 1 # Increment coverage counter

            ##############################
            ## 2.  Get Methylation positions in read. 
            ##############################

            ml_scores = ml_tag.tolist()

            mod_df = pd.DataFrame({'pos': [], 'mod': []})
            mod_ix = 1
            for entry in mm_tag.split(";"): 
                if "," not in entry:
                    continue
                mod_base, positions = entry.split(",", 1)  # Example: "A+m?,5,10,15" 
                abs_positions = [] ## Store absolute positions in MM tag. 
                pos = 0  # Relative to read start
                for rel_pos in map(int, positions.split(",")):
                    pos += rel_pos
                    abs_positions.append(pos)
                    mod_ix += 1
                ## Add rows to mod_df: Pos, moditification type. 
                mod_df = pd.concat([mod_df, pd.DataFrame({'pos': abs_positions, 'mod': [mod_base] * len(abs_positions)})], ignore_index=True)
                
            ## Add score column to mod_df
            mod_df['score'] = list(ml_tag)

            ## Modify mod column: Keep only the first 3 letters per element in the column.
            #mod_df['mod'] = mod_df['mod'].str[:3]
            mod_df['mod'] = mod_df['mod'].astype(str).str[:3]

            ## Count frequency of each modification type: column 'mod'
            mod_counts = mod_df['mod'].value_counts()
            #print(f"mod_counts: {mod_counts}")
            

            ## Check if mod T-a is present in the read.
            ## If mod T-a is present, filter to keep only positions with mod == 'A+' and score >= 160
            flag = (mod_df['mod'] == 'T-a').sum() > 0 # Flag -> If we have at least 1 entry from T-a.
            if flag > 0:
                at_sites += t_sites
            #print(f"at_sites: {at_sites}, t_sites: {t_sites}, flag: {flag}")

            ## Filter mod_df to keep only positions within the region of interest, mod == 'A+' and score >= 160
            # mod_df = mod_df[(mod_df['mod'] == 'A+a') & (mod_df['score'] >= 160)]
            mod_df = mod_df[(mod_df['mod'].isin(['A+a', 'T-a'])) & (mod_df['score'] >= 160)]

            ## Print the number of entries in mod_df.
            #if mod_df.empty:
            #    print(f"No valid methylation positions found for {chrom}:{start}-{end}")

            ## Number of entries in positions of interest. 
            met_A = sum ( (mod_df['pos'] >= read_start_pos) & (mod_df['pos'] < read_end_pos) )
            ## Append to methylation_prop list.
            ## Only if at_sites > 0. 
            if at_sites > 0:
                methylation_prop.append(met_A / at_sites)

        ##############################
        ## END OF FOR LOOP PER READ: SUMMARY PER REGION IN BED.
        ##############################
        # Convert methylation table to data.frame
        #methylation_df = pd.DataFrame.from_dict(
        #    {pos: np.mean(scores) for pos, scores in methylation_matrix.items()},
        #    orient="index",
        #    columns=["Avg_6mA_Methylation_Score"]
        #).sort_index()
        
        ## Filter to keep only positions with avg. score above threshold (> 160)
        # methylation_df = methylation_df[methylation_df["Avg_6mA_Methylation_Score"] > 160]
        ## Average methylation: Number of rows / length (end - start + 1) -> Updated to length of methylation_df / at_sites
        # avg_methylation = len(methylation_df) / (end - start + 1)
        if len(methylation_prop) > 0:
            avg_methylation = np.mean(methylation_prop)
        else:
            avg_methylation = 0  # Default to 0 if no methylation data is available

        ## Print to output file: Chrom, start, end, strand, avg_methylation & coverage_counter. 
        with open(output_file, "a") as f:
            f.write(f"{chrom}\t{start}\t{end}\t{strand}\t{coord}\t{avg_methylation}\t{cov_counter}\n")
            
    print(f"Processing completed. Output saved to: {output_file}")
    bam.close()


## MAIN
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process a BED-like file and BAM file for coverage and adenine methylation analysis.")
    parser.add_argument("-b", "--bed", required=True, help="Path to the BED-like input file.")
    parser.add_argument("-bam", "--bam", required=True, help="Path to the BAM file (indexed).")
    parser.add_argument("-o", "--output", required=True, help="Path to save the output TSV file.")

    args = parser.parse_args()
    ## Check that the input files exist. 
    if not os.path.exists(args.bed):
        raise FileNotFoundError(f"BED-like file not found: {args.bed}")
    if not os.path.exists(args.bam):
        raise FileNotFoundError(f"BAM file not found: {args.bam}")

    ## Function
    methylation_per_region(args.bed, args.bam, args.output)