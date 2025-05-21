## Python script to parse HiFi BAM and generate per read statistics. 
## Author: Diana Vera Cruz
## Last updated: 2025-01-24

## Run: python qc_by_read.py --bam_file path/to/input.bam --output_file path/to/output.tsv --adenine_threshold 160 --cytosine_threshold 180

## Libraries
import pysam
import random
import os
import argparse
import statistics

## Function
## Function to obtain stats per read: Length, number of methylated A and C. 
def stats_by_read(
    bam_file: str,
    output_file: str,
    adenine_threshold: int = 160,
    cytosine_threshold: int = 180
) -> None:
    """
    Summary stats per read in a BAM file: Read metadata + methylation, including percentages relative to total adenines/cytosines.

    Args:
        bam_file (str): Path to the input BAM file.
        output_file (str): Path to the output TSV file.
        adenine_threshold (int): Threshold for adenine methylation probability (default: 160).
        cytosine_threshold (int): Threshold for cytosine methylation probability (default: 180).
    """
    # Open BAM file and output file
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_file, "w") as out:
        # Write header
        out.write(
            "ReadID\tLength\tcontig\tStart\tEnd\tStrand\tHaplotype\tMapQ\tMedianQV\tclass\t"
            "MedianNucleosomeLength\tNumNucleosomes\tNumMethylatedA\tNumMethylatedC\tFractionMethylatedA\t"
            "FractionMethylatedC\tPercentageMethylatedA\tPercentageMethylatedC\n"
        )

        ## Start counter to check total reads processed.
        total_reads=0
        unmapped_reads=0
        noseq_reads=0
        notag_reads=0

        ## MAIN for loop to process reads.
        for read in bam.fetch(until_eof=True):
            ## Increase total reads counter. 
            total_reads+=1

            read_id = read.query_name
            read_length = read.query_length
            if read.query_qualities:
                median_QV = statistics.median(read.query_qualities)
            else:
                median_QV = 0

            # Include unmapped reads, but just to count them.
            if read.is_unmapped:
                unmapped_reads+=1
                out.write(
                    f"{read_id}\t{read_length}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tUnmapped\t"
                    f"0.0000\t0\t0\t0\t0.0000\t0.0000\t0.0000\t0.0000\n"
                )
                continue

            # Read metadata
            contig = read.reference_name
            strand = "-" if read.is_reverse else "+"
            start = read.reference_start
            end = read.reference_end
            mapq = read.mapping_quality
            # median_QV = statistics.median(read.query_qualities)

            # Get the haplotype tag if available
            haplotype = read.get_tag("HP") if read.has_tag("HP") else "NA"

            # Check if query_sequence is available
            if read.query_sequence is None:
                # If query_sequence is missing, skip the read
                noseq_reads+=1
                out.write(
                    f"{read_id}\t{read_length}\t{contig}\t{start}\t{end}\t{strand}\t{haplotype}\t{mapq}\t{median_QV}\tNoSequence\t"
                    f"0.0000\t0\t0\t0\t0.0000\t0.0000\t0.0000\t0.0000\n"
                )
                continue

            # Count total adenines, tiamines and CG in the read sequence
            query_sequence = read.query_sequence
            #total_adenines = query_sequence.count("A") + query_sequence.count("T")
            total_adenines = query_sequence.count("A")
            total_tiamines = query_sequence.count("T")

            ## CG sites. 
            total_cg = query_sequence.count("CG")

            # Initialize methylation counters
            num_methylated_a = 0
            num_methylated_ta = 0 ## This is for PacBio 6mA methylation in T context.
            num_methylated_c = 0

            total_a_positions = 0
            total_ta_positions = 0 ## This is for PacBio 6mA methylation in T context.
            total_c_positions = 0
            flag = 0 ## Flag to determine whether the T-a methylation is present (PacBio) or not (Nanopore).

            ## Extract Nucleosomes info: Tags from ft: ns-positions. 
            if read.has_tag("ns"):
                nuc_positions = read.get_tag("ns")
                ## Check the type of nuc_positions is not a single value or empty.
                ## Let's do the inverse, check if nuc_positions is an array, a
                if nuc_positions and len(nuc_positions) >= 2:
                    ## Check if the length is an even number: If so, remove last value.
                    if len(nuc_positions) % 2 != 0:
                        nuc_positions =nuc_positions[:-1]
                    ## Get substraction of even positions / odd positions.
                    nuc_ln = [x - y for x, y in zip(nuc_positions[1::2], nuc_positions[::2])]
                    ## Get the median of the nucleosome length.
                    median_nuc = statistics.median(nuc_ln)
                    ## Get the total number of nucleosomes.
                    num_nuc = len(nuc_ln)
                else:
                    num_nuc = "NA"
                    median_nuc = "NA"
            else:
                ## Set default values to NA if no nucleosome tag is found.
                num_nuc = "NA"
                median_nuc = "NA"

            # Extract MM and ML tags if available
            if read.has_tag("MM") and read.has_tag("ML"):
                mm_tag = read.get_tag("MM")
                ml_tag = list(read.get_tag("ML"))  # ML is a list of probabilities

                # Parse MM tag for methylation data
                methylation_data = mm_tag.split(";")

                # Process adenine methylation (A+a)
                adenine_methylation = [
                    entry for entry in methylation_data if entry.startswith("A+a")
                ]
                if adenine_methylation:
                    a_context_data = adenine_methylation[0]  # e.g., "A+a,1,10,20"
                    a_positions = a_context_data.split(",")[1:]
                    a_positions = [int(pos) for pos in a_positions if pos.isdigit()]

                    # Count methylated adenines above the threshold
                    for pos in a_positions:
                        if pos < len(ml_tag) and ml_tag[pos] > adenine_threshold:
                            num_methylated_a += 1
                        total_a_positions += 1
                

                # Process adenine methylation (T-a) if present -> Adenine methylation in T context for PacBio.  
                adenine_pbTa_methylation = [
                    entry for entry in methylation_data if entry.startswith("T-a")
                ]
                if adenine_pbTa_methylation:
                    flag = 1 ## Assign flag to 1 if T-a methylation is present.
                    t_context_data = adenine_pbTa_methylation[0]  # e.g., "A+a,1,10,20"
                    t_positions = t_context_data.split(",")[1:]
                    t_positions = [int(pos) for pos in t_positions if pos.isdigit()]

                    # Count methylated adenines above the threshold
                    for pos in t_positions:
                        if pos < len(ml_tag) and ml_tag[pos] > adenine_threshold:
                            num_methylated_ta += 1
                        total_ta_positions += 1
                

                # Process cytosine methylation (C+m)
                cytosine_methylation = [
                    entry for entry in methylation_data if entry.startswith("C+m")
                ]
                if cytosine_methylation:
                    c_context_data = cytosine_methylation[0]  # e.g., "C+m,5,15,25"
                    c_positions = c_context_data.split(",")[1:]
                    c_positions = [int(pos) for pos in c_positions if pos.isdigit()]

                    # Count methylated cytosines above the threshold
                    for pos in c_positions:
                        if pos < len(ml_tag) and ml_tag[pos] > cytosine_threshold:
                            num_methylated_c += 1
                        total_c_positions += 1
            else:
                # If MM and ML tags are missing, skip the read
                out.write(
                    f"{read_id}\t{read_length}\t{contig}\t{start}\t{end}\t{strand}\t{haplotype}\t{mapq}\t{median_QV}\tNoMM-MLTag\t"
                    f"{median_nuc}\t{num_nuc}\t0\t0\t0.0000\t0.0000\t0.0000\t0.0000\n"
                )
                continue


            ## Update num_methylated_a if flag is 1 (PacBio T-a methylation)
            ## If pacbio: Add Ta methylation to A methylation. 
            if flag == 1:
                total_adenines += total_tiamines
                num_methylated_a += num_methylated_ta
                total_a_positions += total_ta_positions

            # Calculate fractions of methylated sites (based on positions reported in MM tag)
            fraction_methylated_a = (
                num_methylated_a / total_a_positions if total_a_positions > 0 else 0
            )
            fraction_methylated_c = (
                num_methylated_c / total_c_positions if total_c_positions > 0 else 0
            )

            # Calculate percentages of methylation (relative to total bases in the read)
            percentage_methylated_a = (
                (num_methylated_a / total_adenines) * 100 if total_adenines > 0 else 0
            )
            percentage_methylated_c = (
                (num_methylated_c / total_cg) * 100 if total_cg > 0 else 0
            )

            # Write the results for the read
            out.write(
                f"{read_id}\t{read_length}\t{contig}\t{start}\t{end}\t{strand}\t{haplotype}\t{mapq}\t{median_QV}\tComplete\t" 
                f"{median_nuc}\t{num_nuc}\t{num_methylated_a}\t{num_methylated_c}\t"
                f"{fraction_methylated_a:.4f}\t{fraction_methylated_c:.4f}\t"
                f"{percentage_methylated_a:.2f}\t{percentage_methylated_c:.2f}\n"
            )
        ## Print counters. 
        print(total_reads)
        print(unmapped_reads)
        print(noseq_reads)
        print(notag_reads)


## MAIN 

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Run Statistics by read in a HiFi BAM file")

    # Define expected parameters
    parser.add_argument("--bam_file", type=str, required=True, help="Path to the input BAM file.")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output TSV file.")
    parser.add_argument("--adenine_threshold", type=int, default=160, help="Adenine methylation score threshold. Default 160")
    parser.add_argument("--cytosine_threshold", type=int, default=180, help="Cytosine methylation score threshold. Default 180")

    # Parse arguments
    args = parser.parse_args()
    ## Check first the bam file exists.
    if not os.path.exists(args.bam_file):
        raise FileNotFoundError(f"File not found: {args.bam_file}")
    ## Print first current time/date
    #print("Start time: ", datetime.now())
    ## Run main function
    stats_by_read(
        bam_file=args.bam_file,
        output_file=args.output_file,
        adenine_threshold=args.adenine_threshold,
        cytosine_threshold=args.cytosine_threshold
    )

