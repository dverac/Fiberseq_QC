##  FIRE output parsing: Extracting nucleosomes + storing the lengths per each. 
# Base file: fire_dir/trackHub/bb/fire-fiber-decorators.bb
## Since the fiber-calls only include FIREs, the nucleosomes calls can only be extracted from the trackHub. 
import pandas as pd
import gzip


def extract_nuc(fire_file):
    """
    Extract nucleosome calls from a FIRE file.
    Args:
        fire_file (str): Path to the FIRE file: fire_dir/trackHub/bb/fire-fiber-decorators.bb
    Returns:
        Dictionary with nucleosome lenghts. 
    """


# Define input file path
input_file = "your_file.bed.gz"  # Change this to your actual gzipped BED file

# Initialize lists to store values for aggregation
col10_values = []
region_lengths = []

# Read the gzip-compressed BED file in chunks
chunksize = 100000  # Adjust chunk size as needed for memory efficiency

with gzip.open(input_file, "rt") as f:  # Open as a text file
    for chunk in pd.read_csv(f, sep="\t", header=None, chunksize=chunksize, dtype={8: str}):
        # Filter rows where column 9 (index 8) is "147,112,219"
        filtered_chunk = chunk[chunk[8] == "147,112,219"]

        # Append column 10 (index 9) values
        col10_values.extend(filtered_chunk[9].astype(float).tolist())

        # Compute region length (column 3 - column 2) and store
        region_lengths.extend((filtered_chunk[2] - filtered_chunk[1]).tolist())

# Compute mean and median
mean_col10 = sum(col10_values) / len(col10_values) if col10_values else 0
median_col10 = pd.Series(col10_values).median() if col10_values else 0

mean_region_length = sum(region_lengths) / len(region_lengths) if region_lengths else 0
median_region_length = pd.Series(region_lengths).median() if region_lengths else 0

# Print results
print(f"Mean of column 10: {mean_col10}")
print(f"Median of column 10: {median_col10}")
print(f"Mean of region length: {mean_region_length}")
print(f"Median of region length: {median_region_length}")




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
    process_bed_bam(args.bed, args.bam, args.output)