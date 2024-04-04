import sys
import os
import glob
import argparse
import pandas as pd


def process_read_counts(READCOUNTS_file_path):
    # Function to process a READ_COUNTS.txt file
        # returns the percent of reads matching the ref genomes (perc_match) and the total reads passing filtering (total_reads)
    df = pd.read_csv(READCOUNTS_file_path, sep='\t').set_index('Record').fillna(0)
    try:
        match = df.loc['3-match', 'Reads']
        match = match.item()
    except:
        match = 0

    try:
        total_reads = df.loc['2-passQC', 'Reads']
        total_reads = total_reads.item()
    except:
        total_reads = 0

    perc_match = match / total_reads * 100
    
    return perc_match, total_reads

def process_nfflu_subfolder(subfolder_path):
    # Function to process each sample subdirectory
        # returns a dict with sample_id as well as coverage stats

    sample_name = os.path.basename(subfolder_path)
    
    # Initialise dict to contain sample values
        # Keys in this dict will be headers in the final output. They can be rearranged, except sample_id, which must remain first. 
        # Header names for segments here are simplified and homogenised (for example, both 'A_HA_H1' and 'A_HA_H3' will become simply 'HA')
        # Headers of the segments in this dict must match, or be substrings of, the segment names in the nf-flu output (e.g.:'A_HA_H1' )
         
    sample_dict = {
        'sample_id':'',
        'total_reads':0,
        '%_match':0,
        'total_coverage':0,
        'PB2':0,
        'PB1':0,
        'PA':0,
        'HA':0,
        'NP':0,
        'NA':0,
        'M':0,
        'NS':0
        }
    
    sample_dict_keys = list(sample_dict.keys())
    sample_dict['sample_id'] = sample_name

    # Define the path to the 'tables' folder within the current subfolder
    tables_folder = os.path.join(subfolder_path, 'tables')
    perc_match = 0
    total_reads = 0

    # Check if the 'tables' folder exists and if so, process
    if os.path.exists(tables_folder):
        # Process read counts data and set in dict
        read_counts_filepath = os.path.join(tables_folder, 'READ_COUNTS.txt')
        perc_match, total_reads = process_read_counts(read_counts_filepath)
        
        sample_dict['total_reads'] = total_reads
        sample_dict['%_match'] = perc_match
        
        # Iterate over all segment coverage files in the 'tables' subfolder and obtain the mean coverage for each segment, and add to sample dict
        # Mean coverage is the mean of the coverage at each genome position in the segment
        # Also calculates total_coverage, which is the mean of the coverage of each segment
        totals = [] 
        for file_path in glob.glob(os.path.join(tables_folder,'*-coverage.txt')):

            # Obtain the segment name from filename
            file_name = os.path.basename(file_path)
            segment_name = file_name.split('-')[0]
            
            # Rename the segment to match sample_dict keys (aka simplified segment names such as: 'A_HA_H1' or 'A_HA_H3' --> 'HA')
            for header in sample_dict_keys:
                if header in segment_name:
                    segment_name = header

            # Read the coverage file for given segment and append mean to dict
            df = pd.read_csv(file_path, sep='\t')
            mean = df['Coverage Depth'].mean()
            sample_dict[segment_name] = mean
            totals.append(mean)

        if len(totals) > 0:            
            sample_dict['total_coverage'] = sum(totals) / 8

    return sample_dict

def main(input, output_dir, output_basename, quiet):

    # Create list of subdirectories within the nf-flu results directory.
    IRMA_FOLDER = os.path.join(input, 'irma')
    subdirs = [ os.path.join(IRMA_FOLDER, item) for item in os.listdir(IRMA_FOLDER) if os.path.isdir(os.path.join(IRMA_FOLDER, item)) ]

    # Apply process_nfflu_subfolder function to each subdir (aka: each sample) and append to list
    results = []
    for subdir in subdirs:
        results.append(process_nfflu_subfolder(subdir))
    
    header_convert_dict= {
        'PB2':'PB2_coverage',
        'PB1':'PB1_coverage',
        'PA':'PA_coverage',
        'HA':'HA_coverage',
        'NP':'NP_coverage',
        'NA':'NA_coverage',
        'M':'M_coverage',
        'NS':'NS_coverage'
        }

    # Create a DataFrame from the results and replace any NaNs with 0
    coverage_stats_df = pd.DataFrame(results).fillna(0)
    coverage_stats_df = coverage_stats_df.rename(columns=header_convert_dict)

    # Convert df columns to int (except first col, which is str). This is done simply to remove any decimals.
    coverage_stats_df_headers = list(coverage_stats_df.keys())
    for header in coverage_stats_df_headers[1:]:
        coverage_stats_df[header] = coverage_stats_df[header].astype(int)

    # Write to file
    coverage_stats_df.to_csv(os.path.join(output_dir, output_basename), index=False)
    if not quiet:
        print(f'Wrote coverage stats for {coverage_stats_df.shape[0]} samples.\nFile saved here: {os.path.join(output_dir, output_basename)}')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='Collate some coverage statistics from a specified nf-flu analysis into a single csv file.\n\n \
                                     The output csv file will contain the following information:\n \
                                     \t\t\t"total_reads": the total reads which passed QC\n \
                                     \t\t\t"%_match": the number of reads which have matched against an Influenza reference\n \
                                     \t\t\t"total_coverage": the average of the coverage across all 8 segments\n \
                                     \t\t\t"PB2_coverage": the average coverage across the PB2 segment\n \
                                     \t\t\t"PB1_coverage": the average coverage across the PB1 segment\n \
                                     \t\t\t"PA_coverage": the average coverage across the PA segment\n \
                                     \t\t\t"HA_coverage": the average coverage across the HA segment\n \
                                     \t\t\t"NP_coverage": the average coverage across the NP segment\n \
                                     \t\t\t"NA_coverage": the average coverage across the NA segment\n \
                                     \t\t\t"M_coverage": the average coverage across the M segment\n \
                                     \t\t\t"NS_coverage": the average coverage across the NS segment\n\n \
                                     "This script expects that the nf-flu analysis results folder contains the "irma" subfolder and its \n \
                                     individual sample subfolders as they were after nf-flu analysis.')

    parser.add_argument('-i','--input', help='Path to the nf-flu output folder. This folder should contain the "irma" subfolder from the nf-flu analysis.', required=True)
    parser.add_argument('-o','--output_dir', help='Path specifiying where to save the output file. If none is given, output file will be saved in path specified \
                        in the input argument (e.g. the nf-flu output folder).')
    parser.add_argument('-b','--output_basename', help='Basename of the output filename. If none is given, default will be: "nf-flu_coverage_stats.csv"')
    parser.add_argument('-q','--quiet', action="store_true", help='Adding this flag will disable printing of stats to screen.')


    args = parser.parse_args()
    
    # Check if input path exists
    if not os.path.exists(args.input):
        sys.exit(f"Could not find the {args.input} folder, please check and try again.")

    # Set the output directory, and make one if it doesn't exist
    if args.output_dir is None:
        args.output_dir = args.input

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Set the basename
    if args.output_basename is None:
        args.output_basename = 'nf-flu_coverage_stats.csv'

    main(args.input, args.output_dir, args.output_basename, args.quiet)
