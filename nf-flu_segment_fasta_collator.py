import os
import sys
import argparse
import glob


def make_fasta_segments_multifasta(segments, sample_paths_list, flu_type, output_dir, output_basename, quiet):
        for segment in segments:
            fasta_list = []

            if flu_type == 'mix':
                search_pattern = f'*{segment.upper()}*.fasta'
                outfile_prefix = ''
            else:
                search_pattern = f'*{flu_type.upper()}_{segment.upper()}*.fasta'
                outfile_prefix = f'{flu_type.upper()}_' 
            
            for sample_path in sample_paths_list:
                sample_name = os.path.basename(sample_path)

                fasta_paths = glob.glob(os.path.join(sample_path, search_pattern))

                if len(fasta_paths) > 0:
                    for fasta_path in fasta_paths:
                        with open(fasta_path, 'r') as f:
                            fasta_list.append(f'>{sample_name}\n')
                            fasta_list.append(f.readlines()[1])
                """ else:
                    if not quiet:
                        print(f'No fasta files found for Influenza {flu_type.upper()} {segment.upper()} segment in sample {sample_name}.') """
    
            output_filename = f'{outfile_prefix}{segment.upper()}{output_basename}'
            with open(os.path.join(output_dir, output_filename), 'w') as fp:
                fp.writelines(fasta_list)
            
            # Print some stats to screen if -q/--quiet not provided at command line.
            if not quiet:

                plural = ''
                if (len(fasta_list) / 2) > 1 or (len(fasta_list) / 2) == 0:
                    plural = 's'

                if flu_type == 'mix':
                    print(f'Added {len(fasta_list) / 2:.0f} {segment.upper()} segment{plural} to file: {os.path.join(output_dir, output_filename)}') 
                else:
                    print(f'Added {len(fasta_list) / 2:.0f} Influenza {flu_type.upper()} {segment.upper()} segment{plural} to file: {os.path.join(output_dir, output_filename)}') 




def main(input_path, flu_type, segments, output_dir, output_basename, quiet):

    # Iterate over all subfolders in the root folder
    IRMA_FOLDER = os.path.join(input_path,'irma')

    #Make a list of all subdirectories in the irma folder. Each subdir represents a unique sample
    sample_paths_list = [ os.path.join(IRMA_FOLDER, item) for item in os.listdir(IRMA_FOLDER) if os.path.isdir(os.path.join(IRMA_FOLDER, item)) ]

    if flu_type == "separate":
        make_fasta_segments_multifasta(segments, sample_paths_list, 'a', output_dir, output_basename, quiet)
        make_fasta_segments_multifasta(segments, sample_paths_list, 'b', output_dir, output_basename, quiet)
    else:
        make_fasta_segments_multifasta(segments, sample_paths_list, flu_type, output_dir, output_basename, quiet)
               


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collate all fasta segments from a specific nf-flu analysis into either a single multifasta file \
                                     or separate multifasta files.\n\
                                     This script expects that the nf-flu analysis results folder contains the "irma" subfolder and its \
                                     individual sample subfolders as they were after nf-flu analysis.')

    parser.add_argument('-i','--input', help='Path to the nf-flu output folder. This folder should contain the "irma" subfolder from the nf-flu analysis.', required=True)
    parser.add_argument('-t','--types', help='Specify which flu types you wish to collate. Options are: "a", "b", "mix", or "separate" (either lowercase or UPPERCASE). \
                        Default is "separate", which will generate separate multifastas for both Influenza A and Influenza B. \
                        "mix" will collate both types into the same multifasta. "a" or "b" will collate only the segments of the specified Influenza type.')
    parser.add_argument('-s','--segments', nargs='+', help='Specify which segments you wish to collate into a multifasta. \
                        Options are: "pb2","pb1","pa","ha","np","na","m","ns", or "all".\
                         A separate multifasta file will be created for each segment specified.\
                         If not specified, will default to "ha".')
    parser.add_argument('-o','--output_dir', help='Path specifiying where to save the output file. If none is given, output file will be saved in path specified \
                        in the input argument (e.g. the nf-flu output folder).')
    parser.add_argument('-b','--output_basename', help='Basename of the output filename. If none is given, default will be "_segments.fasta", therefore output file will \
                        look like: "[type]_[segment]_segments.fasta", unless "mix" is set for the type, in which case output file will be: [segment]_segments.fasta ')
    parser.add_argument('-q','--quiet', action="store_true", help='Adding this flag will disable printing of stats to screen.')



    args = parser.parse_args()

    # Check if input path exists
    if not os.path.exists(args.input):
        sys.exit(f"Could not find the {args.input} folder, please check and try again.")

    # Set the type
    if args.types is None:
        flu_type = 'separate'
    else:
        flu_type = args.types.lower()
        type_options = ['a','b','mix','separate']
      
        while flu_type not in type_options:
            input(f'You must specify a valid flu type to collate, options are: "a","b", "mix", or "separate". You entered: {args.types}\nPlease try again: ')
    
    # Set the segments to be collated
    if args.segments is None:
        segments = ['ha']
    
    else:
        segments = [s.lower() for s in args.segments]
    
        if segments == ['all']:
            segments = ['pb2','pb1','pa','ha','np','na','m','ns']
    
        else:
            segments = [s for s in segments if s != 'all']
            segment_options = ['pb2','pb1','pa','ha','np','na','m','ns']
            for s in segments:
                if s not in segment_options:
                    sys.exit(f'You entered an invalid segment: {s}\nOptions are: "pb2","pb1","pa","ha","np","na","m","ns", or "all".\nPlease try again.')


    
    
    if not isinstance(segments, list):
        segments = [segments]
        

    
    # Set the output directory, and make one if it doesn't exist
    if args.output_dir is None:
        args.output_dir = args.input

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Set the basename
    if args.output_basename is None:
        args.output_basename = '_segments.fasta'
    

    # Run the main function
    main(args.input, flu_type, segments, args.output_dir, args.output_basename, args.quiet)
