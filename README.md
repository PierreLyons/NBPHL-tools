# NBPHL-tools
An initial repository for keeping and sharing simple scripts and tools mainly used to manipulate and organise pathogen WGS pipeline output data.


## nf-flu_cov_stats_collator.py

This script will go through the output of a nf-flu run and collate some coverage stats into a csv table. 

Usage example: python nf-flu_cov_stats_collator.py -i path/to/nf-flu/output

To see other options: python nf-flu_cov_stats_collator.py -h

## nf-flu_segment_fasta_collator.py

This script will go through the output of a nf-flu run and collate segment fasta files into one or more multifastas. These can then be used as input for nextclade. 

Usage example: python nf-flu_segment_fasta_collator.py -i path/to/nf-flu/output

To see other options: python nf-flu_segment_fasta_collator.py -h
