# Bio-scripts
Some rough scripts to modify fasta datasets and to parse maf files

filter-reads.py can use data from a maf file to create a dataset of reads that have alignments. 
It can also be used to filter reads based on total sequence length.

alignment_stats.py will parse a maf file to record where alignments begin and end. It can also be used to filter alignments based on the total number of alignments found per read and based on read lengths.
