# README #

C++ code for barcoding-based consensus building tool.

Usage: ./build-consensus -i <input_bam> -o <output_prefix> -l <interval_list> -r <ref_fasta> -d <dbsnp_file> \

 [-b min_base_qual] [-m min_map_qual] [-U UMI Tag] [-e max_umi_edits] [-n consensus_unique_count] \
 
 [-c min_consensus_size] [-s use_singleton_reads] [-p padding_length] [-v]

Parameter Description:

 -i Path to input bam file
 
 -o Prefix to output files including directory path (Will be used for consensus bam and coverage file)
 
 -l Interval file in bed format

 -r Path to reference fasta file
 
 -d Path to dbSNP file
 
 -b Minimum base quality to consider for consensus building [Default 20]
 
 -m Minimum mapping quality to consider for consensus building [Default 30]
 
 -U UMI Tag name [Default BC]
 
 -e Max edits permitted in UMI for consensus building [Default 2]
 
 -n Number of consensus reads in a family to be marked as unique. Used to downweight singleton reads [Default 1]
 
 -c Minimum number of reads in a family to call consensus [Default 2]
 
 -s Use singleton reads towards coverage and variant calling [Default true]
 
 -v Call Variants from consensus building. Switched off for now [Default false]
 
 -p Pad intervals for mutation calling [Default 0]
 


The consensus building tool accepts a fully aligned and BQSR bam file as input and the interval list, and prepares consensus reads for further variant calling. The output bam file can be fed to a variant caller such as MuTect to call mutations on a clean set of consensus reads. Usage is specified above.
