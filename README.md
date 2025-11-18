# Alport_Syndome_project
Alport Symdrome Chinese cohort project focused on well-known causal genes COL4A3/4/5
Author: Zhen Y

### contents
#### test_data: 
folder contains genotype-phenotype data for AS patients with hearing loss, cysts, and ESKD;
#### pheno-geno: 
folder contains scripts to perform genotype-phenotype association analysis between genotypes and hearing loss, cysts, and ESKD. Use data from test_data as input.
#### TGS_Figure_QC_R_plot.R
R script to visualize ONT nanopore WGS data quality, including read length distribution and reads quality distribution
#### extractMatrix_bedGraph.py
python script use to extract individual signals of overlap genomic intervals from multiple bedGraph files; ussage:
extractMatrix_bedGraph.py --path /input/dir --pattern bedGraph_file_prefix --output output_file
#### extract_clipped_seq.py
extract insersion sequence from regional TGS bam file
Usage: python extract_clipped_seq.py insertion_bearing_regional.bam --min-length 600 --fasta -o inserted_reads.fa
