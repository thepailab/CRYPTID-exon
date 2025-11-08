# CRYPTID-exon
# Overview

CRYPTID-exon (CRYPTic IDentification of exons) is a python framework that predicts exon boundaries using local changes in read coverage in high-throughput short-read RNA-seq data around empirically derived splice sites.
Using mapped reads (BAM file) and a text file containing coordinates of splice junctions as inputs, CRYPTID-exon proceeds with these steps:
- Get read coverage within a user-defined window surrounding 3'splice sites and 5' splice sites
- Obtain putative exon boundaries by predicting changepoints within the read coverage window
- Generate a list of high-confidence exons by matching predictions from 3'splice sites and 5' splice sites

# Prerequisites (versions used for development)
- Python 3.11.0
- Bedtools 2.28.0
- Pandas 2.0.0
- Numpy 1.26.4
- Ruptures 1.1.9
  
# Installation
Download and install the prerequisite packages. Download cryptidEX_main.py

# Usage 
```
cryptidExon_main.py -splice_file .csv -bam .bam -output_dir output_path

Arguments:
-splice_file Output file from CRYPTID-SS or text file containing splice junctions (.csv). If providing a custom
             file, first 5 columns should follow this format:  unique identifier of junctions, chromosome, 5'splice site coordinate, 3' splice site coordinate, strand
-bam path to BAM file containing mapped reads. 
-output_dir path to the directory where output files will be generated 

Optional arguments:
-intron_ window  size of coverage window extended into the intron from exon boundary. Default is 100 nt
-exon_ window  size of coverage window extended into the intron from exon boundary. Default is 200 nt


```
