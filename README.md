# Comparative Meta RNA-Seq Data Standardized Analysis Pipeline
Comparative Meta RNA-Seq Data Standardized Analysis Pipeline (CMRP) is a processing frame for the standardized analysis of Meta RNA-Seq raw data from wide-ranged species.

## Motivation
The development of RNA sequencing (RNA-Seq) techniques has provided efficient methods to analyze whole transcriptomes, allowing for a full detection of the transcriptomic response of organisms under various conditions. Comparative transcriptomics analysis insights of comparaing gene expression across species, while meta-analysis provides methods of comparative analysis cross different studies. Here we defined a pipeline, CMRP, using the raw data of RNA-Seq datasets from published papers, processed with constant protocols for various species, and finally conducted with two integration methods of differential expression analysis.  

Created by Yiru Sheng, Andreas Heyland, and Ayesha Ali at the University of Guelph, Canada. Released under the terms of the General Public License, version 3.0 (GPLv3).

## Installation
See the [installation]() documentation.


CMRP
Trinity assembles transcript sequences from Illumina RNA-Seq data.

Download Trinity here.

Build Trinity by typing 'make' in the base installation directory.

Assemble RNA-Seq data like so:

 Trinity --seqType fq --left reads_1.fq --right reads_2.fq --CPU 6 --max_memory 20G
Find assembled transcripts as: 'trinity_out_dir/Trinity.fasta'

Use the documentation links in the right-sidebar to navigate this documentation, and contact our Google group for technical support.
