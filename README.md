# Comparative Meta RNA-Seq Data Standardized Analysis Pipeline
Comparative Meta RNA-Seq Data Standardized Analysis Pipeline (CMRP) is a processing frame for the standardized analysis of Meta RNA-Seq raw data from wide-ranged species.

## Motivation
The development of RNA sequencing (RNA-Seq) techniques has provided efficient methods to analyze whole transcriptomes, allowing for a full detection of the transcriptomic response of organisms under various conditions. Comparative transcriptomics analysis insights of comparaing gene expression across species, while meta-analysis provides methods of comparative analysis cross different studies. Here we defined a pipeline, CMRP, using the raw data of RNA-Seq datasets from published papers, processed with constant protocols for various species, and finally conducted with two integration methods of differential expression analysis.  

Created by Yiru Sheng, Andreas Heyland, and Ayesha Ali at the University of Guelph, Canada. Released under the terms of the General Public License, version 3.0 (GPLv3).

## Installation
See the [installation](https://github.com/rubysheng/CMRP/blob/bombina/doc/Install.md) documentation.

## Usage
To use CMRP, run each section of code step by step (recommended). For example:

  $ source script/section1.1_download_datasets.sh

For quick usage, run the automatically processing command:

  $ source CMRP.sh

More documents about the usage of all sections, please check [here](https://github.com/rubysheng/CMRP/blob/bombina/doc/Usage.md).
