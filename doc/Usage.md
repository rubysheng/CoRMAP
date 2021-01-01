# Usage the CMRP

Use the documentation links below to navigate this documentation.

## Table of contents

-   [Usage the CMRP](#usage-the-cmrp)
    -   [Table of contents](#table-of-contents)
    -   [Overview of CMRP](#overview-of-cmrp)
    -   [Prepare](#prepare)
        -   [Set up the structure of folder](#set-up-the-structure-of-folder)
    -   [Download datasets](#download-datasets)
    -   [Clean data](#clean-data)
    -   [Obtain Assembly](#obtain-assembly)
    -   [Quantify gene expression](#quantify-gene-expression)
    -   [Predict peptides sequences](#predict-peptides-sequences)
    -   [Search Orthologs](#search-orthologs)
        -   [Integrate across taxonomy groups](#integrate-across-taxonomy-groups)
        -   [Annotate orthologous groups](#annotate-orthologous-groups)
        -   [Compare expression of orthologous groups](#compare-expression-of-orthologous-groups)

## Overview of CMRP

The overview of processing workflow in CMRP:

![workflow](https://github.com/rubysheng/CMRP/blob/bombina/pics/pipeline.png)

[Back to top](#usage-the-cmrp)

## Prepare

### Set up the structure of folder

For the data management, create those suggested folder structure by using the script ___script/section0.1_setup.sh___ with a list of meta-dataset names (check ___sample/list_of_projectnames___) or manage your input datasets as the following structure:

![folder-structure](https://github.com/rubysheng/CMRP/blob/bombina/pics/dir_structure.png)

[Back to top](#usage-the-cmrp)

## Download datasets

CMRP includes a small utility ___script/section0.3_download.sh___ for downloading raw data from published RNA-Seq experiments by their database accession numbers (check ___sample/list_of_projectnames___).

Datasets containing multiple runs were downloaded from the Sequence Read Archive (SRA) through a list of FTP links, obtained from the European Nucleotide Archive (ENA).

[Ascp](https://www.ncbi.nlm.nih.gov/books/NBK158898/) is a tool allowing to quickly download a batch of data by links.

The example input text files for `Ascp` were provided under the ___sample/download_input/___ folder.

[Back to top](#usage-the-cmrp)

## Clean data

For every single dataset, go to its own directory under ___./data/___, and run with the pipeline.

Raw data and cleaned data quality checks were performed using the [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc) and [MultiQC](https://multiqc.info/). `MultiQC` is useful to generate the combined reports of quality check for several runs.

RNA-Seq data quality control using [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) included low-quality base calls (using default Phred score 20) trimming, adaptor auto-detection and trimming, and short reads filtering (the default minimum length cut-off was 20 bp).

Check the example output [reports](https://github.com/rubysheng/CMRP/blob/bombina/sample/output/qc_report) after quality control.

[Back to top](#usage-the-cmrp)

## Obtain Assembly

To reduce the number of input files and decrease the computational complexity and processing time, data normalization was performed using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) separately from and before de novo assembly. The targeted maximum coverage for reads was 50, while the Kmer size and the maximum percent of the mean for the standard deviation of Kmer coverage across reads were both using the default value, 25 and 200, respectively.

De novo assembly with Trinity was carried out in [Galaxy](http://usegalaxy.org) using default parameter settings. Assemblies were assessed by basic contig N50 statistics, which specifies the length of the shortest contig that can cover 50% of the total genome length.  


## Quantify gene expression

Transcript mapping and quantifying back to the assembly were performed using the plugin in `Trinity`. RNA-seq by Expectation-Maximization (RSEM) was used to estimate the contigs abundance in an alignment-based algorithm with corrections for transcript lengths. Briefly, RESM assumes that each sequence read and read length are observed (one end for single reads and both ends for paired reads) and are generated according to a directed acyclic graph that includes parent sequence, length, start position and orientation.  A Bayesian version of the expectation-maximization algorithm is used to obtain maximum likelihood estimates of the model parameters, transcript fractions and posterior mean estimates of abundances. The resulting gene expression matrices were then normalized to Transcripts Per Kilobase Million (TPM) to make them comparable across samples. The original TPM matrix generated from the Trinity plugin included all gene names in each dataset. For each dataset, the expression matrix was transformed into a long format with gene names, sample names and TPM values.

[Back to top](#usage-the-cmrp)

## Predict peptides sequences

[TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki) was used to predict coding regions from transcripts, and consequently provided the amino acid sequences used in further analysis.

[Back to top](#usage-the-cmrp)

## Search Orthologs

Orthologous gene groups were clustered using an automated pipeline for [OrthoMCL](https://github.com/stajichlab/OrthoMCL), called [Orthomcl-pipeline](https://github.com/apetkau/orthomcl-pipeline).

In this tool, all-versus-all Blast search with Opening Reading Frames provided the score of pairwise sequence similarity. Blast used the default setting except for E-value cut-off, which was different based on the diversity of species (or higher-level taxonomy groups). Then all matching pairs were filtered with “percent match length” (>= 50% left) and linked across or within the objects (species or class) for identifying the orthologs or in-paralogs pairs. Finally, the clustering of orthologs gene groups used Markov Clustering (MCL) algorithms.

### Integrate across taxonomy groups

To allow analyzing data from a wide range of species, orthologs searching was required to identify the orthologous gene groups. For large amount of orthologs searching among distantly related species, datasets would be split into several groups based on the taxonomy group of their experiment objects.

Orthologs searching was performed at two levels- clustering intraclass and clustering interclass. The general process involved 1) in-parallel orthologs searching within each class; 2) filtering clustered genes from those OGGs that contained orthologous genes from all species in one class; 3) integrating filtered genes from 5 classes for interclass orthologs searching; and 4) filtering clustered genes from those OGGs that contained orthologous genes from all classes.

If the input datasets had to conduct with the interclass orthologs searching followed by some filtering, the cross-class OGGs would be obtained, defined as orthologs groups (or cOGGs).

[Back to top](#usage-the-cmrp)

### Annotate orthologous groups

From each orthologs group (or cOGG), only all clustered Mammalia sequences would be extracted for annotation because they had much more well-developed and high-quality gene annotation database. Annotation of Mammalia proteins (using BlastP) and transcripts (using BlastX) were against the UniProt database obtained from [Trinotate](http://trinotate.github.io).

The E-value cut-off for both were 10-3; BlastP only kept the most common target, while BlastX returned the top 5 hits; other parameters were set as defaults except for the parallel threads using the maximum cores. Protein sequences were also searched against a protein profile HMM database.

Trinotate generated the final report. After annotation, only the most frequent gene symbol from BlastX and annotated with Gene Ontology (GO) terms in biological processes for each orthologs group (or cOGG) was selected as the representative annotation identifier. When duplicate identifiers occurred, the second most frequent gene symbol was applied to the orthologs group (or cOGG) with the most frequency of its second-most gene symbol.

[Back to top](#usage-the-cmrp)

### Compare expression of orthologous groups

[Back to top](#usage-the-cmrp)
