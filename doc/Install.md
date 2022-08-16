# Installing the CoRMAP

Installing the Comparative Transcriptomics Analysis Pipeline (CoRMAP) can be accomplished by downloading the code with the following command and then following the steps below.

    $ git clone https://github.com/rubysheng/CoRMAP.git

<!-- TOC -->

- [Installing the CoRMAP](#installing-the-cormap)
    - [Requirement](#requirement)
        - [Ascp](#ascp)
        - [Fastqc](#fastqc)
        - [Multiqc](#multiqc)
        - [Trim Galore](#trim-galore)
            - [Cutadapt](#cutadapt)
            - [Fastqc](#fastqc-1)
        - [Trinity](#trinity)
            - [Bowtie2](#bowtie2)
            - [Jellyfish](#jellyfish)
            - [Salmon](#salmon)
            - [Samtools](#samtools)
            - [Rsem-calculate-expression](#rsem-calculate-expression)
        - [Transdecoder](#transdecoder)
        - [Trinotate](#trinotate)
            - [BlastX](#blastx)
            - [Update_blastdb.pl](#update_blastdbpl)
        - [Bioawk](#bioawk)
        - [STAR](#star)
        - [Orthomcl-pipeline](#orthomcl-pipeline)
    - [Notes](#notes)

<!-- /TOC -->

## Requirement

### Ascp

Description: for the fast datasets download

```bash
####################  VERSION LOG  #############################
# Aspera Connect version 3.9.6.173386
# ascp version 3.9.1.168302
# Operating System: Linux
# FIPS 140-2-validated crypto ready to configure
# AES-NI Supported
# License max rate=(unlimited), account no.=1, license no.=1
################################################################
```

**Guidience to install:**

1.  open the the website to download the specific version of [aspera connnect](http://downloads.asperasoft.com/en/downloads/8?list)

        $ wget -qO- https://download.asperasoft.com/download/sw/connect/3.9.6/ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz | tar xvz

2.  installation

        $ chmod 777 ibm-aspera-connect-\*.sh
        $ ./ibm-aspera-connect-\*.sh

3.  add it to the path for using it more simplely in the future

        $ export PATH=$PATH:~/.aspera/connect/bin/
        $ echo 'export PATH=$PATH:~/.aspera/connect/bin/' >> ~/.bashrc

### Fastqc

Description: for RNA sequencing data quality check

```bash
####################  VERSION LOG  #############################
# FastQC v0.11.8
################################################################
```

**Guidience to install:**

1.  use Miniconda3 to install

        $ conda install -c bioconda fastqc

### Multiqc

Description: to combine multiple quality check reports to one

```bash
####################  VERSION LOG  #############################
# /path/to/miniconda3/envs/multiqc/lib/python3.6/site-packages/multiqc/utils/config.py:45: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
#   configs = yaml.load(f)
# /path/to/miniconda3/envs/multiqc/lib/python3.6/site-packages/multiqc/utils/config.py:51: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
#   sp = yaml.load(f)
# multiqc, version 1.7
################################################################
```

**Guidience to install:**

1.  use Miniconda3 to install

        $ conda config --add channels conda-forge
        $ conda config --add channels bioconda
        $ conda create -n multiqc multiqc

### Trim Galore

Description: for RNA sequencing data quality control

```bash
####################  VERSION LOG  #############################
#   Quality-/Adapter-/RRBS-/Speciality-Trimming
#              [powered by Cutadapt]
#                  version 0.6.4
#
#              Last update: 31 07 2019
################################################################
```

**Guidience to install:**

1.  use package downloaded from GitHub

        $ TRIMGALORE_HOME/trim_galore

#### Cutadapt

```bash
####################  VERSION LOG  #############################
# 2.5
################################################################
```

#### Fastqc

```bash
####################  VERSION LOG  #############################
# FastQC v0.11.8
################################################################
```

### Trinity

Description: normalization, de novo assembly, expression matrix generating, etc.

```bash
####################  VERSION LOG  #############################
# Trinity version: Trinity-v2.8.5
# ** NOTE: Latest version of Trinity is v2.8.6, and can be obtained at:
# 	https://github.com/trinityrnaseq/trinityrnaseq/releases
################################################################
```

**Guidience to install:**

1.  check the specific version of [release](https://github.com/trinityrnaseq/trinityrnaseq/releases) and download

        $ wget -qO- https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.8.5.tar.gz | tar xvz

2.  follow the guidence of [Trinityseq on GitHub](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity)

#### Bowtie2

```bash
####################  VERSION LOG  #############################
# /home/lewis/Bowtie2/bowtie2-2.3.5.1-linux-x86_64/bowtie2-align-s version 2.3.5.1
# 64-bit
# Built on
# Wed Apr 17 02:50:12 UTC 2019
# Compiler: gcc version 7.3.1 20180303 (Red Hat 7.3.1-5) (GCC)
# Options: -O3 -m64 -msse2 -funroll-loops -g3 -g -O2 -fvisibility=hidden -I/hbb_exe_gc_hardened/include -ffunction-sections -fdata-sections -fstack-protector -D_FORTIFY_SOURCE=2 -fPIE -std=c++98 -DPOPCNT_CAPABILITY -DWITH_TBB -DNO_SPINLOCK -DWITH_QUEUELOCK=1
# Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
################################################################
```

#### Jellyfish

```bash
####################  VERSION LOG  #############################
# jellyfish 2.2.4
################################################################
```

#### Salmon

```bash
####################  VERSION LOG  #############################
# salmon 0.14.1
#   installation guide:
#     conda config --add channels conda-forge
#     conda config --add channels bioconda
#     conda create -n salmon salmon
################################################################
```

#### Samtools

```bash
####################  VERSION LOG  #############################
# samtools 1.9
# Using htslib 1.9
# Copyright (C) 2018 Genome Research Ltd.
################################################################
```

#### Rsem-calculate-expression

```bash
####################  VERSION LOG  #############################
# Current version: RSEM v1.3.1
################################################################
```

### Transdecoder

Description: annotation

**Guidience to install:**

1.  use Miniconda3 to install

        $ conda install -c bioconda transdecoder

### Trinotate

Description: annotation

```bash
####################  VERSION LOG  #############################
# FastQC v0.11.8
################################################################
```

**Guidience to install:**

1.  download from GitHub

        $ wget -c https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.2.1.tar.gz -O -| tar -xz

2.  set the path to TRINOTATE_HOME

        $ echo 'export TRINOTATE_HOME=/path/to/Trinotate-Trinotate-v3.2.1' >> ~/.bashrc

3.  if there are some modules missing in `@INC`, use miniconda3 to install them

        $ conda install -c bioconda perl-dbi
        $ conda install -c bioconda perl-dbd-sqlite

#### BlastX

```bash
####################  VERSION LOG  #############################
# blastx: 2.9.0+
# Package: blast 2.9.0, build May 31 2019 20:53:30
################################################################
```

#### Update_blastdb.pl

```bash
####################  VERSION LOG  #############################
# /path/to/miniconda3/bin/update_blastdb.pl version 581818
################################################################
```

### Bioawk

Description: extraction

**Guidience to install:**

1.  use Miniconda3 to install
        $ conda install -c bioconda bioawk

### STAR

Description: RNA-Seq read alignment tool to generate the bam file

**Guidience to install:**

1.  use Miniconda3 to install
        $ conda install -c bioconda star

### Orthomcl-pipeline

Description: search orthologs gene groups

**Guidience to install:**

1.  follow the instruction for [Orthomcl-pipeline](https://github.com/apetkau/orthomcl-pipeline/blob/master/INSTALL.md)

## Notes

Errors may occur while running software included in CoRMAP due to using different version of software.

After installation, please check the path of those software have been added to the `$PATH` environment variable.
