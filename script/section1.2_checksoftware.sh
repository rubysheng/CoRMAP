#!/bin/bash
#title          :section1.2_checksoftware.sh
#description    :check if all software packages are well installed.
#author         :Ruby(Yiru) Sheng
#date           :20191016
#version        :1.0
#usage          :./section1.2_checksoftware.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

source ${scriptloc}/section1.1_setup_directory.sh
# cd ${maindirectory}/applications
#
# # all needed software packages:
#
# # 1. “ascp”:
#
# # Description: for the fast datasets download
# ####################  VERSION LOG  #############################
# # Aspera Connect version 3.9.6.173386
# # ascp version 3.9.1.168302
# # Operating System: Linux
# # FIPS 140-2-validated crypto ready to configure
# # AES-NI Supported
# # License max rate=(unlimited), account no.=1, license no.=1
# ################################################################
#
# # Guidience to install:
# #   open the the website to find lastest version of 'aspera connnect' http://downloads.asperasoft.com/en/downloads/8?list
# wget -qO- https://download.asperasoft.com/download/sw/connect/3.9.6/ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz | tar xvz
#
# #   run it
# chmod 777 ibm-aspera-connect-*.sh
# ./ibm-aspera-connect-*.sh
#
# #  add it to the path for using it more simplely in the future
# export PATH=$PATH:~/.aspera/connect/bin/
# echo 'export PATH=$PATH:~/.aspera/connect/bin/' >> ~/.bashrc
#
# # 2. “fastqc”:
#
# # Description: for RNA sequencing data quality check
# ####################  VERSION LOG  #############################
# # FastQC v0.11.8
# ################################################################
#
# # Guidience to install:
# # use Miniconda3 to install
# conda install -c bioconda fastqc
#
# # 3. “multiqc”:
#
# # Description: to combine multiple quality check reports to one
# ####################  VERSION LOG  #############################
# # /home/lewis/miniconda3/envs/multiqc/lib/python3.6/site-packages/multiqc/utils/config.py:45: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
# #   configs = yaml.load(f)
# # /home/lewis/miniconda3/envs/multiqc/lib/python3.6/site-packages/multiqc/utils/config.py:51: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
# #   sp = yaml.load(f)
# # multiqc, version 1.7
# ################################################################
#
# # Guidience to install:
# # use Miniconda3 to install
# conda config --add channels conda-forge
# conda config --add channels bioconda
# conda create -n multiqc multiqc
#
# # 4. “trim_galore”
#
# # Description: for RNA sequencing data quality control
# ####################  VERSION LOG  #############################
# #   Quality-/Adapter-/RRBS-/Speciality-Trimming
# #              [powered by Cutadapt]
# #                  version 0.6.4
# #
# #              Last update: 31 07 2019
# ################################################################
#
# # Guidience to install:
# # use Miniconda3 to install
# conda install -c bioconda trim_galore
#
# # 5. “trinity”
#
# # Description: normalization, de novo assembly, expression matrix generating, etc.
# ####################  VERSION LOG  #############################
# # Trinity version: Trinity-v2.8.5
# # ** NOTE: Latest version of Trinity is v2.8.6, and can be obtained at:
# # 	https://github.com/trinityrnaseq/trinityrnaseq/releases
# ################################################################
#
# # Guidience to install:
# # check the lastest release https://github.com/trinityrnaseq/trinityrnaseq/releases and download
# wget -qO- https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.8.5.tar.gz | tar xvz
# # follow the guidence of Trinityseq/github https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity
#
# # related utilities:
# # bowtie2
# ####################  VERSION LOG  #############################
# # /home/lewis/Bowtie2/bowtie2-2.3.5.1-linux-x86_64/bowtie2-align-s version 2.3.5.1
# # 64-bit
# # Built on
# # Wed Apr 17 02:50:12 UTC 2019
# # Compiler: gcc version 7.3.1 20180303 (Red Hat 7.3.1-5) (GCC)
# # Options: -O3 -m64 -msse2 -funroll-loops -g3 -g -O2 -fvisibility=hidden -I/hbb_exe_gc_hardened/include -ffunction-sections -fdata-sections -fstack-protector -D_FORTIFY_SOURCE=2 -fPIE -std=c++98 -DPOPCNT_CAPABILITY -DWITH_TBB -DNO_SPINLOCK -DWITH_QUEUELOCK=1
# # Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
# ################################################################
# #
# # jellyfish
# ####################  VERSION LOG  #############################
# # jellyfish 2.2.4
# ################################################################
# #
# # salmon
# ####################  VERSION LOG  #############################
# # salmon 0.14.1
# #   installation guide:
# #     conda config --add channels conda-forge
# #     conda config --add channels bioconda
# #     conda create -n salmon salmon
# ################################################################
# #
# # samtools
# ####################  VERSION LOG  #############################
# # samtools 1.9
# # Using htslib 1.9
# # Copyright (C) 2018 Genome Research Ltd.
# ################################################################
# #
# # rsem-calculate-expression
# ####################  VERSION LOG  #############################
# # Current version: RSEM v1.3.1
# ################################################################
#
#
# # 6. “transdecoder”
#
# # Description: annotation
# ####################  VERSION LOG  #############################
# # FastQC v0.11.8
# ################################################################
#
# # Guidience to install:
# # use Miniconda3 to install
# conda install -c bioconda transdecoder
#
# # 7. “trinotate”
#
# # Description: annotation
# ####################  VERSION LOG  #############################
# # FastQC v0.11.8
# ################################################################
#
# # Guidience to install:
# # use Miniconda3 to install
# conda install -c bioconda trinotate
