#!/bin/bash
#title          :section0.2_checksoftware.sh
#description    :check if all software packages are well installed.
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/section1.2_checksoftware.sh
#bash_version   :4.4.19(1)-release
#=======================================================================================================================

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
# # use package downloaded from github
# $TRIMGALORE_HOME/trim_galore
#
# # related utilities:
# cutadapt
# ####################  VERSION LOG  #############################
# 2.5
# ################################################################
#
# fastqc
# ####################  VERSION LOG  #############################
# FastQC v0.11.8
# ################################################################


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
# # download from github
# wget -c https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.2.1.tar.gz -O -| tar -xz
# # set the path to TRINOTATE_HOME
# echo 'export TRINOTATE_HOME=/media/lewis/Seagate_Backup_Plus_Drive/ruby/applications/Trinotate-Trinotate-v3.2.1' >> ~/.bashrc
# # if there are some modules missing in @INC, use miniconda3 to install them
# conda install -c bioconda perl-dbi
# conda install -c bioconda perl-dbd-sqlite
#
# # related utilities:
# # blastx
# ####################  VERSION LOG  #############################
# blastx: 2.9.0+
# Package: blast 2.9.0, build May 31 2019 20:53:30
# ################################################################
#
# # update_blastdb.pl --version
# ####################  VERSION LOG  #############################
# # /home/lewis/miniconda3/bin/update_blastdb.pl version 581818
# ################################################################
#
#
# # 8. "bioawk"
#
# # Description: extraction
#
# # Guidience to install:
# # use Miniconda3 to install
# conda install -c bioconda bioawk
#
# 9. "STAR"
#
# # Description: RNA-Seq read alignment tool to generate the bam file
#
# # Guidience to install:
# # use Miniconda3 to install
# conda install -c bioconda star





echo "WHERE IS ascp: `which ascp` "
echo "WHERE IS fastqc: `which fastqc`"
echo "WHERE IS multiqc: `which multiqc`"
echo "WHERE IS trim_galore: `which trim_galore`"
echo "WHERE IS Trinity: `which $TRINITY_HOME/Trinity`"
echo "WHERE IS bowtie2: `which bowtie2`"
echo "WHERE IS jellyfish: `which jellyfish`"
echo "WHERE IS salmon: `which salmon`"
echo "WHERE IS samtools: `which samtools`"
echo "WHERE IS transdecoder: `which transdecoder`"
echo "WHERE IS trinotate: `which trinotate`"
echo "If you miss some software packages, please go and check the INSTALLATION_GUIDENCE.txt"
echo "Also, You can check the version of applications I used ;)"
