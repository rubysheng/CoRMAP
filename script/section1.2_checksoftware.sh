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

# all needed software packages:
# 1. “ascp”:
# Description: for the fast datasets download
################################################################
# Aspera Connect version 3.9.6.173386
# ascp version 3.9.1.168302
# Operating System: Linux
# FIPS 140-2-validated crypto ready to configure
# AES-NI Supported
# License max rate=(unlimited), account no.=1, license no.=1
################################################################
ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.sh
# open the the website to find lastest version of 'aspera connnect' http://downloads.asperasoft.com/en/downloads/8?list
wget -qO- https://download.asperasoft.com/download/sw/connect/3.9.6/ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz | tar xvz

## run it
chmod +x ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh
./ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh

## add it to the path now and in the future
export PATH=$PATH:~/.aspera/connect/bin/
echo 'export PATH=$PATH:~/.aspera/connect/bin/' >> ~/.bash_profile

# “fastqc”,
# “multiqc”,
# “trim_galore”,
# “trinity”,
# “transdecoder”,
# “trinotate”
