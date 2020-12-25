#!/bin/bash
#title          :change_dir_samplefile.sh
#description    :change the dir name in all new sample files
#author         :Ruby(Yiru) Sheng
#date           :20191030
#version        :1.0
#usage          :./change_dir_samplefile.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================
sed "s/SE/SR/g" ./data/PRJNA201010/sample_file_1.txt
sed "s/SE/SR/g" ./data/PRJNA201386/sample_file_2.txt
sed "s/SE/SR/g" ./data/PRJNA390522/sample_file_8.txt
sed "s/SE/SR/g" ./data/PRJNA450614/sample_file_10.txt
sed "s/SE/SR/g" ./data/PRJNA529794/sample_file_12.txt
