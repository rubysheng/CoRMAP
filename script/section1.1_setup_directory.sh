#!/bin/bash
#title          :section1.1_setup_directory.sh
#description    :Here is the very beginning to setup the main work directory, and build up folder structure for data analysis
#author         :Ruby(Yiru) Sheng
#date           :20191016
#version        :1.1
#usage          :./section1.1_setup_directory.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

# set the main work directory that contains all datasets
maindirectory=$(pwd)
scriptloc='comparative-transcriptomic-analysis-pip/script/'

mkdir applications datasets
