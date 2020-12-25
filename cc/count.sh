#!/bin/bash
mkdir ./processed_data
find . -type f -name *renamed.fq.gz -exec mv {} ./processed_data \;
find . -name sample_file_*.txt -exec mv {} ./processed_data \;
pro_name=`basename $(pwd)`
sshpass -f "/path/to/passwordfile" scp -r ./processed_data ruby@graham.computecanada.ca:
