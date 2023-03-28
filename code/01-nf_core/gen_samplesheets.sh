#!/bin/bash

wget -L https://raw.githubusercontent.com/nf-core/rnaseq/master/bin/fastq_dir_to_samplesheet.py

# for each folder in /mnt/isilon/cscb/data/BaseSpace/HTS019-369865496/HTS019.1-369865496/BCL_Convert_10_20_2022_10_37_35-620098007/, run the following command

# variable set to 1
q=0

for i in /mnt/isilon/cscb/data/BaseSpace/HTS019-369865496/HTS019.1-369865496/BCL_Convert_10_20_2022_10_37_35-620098007/*; 
do 
    q=$((q+1));
    python fastq_dir_to_samplesheet.py $i samplesheet_$q.csv \
    --strandedness forward; 
done

# concatenate all samplesheets into one file using a for loop ranging from 1 - 16
for i in {1..16};
do
    # if i > 1, omit first line
    if [ $i -gt 1 ]; then
        tail -n +2 samplesheet_$i.csv >> samplesheet_main.csv;
    else
        cat samplesheet_$i.csv >> samplesheet_main.csv;
    fi
done