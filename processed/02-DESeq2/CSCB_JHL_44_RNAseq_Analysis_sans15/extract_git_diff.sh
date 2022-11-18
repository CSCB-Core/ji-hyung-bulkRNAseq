#!/bin/bash

# 'dc' is the ear
# git diff dc915f7cda6b63921766b3a95c02f4607517876b 5b24014bb4e08b91c7026afe4e12d4ebda363215 genomeKOvsWT.conditioncondition2.csv > gitdiff_genomeKOvsWT.conditioncondition2.csv

# # add comma after every '-' and '+' in the diff file
# sed -i 's/^-/-,/g' gitdiff_genomeKOvsWT.conditioncondition2.csv

# for every file that matches the patter *conditioncondition*.csv
for file in genome*conditioncondition*.csv
do
    # get the file name without the extension
    filename=$(basename "$file")
    extension="${filename##*.}"
    filename="${filename%.*}"
    # get the git diff for the file
    git diff dc915f7cda6b63921766b3a95c02f4607517876b 5b24014bb4e08b91c7026afe4e12d4ebda363215  $file > gitdiff_$filename.csv
    # add comma after every '-' and '+' in the diff file
    sed -i 's/^-/-,/g' gitdiff_$filename.csv
    sed -i 's/^+/+,/g' gitdiff_$filename.csv
done
