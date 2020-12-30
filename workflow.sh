#!/bin/bash

# Create folder final and copy data into folder /data
cd final
mkdir data
cp /data-shared/vcf_examples/luscinia_vars.vcf.gz data/

# Set up the INPUT file
# INPUT= /data/luscinia_vars.vcf

# Look at the data
# zcat <$INPUT | less -S


# Unzip the file 
gunzip data/luscinia_vars.vcf.gz

# Set up the INPUT file - unzipped
INPUT=data/luscinia_vars.vcf

# Get rid of comment lines starting with ## and visually inspect
# <$INPUT grep -v '##' | less -S

# Select the columns 1 to 6 from orig vcf file
<$INPUT grep -v '^##' | cut -f1-6  >data/01-data.tsv

# Create new columns containing DP, INDEL and AF1 from the column INFO
<$INPUT egrep -o 'DP=[^;]*' | sed 's/DP=//' > data/01-data-DP.tsv
<$INPUT awk '{if($0 ~ /INDEL/) print "INDEL"; else print "SNP"}' > data/01-data-is-INDEL.tsv
<$INPUT egrep -o 'AF1=[^;]*' | sed 's/AF1=//' > data/01-data-AF.tsv

# Create column names for the new columns
sed  -i '1i DP' data/01-data-DP.tsv 
sed  -i '1i ISINDEL' data/01-data-is-INDEL.tsv 
sed  -i '1i AF' data/01-data-AF.tsv 

# Paste all columns to the final data file
paste data/01-data.tsv data/01-data-DP.tsv data/01-data-is-INDEL.tsv data/01-data-AF.tsv >data/02-data-with-DP-is-INDEL.tsv
