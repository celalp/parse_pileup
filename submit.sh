
#!/bin/bash

################################################
################################################
# This is a wrapper script to submit single.sh 
# this script takes a comma separated file that
# includes sample id, and absolute paths of 
# fastq files and submits a job with it's own
# error and log files. The results will be 
# specified in the $outdir variable.
# see Readme.md for more info
###############################################

files_file=$1
outdir=$2
submission=$3

mkdir -p $outdir

while read line
 do
 sample=$(echo $line | cut -d "," -f 1)
 fastq1=$(echo $line | cut -d "," -f 2)
 fastq2=$(echo $line | cut -d "," -f 3)
 echo "$sample is being submitted results will be in $outdir"
 qsub -v sample=$sample,outdir=$outdir,fastq1=$fastq1,fastq2=$fastq2 -e $outdir/$sample.err -o $outdir/$sample.out -N $sample single.sh
done < $files_file > $submission
