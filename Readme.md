# Parse Pileups

## Overview
This repository consists of 3 scripts that together takes fastq files aligns them to the mitochondrial genome (provided) 
and then parses the resulting pileup file. 

## Installation 

All you need to to do is clone this git repository to a local folder using and run the first time setup script

```bash
git clone https://github.com/celalp/parse_pileup
cd parse_pileup
bash prepare_submit.sh
```

## Usage

Before you start submitting samples you need to prepare a samples file. This is a three column csv without headers where
each row represents a single sample. Columns are as follows strictly in this order. 

1) Sample id
2) full path of the first fastq read
3) full path of the second fastq read, if this is a single end library then it has to be `None`

There is an example of this this file called `example.csv`

After the sample file is created you can use the `submit.sh` script to submit. This script also has arguments they need 
to be used in the following strict order. 

1) Samples file
2) directory to place the outputs if it does not exist will be created
3) script output this file will contain the job ids of each submitted job

You can then submit jobs to the cluster using:

```bash
bash submit.sh sample_file results_directory output_file
```

## Output

If the pipeline finishes successfully under the `results_directory` there will be several folders. Sorted and indexed 
bam files will be under the `alignment` folder. Rotated and original pileup files will be under `pileup` folder and 
finally the procesed pileups and a separate fasta will be under `processed` folder. If you wish to include the fasta 
sequence within the report you can delete the -f option from the `parse_pileup.py`. 

If there is an error each of the steps (alignment, pileup, processing) have their own log files. Please check these files 
for any errors. These log files have the structure `<samplename>.<step>.log`. For the pileup files any warnings regarding
whether therea are any regions that have 0 coverage will also be under `<samplename>.process.log`. 

Additionally there will be general .out and .err files for each sample within the output directory. These 
also might give clues about what went wrong.    

If you have any questions email me at alper.celik@sickkids.ca