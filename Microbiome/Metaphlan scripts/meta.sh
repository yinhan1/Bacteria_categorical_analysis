#!/bin/bash 


folder='CD_paper'     
replicate='*_B_R2'


# -----------------    set up directories    ------------------ #

# home directory 
home='/Volumes/Han_SSD'

# working directories for input and output data
in_dir=${home}/${folder}
out_dir=${home}/output/${folder}



# -----------------    metaphlan functions    -----------------  # 

# activate conda environment
source activate mpa 

# function to run metaphlan script
function run_metaphlan () 
{ 
	metaphlan ${in_dir}/$1.fastq.gz --input_type fastq --no_map -o ${out_dir}/$1.txt 
}



# ---------------   for-loop over the folder   ---------------  # 

# print start time 
start_time=$(date)
printf '\n======================   SCRIPT STARTS @ %s  ======================\n\n\n' "$start_time"


for path in ${in_dir}/${replicate}.fastq.gz; do
    
    # extract file names without extensions (fastq.gz)
	file=$(basename $path)
	filename=${file%%.*}

	# apply metaphlan script
	run_metaphlan $filename

	# create timestamp when job is done 
	timestamp=$(date)

	# print timestamp when job is done 
	printf '\n--------------------   %s COMPLETED @ %s  -------------------\n\n\n' "$filename" "$timestamp"

done 

# Lengths of sequence and quality values differs for UNC10-SN254:591:HAA11ADXX:1:2203:7086:40592 2:N:0:ACGATA (100 and 2305).

