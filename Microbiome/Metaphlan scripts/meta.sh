#!/bin/bash 


folder='CD_paper'

home='/Volumes/Han_SSD'
in_dir=${home}/${folder}
out_dir=${home}/output/${folder}

source activate mpa 


function run_metaphlan () 
{
	metaphlan ${in_dir}/$1.fastq.gz --input_type fastq --no_map -o ${out_dir}/$1.txt
}

# run_metaphlan 'test'

for path in ${in_dir}/C*_A_R2.fastq.gz; do

	file=$(basename $path)
	filename=${file%%.*}

	run_metaphlan $filename
	
	printf '\t--------------------   %s DONE   -------------------\n\n\n' "$filename"

done 

