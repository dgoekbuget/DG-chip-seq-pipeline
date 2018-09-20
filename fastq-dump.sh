#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						             # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						               # rerun job if necessary
#$ -N 'DG-fastqdmp'       	       # give name to job
#$ -V                              # exports all environmental variables to qsub

workdir="$1"
files=($(awk '{print $1}' $2 ))
infile="${files[$SGE_TASK_ID]}"  # Pick one item from that array

#FASTQ dump
fastq-dump --gzip --outdir $workdir $infile
