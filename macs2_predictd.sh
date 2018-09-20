#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						             # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						               # rerun job if necessary
#$ -N 'DG-crosscor'	       # give name to job
#$ -V                              # exports all environmental variables to qsub

workdir="$1"
files=($(awk '{print $1}' $2 ))
infile="$workdir"/bam/"${files[$SGE_TASK_ID]}".sorted.dedup.bam  # Pick one item from that array

#Create required directories

mkdir -p "$workdir"/macs2/qc

#MACS2 peak finding
source activate python2

macs2 predictd -i "$infile" -f BAM -g mm -m 5 50 --outdir "$workdir"/macs2/qc --rfile $(basename "${infile%.sorted.dedup.bam}"_model.r)

source deactivate
