#$ -S /bin/bash
#$ -j y   						                                 #standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	                     #location of standard output
#$ -r y  						                                   #rerun job if necessary
#$ -N 'DG-ChIPseq-initiation'	                           #give name to job
#$ -t 1                                                #no. jobs
#$ -V                                                  # exports all environmental variables to qsub

#Read samplesheet and scripts
samples="/mnt/iscsi_speed/blelloch/Deniz/test.tab"
scripts="/mnt/iscsi_speed/blelloch/Deniz/scripts/DG-chipseq-pipeline"

#Extract parameters
workdir=""

#Run mapping ($1=wordir, $2=GENOME, $3="--nodedup")
qsub -t 1-"${nsamples}" "$scripts"/mapping.sh "$workdir" "$GENOME"

#Run MACS2

qsub -hold_jid "DG-ChIPseq-mapping" "$scripts"/macs2.sh -t -c

#Convert to .bw
