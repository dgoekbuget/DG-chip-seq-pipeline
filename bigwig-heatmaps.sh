#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						                                                         # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						                                                           # rerun job if necessary
#S -t 1
#$ -N 'DG-macs-atac'	                                                               # give name to job
#$ -V                                                                          # exports all environmental variables to qsub

#Generate matrix file and heatmap for regions defined in .bed file and .bw files of interest

bws=("$@")
workdir="/mnt/iscsi_speed/blelloch/deniz/analysis"
ns=${#bws[@]}
peaks="$workdir"/out.sorted.bed
scripts="/mnt/iscsi_speed/blelloch/deniz/scripts"
heatwidth=1000

#Compute matrix for each file
qsub -t "$ns" "$scripts"/bw2matrix.sh $workdir $bws $peaks $heatwidth

#qsub -hold_jid "bw2matrix" "$scripts"/matrix-cbind.sh
