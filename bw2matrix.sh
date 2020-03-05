#$ -S /bin/bash
#$ -M deniz.goekbuget@ucsf.edu
#$ -m besa
#$ -l mem_free=32G                             #star needs > 32G
#$ -j y                                        # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs                 # location of standard output
#$ -r y                                        # rerun job if necessary
#$ -N 'bw2matrix'                            # give name to job
#$ -V                                          # exports all environmental variables to qsub

workdir=$1
bws=$2
peaks=$3
heatwidth=$4
infile=${bws[$SGE_TASK_ID]}
outfile=${$(basename $infile)%.bw}.matrix
bl="/mnt/iscsi_speed/blelloch/deniz/genomes/blacklist/mm10.blacklist.bed"

computeMatrix reference-point -R $peaks -b $heatwidth -a $heatwidth -bl $bl -p "max" -out $outfile -S $infile
