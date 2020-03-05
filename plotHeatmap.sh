#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y                                                                        # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs     # location of standard output
#$ -r y                                                                        # rerun job if necessary
#$ -N 'DG-matrix'                                                                 # give name to job
#$ -V
#$ -t 1                                                                       # exports all environmental variables to qsub

#Input .bed file of regions
infile="$1"

#Plot heatmap

plotHeatmap -m "$infile" -o ${infile%.matrix}.pdf --samplesLabel "FOXD3 n" "OCT4 n" "KLF4 n" "NANOG n" "ATAC n" --sortUsing sum --sortUsingSamples 1 --outFileSortedRegions ${infile%.matrix}.SortedRegions.bed --outFileNameMatrix ${infile%.matrix}.tab --colorList '#ffffff,blue,#000000' --missingDataColor 1
