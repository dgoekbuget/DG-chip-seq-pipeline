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
#optionally specify number of clusters for k-means
if [ -z "$2" ]; then
	km=10
else
	km="$2"
fi
#Output
outfile=${infile%.bed}.matrix
#Blacklist .bed file to be excluded
blacklist="/mnt/iscsi_speed/blelloch/deniz/genomes/blacklist/mm10.blacklist.bed"

#Bigwig files to be included in heatmap
foxd3c1="/mnt/iscsi_speed/blelloch/deniz/dgrb02_foxd3/naive/bw/"
foxd3c2="/mnt/iscsi_speed/blelloch/deniz/dgrb02_foxd3/naive/bw/"
foxd3c3="/mnt/iscsi_speed/blelloch/deniz/dgrb02_foxd3/naive/bw/"
foxd3t1="/mnt/iscsi_speed/blelloch/deniz/dgrb02_foxd3/naive/bw/"
foxd3t2="/mnt/iscsi_speed/blelloch/deniz/dgrb02_foxd3/naive/bw/"
foxd3t3="/mnt/iscsi_speed/blelloch/deniz/dgrb02_foxd3/naive/bw/"


if computeMatrix reference-point -R "$infile" -b 3000 -a 3000 -bs 20 \
  -bl "$blacklist" -p "max" -out "$outfile" \
    -S $foxd3 $sall4 $oct4 $p300 $h3k27ac $h3k4me1 $h3k4me3 $mll34 $atac $smc1 $pol2 "$GRO_F" "$GRO_R"; then
   plotHeatmap -m "$outfile" \
    -o ${outfile%.matrix}.eps \
      --kmeans $km --colorList '#ffffff,blue,#000000' \
        --samplesLabel "FOXD3" "SALL4" "OCT4" "P300" "H3K27ac" "H3K4me1" "H3K4me3" "MLL3/4" "ATAC" "SMC1" "POL2" "GRO-F" "GRO-R" \
					--outFileSortedRegions ${outfile%.matrix}.SortedRegions.bed
fi
