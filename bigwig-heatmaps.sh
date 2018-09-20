#Generate matrix file and heatmap for regions defined in .bed file and .bw files of interest

workdir="$1"
files=($(awk '{print $1}' $2 ))
peaks="$workdir"/"$3"

#Compute matrix
computeMatrix reference-point -R $peaks -b 1000 -a 1000 -bl ../../genomes/blacklist/mm10.blacklist.bed -p "max" -out output -S ../../dgrb01_foxd3/macs2/dgrb14_1e3_treat_pileup_FE.bw
