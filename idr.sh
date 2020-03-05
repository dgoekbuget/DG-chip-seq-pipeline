#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						                                                         # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						                                                           # rerun job if necessary
#S -t 1
#$ -N 'DG-MACS2'	                                                               # give name to job
#$ -V                                                                          # exports all environmental variables to qsub

#Read treatment and control files by quoting paths to treatment (-t) and control (-c) files.

treatment=()
control=()
while getopts ":t:c:w:f:" opt; do
  case $opt in
  t)
    treatment+=("${OPTARG}")
  ;;
  c)
    control+=("${OPTARG}")
  ;;
  w)
    workdir="$OPTARG"
  ;;
  f)
    factor="$OPTARG"
  ;;
  esac
done

echo "Treatments are: "${treatment[@]}""
echo "Controls are: "${control[@]}""
echo "Workdir is: "$workdir""

#Outputname
outname=$(basename "${treatment[0]}")

#Create required directories

mkdir -p "$workdir"/macs2

mkdir -p "$workdir"/temp

mkdir -p "$workdir"/macs2/"$factor"

mkdir -p "$workdir"/macs2/"$factor"/idr

#idr filter

idr --samples samples1.narrowPeak sample2.narrowpeak --output-file-type "bed" --plot --output-file DG14vsDG22.txt
