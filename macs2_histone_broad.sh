#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						                                                         # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						                                                           # rerun job if necessary
#$ -N 'DG-MACS2'	                                                               # give name to job
#$ -V                                                                          # exports all environmental variables to qsub

#Read treatment and control files by quoting paths to treatment (-t) and control (-c) files.

treatment=()
control=()
while getopts ":t:c:w:" opt; do
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
  esac
done

#Outputname
outname="${treatment[0]}"

#Create required directories

if [[ ! -d "$workdir"/macs2 ]]; then
  mkdir "$workdir"/macs2
fi

if [[ ! -d "$workdir"/temp ]]; then
  mkdir "$workdir"/temp
fi

if [[ ! -d "$workdir"/macs2/"${outname%.sorted.dedup.bam}" ]]; then
  mkdir "$workdir"/macs2/"${outname%.sorted.dedup.bam}"
fi

#MACS2 peak finding
source activate python2

macs2 callpeak -t "${treatment[@]}" -c "${control[@]}" -n "${outname%.sorted.dedup.bam}" --outdir "$workdir"/macs2/"${outname%.sorted.dedup.bam}" -f BAM -g mm --broad -p 1e-2 --nomodel --shift 0 --extsize 200 --keep-dup all --tempdir "$workdir"/temp --verbose 0

source deactivate

#Remove temporary files
if [[ "$?" -eq 0 ]];then
  rm "$workdir"/temp/*
  rmdir "$workdir"/temp
fi

#Generate crosscorrelation PDFs
Rscript "$workdir"/macs2/${outname%.sorted.dedup.bam}"_model.R
