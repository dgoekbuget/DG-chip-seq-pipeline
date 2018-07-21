#$ -S /bin/bash
#$ -l mem_free=20G
#$ -j y   						                                                         # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						                                                           # rerun job if necessary
#$ -N 'DG-ChIPseq-MACS2'	                                                               # give name to job
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

#Create required directories

if [[ ! -d "$workdir"/macs2 ]]; then
  mkdir "$workdir"/macs2
fi

if [[ ! -d "$workdir"/temp ]]; then
  mkdir "$workdir"/temp
fi

#MACS2 peak finding
source activate python2

macs2 callpeak -t "${treatment[@]}" -c "${control[@]}" -n "${treatment[0]}" --outdir "$workdir"/macs2 -f BAM -g mm -B -p 0.01 --tempdir "$workdir"/temp --verbose 0
