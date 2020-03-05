#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						                                                         # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						                                                           # rerun job if necessary
#S -t 1
#$ -N 'DG-MACS2'	                                                               # give name to job
#$ -V                                                                          # exports all environmental variables to qsub

#Read treatment and control files by quoting paths to treatment (-t) and control (-c) files.

chromsize="/mnt/iscsi_speed/blelloch/deniz/genomes/mm10.chrom.sizes"

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

#MACS2 peak finding
source activate python2

macs2 callpeak -t "${treatment[@]}" -c "${control[@]}" -n "${outname%.sorted.dedup.bam}" --outdir "$workdir"/macs2/"$factor" --nomodel --extsize 200 -f BAM -g mm -B --SPMR -p 0.01 --tempdir "$workdir"/temp

if [[ "$?" -eq 0 ]];then
  echo "Peak Calling Finished. Starting FE calculation..."
  macs2 bdgcmp -t "$workdir"/macs2/"$factor"/"${outname%.sorted.dedup.bam}"_treat_pileup.bdg -c "$workdir"/macs2/"$factor"/"${outname%.sorted.dedup.bam}"_control_lambda.bdg -o "${outname%.sorted.dedup.bam}"_FE.bdg -m FE --outdir "$workdir"/macs2/"$factor"
fi

conda deactivate
source activate chip

if [[ "$?" -eq 0 ]];then
  echo "FE calculated. Sorting bedgraph..."
	bedtools sort -i "$workdir"/macs2/"$factor"/"${outname%.sorted.dedup.bam}"_FE.bdg > "$workdir"/macs2/"$factor"/${outname%.sorted.dedup.bam}_FE.sorted.bdg
elif [[ "$?" -eq 1 ]];then
	echo "Creating bedGraph failed"
  exit 1
fi

if [[ "$?" -eq 0 ]];then
  echo "Bedgraph sorted. Converting to bigwig..."
  bedGraphToBigWig "$workdir"/macs2/"$factor"/"${outname%.sorted.dedup.bam}"_FE.sorted.bdg "$chromsize" "$workdir"/macs2/"$factor"/"${outname%.sorted.dedup.bam}"_FE.bw
fi

if [[ "$?" -eq 0 ]];then
  echo "Bigwig generated. Removing intermediate files..."
  rm "$workdir"/macs2/"$factor"/"${outname%.sorted.dedup.bam}"_FE.bdg "$workdir"/macs2/"$factor"/"${outname%.sorted.dedup.bam}"_FE.sorted.bdg
fi

conda deactivate
