#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						                                                         # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						                                                           # rerun job if necessary
#S -t 1
#$ -N 'DG-MACS2-1rep'	                                                               # give name to job
#$ -V                                                                          # exports all environmental variables to qsub

#Read treatment and control files by quoting paths to treatment (-t) and control (-c) files.

chromsize="/mnt/iscsi_speed/blelloch/deniz/genomes/mm10.chrom.sizes"
workdir="$1"
samplesheet="$2"
samples="$workdir"/"$samplesheet"
#Remove ^Ms and add newline to EOF
{ cat $samples; echo; } > "$samples"_temp
mv "$samples"_temp "$samples"

files=($(awk '{print $1}' $samples ))
infile="${files[$SGE_TASK_ID]}".sorted.dedup.bam  # Pick one item from that array

#Create required directories

mkdir -p "$workdir"/macs2

mkdir -p "$workdir"/temp

mkdir -p "$workdir"/macs2/raw/

#MACS2 peak finding
source activate python2

macs2 callpeak -t "$workdir"/"$infile" -n "${infile%.sorted.dedup.bam}" --outdir "$workdir"/macs2/raw --nomodel --extsize 200 -f BAM -g mm -B --SPMR -p 0.01 --tempdir "$workdir"/temp

conda deactivate

source activate chip

if [[ "$?" -eq 0 ]];then
  echo "Peaks called. Sorting bedgraph..."
	bedtools sort -i "$workdir"/macs2/raw/"${infile%.sorted.dedup.bam}"_treat_pileup.bdg > "$workdir"/macs2/raw/"${infile%.sorted.dedup.bam}".pileup.sorted.bdg
elif [[ "$?" -eq 1 ]];then
	echo "Creating bedGraph failed"
  exit 1
fi

if [[ "$?" -eq 0 ]];then
  echo "Bedgraph sorted. Converting to bigwig..."
  bedGraphToBigWig "$workdir"/macs2/raw/"${infile%.sorted.dedup.bam}".pileup.sorted.bdg "$chromsize" "$workdir"/macs2/raw/"${infile%.sorted.dedup.bam}".pileup.bw
fi

if [[ "$?" -eq 0 ]];then
  echo "Bigwig generated. Removing intermediate files..."
  rm "$workdir"/macs2/raw/"${infile%.sorted.dedup.bam}".pileup.sorted.bdg
fi

conda deactivate
