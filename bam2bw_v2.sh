#$ -S /bin/bash
#$ -M deniz.goekbuget@ucsf.edu
#$ -m besa
#$ -l mem_free=32G                             #star needs > 32G
#$ -j y                                        # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs                 # location of standard output
#$ -r y                                        # rerun job if necessary
#$ -N 'DG-bam2bw'                            # give name to job
#$ -V                                          # exports all environmental variables to qsub

workdir=$1
files=($(awk '{print $1}' $2 ))
infile="$workdir"/bam/"${files[$SGE_TASK_ID]}".sorted.dedup.bam
chromsize="/mnt/iscsi_speed/blelloch/deniz/genomes/mm10.chrom.sizes"

#Make outputdir
if [[ ! -d "$workdir"/bw ]]; then
  mkdir "$workdir"/bw
fi

# Calculate coverage with bedtools
reads=$(samtools view -c "$infile")
scale=$(bc -l <<< "1000000/${reads}") #expr only works with integers
rundate=$(date +%Y-%m-%d)

# generate bedgraph of coverage
bedtools genomecov -scale "${scale}" -bga -split -ibam "$infile" > "${infile%.sorted.dedup.bam}".dedup.unsorted.bedGraph

# sort bedGraph
if [[ "$?" -eq 0 ]];then
  echo "bedGraph created"
	bedtools sort -i "${infile%.sorted.dedup.bam}".dedup.unsorted.bedGraph > "${infile%.sorted.dedup.bam}".dedup.sorted.bedGraph
elif [[ "$?" -eq 1 ]];then
	echo "Creating bedGraph failed"
  exit 1
fi

# Convert bedgraph to bigwig
if [[ "$?" -eq 0 ]]; then
  bedGraphToBigWig "${infile%.sorted.dedup.bam}".dedup.sorted.bedGraph "$chromsize" "${infile%.sorted.dedup.bam}".sorted.bw
elif [[ "$?" -eq 1 ]]; then
  echo "sorting bedGraph failed"
  exit 1
fi

#Remove intermediate files and move bw files
if [[ "$?" -eq 0 ]]; then
  rm "${infile%.sorted.dedup.bam}".dedup.unsorted.bedGraph "${infile%.sorted.dedup.bam}".dedup.sorted.bedGraph
	mv "${infile%.sorted.dedup.bam}".sorted.bw "$workdir"/bw
fi
