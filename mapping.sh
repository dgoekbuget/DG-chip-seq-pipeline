#$ -S /bin/bash
#$ -l mem_free=20G
#$ -j y   						             # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						               # rerun job if necessary
#$ -N 'DG-ChIPseq-mapping'	       # give name to job
#$ -V                              # exports all environmental variables to qsub

i=$(expr $SGE_TASK_ID - 1)
workdir="$1"
files=( "$workdir"/*.fastq.gz) # Expand the glob, store it in an array (Charles Duffy)
infile="${files[$i]}"  # Pick one item from that array
GENOME="$2"

#Bowtie
bowtie2 -p 4 -q -D 15 -R 10 -L 22 -i S,1,1.15 -x $GENOME "$infile" -S "${infile%.fastq.gz}.sam"

#Convert sam to bam
if [[ "$?" -eq 0 ]];then
  echo "Bowtie done"
	samtools view -bS "${infile%.fastq.gz}.sam" > "${infile%.fastq.gz}.sam.bam"
elif [[ "$?" -eq 1 ]];then
	echo "Bowtie2 failed"
  exit 1
fi

#Remove mapped reads of low quality
if [[ "$?" -eq 0 ]];then
	echo"sam2bam done"
  samtools view -h -q 31 -F 4 -b "${infile%.fastq.gz}.sam.bam" > "${infile%.fastq.gz}.mapq31.bam"
elif [[ "$?" -eq 1 ]];then
	echo "Removing mapped reads of low quality failed"
  exit 1
fi

#Sort bam file
if [[ "$?" -eq 0 ]];then
	echo "lowQ reads removed"
  samtools sort "${infile%.fastq.gz}.mapq31.bam" > "${infile%.fastq.gz}.sorted.bam"
elif [[ "$?" -eq 1 ]];then
	echo "Sorting of bam file failed"
  exit 1
fi

#Deduplicate
if [[ "$3" -eq "--nodedup" ]];then
  echo "no deduplication performed"
elif [[ "$?" -eq 0 ]];then
  echo "bam file sorted"
	samtools rmdup -s "${infile%.fastq.gz}.mapq31.sorted.bam" "${infile%.fastq.gz}.sorted.dedup.bam"
elif [[ "$?" -eq 1 ]];then
  echo "Deduplication failed"
  exit 1
fi

#QC
##Fastqc
if [[ "$?" -eq 0 ]]; then
  zcat "$infile" | fastqc stdin -o "$workdir"
fi

#Run multiqc
if [[ "$?" -eq 0 ]];then
	multiqc "$workdir"
fi

#Remove intermediate files and tidy up
if [[ "$?" -eq 0 ]];then
  rm "${infile%.fastq.gz}.sam" "${infile%.fastq.gz}.mapq31.bam"
  if [[ ! -d "$workdir"/bamfiles ]]; then
    mkdir "$workdir"/bamfiles
    mv "$workdir"/*.bam "$workdir"/bamfiles
  fi
elif [[ "$?" -eq 1 ]];then
	echo "bowtie2 failed"
  exit 1
fi
