#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						             # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						               # rerun job if necessary
#$ -N 'DG-mapping'	       # give name to job
#$ -V                              # exports all environmental variables to qsub

workdir="$1"
GENOME="$2"
files=($(awk '{print $1}' $3 ))
infile="$workdir"/"${files[$SGE_TASK_ID]}".fastq.gz  # Pick one item from that array

#Bowtie
bowtie2 -p 4 -q -D 15 -R 10 -L 22 -i S,1,1.15 -x $GENOME "$infile" -S "${infile%.fastq.gz}".sam

#Convert sam to bam
if [[ "$?" -eq 0 ]];then
  echo "Bowtie done"
	samtools view -bS "${infile%.fastq.gz}".sam > "${infile%.fastq.gz}".sam.bam
elif [[ "$?" -eq 1 ]];then
	echo "Bowtie2 job $JOB_ID failed" >> "$workdir"/log
  exit 1
fi

#Remove mapped reads of low quality
if [[ "$?" -eq 0 ]];then
	echo"sam2bam done"
  samtools view -h -q 31 -F 4 -b "${infile%.fastq.gz}".sam.bam > "${infile%.fastq.gz}".mapq31.bam
elif [[ "$?" -eq 1 ]];then
	echo "sam2bam job $JOB_ID failed"
  exit 1
fi

#Sort bam file
if [[ "$?" -eq 0 ]];then
	echo "lowQ reads removed"
  samtools sort "${infile%.fastq.gz}".mapq31.bam > "${infile%.fastq.gz}".sorted.bam
elif [[ "$?" -eq 1 ]];then
	echo "Removal of lowQ reads failed" >> "$workdir"/log
  exit 1
fi

#Deduplicate
if [[ "$?" -eq 0 ]];then
  echo "bam file sorted"
	samtools rmdup -s "${infile%.fastq.gz}".sorted.bam "${infile%.fastq.gz}".sorted.dedup.bam
elif [[ "$?" -eq 1 ]];then
  echo "Sorting job $JOB_ID failed" >> "$workdir"/log
  exit 1
fi

#Index Bam files

if [[ "$?" -eq 0 ]]; then
  echo "deduplication successful"
  samtools index "${infile%.fastq.gz}.sorted.dedup.bam"
elif [[ "$?" -eq 1 ]];then
  echo "Deduplication job $JOB_ID failed" >> "$workdir"/log
  exit 1
fi

#QC
##Fastqc
if [[ "$?" -eq 0 ]]; then
  if [[ ! -d "$workdir"/qc ]]; then
    mkdir "$workdir"/qc
  fi
  fastqc "$infile" -o "$workdir"/qc
fi

#Run multiqc
if [[ "$?" -eq 0 ]];then
	multiqc "$workdir"/qc -o "$workdir"/qc
fi

#Remove intermediate files and tidy up
if [[ "$?" -eq 0 ]];then
  rm "${infile%.fastq.gz}".sam "${infile%.fastq.gz}".mapq31.bam "${infile%.fastq.gz}".sorted.bam "${infile%.fastq.gz}".sam.bam
  #Create directories
  if [[ ! -d "$workdir"/fastq ]]; then
    mkdir "$workdir"/fastq
  fi
  if [[ ! -d "$workdir"/bam ]]; then
    mkdir "$workdir"/bam
  fi
  #Move output files in directories
  mv "${infile%.fastq.gz}".sorted.dedup.bam "$workdir"/bam
  mv "${infile%.fastq.gz}".sorted.dedup.bam.bai "$workdir"/bam
  mv "$infile" "$workdir"/fastq
elif [[ "$?" -eq 1 ]];then
	echo "multiqc failed"
  exit 1
fi
