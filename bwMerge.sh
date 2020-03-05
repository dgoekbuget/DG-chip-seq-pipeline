#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y                                                                        # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs     # location of standard output
#$ -r y                                                                        # rerun job if necessary
#$ -N 'DG-bwMerge'                                                                 # give name to job
#$ -V
#$ -t 1

workdir=$1
i=$(expr $SGE_TASK_ID - 1)
files=("$workdir"/*.bw) # Expand the glob, store it in an array (Charles Duffy)
chromsize="/mnt/iscsi_speed/blelloch/deniz/genomes/mm10.chrom.sizes"
outfile=$(basename ${files[0]})

if bigWigMerge "${files[@]}" "$workdir"/out.bedGraph; then
 if sort -k1,1 -k2,2n "$workdir"/out.bedGraph > "$workdir"/out.sorted.bedGraph; then
   if bedGraphToBigWig "$workdir"/out.sorted.bedGraph $chromsize "$workdir"/${outfile%.bw}.merged.bw; then
   rm "$workdir"/out.bedGraph "$workdir"/out.sorted.bedGraph
   fi
 fi
fi
