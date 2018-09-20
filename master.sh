#$ -S /bin/bash
#$ -j y   						                                 #standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	                     #location of standard output
#$ -r y  						                                   #rerun job if necessary
#$ -N 'DG-ChIPseq'        	                           #give name to job
#$ -t 1                                                #no. jobs
#$ -V                                                  # exports all environmental variables to qsub

#Read samplesheet and scripts
##samplesheet needs to have one set of controls for all experimentals in tab separated file

workdir="/mnt/iscsi_speed/blelloch/deniz/testdata"
samplesheet="test_p300.txt"
samples="$workdir"/"$samplesheet"
scripts="/mnt/iscsi_speed/blelloch/deniz/scripts/DG-chipseq-pipeline"
GENOME="/mnt/iscsi_speed/blelloch/deniz/genomes/mm10/mm10"

#Remove ^Ms and add newline to EOF
{ tr '\r' '\n' < $samples; echo; } > "$samples"_temp
mv "$samples"_temp "$samples"

#Extract variables
nsamples=$(expr $(wc -l < "$samples") - 1)
treatments=($(awk '{ if(($6 == 1)) { print $2} }' "$samples" | uniq))
ctrs=($(awk '{ if(($6 == 0)) { print $1} }' "$samples"))

nctr=$(echo ${#ctrs[@]})

echo "Treatment conditions are: $treatments"

#Run mapping ($1=wordir, $2=GENOME, $3="--nodedup")
qsub -t 1-"$nsamples" "$scripts"/mapping-v2.sh "$workdir" "$GENOME" "$samples"

#Convert to .bw
qsub -t 1-"$nsamples" -hold_jid "DG-ChIPseq-mapping" "$scripts"/bam2bw-v2.sh "$workdir" "$samples"

#Cross-crosscorrelation
qsub -t 1-"$nsamples" -hold_jid "DG-bam2bw" "$scripts"/macs2_predictd.sh "$workdir" "$samples"

#Run MACS2 for each treatment vs controls

#Define merged control string
for (( i = 0; i < "$nctr"; i++ )); do
  optList1+=(-c $workdir/bam/${ctrs[$i]}.sorted.dedup.bam)
done

#Define each treatment and run MACS2 for each replicate of each treatment vs all ctrs (only if there are controls)
if [[ "$nctr" -ne 0 ]]; then
  for i in "${treatments[@]}"; do
    treat=($(awk -v i=$i '{ if($2 == i && $6 == 1) { print $1} }' "$samples")) #-v passes bash
    ntreat=$(echo ${#treat[@]})                                                #variable to awk
    for (( j = 0; j < "$ntreat"; j++ )); do
      optList2="-t ${workdir}/bam/${treat[$j]}.sorted.dedup.bam"
      qsub -hold_jid "DG-bam2bw" $scripts/macs2.sh $optList2 ${optList1[@]} -w $workdir -f $i
      echo "MACS2 call: qsub -hold_jid "DG-bam2bw","DG-crosscor" $scripts/macs2.sh $optList2 ${optList1[@]} -w $workdir"
    done
  done
fi
