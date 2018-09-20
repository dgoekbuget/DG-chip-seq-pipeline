#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						                                                         # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						                                                           # rerun job if necessary
#$ -N 'centrimo_meme'	                                                               # give name to job
#$ -t 1                                                                      # adjust depending of no jobs needed
#$ -V                                                                          # exports all environmental variables to qsub

#Runs centrimo meme motif algorithm on mouse PWMs converted in .meme format using pwm2meme.sh
#Takes .fa file with equal lengths as $1 and working directory as $2

genome="/mnt/iscsi_speed/blelloch/deniz/genomes/mm10_UCSC.fa"
pwm="/mnt/iscsi_speed/blelloch/deniz/pwm/Mus_musculus_2018_05_04_6-16_pm/pwms_all_motifs"
workdir="$2"
RUNDATE=$(date +%Y-%m-%d)
infile="$1"
outdir=$(basename $infile)

#Run CentriMo

centrimo -o "$workdir"/centrimo_"${outdir%.fa}"_"${RUNDATE}" "$infile" "$pwm"/cis_db_mm_pwm.meme

if [[ "$?" -eq 0 ]]; then
	#Name motifs with file names
	for i in "${files[@]}"; do
		sed -i 's/MOTIF [0-9]*/${i%.txt}/' "$pwm"/cis_db_mm_pwm.meme
	done
else
	echo "CentriMo failed"
fi

if [[ "$?" -eq 0 ]]; then
  echo "centrimo_meme finished successfully!"
else
  echo "Text editing step failed."
fi
