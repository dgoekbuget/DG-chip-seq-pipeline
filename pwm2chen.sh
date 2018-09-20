#$ -S /bin/bash
#$ -l mem_free=20G
#$ -j y   						                                                         # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						                                                           # rerun job if necessary
#$ -N 'pwm2meme'	                                                               # give name to job
#$ -t 1                                                                      # adjust depending of no jobs needed
#$ -V                                                                          # exports all environmental variables to qsub

genome="/mnt/iscsi_speed/blelloch/deniz/genomes/mm10_UCSC.fa"
pwm="/mnt/iscsi_speed/blelloch/deniz/pwm/Mus_musculus_2018_05_04_6-16_pm/pwms_all_motifs"
RUNDATE=$(date +%Y-%m-%d)

#Compute background markov model for whole genome

if [[ ! -f "${genome%.fa}_background" ]]; then
	fasta-get-markov "$genome" > "${genome%.fa}_background"
fi

#Convert PWMs to Chen Matrix format
files=( "$pwm"/*.txt)

if [[ -f "$pwm"/pwms_all_motifs.chen ]]; then
	rm "$pwm"/pwms_all_motifs.chen
fi

touch "$pwm"/pwms_all_motifs.chen

for i in "${files[@]}"; do
	echo ">${i%.txt}" | sed "s?${pwm}/??" >> "$pwm"/pwms_all_motifs.chen
	cut -f2- "$i" | sed '1d' >> "$pwm"/pwms_all_motifs.chen
done

if [[ "$?" -eq 0 ]]; then
  cat "$pwm"/pwms_all_motifs.chen | chen2meme -bg "${genome%.fa}_background" > "$pwm"/cis_db_mm_pwm.meme
else
  echo "conversion to matrix failed!"
fi

if [[ "$?" -eq 0 ]]; then
  echo "pwm2meme finished successfully!"
else
  echo "chen2meme failed!"
fi
