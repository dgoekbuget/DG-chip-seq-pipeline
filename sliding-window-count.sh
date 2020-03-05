#$ -S /bin/bash
#$ -l mem_free=32G
#$ -j y   						                                                         # standard error and output merged
#$ -o /mnt/blelloch/Deniz/logs  	 # location of standard output
#$ -r y  						                                                           # rerun job if necessary
#S -t 1
#$ -N 'DG-sliding'	                                                               # give name to job
#$ -V                                                                          # exports all environmental variables to qsub

#Required arguments: each treatment (-t) and each control (-c) bed (tab-separated) file names and path to workdir (-w)

chromsize="/mnt/iscsi_speed/blelloch/deniz/genomes/mm10.chrom.sizes"

treats=()
controls=()
while getopts ":t:c:w:" opt; do
  case $opt in
  t)
    treats+=("${OPTARG}")
  ;;
  c)
    controls+=("${OPTARG}")
  ;;
  w)
    workdir="$OPTARG"
  ;;
  *)
    echo "Option not recognized. Enter each treatment (-t), control (-c) and work directory (-w)"
    exit 1
  ;;
  esac
done

echo Treatments are: "${treats[@]}"
echo Controls are: "${controls[@]}"
echo Workdir is: "$workdir"

#Define sliding parameters
windowsize=10000
slide=500

#Define samples sizes
ntreats=${#treats[@]}
nctr=${#controls[@]}
echo No. treatments: "$ntreats"
echo No. controls: "$nctr"

#Generate output folder
mkdir -p "$workdir"/sliding-window-"$windowsize"-"$slide"
mkdir -p "$workdir"/sliding-window-"$windowsize"-"$slide"/temp
temp="$workdir"/sliding-window-"$windowsize"-"$slide"/temp
outdir="$workdir"/sliding-window-"$windowsize"-"$slide"

#Define regions (ignoring last less than $windowsize bases)
if [ ! -f "$outdir"/chr_windows.bed ]; then
  if while IFS=$'\t' read -r chr length; do
    start=1
    end=$(( windowsize + 1 ))
      for (( i = 0; i < $(( length / slide )); i++ )); do
        echo -e "$chr\t$start\t$end"
        start=$(( start + 500 ))
        end=$(( end + 500 ))
      done
    done <"$chromsize" >> "$outdir"/chr_windows.bed
  then
    echo "Slided chromosome windows generated"
  else
    exit 1
    echo "Slided chromosome window generation failed."
  fi
else
  echo "Slided chromosome window generation skipped. File already exists."
fi

#Count reads mapping to each window
##Normalize treatments and controls

if ! for i in "${treats[@]}"; do
    out=$(basename "$i")
    reads=$(samtools view -c "$workdir"/"$i")
    bedtools intersect -c -a "$outdir"/chr_windows.bed -b "$workdir"/"$i" | awk -v OFS="\t" -v reads="$reads" '{print $1,$2,$3,($4+0.25*reads/1000000)/reads}' > "$temp"/"${out%.bam}".treat
  done
then
  exit 1
  echo "Normalizing treatments failed."
fi

if ! for i in "${controls[@]}"; do
    out=$(basename "$i")
    reads=$(samtools view -c "$workdir"/"$i")
    bedtools intersect -c -a "$outdir"/chr_windows.bed -b "$workdir"/"$i" | awk -v OFS="\t" -v reads="$reads" '{print $1,$2,$3,($4+0.25*reads/1000000)/reads}' > "$temp"/"${out%.bam}".control
  done
then
  exit 1
  echo "Normalizing controls failed."
fi

##Average treatments and controls

args=()
for i in "$temp"/*.treat; do
    args+=("$i")
done

if ! paste "${args[@]}" | awk -v N="$ntreats" '{ avg=0; for (i=4; i<=NF; i+=4) avg+=$i/N; print $1,$2,$3,avg }' > "${args[0]}".average
then
  exit 1
  echo "Averaging treatments failed."
fi

args=()
for i in "$temp"/*.control; do
    args+=("$i")
done

if ! paste "${args[@]}" | awk -v N="$nctr" '{ avg=0; for (i=4; i<=NF; i+=4) avg+=$i/N; print $1,$2,$3,avg }' > "${args[0]}".average
then
  exit 1
  echo "Averaging controls failed."
fi

##Calculate fold enrichment over controls

if ! paste "$temp"/"${treats[0]%.bam}".treat.average "$temp"/"${controls[0]%.bam}".control.average | awk '{ print $1,$2,$3,$4/$8}' > "$outdir"/"${treats[0]%.bam}"."$windowsize"bp.s"$slide".bed
then
  exit 1
  echo "Calculating enrichment over control failed."s
fi

rm "$temp"/*.treat
rm "$temp"/*.control
rm "$temp"/*.treat.average
rm "$temp"/*.control.average
if ! rmdir "$temp"; echo "Script finished successfully."
then
  exit 1
fi
