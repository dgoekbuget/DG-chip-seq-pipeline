#Converts .SortedRegions.bed output from deeptools into .bed file with nearest gene, and the respective logFC and p-val in naive to epilc

peakfile="/path/to/file.bed"
tss="/mnt/iscsi_speed/blelloch/deniz/genomes/TSS.mouse.GRCm38-2.sorted.bed"

#Catalogue of peaks and nearest TSS
bedtools closest -d -a <(sort -k 1,1 -k2,2n $peakfile) -b $tss > DG1420-idr005.noTSS.nearest.bed

#Add nearest gene to .SortedRegions.bed file
awk 'FNR==NR {a[$4]=$4; b[$4]=$6; c[$4]=$9 ; next} {if (!a[$4]) {a[$1]="NA"; b[$1]="NA"; c[$1]="NA"} printf "%s %s %s\n", $0, b[$4], c[$4]} ' SRR1202455_summits.noTSS.nearest.bed SRR1202455_summits.noTSS.nearest.SortedRegions.bed

#Add FC and p-Value naive vs EpiLC to .nearest.bed file
awk 'FNR==NR {a[$8]=$8; b[$8]=$2; c[$8]=$6; d[$8]=$9 ; next} {if (!a[$7]) {a[$7]="NA"; b[$7]="NA"; c[$7]="NA"; d[$7]="NA"} printf "%s %s %s %s\n", $0, b[$7], c[$7], d[$7]} ' ../../../../jake/rna-seq_for_deniz.tab DG1420-idr005.noTSS.nearest.bed

#Add FC and p-Value naive vs EpiLC to .SortedRegions.nearest.bed file
awk 'FNR==NR {a[$2]=$1; b[$2]=$2; c[$2]=$10; d[$2]=$11; e[$2]=$12; f[$2]=$13 ; next} {if (!a[$1] && !b[$2]) {c[$2]="NA"; d[$2]="NA"; e[$2]="NA"; f[$2]="NA"} printf "%s %s %s %s %s\n", $0, c[$2], d[$2], e[$2], f[$2]} ' DG1420-idr005.noTSS.nearest.rnaseq.bed DG1420-idr005.noTSS.SortedRegions.bed 
