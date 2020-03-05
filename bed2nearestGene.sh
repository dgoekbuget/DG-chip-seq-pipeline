awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.vM18.basic.annotation.gtf | gtf2bed - > gencode.vM18.basic.annotation.bed

#find closest TSS (with overlaps) and report distance. Keeps gene name column, strand, distance.
bedtools closest -d -a SRR1202455_summits.noTSS.bed -b /mnt/iscsi_speed/blelloch/deniz/genomes/TSS.mouse.GRCm38_2.bed | cut -f 1-5,9-12 >
