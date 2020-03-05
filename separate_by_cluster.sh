
nclust=$(awk '{ print $13 }' $1 | uniq | wc -l)
nclust=$(( nclust - 1 ))

for (( i = 1; i <= $nclust; i++ )); do
  awk -v i=$i -v OFS="\t" '{ if($13 == "cluster_"i) { print $0 }}' "$1" > ${1%.bed}.cluster_"$i".bed
done
