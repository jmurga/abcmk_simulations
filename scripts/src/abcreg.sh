seq ${1} ${2} | parallel -u -j ${} "/home/jmurga/ABCreg/srcreg -d alpha_${4}.tsv.gz -p ${4}_1.tsv.gz -P 3 -S 100 -L -t 0.01 -b test_{}"