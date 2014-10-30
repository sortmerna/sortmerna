#!/bin/bash

rootdir="/home/evko1434/ARB_SILVA_databases"
reference_db=(bacteria_16S_119_clean-T.fasta archaea_16S_119_clean-T.fasta eukarya_18S_119_clean-T.fasta bacteria_23S_clean-T.fasta archaea_23S_clean-T.fasta eukarya_28S_clean-T.fasta RF00001-T.fasta.txt RF00002-T.fasta.txt)
id=(0.90 0.95 0.95 0.98 0.98 0.98 0.98 0.98)

#for i in 0;
for ((i=0;i<${#reference_db[@]};i++));
do
    echo "sumaclust -l -p 40 -t ${id[i]} -F $rootdir/step4-c/${reference_db[i]} $rootdir/${reference_db[i]}"
    echo "sumaclust -l -p 40 -t ${id[i]} -F $rootdir/step4-c/${reference_db[i]} $rootdir/${reference_db[i]}" | qsub -l nodes=1:ppn=40 -q mem4gbq -k oe -N ${db}_sumaclust; sleep 2
done