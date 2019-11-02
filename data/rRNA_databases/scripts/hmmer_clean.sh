#!/bin/bash

# usage: run hmmer_clean.sh on rRNA databases (in FASTA format) to remove contaminant rRNA (ex. 23S rRNA found in 16S database)
# purpose: to remove sequences spanning multiple genes (+ITS)
# ex.
#        18S           ITS1     5.8S    ITS2                 28S
# |----------------||--------||-----||---------||----------------------------|
# {*******************************}
#          sequence spans 18s + 5.8S
#
# requirements: remove_partial_seqs.py (included with this bash script) and filter_fasta.py (part of QIIME's scripts)

rootdir=/home/evko1434/ARB_SILVA_databases
outdir=$rootdir/hmm_output
declare -A databases=( ["bacteria_16S_119"]="bac_lsu,bac_tsu" ["archaea_16S_119"]="arc_lsu,arc_tsu" ["eukarya_18S_119"]="euk_lsu,euk_tsu" ["bacteria_23S"]="bac_ssu,bac_tsu" ["archaea_23S"]="arc_ssu,arc_tsu" ["eukarya_28S"]="euk_ssu,euk_tsu" )
hmmer=/home/evko1434/hmmer-3.1b1-linux-intel-x86_64/binaries/nhmmer
qsub_param="-l nodes=1:ppn=20 -q mem4gbq -k oe"
mkdir $outdir

for db in "${!databases[@]}"
do
    hmms=(${databases["$db"]//,/ })

    echo "source /home/evko1434/.bashrc; $hmmer -o $outdir/${db}_vs_${hmms[0]}.txt --noali --tblout $outdir/${db}_vs_${hmms[0]}_table.txt --acc -E 1e-5 --cpu 20 ${rootdir}/HMM3/${hmms[0]}.hmm $rootdir/${db}.fasta; $hmmer -o $outdir/${db}_vs_${hmms[1]}.txt --noali --tblout $outdir/${db}_vs_${hmms[1]}_table.txt --acc -E 1e-5 --cpu 20 ${rootdir}/HMM3/${hmms[1]}.hmm $rootdir/${db}.fasta; python $rootdir/remove_partial_seqs.py $outdir/${db}_vs_${hmms[0]}_table.txt $outdir/${db}_partial_${hmms[0]}_to_remove.txt; python $rootdir/remove_partial_seqs.py $outdir/${db}_vs_${hmms[1]}_table.txt $outdir/${db}_partial_${hmms[1]}_to_remove.txt; cat $outdir/${db}_partial_${hmms[0]}_to_remove.txt $outdir/${db}_partial_${hmms[1]}_to_remove.txt > $outdir/${db}_partial_seqs_to_remove.txt; filter_fasta.py -f $rootdir/${db}.fasta -s $outdir/${db}_partial_seqs_to_remove.txt -n -o $rootdir/${db}_clean.fasta; python $rootdir/edit_U_to_T_rna.py $rootdir/${db}_clean.fasta $rootdir/${db}_clean-T.fasta" | qsub $qsub_param -N ${db}_clean; sleep 2
done






