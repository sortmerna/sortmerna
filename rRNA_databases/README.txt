Date: 29/11/2014
Contact: Jenya Kopylov (jenya.kopylov@gmail.com)

-------
Content
-------

Use HMMER & SumaClust to generate representative rRNA database for SortMeRNA metatranscriptomic read filtering.
Scripts used to generate these databases are in the folder ./rRNA_databases/scripts/ . They may need to be
modified for your specific system in order to work (if you wish to re-build the databases).

-------
Methods
-------

1. Use the ARB package [1] to extract FASTA files for 16S bacteria, 16S archaea and 18S eukarya using SSURef_NR99_119_SILVA_14_07_14_opt.arb
and 23S bacteria, 23S archaea and 28S eukarya using LSURef_119_SILVA_15_07_14_opt.arb.

2. RFAM 5S and 5.8S databases [2] will be updated at the next release of SortMeRNA (2.1) 

3. Remove partial sequences from all databases using HMMER 3.1b1 [3] (meta_RNA [4] HMM profiles used)

>> bash ./scripts/hmmer_clean.sh

4. Use SumaClust version 1.0.00 [5] to cluster sequences at various % id using the FASTA files generated in (3) 

>> bash ./scripts/generate_databases.sh

5. Extract only the cluster centers

ex.

awk '/^>/ {printf("\n%s\t",$0);next;} {printf("%s",$0);} END {printf("\n");}' sumaclust_output_bacteria_16S_119-T.fasta \
           | egrep -v '^$' \
           | grep "cluster_center=True;" \
           | tr "\t" "\n" \
           > silva-bac-16s-id90.fasta 

---------
Citations
---------

[1] Wolfgang Ludwig et al., "ARB: a software environment for sequence data". Nucleic Acids Research (2004) 32:1363-1371
[2] Sarah W. Burge et al., "Rfam 11.0: 10 years of RNA families". Nucleic Acids Research (2012) doi: 10.1093/nar/gks1005
[3] R.D. Finn, J. Clements and S.R. Eddy, "HMMER web server: interactive sequence similarity searching". Nucleic Acids Research (2011) Web Server Issue 39:W29-W37
[4] Ying Huang, Paul Gilna and Weizhong Li, "Identification of ribosomal RNA genes in metagenomic fragments". Bioinformatics (2009) 25:1338-1340
[5] Céline Mercier et al., SUMATRA and SUMACLUST: fast and exact comparison and clustering of full-length barcode sequences. (submitted, 2014)
