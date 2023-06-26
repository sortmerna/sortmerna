Sample execution statistics v4.0.1
==================================

Statistics calculation on 2 reads files of ~50M reads each::
   
   [calculate:107] Starting statistics calculation on file: '/home/reads/NG-16660_IS1_lib274081_6037_6_f_qc.fastq.gz'  ...   [inflatez:169] xINFO: infl
   [calculate:225] Done statistics on file. Elapsed time: 185.09 sec. all_reads_count= 50207959
   [calculate:107] Starting statistics calculation on file: '/home/reads/NG-16660_IS1_lib274081_6037_6_r_qc.fastq.gz'  ...   [inflatez:169] xINFO: infl
   [calculate:225] Done statistics on file. Elapsed time: 186.63 sec. all_reads_count= 100415918

Run statistics:
- 8 databases
2 read files of ~100M reads total
8 threads (a high-end laptop as of 202001)

::

   num reads: 100,415,918

   ref                           hash                  size           sec        min     hr
   silva-bac-16s-id90.fasta      15734375058464002811  19,437,013     19589.84   326.48  5.44
   silva-bac-23s-id98.fasta      17299952793705614139  12,911,743      7313.21   121.88  2.03
   silva-arc-16s-id95.fasta      3436099190853847617    3,893,959      3047.73    50.78  0.87
   silva-arc-23s-id98.fasta      3400685301612210653      752,022       370.21     6.17  0.1
   silva-euk-18s-id95.fasta      2700646386527218729   13,259,584     11259.34   187.66  3.13
   silva-euk-28s-id98.fasta      1845323523482939374   14,945,070      4182.19    69.70  1.16
   rfam-5s-database-id98.fasta   13019673092862722585   8,525,326      3263.54    54.39  0.90
   rfam-5.8s-database-id98.fasta 2169995244134016533    2,280,449      3259.41    54.32  0.90
                                                       Total time (hr) for alignment:    14.5
