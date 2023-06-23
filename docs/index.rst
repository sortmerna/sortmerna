User guide
==========

**SortMeRNA** is a local sequence alignment tool for filtering, mapping and OTU clustering. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of NGS reads. The main application of SortMeRNA is filtering rRNA from metatranscriptomic data. Additional applications include OTU-picking and taxonomy assignation available through `QIIME v1.9+
<http://qiime.org>`_. SortMeRNA takes as input one or two (paired) reads file(s) (fasta/fasta.gz/fastq/fastq.gz), and one or multiple rRNA database file(s), and sorts apart aligned and rejected reads into two files. SortMeRNA works with Illumina, 454, Ion Torrent and PacBio data, and can produce SAM and BLAST-like alignments.

.. note::
   
   This project is under active development.

Basic usage
-----------

The only required options are :code:`--ref` and :code:`--reads`. Options (any) can be specified usig a single dash e.g. :code:`-ref` and :code:`-reads`. Both plain :code:`fasta/fastq` and archived :code:`fasta.gz/fastq.gz` files are accepted. File extensions :code:`.fastq`, :code:`.fastq.gz`, :code:`.fq`, :code:`.fq.gz`, :code:`.fasta`, ... are optional. The format and compression are automatically recognized. Relative paths are accepted.

Example 1
#########

Single reference and single reads file::

   sortmerna --ref REF_PATH --reads READS_PATH

Example 2
#########

For multiple references use multiple :code:`--ref`::

   sortmerna --ref REF_PATH_1 --ref REF_PATH_2 --ref REF_PATH_3 --reads READS_PATH

Example 3
#########

For Paired reads use :code:`--reads` twice::

   sortmerna --ref REF_PATH_1 --ref REF_PATH_2 --ref REF_PATH_3 --reads READS_PATH_1 --reads READS_PATH_2

B
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

Sample execution trace v4.0.1
=============================

::
   
   biocodz@ubuntu16:~/sortmerna$ python scripts/run.py --name t17
   
   Current dir: /home/biocodz/sortmerna/scripts
   Using Environment configuration file: /home/biocodz/sortmerna/scripts/env.yaml
   Using Build configuration template: /home/biocodz/sortmerna/scripts/test.jinja.yaml
   Removing KVDB dir: /home/biocodz/sortmerna/run/kvdb
   Removing OUT_DIR: /home/biocodz/sortmerna/run/out
   Running t17: test_indexing
   
   [run] Running: /home/biocodz/sortmerna/dist/bin/sortmerna -ref data/rRNA_databases/silva-euk-28s-id98.fasta -ref data/rRNA_databases/silva-euk-18s-id95.fasta -ref data/rRNA_databases/silva-bac-23s-id98.fasta -ref data/rRNA_databases/silva-bac-16s-id90.fasta -ref data/rRNA_databases/silva-arc-23s-id98.fasta -ref data/rRNA_databases/silva-arc-16s-id95.fasta -ref data/rRNA_databases/rfam-5s-database-id98.fasta -ref data/rRNA_databases/rfam-5.8s-database-id98.fasta -reads data/set4_mate_pairs_metatranscriptomics_1.fastq.gz -reads data/set4_mate_pairs_metatranscriptomics_2.fastq.gz -num_alignments 1 -v -workdir run in /home/biocodz/sortmerna
   
   [process:1215] === Options processing starts ... ===
   
   Found value: /home/biocodz/sortmerna/dist/bin/sortmerna
   Found flag: -ref
   Found value: data/rRNA_databases/silva-euk-28s-id98.fasta of previous flag: -ref
   Found flag: -ref
   Found value: data/rRNA_databases/silva-euk-18s-id95.fasta of previous flag: -ref
   Found flag: -ref
   Found value: data/rRNA_databases/silva-bac-23s-id98.fasta of previous flag: -ref
   Found flag: -ref
   Found value: data/rRNA_databases/silva-bac-16s-id90.fasta of previous flag: -ref
   Found flag: -ref
   Found value: data/rRNA_databases/silva-arc-23s-id98.fasta of previous flag: -ref
   Found flag: -ref
   Found value: data/rRNA_databases/silva-arc-16s-id95.fasta of previous flag: -ref
   Found flag: -ref
   Found value: data/rRNA_databases/rfam-5s-database-id98.fasta of previous flag: -ref
   Found flag: -ref
   Found value: data/rRNA_databases/rfam-5.8s-database-id98.fasta of previous flag: -ref
   Found flag: -reads
   Found value: data/set4_mate_pairs_metatranscriptomics_1.fastq.gz of previous flag: -reads
   Found flag: -reads
   Found value: data/set4_mate_pairs_metatranscriptomics_2.fastq.gz of previous flag: -reads
   Found flag: -num_alignments
   Found value: 1 of previous flag: -num_alignments
   Found flag: -v
   Previous flag: -v is Boolean. Setting to True
   Found flag: -workdir
   Found value: run of previous flag: -workdir
   [opt_workdir:1027] Using WORKDIR ["/home/biocodz/sortmerna/run" as specified
   [process:1298] Processing option: num_alignments with value: 1
   [process:1298] Processing option: reads with value: data/set4_mate_pairs_metatranscriptomics_1.fastq.gz
   [opt_reads:74] Processing reads file [1] out of total [2] files
   [process:1298] Processing option: reads with value: data/set4_mate_pairs_metatranscriptomics_2.fastq.gz
   [opt_reads:74] Processing reads file [2] out of total [2] files
   [process:1298] Processing option: ref with value: data/rRNA_databases/silva-euk-28s-id98.fasta
   [opt_ref:168] Processing reference [1] out of total [8] references
   [opt_ref:223] File ["/home/biocodz/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta"] exists and is readable
   [process:1298] Processing option: ref with value: data/rRNA_databases/silva-euk-18s-id95.fasta
   [opt_ref:168] Processing reference [2] out of total [8] references
   [opt_ref:223] File ["/home/biocodz/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta"] exists and is readable
   [process:1298] Processing option: ref with value: data/rRNA_databases/silva-bac-23s-id98.fasta
   [opt_ref:168] Processing reference [3] out of total [8] references
   [opt_ref:223] File ["/home/biocodz/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta"] exists and is readable
   [process:1298] Processing option: ref with value: data/rRNA_databases/silva-bac-16s-id90.fasta
   [opt_ref:168] Processing reference [4] out of total [8] references
   [opt_ref:223] File ["/home/biocodz/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta"] exists and is readable
   [process:1298] Processing option: ref with value: data/rRNA_databases/silva-arc-23s-id98.fasta
   [opt_ref:168] Processing reference [5] out of total [8] references
   [opt_ref:223] File ["/home/biocodz/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta"] exists and is readable
   [process:1298] Processing option: ref with value: data/rRNA_databases/silva-arc-16s-id95.fasta
   [opt_ref:168] Processing reference [6] out of total [8] references
   [opt_ref:223] File ["/home/biocodz/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta"] exists and is readable
   [process:1298] Processing option: ref with value: data/rRNA_databases/rfam-5s-database-id98.fasta
   [opt_ref:168] Processing reference [7] out of total [8] references
   [opt_ref:223] File ["/home/biocodz/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta"] exists and is readable
   [process:1298] Processing option: ref with value: data/rRNA_databases/rfam-5.8s-database-id98.fasta
   [opt_ref:168] Processing reference [8] out of total [8] references
   [opt_ref:223] File ["/home/biocodz/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta"] exists and is readable
   [process:1298] Processing option: v with value:
   
   [process:1318] === Options processing done ===
   
   [test_kvdb_path:1145] Key-value DB location ("/home/biocodz/sortmerna/run/kvdb")
   [test_kvdb_path:1182] Database ("/home/biocodz/sortmerna/run/kvdb") will be created
   [validate:1340] No output format has been chosen (fastx/sam/blast/otu_map). Using default 'blast'
   
     Program:      SortMeRNA version 4.0.1
     Copyright:    2016-2019 Clarity Genomics BVBA:
                   Turnhoutseweg 30, 2340 Beerse, Belgium
                   2014-2016 Knight Lab:
                   Department of Pediatrics, UCSD, La Jolla
                   2012-2014 Bonsai Bioinformatics Research Group:
                   LIFL, University Lille 1, CNRS UMR 8022, INRIA Nord-Europe
     Disclaimer:   SortMeRNA comes with ABSOLUTELY NO WARRANTY; without even the
                   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
                   See the GNU Lesser General Public License for more details.
     Contributors: Jenya Kopylova   jenya.kopylov@gmail.com
                   Laurent Noé      laurent.noe@lifl.fr
                   Pierre Pericard  pierre.pericard@lifl.fr
                   Daniel McDonald  wasade@gmail.com
                   Mikaël Salson    mikael.salson@lifl.fr
                   Hélène Touzet    helene.touzet@lifl.fr
                   Rob Knight       robknight@ucsd.edu
   
   [main:72] Running command:
   /home/biocodz/sortmerna/dist/bin/sortmerna -ref data/rRNA_databases/silva-euk-28s-id98.fasta -ref data/rRNA_databases/silva-euk-18s-id95.fasta -ref data/rRNA_databases/silva-bac-23s-id98.fasta -ref data/rRNA_databases/silva-bac-16s-id90.fasta -ref data/rRNA_databases/silva-arc-23s-id98.fasta -ref data/rRNA_databases/silva-arc-16s-id95.fasta -ref data/rRNA_databases/rfam-5s-database-id98.fasta -ref data/rRNA_databases/rfam-5.8s-database-id98.fasta -reads data/set4_mate_pairs_metatranscriptomics_1.fastq.gz -reads data/set4_mate_pairs_metatranscriptomics_2.fastq.gz -num_alignments 1 -v -workdir run
   [Index:83] Index file [run/idx/1845323523482939374.bursttrie_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/1845323523482939374.pos_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/1845323523482939374.kmer_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/1845323523482939374.stats] already exists and is not empty.
   [Index:83] Index file [run/idx/2700646386527218729.bursttrie_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/2700646386527218729.pos_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/2700646386527218729.kmer_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/2700646386527218729.stats] already exists and is not empty.
   [Index:83] Index file [run/idx/17299952793705614139.bursttrie_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/17299952793705614139.pos_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/17299952793705614139.kmer_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/17299952793705614139.stats] already exists and is not empty.
   [Index:83] Index file [run/idx/15734375058464002811.bursttrie_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/15734375058464002811.pos_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/15734375058464002811.kmer_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/15734375058464002811.stats] already exists and is not empty.
   [Index:83] Index file [run/idx/3400685301612210653.bursttrie_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/3400685301612210653.pos_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/3400685301612210653.kmer_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/3400685301612210653.stats] already exists and is not empty.
   [Index:83] Index file [run/idx/3436099190853847617.bursttrie_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/3436099190853847617.pos_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/3436099190853847617.kmer_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/3436099190853847617.stats] already exists and is not empty.
   [Index:83] Index file [run/idx/13019673092862722585.bursttrie_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/13019673092862722585.pos_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/13019673092862722585.kmer_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/13019673092862722585.stats] already exists and is not empty.
   [Index:83] Index file [run/idx/2169995244134016533.bursttrie_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/2169995244134016533.pos_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/2169995244134016533.kmer_0.dat] already exists and is not empty.
   [Index:83] Index file [run/idx/2169995244134016533.stats] already exists and is not empty.
   [Index:104] Found 32 non-empty index files. Skipping indexing.
   [Index:105] TODO: a better validation using an index descriptor to decide on indexing
   [calculate:107] Starting statistics calculation on file: 'data/set4_mate_pairs_metatranscriptomics_1.fastq.gz'  ...   [inflatez:169] xINFO: inflateEnd status is 0
   [calculate:225] Done statistics on file. Elapsed time: 0.02 sec. all_reads_count= 5000
   [calculate:107] Starting statistics calculation on file: 'data/set4_mate_pairs_metatranscriptomics_2.fastq.gz'  ...   [inflatez:169] xINFO: inflateEnd status is 0
   [calculate:225] Done statistics on file. Elapsed time: 0.02 sec. all_reads_count= 10000
   [store_to_db:421] Stored Reads statistics to DB:
       min_read_len= 100 max_read_len= 100 all_reads_count= 10000 all_reads_len= 1000000 total_reads_mapped= 0 total_reads_mapped_cov= 0 reads_matched_per_db= TODO is_total_reads_mapped_cov= 0 is_stats_calc= 0
   
   
   [align:358] ==== Starting alignment ====
   
   [align:368] Using default number of Processor threads equals num CPU cores: 8
   Number of cores: 8 Read threads:  1 Write threads: 1 Processor threads: 8
   [ThreadPool:36] initialized Pool with: [10] threads
   
   [ReadsQueue:57] [read_queue] created with [1] Pushers
   [ReadsQueue:57] [write_queue] created with [8] Pushers
   [Refstats:32] Index Statistics calculation Start ...[Refstats:42] Done. Time elapsed: 4.00 sec
   
   [align:408] Loading index 0 part 1/1 ... done [2.69] sec
   [align:421] Loading references  ... done [0.12] sec
   [write:19] Writer writer_0 thread 139949152134912 started
   Processor proc_1 thread 139949277959936 started
   Processor proc_2 thread 139949286352640 started
   Processor proc_3 thread 139949982586624 started
   Processor proc_4 thread 139949261174528 started
   Processor proc_0 thread 139949160527616 started
   Processor proc_5 thread 139949135349504 started
   Processor proc_7 thread 139949294745344 started
   Processor proc_6 thread 139949143742208 started
   [run:70] thread: 139949269567232 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949269567232] pushers: [0]
   [run:113] thread: 139949269567232 done. Elapsed time: 0.23 sec Reads added: 10000 Num aligned reads (passing E-value): 0 readQueue.size: 25
   [threadEntry:108] number of running_threads= 9 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949135349504] pushers: [7]
   [decrPushers:170] id: [write_queue] thread: [139949277959936] pushers: [4]
   [decrPushers:170] id: [write_queue] thread: [139949261174528] pushers: [6]
   [decrPushers:170] id: [write_queue] thread: [139949160527616] pushers: [2]
   [run:95] Processor proc_0 thread 139949160527616 done. Processed 996 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 0
   [threadEntry:108] number of running_threads= 8 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949286352640] pushers: [3]
   [run:95] Processor proc_2 thread 139949286352640 done. Processed 1076 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 0
   [threadEntry:108] number of running_threads= 7 jobs queue is empty= 1
   [run:95] Processor proc_1 thread 139949277959936 done. Processed 1199 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 0
   [write:50] writer_0 thread 139949152134912 done. Elapsed time: 0.23 s Reads written: 10000 Num aligned reads (passing E-value):0
   [threadEntry:108] number of running_threads= 5 jobs queue is empty= 1
   [run:95] Processor proc_4 thread 139949261174528 done. Processed 1242 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 0
   [threadEntry:108] number of running_threads= 4 jobs queue is empty= 1
   [threadEntry:108] number of running_threads= 6 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949294745344] pushers: [1]
   [run:95] Processor proc_7 thread 139949294745344 done. Processed 1022 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 0
   [threadEntry:108] number of running_threads= 3 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949143742208] pushers: [5]
   [run:95] Processor proc_6 thread 139949143742208 done. Processed 1438 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 0
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [run:95] Processor proc_5 thread 139949135349504 done. Processed 1046 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 0
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949982586624] pushers: [0]
   [run:95] Processor proc_3 thread 139949982586624 done. Processed 1981 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [8]
   [reset:137] [read_queue] pushers: [1]
   [align:461] Done index 0 Part: 1 Time: 0.36 sec
   
   [align:408] Loading index 1 part 1/1 ... done [2.32] sec
   [align:421] Loading references  ... done [0.12] sec
   [write:19] Writer writer_0 thread 139949261174528 started
   Processor proc_4 thread 139949135349504 started
   Processor proc_2 thread 139949277959936 started
   Processor proc_3 thread 139949160527616 started
   Processor proc_7 thread 139949982586624 started
   Processor proc_5 thread 139949294745344 started
   Processor proc_6 thread 139949286352640 started
   Processor proc_1 thread 139949269567232 started
   Processor proc_0 thread 139949152134912 started
   [run:70] thread: 139949143742208 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949143742208] pushers: [0]
   [run:113] thread: 139949143742208 done. Elapsed time: 0.77 sec Reads added: 10000 Num aligned reads (passing E-value): 0 readQueue.size: 100
   [threadEntry:108] number of running_threads= 9 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949294745344] pushers: [7]
   [run:95] Processor proc_5 thread 139949294745344 done. Processed 1188 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 169
   [threadEntry:108] number of running_threads= 8 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949135349504] pushers: [6]
   [run:95] Processor proc_4 thread 139949135349504 done. Processed 2020 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 266
   [threadEntry:108] number of running_threads= 7 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949160527616] pushers: [5]
   [run:95] Processor proc_3 thread 139949160527616 done. Processed 1129 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 147
   [threadEntry:108] number of running_threads= 6 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949277959936] pushers: [4]
   [run:95] Processor proc_2 thread 139949277959936 done. Processed 509 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 74
   [threadEntry:108] number of running_threads= 5 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949269567232] pushers: [3]
   [run:95] Processor proc_1 thread 139949269567232 done. Processed 1279 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 184
   [threadEntry:108] number of running_threads= 4 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [2]
   [run:95] Processor proc_0 thread 139949152134912 done. Processed 1431 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 214
   [threadEntry:108] number of running_threads= 3 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949982586624] pushers: [1]
   [run:95] Processor proc_7 thread 139949982586624 done. Processed 1518 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 229
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949286352640] pushers: [0]
   [run:95] Processor proc_6 thread 139949286352640 done. Processed 926 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 136
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949261174528 done. Elapsed time: 0.84 s Reads written: 10000 Num aligned reads (passing E-value):1419
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [8]
   [reset:137] [read_queue] pushers: [1]
   [align:461] Done index 1 Part: 1 Time: 0.96 sec
   
   [align:408] Loading index 2 part 1/1 ... done [2.48] sec
   [align:421] Loading references  ... done [0.12] sec
   [write:19] Writer writer_0 thread 139949294745344 started
   [run:70] thread: 139949135349504 started
   Processor proc_0 thread 139949143742208 started
   Processor proc_1 thread 139949269567232 started
   Processor proc_4 thread 139949286352640 started
   Processor proc_5 thread 139949261174528 started
   Processor proc_6 thread 139949152134912 started
   Processor proc_2 thread 139949982586624 started
   Processor proc_7 thread 139949277959936 started
   Processor proc_3 thread 139949160527616 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949135349504] pushers: [0]
   [run:113] thread: 139949135349504 done. Elapsed time: 0.18 sec Reads added: 10000 Num aligned reads (passing E-value): 1419 readQueue.size: 2
   [threadEntry:108] number of running_threads= 9 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949160527616] pushers: [7]
   [run:95] Processor proc_3 thread 139949160527616 done. Processed 831 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 137
   [threadEntry:108] number of running_threads= 8 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949269567232] pushers: [6]
   [run:95] Processor proc_1 thread 139949269567232 done. Processed 1588 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 227
   [threadEntry:108] number of running_threads= 7 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [5]
   [run:95] Processor proc_6 thread 139949152134912 done. Processed 1584 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 231
   [threadEntry:108] number of running_threads= 6 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949286352640] pushers: [4]
   [run:95] Processor proc_4 thread 139949286352640 done. Processed 1642 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 226
   [threadEntry:108] number of running_threads= 5 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949277959936] pushers: [3]
   [run:95] Processor proc_7 thread 139949277959936 done. Processed 797 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 107
   [threadEntry:108] number of running_threads= 4 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949143742208] pushers: [1]
   [run:95] Processor proc_0 thread 139949143742208 done. Processed 980 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 143
   [threadEntry:108] number of running_threads= 3 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949261174528] pushers: [1]
   [run:95] Processor proc_5 thread 139949261174528 done. Processed 970 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 144
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949982586624] pushers: [0]
   [run:95] Processor proc_2 thread 139949982586624 done. Processed 1608 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 204
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949294745344 done. Elapsed time: 0.18 s Reads written: 10000 Num aligned reads (passing E-value):1419
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [8]
   [reset:137] [read_queue] pushers: [1]
   [align:461] Done index 2 Part: 1 Time: 0.30 sec
   
   [align:408] Loading index 3 part 1/1 ... done [2.70] sec
   [align:421] Loading references  ... done [0.11] sec
   [run:70] thread: 139949135349504 started
   Processor proc_0 thread 139949143742208 started
   Processor proc_5 thread 139949160527616 started
   [write:19] Writer writer_0 thread 139949152134912 started
   Processor proc_3 thread 139949982586624 started
   Processor proc_4 thread 139949294745344 started
   Processor proc_1 thread 139949269567232 started
   Processor proc_2 thread 139949286352640 started
   Processor proc_6 thread 139949277959936 started
   Processor proc_7 thread 139949261174528 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949135349504] pushers: [0]
   [run:113] thread: 139949135349504 done. Elapsed time: 1.82 sec Reads added: 10000 Num aligned reads (passing E-value): 1419 readQueue.size: 100
   [threadEntry:108] number of running_threads= 9 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949269567232] pushers: [7]
   [run:95] Processor proc_1 thread 139949269567232 done. Processed 1229 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 730
   [threadEntry:108] number of running_threads= 8 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949261174528] pushers: [6]
   [run:95] Processor proc_7 thread 139949261174528 done. Processed 1404 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 857
   [threadEntry:108] number of running_threads= 7 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949982586624] pushers: [5]
   [run:95] Processor proc_3 thread 139949982586624 done. Processed 1376 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 843
   [threadEntry:108] number of running_threads= 6 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949143742208] pushers: [4]
   [run:95] Processor proc_0 thread 139949143742208 done. Processed 1236 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 746
   [threadEntry:108] number of running_threads= 5 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949294745344] pushers: [3]
   [run:95] Processor proc_4 thread 139949294745344 done. Processed 1192 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 726
   [threadEntry:108] number of running_threads= 4 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949160527616] pushers: [2]
   [run:95] Processor proc_5 thread 139949160527616 done. Processed 1322 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 761
   [threadEntry:108] number of running_threads= 3 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949286352640] pushers: [1]
   [run:95] Processor proc_2 thread 139949286352640 done. Processed 1141 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 623
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949277959936] pushers: [0]
   [write:50] writer_0 thread 139949152134912 done. Elapsed time: 1.85 s Reads written: 10000 Num aligned reads (passing E-value):5942
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:95] Processor proc_6 thread 139949277959936 done. Processed 1100 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 656
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [8]
   [reset:137] [read_queue] pushers: [1]
   [align:461] Done index 3 Part: 1 Time: 2.02 sec
   
   [align:408] Loading index 4 part 1/1 ... done [0.31] sec
   [align:421] Loading references  ... done [0.01] sec
   [run:70] thread: 139949135349504 started
   [write:19] Writer writer_0 thread 139949286352640 started
   Processor proc_0 thread 139949294745344 started
   Processor proc_1 thread 139949160527616 started
   Processor proc_2 thread 139949269567232 started
   Processor proc_4 thread 139949143742208 started
   Processor proc_7 thread 139949982586624 started
   Processor proc_6 thread 139949277959936 started
   Processor proc_3 thread 139949152134912 started
   Processor proc_5 thread 139949261174528 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949135349504] pushers: [0]
   [run:113] thread: 139949135349504 done. Elapsed time: 0.40 sec Reads added: 10000 Num aligned reads (passing E-value): 5942 readQueue.size: 0
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [6]
   [decrPushers:170] id: [write_queue] thread: [139949143742208] pushers: [4]
   [decrPushers:170] id: [write_queue] thread: [139949277959936] pushers: [3]
   [threadEntry:108] number of running_threads= 9 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949294745344] pushers: [5]
   [run:95] Processor proc_0 thread 139949294745344 done. Processed 1470 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 868
   [threadEntry:108] number of running_threads= 8 jobs queue is empty= 1
   [run:95] Processor proc_4 thread 139949143742208 done. Processed 1465 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 850
   [threadEntry:108] number of running_threads= 7 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949269567232] pushers: [5]
   [run:95] Processor proc_2 thread 139949269567232 done. Processed 1441 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 853
   [threadEntry:108] number of running_threads= 6 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949982586624] pushers: [0]
   [run:95] Processor proc_7 thread 139949982586624 done. Processed 1296 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 784
   [threadEntry:108] number of running_threads= 5 jobs queue is empty= 1
   [write:50] writer_0 thread 139949286352640 done. Elapsed time: 0.40 s Reads written: 10000 Num aligned reads (passing E-value):5942
   [threadEntry:108] number of running_threads= 4 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949160527616] pushers: [1]
   [run:95] Processor proc_1 thread 139949160527616 done. Processed 1441 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 863
   [threadEntry:108] number of running_threads= 3 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949261174528] pushers: [1]
   [run:95] Processor proc_5 thread 139949261174528 done. Processed 837 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 481
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [run:95] Processor proc_6 thread 139949277959936 done. Processed 996 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 589
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:95] Processor proc_3 thread 139949152134912 done. Processed 1054 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 654
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [8]
   [reset:137] [read_queue] pushers: [1]
   [align:461] Done index 4 Part: 1 Time: 0.42 sec
   
   [align:408] Loading index 5 part 1/1 ... done [0.65] sec
   [align:421] Loading references  ... done [0.02] sec
   Processor proc_0 thread 139949143742208 started
   Processor proc_1 thread 139949135349504 started
   Processor proc_2 thread 139949286352640 started
   [write:19] Writer writer_0 thread 139949294745344 started
   Processor proc_3 thread 139949261174528 started
   Processor proc_4 thread 139949152134912 started
   Processor proc_5 thread 139949982586624 started
   Processor proc_6 thread 139949277959936 started
   [run:70] thread: 139949269567232 started
   Processor proc_7 thread 139949160527616 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949269567232] pushers: [0]
   [run:113] thread: 139949269567232 done. Elapsed time: 0.18 sec Reads added: 10000 Num aligned reads (passing E-value): 5942 readQueue.size: 29
   [threadEntry:108] number of running_threads= 9 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949143742208] pushers: [7]
   [run:95] Processor proc_0 thread 139949143742208 done. Processed 1596 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 960
   [threadEntry:108] number of running_threads= 8 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949277959936] pushers: [5]
   [decrPushers:170] id: [write_queue] thread: [139949261174528] pushers: [3]
   [run:95] Processor proc_3 thread 139949261174528 done. Processed 964 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 593
   [threadEntry:108] number of running_threads= 7 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949135349504] pushers: [4]
   [run:95] Processor proc_1 thread 139949135349504 done. Processed 808 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 455
   [threadEntry:108] number of running_threads= 6 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949160527616] pushers: [2]
   [run:95] Processor proc_7 thread 139949160527616 done. Processed 1395 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 831
   [threadEntry:108] number of running_threads= 5 jobs queue is empty= 1
   [run:95] Processor proc_6 thread 139949277959936 done. Processed 922 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 543
   [threadEntry:108] number of running_threads= 4 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [6]
   [run:95] Processor proc_4 thread 139949152134912 done. Processed 961 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 561
   [threadEntry:108] number of running_threads= 3 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949982586624] pushers: [0]
   [run:95] Processor proc_5 thread 139949982586624 done. Processed 1581 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 968
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949286352640] pushers: [1]
   [run:95] Processor proc_2 thread 139949286352640 done. Processed 1773 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 1033
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949294745344 done. Elapsed time: 0.18 s Reads written: 10000 Num aligned reads (passing E-value):5944
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [8]
   [reset:137] [read_queue] pushers: [1]
   [align:461] Done index 5 Part: 1 Time: 0.22 sec
   
   [align:408] Loading index 6 part 1/1 ... done [0.84] sec
   [align:421] Loading references  ... done [0.06] sec
   [write:19] Writer writer_0 thread 139949160527616 started
   Processor proc_0 thread 139949277959936 started
   Processor proc_1 thread 139949982586624 started
   Processor proc_5 thread 139949261174528 started
   Processor proc_2 thread 139949143742208 started
   Processor proc_6 thread 139949294745344 started
   Processor proc_7 thread 139949269567232 started
   Processor proc_3 thread 139949135349504 started
   [run:70] thread: 139949286352640 started
   Processor proc_4 thread 139949152134912 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949286352640] pushers: [0]
   [decrPushers:170] id: [write_queue] thread: [139949982586624] pushers: [7]
   [run:95] Processor proc_1 thread 139949982586624 done. Processed 1369 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 818
   [threadEntry:108] number of running_threads= 9 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949143742208] pushers: [3]
   [decrPushers:170] id: [write_queue] thread: [139949261174528] pushers: [2]
   [run:95] Processor proc_5 thread 139949261174528 done. Processed 949 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 525
   [threadEntry:108] number of running_threads= 8 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949135349504] pushers: [1]
   [run:113] thread: 139949286352640 done. Elapsed time: 0.30 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 0
   [threadEntry:108] number of running_threads= 7 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [0]
   [run:95] Processor proc_4 thread 139949152134912 done. Processed 1459 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 907
   [threadEntry:108] number of running_threads= 6 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949277959936] pushers: [6]
   [decrPushers:170] id: [write_queue] thread: [139949269567232] pushers: [5]
   [run:95] Processor proc_7 thread 139949269567232 done. Processed 1362 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 835
   [threadEntry:108] number of running_threads= 5 jobs queue is empty= 1
   [run:95] Processor proc_3 thread 139949135349504 done. Processed 952 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 547
   [threadEntry:108] number of running_threads= 4 jobs queue is empty= 1
   [run:95] Processor proc_0 thread 139949277959936 done. Processed 1265 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 737
   [threadEntry:108] number of running_threads= 3 jobs queue is empty= 1
   [run:95] Processor proc_2 thread 139949143742208 done. Processed 1326 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 798
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949294745344] pushers: [4]
   [run:95] Processor proc_6 thread 139949294745344 done. Processed 1318 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 777
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949160527616 done. Elapsed time: 0.30 s Reads written: 10000 Num aligned reads (passing E-value):5944
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [8]
   [reset:137] [read_queue] pushers: [1]
   [align:461] Done index 6 Part: 1 Time: 0.36 sec
   
   [align:408] Loading index 7 part 1/1 ... done [0.29] sec
   [align:421] Loading references  ... done [0.01] sec
   [run:70] thread: 139949286352640 started
   [write:19] Writer writer_0 thread 139949982586624 started
   Processor proc_0 thread 139949261174528 started
   Processor proc_1 thread 139949152134912 started
   Processor proc_2 thread 139949135349504 started
   Processor proc_3 thread 139949269567232 started
   Processor proc_5 thread 139949277959936 started
   Processor proc_4 thread 139949143742208 started
   Processor proc_6 thread 139949160527616 started
   Processor proc_7 thread 139949294745344 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949286352640] pushers: [0]
   [decrPushers:170] id: [write_queue] thread: [139949294745344] pushers: [7]
   [run:95] Processor proc_7 thread 139949294745344 done. Processed 1364 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 791
   [run:113] thread: 139949286352640 done. Elapsed time: 0.29 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 0
   [decrPushers:170] id: [write_queue] thread: [139949143742208] pushers: [6]
   [run:95] Processor proc_4 thread 139949143742208 done. Processed 1474 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 839
   [threadEntry:108] number of running_threads= 8 jobs queue is empty= 1
   [threadEntry:108] number of running_threads= 9 jobs queue is empty= 1
   [threadEntry:108] number of running_threads= 7 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949135349504] pushers: [3]
   [run:95] Processor proc_2 thread 139949135349504 done. Processed 1082 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 627
   [threadEntry:108] number of running_threads= 6 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949261174528] pushers: [4]
   [run:95] Processor proc_0 thread 139949261174528 done. Processed 1142 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 697
   [decrPushers:170] id: [write_queue] thread: [139949160527616] pushers: [1]
   [run:95] Processor proc_6 thread 139949160527616 done. Processed 828 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 487
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [5]
   [decrPushers:170] id: [write_queue] thread: [139949269567232] pushers: [2]
   [threadEntry:108] number of running_threads= 5 jobs queue is empty= 1
   [run:95] Processor proc_3 thread 139949269567232 done. Processed 1436 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 884
   [threadEntry:108] number of running_threads= 3 jobs queue is empty= 1
   [run:95] Processor proc_1 thread 139949152134912 done. Processed 1519 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 909
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949277959936] pushers: [0]
   [write:50] writer_0 thread 139949982586624 done. Elapsed time: 0.29 s Reads written: 10000 Num aligned reads (passing E-value):5944
   [run:95] Processor proc_5 thread 139949277959936 done. Processed 1155 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 710
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [threadEntry:108] number of running_threads= 4 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [8]
   [reset:137] [read_queue] pushers: [1]
   [align:461] Done index 7 Part: 1 Time: 0.31 sec
   
   [align:468] ==== Done alignment ====
   
   [store_to_db:421] Stored Reads statistics to DB:
       min_read_len= 100 max_read_len= 100 all_reads_count= 10000 all_reads_len= 1000000 total_reads_mapped= 5944 total_reads_mapped_cov= 5944 reads_matched_per_db= TODO is_total_reads_mapped_cov= 1 is_stats_calc= 0
   
   [~ReadsQueue:68] Destructor called on write_queue  recs.size= 0 pushed: 80000  popped: 80000
   [~ReadsQueue:68] Destructor called on read_queue  recs.size= 0 pushed: 80000  popped: 80000
   Thread  139949143742208 job done
   Thread  139949294745344 job done
   Thread  139949286352640 job done
   Thread  139949135349504 job done
   Thread  139949261174528 job done
   Thread  139949269567232 job done
   Thread  139949152134912 job done
   Thread  139949982586624 job done
   Thread  139949277959936 job done
   Thread  139949160527616 job done
   
   [postProcess:206] ==== Starting Post-processing (alignment statistics report) ====
   
   [ThreadPool:36] initialized Pool with: [3] threads
   
   [ReadsQueue:57] [read_queue] created with [1] Pushers
   [ReadsQueue:57] [write_queue] created with [1] Pushers
   [postProcess:217] Restored Readstats from DB:
       min_read_len= 100 max_read_len= 100 all_reads_count= 10000 all_reads_len= 1000000 total_reads_mapped= 5944 total_reads_mapped_cov= 5944 reads_matched_per_db= TODO is_total_reads_mapped_cov= 1 is_stats_calc= 0
   
   [Refstats:32] Index Statistics calculation Start ...[Refstats:42] Done. Time elapsed: 3.77 sec
   
   [postProcess:236] Loading reference 0 part 1/1  ... done [0.02 sec]
   [run:70] thread: 139949152134912 started
   [write:19] Writer writer_0 thread 139949143742208 started
   [run:111] PostProcessor postproc_0 thread 139949135349504 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949152134912] pushers: [0]
   [decrPushers:170] id: [write_queue] thread: [139949135349504] pushers: [0]
   [run:141] postproc_0 thread 139949135349504 done. Processed 10000 reads. count_reads_aligned: 5944
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [write:50] writer_0 thread 139949143742208 done. Elapsed time: 0.16 s Reads written: 0 Num aligned reads (passing E-value):0
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:113] thread: 139949152134912 done. Elapsed time: 0.16 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [read_queue] pushers: [1]
   [reset:137] [write_queue] pushers: [1]
   [postProcess:278] Done reference 0 Part: 1 Time: 0.16 sec
   
   [postProcess:236] Loading reference 1 part 1/1  ... done [0.02 sec]
   [run:70] thread: 139949143742208 started
   [write:19] Writer writer_0 thread 139949135349504 started
   [run:111] PostProcessor postproc_0 thread 139949152134912 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949143742208] pushers: [0]
   [run:113] thread: 139949143742208 done. Elapsed time: 0.14 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 2
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [0]
   [run:141] postproc_0 thread 139949152134912 done. Processed 10000 reads. count_reads_aligned: 5944
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949135349504 done. Elapsed time: 0.14 s Reads written: 0 Num aligned reads (passing E-value):0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [read_queue] pushers: [1]
   [reset:137] [write_queue] pushers: [1]
   [postProcess:278] Done reference 1 Part: 1 Time: 0.14 sec
   
   [postProcess:236] Loading reference 2 part 1/1  ... done [0.01 sec]
   [write:19] Writer writer_0 thread 139949135349504 started
   [run:70] thread: 139949143742208 started
   [run:111] PostProcessor postproc_0 thread 139949152134912 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949143742208] pushers: [0]
   [run:113] thread: 139949143742208 done. Elapsed time: 0.14 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 2
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [0]
   [run:141] postproc_0 thread 139949152134912 done. Processed 10000 reads. count_reads_aligned: 5944
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949135349504 done. Elapsed time: 0.14 s Reads written: 0 Num aligned reads (passing E-value):0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [read_queue] pushers: [1]
   [reset:137] [write_queue] pushers: [1]
   [postProcess:278] Done reference 2 Part: 1 Time: 0.14 sec
   
   [postProcess:236] Loading reference 3 part 1/1  ... done [0.04 sec]
   [run:111] PostProcessor postproc_0 thread 139949135349504 started
   [run:70] thread: 139949143742208 started
   [write:19] Writer writer_0 thread 139949152134912 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949143742208] pushers: [0]
   [run:113] thread: 139949143742208 done. Elapsed time: 0.13 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 2
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949135349504] pushers: [0]
   [run:141] postproc_0 thread 139949135349504 done. Processed 10000 reads. count_reads_aligned: 5944
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949152134912 done. Elapsed time: 0.13 s Reads written: 0 Num aligned reads (passing E-value):0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [read_queue] pushers: [1]
   [reset:137] [write_queue] pushers: [1]
   [postProcess:278] Done reference 3 Part: 1 Time: 0.13 sec
   
   [postProcess:236] Loading reference 4 part 1/1  ... done [0.00 sec]
   [write:19] Writer writer_0 thread 139949143742208 started
   [run:70] thread: 139949135349504 started
   [run:111] PostProcessor postproc_0 thread 139949152134912 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949135349504] pushers: [0]
   [run:113] thread: 139949135349504 done. Elapsed time: 0.13 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 2
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [0]
   [run:141] postproc_0 thread 139949152134912 done. Processed 10000 reads. count_reads_aligned: 5944
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949143742208 done. Elapsed time: 0.13 s Reads written: 0 Num aligned reads (passing E-value):0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [read_queue] pushers: [1]
   [reset:137] [write_queue] pushers: [1]
   [postProcess:278] Done reference 4 Part: 1 Time: 0.13 sec
   
   [postProcess:236] Loading reference 5 part 1/1  ... done [0.01 sec]
   [run:70] thread: 139949135349504 started
   [write:19] Writer writer_0 thread 139949152134912 started
   [run:111] PostProcessor postproc_0 thread 139949143742208 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949135349504] pushers: [0]
   [run:113] thread: 139949135349504 done. Elapsed time: 0.13 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 1
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949143742208] pushers: [0]
   [run:141] postproc_0 thread 139949143742208 done. Processed 10000 reads. count_reads_aligned: 5944
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949152134912 done. Elapsed time: 0.13 s Reads written: 0 Num aligned reads (passing E-value):0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [read_queue] pushers: [1]
   [reset:137] [write_queue] pushers: [1]
   [postProcess:278] Done reference 5 Part: 1 Time: 0.13 sec
   
   [postProcess:236] Loading reference 6 part 1/1  ... done [0.05 sec]
   [run:70] thread: 139949135349504 started
   [write:19] Writer writer_0 thread 139949143742208 started
   [run:111] PostProcessor postproc_0 thread 139949152134912 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949135349504] pushers: [0]
   [run:113] thread: 139949135349504 done. Elapsed time: 0.13 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 2
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949152134912] pushers: [0]
   [run:141] postproc_0 thread 139949152134912 done. Processed 10000 reads. count_reads_aligned: 5944
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949143742208 done. Elapsed time: 0.13 s Reads written: 0 Num aligned reads (passing E-value):0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [read_queue] pushers: [1]
   [reset:137] [write_queue] pushers: [1]
   [postProcess:278] Done reference 6 Part: 1 Time: 0.13 sec
   
   [postProcess:236] Loading reference 7 part 1/1  ... done [0.02 sec]
   [write:19] Writer writer_0 thread 139949152134912 started
   [run:111] PostProcessor postproc_0 thread 139949143742208 started
   [run:70] thread: 139949135349504 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949135349504] pushers: [0]
   [run:113] thread: 139949135349504 done. Elapsed time: 0.13 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 2
   [threadEntry:108] number of running_threads= 2 jobs queue is empty= 1
   [decrPushers:170] id: [write_queue] thread: [139949143742208] pushers: [0]
   [run:141] postproc_0 thread 139949143742208 done. Processed 10000 reads. count_reads_aligned: 5944
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [write:50] writer_0 thread 139949152134912 done. Elapsed time: 0.13 s Reads written: 0 Num aligned reads (passing E-value):0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [read_queue] pushers: [1]
   [reset:137] [write_queue] pushers: [1]
   [postProcess:278] Done reference 7 Part: 1 Time: 0.13 sec
   [postProcess:287] total_reads_denovo_clustering = 0
   [store_to_db:421] Stored Reads statistics to DB:
       min_read_len= 100 max_read_len= 100 all_reads_count= 10000 all_reads_len= 1000000 total_reads_mapped= 5944 total_reads_mapped_cov= 5944 reads_matched_per_db= TODO is_total_reads_mapped_cov= 1 is_stats_calc= 1
   
   [writeLog:779] Using Log file: run/out/aligned.log
   
   [postProcess:304] ==== Done Post-processing (alignment statistics report) ====
   
   [~ReadsQueue:68] Destructor called on write_queue  recs.size= 0 pushed: 0  popped: 0
   [~ReadsQueue:68] Destructor called on read_queue  recs.size= 0 pushed: 80000  popped: 80000
   Thread  139949135349504 job done
   Thread  139949143742208 job done
   Thread  139949152134912 job done
   
   [generateReports:938] === Report generation starts. Thread: 19971136 ===
   
   [ThreadPool:36] initialized Pool with: [2] threads
   
   [generateReports:946] Restored Readstats from DB: 1
   [ReadsQueue:57] [read_queue] created with [1] Pushers
   [ReadsQueue:57] [write_queue] created with [1] Pushers
   [Refstats:32] Index Statistics calculation Start ...[Refstats:42] Done. Time elapsed: 3.72 sec
   [generateReports:946] Restored Readstats from DB: 1
   
   [generateReports:964] Loading reference 0 part 1/1  ... done [0.02 sec]
   [run:154] Report Processor report_proc_0 thread 139949152134912 started
   [run:70] thread: 139949143742208 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949143742208] pushers: [0]
   [run:113] thread: 139949143742208 done. Elapsed time: 0.12 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 0
   [run:191] Report Processor report_proc_0 thread 139949152134912 done. Processed 10000 reads
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [1]
   [reset:137] [read_queue] pushers: [1]
   [generateReports:994] Done reference 0 Part: 1 Time: 0.13 sec
   [generateReports:994] Done reference 0 Part: 1 Time: 0.13 sec
   
   [generateReports:964] Loading reference 1 part 1/1  ... done [0.01 sec]
   [run:154] Report Processor report_proc_0 thread 139949152134912 started
   [run:70] thread: 139949143742208 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949143742208] pushers: [0]
   [run:113] thread: 139949143742208 done. Elapsed time: 0.13 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 1
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:191] Report Processor report_proc_0 thread 139949152134912 done. Processed 10000 reads
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [1]
   [reset:137] [read_queue] pushers: [1]
   [generateReports:994] Done reference 1 Part: 1 Time: 0.13 sec
   [generateReports:994] Done reference 1 Part: 1 Time: 0.13 sec
   
   [generateReports:964] Loading reference 2 part 1/1  ... done [0.01 sec]
   [run:70] thread: 139949152134912 started
   [run:154] Report Processor report_proc_0 thread 139949143742208 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949152134912] pushers: [0]
   [run:113] thread: 139949152134912 done. Elapsed time: 0.13 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 2
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:191] Report Processor report_proc_0 thread 139949143742208 done. Processed 10000 reads
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [1]
   [reset:137] [read_queue] pushers: [1]
   [generateReports:994] Done reference 2 Part: 1 Time: 0.13 sec
   [generateReports:994] Done reference 2 Part: 1 Time: 0.13 sec
   
   [generateReports:964] Loading reference 3 part 1/1  ... done [0.02 sec]
   [run:154] Report Processor report_proc_0 thread 139949152134912 started
   [run:70] thread: 139949143742208 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949143742208] pushers: [0]
   [run:191] Report Processor report_proc_0 thread 139949152134912 done. Processed 10000 reads
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:113] thread: 139949143742208 done. Elapsed time: 0.12 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [1]
   [reset:137] [read_queue] pushers: [1]
   [generateReports:994] Done reference 3 Part: 1 Time: 0.13 sec
   [generateReports:994] Done reference 3 Part: 1 Time: 0.13 sec
   
   [generateReports:964] Loading reference 4 part 1/1  ... done [0.00 sec]
   [run:70] thread: 139949152134912 started
   [run:154] Report Processor report_proc_0 thread 139949143742208 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949152134912] pushers: [0]
   [run:113] thread: 139949152134912 done. Elapsed time: 0.13 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 0
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:191] Report Processor report_proc_0 thread 139949143742208 done. Processed 10000 reads
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [1]
   [reset:137] [read_queue] pushers: [1]
   [generateReports:994] Done reference 4 Part: 1 Time: 0.13 sec
   [generateReports:994] Done reference 4 Part: 1 Time: 0.13 sec
   
   [generateReports:964] Loading reference 5 part 1/1  ... done [0.01 sec]
   [run:154] Report Processor report_proc_0 thread 139949143742208 started
   [run:70] thread: 139949152134912 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949152134912] pushers: [0]
   [run:191] Report Processor report_proc_0 thread 139949143742208 done. Processed 10000 reads
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:113] thread: 139949152134912 done. Elapsed time: 0.12 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 0
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [1]
   [reset:137] [read_queue] pushers: [1]
   [generateReports:994] Done reference 5 Part: 1 Time: 0.12 sec
   [generateReports:994] Done reference 5 Part: 1 Time: 0.12 sec
   
   [generateReports:964] Loading reference 6 part 1/1  ... done [0.05 sec]
   [run:70] thread: 139949152134912 started
   [run:154] Report Processor report_proc_0 thread 139949143742208 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949152134912] pushers: [0]
   [run:113] thread: 139949152134912 done. Elapsed time: 0.13 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 2
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:191] Report Processor report_proc_0 thread 139949143742208 done. Processed 10000 reads
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [1]
   [reset:137] [read_queue] pushers: [1]
   [generateReports:994] Done reference 6 Part: 1 Time: 0.13 sec
   [generateReports:994] Done reference 6 Part: 1 Time: 0.13 sec
   
   [generateReports:964] Loading reference 7 part 1/1  ... done [0.01 sec]
   [run:70] thread: 139949152134912 started
   [run:154] Report Processor report_proc_0 thread 139949143742208 started
   [inflatez:169] xINFO: inflateEnd status is 0
   [inflatez:169] xINFO: inflateEnd status is 0
   [decrPushers:170] id: [read_queue] thread: [139949152134912] pushers: [0]
   [run:113] thread: 139949152134912 done. Elapsed time: 0.12 sec Reads added: 10000 Num aligned reads (passing E-value): 5944 readQueue.size: 2
   [threadEntry:108] number of running_threads= 1 jobs queue is empty= 1
   [run:191] Report Processor report_proc_0 thread 139949143742208 done. Processed 10000 reads
   [threadEntry:108] number of running_threads= 0 jobs queue is empty= 1
   [reset:137] [write_queue] pushers: [1]
   [reset:137] [read_queue] pushers: [1]
   [generateReports:994] Done reference 7 Part: 1 Time: 0.12 sec
   
   [generateReports:1001] === Done Reports generation ===
   
   [~ReadsQueue:68] Destructor called on write_queue  recs.size= 0 pushed: 0  popped: 0
   [~ReadsQueue:68] Destructor called on read_queue  recs.size= 0 pushed: 80000  popped: 80000
   Thread  139949152134912 job done
   Thread  139949143742208 job done
   [closefiles:766] Flushed and closed
   [run] Run time: 31.779393434524536

Sample execution trace v4.3.2
=============================

::

   python /media/sf_a01_code/sortmerna/scripts/run.py --name t41 --envn LNX_VBox_Ubuntu_1604
   Current dir: /media/sf_a01_code/sortmerna/scripts
   [run.py:__main__] Using Environment configuration file: /media/sf_a01_code/sortmerna/scripts/env.jinja.yaml
   [run.py:__main__] Using Build configuration template: /media/sf_a01_code/sortmerna/scripts/test.jinja.yaml
   [run.py:__main__] using /home/biocodz/sortmerna/dist/bin/sortmerna
   [process_smr_opts] '-workdir' option was provided. Using workdir: [/home/biocodz/sortmerna/run]
   [process_smr_opts] '-workdir' option was provided. Using workdir: [/home/biocodz/sortmerna/run]
   [run.py:__main__] Removing KVDB dir: /home/biocodz/sortmerna/run/kvdb
   [run.py:__main__] Removing Aligned Output: /home/biocodz/sortmerna/run/out
   [run.py:__main__] Running t41: issue 231 8x1M
   [run] Running: /home/biocodz/sortmerna/dist/bin/sortmerna -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta -reads /media/sf_a01_data/bio/reads/rna/SRR1635864_1_2M.fq.gz -reads /media/sf_a01_data/bio/reads/rna/SRR1635864_2_2M.fq.gz -fastx -blast 1 cigar qcov -out2 -sout -other -v -threads 8 -index 2 -workdir /home/biocodz/sortmerna/run in /media/sf_a01_code/sortmerna
   [process:1372] === Options processing starts ... ===
   Found value: /home/biocodz/sortmerna/dist/bin/sortmerna
   Found flag: -ref
   Found value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta of previous flag: -ref
   Found flag: -ref
   Found value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta of previous flag: -ref
   Found flag: -ref
   Found value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta of previous flag: -ref
   Found flag: -ref
   Found value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta of previous flag: -ref
   Found flag: -ref
   Found value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta of previous flag: -ref
   Found flag: -ref
   Found value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta of previous flag: -ref
   Found flag: -ref
   Found value: /media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta of previous flag: -ref
   Found flag: -ref
   Found value: /media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta of previous flag: -ref
   Found flag: -reads
   Found value: /media/sf_a01_data/bio/reads/rna/SRR1635864_1_2M.fq.gz of previous flag: -reads
   Found flag: -reads
   Found value: /media/sf_a01_data/bio/reads/rna/SRR1635864_2_2M.fq.gz of previous flag: -reads
   Found flag: -fastx
   Previous flag: -fastx is Boolean. Setting to True
   Found flag: -blast
   Found value: 1 cigar qcov of previous flag: -blast
   Found flag: -out2
   Previous flag: -out2 is Boolean. Setting to True
   Found flag: -sout
   Previous flag: -sout is Boolean. Setting to True
   Found flag: -other
   Previous flag: -other is Boolean. Setting to True
   Found flag: -v
   Previous flag: -v is Boolean. Setting to True
   Found flag: -threads
   Found value: 8 of previous flag: -threads
   Found flag: -index
   Found value: 2 of previous flag: -index
   Found flag: -workdir
   Found value: /home/biocodz/sortmerna/run of previous flag: -workdir
   [opt_workdir:990] Using WORKDIR: "/home/biocodz/sortmerna/run" as specified
   [process:1456] Processing option: blast with value: 1 cigar qcov
   [process:1456] Processing option: fastx with value:
   [process:1456] Processing option: index with value: 2
   [opt_index:1157] using 'index' with specified value 2
   [process:1456] Processing option: other with value:
   [opt_other:267] other was specified without argument. Will use default Directory and Prefix for the non-aligned output.
   [process:1456] Processing option: out2 with value:
   [process:1456] Processing option: reads with value: /media/sf_a01_data/bio/reads/rna/SRR1635864_1_2M.fq.gz
   [opt_reads:97] Processing reads file [1] out of total [2] files
   [process:1456] Processing option: reads with value: /media/sf_a01_data/bio/reads/rna/SRR1635864_2_2M.fq.gz
   [opt_reads:97] Processing reads file [2] out of total [2] files
   [process:1456] Processing option: ref with value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta
   [opt_ref:157] Processing reference [1] out of total [8] references
   [opt_ref:205] File "/media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta" exists and is readable
   [process:1456] Processing option: ref with value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta
   [opt_ref:157] Processing reference [2] out of total [8] references
   [opt_ref:205] File "/media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta" exists and is readable
   [process:1456] Processing option: ref with value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta
   [opt_ref:157] Processing reference [3] out of total [8] references
   [opt_ref:205] File "/media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta" exists and is readable
   [process:1456] Processing option: ref with value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta
   [opt_ref:157] Processing reference [4] out of total [8] references
   [opt_ref:205] File "/media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta" exists and is readable
   [process:1456] Processing option: ref with value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta
   [opt_ref:157] Processing reference [5] out of total [8] references
   [opt_ref:205] File "/media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta" exists and is readable
   [process:1456] Processing option: ref with value: /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta
   [opt_ref:157] Processing reference [6] out of total [8] references
   [opt_ref:205] File "/media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta" exists and is readable
   [process:1456] Processing option: ref with value: /media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta
   [opt_ref:157] Processing reference [7] out of total [8] references
   [opt_ref:205] File "/media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta" exists and is readable
   [process:1456] Processing option: ref with value: /media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta
   [opt_ref:157] Processing reference [8] out of total [8] references
   [opt_ref:205] File "/media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta" exists and is readable
   [process:1456] Processing option: sout with value:
   [process:1456] Processing option: threads with value: 8
   [process:1456] Processing option: v with value:
   [process:1476] === Options processing done ===
   [process:1477] Alignment type: [best:1 num_alignments:1 min_lis:2 seeds:2]
   [validate_kvdbdir:1223] Key-value DB location "/home/biocodz/sortmerna/run/kvdb"
   [validate_kvdbdir:1259] Creating KVDB directory: "/home/biocodz/sortmerna/run/kvdb"
   [validate_idxdir:1189] Using index directory: "/home/biocodz/sortmerna/run/idx"
   [validate_idxdir:1205] IDX directory: "/home/biocodz/sortmerna/run/idx" exists and is not empty
   [validate_readb_dir:1281] Using split reads directory : "/home/biocodz/sortmerna/run/readb"
   [validate_readb_dir:1297] split reads directory : "/home/biocodz/sortmerna/run/readb" exists and is not empty
   [validate_aligned_pfx:1310] Checking output directory: "/home/biocodz/sortmerna/run/out"
   [main:62] Running command:
   /home/biocodz/sortmerna/dist/bin/sortmerna -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta -ref /media/sf_a01_code/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta -reads /media/sf_a01_data/bio/reads/rna/SRR1635864_1_2M.fq.gz -reads /media/sf_a01_data/bio/reads/rna/SRR1635864_2_2M.fq.gz -fastx -blast 1 cigar qcov -out2 -sout -other -v -threads 8 -index 2 -workdir /home/biocodz/sortmerna/run
   [Index:102] Found 32 non-empty index files. Skipping indexing.
   [init:108] Readfeed init started
   [define_format:885] file: "/media/sf_a01_data/bio/reads/rna/SRR1635864_1_2M.fq.gz" is FASTQ gzipped
   [define_format:885] file: "/media/sf_a01_data/bio/reads/rna/SRR1635864_2_2M.fq.gz" is FASTQ gzipped
   [count_reads:919] started count  ...
   [next:311] EOF FWD reached. Total reads: 500000
   [next:311] EOF REV reached. Total reads: 500000
   [count_reads:949] done count. Elapsed time: 5.73545 sec. Total reads: 1000000
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/fwd_0.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/rev_0.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/fwd_1.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/rev_1.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/fwd_2.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/rev_2.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/fwd_3.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/rev_3.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/fwd_4.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/rev_4.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/fwd_5.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/rev_5.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/fwd_6.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/rev_6.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/fwd_7.fq.gz
   [init_split_files:971] added file: /home/biocodz/sortmerna/run/readb/rev_7.fq.gz
   [is_split_ready:723] found existing readfeed descriptor /home/biocodz/sortmerna/run/readb/readfeed
   [split:583] start splitting. Using number of splits equals number of processing threads: 8
   [clean:1102] found descriptor /home/biocodz/sortmerna/run/readb/readfeed
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/fwd_0.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/rev_0.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/fwd_1.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/rev_1.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/fwd_2.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/rev_2.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/fwd_3.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/rev_3.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/fwd_4.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/rev_4.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/fwd_5.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/rev_5.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/fwd_6.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/rev_6.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/fwd_7.fq.gz
   [clean:1146] removing split file: /home/biocodz/sortmerna/run/readb/rev_7.fq.gz
   [next:311] EOF FWD reached. Total reads: 500000
   [next:311] EOF REV reached. Total reads: 500000
   [split:694] Done splitting. Reads count: 1000000 Runtime sec: 55.3599
   
   [init:135] Readfeed init done in sec [61.114]
   [store_to_db:292] Stored Reads statistics to DB:
        all_reads_count= 1000000 all_reads_len= 99050688 min_read_len= 101 max_read_len= 55 total_aligned= 0 total_aligned_id= 0 total_aligned_cov= 0 total_aligned_id_cov= 0 total_denovo= 0 num_short= 0 reads_matched_per_db= TODO is_stats_calc= 0 is_total_reads_mapped_cov= 0
   
   [align:143] ==== Starting alignment ====
   [align:146] Number of cores: 8
   [align:163] Using number of Processor threads: 8
   [Refstats:60] Index Statistics calculation starts ... done in: 4.26797 sec
   [align:185] Loading index: 0 part: 1/1 Memory KB: 9 ...
   [align:190] done in [0.937722] sec Memory KB: 407
   [align:193] Loading references ...
   [align:197] done in [0.70624] sec. Memory KB: 422
   [align2:70] Processor 0 thread 140010330240768 started
   [align2:70] Processor 1 thread 140010338633472 started
   [align2:70] Processor 7 thread 140010147415808 started
   [align2:70] Processor 4 thread 140010355418880 started
   [align2:70] Processor 5 thread 140010347026176 started
   [align2:70] Processor 6 thread 140010155808512 started
   [align2:70] Processor 2 thread 140010917435136 started
   [align2:70] Processor 3 thread 140010363811584 started
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 5 thread 140010347026176 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 551 Runtime sec: 21.6408
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 0 thread 140010330240768 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 659 Runtime sec: 21.6978
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 4 thread 140010355418880 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 656 Runtime sec: 21.6984
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 2 thread 140010917435136 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 650 Runtime sec: 21.7438
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 1 thread 140010338633472 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 559 Runtime sec: 21.9022
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 7 thread 140010147415808 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 649 Runtime sec: 21.9709
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 3 thread 140010363811584 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 730 Runtime sec: 22.0211
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 6 thread 140010155808512 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 667 Runtime sec: 22.0475
   [align:220] done index: 0 part: 1 in 22.0595 sec Memory KB: 427
   [align:227] Index and References unloaded in 0.197478 sec. Memory KB: 427
   [align:185] Loading index: 1 part: 1/1 Memory KB: 427 ...
   [align:190] done in [0.92918] sec Memory KB: 427
   [align:193] Loading references ...
   [align:197] done in [0.692128] sec. Memory KB: 427
   [align2:70] Processor 0 thread 140010147415808 started
   [align2:70] Processor 1 thread 140010155808512 started
   [align2:70] Processor 2 thread 140010347026176 started
   [align2:70] Processor 3 thread 140010355418880 started
   [align2:70] Processor 4 thread 140010917435136 started
   [align2:70] Processor 5 thread 140010363811584 started
   [align2:70] Processor 6 thread 140010338633472 started
   [align2:70] Processor 7 thread 140010330240768 started
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 1 thread 140010155808512 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 590 Runtime sec: 21.1139
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 5 thread 140010363811584 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 576 Runtime sec: 21.6393
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 7 thread 140010330240768 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 680 Runtime sec: 21.7497
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 2 thread 140010347026176 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 677 Runtime sec: 21.7859
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 3 thread 140010355418880 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 769 Runtime sec: 21.9126
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 4 thread 140010917435136 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 702 Runtime sec: 21.97
   
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 0 thread 140010147415808 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 701 Runtime sec: 22.0073
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 6 thread 140010338633472 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 709 Runtime sec: 23.2707
   [align:220] done index: 1 part: 1 in 23.2729 sec Memory KB: 429
   [align:227] Index and References unloaded in 0.170446 sec. Memory KB: 429
   [align:185] Loading index: 2 part: 1/1 Memory KB: 429 ...
   [align:190] done in [0.797775] sec Memory KB: 429
   [align:193] Loading references ...
   [align:197] done in [0.672418] sec. Memory KB: 429
   [align2:70] Processor 0 thread 140010330240768 started
   [align2:70] Processor 1 thread 140010338633472 started
   [align2:70] Processor 2 thread 140010363811584 started
   [align2:70] Processor 3 thread 140010917435136 started
   [align2:70] Processor 4 thread 140010355418880 started
   [align2:70] Processor 5 thread 140010347026176 started
   [align2:70] Processor 6 thread 140010155808512 started
   [align2:70] Processor 7 thread 140010147415808 started
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 5 thread 140010347026176 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 4344 Runtime sec: 23.7947
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 2 thread 140010363811584 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 4550 Runtime sec: 24.0882
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 1 thread 140010338633472 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 4313 Runtime sec: 24.5213
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 0 thread 140010330240768 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 4526 Runtime sec: 24.5237
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 6 thread 140010155808512 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 4659 Runtime sec: 24.6418
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 7 thread 140010147415808 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 4527 Runtime sec: 24.6815
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 4 thread 140010355418880 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 4716 Runtime sec: 24.9354
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 3 thread 140010917435136 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 4652 Runtime sec: 24.9785
   [align:220] done index: 2 part: 1 in 24.9854 sec Memory KB: 453
   [align:227] Index and References unloaded in 0.165393 sec. Memory KB: 453
   [align:185] Loading index: 3 part: 1/1 Memory KB: 453 ...
   [align:190] done in [1.09961] sec Memory KB: 528
   [align:193] Loading references ...
   [align:197] done in [0.959064] sec. Memory KB: 547
   [align2:70] Processor 0 thread 140010147415808 started
   [align2:70] Processor 1 thread 140010155808512 started
   [align2:70] Processor 2 thread 140010347026176 started
   [align2:70] Processor 3 thread 140010355418880 started
   [align2:70] Processor 4 thread 140010917435136 started
   [align2:70] Processor 5 thread 140010363811584 started
   [align2:70] Processor 6 thread 140010338633472 started
   [align2:70] Processor 7 thread 140010330240768 started
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 0 thread 140010147415808 done. Processed 124991 reads. Skipped already processed: 5 reads Aligned reads (passing E-value): 6149 Runtime sec: 26.0607
   [next:433] EOF REV reached. Total reads: 62500
   [next:433] EOF FWD reached. Total reads: 62500
   [align2:133] Processor 3 thread 140010355418880 done. Processed 124997 reads. Skipped already processed: 3 reads Aligned reads (passing E-value): 6341 Runtime sec: 27.4329
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 6 thread 140010338633472 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6333 Runtime sec: 27.4787
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 1 thread 140010155808512 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 5921 Runtime sec: 27.7387
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 4 thread 140010917435136 done. Processed 124993 reads. Skipped already processed: 6 reads Aligned reads (passing E-value): 6326 Runtime sec: 27.956
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 7 thread 140010330240768 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6202 Runtime sec: 28.2308
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 5 thread 140010363811584 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 5950 Runtime sec: 28.3352
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 2 thread 140010347026176 done. Processed 124995 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6285 Runtime sec: 28.5799
   [align:220] done index: 3 part: 1 in 28.5808 sec Memory KB: 573
   [align:227] Index and References unloaded in 0.221517 sec. Memory KB: 573
   [align:185] Loading index: 4 part: 1/1 Memory KB: 573 ...
   [align:190] done in [0.252196] sec Memory KB: 573
   [align:193] Loading references ...
   [align:197] done in [0.0478474] sec. Memory KB: 573
   [align2:70] Processor 0 thread 140010330240768 started
   [align2:70] Processor 1 thread 140010338633472 started
   [align2:70] Processor 7 thread 140010221201152 started
   [align2:70] Processor 6 thread 140010229593856 started
   [align2:70] Processor 2 thread 140010363811584 started
   [align2:70] Processor 3 thread 140010917435136 started
   [align2:70] Processor 4 thread 140010355418880 started
   [align2:70] Processor 5 thread 140010347026176 started
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 1 thread 140010338633472 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 5921 Runtime sec: 10.3744
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 2 thread 140010363811584 done. Processed 124995 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6285 Runtime sec: 10.457
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 5 thread 140010347026176 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 5950 Runtime sec: 10.4737
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 6 thread 140010229593856 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6333 Runtime sec: 10.5179
   [next:433] EOF REV reached. Total reads: 62500
   [next:433] EOF FWD reached. Total reads: 62500
   [align2:133] Processor 3 thread 140010917435136 done. Processed 124997 reads. Skipped already processed: 3 reads Aligned reads (passing E-value): 6341 Runtime sec: 10.5141
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 4 thread 140010355418880 done. Processed 124993 reads. Skipped already processed: 6 reads Aligned reads (passing E-value): 6326 Runtime sec: 10.5275
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 7 thread 140010221201152 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6202 Runtime sec: 10.5978
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 0 thread 140010330240768 done. Processed 124991 reads. Skipped already processed: 5 reads Aligned reads (passing E-value): 6149 Runtime sec: 10.8378
   [align:220] done index: 4 part: 1 in 10.9213 sec Memory KB: 573
   [align:227] Index and References unloaded in 0.109545 sec. Memory KB: 573
   [align:185] Loading index: 5 part: 1/1 Memory KB: 573 ...
   [align:190] done in [0.469406] sec Memory KB: 573
   [align:193] Loading references ...
   [align:197] done in [0.244222] sec. Memory KB: 573
   [align2:70] Processor 0 thread 140010221201152 started
   [align2:70] Processor 1 thread 140010229593856 started
   [align2:70] Processor 7 thread 140010330240768 started
   [align2:70] Processor 2 thread 140010347026176 started
   [align2:70] Processor 3 thread 140010355418880 started
   [align2:70] Processor 4 thread 140010917435136 started
   [align2:70] Processor 5 thread 140010363811584 started
   [align2:70] Processor 6 thread 140010338633472 started
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 6 thread 140010338633472 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6334 Runtime sec: 11.6236
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 7 thread 140010330240768 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6202 Runtime sec: 11.9621
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 4 thread 140010917435136 done. Processed 124993 reads. Skipped already processed: 6 reads Aligned reads (passing E-value): 6327 Runtime sec: 12.446
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 5 thread 140010363811584 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 5950 Runtime sec: 13.5967
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 1 thread 140010229593856 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 5922 Runtime sec: 13.7651
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 2 thread 140010347026176 done. Processed 124995 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6285 Runtime sec: 13.8438
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 0 thread 140010221201152 done. Processed 124991 reads. Skipped already processed: 5 reads Aligned reads (passing E-value): 6151 Runtime sec: 13.992
   [next:433] EOF REV reached. Total reads: 62500
   [next:433] EOF FWD reached. Total reads: 62500
   [align2:133] Processor 3 thread 140010355418880 done. Processed 124997 reads. Skipped already processed: 3 reads Aligned reads (passing E-value): 6344 Runtime sec: 14.0103
   [align:220] done index: 5 part: 1 in 14.0488 sec Memory KB: 574
   [align:227] Index and References unloaded in 0.0635204 sec. Memory KB: 574
   [align:185] Loading index: 6 part: 1/1 Memory KB: 574 ...
   [align:190] done in [0.545306] sec Memory KB: 574
   [align:193] Loading references ...
   [align:197] done in [0.473546] sec. Memory KB: 567
   [align2:70] Processor 0 thread 140010330240768 started
   [align2:70] Processor 1 thread 140010338633472 started
   [align2:70] Processor 2 thread 140010363811584 started
   [align2:70] Processor 3 thread 140010917435136 started
   [align2:70] Processor 4 thread 140010355418880 started
   [align2:70] Processor 5 thread 140010347026176 started
   [align2:70] Processor 6 thread 140010229593856 started
   [align2:70] Processor 7 thread 140010221201152 started
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 7 thread 140010221201152 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6250 Runtime sec: 15.0544
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 1 thread 140010338633472 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 5967 Runtime sec: 15.0637
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 6 thread 140010229593856 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6378 Runtime sec: 15.0709
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 4 thread 140010355418880 done. Processed 124993 reads. Skipped already processed: 6 reads Aligned reads (passing E-value): 6369 Runtime sec: 15.088
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 0 thread 140010330240768 done. Processed 124991 reads. Skipped already processed: 5 reads Aligned reads (passing E-value): 6189 Runtime sec: 15.1061
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 2 thread 140010363811584 done. Processed 124995 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6344 Runtime sec: 15.1101
   [next:433] EOF REV reached. Total reads: 62500
   [next:433] EOF FWD reached. Total reads: 62500
   [align2:133] Processor 3 thread 140010917435136 done. Processed 124997 reads. Skipped already processed: 3 reads Aligned reads (passing E-value): 6382 Runtime sec: 15.1477
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 5 thread 140010347026176 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 6000 Runtime sec: 15.4381
   [align:220] done index: 6 part: 1 in 15.4438 sec Memory KB: 567
   [align:227] Index and References unloaded in 0.0704959 sec. Memory KB: 555
   [align:185] Loading index: 7 part: 1/1 Memory KB: 555 ...
   [align:190] done in [0.228702] sec Memory KB: 555
   [align:193] Loading references ...
   [align:197] done in [0.147236] sec. Memory KB: 555
   [align2:70] Processor 0 thread 140010221201152 started
   [align2:70] Processor 1 thread 140010229593856 started
   [align2:70] Processor 7 thread 140010330240768 started
   [align2:70] Processor 2 thread 140010347026176 started
   [align2:70] Processor 3 thread 140010355418880 started
   [align2:70] Processor 6 thread 140010338633472 started
   [align2:70] Processor 4 thread 140010917435136 started
   [align2:70] Processor 5 thread 140010363811584 started
   [next:433] EOF REV reached. Total reads: 62500
   [next:433] EOF FWD reached. Total reads: 62500
   [align2:133] Processor 3 thread 140010355418880 done. Processed 124997 reads. Skipped already processed: 3 reads Aligned reads (passing E-value): 6382 Runtime sec: 9.78541
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 1 thread 140010229593856 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 5967 Runtime sec: 9.83089
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 6 thread 140010338633472 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6378 Runtime sec: 9.8955
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 7 thread 140010330240768 done. Processed 124996 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6250 Runtime sec: 9.92369
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 4 thread 140010917435136 done. Processed 124993 reads. Skipped already processed: 6 reads Aligned reads (passing E-value): 6369 Runtime sec: 9.94304
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 5 thread 140010363811584 done. Processed 125000 reads. Skipped already processed: 0 reads Aligned reads (passing E-value): 6000 Runtime sec: 10.0105
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 2 thread 140010347026176 done. Processed 124995 reads. Skipped already processed: 4 reads Aligned reads (passing E-value): 6344 Runtime sec: 10.1247
   [next:433] EOF REV reached. Total reads: 62500
   [align2:133] Processor 0 thread 140010221201152 done. Processed 124991 reads. Skipped already processed: 5 reads Aligned reads (passing E-value): 6189 Runtime sec: 10.675
   [align:220] done index: 7 part: 1 in 10.6762 sec Memory KB: 556
   [align:227] Index and References unloaded in 0.0470631 sec. Memory KB: 556
   [align:237] ==== Done alignment in 160.304 sec ====
   
   [store_to_db:292] Stored Reads statistics to DB:
        all_reads_count= 1000000 all_reads_len= 99050688 min_read_len= 101 max_read_len= 55 total_aligned= 49909 total_aligned_id= 0 total_aligned_cov= 0 total_aligned_id_cov= 0 total_denovo= 0 num_short= 0 reads_matched_per_db= TODO is_stats_calc= 0 is_total_reads_mapped_cov= 0
   
   [writeSummary:179] ==== Starting summary of alignment statistics ====
   [Refstats:60] Index Statistics calculation starts ... done in: 4.45322 sec
   [write:62] Using summary file: /home/biocodz/sortmerna/run/out/aligned.log
   [writeSummary:185] ==== Done summary in sec [4.45347] ====
   
   [writeReports:160] === Report generation starts ===
   [writeReports:175] Restored Readstats from DB: 1
   [Refstats:60] Index Statistics calculation starts ... done in: 4.43634 sec
   [validate_out_type:139] Output type:
   1-file  2-files  paired  paired_in  paired_out  out2  sout  other  otype
              +        +                              +     +     +     66
   [set_num_out:162] num_out: 4
   [init:60] num_out: 4
   [validate_out_type:139] Output type:
   1-file  2-files  paired  paired_in  paired_out  out2  sout  other  otype
              +        +                              +     +     +     66
   [set_num_out:162] num_out: 4
   [init:60] num_out: 4
   [writeReports:190] loading reference 0 part 1/1 ... done in 0.791885 sec
   [report:93] Report Processor: 0 thread: 140010330240768 started. Memory KB: 556
   [report:93] Report Processor: 1 thread: 140010338633472 started. Memory KB: 556
   [report:93] Report Processor: 2 thread: 140010363811584 started. Memory KB: 556
   [report:93] Report Processor: 3 thread: 140010917435136 started. Memory KB: 556
   [report:93] Report Processor: 4 thread: 140010355418880 started. Memory KB: 556
   [report:93] Report Processor: 5 thread: 140010347026176 started. Memory KB: 556
   [report:93] Report Processor: 6 thread: 140010229593856 started. Memory KB: 556
   [report:93] Report Processor: 7 thread: 140010221201152 started. Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 1 thread: 140010338633472 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 4 thread: 140010355418880 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 5 thread: 140010347026176 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 2 thread: 140010363811584 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 0 thread: 140010330240768 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 6 thread: 140010229593856 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 7 thread: 140010221201152 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 3 thread: 140010917435136 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [writeReports:211] done reference 0 part: 1 in 24.8879 sec
   [writeReports:217] references unloaded in 0.00277428 sec Memory KB: 556
   [writeReports:190] loading reference 1 part 1/1 ... done in 0.827815 sec
   [report:93] Report Processor: 0 thread: 140010221201152 started. Memory KB: 556
   [report:93] Report Processor: 1 thread: 140010229593856 started. Memory KB: 556
   [report:93] Report Processor: 2 thread: 140010347026176 started. Memory KB: 556
   [report:93] Report Processor: 3 thread: 140010355418880 started. Memory KB: 556
   [report:93] Report Processor: 4 thread: 140010917435136 started. Memory KB: 556
   [report:93] Report Processor: 5 thread: 140010363811584 started. Memory KB: 556
   [report:93] Report Processor: 6 thread: 140010338633472 started. Memory KB: 556
   [report:93] Report Processor: 7 thread: 140010330240768 started. Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 1 thread: 140010229593856 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 2 thread: 140010347026176 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 0 thread: 140010221201152 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 6 thread: 140010338633472 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 5 thread: 140010363811584 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 7 thread: 140010330240768 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 3 thread: 140010355418880 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 4 thread: 140010917435136 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [writeReports:211] done reference 1 part: 1 in 4.43556 sec
   [writeReports:217] references unloaded in 0.00518327 sec Memory KB: 556
   [writeReports:190] loading reference 2 part 1/1 ... done in 0.719366 sec
   [report:93] Report Processor: 0 thread: 140010330240768 started. Memory KB: 556
   [report:93] Report Processor: 1 thread: 140010338633472 started. Memory KB: 556
   [report:93] Report Processor: 7 thread: 140010221201152 started. Memory KB: 556
   [report:93] Report Processor: 6 thread: 140010229593856 started. Memory KB: 556
   [report:93] Report Processor: 5 thread: 140010347026176 started. Memory KB: 556
   [report:93] Report Processor: 2 thread: 140010363811584 started. Memory KB: 556
   [report:93] Report Processor: 4 thread: 140010355418880 started. Memory KB: 556
   [report:93] Report Processor: 3 thread: 140010917435136 started. Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 1 thread: 140010338633472 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 4 thread: 140010355418880 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 5 thread: 140010347026176 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 2 thread: 140010363811584 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 3 thread: 140010917435136 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 6 thread: 140010229593856 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 7 thread: 140010221201152 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 0 thread: 140010330240768 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [writeReports:211] done reference 2 part: 1 in 4.91457 sec
   [writeReports:217] references unloaded in 0.00132068 sec Memory KB: 556
   [writeReports:190] loading reference 3 part 1/1 ... done in 1.05361 sec
   [report:93] Report Processor: 0 thread: 140010221201152 started. Memory KB: 556
   [report:93] Report Processor: 1 thread: 140010229593856 started. Memory KB: 556
   [report:93] Report Processor: 2 thread: 140010347026176 started. Memory KB: 556
   [report:93] Report Processor: 3 thread: 140010355418880 started. Memory KB: 556
   [report:93] Report Processor: 4 thread: 140010917435136 started. Memory KB: 556
   [report:93] Report Processor: 5 thread: 140010363811584 started. Memory KB: 556
   [report:93] Report Processor: 6 thread: 140010338633472 started. Memory KB: 556
   [report:93] Report Processor: 7 thread: 140010330240768 started. Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 1 thread: 140010229593856 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 2 thread: 140010347026176 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 5 thread: 140010363811584 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 4 thread: 140010917435136 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 3 thread: 140010355418880 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 0 thread: 140010221201152 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 6 thread: 140010338633472 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 7 thread: 140010330240768 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [writeReports:211] done reference 3 part: 1 in 4.36404 sec
   [writeReports:217] references unloaded in 0.00403434 sec Memory KB: 556
   [writeReports:190] loading reference 4 part 1/1 ... done in 0.0597788 sec
   [report:93] Report Processor: 0 thread: 140010330240768 started. Memory KB: 556
   [report:93] Report Processor: 1 thread: 140010338633472 started. Memory KB: 556
   [report:93] Report Processor: 2 thread: 140010363811584 started. Memory KB: 556
   [report:93] Report Processor: 6 thread: 140010229593856 started. Memory KB: 556
   [report:93] Report Processor: 7 thread: 140010221201152 started. Memory KB: 556
   [report:93] Report Processor: 3 thread: 140010917435136 started. Memory KB: 556
   [report:93] Report Processor: 5 thread: 140010347026176 started. Memory KB: 556
   [report:93] Report Processor: 4 thread: 140010355418880 started. Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 0 thread: 140010330240768 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 7 thread: 140010221201152 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 6 thread: 140010229593856 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 1 thread: 140010338633472 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 2 thread: 140010363811584 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 3 thread: 140010917435136 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 4 thread: 140010355418880 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 5 thread: 140010347026176 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [writeReports:211] done reference 4 part: 1 in 4.04392 sec
   [writeReports:217] references unloaded in 4.0988e-05 sec Memory KB: 556
   [writeReports:190] loading reference 5 part 1/1 ... done in 0.229265 sec
   [report:93] Report Processor: 0 thread: 140010221201152 started. Memory KB: 556
   [report:93] Report Processor: 6 thread: 140010338633472 started. Memory KB: 556
   [report:93] Report Processor: 3 thread: 140010355418880 started. Memory KB: 556
   [report:93] Report Processor: 5 thread: 140010363811584 started. Memory KB: 556
   [report:93] Report Processor: 7 thread: 140010330240768 started. Memory KB: 556
   [report:93] Report Processor: 4 thread: 140010917435136 started. Memory KB: 556
   [report:93] Report Processor: 1 thread: 140010229593856 started. Memory KB: 556
   [report:93] Report Processor: 2 thread: 140010347026176 started. Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 5 thread: 140010363811584 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 3 thread: 140010355418880 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 4 thread: 140010917435136 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 6 thread: 140010338633472 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 0 thread: 140010221201152 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 7 thread: 140010330240768 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 2 thread: 140010347026176 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 1 thread: 140010229593856 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [writeReports:211] done reference 5 part: 1 in 3.89069 sec
   [writeReports:217] references unloaded in 0.00141204 sec Memory KB: 556
   [writeReports:190] loading reference 6 part 1/1 ... done in 0.494027 sec
   [report:93] Report Processor: 0 thread: 140010330240768 started. Memory KB: 556
   [report:93] Report Processor: 4 thread: 140010355418880 started. Memory KB: 556
   [report:93] Report Processor: 1 thread: 140010338633472 started. Memory KB: 556
   [report:93] Report Processor: 5 thread: 140010347026176 started. Memory KB: 556
   [report:93] Report Processor: 6 thread: 140010229593856 started. Memory KB: 556
   [report:93] Report Processor: 2 thread: 140010363811584 started. Memory KB: 556
   [report:93] Report Processor: 7 thread: 140010221201152 started. Memory KB: 556
   [report:93] Report Processor: 3 thread: 140010917435136 started. Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 0 thread: 140010330240768 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 1 thread: 140010338633472 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 3 thread: 140010917435136 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 7 thread: 140010221201152 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 5 thread: 140010347026176 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 6 thread: 140010229593856 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 2 thread: 140010363811584 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 4 thread: 140010355418880 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [writeReports:211] done reference 6 part: 1 in 4.31935 sec
   [writeReports:217] references unloaded in 0.00583015 sec Memory KB: 556
   [writeReports:190] loading reference 7 part 1/1 ... done in 0.151042 sec
   [report:93] Report Processor: 0 thread: 140010221201152 started. Memory KB: 556
   [report:93] Report Processor: 1 thread: 140010229593856 started. Memory KB: 556
   [report:93] Report Processor: 2 thread: 140010347026176 started. Memory KB: 556
   [report:93] Report Processor: 3 thread: 140010355418880 started. Memory KB: 556
   [report:93] Report Processor: 7 thread: 140010330240768 started. Memory KB: 556
   [report:93] Report Processor: 4 thread: 140010917435136 started. Memory KB: 556
   [report:93] Report Processor: 6 thread: 140010338633472 started. Memory KB: 556
   [report:93] Report Processor: 5 thread: 140010363811584 started. Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 4 thread: 140010917435136 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 7 thread: 140010330240768 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 6 thread: 140010338633472 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 0 thread: 140010221201152 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 2 thread: 140010347026176 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 5 thread: 140010363811584 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 3 thread: 140010355418880 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [next:433] EOF FWD reached. Total reads: 62500
   [next:433] EOF REV reached. Total reads: 62500
   [report:152] Report processor: 1 thread: 140010229593856 done. Processed reads: 125000 Invalid reads: 0 Memory KB: 556
   [writeReports:211] done reference 7 part: 1 in 4.29868 sec
   [writeReports:217] references unloaded in 0.0028677 sec Memory KB: 556
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_fwd_1.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_fwd_1.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_fwd_1.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_fwd_2.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_fwd_2.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_fwd_2.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_fwd_3.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_fwd_3.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_fwd_3.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_fwd_4.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_fwd_4.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_fwd_4.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_fwd_5.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_fwd_5.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_fwd_5.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_fwd_6.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_fwd_6.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_fwd_6.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_fwd_7.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_fwd_7.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_fwd_7.fq.gz
   [strip_path_sfx:154] moving /home/biocodz/sortmerna/run/out/aligned_paired_fwd_0.fq.gz -> "/home/biocodz/sortmerna/run/out/aligned_paired_fwd.fq.gz"
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_rev_1.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_rev_1.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_rev_1.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_rev_2.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_rev_2.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_rev_2.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_rev_3.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_rev_3.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_rev_3.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_rev_4.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_rev_4.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_rev_4.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_rev_5.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_rev_5.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_rev_5.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_rev_6.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_rev_6.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_rev_6.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_paired_rev_7.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_paired_rev_7.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_paired_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_paired_rev_7.fq.gz
   [strip_path_sfx:154] moving /home/biocodz/sortmerna/run/out/aligned_paired_rev_0.fq.gz -> "/home/biocodz/sortmerna/run/out/aligned_paired_rev.fq.gz"
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_1.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_1.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_1.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_2.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_2.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_2.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_3.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_3.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_3.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_4.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_4.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_4.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_5.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_5.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_5.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_6.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_6.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_6.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_7.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_7.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_7.fq.gz
   [strip_path_sfx:154] moving /home/biocodz/sortmerna/run/out/aligned_singleton_fwd_0.fq.gz -> "/home/biocodz/sortmerna/run/out/aligned_singleton_fwd.fq.gz"
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_rev_1.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_rev_1.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_rev_1.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_rev_2.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_rev_2.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_rev_2.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_rev_3.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_rev_3.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_rev_3.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_rev_4.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_rev_4.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_rev_4.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_rev_5.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_rev_5.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_rev_5.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_rev_6.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_rev_6.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_rev_6.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/aligned_singleton_rev_7.fq.gz for reading.
   [merge:155] merged /home/biocodz/sortmerna/run/out/aligned_singleton_rev_7.fq.gz -> /home/biocodz/sortmerna/run/out/aligned_singleton_rev_0.fq.gz
   [merge:158] deleted /home/biocodz/sortmerna/run/out/aligned_singleton_rev_7.fq.gz
   [strip_path_sfx:154] moving /home/biocodz/sortmerna/run/out/aligned_singleton_rev_0.fq.gz -> "/home/biocodz/sortmerna/run/out/aligned_singleton_rev.fq.gz"
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_fwd_1.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_fwd_1.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_fwd_1.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_fwd_2.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_fwd_2.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_fwd_2.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_fwd_3.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_fwd_3.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_fwd_3.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_fwd_4.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_fwd_4.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_fwd_4.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_fwd_5.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_fwd_5.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_fwd_5.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_fwd_6.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_fwd_6.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_fwd_6.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_fwd_7.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_fwd_7.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_fwd_7.fq.gz
   [strip_path_sfx:154] moving /home/biocodz/sortmerna/run/out/other_paired_fwd_0.fq.gz -> "/home/biocodz/sortmerna/run/out/other_paired_fwd.fq.gz"
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_rev_1.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_rev_1.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_rev_1.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_rev_2.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_rev_2.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_rev_2.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_rev_3.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_rev_3.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_rev_3.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_rev_4.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_rev_4.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_rev_4.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_rev_5.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_rev_5.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_rev_5.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_rev_6.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_rev_6.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_rev_6.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_paired_rev_7.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_paired_rev_7.fq.gz -> /home/biocodz/sortmerna/run/out/other_paired_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_paired_rev_7.fq.gz
   [strip_path_sfx:154] moving /home/biocodz/sortmerna/run/out/other_paired_rev_0.fq.gz -> "/home/biocodz/sortmerna/run/out/other_paired_rev.fq.gz"
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_fwd_1.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_fwd_1.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_fwd_1.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_fwd_2.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_fwd_2.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_fwd_2.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_fwd_3.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_fwd_3.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_fwd_3.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_fwd_4.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_fwd_4.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_fwd_4.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_fwd_5.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_fwd_5.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_fwd_5.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_fwd_6.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_fwd_6.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_fwd_6.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_fwd_7.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_fwd_7.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_fwd_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_fwd_7.fq.gz
   [strip_path_sfx:154] moving /home/biocodz/sortmerna/run/out/other_singleton_fwd_0.fq.gz -> "/home/biocodz/sortmerna/run/out/other_singleton_fwd.fq.gz"
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_rev_1.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_rev_1.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_rev_1.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_rev_2.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_rev_2.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_rev_2.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_rev_3.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_rev_3.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_rev_3.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_rev_4.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_rev_4.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_rev_4.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_rev_5.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_rev_5.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_rev_5.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_rev_6.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_rev_6.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_rev_6.fq.gz
   [openfr:114] Opened output file /home/biocodz/sortmerna/run/out/other_singleton_rev_7.fq.gz for reading.
   [merge:130] merged /home/biocodz/sortmerna/run/out/other_singleton_rev_7.fq.gz -> /home/biocodz/sortmerna/run/out/other_singleton_rev_0.fq.gz
   [merge:133] deleted /home/biocodz/sortmerna/run/out/other_singleton_rev_7.fq.gz
   [strip_path_sfx:154] moving /home/biocodz/sortmerna/run/out/other_singleton_rev_0.fq.gz -> "/home/biocodz/sortmerna/run/out/other_singleton_rev.fq.gz"
   [merge:67] merged /home/biocodz/sortmerna/run/out/aligned_1.blast.gz -> /home/biocodz/sortmerna/run/out/aligned_0.blast.gz
   [merge:70] deleted /home/biocodz/sortmerna/run/out/aligned_1.blast.gz
   [merge:67] merged /home/biocodz/sortmerna/run/out/aligned_2.blast.gz -> /home/biocodz/sortmerna/run/out/aligned_0.blast.gz
   [merge:70] deleted /home/biocodz/sortmerna/run/out/aligned_2.blast.gz
   [merge:67] merged /home/biocodz/sortmerna/run/out/aligned_3.blast.gz -> /home/biocodz/sortmerna/run/out/aligned_0.blast.gz
   [merge:70] deleted /home/biocodz/sortmerna/run/out/aligned_3.blast.gz
   [merge:67] merged /home/biocodz/sortmerna/run/out/aligned_4.blast.gz -> /home/biocodz/sortmerna/run/out/aligned_0.blast.gz
   [merge:70] deleted /home/biocodz/sortmerna/run/out/aligned_4.blast.gz
   [merge:67] merged /home/biocodz/sortmerna/run/out/aligned_5.blast.gz -> /home/biocodz/sortmerna/run/out/aligned_0.blast.gz
   [merge:70] deleted /home/biocodz/sortmerna/run/out/aligned_5.blast.gz
   [merge:67] merged /home/biocodz/sortmerna/run/out/aligned_6.blast.gz -> /home/biocodz/sortmerna/run/out/aligned_0.blast.gz
   [merge:70] deleted /home/biocodz/sortmerna/run/out/aligned_6.blast.gz
   [merge:67] merged /home/biocodz/sortmerna/run/out/aligned_7.blast.gz -> /home/biocodz/sortmerna/run/out/aligned_0.blast.gz
   [merge:70] deleted /home/biocodz/sortmerna/run/out/aligned_7.blast.gz
   [strip_path_sfx:154] moving /home/biocodz/sortmerna/run/out/aligned_0.blast.gz -> "/home/biocodz/sortmerna/run/out/aligned.blast.gz"
   [writeReports:259] === done Reports in 64.4271 sec ===
   
   [run] Run time: 294.6793098449707
   Testing num_reads: 1000000 Expected: 1000000
   Testing num_hits: 49909 Expected: 49909
   Testing num_fail: 950091 Expected: 950091
   processing Blast file: /home/biocodz/sortmerna/run/out/aligned.blast.gz
   [process_blast] TODO: implement gz processing

Choosing parameters for filtering and read mapping
==================================================
   
Users have the option to output sequence alignments for their matching rRNA reads in the SAM or BLAST-like formats. Depending on the desired quality of alignments, different parameters must be set. Table 1 presents a guide to setting parameters for most use cases. In all cases, output alignments are always guaranteed to reach the threshold E-value score (default E-value=1). An E-value of 1 signifies that one random alignment is expected for aligning all reads against the reference database. The E-value is computed for the entire search space, not per read. 
  
+------------------+--------------------+-----------------------------------------------------------------------------------------------+
| Option           | Speed              | Description                                                                                   |
+==================+====================+===============================================================================================+
|                  | Very fast for INT=1| Output the first alignment passing E-value threshold (best choice if only filtering is needed)|
|                  |                    |                                                                                               |
|                  +--------------------+-----------------------------------------------------------------------------------------------+
| --num-alignment  | Speed decreases    |                                                                                               |
| INT              | for higher value   | Higher INT signifies more alignments will be made & output                                    |
|                  | INT                |                                                                                               |
|                  +--------------------+-----------------------------------------------------------------------------------------------+
|                  | Very slow for INT=0| All alignments reaching the E-value threshold are reported (this option is not suggested for  |
|                  |                    | high similarity rRNA databases, due to many possible alignments per read causing a very       |
|                  |                    | large file output)                                                                            |
+------------------+--------------------+-----------------------------------------------------------------------------------------------+
|                  | Fast for INT=1     | Only one high-candidate reference sequence will be searched for alignments (determined        |
|                  |                    | heuristically using a Longest Increasing Sub-sequence of seed matches). The single best       |
|                  |                    | alignment of those will be reported                                                           | 
|                  |                    |                                                                                               |
|                  +--------------------+-----------------------------------------------------------------------------------------------+
| --best INT       | Speed decreases    |                                                                                               |
|                  | for higher value   | Higher INT signies more alignments will be made, though only the best one will be reported    |
|                  | INT                |                                                                                               |
|                  +--------------------+-----------------------------------------------------------------------------------------------+
|                  | Very slow for INT=0| All high-candidate reference sequences will be searched for alignments, though only the best  | 
|                  |                    | one will be reported                                                                          |
|                  |                    |                                                                                               |
|                  |                    |                                                                                               |
+------------------+--------------------+-----------------------------------------------------------------------------------------------+

Contents
--------

.. toctree::

   installation