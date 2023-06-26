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
