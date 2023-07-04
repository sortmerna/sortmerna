Execution trace
===============

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
