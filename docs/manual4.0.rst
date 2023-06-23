Introduction
============

SortMeRNA is a local sequence alignment tool for filtering, mapping and OTU clustering. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of NGS reads. The main application of SortMeRNA is filtering rRNA from metatranscriptomic data. Additional applications include OTU-picking and taxonomy assignation available through QIIME v1.9+ http://qiime.org. SortMeRNA takes as input one or two (paired) reads file(s) (fasta/fasta.gz/fastq/fastq.gz), and one or multiple rRNA database file(s), and sorts apart aligned and rejected reads into two files. SortMeRNA works with Illumina, 454, Ion Torrent and PacBio data, and can produce SAM and BLAST-like alignments.

.. image:: img/smr-3-processing_pipeline.svg

::

    Figure 1. Processing smr-3-processing_pipeline

Databases
---------

+-----------------------+---+----------------+----------------------+---------------+
|representative database|%id|#seq (clustered)|origin                |#seq (original)|
+=======================+===+================+======================+===============+
|silva-bac-16s-id90     |90 |12798           |SILVA SSU Ref NR v.119|464618         | 
+-----------------------+---|----------------+----------------------+---------------+
|silva-arc-16s-id95     |95 |3193            |SILVA SSU Ref NR v.119|18797          |
+-----------------------+---|----------------+----------------------+---------------+
|silva-euk-18s-id95     |95 |7348            |SILVA SSU Ref NR v.119|51553          |
+-----------------------+---|----------------+----------------------+---------------+
|silva-bac-23s-id98     |98 |4488            |SILVA LSU Ref v.119   |43822          |
+-----------------------+---|----------------+----------------------+---------------+
|silva-arc-23s-id98     |98 |251             |SILVA LSU Ref v.119   |629            |
+-----------------------+---|----------------+----------------------+---------------+
|silva-euk-28s-id98     |98 |4935            |SILVA LSU Ref v.119   |13095          |
+-----------------------+---|----------------+----------------------+---------------+
|rfam-5s-id98 98        |98 |59513           |RFAM                  |116760         |
+-----------------------+---|----------------+----------------------+---------------+
|rfam-5.8s-id98         |98 |13034           |RFAM                  |225185         |
+-----------------------+---|----------------+----------------------+---------------+

HMMER 3.1b1 and SumaClust v1.0.00 were used to reduce the size of the original databases to the similarity listed in column 2 (%id) of the table above (see /sortmerna/rRNA databases/README.txt for a list of complete steps).

These representative databases were specifically made for fast fltering of rRNA. Approximately the same number of rRNA will be fltered using silva-bac-16s-id90 (12802 rRNA) as using Greengenes 97% (99322 rRNA), but the former will run signifcantly faster.

id %: members of the cluster must have identity at least this % id with the representative sequence.

Remark: The user must first index the fasta database by using the command indexdb and then filter/map reads against the database using the command sortmerna. 

Usage 
-----

Highlights:
- 4.0 has a single executable :code:`sortmerna` (no more :code:`indexdb`). Indexing is controlled using the :code:`sortmerna indexing options.
- The only required options are :code:`-ref` and :code:`-reads`
- Multiple :code:`ref` accepted
- Two :code:`-reads` options accepted for paired reads files
- Options (any) can be specified using a double or single dasg e.g :code:`-ref`/:code:`--ref` and :code:`-reads`/:code:`--reads`
- Both plain :code:`fasta`/:code:`fastq` anc archived :code:`fasta.gz`/:code:`fastq.gz` files are accepted
- Relative paths are accepted
- Work Directory:
    - Sortmerna uses a :code:`work directory` to store
        - Reference index Database
        - Key-value Database for alignment results
        - Output

Overview
--------

Work directory
##############

SortMeRNA uses a Work directory with the default location :code:`$USER/sortmerna/run/`. The work directory has the following structure::

    $USER/sortmerna/run/
    kvdb/       key-value datastore for alignment results
    idx/        index database
    out/        output files like aligned.blast

    On each new run SortMeRNA checks the :code:`idx` directory for already existing index files. If index files exist, the indexing won't be repeated, and the execution trace will show::

    [Index:104] Found 32 non-empty index files. Skipping indexing.

    if for whatever reason you'd like to run the indexing from scratch, just pourge the :code:`idx`::

    rm -rf $HOME/sortmerna/run/idx/

Sortmerna also checks :code:`kvdb` directory. This directory has to be Empty prior each new run. Unlike the index database, the alignment results are not reused, unless only report generation task is performed using the :code:`-task` option i.e. :code:`-task 1` or :code:`-task 2`. 

Purge :code:`kvdb` directory prior to each new run::

    rm -rf $HOME/sortmerna/run/kvdb/

For reference, a sample trace of a successful execution is provided. The trace is from a standard Sortmerna test found in the distribution.

Options
-------

+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|Option             |argument|description                                                                                                              |default              |
+===================+========+=========================================================================================================================+=====================+
|-ref               |PATH    |Reference file (FASTA) absolute or relative path. Use mutliple times, once per a reference file                          |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|-reads             |PATH    |Raw reads file (FASTA/FASTQ/FASTA.GZ/FASTQ.GZ). Use twice for files with paired reads                                    |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-workdir]         |PATH    |Working directory for storing the Reference index, Key-value database, Output                                            |USRDIR/sortmerna/run/| 
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-fastx]           |Boolean |Output aligned reads into FASTA/FASTQ file                                                                               |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-SQ]              |Boolean |Add SQ tags to the SAM file                                                                                              |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-blast]           |String  |output alignments in various Blast-like formats::                                                                        |                     |
|                   |        |   '0'                    - pairwise                                                                                     |                     |
|                   |        |   '1'                    - tabular(Blast - m 8 format)                                                                  |                     |
|                   |        |   '1 cigar'              - tabular + column for CIGAR                                                                   |                     |
|                   |        |   '1 cigar qcov'         - tabular + columns for CIGAR                                                                  |                     |
|                   |        |                            and query coverage                                                                           |                     |
|                   |        |   '1 cigar qcov qstrand' - tabular + columns for CIGAR,                                                                 |                     |
|                   |        |                            query coverage and strand                                                                    |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-other]           |Boolean |Create Non-aligned reads output file. Must be used with :code:`fastx`                                                    |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-num_alignments]  |        | Positive integer (>=0). Report first INT alignments per read reaching E-value. If Int = 0, all alignments will be output|False                |  
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-best]            |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-min_lis]         |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-print_all_reads] |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-paired_in]       |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-paired_out]      |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-match]           |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-mismatch]        |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-gap_open]        |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-gap_ext]         |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-a]               |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-e]               |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-F]               |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-N]               |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-R]               |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-id]              |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+
|[-coverage]        |        |                                                                                                                         |                     |
+-------------------+--------+-------------------------------------------------------------------------------------------------------------------------+---------------------+


List the version::

    sortmerna --version

List the hekp::

    sortmerna -h

Choosing parameters for filtering and read mapping
--------------------------------------------------
   
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

Table 1: Alignment parameter guide

Example: multiple databases and the fastes alignment option
##########################################################

::

    sortmerna -ref $SMR_HOME/data/silva-bac-16s-database-id85.fasta \
    -ref $SMR_HOME/data/silva-arc-16s-database-id95.fasta \
    -reads $SMR_HOME/data/set2_environmental_study_550_amplicon.fasta \
    -sam -fastx -blast 1 -num_alignments 1 -v

Filtering paired-end reads
##########################

The reads files are considred :code:`paired` if two reads files are provided (using :code:`reads` twice) e.g::

    sortmerna -ref Ref_1 -reads File_1 -reads File_2 ...    

The following alignment situations are possible for paired reads: 

+---+-----------------------------------------------------------------------------------------+-------------------------------------------------------------+
|A_0|Neither of the two paired reads is aligned                                               |No output is generated                                       |
+---+-----------------------------------------------------------------------------------------+-------------------------------------------------------------+
|A_1|A single out of a pair of reads is aligned (alignment score is above the given threshold)|see Table 3                                                  |
+---+-----------------------------------------------------------------------------------------+-------------------------------------------------------------+
|A_2|Both paired reads are aligned                                                            |both paired reads written into the :code:`aligned.fasta` file|
+---+-----------------------------------------------------------------------------------------+-------------------------------------------------------------+

::

    Table 2

In the A_1 case the output can be controlled using the following options:
- fastx
- other
- paired_in - put both reads into the 'aligned.fasta' file
- paired_out - put both reads into the 'other.fasta' file

The table 3 shows the output generated depending on selected options. 

+------+--------------------------------+-------------------------------------------------+
|case #|options selected                |Number of reads per pair written into output file|
+======+================================+=================================================+
|-     |fastx|other|paired_in|paired_out|aligned.fasta|other.fasta                        |
+------+-----+-----+---------+----------+-------------+-----------------------------------+
|1     |true |false|false    |false     |1            |0                                  |
+------+-----+-----+---------+----------+-------------+-----------------------------------+
|2     |true |false|true     |false     |2            |0                                  |
+------+-----+-----+---------+----------+-------------+-----------------------------------+
|3     |true |false|false    |true      |0            |0                                  |
+------+-----+-----+---------+----------+-------------+-----------------------------------+
|4     |true |true |true     |false     |2            |0                                  |
+------+-----+-----+---------+----------+-------------+-----------------------------------+
|5     |true |true |false    |true      |0            |2                                  |
+------+-----+-----+---------+----------+-------------+-----------------------------------+

::

    Table 3

Case 1 results in the splitting of some paired reads in the output files and not optimal for users who require paired order of the reads for downstream analyses.

The option 'paired_in' is optimal for users who want all reads in the :code:`other.fasta` file to be non-rRNA. However, there are small chances that reads which are non-rRNA will also be put into the :code:`aligned.fasta` file.

The option 'paired_out' is optimal for users who want only rRNA reads in the :code:`aligned.fasta` file. However, there are small chances that reads which are rRNA will also be put into the 'other.fasta' file.

If neither of these two options is used, then aligned and non-aligned reads will be properly output to the :code:`aligned.fasta` and :code:`other.fasta` files, possibly breaking the order for a set of paired reads between two output files.

It's important to note that regardless of the options used, the :code:`aligned.log` file will always report the true number of reads classified as rRNA (not the number of reads in the :code:`aligned.fasta` file). 
Mapping reads for classification
################################

Although SortMeRNA is very sensitive with the small rRNA databases distributed with the source code, these databases are not optimal for classification since often alignments with 75-90% identity will be returned (there are only several thousand rRNA in most of the databases, compared to the original SILVA or Greengenes databases containing millions of rRNA). Classification at the species level generally considers alignments at 97% and above, so it is suggested to use a larger database is species classification is the main goal.

Moreover, SortMeRNA is a local alignment tool, so it's also important to look at the query coverage % for each alignment. In the SAM output format, neither % id or query coverage are reported. If the user wishes for these values, then the Blast tabular format with CIGAR + query coverage option (--blast '1 cigar qcov') is the way to go. 

OTU-picking
-----------

SortMeRNA is implemented in QIIME's closed-reference and open-reference OTU-picking workflows. The readers are referred to QIIME's tutorials for an in-depth discussion of these methods http://qiime.org/tutorials/otu_picking.html. 

Advanced
--------

::

    -num seeds INT



The threshold number of seeds required to match in the primary seed-search filter before moving on to the secondary seed-cluster filter. More specically, the threshold number of seeds required before searching for a longest increasing subsequence (LIS) of the seeds' positions between the read and the closest matching reference sequence. By default, this is set to 2 seeds. 

::

    passes INT,INT,INT

In the primary seed-search filter, SortMeRNA moves a seed of length L (parameter of indexdb rna) across the read using three passes. If at the end of each pass a threshold number of seeds (defined by --num seeds) did not match to the reference database, SortMeRNA attempts to find more seeds by decreasing the interval at which the seed is placed along the read by using another pass. In default mode, these intervals are set to L,L/2,3 for Pass 1, 2 and 3, respectively. Usually, if the read is highly similar to the reference database, a threshold number of seeds will be found in the first pass. 

::

    -edges INT(%)



The number (or percentage if followed by %) of nucleotides to add to each edge of the alignment region on the reference sequence before performing Smith-Waterman alignment. By default, this is set to 4 nucleotides. 

:: 

    -full_search flag

During the index traversal, if a seed match is found with 0-errors, SortMeRNA will stop searching for further 1-error matches. This heuristic is based upon the assumption that 0-error matches are more signicant than 1-error matches. By turning it off using the --full_search flag, the sensitivity may increase (often by less than 1%) but with up to four-fold decrease in speed. 

::

    -pid FLAG

The pid of the running sortmerna process will be added to the output files in order to avoid over-writing output if the same --aligned STRING base name is provided for different runs. 

Help 
----

Any issues or bug reports should be reported to https://github.com/biocore/sortmerna/issues. Comments and suggestions are also always appreciated! 

Citation
--------

If you use SortMeRNA please cite,

Kopylova E., Noe L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611.

Copyright (C) 2016-2019 Clarity Genomics BVBA Turnhoutseweg 30, 2340 Beerse, Belgium http://www.clarity-genomics.com

Copyright (C) 2014-2016 Knight Lab Department of Pediatrics, UCSD School of Medicine, La Jolla, California, USA https://knightlab.colorado.edu

Copyright (C) 2012-2014 Bonsai Bioinformatics Research Group (LIFL - Université Lille 1), CNRS UMR 8022, INRIA Nord-Europe, France http://bioinfo.lifl.fr/RNA/sortmerna/