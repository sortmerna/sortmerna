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