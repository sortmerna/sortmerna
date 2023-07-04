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

Help
----

Any issues or bug reports should be reported to https://github.com/biocore/sortmerna/issues. Comments and suggestions are also always appreciated!

Citation
--------

If you use SortMeRNA please cite,

Kopylova E., Noe L. and Touzet H., “SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data”, Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611.

Copyright (C) 2016-2019 Clarity Genomics BVBA Turnhoutseweg 30, 2340 Beerse, Belgium http://www.clarity-genomics.com

Copyright (C) 2014-2016 Knight Lab Department of Pediatrics, UCSD School of Medicine, La Jolla, California, USA https://knightlab.colorado.edu

Copyright (C) 2012-2014 Bonsai Bioinformatics Research Group (LIFL - Université Lille 1), CNRS UMR 8022, INRIA Nord-Europe, France http://bioinfo.lifl.fr/RNA/sortmerna/

Contents
--------

.. toctree::

   installation
   building
   usage
   databases
   trace4.3.2.rst
   contributors
   citation