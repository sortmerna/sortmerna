==========
User guide
==========

**SortMeRNA** is a local sequence alignment tool for filtering, mapping and OTU clustering. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of NGS reads. The main application of SortMeRNA is filtering rRNA from metatranscriptomic data. Additional applications include OTU-picking and taxonomy assignation available through `QIIME v1.9+
<http://qiime.org>`_. SortMeRNA takes as input one or two (paired) reads file(s) (fasta/fasta.gz/fastq/fastq.gz), and one or multiple rRNA database file(s), and sorts apart aligned and rejected reads into two files. SortMeRNA works with Illumina, 454, Ion Torrent and PacBio data, and can produce SAM and BLAST-like alignments.

.. note::
   
   This project is under active development.

Instalation
==========

Using Conda
-----------

Install Conda::

   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh

As per the `Bioconda guidelines
<https://bioconda.github.io/>`_, add the following conda channels::

   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict

Then install SortMeRNA::

   conda search sortmerna
     Loading channels: done
     # Name                       Version           Build  Channel
     sortmerna                        2.0               0  bioconda
     ...
     sortmerna                      4.3.4               0  bioconda
     ...
     sortmerna                      4.3.6               0  bioconda

   # create a new environment and install SortMeRNA in it
   conda create --name sortmerna_env
   conda activate sortmerna_env
   conda install sortmerna
   which sortmerna
     /home/biocodz/miniconda3/envs/sortmerna_env/bin/sortmerna

   # test the installation
   sortmerna --version
     SortMeRNA version 4.3.6
     Build Date: Aug 16 2022
     sortmerna_build_git_sha:@db8c1983765f61986b46ee686734749eda235dcc@
     sortmerna_build_git_date:@2022/08/16 11:42:59@

   # view help
   sortmerna -h

Basic usage
===========

Example 1
---------

Single reference and single reads file::

   sortmerna --ref REF_PATH --reads READS_PATH

Download latest SortMeRNA databases (v4.3.6)::

   wget https://github.com/biocore/sortmerna/releases/download/v4.3.6/database.tar.gz
   mkdir rRNA_databases_v4.3.6
   tar -xvf database.tar.gz -C rRNA_databases_v4.3.6

Download test reads file::

   wget https://github.com/biocore/sortmerna/blob/master/data/test_read.fasta

Run with single reads file::

   sortmerna --ref rRNA_databases_v4.3.6/smr_v4.3_default_db.fasta --reads test_read.fasta

.. note::
   
   Index file for ``rRNA_databases_v4.3.6/smr_v4.3_default_db.fasta`` created.
