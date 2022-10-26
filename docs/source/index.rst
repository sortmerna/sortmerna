==========
User guide
==========

**SortMeRNA** is a local sequence alignment tool for filtering, mapping and clustering.

.. note::
   
   This project is under active development.

Basic usage
===========

Example 1
---------

Single reference and single reads file::

   sortmerna --ref REF_PATH --reads READS_PATH

Download latest SortMeRNA databases (v4.3.4)::

   wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
   mkdir rRNA_databases_v4.3.4
   tar -xvf database.tar.gz -C rRNA_databases_v4.3.4

Download test reads file::

   wget https://github.com/biocore/sortmerna/blob/master/data/test_read.fasta

Run with single reads file::

   sortmerna --ref rRNA_databases_v4.3.4/smr_v4.3_default_db.fasta --reads test_read.fasta
