==========
User guide
==========

**SortMeRNA** is a local sequence alignment tool for filtering, mapping and OTU clustering. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of NGS reads. The main application of SortMeRNA is filtering rRNA from metatranscriptomic data. Additional applications include OTU-picking and taxonomy assignation available through `QIIME v1.9+
<http://qiime.org>`_. SortMeRNA takes as input one or two (paired) reads file(s) (fasta/fasta.gz/fastq/fastq.gz), and one or multiple rRNA database file(s), and sorts apart aligned and rejected reads into two files. SortMeRNA works with Illumina, 454, Ion Torrent and PacBio data, and can produce SAM and BLAST-like alignments.

.. note::
   
   This project is under active development.

Installation
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

The only required options are :code:`--ref` and :code:`--reads`. Options (any) can be specified usig a single dash e.g. :code:`-ref` and :code:`-reads`. Both plain :code:`fasta/fastq` and archived :code:`fasta.gz/fastq.gz` files are accepted. File extensions :code:`.fastq`, :code:`.fastq.gz`, :code:`.fq`, :code:`.fq.gz`, :code:`.fasta`, ... are optional. The format and compression are automatically recognized. Relative paths are accepted.

Example 1
---------

Single reference and single reads file::

   sortmerna --ref REF_PATH --reads READS_PATH

Example 2
---------

For multiple references use multiple :code:`--ref`::

   sortmerna --ref REF_PATH_1 --ref REF_PATH_2 --ref REF_PATH_3 --reads READS_PATH

Example 3
---------

For Paired reads use :code:`--reads` twice::

   sortmerna --ref REF_PATH_1 --ref REF_PATH_2 --ref REF_PATH_3 --reads READS_PATH_1 --reads READS_PATH_2

Reference Databases
===================

Download latest SortMeRNA databases (v4.3.6)::

   wget https://github.com/biocore/sortmerna/releases/download/v4.3.6/database.tar.gz
   mkdir rRNA_databases_v4.3.6
   tar -xvf database.tar.gz -C rRNA_databases_v4.3.6

Building
========

General Notes
-------------

General build steps:

1. Prepare the build environment
2. Get the sources from GitHub or from GitHub releases
3. Build

SortmeRNA-4 is C++17 compliant, and requires a compiler that supports filesystem standard library. Currently only GCC-9 and Windows SDK meet that requirement.

Ready to use GCC-9 distribution is only available on Debian (Ubuntu). It is not available on Centos, i.e. no Devtoolset-9 yet, and building GCC-9 on Centos is not trivial.

Clang has no support for the filesystem yet. LLVM 9.0 is due to be released later this year. At that time we'll get back to supporting builds with Clang on Linux and OSX.

The build is performed using a Python script provided with the Sortmerna distribution. The script uses a configuration file env.yaml, which can be Optionally modified to customize the build.

Building on Linux
-----------------

Quick Build
^^^^^^^^^^^

If you already have GCC, CMake and Conda, the following command will compile SortMeRNA local sources in the current working directory on a linux machine::

   python scripts/build.py --name all --local-linux

Install GCC 9
^^^^^^^^^^^^^

This is for Debian distros (Ubuntu)::

   sudo add-apt-repository ppa:ubuntu-toolchain-r/test
   sudo apt update
   sudo apt -y install gcc-9 g++-9
   sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-9
   sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 9
   sudo update-alternatives --install /usr/bin/cpp cpp-bin /usr/bin/cpp-9 60
   
   # select gcc-9
   sudo update-alternatives --config gcc
   # select cpp-9
   sudo update-alternatives --config cpp-bin
   
   gcc --version
     gcc (Ubuntu 9.2.1-17ubuntu1~16.04) 9.2.1 20191102

Get SortMeRNA sources
^^^^^^^^^^^^^^^^^^^^^

The sources can be placed in any directory, but here we use the user's Home directory::

   # clone the repository
   git clone https://github.com/biocore/sortmerna.git
   
   # alternatively get the release sources
   wget https://github.com/biocore/sortmerna/archive/v4.0.0.tar.gz
   
   tar xzf v4.0.0.tar.gz
   
   pushd sortmerna
   
   # If you need a particular release (tag)
   git checkout v4.0.0

Install Conda
^^^^^^^^^^^^^

Use the :code:`build.py` python script provided with Sortmerna distro. The following installs Conda, and the python packages :code:`pyyaml`, and :code:`jinja2` in the User's Home directory::
   
   SMR_HOME=$HOME/sortmerna
   python $SMR_HOME/scripts/build.py --name cmake
   
   ls -lrt
   drwxrwxr-x 15 biocodz biocodz     4096 Nov 18 09:43 miniconda3
   
   # add Conda binaries to the PATH
   export PATH=$HOME/miniconda3/bin:$PATH

Install CMake
^^^^^^^^^^^^^

The following installs CMake in user's home directory::

   SMR_HOME=$HOME/sortmerna
   python $SMR_HOME/scripts/build.py --name cmake
     [cmake_install] Installed CMake /home/biocodz/cmake-3.15.5-Linux-x86_64/bin/cmake
   
   # add cmake to PATH
   export PATH=$HOME/cmake-3.15.5-Linux-x86_64/bin:$PATH

Build
^^^^^

All required third party libraries will be checked and installed automatically (in User directory by default) The default build won't interfere with any existing system installation. By default the build produces statically linked executable i.e. portable.

::

   SMR_HOME=$HOME/sortmerna
   
   # modify configuration (optional)
   vi $SMR_HOME/scripts/env.yaml
   
   # run the build
   python $SMR_HOME/scripts/build.py --name all [--env $SMR_HOME/script/my_env.yaml]
   
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
