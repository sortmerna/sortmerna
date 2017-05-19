sortmerna
=========

[![Build Status](https://travis-ci.org/biocore/sortmerna.png?branch=master)](https://travis-ci.org/biocore/sortmerna)

SortMeRNA is a local sequence alignment tool for filtering, mapping and clustering.

The core algorithm is based on approximate seeds and allows for sensitive analysis of NGS reads.
The main application of SortMeRNA is filtering rRNA from metatranscriptomic data.
SortMeRNA takes as input a file of reads (fasta or fastq format) and one or multiple
rRNA database file(s), and sorts apart aligned and rejected reads into two files specified by the user.
Additional applications include clustering and taxonomy assignation available through QIIME v1.9.1
(http://qiime.org). SortMeRNA works with Illumina, Ion Torrent and PacBio data, and can produce SAM and
BLAST-like alignments.

Visit http://bioinfo.lifl.fr/RNA/sortmerna/ for more information.


# Table of Contents
* [Support](#support)
* [Documentation](#documentation)
* [Getting Started](#getting-started)
* [Compilation](#sortmerna-compilation)
	* [Linux OS](#linux-os)
	* [Mac OS](#mac-os)
* [Install compilers, ZLIB and autoconf](#install-compilers-zlib-and-autoconf)
	* [Clang for Mac OS](#clang-for-mac-os)
	* [GCC and Zlib though MacPorts](#gcc-and-zlib-though-macports)
	* [autoconf](#autoconf)
* [Tests](#tests)
* [Third-party libraries](#third-party-libraries)
* [Wrappers and packages](#wrappers-and-packages)
	* [Galaxy](#galaxy)
	* [Debian](#debian)
	* [GNU Guix](#gnu-guix)
	* [QIIME](#qiime)
* [Taxonomies](#taxonomies)
* [Citation](#citation)
* [Contributors](#contributors)


# Support
For questions and comments, please use the SortMeRNA [forum](https://groups.google.com/forum/#!forum/sortmerna).

  
# Documentation

If you have [Doxygen](http://www.stack.nl/~dimitri/doxygen/) installed, you can generate the documentation
by modifying the following lines in ```doxygen_configure.txt```:

```
INPUT = /path/to/sortmerna/include /path/to/sortmerna/src
IMAGE_PATH = /path/to/sortmerna/algorithm
```

and running the following command:

```
doxygen doxygen_configure.txt
```

This command will generate a folder ```html``` in the directory from which the
command was run.


# Getting Started

There are 3 methods to install SortMeRNA:

1. [GitHub repository](https://github.com/biocore/sortmerna) development version (master branch)
...* [Installation instructions](#sortmerna-compilation)
2. [GitHub releases](https://github.com/biocore/sortmerna/releases) (tar balls)
...* [Installation instructions Linux](#linux-os)
...* [Installation instructions Mac OS](#mac-os)
3. [BioInfo releases](http://bioinfo.lifl.fr/RNA/sortmerna/) (tar balls including compiled binaries)

Option (3) is the most straight-forward, as it does not require ```autoconf``` and provides
access to pre-compiled binaries to various OS.

# SortMeRNA Compilation

NOTE: You will require ```autoconf``` to build from the git cloned
repository or from source code in the `Source code` tar
balls under release GitHub `Downloads`.

(0) Prepare your build system for compilation:

```bash
bash autogen.sh
```

## Linux OS

(1) Check your GCC compiler is version 4.0 or above:

```bash
gcc --version
```

(2) Run configure and make scripts:

```bash
./configure LDFLAGS="-L/usr/lib64"
make
```

(3) To install:

```bash
make install
```

You can define an alternative installation directory by
specifying ```--prefix=/path/to/installation/dir``` to ```configure```.


## Mac OS

(1) Check the version of your C/C++ compiler:

```bash
gcc --version
```

(2a) If the compiler is Clang or GCC, proceed to run configure and make scripts:

```bash
./configure
make
```

Note: If the compiler is Clang, you will not have access to multithreading.

(2b) If the compiler is LLVM-GCC, you will need to change it
(see [Deprecation and Removal Notice](https://developer.apple.com/library/ios/documentation/DeveloperTools/Conceptual/WhatsNewXcode/Articles/xcode_5_0.html)).

To set your compiler to Clang (see [instructions](#set-clang-compiler-for-mac-os))
or the original GCC compiler (see [instructions](#set-gcc-compiler-for-mac-os)).

(3) To install:

```bash
make install
```

### Set Clang compiler for Mac OS

(1) Check if you have Clang installed:

```bash
clang --version
```

(2a) If Clang is installed, set your compiler to Clang:

```bash
export CC=clang
export CXX=clang++
```

(2b) If Clang is not installed, see [Clang for Mac OS](#clang-for-mac-os)
for installation instructions.

### Set GCC compiler for Mac OS

(1) Check if you have GCC installed:

```bash
gcc --version
```

(2a) If GCC is installed, set your compiler to GCC:

```bash
export CC=gcc-mp-4.8
export CXX=g++-mp-4.8
```

(2b) If GCC is not installed, see [Install GCC through MacPorts](#gcc-and-zlib-though-macPorts)
for installation instructions.

(3) Next, if you would like zlib support (reading compressed .zip and .gz FASTA/FASTQ files), Zlib
should also be installed via MacPorts. See section [Install GCC and Zlib though MacPorts](#install-gcc-and-zlib-though-macports)
for installation instructions.

(4a) Assuming you have Zlib installed, run configure and make scripts
(if compression feature wanted):

```bash
./configure --with-zlib="/opt/local"
make
```

(4b) Otherwise (if option to read compressed files is not wanted):

```bash
./configure --without-zlib
make
```

You can define an alternative installation directory by
specifying ```--prefix=/path/to/installation/dir``` to ```configure```.


Install compilers, ZLIB and autoconf
====================================

NOTE: the Clang compiler on Mac (distributed through Xcode) does not support OpenMP (multithreading).
A preliminary implementation of OpenMP for Clang has been made at "http://clang-omp.github.io"
though has not been yet incorporated into the Clang mainline. The user may follow the
steps outlined in the above link to install the version of Clang with multithreading support, 
though this version has not yet been tested with SortMeRNA. Otherwise, the user is 
recommended to install the original GCC compiler via MacPorts (contains full multithreading support).

Clang for Mac OS
----------------

Installing Xcode (free through the App Store) and Xcode command line tools will automatically 
install the latest version of Clang supported with Xcode. 

After installing Xcode, the Xcode command line tools may be installed via:

Xcode -> Preferences -> Downloads

Under "Components", click to install "Command Line Tools"


GCC and Zlib though MacPorts
----------------------------

Assuming you have MacPorts installed, type:

```bash
sudo port selfupdate
sudo port install gcc48
sudo port install zlib
```

After the installation, you should find the compiler installed in /opt/local/bin/gcc-mp-4.8 and /opt/local/bin/g++-mp-4.8
as well as Zlib in /opt/local/lib/libz.dylib and /opt/local/include/zlib.h .

autoconf
--------

```bash
wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz
tar -zxvf autoconf-2.69.tar.gz
cd autoconf-2.69
./configure
make
make install
```

You can define an alternative installation directory by
specifying ```--prefix=/path/to/installation/dir``` to ```configure```
(before calling ```make```).

If installing in a directory other than those listed in $PATH,
add the installation directory to $PATH:

```bash
export PATH=$PATH:/path/to/installation/dir
```

Tests
=====

Usage tests can be run with the following command:
```
python ./tests/test_sortmerna.py
python ./tests/test_sortmerna_zlib.py
```
Make sure the ```data``` folder is in the same directory as ```test_sortmerna.py```

Users require [scikit-bio](https://github.com/biocore/scikit-bio) 0.5.0 to run the tests.


Third-party libraries
=====================
Various features in SortMeRNA are dependent on third-party libraries, including:
* [ALP](http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/software/program.html?uid=6): computes statistical parameters for Gumbel distribution (K and Lambda)
* [CMPH](http://cmph.sourceforge.net): C Minimal Perfect Hashing Library
* [KSEQ](http://lh3lh3.users.sourceforge.net/parsefastq.shtml): FASTA/FASTQ parser (including compressed files)
* [PARASAIL](https://github.com/jeffdaily/parasail): Pairwise Sequence Alignment Library

Wrappers and Packages
=====================

Galaxy
------

Thanks to Björn Grüning and Nicola Soranzo, an up-to-date Galaxy wrapper exists for SortMeRNA.
Please visit Björn's [github page](https://github.com/bgruening/galaxytools/tree/master/tools/rna_tools/sortmerna) for installation.

Debian
------

Thanks to the [Debian Med](https://www.debian.org/devel/debian-med/) team, SortMeRNA 2.0 is now a package in Debian.
Thanks to Andreas Tille for the sortmerna and indexdb_rna man pages (version 2.0).
These have been updated for 2.1 in the master repository.

GNU Guix
--------

Thanks to Ben Woodcroft for adding SortMeRNA 2.1 to GNU Guix, find the package [here](https://www.gnu.org/software/guix/packages/).

QIIME
-----

SortMeRNA 2.0 can be used in [QIIME](http://qiime.org)'s [pick_closed_reference_otus.py](http://qiime.org/scripts/pick_closed_reference_otus.html),
[pick_open_reference_otus.py](http://qiime.org/scripts/pick_open_reference_otus.html) and [assign_taxonomy.py](http://qiime.org/scripts/assign_taxonomy.html) scripts.

Note: At the moment, only 2.0 is compatible with QIIME.

Taxonomies
==========

The folder `rRNA_databases/silva_ids_acc_tax.tar.gz` contains SILVA taxonomy strings (extracted from XML file generated by ARB)
for each of the reference sequences in the representative databases. The format of the files is three tab-separated columns,
the first being the reference sequence ID, the second being the accession number and the final column is the taxonomy.

Citation
========

If you use SortMeRNA, please cite:
Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611.

Contributors
============
See [AUTHORS](./AUTHORS) for a list of contributors to this project.

