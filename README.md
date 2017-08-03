SortMeRNA
=========

[![Build Status](https://travis-ci.org/biocore/sortmerna.png?branch=master)](https://travis-ci.org/biocore/sortmerna)

SortMeRNA is a local sequence alignment tool for filtering, mapping and clustering high throughput sequences.

SortMeRNA is supported for Linux, Mac and Windows.

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
		* [Install GCC through MacPorts](#install-gcc-though-macports)
		* [Install Clang for Mac OS](#install-clang-for-mac-os)
	* [Windows OS](#windows-os)
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

SortMeRNA can be built and run on Linux, Mac and Windows:

* [GitHub development version](#sortmerna-compilation) master branch
* [GitHub releases](https://github.com/biocore/sortmerna/releases) (tar balls, zip)
	* [Linux](#linux-os)
	* [Mac OS](#mac-os)
	* [Windows OS](#windows-os)


# SortMeRNA Compilation

CMake is used for build files generation and should be installed prior the build (version 3.0 or above required).
CMake distributions are available for all major operating systems.
Please visit [CMake project website](https://cmake.org/) for download and installation instructions.

## Linux OS

(1) Check your GCC compiler is version 4.0 or above:

```bash
gcc --version
```

(2) Set environment variables CC and CXX to use gcc:
```bash
export CC=gcc
export CXX=g++
```

(3) Clone sortmerna repository, or download a non-binary release:

```bash
git clone https://github.com/biocore/sortmerna.git
```

(4) Set SortMeRNA home path:

```bash
pwd
/path/to/sortmerna
SMR_HOME=/path/to/sortmerna
```

(5) Generate the build files:

```bash
mkdir -p $SMR_HOME/build/Release
pushd $SMR_HOME/build/Release
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release ../..
```
`$SMR_HOME` is the top directory where sortmerna code (e.g. git repo) is located (see step (4) on how to set).

The above commands will perform necessary system check-ups, dependencies, and generate Makefile.

(6) Compile and build executables:

```bash
make
```

The binaries are created in `$SMR_HOME/build/Release/src/indexdb` and `$SMR_HOME/build/Release/src/sortmerna`
Simply add the build binaries to the PATH:

```bash
export PATH="$SMR_HOME/build/Release/src/indexdb:$SMR_HOME/build/Release/src/sortmerna:$PATH"
```

## Mac OS

For Mac OS, you can use compiler g++/gcc (supports multithreading) or Clang (does NOT support multithreading).

(1) Check if you have GCC installed:

```bash
gcc --version
```

(2a) If installed, set your compiler to GCC:

```bash
export CC=gcc-mp-4.8
export CXX=g++-mp-4.8
```

(2b) If GCC is not installed, see [Install GCC through MacPorts](#install-gcc-though-macPorts)
for installation instructions. Or, use the Clang compiler (though multithreading will NOT be supported,
see [Install Clang for Mac OS](#install-clang-for-mac-os)):

```bash
clang --version
export CC=clang
export CXX=clang++
```

(3) Follow steps (3)-(6) above for [Linux OS](#linux-os).

### Install GCC through MacPorts

Assuming you have MacPorts installed, type:

```bash
sudo port selfupdate
sudo port install gcc48
```

After the installation, you should find the compiler installed in /opt/local/bin/gcc-mp-4.8 and /opt/local/bin/g++-mp-4.8

### Install Clang for Mac OS

Installing Xcode (free through the App Store) and Xcode command line tools will automatically 
install the latest version of Clang supported with Xcode. 

After installing Xcode, the Xcode command line tools may be installed via:

Xcode -> Preferences -> Downloads

Under "Components", click to install "Command Line Tools"

## Windows OS

MS Visual Studio Community edition and CMake for Windows are required for building SortMeRNA.
Download and Install VS Community edition from [Visual Studio community website](https://www.visualstudio.com/vs/community/)
The following assumes `Visual Studio 14 2015`.

Open Win CMD (command shell)

```bash
mkdir %SMR_HOME%\build
pushd %SMR_HOME%\build
cmake -G "Visual Studio 14 2015 Win64" ..
```

The above generates VS project files in `%SMR_HOME%\build\` directory. It also downloads required 3rd party source packages like `zlib` (in `%SMR_HOME%\3rdparty\`).
`%SMR_HOME%` is the top directory where SortMeRNA source distribution (e.g. Git repo) is installed.

Start Visual Studio and open Sortmerna solution
`File -> Open -> Project/Solution .. open %SMR_HOME%\build\sortmerna.sln`

Select desired build type: `Release | Debug | RelWithDebInfo | MinSizeRel`.
In Solution explorer right-click `ALL_BUILD' and select `build` in pop-up menu.

Depending on the build type the binaries are generated in 
`%SMR_HOME%\build\src\sortmerna\Release` (or `Debug | RelWithDebInfo | MinSizeRel`).

Add sortmerna executables to PATH

```bash
set PATH=%SMR_HOME%\build\src\indexdb\Release;%SMR_HOME%\build\src\sortmerna\Release;%PATH%
```


Tests
=====

Python code is provided for running tests in $SRM_HOME/tests (%SRM_HOME%\tests) and requires Python 3.5 or higher.

Tests can be run with the following command:

```bash
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
* [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library): SSW Library for SIMD Smith-Waterman

Wrappers and Packages
=====================

Galaxy
------

Thanks to Björn Grüning and Nicola Soranzo, a Galaxy wrapper exists for SortMeRNA.
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

