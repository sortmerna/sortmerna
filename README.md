sortmerna
=========

[![Build Status](https://travis-ci.org/biocore/sortmerna.png?branch=master)](https://travis-ci.org/biocore/sortmerna)

This is ongoing development of SortMeRNA. See 'release' on this page for
releases or visit http://bioinfo.lifl.fr/RNA/sortmerna/ for more information.

WARNING: Do not run `_autotools_batch_.sh` unless you have Autotools installed. 
This script is for producing a distribution version of the code.

NOTE: the Clang compiler on Mac (distributed through Xcode) does not support multithreading.
A preliminary implementation of OpenMP for Clang has been made at "http://clang-omp.github.io"
though has not been yet incorporated into the Clang mainline. The user may follow the
steps outlined in the above link to install the version of Clang with multithreading support, 
though this version has not yet been tested with SortMeRNA. Otherwise, the user is 
recommended to install the original GCC compiler via MacPorts (contains full multithreading support).
  
Documentation:
--------------

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


To compile on Linux OS:
-----------------------

(1) Check your GCC compiler is version 4.0 or higher:

```bash
gcc --version
```

(2) Run configure and make scripts:

```bash
bash ./build.sh
```

Note: `make install` is not called in this script. However, any arguments given
to `build.sh` will be passed to the configure script. If you plan on calling `make install`
afterwards, then you can set your installation directory using
`build.sh --prefix=/path/to/installation/dir`. Otherwise, you can simply copy the
binaries `sortmerna` and `indexdb_rna` to your installation directory after
calling `build.sh`.



To compile on Mac OS:
---------------------

(1) Check the version of your C/C++ compiler:

```bash
gcc --version
```

If the compiler is LLVM-GCC (deprecated, see [Deprecation and Removal Notice](https://developer.apple.com/library/ios/documentation/DeveloperTools/Conceptual/WhatsNewXcode/Articles/xcode_5_0.html)), 
then you must set it to Clang or the original GCC compiler (installable via MacPorts).

(2a) To set your compiler to Clang (check you have it "clang --version", if not then 
see "Install Clang for Mac OS" below):

```bash
export CC=clang
export CXX=clang++
```

(2b) To set your compiler to the original GCC (check you have it "gcc-mp-4.8 --version", 
if not then see "Install GCC through MacPorts" below. Note the GCC compiler comes in many 
versions, 4.8 is one of the latest):

```bash
export CC=gcc-mp-4.8
export CXX=g++-mp-4.8
```

(3) Run configure and make scripts:

```bash
bash ./build.sh
```

Note: `make install` is not called in this script. However, any arguments given
to `build.sh` will be passed to the configure script. If you plan on calling `make install`
afterwards, then you can set your installation directory using
`build.sh --prefix=/path/to/installation/dir`. Otherwise, you can simply copy the
binaries `sortmerna` and `indexdb_rna` to your installation directory after
calling `build.sh`.


Install Clang for Mac OS 
------------------------

Installing Xcode (free through the App Store) and Xcode command line tools will automatically 
install the latest version of Clang supported with Xcode. 

After installing Xcode, the Xcode command line tools may be installed via:

Xcode -> Preferences -> Downloads

Under "Components", click to install "Command Line Tools"


Install GCC though MacPorts
---------------------------

Assuming you have MacPorts installed, type:

```bash
sudo port selfupdate
sudo port install gcc48
```

After the installation, you should find the compiler installed in /opt/local/bin/gcc-mp-4.8 and /opt/local/bin/g++-mp-4.8 .


Tests
-----

Usage tests can be run with the following command:
```
python ./tests/test_sortmerna.py
```
Make sure the ```data``` folder is in the same directory as ```test_sortmerna.py```


Galaxy Wrapper
--------------

Currently a Galaxy wrapper exists for SortMeRNA 2.0 (thanks to Björn Grüning and Nicola Soranzo,
see this [PR](https://github.com/bgruening/galaxytools/pull/202) for more details).
Please visit "https://github.com/bgruening/galaxytools/tree/master/tools/rna_tools/sortmerna" for installation.

Debian package
--------------

Thanks to the [Debian Med](https://www.debian.org/devel/debian-med/) team, SortMeRNA 2.0 is now a package in Debian.
Thanks to Andreas Tille for the sortmerna and indexdb_rna man pages (version 2.0).
These have been updated for 2.1 in the master repository.

GNU Guix package
----------------

Thanks to Ben Woodcroft for adding SortMeRNA 2.1 to GNU Guix, find the package [here](https://www.gnu.org/software/guix/packages/).

Running in QIIME
----------------

SortMeRNA 2.0 can be used in [QIIME](http://qiime.org)'s [pick_closed_reference_otus.py](http://qiime.org/scripts/pick_closed_reference_otus.html),
[pick_open_reference_otus.py](http://qiime.org/scripts/pick_open_reference_otus.html) and [assign_taxonomy.py](http://qiime.org/scripts/assign_taxonomy.html) scripts.

Note: At the moment, only 2.0 is compatible with QIIME.

Taxonomies
----------

The folder `rRNA_databases/silva_ids_acc_tax.tar.gz` contains SILVA taxonomy strings (extracted from XML file generated by ARB)
for each of the reference sequences in the representative databases. The format of the files is three tab-separated columns,
the first being the reference sequence ID, the second being the accession number and the final column is the taxonomy.

Citation
--------

If you use SortMeRNA, please cite:
Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611.



