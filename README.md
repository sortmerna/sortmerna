# sortmerna

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

* [Getting Started](#getting-started)
	* [Attention the new sortmerna '-d' option](#beware-issue-171)
	* [DNANexus cloud](#dnanexus-cloud)
	* [Using GitHub release binaries](#using-github-release-binaries)
	* [Building from source code](#building-from-source-code)
		* [General notes](#general-notes)
		* [Obtaining the source code](#obtaining-the-source-code)
			* [Download GitHub Release](#download-github-release)
			* [Clone the GitHub Repository and use the head or a release](#clone-the-github-repository-and-use-the-head-or-a-release)
		* [building on Linux OS](#building-on-linux-os)
		* [Building on Mac OS](#building-on-mac-os)
			* [Install Clang on Mac](#install-clang-on-mac)
			* [Configure shell to use Clang compiler on Mac](#configure-shell-to-use-clang-compiler-on-mac)
			* [Configure shell to use GCC compiler on Mac](#configure-shell-to-use-gcc-compiler-on-Mac)
		* [Building on Windows OS](#building-on-windows-os)
* [Running tests](#running)
* [User Manual](#user-manual)
* [Third-party libraries](#third-party-libraries)
* [Taxonomies](#taxonomies)
* [Citation](#citation)
* [Contributors](#contributors)
* [Support](#support)
* [Documentation](#documentation)
* [References](#references)


# Getting Started

SortMeRNA can be built and run on Windows, Linux, and Mac.

The following methods can be used for building/running SortMeRNA:

## Beware Issue 171

Beware that the directory specified using `-d` options is completely wiped out on each `sortmerna` run. Refer to the [Issue 171](https://github.com/biocore/sortmerna/issues/171)

## DNANexus cloud

Arguably the easiest way to use Sortmerna.

* Ready-to-use `sortmerna` application is available on DNANexus cloud. Just upload your data and run the application.
* Sortmerna distribution also contains DNANexus applets for building, testing and running Sortmerna. See [dnanexus](https://github.com/biocore/sortmerna/blob/master/dnanexus/README.md) subdirectory for detailed instructions.

Note that DNANexus trial accounts can be used.

## Using GitHub release binaries
Visit [Sortmerna GitHub Releases](https://github.com/biocore/sortmerna/releases)

The file [sortmerna-3.0.0-linux.tar.gz](https://github.com/biocore/sortmerna/releases/download/v3.0.0/sortmerna-3.0.0-linux.tar.gz) contains 3 binaries
`indexdb`, `libstdc++.so.6`, `sortmerna`

`libstdc++.so.6` is a dependency of `sortmerna`. It is needed on machines that have no GCC installed.

The binaries were built on DNANexus Ubuntu 16.04 VM using [sortmerna-3.build.on.u16asset](https://github.com/biocore/sortmerna/tree/master/dnanexus/applets/sortmerna-3.build.on.u16asset) and the [sortmerna-3.asset](https://github.com/biocore/sortmerna/tree/master/dnanexus/assets/sortmerna-3.asset)

The `sortmerna` binary was patched using `patchelf` utility to search for `libstdc++.so.6` in the same directory where `sortmerna` is located.

Put all the binaries into the same location e.g. `SORTMERNA_HOME/bin`, and add the `bin` to the PATH.

Below are all the necessary `bash` commands:

```
# get binary distro
wget https://github.com/biocore/sortmerna/releases/download/v3.0.0/sortmerna-3.0.0-linux.tar.gz

# get sources
wget https://github.com/biocore/sortmerna/archive/v3.0.0.tar.gz

# list the source tar content
less v3.0.0.tar.gz
drwxrwxr-x root/root         0 2018-10-29 08:13 sortmerna-3.0.0/
-rw-rw-r-- root/root       119 2018-10-29 08:13 sortmerna-3.0.0/.gitignore
-rw-rw-r-- root/root      4909 2018-10-29 08:13 sortmerna-3.0.0/.travis.yml
drwxrwxr-x root/root         0 2018-10-29 08:13 sortmerna-3.0.0/3rdparty/
drwxrwxr-x root/root         0 2018-10-29 08:13 sortmerna-3.0.0/3rdparty/alp/
-rw-rw-r-- root/root      1344 2018-10-29 08:13 sortmerna-3.0.0/3rdparty/alp/CMakeLists.txt
...

# extract sources
tar xzf v3.0.0.tar.gz

# create bin/
mkdir sortmerna-3.0.0/bin

# extract binaries into bin/
tar xzf sortmerna-3.0.0-linux.tar.gz -C sortmerna-3.0.0/bin/

# list the binaries
ls -lrt sortmerna-3.0.0/bin/
-rwxrwxr-x 1 biocodz biocodz  174560 Oct 30 06:48 indexdb
-rw-rw-r-- 1 biocodz biocodz 1607120 Oct 30 06:48 libstdc++.so.6
-rwxrwxr-x 1 biocodz biocodz  849928 Oct 30 06:56 sortmerna

# set PATH
export PATH=$PWD/sortmerna-3.0.0/bin:$PATH

# test the installation
sortmerna --version
SortMeRNA version 3.0.0
Build Date: Oct 29 2018
sortmerna_build_git_sha:@bc47297ce3e8580ef32d161255bfcbe1f718c322@
sortmerna_build_git_date:@2018/10/29 12:35:39@

# view help
sortmerna -h
```

## Building from source code

### General notes

The build was tested using the following environments:
1. Linux
	* Ubuntu 14.04 Trusty with GCC 7.3.0
	* Ubuntu 16.04 Xenial with GCC 7.3.0
	* Ubuntu 16.04 Xenial with GCC 5.4.0 (default)
	* Centos 6.6 with GCC 7.3.0
2. Windows
	* Windows 10 with Visual Studio 15 2017 Win64
3. MAC
	* macOS 10.13 High Sierra (64-bit) with AppleClang 9.0.0.9000039

Getting latest GCC on _old_ Linux distros requires either installing GCC from PPA (Ubuntu), or building from sources - a lengthy process (around 10 hours on Centos VM running on VirtualBox Windows 10 host).

`CMake` is necessary for building. The distributions are available for all major operating systems. CMake can be easily built _if_ not available through a standard packager. Please visit [CMake project website](https://cmake.org/) for download and installation instructions.

The following libraries have to be installed using a packager or to be built
* `ZLib` (zlib1g-dev)
* `RocksDB` (librocksdb-dev)
* `RapidJson` (rapidjson-dev)
	
`Git` has to be installed _if_ building from the GitHub repository sources.

The following Flags can be used when generating the build files using CMake (`-D<FLAG>=VALUE`):
* `ROCKSDB_INCLUDE_DIR` (path to RocksDB include directory)
* `ROCKSDB_LIB_DEBUG` (path to RocksDB library for Debug)
* `ROCKSDB_LIB_RELEASE` (path to RocksDB library for Release)
* `ZLIB_LIB_DEBUG` (path to ZLib debug library location. Use if location is custom)
* `ZLIB_LIB_RELEASE` (path to ZLib release library locations. Use if location is custom)
* `SRC_ZLIB` (download Zlib sources. Use if ZLib is to be built)
* `SRC_ROCKSDB` (download RocksDB sources. Use if RocksDB is to be built)
* `SRC_RAPIDJSON` (download RapidJson sources. Use if 'apt install rapidjson' not available)
* `SET_ROCKSDB` (set to 1 to indicate RocksDB was built from sources. Not nesessary if RocksDB is installed using packager)
* `SET_ZLIB` (set to 1 to indicate ZLib was built from sources.)

The above flags can be ignored if the dependencies (`zlib`, `rocksdb`, `rapidjson`) are installed using a standard packager like `apt` (on Linux) or `homebrew` (on Mac)

General steps for building Sortmerna from the sources are as follows:
1. Prepare the build environment
2. Get the sources from GitHub or from GitHub releases
3. Build

The minute details on setting the build environment and performing the build can be found in [dnanexus folder](https://github.com/biocore/sortmerna/tree/master/dnanexus), which includes code for the build automation.

The build environment (step 1) for Ubuntu 16.04 is prepared using [sortmerna-3.asset](https://github.com/biocore/sortmerna/blob/master/dnanexus/assets/sortmerna-3.asset)
* The dependencies available as standard distros are prepared in [dxasset.json](https://github.com/biocore/sortmerna/blob/master/dnanexus/assets/sortmerna-3.asset/dxasset.json) ( see this [line](https://github.com/biocore/sortmerna/blob/master/dnanexus/assets/sortmerna-3.asset/dxasset.json#L10)).
* Non-standard packages (GCC-7 and CMake) are prepared in [Makefile](https://github.com/biocore/sortmerna/blob/master/dnanexus/assets/sortmerna-3.asset/Makefile)

The build (steps 2, 3) is performed using the applet [sortmerna-3.build.on.u16asset](https://github.com/biocore/sortmerna/tree/master/dnanexus/applets/sortmerna-3.build.on.u16asset)
* Git clone is called [here](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna-3.build.on.u16asset/src/build.sh#L37)
* CMake called to generate the build files [here](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna-3.build.on.u16asset/src/build.sh#L40)
* Make is called [here](https://github.com/biocore/sortmerna/blob/master/dnanexus/applets/sortmerna-3.build.on.u16asset/src/build.sh#L45)

### Obtaining the source code

The source code can be obtained using the following methods

1. Download GitHub Release
2. Clone the GitHub repository and use either master head branch (development) or checkout a release (tag)

#### Download GitHub Release

[GitHub releases](https://github.com/biocore/sortmerna/releases)

TODO

#### Clone the GitHub Repository and use the head or a release

```
# clone the repository
git clone https://github.com/biocore/sortmerna.git

pushd sortmerna

# If you need a particular release (tag)
git checkout v3.0.0
```

### building on Linux OS

(1) Install GCC if not already installed. SortmeRNA is C++14 compliant, so the GCC needs to be fairly new e.g. 5.4.0 works OK.

```bash
gcc --version
	gcc (Ubuntu 5.4.0-6ubuntu1~16.04.4) 5.4.0 20160609
```

(2) Install pre-requisites (CMake, Git, Zlib, RocksDB, RapidJson)

```
sudo apt update
sudo apt install cmake
sudo apt install git
sudo apt install zlib1g-dev librocksdb-dev rapidjson-dev
```
	
If the dependencies cannot be installed using a package manager, they need to be built (read below).
	
(3) Clone the Git repository

```
git clone https://github.com/biocore/sortmerna.git
```
	
(4) Generate the build files using CMake:

```bash
mkdir -p $SMR_HOME/build/Release
pushd $SMR_HOME/build/Release
```

(4.1) If all the dependencies are available on the system

```bash
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release ../..
```
	
(4.2) If RocksDB and RapidJson have to be installed from sources (see the flags description above)
	
```bash
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DSRC_ROCKSDB=1 -DSRC_RAPIDJSON=1 -DSET_ROCKSDB=1 ../..
```

The above will download RocksDB and RapidJson into default locations ($SMR_HOME/3rdparty/rocksdb) and ($SMR_HOME/3rdparty/rapidjson) correspondingly.

OR with custom values for RocksDB include/lib
	
```bash
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DSRC_ROCKSDB=1 -DSRC_RAPIDJSON=1 -DSET_ROCKSDB=1 -DROCKSDB_INCLUDE_DIR=$SOME_DIR/rocksdb/include -DROCKSDB_LIB_RELEASE=$SOME_DIR/rocksdb/build/Release ../..
```

NOTE: `$SMR_HOME` is the top directory where sortmerna code (e.g. git repo) is located.

Other compiler/linker flags that might be necessary depending on the system:

* -DEXTRA_CXX_FLAGS_RELEASE="-lrt" (had to use this on Centos 6.6 + GCC 7.3.0)

The above commands will perform necessary system check-ups, dependencies, and generate Makefile.

(5) Compile and build executables:

(5.1) If RocksDB needs to be built

```bash
mdir -p SMR_HOME/3rdparty/rocksdb/build/Release
pushd SMR_HOME/3rdparty/rocksdb/build/Release
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DPORTABLE=1 -DWITH_ZLIB=1 -DWITH_TESTS=0 -DWITH_TOOLS=0 ../..
make
popd
```

(5.2)

Build SortMeRNA

```bash
make
```

The binaries are created in `$SMR_HOME/build/Release/src/indexdb` and `$SMR_HOME/build/Release/src/sortmerna`

```
# add the build binaries to the PATH
export PATH="$SMR_HOME/build/Release/src/indexdb:$SMR_HOME/build/Release/src/sortmerna:$PATH"
```

### Building on Mac OS

We tested the build on macOS 10.13 High Sierra (64-bit).
We recommend the Homebrew - an excellent packager for Mac [1], which has all the latest packages required to build SortmeRNA.
The build can be performed using either Clang or GCC.

(1) Install Homebrew:

```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" [1]

brew --version
brew help
```
	
(2) Install pre-requisites (CMake, Git, Zlib, RocksDB, RapidJson)

```bash
brew install cmake
brew install git
brew install zlib
brew install rocksdb
brew install rapidjson
```

(3) Clone the GIt repository

```
git clone https://github.com/biocore/sortmerna.git
```

(4) Generate the build files:

```bash
mkdir -p $SMR_HOME/build/Release
pushd $SMR_HOME/build/Release
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DEXTRA_CXX_FLAGS_RELEASE="-pthread" ../..
	-- The CXX compiler identification is AppleClang 9.0.0.9000039
	-- The C compiler identification is AppleClang 9.0.0.9000039
	-- Check for working CXX compiler: /Library/Developer/CommandLineTools/usr/bin/c++
	-- Check for working CXX compiler: /Library/Developer/CommandLineTools/usr/bin/c++ -- works
	-- Detecting CXX compiler ABI info
	-- Detecting CXX compiler ABI info - done
	-- Detecting CXX compile features
	-- Detecting CXX compile features - done
	-- Check for working C compiler: /Library/Developer/CommandLineTools/usr/bin/cc
	-- Check for working C compiler: /Library/Developer/CommandLineTools/usr/bin/cc -- works
	-- Detecting C compiler ABI info
	-- Detecting C compiler ABI info - done
	-- Detecting C compile features
	-- Detecting C compile features - done
	CMAKE_CXX_COMPILER_ID = AppleClang
	CMAKE_CONFIGURATION_TYPES =
	CMAKE_CXX_FLAGS_RELEASE: -O3 -DNDEBUG
	EXTRA_CXX_FLAGS_RELEASE: -pthread
	Cloning into 'concurrentqueue'...
	Checking out files: 100% (1613/1613), done.
	-- Configuring done
	-- Generating done
	-- Build files have been written to: /Users/bc/sortmerna/build/Release
```

Note: `$SMR_HOME` is the top directory where sortmerna code (e.g. git repo) is located.

CMake will perform necessary system check-ups, dependencies, and generate Makefile.

(5) Compile and build executables:

```bash
make
```

The binaries are created in `$SMR_HOME/build/Release/src/indexdb` and `$SMR_HOME/build/Release/src/sortmerna`

Simply add the build binaries to the PATH e.g.

```
export PATH="$SMR_HOME/build/Release/src/indexdb:$SMR_HOME/build/Release/src/sortmerna:$PATH"
```

#### Install Clang on Mac

Installing Xcode (free through the App Store) and Xcode command line tools will automatically 
install the latest version of Clang supported with Xcode. 

After installing Xcode, the Xcode command line tools may be installed via:

Xcode -> Preferences -> Downloads

Under "Components", click to install "Command Line Tools"

#### Configure shell to use Clang compiler on Mac

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

#### Configure shell to use GCC compiler on Mac

(1) Check if you have GCC installed:

```bash
gcc --version
```

(2a) If GCC is installed, set your compiler to GCC:

```bash
export CC=gcc-mp-5.4
export CXX=g++-mp-5.4
```

(2b) If GCC is not installed, it can be installed through Homebrew or MacPorts.

```
brew tap homebrew/versions
brew install [flags] gcc54
```

To list available flags
```
brew options gcc54
```

### Building on Windows OS

MS Visual Studio Community edition and CMake for Windows are required for building SortMeRNA.

We tested the build using `Visual Studio 15 2017 Win64` and `Visual Studio 14 2015 Win64`

(1) Download and Install VS Community edition from [Visual Studio community website](https://www.visualstudio.com/vs/community/)

(2) Install CMake

CMake can be installed using either Windows Installer or binaries from archive.
Download binary distributions from [here](https://cmake.org/download/)

If you choose portable binaries (not the installer) e.g. cmake-3.11.0-rc1-win64-x64.zip,
just download and extract the archive in a directory of your choice e.g.

```
C:\libs\cmake-3.11.0-rc1-win64-x64\
	bin\
	doc\
	man\
	share\
```

The `bin` directory above contains `cmake.exe` and `cmake-gui.exe`. Add the `bin` directory to your PATH
Start `cmd` and

```
set PATH=C:\libs\cmake-3.11.0-rc1-win64-x64\bin;%PATH%
cmake --version
cmake-gui
```

(3) Install Git for Windows

Download binary distribution either portable or the installer from [here](https://git-scm.com/download/win)

The portable distribution is a self-extracting archive that can be installed in a directory of your choice e.g.

```
C:\libs\git-2.16.2-64\
	bin\
	cmd\
	dev\
	etc\
	mingw64\
	tmp\
	usr\
```
You can use either `bash.exe` or native Windows CMD `cmd.exe`.

If you choose to work with CMD, add the following to your path:

```
set GIT_HOME=C:\libs\git-2.16.2-64
set PATH=%GIT_HOME%\bin;%GIT_HOME%\usr\bin;%GIT_HOME%\mingw64\bin;%PATH%

git --version
```

(4) Clone the GIt repository

```
git clone https://github.com/biocore/sortmerna.git
```

(5) Prepare the build files:

On Windows we recommend using the `cmake-gui` utility.
Either navigate to CMake installation directory (using Windows Explorer) and double-click
`cmake-gui`, or launch it from command line as shown below:

```
set PATH=C:\libs\cmake-3.11.0-rc1-win64-x64\bin;%PATH%
cmake-gui
```

In the CMake GUI
 - click `Browse source` button and navigate to the directory where Sortmerna sources are located (SMR_HOME).
 - click `Browse Build` and navigate to the directory where to build the binaries e.g. %SMR_HOME%\build
 - at the prompt select the Generator from the list e.g. "Visual Studio 15 2017 Win64"
 - click `Configure`
 - Set the following variables:
   - ZLIB_INCLUDE_DIR=%SMR_HOME%/3rdparty/zlib
   - ZLIB_LIB_DEBUG=%SMR_HOME%/3rdparty/zlib/build/Debug
   - ZLIB_LIB_RELEASE=%SMR_HOME%/3rdparty/zlib/build/Release
   - ROCKSDB_INCLUDE_DIR=%SMR_HOME%/3rdparty/rocksdb/include
   - ROCKSDB_LIB_DEBUG=%SMR_HOME%/3rdparty/rocksdb/build/Debug
   - ROCKSDB_LIB_RELEASE=%SMR_HOME%/3rdparty/rocksdb/build/Release
 - click `Configure` again
 - click `Generate` if all variables were set OK (no red background)

The `Generate` generates VS project files in `%SMR_HOME%\build\` directory.
`%SMR_HOME%` is the top directory where SortMeRNA source distribution (e.g. Git repo) is installed.

(6) Configure and build Zlib library

When Cmake-gui `Configure` is run it downloads required 3rd party source packages into `%SMR_HOME%\3rdparty\` directory.

In Cmake-gui:
- click `Browse Source...` and select `%SMR_HOME%\3rdparty\zlib\`
- click `Browse Build...` and select `%SMR_HOME%\3rdparty\zlib\build\` (confirm to create the `build` directory if not already exists)
- click `Configure` and set the required variables or accept defaults
- click `Generate`

In Visual Studio
- `File -> Open -> Project/Solution` and select `%SMR_HOME%\3rdparty\zlib\build\zlib.sln`
- In Solution Explorer right-click `ALL_BUILD` and select `build` from drop-down menu

(7) COnfigure and build RockDB library

In Cmake-gui:
- click `Browse Source...` and select `%SMR_HOME%\3rdparty\rocksdb\`
- click `Browse Build...` and select `%SMR_HOME%\3rdparty\rocksdb\build\` (confirm to create the `build` directory if not already exists)
- click `Configure` and set the following variables:
  - Ungrouped Entries
    - PORTABLE (checkbox)
	- GIT_EXECUTABLE (select path to `git.exe` e.g. `C:/libs/git-2.16.2-64/bin/git.exe`
  - WITH
    - WITH_MD_LIBRARY
	- WITH_ZLIB
	- Accept defaults for the rest
- click `Generate`

In Visual Studio
- `File -> Open -> Project/Solution` and select `%SMR_HOME%\3rdparty\rocksdb\build\rocksdb.sln`
- In Solution Explorer right-click `ALL_BUILD` and select `build` from drop-down menu

(8) Build SormeRNA

In Visual Studio:
- `File -> Open -> Project/Solution .. open %SMR_HOME%\build\sortmerna.sln`
- Select desired build type: `Release | Debug | RelWithDebInfo | MinSizeRel`.
- In Solution explorer right-click `ALL_BUILD' and select `build` in drop-down menu.

Depending on the build type the binaries are generated in 
`%SMR_HOME%\build\src\sortmerna\Release` (or `Debug | RelWithDebInfo | MinSizeRel`).

(9) Add sortmerna executables to PATH

```
set PATH=%SMR_HOME%\build\src\indexdb\Release;%SMR_HOME%\build\src\sortmerna\Release;%PATH%
```

# Running

Python code is provided for running integration tests in $SRM_HOME/tests (%SRM_HOME%\tests) and requires Python 3.5 or higher.

Tests can be run with the following command:

```
python ./tests/test_sortmerna.py
```

OR individual tests

```
python ./tests/test_sortmerna.py SortmernaTests.test_simulated_amplicon_generic_buffer
```

Tests on compressed data files

```
python ./tests/test_sortmerna_zlib.py
```

Users require [scikit-bio](https://github.com/biocore/scikit-bio) 0.5.0 to run the tests.

# User Manual

User manual is available in [docs/web folder](https://github.com/biocore/sortmerna/tree/master/docs/web).
The manual is written as a single web page using HTML, CSS, and JS (very minimal for changing color theme light/dark), so cannot be viewed directly on GitHub.
Clone the repository to your local machine and open [index.html](https://github.com/biocore/sortmerna/blob/master/docs/web/index.html) in the web browser.
In case you prefer PDF, any decent browser can print web pages to PDF.
Please, note the manual was tested so far only using Chrome on FHD display (1920 x 1080).

# Third-party libraries

Various features in SortMeRNA are dependent on third-party libraries, including:
* [ALP](http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/software/program.html?uid=6): computes statistical parameters for Gumbel distribution (K and Lambda)
* [CMPH](http://cmph.sourceforge.net): C Minimal Perfect Hashing Library
* [Zlib](https://github.com/madler/zlib): reading compressed Reads files
* [RocksDB](https://github.com/facebook/rocksdb): storage for SortmeRNA alignment results
* [RapidJson](https://github.com/Tencent/rapidjson): serialization of Reads objects to store in RocksDB
* [Concurrent Queue](https://github.com/cameron314/concurrentqueue): Lockless buffer for Reads accessed from multiple processing threads

# Taxonomies

The folder `rRNA_databases/silva_ids_acc_tax.tar.gz` contains SILVA taxonomy strings (extracted from XML file generated by ARB)
for each of the reference sequences in the representative databases. The format of the files is three tab-separated columns,
the first being the reference sequence ID, the second being the accession number and the final column is the taxonomy.

# Citation

If you use SortMeRNA, please cite:
Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611.

# Contributors

See [AUTHORS](./AUTHORS) for a list of contributors to this project.

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

This command will generate a folder `html` in the directory from which the
command was run.

# References

1. Homebrew 
	- [home](https://brew.sh/)
	- [github](https://github.com/Homebrew)
	
