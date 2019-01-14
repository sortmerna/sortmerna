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
	* [Using GitHub release binaries](#using-github-release-binaries)
	* [Building from source code](#building-from-source-code)
		* [General notes](#general-notes)
		* [Obtaining the source code](#obtaining-the-source-code)
			* [Download GitHub Release](#download-github-release)
			* [Clone the GitHub Repository and use the head or a release](#clone-the-github-repository-and-use-the-head-or-a-release)
		* [building on Linux OS](#building-on-linux-os)
		* [Building on MacOS](#building-on-macos)
		* [Building on Windows OS](#building-on-windows-os)
	* [DNANexus cloud](#dnanexus-cloud)
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

SortMeRNA can be run/built on Windows, Linux, and Mac.

The following methods can be used:

## Using GitHub release binaries

Visit [Sortmerna GitHub Releases](https://github.com/biocore/sortmerna/releases)

Any of the following 2 files can be used for installation:
- `sortmerna-<VERSION>-linux.sh` (convenience script that embeds `sortmerna-<VERSION>-linux.tar.gz`)
- `sortmerna-<VERSION>-linux.tar.gz`

Installation commands in `bash`:

```
# get binary distro (substitute the release version as desired)
wget https://github.com/biocore/sortmerna/releases/download/v3.0.3/sortmerna-3.0.3-linux.sh

# view the installer usage
chmod +x ./sortmerna-3.0.3-linux.sh
./sortmerna-3.0.3-linux.sh --help
Usage: dist/sortmerna-3.0.3-Linux.sh [options]
Options: [defaults in brackets after descriptions]
  --help            print this message
  --version         print cmake installer version
  --prefix=dir      directory in which to install
  --include-subdir  include the sortmerna-3.0.3-Linux subdirectory
  --exclude-subdir  exclude the sortmerna-3.0.3-Linux subdirectory
  --skip-license    accept license

./sortmerna-3.0.3-Linux.sh --prefix=$HOME/sortmerna --exclude-subdir --skip-license
sortmerna Installer Version: 3.0.3, Copyright (c) Clarity Genomics
This is a self-extracting archive.
The archive will be extracted to: /home/biocodz/sortmerna

Using target directory: /home/biocodz/sortmerna
Extracting, please wait...

Unpacking finished successfully

ls -lrt /home/biocodz/sortmerna/bin/
indexdb
sortmerna

# set PATH
export PATH=$HOME/sortmerna/bin:$PATH

# test the installation
sortmerna --version
SortMeRNA version 3.0.3
Build Date: Dec 11 2019
sortmerna_build_git_sha:@bc47297ce3e8580ef32d161255bfcbe1f718c322@
sortmerna_build_git_date:@2019/01/11 12:35:39@

# view help
sortmerna -h
```

## Building from source code

### General notes

General steps for building Sortmerna from the sources are as follows:
1. Prepare the build environment
2. Get the sources from GitHub or from GitHub releases
3. Build

`HINT`: The minute details on setting the build environment and performing the build can be found in [travis.yml](https://github.com/biocore/sortmerna/blob/master/.travis.yml), which is a [Travis CI](https://travis-ci.org) configuration file used for Sortmerna automated builds.

The build was tested using the following environments:
1. Linux
	* Ubuntu 18.04 Bionic with GCC 7.3.0 or Clang 7.0.0
	* Ubuntu 16.04 Xenial with GCC 7.3.0 or Cland 7.0.0
	* Ubuntu 16.04 Xenial with GCC 5.4.0 (default)
	* Ubuntu 14.04 Trusty with GCC 7.3.0 or Cland 7.0.0
	* Centos 6.6 with GCC 7.3.0
2. Windows
	* Windows 10 with Visual Studio 15 2017 Win64
3. MacOS
	* macOS 10.13 High Sierra (64-bit) with AppleClang 9.0.0.9000039
	* macOS 10.12.6 (64-bit) with clang-900.0.39.2

Getting latest GCC on _old_ Linux distros requires either installing GCC from PPA (Ubuntu), or building from sources - a lengthy process (around 10 hours on Centos VM running on VirtualBox Windows 10 host).

`CMake-3.13` of higher is necessary for building. The distributions are available for all major operating systems. CMake can be easily built _if_ not available through a standard packager. Please visit [CMake project website](https://cmake.org/) for download and installation instructions.

The following libraries have to be installed using a packager or to be built
* `ZLib` (zlib1g-dev)
* `RocksDB` (librocksdb-dev)
* `RapidJson` (rapidjson-dev)
	
`Git` has to be installed _if_ building from the GitHub repository sources.

The following Flags can be used when generating the build files using CMake (`-D<FLAG>=VALUE`):
* `CMAKE_INSTALL_PREFIX` (directory path to use for installation)
* `ZLIB_ROOT` (ZLIB installation directory)
* `ZLIB_LIBRARY_DEBUG` (path to ZLib debug library location. Use if location is custom)
* `ZLIB_LIBRARY_RELEASE` (path to ZLib release library locations. Use if location is custom)
* `ROCKSDB_SRC` (RocksDB source root directory)
* `ROCKSDB_HOME` (RocksDB installation directory)
* `RAPIDJSON_HOME` (RapidJSON installation directory)

The above flags can be ignored if the dependencies (`zlib`, `rocksdb`, `rapidjson`) are installed using a standard packager like `apt` (on Linux) or `homebrew` (on Mac)


### Obtaining the source code

The source code can be obtained using `any` of the following:

1. Download archived sources from GitHub Releases
2. Clone the GitHub repository and use either master head branch (development) or checkout a release (tag)

#### Download GitHub Release

Visit [GitHub releases](https://github.com/biocore/sortmerna/releases) and click `Source code` link to download the archive, or alternatively on command line:

```
wget https://github.com/biocore/sortmerna/archive/v3.0.3.tar.gz
```

#### Clone the GitHub Repository and use the HEAD or a release

```
# clone the repository
git clone https://github.com/biocore/sortmerna.git

pushd sortmerna

# If you need a particular release (tag)
git checkout v3.0.3
```

### building on Linux OS

(1) Install GCC or/and LLVM Clang if not already installed. SortmeRNA is C++14 compliant, so the compiler needs to be fairly new e.g. GCC 5.4.0 works OK.

```bash
gcc --version
	gcc (Ubuntu 5.4.0-6ubuntu1~16.04.4) 5.4.0 20160609

# for Clang
cc --version
	clang version 7.0.1-svn348686-1~exp1~20181211132844.50 (branches/release_70)
	Target: x86_64-pc-linux-gnu
	Thread model: posix
	InstalledDir: /usr/bin
```

(2) Install pre-requisites (CMake, Git, Zlib, RocksDB, RapidJson)

```
sudo apt update
sudo apt install cmake
sudo apt install git
sudo apt install zlib1g-dev librocksdb-dev rapidjson-dev
```
	
If the dependencies cannot be installed using a package manager, they need to be built (read below).

(3) Build Zlib (fast - under a minute)
```
git clone https://github.com/madler/zlib.git
mkdir -p zlib/build/Release
pushd zlib/build/Release
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../../dist ../.. 
cmake --build . 
cmake --build . --target install
ls -lrt ../../dist/lib/ 
popd
```

(4) Build RocksDB (may take ~20min)

```bash
git clone https://github.com/facebook/rocksdb.git;
pushd rocksdb

# checkout the latest release
ROCKSDB_RELEASE=v5.17.2; git checkout tags/${ROCKSDB_RELEASE}

# modify thirdparty.inc to point to correct ZLIB installation. 
# For example:
vi thirdparty.inc
...
set(ZLIB_HOME ${HOME}/zlib/dist)
set(ZLIB_INCLUDE ${ZLIB_HOME}/include)
set(ZLIB_LIB_DEBUG ${ZLIB_HOME}/lib/libz.a)
set(ZLIB_LIB_RELEASE ${ZLIB_HOME}/lib/libz.a)
...

mkdir -p build/Release
pushd build/Release;
cmake -G "Unix Makefiles" \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX=../../dist \
    -DWITH_ZLIB=1 \
	-DWITH_GFLAGS=0 \
	-DPORTABLE=1 \
	-DWITH_TESTS=0 \
	-DWITH_TOOLS=0 ../..

cmake --build . ;
cmake --build . --target install;
popd;
# check the distro
ls -l rocksdb/dist
	drwxrwxr-x 3 biocodz biocodz 4096 Dec 27 13:58 include
	drwxrwxr-x 3 biocodz biocodz 4096 Dec 27 13:58 lib

ls -l dist/lib
	/home/biocodz/rocksdb/dist/lib/librocksdb.a
	/home/biocodz/rocksdb/dist/lib/librocksdb.so.5.18.0
popd
```

(5) Build RapidJson (fast - under a minute)

```
git clone https://github.com/Tencent/rapidjson.git
mkdir -p rapidjson/build/Release
pushd rapidjson/build/Release

cmake -G "Unix Makefiles" \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX=../../dist \
    -DRAPIDJSON_BUILD_EXAMPLES=OFF \
	-DRAPIDJSON_BUILD_DOC=OFF ../..

cmake --build .
cmake --build . --target install
popd
# check the distro (only 'include/' is necessary)
ls -l ~/rapidjson/dist/
	drwxrwxr-x 3 biocodz biocodz 4096 Dec 27 14:52 include
	drwxrwxr-x 4 biocodz biocodz 4096 Dec 27 14:52 lib
	drwxrwxr-x 3 biocodz biocodz 4096 Dec 27 14:52 share
```
	
(6) Build Sortmerna

```
git clone https://github.com/biocore/sortmerna.git

mkdir -p $HOME/sortmerna/build/Release
pushd $HOME/sortmerna/build/Release

# If alternative toolchains (compilers) are installed on the system

# to use Clang
rm CMakeCache.txt
cmake -E env CC=/usr/bin/clang-7 CXX=/usr/bin/clang++-7 \
cmake -G "Unix Makefiles" \
	-DCMAKE_BUILD_TYPE=Release \
	... # other flags as necessary (see above)
	../..

# to use GCC
rm CMakeCache.txt
CC=gcc-7 CXX=g++-7 cmake \
cmake -G "Unix Makefiles" \
	-DCMAKE_BUILD_TYPE=Release \
	... # other flags as necessary (see above)
	../..

# If all the dependencies are available on the system
cmake -G "Unix Makefiles" \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX=../../dist \
	-DCPACK_BINARY_TGZ=ON ../..

##				OR
#  If ZLib, RocksDB, and RapidJson were built from sources (see above)

cmake -G "Unix Makefiles" \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX=../../dist \
	-DCPACK_BINARY_TGZ=ON \
	-DZLIB_ROOT=~/zlib/dist \
	-DZLIB_LIBRARY_RELEASE=~/zlib/dist/lib/libz.a \
	-DZLIB_LIBRARY_DEBUG=~/zlib/dist/lib/libzd.a \
	-DROCKSDB_HOME=~/rocksdb/dist \
	-DROCKSDB_SRC=~/rocksdb \
	-DRAPIDJSON_HOME=~/rapidjson/dist ../..

cmake --build .
cmake --build . --target install
cmake --build . --target package

popd

# check the distro
ls -l ~/sortmerna/dist
	drwxrwxr-x 2 biocodz biocodz    4096 Jan  2 18:09 bin
	-rwxrwxrwx 1 biocodz biocodz 2074449 Jan  2 18:09 sortmerna-3.0.3-Linux.sh
	-rw-rw-r-- 1 biocodz biocodz 2070507 Jan  2 18:09 sortmerna-3.0.3-Linux.tar.gz
	-rw-rw-r-- 1 biocodz biocodz 2877845 Jan  2 18:09 sortmerna-3.0.3-Linux.tar.Z

# add the build binaries to the PATH
export PATH="$HOME/sortmerna/dist/bin/indexdb:$HOME/sortmerna/dist/bin/sortmerna:$PATH"

# test
sortmerna --version
	SortMeRNA version 3.0.3
	Build Date: Jan  9 2019
	sortmerna_build_git_sha:@2a619c87cdcb6e98c5460630182f161359cfd1d1@
	sortmerna_build_git_date:@2019/01/09 18:04:04@

sortmerna -h
```

Other compiler/linker flags that might be necessary depending on the system:

* -DEXTRA_CXX_FLAGS_RELEASE="-lrt" (had to use this on Centos 6.6 + GCC 7.3.0)


### Building on MacOS

The build can be performed using either `Clang` or `GCC`.

For third-party software installation we recommend the Homebrew - an excellent packager for Mac [1], which has all the latest packages. 

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

(3) Install `Clang` on Mac (if necessary)

Installing Xcode (free through the App Store) and Xcode command line tools will automatically 
install the latest version of Clang supported with Xcode. 

After installing Xcode, the Xcode command line tools may be installed via:

Xcode -> Preferences -> Downloads

Under "Components", click to install "Command Line Tools"

(4) Install `GCC` (if necessary) using Homebrew.

```
# list available flags
brew options gcc7
brew tap homebrew/versions
brew install [flags] gcc7
```

The rest is the same as for Linux (see above starting with step 3)
For older versions of OSX you may need to use `-DEXTRA_CXX_FLAGS_RELEASE="-pthread"`.


### Building on Windows OS

MS Visual Studio Community edition and CMake for Windows are required for building SortMeRNA.

(1) Download and Install VS Community edition from [Visual Studio community website](https://www.visualstudio.com/vs/community/)

(2) Install CMake

CMake can be installed using either Windows Installer or binaries from archive.
Download binary distributions from [here](https://cmake.org/download/)

If you choose `portable` binaries (not the installer) e.g. `cmake-3.13.1-win64-x644.zip`,
just download and extract the archive in a directory of your choice e.g.

```
C:\libs\cmake-3.13.1-win64-x64\
	bin\
	doc\
	man\
	share\
```

The `bin` directory above contains `cmake.exe` and `cmake-gui.exe`. Add the `bin` directory to your PATH
Start `cmd` and

```
set PATH=C:\libs\cmake-3.13.1-win64-x64\bin;%PATH%
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
You can use either `bash.exe` provided with the Git distro, or native Windows CMD `cmd.exe`.

If you choose to work with `cmd.exe`, add the following to your path:

```
set GIT_HOME=C:\libs\git-2.16.2-64
set PATH=%GIT_HOME%\bin;%GIT_HOME%\usr\bin;%GIT_HOME%\mingw64\bin;%PATH%

git --version
```

(4) Configure and build Zlib library

Because we use CMake as the build tool, the steps are the same a for Linux (see Linux instructions above).
For the toolchain specify `cmake -G "Visual Studio 15 2017 Win64"`.
No need to specify `-DCMAKE_BUILD_TYPE` because Visual studio is a multi-type configuration in CMake parlant.

```
git clone https://github.com/madler/zlib.git
```

If using command line the steps are the same as on Linux (see above)

If using CMake GUI:
- click `Browse Source...` and select `%ZLIB_HOME%` i.e. the root directory where the ZLib sources are.
- click `Browse Build...` and select `%ZLIB_HOME%\build\` (confirm to create the `build` directory if not already exists)
- click `Configure` and set the required variables or accept defaults
- click `Generate`

In Visual Studio
- `File -> Open -> Project/Solution` and select `%ZLIB_HOME%\build\zlib.sln`
- In Solution Explorer right-click `ALL_BUILD` and select `build` from drop-down menu

(5) Configure and build RockDB library

See instructions for Linux if using command line.

If using CMake GUI:
- click `Browse Source...` and select `%ROCKSDB_SRC%"`
- click `Browse Build...` and select `%ROCKSDB_SRC%\build\` (confirm to create the `build` directory if not already exists)
- click `Configure` and set the following variables:
  - Ungrouped Entries
    - PORTABLE (check)
	- GIT_EXECUTABLE (select path to `git.exe` e.g. `C:/libs/git-2.16.2-64/bin/git.exe`
  - WITH
    - WITH_MD_LIBRARY=1
	- WITH_ZLIB=1
	- WITH_TESTS=0
	- WITH_TOOLS=0
	- Accept defaults for the rest
  - CMAKE
	- CMAKE_INSTALL_PREFIX=../../dist
- click `Generate`

In Visual Studio
- `File -> Open -> Project/Solution` and select `%ROCKSDB_SRC%\build\rocksdb.sln`
- In Solution Explorer right-click `ALL_BUILD` and select `build` from drop-down menu

(6) Configure and build RapidJson

See instructions for Linux if using command line.

(7) Configure and build Dirent for Windows [2]

This is optional. If CMake won't find the Dirent when configuring Sortemrna, it will download it from GitHub [2] into `sortmerna/3rdparty/` folder. Only a single header file is used, so no build is required.

(8) Configure and build Sortmerna

Using standard `cmd` Windows shell:

`SMR_HOME` is the top directory where SortMeRNA source distribution (e.g. Git repo) is installed.

```
git clone https://github.com/biocore/sortmerna.git

set SMR_HOME=C:/projects/sortmerna 

mkdir %SMR_HOME%/build
pushd %SMR_HOME%/build

cmake -G "Visual Studio 15 2017 Win64" \
	-DCMAKE_INSTALL_PREFIX=%SMR_HOME%/dist \
	-DCPACK_BINARY_7Z=ON \
	-DCPACK_SOURCE_7Z=ON \
	-DCPACK_SOURCE_ZIP=OFF \
	-DWITH_MD_LIBRARY=ON \
	-DZLIB_ROOT="%LIBS_HOME%/zlib/dist" \
	-DZLIB_LIBRARY_RELEASE="%ZLIB_HOME%/dist/lib/zlibstatic.lib" \
	-DZLIB_LIBRARY_DEBUG="%ZLIB_HOME%/dist/lib/zlibstaticd.lib" \
	-DROCKSDB_SRC="%ROCKSDB_SRC%" \
	-DROCKSDB_HOME="%ROCKSDB_SRC%/dist" \
	-DRAPIDJSON_HOME="%RAPIDJSON_SRC/dist" \
	-DDIRENTWIN_HOME=%DIRENTWIN_SRC%/dist ..

# possible builds: Debug | Release | RelWithDebInfo | MinSizeRel
# Debug - default
cmake --build .
#     OR Release
cmake --build . --config Release

# generate distro
cmake --build . --target install
cmake --build . --target package

# test
set PATH=%SMR_HOME%\dist\bin;%PATH%
sortmerna --version
```

If using CMake GUI

To launch `CMake GUI` navigate to CMake installation directory (using Windows Explorer) and double-click
`cmake-gui`, or launch it from command line:

```
set PATH=C:\libs\cmake-3.13.1-win64-x64\bin;%PATH%
cmake-gui
```

In GUI:
 - click `Browse source` button and navigate to the directory where Sortmerna sources are located (SMR_HOME).
 - click `Browse Build` and navigate to the directory where to build the binaries e.g. `%SMR_HOME%\build`
 - at the prompt select the Generator from the list e.g. "Visual Studio 15 2017 Win64"
 - click `Configure`

 - Set the same variables as shown above for command line i.e.:
	- CMAKE_INSTALL_PREFIX
	- CPACK_BINARY_7Z
	- CPACK_SOURCE_7Z
	- CPACK_SOURCE_ZIP
	- WITH_MD_LIBRARY
	- ZLIB_ROOT
	- ZLIB_LIBRARY_RELEASE
	- ZLIB_LIBRARY_DEBUG
	- ROCKSDB_SRC
	- ROCKSDB_HOME
	- RAPIDJSON_HOME
	- DIRENTWIN_HOME

 - click `Configure` again
 - click `Generate` if all variables are OK (no red background)

The `Generate` generates VS project files in `%SMR_HOME%\build\` directory.
Start Visual Studio and:
- `File -> Open -> Project/Solution .. open %SMR_HOME%\build\sortmerna.sln`
- Select desired build type: `Release | Debug | RelWithDebInfo | MinSizeRel`.
- In Solution explorer right-click `ALL_BUILD' and select `build` in drop-down menu.
- Execute target 'INSTALL'
- Execute target 'PACKAGE'

## DNANexus cloud

So far (since 20181117) we were not able to make Sortmerna application public due DNANexus lack of interest or resources.
If you need the access to the application, please, contact `biocodz at protonmail dot com` to be added to the list of authorized users.

* Ready-to-use `sortmerna` application is available on DNANexus cloud. Just upload your data and run the application.
* Sortmerna distribution also contains DNANexus applets for building, testing and running Sortmerna. See [dnanexus](https://github.com/biocore/sortmerna/blob/master/dnanexus/README.md) subdirectory for detailed instructions.

Note that DNANexus trial accounts can be used.

# Running

Python code is provided for running integration tests in $SRM_HOME/tests (%SRM_HOME%\tests) and requires Python 3.

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

Users require [scikit-bio](https://github.com/biocore/scikit-bio) to run the tests.

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
2. [Dirent for Windows](https://github.com/tronkko/dirent.git)
	
