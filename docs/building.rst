Building
========

General Notes
-------------

General build steps:

1. Prepare the build environment
2. Get the sources from GitHub or from GitHub releases
3. Build

SortmeRNA-4 is C++17 compliant.

The build uses CMake controlled from a Python script provided with the Sortmerna distribution.

Currently the builds are performed on Linux with GCC, on MacOS with AppleCland, and Windows with native windows compilers (available with Visual Studio Community edition).

Build instructions are detailed in this github workflow: https://github.com/sortmerna/sortmerna/blob/master/.github/workflows/build_multi_platform.yml


Building on Linux
-----------------

Quick Build
###########

The necessary build pre-requisits are::

   GCC (>=9, <=11.4.0)
   Conda
   CMake

If GCC and Conda are available, the following commands will build SortMeRNA::

   # create and activate conda environment
   conda create -n sortmerna -c conda-forge pyyaml jinja2 requests cmake
   conda activate sortmerna
   # clone git repo
   git clone https://github.com/sortmerna/sortmerna.git
   pushd sortmerna
   # build
   python setup.py -n all
   # test
   dist/bin/sortmerna -h

The above builds all required dependencies as listed in 'sortmerna/3rdparty.jinja'. Those include RocksDB and ZLib libraries.
The built artifacts are output into 'sortmerna/dist' directory

Install GCC
###########

For Ubuntu 20.04 and later::

   # check the OS release
   lsb_release -a
      Ubuntu 22.04.3 LTS

   # check versions installed/available
   apt policy gcc
      Installed: 4:11.2.0-1ubuntu1
      Candidate: 4:11.2.0-1ubuntu1

   # install if necessary
   sudo apt install gcc

For older Debian distros (Ubuntu) GCC can be installed from PPA::

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

Install Conda
#############

Using official installer::
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh

Get SortMeRNA sources
#####################

The sources can be placed anywhere, but here the user's Home directory is assumed::

   # clone the repository
   git clone https://github.com/biocore/sortmerna.git
   
   pushd sortmerna
   
   # If you need a particular release (tag)
   git checkout v4.3.7
   
   # alternatively get the release sources
   wget https://github.com/biocore/sortmerna/archive/v4.3.7.tar.gz
   
   tar xzf v4.3.7.tar.gz
