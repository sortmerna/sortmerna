# Building from source code

# Contents

* [General notes](#general-notes)
* [Building on Linux](#building-on-linux)
	* [Install GCC 9](#install-gcc-9)
	* [Install Devtoolset 9](#install-devtoolset-9)
	* [Install Conda](#install-conda)
	* [get Sortmerna sources](#get-sortmerna-sources)
	* [build](#build)
* [Building on Windows](#building-on-windows)
	* [Install Visual Studio](#install-visual-studio)
	* [Install Conda](#install-conda)
	* [get Sortmerna sources](#get-sortmerna-sources)
	* [build](#build)
* [building on Max](#building-on-mac)
* [Third-party libraries](#third-party-libraries)

## General notes

General build steps:
1. Prepare the build environment
2. Get the sources from GitHub or from GitHub releases
3. Build

SortmeRNA-4 is C++17 compliant, and requires a compiler that supports `filesystem` standard library.
Currently only GCC-9 and Windows SDK meet that requirement.

Ready to use GCC-9 distribution is only available on Debian (Ubuntu). It is not available on Centos, i.e. no Devtoolset-9 yet, and building GCC-9 on Centos is not trivial.

Clang has no support for the `filesystem` yet. LLVM 9.0 is due to be released later this year. At that time we'll get back to supporting builds with Clang on Linux and OSX.

The build is performed using a Python script provided with the Sortmerna distribution.
The script uses a configuration file `env.yaml`, which can be Optionally modified to customize the build.

## Building on Linux

### Install GCC 9

This is for Debian distros (Ubuntu)

```
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
```

### Install Devtoolset 9

TODO. 

This section is for RHEL (Centos) distros. Waiting for Devtoolset-9 to become available.

### Install Conda

The following will install Conda in the User's Home directory.

```
pushd $HOME
wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
export PATH=$HOME/miniconda3/bin:$PATH
pip install -U pip
pip install pyyaml
pip install jinja2
```

### get Sortmerna sources

The sources can be placed in any directory, but here we use the user's Home directory
```
# clone the repository
git clone https://github.com/biocore/sortmerna.git

# alternatively get the release sources
wget https://github.com/biocore/sortmerna/archive/v4.0.0.tar.gz

tar xzf v4.0.0.tar.gz

pushd sortmerna

# If you need a particular release (tag)
git checkout v4.0.0
```

### install CMake

Use the Python script provided with the Sortmerna distro

```
# 
SMR_HOME=$HOME/sortmerna
python $SMR_HOME/scripts/build.py --name cmake_install
  [cmake_install] Installed CMake /home/biocodz/cmake-3.15.5-Linux-x86_64/bin/cmake

# add cmake to PATH
export PATH=$HOME/cmake-3.15.5-Linux-x86_64/bin:$PATH
```

### build

Use python script provided with the Sortmerna distro.
By default the build is performed in the User directory.
All required third party libraries will be checked and installed automatically (in User directory by default)
The default build won't interfere with any existing system installation.

By default the build produces statically linked executable i.e. portable.

```
# navigate to the Sortmerna source directory SMR_HOME
SMR_HOME=$HOME/sortmerna

# modify configuration (optional)
vi $SMR_HOME/scripts/env.yaml

# run the build
python $SMR_HOME/scripts/build.py --name all [--env $SMR_HOME/script/my_env.yaml]
```

## Building on Windows

### Install Visual Studio

Download and Install VS Community edition from [Visual Studio community website](https://www.visualstudio.com/vs/community/)

### Install Conda

TODO

### get Sortmerna sources

TODO

### build

TODO

## Building on Mac

TODO

## Third-party libraries

Various features in SortMeRNA are dependent on third-party libraries, including:
* [ALP](http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/software/program.html?uid=6): computes statistical parameters for Gumbel distribution (K and Lambda)
* [CMPH](http://cmph.sourceforge.net): C Minimal Perfect Hashing Library
* [Zlib](https://github.com/madler/zlib): reading compressed Reads files
* [RocksDB](https://github.com/facebook/rocksdb): storage for SortmeRNA alignment results
* [RapidJson](https://github.com/Tencent/rapidjson): serialization of Reads objects