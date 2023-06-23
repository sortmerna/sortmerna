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
###########

If you already have GCC, CMake and Conda, the following command will compile SortMeRNA local sources in the current working directory on a linux machine::

   python scripts/build.py --name all --local-linux

Install GCC 9
#############

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
#####################

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
#############

Use the :code:`build.py` python script provided with Sortmerna distro. The following installs Conda, and the python packages :code:`pyyaml`, and :code:`jinja2` in the User's Home directory::
   
   SMR_HOME=$HOME/sortmerna
   python $SMR_HOME/scripts/build.py --name cmake
   
   ls -lrt
   drwxrwxr-x 15 biocodz biocodz     4096 Nov 18 09:43 miniconda3
   
   # add Conda binaries to the PATH
   export PATH=$HOME/miniconda3/bin:$PATH

Install CMake
#############

The following installs CMake in user's home directory::

   SMR_HOME=$HOME/sortmerna
   python $SMR_HOME/scripts/build.py --name cmake
     [cmake_install] Installed CMake /home/biocodz/cmake-3.15.5-Linux-x86_64/bin/cmake
   
   # add cmake to PATH
   export PATH=$HOME/cmake-3.15.5-Linux-x86_64/bin:$PATH

Build
#####

All required third party libraries will be checked and installed automatically (in User directory by default) The default build won't interfere with any existing system installation. By default the build produces statically linked executable i.e. portable.

::

   SMR_HOME=$HOME/sortmerna
   
   # modify configuration (optional)
   vi $SMR_HOME/scripts/env.yaml
   
   # run the build
   python $SMR_HOME/scripts/build.py --name all [--env $SMR_HOME/script/my_env.yaml]
