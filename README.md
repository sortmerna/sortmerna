sortmerna
=========

This is ongoing development of SortMeRNA. For the latest stable release, see http://bioinfo.lifl.fr/RNA/sortmerna/.


WARNING: Do not attempt to run _create_dist_mac.sh unless you have Autotools installed. 
This script is for producing a distribution version of the code.

NOTE: the Clang compiler on Mac (distributed through Xcode) does not support multithreading.
A preliminary implementation of OpenMP for Clang has been made at "http://clang-omp.github.io"
though has not been yet incorporated into the Clang mainline. The user may follow the
steps outlined in the above link to install the version of Clang with multithreading support, 
though this version has not yet been tested with SortMeRNA. Otherwise, the user is 
recommended to install the original GCC compiler via MacPorts (contains full multithreading support).
  


To compile on Linux OS:
-----------------------

(1) Check your GCC compiler is version 4.0 or higher:

```bash
gcc --version
```

(2) Compile the code:

```bash
./configure
make
```



To compile on Mac OS:
---------------------

(1) Check the version of your C/C++ compiler:

```bash
gcc --version
```

If the compiler is LLVM-GCC (deprecated, see note on http://bioinfo.lifl.fr/RNA/sortmerna/), 
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
export CXX=gcc-mp-4.8
```

(3) Compile the code

```bash
./configure
make
```

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


