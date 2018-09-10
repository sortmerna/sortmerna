#!/bin/bash
#==============================================================================
# FILE: build.sh
# Created: Jul 20, 2018
##==============================================================================

main() {

    echo "Upload directory: ${upload_dir}"

    # Bypass the APT caching proxy that is built into the execution environment.
    # It's configured to only allow access to the stock Ubuntu repos.
    sudo mv /etc/apt/apt.conf.d/99dnanexus ./

    # install gcc-7
    #
    sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
    sudo apt update
    sudo apt install -y gcc-7 g++-7
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-7

    # install latest CMake
    sudo apt -y remove cmake
    sudo apt -y purge --auto-remove cmake
    wget https://cmake.org/files/v3.12/cmake-3.12.0.tar.gz
    tar xzf cmake-3.12.0.tar.gz
    pushd cmake-3.12.0/
    ./bootstrap
    make
    sudo make install
    popd

    # ZLib
    sudo apt install -y zlib1g-dev

    # Sortmerna - cmake
    SMR_ROOT=$PWD/sortmerna
    git clone https://github.com/biocore/sortmerna.git
    mkdir -p sortmerna/build/Release
    pushd sortmerna/build/Release
    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DWITH_TESTS=0 -DSRC_ROCKSDB=1 -DSRC_RAPIDJSON=1 -DSET_ROCKSDB=1 -DROCKSDB_INCLUDE_DIR=$SMR_ROOT/3rdparty/rocksdb/include -DROCKSDB_LIB_RELEASE=$SMR_ROOT/3rdparty/rocksdb/build/Release ../..
    popd

    # RocksDB
    mkdir -p sortmerna/3rdparty/rocksdb/build/Release
    pushd sortmerna/3rdparty/rocksdb/build/Release
    cmake -E env CXXFLAGS="-Wno-error=maybe-uninitialized" cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DPORTABLE=1 -DWITH_ZLIB=1 -DWITH_TESTS=0 -DWITH_TOOLS=0 ../..
    make
    popd

    # sortmerna build
    pushd sortmerna/build/Release
    make
    popd

    # upload built artifacts and library dependencies to dnanexus project and add output entries to 'job_output.json'
    sortmerna=$(dx upload sortmerna/build/Release/src/sortmerna/sortmerna --path ${upload_dir} --brief)
    dx-jobutil-add-output sortmerna "$sortmerna" --class=file
    indexdb=$(dx upload sortmerna/build/Release/src/indexdb/indexdb --path ${upload_dir} --brief)
    dx-jobutil-add-output indexdb "$indexdb" --class=file
    libstdcpp=""
    libstdcline=$(ldd sortmerna/build/Release/src/sortmerna/sortmerna | grep libstdc)
    strarr=($libstdcline)
    for i in "${strarr[@]}"
    do
      if [ -e $i]
      then
        echo "Uploading output: $i"
        libstdcpp=$(dx upload $i --path ${upload_dir} --brief)
        dx-jobutil-add-output libstdcpp "$libstdcpp" --class=file
        break
      fi
    done
    if [ -z "$libstdcpp"]
    then
        echo "ERROR: libstdc++.so.6 Not Found"
    fi
} # ~main
