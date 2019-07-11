#!/bin/bash
#===============================================================================
# FILE: build.sh
# Created: Sep 09, 2018 Sun
#
#===============================================================================

function misc {
    sudo apt install -y zlib1g-dev librocksdb-dev rapidjson-dev

    gcc-9 --version
    cmake --version
    apt-cache policy librocksdb-dev
    apt-cache policy zlib1g-dev
    apt-cache policy rapidjson-dev

    SMR_ROOT=$PWD/sortmerna
    git clone https://github.com/biocore/sortmerna.git
    mkdir -p $SMR_ROOT/build/Release
    pushd $SMR_ROOT/build/Release
    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DWITH_TESTS=0 ../..
    # sortmerna build
    make
    popd
}

function build_wsl {
    #
    # @param 1 String  Optional   type of build: Release | Debug
    #
    # Build on Windows Linux (WSL) - sources on Windows side
    #
    # Usage: /mnt/c/a02_projects/sortmerna/build.sh wsl [Release]
    #
    if [ -z ${1+x} ]; then BUILD_TYPE="Release" ; else BUILD_TYPE=$1 ; fi

    BUILD_HOME=/home/biocodz/sortmerna

    ZLIB_ROOT=/home/biocodz/zlib/dist
    ZLIB_LIBRARY_RELEASE=/home/biocodz/zlib/dist/lib/libz.a
    ROCKSDB_HOME=/home/biocodz/rocksdb/dist
    ROCKSDB_SRC=/mnt/c/a02_projects/rocksdb
    RAPIDJSON_HOME=/home/biocodz/rapidjson/dist
    SMR_SRC_DIR=/mnt/c/a02_projects/sortmerna
    CMAKE_INSTALL_PREFIX=$BUILD_HOME/dist

    mkdir -p $BUILD_HOME/build/$BUILD_TYPE
    pushd $BUILD_HOME/build/$BUILD_TYPE

    #CC=gcc-7 CXX=g++-7 \
    cmake -G "Unix Makefiles" \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX \
    -DZLIB_ROOT=$ZLIB_ROOT \
    -DZLIB_LIBRARY_RELEASE=$ZLIB_LIBRARY_RELEASE \
    -DROCKSDB_HOME=$ROCKSDB_HOME \
    -DROCKSDB_SRC=$ROCKSDB_SRC\
    -DRAPIDJSON_HOME=$RAPIDJSON_HOME \
    -DCPACK_BINARY_TGZ=ON $SMR_SRC_DIR

    # install
    cmake --build . --target install

    # test
    $BUILD_HOME/dist/bin/sortmerna -version
    #$BUILD_HOME/dist/bin/sortmerna -h
}

if [ "${1}" == "misc" ]; then misc
elif [ "${1}" == "wsl" ]; then build_wsl $2 $3
else
	echo "Unknown target: [${1}]"
fi