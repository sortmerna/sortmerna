#!/bin/bash
#==============================================================================
# FILE: build.sh
# Created: AUg 23, 2018
##==============================================================================

main() {

    echo "Upload directory: ${upload_root}"

    # get executable name from the Job description to use in output (TODO: is there a better way?)
    tokens=($(dx describe ${DX_JOB_ID} | grep 'Executable name'))
    exe_name=${tokens[2]}
    upload_path="${upload_root}/$exe_name"
    echo "Applet/Executable name: $exe_name" # e.g. sortmerna-3.0-beta.build.on.asset
    echo "Upload will be done into $upload_path"
    stat=$(dx mkdir -p $upload_path)
    echo "Upload path status: $stat"

    # check gcc-7 - part of the asset
    #
    echo "Testing: which gcc-7"
    which gcc-7
    echo "Testing: gcc-7 --version"
    gcc-7 --version

    # check RocksDB, zlib, rapidjson - part of the asset
    echo "Testing: apt-cache policy librocksdb-dev"
    apt-cache policy librocksdb-dev
    echo "Testing: apt-cache policy zlib1g-dev"
    apt-cache policy zlib1g-dev
    echo "Testing: apt-cache policy rapidjson-dev"
    apt-cache policy rapidjson-dev

    # Sortmerna - cmake generate build files
    SMR_ROOT=$PWD/sortmerna
    git clone https://github.com/biocore/sortmerna.git
    mkdir -p sortmerna/build/Release
    pushd sortmerna/build/Release
    cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DWITH_TESTS=0 ../..
    popd

    # sortmerna build
    pushd sortmerna/build/Release
    make
    popd

    # upload built artifacts and library dependencies to dnanexus project and add output entries to 'job_output.json'
    sortmerna=$(dx upload sortmerna/build/Release/src/sortmerna/sortmerna --path ${upload_path}/ --brief)
    dx-jobutil-add-output sortmerna "$sortmerna" --class=file
    indexdb=$(dx upload sortmerna/build/Release/src/indexdb/indexdb --path ${upload_path}/ --brief)
    dx-jobutil-add-output indexdb "$indexdb" --class=file
    libstdcpp=""
    libstdcline=$(ldd sortmerna/build/Release/src/sortmerna/sortmerna | grep libstdc)
    echo "libstdcline: $libstdcline"
    strarr=($libstdcline)
    for i in "${strarr[@]}"
    do
      if [ -e $i ]; then
        echo "Uploading output: $i"
        libstdcpp=$(dx upload $i --path ${upload_path}/ --brief)
        dx-jobutil-add-output libstdcpp "$libstdcpp" --class=file
        break
      fi
    done
    if [ -z "$libstdcpp" ]; then
        echo "ERROR: libstdc++.so.6 Not Found"
    fi
} # ~main
