#!/bin/bash
# sortmerna 3.0-beta.run.tests.u16

main() {

    mkdir $HOME/bin

    # download binaries
    for BIN in "${BINS[@]}"
    do
        echo "dx download -o $HOME/bin/ \"${BIN}\""
        dx download -o $HOME/bin/ "${BIN}"
    done

    chmod u+x $HOME/bin/indexdb
    chmod u+x $HOME/bin/sortmerna
    export PATH=$HOME/bin:$PATH

    echo "ls -lrt $HOME/bin/"
    ls -lrt $HOME/bin/

    # check Python
    $HOME/miniconda/bin/python --version

    # check scikit-bio
    $HOME/miniconda/bin/conda list scikit-bio

    # check patchelf
    echo "apt-cache policy patchelf"
    apt-cache policy patchelf
    echo "which patchelf"
    which patchelf # /usr/local/bin/patchelf

    # check rocksdb
    echo "apt-cache policy librocksdb-dev"
    apt-cache policy librocksdb-dev

    # patch sortmerna to look for libstdc++.so.6 in the ORIGIN directory
    echo "Running: patchelf --set-rpath '$ORIGIN' $HOME/bin/sortmerna"
    patchelf --set-rpath '$ORIGIN' $HOME/bin/sortmerna

    # clone Sortmerna to get tests' data and code
    git clone ${sortmerna_git}

    echo "[INFO] Files in $HOME/"
    find . -type f

    # run tests
    $HOME/miniconda/bin/python $HOME/sortmerna/tests/test_sortmerna.py

    echo "==== DONE ===="
}
