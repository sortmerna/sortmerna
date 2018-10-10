#!/bin/bash
# sortmerna-3.run.tests.u16

main() {

    echo "find $HOME -type f"
    find $HOME -type f

    # prepare binaries and set PATH
    chmod u+x $HOME/bin/indexdb
    chmod u+x $HOME/bin/sortmerna
    export PATH=$HOME/bin:$HOME/miniconda/bin:$PATH

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

    # install Conda Python 3 and bio packages
    echo "PYTHONPATH: $PYTHONPATH"
    unset PYTHONPATH
    echo "wget ${conda_url} -O miniconda.sh"
    wget ${conda_url} -O miniconda.sh
    chmod +x miniconda.sh
    ./miniconda.sh -b -p $HOME/miniconda &> /dev/null
    echo "Running: conda update --yes conda"
    conda update --yes conda
    echo "Running: conda config --add channels ${bio_channel}"
    conda config --add channels ${bio_channel}
    echo "Running: conda create --yes -n env_name python=3.5 numpy scipy"
    conda create --yes -n env_name python=3.5 numpy scipy
    echo "Running: source activate env_name"
    source activate env_name
    echo "Running: pip install scikit-bio==0.2.3"
    pip install -q scikit-bio==0.2.3

    ls -lrt $HOME/bin
    echo "[INFO] PATH: $PATH"

    # clone Sortmerna to get tests' data and code
    git clone ${sortmerna_git}

    echo "[INFO] Files in $HOME/"
    find . -type f

    # run tests
    python $HOME/sortmerna/tests/test_sortmerna.py

    echo "==== DONE ===="
}
