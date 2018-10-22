SHELL=/bin/bash
BIOCHANNEL=https://conda.anaconda.org/biobuilds

.ONESHELL:
all:
	# Bypass the APT caching proxy that is built into the execution environment
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus

	uname −a
	lsb_release −a
	which bash
	bash −−version

	# install Conda Python 3 and bio packages
	#
    #@echo "PYTHONPATH: $PYTHONPATH"
    #unset PYTHONPATH
    #@echo "wget ${conda_url} -O miniconda.sh"
    wget ${conda_url} -O miniconda.sh
    chmod +x miniconda.sh
    ./miniconda.sh -b -p $HOME/miniconda &> /dev/null
    #@echo "Running: conda update --yes conda"
    $HOME/miniconda/bin/conda update --yes conda
    #@echo "Running: conda config --add channels ${bio_channel}"
    $HOME/miniconda/bin/conda config --add channels ${BIOCHANNEL}
    #@echo "Running: conda create --yes -n env_name python=3.5 numpy scipy"
    #conda create --yes -n env_name python=3.5 numpy scipy
    #@echo "Running: source activate env_name"
    #source activate env_name
    #@echo "Running: pip install scikit-bio==0.2.3"
    $HOME/miniconda/bin/conda install -q scikit-bio

    @echo "$HOME/miniconda/bin/python --version"
    $HOME/miniconda/bin/python --version
    echo "find $HOME -type f"
    find $HOME -type f


	# clean-up
	@echo "Before cleanup:"
	df -h
	sudo apt autoremove --purge
	@echo "Done: apt autoremove --purge"
	@echo "After cleanup:"
	df -h