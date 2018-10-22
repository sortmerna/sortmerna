SHELL=/bin/bash
BIO_CHANNEL=https://conda.anaconda.org/biobuilds
CONDA_URL=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

.ONESHELL:
all:
	# Bypass the APT caching proxy that is built into the execution environment
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus

	uname −a
	bash --version
	lsb_release -a

	# install Conda and bio packages
	#
	pushd /opt && \
	unset PYTHONPATH && \
	wget $(CONDA_URL) -O miniconda.sh && \
	chmod +x miniconda.sh && \
	./miniconda.sh -b -p /opt/miniconda &> /dev/null && \
	/opt/miniconda/bin/conda update -y -q conda && \
	/opt/miniconda/bin/conda config --add channels $(BIO_CHANNEL) && \
	/opt/miniconda/bin/conda install -q -y scikit-bio

	/opt/miniconda/bin/python --version
	/opt/miniconda/bin/conda list python
	/opt/miniconda/bin/conda list scikit-bio
	/opt/miniconda/bin/conda list numpy
	/opt/miniconda/bin/conda list scipy

	# clean-up     TODO: Causes problems - leave for now
	#$(info "Before cleanup:")
	#out := $(shell df -h)
	#$(info $(out))
	#sudo apt autoremove --purge
	#@echo "Done: apt autoremove --purge"
	#$(info "After cleanup:")
	#out := $(shell df -h)
	#$(info $(out))