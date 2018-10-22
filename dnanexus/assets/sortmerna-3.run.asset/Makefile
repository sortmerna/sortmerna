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
	wget $(CONDA_URL) -O miniconda.sh && \
	chmod +x miniconda.sh && \
	./miniconda.sh -b -p $(HOME)/miniconda &> /dev/null && \
	$(HOME)/miniconda/bin/conda update -y -q conda && \
	$(HOME)/miniconda/bin/conda config --add channels $(BIO_CHANNEL) && \
	$(HOME)/miniconda/bin/conda install -q -y scikit-bio

	$(HOME)/miniconda/bin/python --version
	$(HOME)/miniconda/bin/conda list python
	$(HOME)/miniconda/bin/conda list scikit-bio
	$(HOME)/miniconda/bin/conda list numpy
	$(HOME)/miniconda/bin/conda list scipy

	# clean-up
	$(info "Before cleanup:")
	out := $(shell df -h)
	$(info $(out))
	sudo apt autoremove --purge
	@echo "Done: apt autoremove --purge"
	$(info "After cleanup:")
	out := $(shell df -h)
	$(info $(out))