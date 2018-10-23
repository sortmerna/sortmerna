SHELL=/bin/bash
CMAKE_TAR=https://cmake.org/files/v3.12/cmake-3.12.0.tar.gz

.ONESHELL:
all:
	# Bypass the APT caching proxy that is built into the execution environment
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus

	# these commands fail to execute. DNANexus bug?
	#uname −a  # uname: extra operand ‘−a’
	#lsb_release −a # lsb_release: error: No arguments are permitted
	#which bash
	#bash −−version # bash: −−version: No such file or directory
	#$(shell uname -a)
	#$(shell lsb_release -a)
	#$(shell bash --version)

	# install gcc-7
	#
	sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
	sudo apt update && \
	sudo apt install -y gcc-7 g++-7 && \
	sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-7
	which gcc-7
	gcc-7 --version

	# install latest CMake 3.12.0
	sudo apt -y remove cmake && sudo apt -y purge --auto-remove cmake
	wget $(CMAKE_TAR) && tar xzf cmake-3.12.0.tar.gz
	pwd
	ls -lrt ./cmake-3.12.0
	pushd ./cmake-3.12.0 && ./bootstrap && sudo make install && popd
	which cmake
	cmake --version

	# clean-up
	@echo "Before cleanup:"
	df -h
	sudo apt autoremove --purge
	@echo "Done: apt autoremove --purge"
	@echo "After cleanup:"
	df -h