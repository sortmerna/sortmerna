Installation
============

Using Conda
-----------

Install Conda::

   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh

As per the `Bioconda guidelines
<https://bioconda.github.io/>`_, add the following conda channels::

   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict

Then install SortMeRNA::

   conda search sortmerna
     Loading channels: done
     # Name                       Version           Build  Channel
     sortmerna                        2.0               0  bioconda
     ...
     sortmerna                      4.3.4               0  bioconda
     ...
     sortmerna                      4.3.6               0  bioconda

   # create a new environment and install SortMeRNA in it
   conda create --name sortmerna_env
   conda activate sortmerna_env
   conda install sortmerna
   which sortmerna
     /home/biocodz/miniconda3/envs/sortmerna_env/bin/sortmerna

   # test the installation
   sortmerna --version
     SortMeRNA version 4.3.6
     Build Date: Aug 16 2022
     sortmerna_build_git_sha:@db8c1983765f61986b46ee686734749eda235dcc@
     sortmerna_build_git_date:@2022/08/16 11:42:59@

   # view help
   sortmerna -h