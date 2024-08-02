# sortmerna

SortMeRNA is a local sequence alignment tool for filtering, mapping and clustering.

The core algorithm is based on approximate seeds and allows for sensitive analysis of NGS reads.
The main application of SortMeRNA is filtering rRNA from metatranscriptomic data.
SortMeRNA takes as input files of reads (fasta, fastq, fasta.gz, fastq.gz) and one or multiple
rRNA database file(s), and sorts apart aligned and rejected reads into two files. SortMeRNA works
with Illumina, Ion Torrent and PacBio data, and can produce SAM and BLAST-like alignments.

SortMeRNA is also available through [QIIME v1.9.1](http://qiime.org) and
the [nf-core RNA-Seq pipeline v.3.9](https://nf-co.re/rnaseq/3.9).

## Table of Contents

- [Getting Started](#getting-started)
  - [Using Conda package](#using-conda-package)
  - [Using GitHub release binaries on Linux](#using-github-release-binaries-on-linux)
  - [Running](#running)
    - [Execution trace](#execution-trace)
- [Building from sources](#building-from-sources)
- [User Manual](#user-manual)
- [Taxonomies](#taxonomies)
- [Citation](#citation)
- [Contributors](#contributors)
- [Support](#support)


## Getting Started

SortMeRNA 4 is C++17 compliant, and mostly uses standard libraries. It uses CMake as the build system, and can be run/built on all major OS including Linux, Windows, and Mac, on AMD64 and ARM64 processors.

### Using Conda package

Install conda - [official docs](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
The conda packages before Sortmerna 4.3.7 were hosted on Bioconda. Starting with 4.3.7 the packages are hosted on conda-forge.
Erroneously an empty 4.3.7 package made its way to Bioconda, and should be ignored until removed (from Bioconda).

Currently the build on conda-forge still waiting to be merged. Until it is ready, the local installation package can be used:
```
# == only for 4.3.7 until ready on conda-forge ==
# download the conda-build package into a directory of your choice e.g. Downloads/
wget https://github.com/sortmerna/sortmerna/releases/download/v4.3.7/sortmerna-4.3.7-conda-linux-64.tar.bz2 -P ~/Downloads/

# create a new environment and install SortMeRNA in it
conda create --name sortmerna
conda activate sortmerna
conda install ~/Downloads/sortmerna-4.3.7-conda-linux-64.tar.bz2

which sortmerna  # check the installed binary e.g. miniforge3/envs/sortmerna/bin/sortmerna 
sortmerna -h
```
For versions older then 4.3.7 per the [Bioconda guidelines](https://bioconda.github.io), add the following conda channels:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict


conda search sortmerna
  Loading channels: done
  # Name                       Version           Build  Channel
  sortmerna                        2.0               0  bioconda
  ...
  sortmerna                      4.3.4               0  bioconda
  ...
  sortmerna                      4.3.6               0  bioconda
  ...
  sortmerna                      4.3.7      hdbdd923_1  bioconda <- (!) ignore - corrupt, see instructions above

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
```

### Using GitHub release binaries on Linux

Visit [Sortmerna GitHub Releases](https://github.com/biocore/sortmerna/releases)

Linux distribution is a Shell script with the embedded installation archive.

Issue the following bash commands:

```
pushd ~

# get the distro
wget https://github.com/biocore/sortmerna/releases/download/v4.3.6/sortmerna-4.3.6-Linux.sh

# view the installer usage
bash sortmerna-4.3.6-Linux.sh --help
    Options: [defaults in brackets after descriptions]
      --help            print this message
      --version         print cmake installer version
      --prefix=dir      directory in which to install
      --include-subdir  include the sortmerna-4.3.6-Linux subdirectory
      --exclude-subdir  exclude the sortmerna-4.3.6-Linux subdirectory
      --skip-license    accept license

# run the installer
bash sortmerna-4.3.6-Linux.sh --skip-license
  sortmerna Installer Version: 4.3.6, Copyright (c) Clarity Genomics
  This is a self-extracting archive.
  The archive will be extracted to: $HOME/sortmerna
  
  Using target directory: /home/biocodz/sortmerna
  Extracting, please wait...
  
  Unpacking finished successfully

# check the installed binaries
ls -lrt /home/biocodz/sortmerna/bin/
sortmerna

# set PATH
export PATH=$HOME/sortmerna/bin:$PATH

# test the installation
sortmerna --version
  SortMeRNA version 4.3.6
  Build Date: Jul 17 2021
  sortmerna_build_git_sha:@921fa40256760ea2d44c49b21eb326afda748d5e@
  sortmerna_build_git_date:@2022/08/16 10:59:31@

# view help
sortmerna -h
```

### Running

* The only required options are `--ref` and `--reads`
* Options (any) can be specified usig a single dash e.g. `-ref` and `-reads`
* Both plain `fasta/fastq` and archived `fasta.gz/fastq.gz` files are accepted
* file extensions `.fastq, .fastq.gz, .fq, .fq.gz, .fasta, ...` are optional. The format and compression are automatically recognized
* Relative paths are accepted

for example

```
# single reference and single reads file
sortmerna --ref REF_PATH --reads READS_PATH

# for multiple references use multiple '--ref'
sortmerna --ref REF_PATH_1 --ref REF_PATH_2 --ref REF_PATH_3 --reads READS_PATH

# for paired reads use '--reads' twice
sortmerna --ref REF_PATH_1 --ref REF_PATH_2 --ref REF_PATH_3 --reads READS_PATH_1 --reads READS_PATH_2

```

More examples can be found in [test.jinja](https://github.com/biocore/sortmerna/blob/master/scripts/test.jinja) and [run.py](https://github.com/biocore/sortmerna/blob/master/scripts/run.py)

#### Execution trace

Here is a [sample execution trace](https://sortmerna.readthedocs.io/en/latest/trace4.3.2.html).  

`IMPORTANT`
- Progressing execution trace showing the number of reads processed so far indicates a normally running program. 
- Non-progressing trace means a problem. Please, kill the process (no waiting for two days), and file an issue [here](https://github.com/biocore/sortmerna/issues)  
- please, provide the execution trace when filing issues.

[Sample execution statistics](https://github.com/biocore/sortmerna/wiki/sample-execution-statistics) are provided to give an idea on what the execution time might be.

## Building from sources

[Build instructions](https://sortmerna.readthedocs.io/en/latest/building.html)

## User Manual

See [Sortmerna Read The Docs project](https://sortmerna.readthedocs.io/en/latest/index.html).

In case you need PDF, any modern browser can print web pages to PDF.

## Taxonomies

The folder `data/rRNA_databases/silva_ids_acc_tax.tar.gz` contains SILVA taxonomy strings (extracted from XML file generated by ARB)
for each of the reference sequences in the representative databases. The format of the files is three tab-separated columns,
the first being the reference sequence ID, the second being the accession number and the final column is the taxonomy.

## Citation

If you use SortMeRNA, please cite:
Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611.

## Contributors

See [AUTHORS](./AUTHORS) for a list of contributors to this project.

## Support

For questions and comments, feel free to file an [issue](https://github.com/sortmerna/sortmerna/issues), or start a [discussion](https://github.com/sortmerna/sortmerna/discussions).
	
