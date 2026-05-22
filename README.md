# sortmerna

SortMeRNA is a local sequence alignment tool for filtering, mapping and clustering.

The core algorithm is based on approximate seeds and allows for sensitive analysis of NGS reads.
The main application of SortMeRNA is filtering rRNA from metatranscriptomic data.
SortMeRNA takes as input files of reads (fasta, fastq, fasta.gz, fastq.gz) and one or multiple
rRNA database file(s), and sorts apart aligned and rejected reads into two files. SortMeRNA works
with Illumina, Ion Torrent and PacBio data, and can produce SAM and BLAST-like alignments.

SortMeRNA is also available through the [nf-core RNA-Seq pipeline v.3.26.0](https://nf-co.re/rnaseq/3.26.0).

## Table of Contents

- [What's new in 7.0.0](#whats-new-in-700)
- [Getting Started](#getting-started)
  - [Using Conda package](#using-conda-package)
  - [Using GitHub release binaries on Linux](#using-github-release-binaries-on-linux)
  - [Running](#running)
    - [Execution trace](#execution-trace)
    - [Read feed modes (gzipped FASTQ)](#read-feed-modes-gzipped-fastq)
- [Building from sources](#building-from-sources)
- [User Manual](#user-manual)
- [Databases](#databases)
- [Taxonomies](#taxonomies)
- [Citation](#citation)
- [Contributors](#contributors)
- [Third-party dependencies](#third-party-dependencies)
- [Support](#support)

## What's new in 7.0.0

### Resume after interruption (headline)

A sortmerna run that is killed mid-alignment now resumes automatically
on the next invocation against the same kvdb. There is no `--resume`
flag: presence of a non-empty kvdb plus matching opts/input fingerprints
is sufficient. If any alignment-relevant option or input file fingerprint
has changed since the previous run, sortmerna aborts with a diagnostic
naming the kvdb so the user can either restore the original inputs or
wipe the kvdb. Internally, progress is persisted at pass boundaries via
atomic kvdb batches; if a pass is interrupted, the partial work is
rolled back on resume before the pass is retried, so end-of-run
statistics and outputs are identical to an uninterrupted run.

### Testing & dev

`scripts/run.py test ... --interrupt [N]` exercises the auto-resume
path by SIGKILLing the binary mid-alignment N times (max 2, default
1) and re-launching, then validating the final output normally.

## Getting Started

SortMeRNA 6 is C++17 compliant. It uses CMake as the build system, and can be run/built on Linux, Mac, and on AMD64 and ARM64 architectures. Support for native Windows was dropped in version 6.0. WSL can be used to run Sortmerna on Windows.

### Using Conda package

Install conda - [official docs](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Sortmerna conda packages are hosted on conda-forge.

```
conda config --add channels conda-forge

conda search sortmerna

 Name      Version Build                         Channel     Subdir  
──────────────────────────────────────────────────────────────────────
 sortmerna 6.0.0   py310hea66570_0 (+  4 builds) conda-forge linux-64
 sortmerna 5.0.0   py310hea66570_0 (+  4 builds) conda-forge linux-64
 sortmerna 4.4.0   py310hea66570_0 (+  4 builds) conda-forge linux-64
 sortmerna 4.3.7   py310hea66570_0 (+ 13 builds) conda-forge linux-64

# create a new environment and install SortMeRNA in it
conda create --name sortmerna_run
conda activate sortmerna_run
conda install sortmerna
which sortmerna
  /home/biocodz/miniconda3/envs/sortmerna_env/bin/sortmerna

# check version
sortmerna --version

# view help
sortmerna -h
```

### Using GitHub release binaries on Linux

Visit [Sortmerna GitHub Releases](https://github.com/sortmerna/sortmerna/releases)

Linux distribution is a Shell script with the embedded installation archive.

Issue the following bash commands:

```
pushd ~

# get the distro
wget https://github.com/sortmerna/sortmerna/releases/download/v6.0.0/sortmerna-6.0.0-Linux.sh

# view the installer usage
bash sortmerna-6.0.0-Linux.sh --help
    Options: [defaults in brackets after descriptions]
      --help            print this message
      --version         print cmake installer version
      --prefix=dir      directory in which to install
      --include-subdir  include the sortmerna-6.0.0-Linux subdirectory
      --exclude-subdir  exclude the sortmerna-6.0.0-Linux subdirectory
      --skip-license    accept license

# run the installer
bash sortmerna-6.0.0-Linux.sh --skip-license
  sortmerna Installer Version: 6.0.0, Copyright (c) Clarity Genomics
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
  SortMeRNA version 6.0.0
  Build Date: May 15 2026
  sortmerna_build_git_sha:@c750937be9a37bfde9a3d1d5157fe185becd384e@
  sortmerna_build_git_date:@2026/05/15 10:41:11@

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

More examples can be found in [test.jinja](https://github.com/sortmerna/sortmerna/blob/master/scripts/test.jinja) and [run.py](https://github.com/sortmerna/sortmerna/blob/master/scripts/run.py)

#### Execution trace

Here is a [sample execution trace](https://sortmerna.readthedocs.io/en/latest/trace4.3.2.html).  

`IMPORTANT`
- Progressing execution trace showing the number of reads processed so far indicates a normally running program. 
- Non-progressing trace means a problem. Please, kill the process (no waiting for two days), and file an issue [here](https://github.com/sortmerna/sortmerna/issues)  
- please, provide the execution trace when filing issues.

[Sample execution statistics](https://github.com/sortmerna/sortmerna/wiki/sample-execution-statistics) are provided to give an idea on what the execution time might be.

#### Read feed modes (gzipped FASTQ)

Since 5.0.0 SortMeRNA processes both indexed and flat files in parallel rather than physically
splitting them. SortMeRNA automatically defaults to indexed in-memory feed (`INDEXED_GZ`) that 
uses [rapidgzip](https://github.com/mxmlnkn/rapidgzip) to assign each worker thread its own 
byte range — no split files are written to disk.

The legacy `--readfeed 1` on-disk splitting still works but is deprecated. 
The previously-documented `run.py split` pre-splitting workflow is no longer necessary
and not advised.

## Building from sources

[Build instructions](https://sortmerna.readthedocs.io/en/latest/building.html)

## User Manual

See [Sortmerna Read The Docs project](https://sortmerna.readthedocs.io/en/latest/index.html).

In case you need PDF, any modern browser can print web pages to PDF.

## Databases

[database.tar.gz](https://github.com/sortmerna/sortmerna/releases/download/v4.3.4/database.tar.gz) provided since release 4.3.4. The tarball contains 4 separate files.

Since the version 4.3.7 the four database files are also provided as separate .gz archives.

We recommend to use smr_v4.3_default_db.fasta.

Original source databases (clustering parameters given below):
* Silva 138 SSURef NR99 (16S, 18S)
* Silva 132 LSURef (23S, 28S)
* RFAM v14.1 (5S, 5.8S)

The difference between the databases is the % ID for clustering the sequences for each kingdom + rRNA component.

Specifically,

* smr_v4.3_fast_db.fasta
  * bac-16S 85%, 5S & 5.8S seeds, rest 90% (benchmark accuracy: 99.888%)
* smr_v4.3_default_db.fasta
  * bac-16S 90%, 5S & 5.8S seeds, rest 95% (benchmark accuracy: 99.899%)
* smr_v4.3_sensitive_db.fasta
  * all 97% (benchmark accuracy: 99.907%)
* smr_v4.3_sensitive_db_rfam_seeds.fasta
  * all 97%, except RFAM database which includes the full seed database sequences

The accuracy (based on sensitivity and selectivity) is very good for all databases, however the "sensitive" databases will run at least 2x slower.

## Taxonomies

The folder `data/rRNA_databases/silva_ids_acc_tax.tar.gz` contains SILVA taxonomy strings (extracted from XML file generated by ARB)
for each of the reference sequences in the representative databases. The format of the files is three tab-separated columns,
the first being the reference sequence ID, the second being the accession number and the final column is the taxonomy.

## Citation

If you use SortMeRNA, please cite:
Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611.

## Contributors

See [AUTHORS](./AUTHORS) for a list of contributors to this project.

Maintained by the team at [Clarity Genomics Inc](https://www.clarity-genomics.com) and [Centre de Recherche en Informatique, Signal et Automatique de Lille](https://www.cristal.univ-lille.fr/).

## Third-party dependencies

Refer to `3rdparty.jinja` for the full list of bundled libraries. In 6.0.0 these include: 
Zlib, RocksDB, Rapidgzip (indexed_bzip2), BBHash (replaces CMPH), Parasail (replaces SSW), 
ALP, and concurrentqueue.

## Support

For questions and comments, feel free to file an [issue](https://github.com/sortmerna/sortmerna/issues), or start a [discussion](https://github.com/sortmerna/sortmerna/discussions).
