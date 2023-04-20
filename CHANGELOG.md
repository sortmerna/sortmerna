# SortMeRNA: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[v4.3.7]] 

### Fixed

- [[#361](https://github.com/sortmerna/sortmerna/issues/361))] - updated README ``--version` output

## [[v4.3.6](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.6)] - 2021-07-21

### Fixed

- [[#312](https://github.com/sortmerna/sortmerna/issues/312)] - fixed an issue where complex output file names where being altered
- [[#328](https://github.com/sortmerna/sortmerna/issues/328)] - fixed an issue where a missing `-read` option wasn't being handle properly

## [[v4.3.5](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.5-pre)] - 2021-07-21

### Fixed

- [[#316](https://github.com/sortmerna/sortmerna/issues/316)] - fixed an issue where `--sq` option wasn't functioning
- [[#321](https://github.com/sortmerna/sortmerna/issues/321)] - fixed an issue where very short fastq records couldn't be processed
- [[#322](https://github.com/sortmerna/sortmerna/issues/322)] - fixed an issue where `--mismatch` option couldn't take negative values
- [[#330](https://github.com/sortmerna/sortmerna/issues/330)] - fixed an issue where `--zip-out` option couldn't take non-integer values

## [[v4.3.4](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.4)] - 2021-07-21

### Fixed

- [[#288](https://github.com/sortmerna/sortmerna/issues/288)] - fixed an issue where gzipped reads were becoming corrupt during splitting

## [[v4.3.3](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.3)] - 2021-05-25

### Added

- new tests added including coverage for recent bugs

## [[v4.3.3-pre](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.3-pre)] - 2021-05-10

### Fixed

- [[#288](https://github.com/sortmerna/sortmerna/issues/288)]
- [[#290](https://github.com/sortmerna/sortmerna/issues/290)]

## [[v4.3.2](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.2)] - 2021-04-02

### Fixed 

- [[#221](https://github.com/sortmerna/sortmerna/issues/221)], [[#263](https://github.com/sortmerna/sortmerna/issues/263)], [[#272](https://github.com/sortmerna/sortmerna/issues/272)] , [[#283](https://github.com/sortmerna/sortmerna/issues/283)] - fixed an issue where processing would get stuck when using reads files made of concatenated gzipped pieces

## [[v4.3.1](https://github.com/sortmerna/sortmerna/releases/tag/v4.3.1)] - 2021-03-28

### Added

- [[#247](https://github.com/sortmerna/sortmerna/issues/247)] - added `sout` and `out2` options to allow separation of paired-end reads output into 'paired aligned' and 'singleton aligned' files
- added possibility to generate compressed output
- allowing for kmer index generation prior running the alignment (new option `-index`)
- `--dbg-level` option to control verbosity of the execution trace

### Changed

- [[#254](https://github.com/sortmerna/sortmerna/issues/254)] - changed default number of threads (`-threads`) to 2
- Reference RNA databases `database.tar.gz`
```
140M smr_v4.3_default_db.fasta
 64M smr_v4.3_fast_db.fasta
399M smr_v4.3_sensitive_db.fasta
379M smr_v4.3_sensitive_db_rfam_seeds.fasta
```

#### Architectural changes


`Split reads` - a new pre-processing step for splitting the original reads files into parts equal the number of processing threads. I.e. if there are 2 reads files and 100 threads are used, then 200 split files are generated prior starting the alignment, so that each thread uses its own files. This added the pre-processing overhead but improved the efficiency of the multithreaded processing. The advantages of the new architecture:
- circumvents the problem of random access to compressed files from multiple threads
- eliminates the use of shared objects concurrently accessed by the processing threads
- allows for performing the alignment on a widely distributed cluster, which might come important when processing terabyte and higher scale input

The new folder `readb` is added to the working directory for holding the split reads. De-facto this is a new database in addition to the kmer index, and the key-value databases. This new DB is planned as a precursor of a more sophisticated DB aimed at in the future for storing and accessing the reads.

### Fixed

- [[#231](https://github.com/sortmerna/sortmerna/issues/231)] - fixed issue where mulriple processing threads were inneficiently idling in context switching

### Removed

- `-best` option - `best` is now the default stategy, added the`--no-best` option to change search staretgy

## [[v4.2.0](https://github.com/sortmerna/sortmerna/releases/tag/v4.2.0)] - 2020-03-09

### Added

- added `tar.gz` to use with Conda recipe
- added `database.tar.gz` with new database files

### Changed

- [[#216](https://github.com/sortmerna/sortmerna/issues/216)] - modified `-workdir`, `-kvdb`, `idx`, `-aligned`, and  `-other` options. The modifications give the user a total control on naming and locations for the output, the key-value DB, and the Index.

## [[v4.1.0](https://github.com/sortmerna/sortmerna/releases/tag/v4.1.0)] - 2020-02-25

### Added

Two new boolean options added for processing paired reads:
- `--paired` to process a single reads file with paired reads
- [[#202](https://github.com/sortmerna/sortmerna/issues/202)] - `--out2` to output paired reads into separate files

Excerpt from help: 
```
    --paired          BOOL        Optional  Indicates a Single reads file with paired reads         False
                                            If a single reads file with paired reads is used,
                                            and neither 'paired_in' nor 'paired_out' are specified,
                                            use this option together with 'out2' to output
                                            FWD and REV reads into separate files

    --out2            BOOL        Optional  Output paired reads into separate files.                False
                                            Must be used with 'fastx'.
                                            Ignored without either of 'paired_in' |
                                            'paired_out' | 'paired' | two 'reads'
```

## [[v4.0.0](https://github.com/sortmerna/sortmerna/releases/tag/v4.0.0)] - 2019-12-02

### Added

- support for accepting two reads files for paired reads
- both plain fasta/fastq and archived fasta.gz/fastq.gz files are automatically recognized
- support for relative file paths

### Changed

- single executable `sortmerna` (no more `indexdb`)
- now builds on C++17 standard

### Fixed

## [[v3.0.3](https://github.com/sortmerna/sortmerna/releases/tag/v3.0.4)] - 2019-10-17

This release fixes a few bugs. The installation file contains a statically linked executable i.e. self-contained, so should be good on any Linux distro (tested on Ubuntu 16.04, 18.04, and Centos 7.7)

## [[v3.0.3](https://github.com/sortmerna/sortmerna/releases/tag/v3.0.3)] - 2019-01-16

### Added

- missing functionality to automatically package and install the binaries after building

### Fixed

- [[#184](https://github.com/sortmerna/sortmerna/issues/184)] - The file `sortmerna-3.0.3-Linux_C6_static.zip` built on Centos 6 (courtesy of [unode](https://gitlab.com/unode) bundles all the dependencies including the system libraries to fix issue 184

## [[v3.0.2](https://github.com/sortmerna/sortmerna/releases/tag/v3.0.2)] - 2018-11-21

### Fixed

- [[#137](https://github.com/sortmerna/sortmerna/issues/137)] - fixed an issue where reads were being mapped more than once on the same reference

### Changed

`sortmerna-3.0.2-linux.tar.gz` contains artifacts built on Ubuntu 16 with RocksDB 4.1. It includes `indexdb`, `sortmerna`, and `libstdc++.so.6`. The `libstdc++.so.6` should normally be available on the Debian flavoured systems and ignored (deleted from the distribution). 
 
The `sortmerna-3.0.2-debian_with_3deps.tar.gz` contains the same `indexdb`, and `sortmerna` executables as above. It also includes 3 shared libraries: `librocksdb.so.4.1`, `libgflags.so.2`, `libsnappy.so.1`. This distribution should work both on Ubuntu 16 and 18 (some basic tests were performed). This archive was prepared to resolve issues [[#173](https://github.com/sortmerna/sortmerna/issues/173)], and [[174](https://github.com/sortmerna/sortmerna/issues/174)] 
 
Note that `sortmerna` binary was patched (RPATH=ORIGIN), so as to look for dependencies in the same folder.

## [[v3.0.1](https://github.com/sortmerna/sortmerna/releases/tag/v3.0.1)] - 2018-11-02

### Changed

- [[#171](https://github.com/sortmerna/sortmerna/issues/171)] - `-d` is optional now. By default Sortmerna will create `kvdb` directory for the Key-Value datastore in the User Home directory. If the `kvdb` directory already exists, the program will prompt the user to empty it, or to provide a different directory using `-d` option
- [[#172](https://github.com/sortmerna/sortmerna/issues/172)] - The multi-threading options are moved to a Developer category. They are still listed in the Help message
- 
## [[v3.0.0](https://github.com/sortmerna/sortmerna/releases/tag/v3.0.0)] - 2018-10-30

### Added

- [KSEQ](http://lh3lh3.users.sourceforge.net/parsefastq.shtml) library
- support for compressed reads input
- support for both mmap and geneic buffer input

### Changed

- Modified architecture to use a concurrent queue for holding Reads. The reads are pushed to the queue by the Reads file Reader, and processed in multiple Processor threads. This allows for using Const memory, only for the Reads being currently processed(aligned). At any time the number of Reads im memory is not more than the number of processing threads
- Modified architecture to store alignment results in a key-value database (RocksDB)
- Separated Reads processing into independent stages: Alignment, Post-Processing, Report generation
- Transfered build system to CMake
- Ported software to Windows
- Changed directory structure to better organise the code, and separated into projects suited for using CMake
- Modified python code to comply with version 3.5
- updated README to include info on SortMeRNA forum
- updated python tests to use Python 3 & scikit-bio version 0.5.0.dev0
- clean up compiler warnings

### Removed 

- Using standard C++ threads instead of OpenMP library. Sortmerna is multi-threaded by design now
- removal of Memory Mapping when processing the Reads file. The Reads are consumed from a stream, put on a queue and immediately processed, so that only a handful of reads is kept in the memory at any time

## [[v2.1b](https://github.com/sortmerna/sortmerna/releases/tag/2.1b)] - 2016-03-04

### Fixed 

- [[#105](https://github.com/sortmerna/sortmerna/pull/105)] - fix issue regarding duplications in include_HEADERS for Galaxy

## [[v2.1](https://github.com/sortmerna/sortmerna/releases/tag/2.1)] - 2016-02-01

### Added

- SILVA associated taxonomy string to representative databases

### Changed

- modified option `--blast INT` to `--blast STRING` to support more fields in BLAST-like tabular output (ex. `--blast '1 cigar qstrand'`)

### Fixed

- [[#70](https://github.com/sortmerna/sortmerna/issues/70)] - fixed issue (-m parameter for sortmerna not working for values greater than 4096); problem was related to large file support (`-D_FILE_OFFSET_BITS=64` flag added to compilation)
- fixed bug that causes incorrect CIGAR string when reference length < read length and based on the LCS, read hangs off end reference (alignment length should be computed based on this setup)
- [[#48](https://github.com/sortmerna/sortmerna/issues/48)] - ixed issue that failed to check all possible temporary directories for writing tmp files
- [[#72](https://github.com/sortmerna/sortmerna/issues/72)] - fixed issue that caused incorrect CIGAR string when reference length < read length and based on the LCS, read hangs off end reference (alignment length should be computed based on this setup)

## [[v2.0](https://github.com/sortmerna/sortmerna/releases/tag/2.0)] - 2014-10-30

### Added

- [affects Installation] added script `build.sh` to call configure, touch commands and make in order to avoid timestamp issues when cloning this repository
- OTU-picking extensions added for closed-reference clustering compatible with QIIME’s v1.9 `pick_otus.py`, `pick_closed_reference_otus.py` and `pick_open_reference_otus.py` scripts
- tests added

### Changed

- representative SILVA databases updated to version 119 (for filtering rRNA from metatranscriptomic data)
- [affects FASTQ paired reads] edited code for splitting FASTQ file
- [affects OTU-picking using `--otu_map --best INT` where INT > 1] changed the read_hits_align_info structure from `map<uint32_t,pair<uint16_t,s_align*> >` to `map<uint32_t, triple_s >` where the structure `triple_s` holds two `uint16_t` variables and an `s_align*` variable. This allows to store an additional integer for giving the index in `s_align` array of the highest scoring alignment. This is necessary if `--otu_map` and `--best [INT]` options are used where INT > 1 since different OTU ids can be used in the OTU-map when multiple alignments score equally as well. To illustrate an example,
	- (a) `--best 1`, the `s_align` array is size 1 and only the single best alignment is stored, being the first encountered alignment if multiple alignments of equal score are found.
	- (b) `--best 4`, the `s_align` array is size 4. Assume the first 2 alignments score 144 (occupying the first 2 slots on `s_align` array) and the next 3 alignments score 197 (2 of these alignments will occupy the final 2 slots, where the 3rd alignment will overwrite the first slot holding 144). We will have a situation like: 197, 144, 197, 197. In order to follow the same principle of (a) where the first encountered alignment of the highest score is output, we need to know that this alignment was in slot 3 (not slot 1).
- [affects multiple split databases] moved the declaration + initialization/deletion of `int32_t *best_x` from outside of `for each index_num` loop to inside the `for each index_part` loop. This is required to maintain similar results when using 1 index part (all database indexed as one part) vs. multiple index parts. The difference occurs because of the following situation:

(a) Database indexed as 1 part:

| candidate sequence | #seed hits |
| --------------------- | ------------- |
| ref1 | 10 |
| ref3 [correct reference] | 9 | 
| ref2 | 8 |
| ref4 | 8 |

(b) Database indexed as 2 parts:

part 1:
		
| candidate sequence | #seed hits |
| -------------- | ------------ |
|ref1 | 10 |
|ref2 | 8 |

part 2:
		
| candidate sequence | #seed hits |
| ------------------- | --------------- |
| ref3 [correct reference] | 9 | 
| ref4 | 8 |
 
If `min_lis_gv` = 2 (best_x[readn] = 2), then ref1 and ref3 will be analyzed in (a) before best_x[readn] = 0 and we stop analysis. However, in (b), if `min_lis_gv` = 2 outside of `for each index_num` loop, only ref1 and ref2 will be analyzed in part 1 at which point best_x[readn] = 0 and sequences in part 2 will not be analyzed. By initializing best_x[readn] = 2 at the start of each index_part, then ref1/ref2 in part1 will be analyzed and ref3/ref4 in part 2, where the correct reference sequence ref3 will be analyzed.

### Fixed

- [affects FASTQ paired reads] fixed the bug regarding `--paired_in` and `--paired_out` output

## [v1.99] - 2014-03-11

### Added 

- SortMeRNA can now perform alignments and output SAM and Blast-like formats (using the SSW Library, see [Zhao M. et al., "SSW Library: An SIMD Smith-Waterman C/C++ Library for Use in Genomic Applications", PLOS ONE, 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138))

### Changed

- Indexing data structures re-written and more optimized for space, requiring considerably less memory than previous versions (integration of the C Minimal Perfect Hashing (CMPH) library, see [http://cmph.sourceforge.net](http://cmph.sourceforge.net))
- Multiple indexes can now be constructed in one command (now indexdb_rna rather than buildtrie)

### Fixed 

- Issues with FASTA/Q output files resolved (thanks to Ali May)

### Removed 

- The `$SORTMERNADIR` environmental variable no longer used

## [v1.9] - 2013-08-30

### Changed 

- updated merge-paired-reads.sh to work on a cluster (thanks to Nicolas Delhomme)

### Fixed

- fixed a bug for naming output log file (thanks to Shaman Narayanasamy)
- the paths for binaries sortmerna and buildtrie have been corrected to work with `make install` for installation directories other than the default `/usr/local`

## [v1.8] - 2013-05-13

### Changed

- modified `merge.sh` and `unmerge-paired-reads.sh`

### Fixed

- fixed a bug to detect last (rRNA) read in fastq files

## [v1.7] - 2013-04-05

### Added

- added `merge_paired_reads.sh` for forward-reverse paired-end reads (see the user manual v-1.7, section 4.2.4)

### Changed

- changed to the usual 'configure, make, make install' (see the user manual v-1.7)

### Fixed

- fixed an integer overflow for mmap calculation for 32-bit systems 

## [v1.6] - 2013-02-26

### Changed

- for taxonomical analysis, the sequence tags in the rRNA databases now follow the format: `>[accession] [taxonomy] [length]`
- changed sysconf library to sysctl for Mac OS 

## [v1.5] - 2013-02-15

### Added 

- opion `-m` for specifying the amount of memory for loading reads
- local timestamp added to `--log` statistics file

### Changed

- reads of length <L (default L=18) are automatically considered as non-rRNA
- SortMeRNA User Manual updated

### Fixed

- error for output of paired reads >1GB resolved

## [v1.4] - 2013-02-06

### Added

- support for Illumina or 454 reads up to 5000 nucleotides
- AUTHORS file added to SortMeRNA directory

## [v1.3] - 2013-01-24

### Added

- support for paired-end reads

### Changed

- I/O file checks modified
- Makefile updated

### Fixed

- L=20 error messages resolved

## [v1.2] - 2013-01-07

### Added

- support for input directory without suggested path (assumes current)

## [v1.1] - 2012-12-20

### Added

- support for input files wthout extensions

## [v1.0] - 2012-10-15

SortMeRNA v1.0 released
