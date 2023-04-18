# SortMeRNA: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[v3.0.0](https://github.com/sortmerna/sortmerna/releases/tag/v3.0.0)] - 2018-02-28

### Changed
- Modified architecture to use C++ threads instead of OpenMP
- Modified architecture to use a concurrent queue for holding Reads. The reads are pushed to the queue by the Reads file Reader, and processed in multiple Processor threads. This allows for using Const memory, only for the Reads being currently processed(aligned). At any time the number of Reads im memory is not more than the number of processing threads. 
- Modified architecture to store alignment results in a key-value database (RocksDB).
- Separated Reads processing into independent stages: Alignment, Post-Processing, Report generation

### Removed 

- OpenMP
- Memory mapping

## [[]()] - 2017-07-14

### Changed 

- Transfered build system to CMake
- Ported software to Windows
- Changed directory structure to better organise the code, and separated into projects suited for using CMake
- Modified python code to comply with version 3.5

## [[]()] - 2017-05-18

### Changed

- updated README to include info on SortMeRNA forum
- updated python tests to use Python 3 & scikit-bio version 0.5.0.dev0

## [[]()] - 2016-02-17

### Added

- [KSEQ](http://lh3lh3.users.sourceforge.net/parsefastq.shtml) library, allow compressed reads input
- support for both mmap and geneic buffer input

### Changed

- clean up compiler warnings

## [[v2.1b](https://github.com/sortmerna/sortmerna/releases/tag/2.1b)] - 2016-03-03

### Fixed 

- fix issue regarding duplications in include_HEADERS for Galaxy (see pull request 105)

## [[v2.1](https://github.com/sortmerna/sortmerna/releases/tag/2.1)] - 2016-02-01

### Added

- SILVA associated taxonomy string to representative databases

### Changed

- modified option --blast INT to --blast STRING to support more fields in BLAST-like tabular output (ex. --blast '1 cigar qstrand')

### Fixed

- [[#70](https://github.com/sortmerna/sortmerna/issues/70)] - update for FASTQ reads
- [[#70](https://github.com/sortmerna/sortmerna/issues/70)] - fixed issue (-m parameter for sortmerna not working for values greater than 4096); problem was related to large file support (-D_FILE_OFFSET_BITS=64 flag added to compilation)
- fixed bug that causes incorrect CIGAR string when reference length < read length and based on the LCS, read hangs off end reference (alignment length should be computed based on this setup)
- [[#48](https://github.com/sortmerna/sortmerna/issues/48)] - resolved issue
- [[#72](https://github.com/sortmerna/sortmerna/issues/72)] - resolved issue

## [[v2.0](https://github.com/sortmerna/sortmerna/releases/tag/2.0)] - 2014-08-30 30 October 2014

### Added

- [affects Installation] added script `build.sh` to call configure, touch commands and make in order to avoid timestamp issues when cloning this repository
- OTU-picking extensions added for closed-reference clustering compatible with QIIME’s v1.9 pick_otus.py, pick_closed_reference_otus.py and pick_open_reference_otus.py scripts

### Changed

- representative SILVA databases updated to version 119 (for filtering rRNA from metatranscriptomic data)
- [affects FASTQ paired reads] edited code for splitting FASTQ file
- [affects OTU-picking using "--otu_map --best INT" where INT > 1] changed the read_hits_align_info structure from map<uint32_t,pair<uint16_t,s_align*> > to map<uint32_t, triple_s > where the structure triple_s holds two uint16_t variables and an s_align* variable. This allows to store an additional integer for giving the index in s_align array of the highest scoring alignment. This is necessary if —otu_map and —best [INT] options are used where INT > 1 since different OTU ids can be used in the OTU-map when multiple alignments score equally as well. To illustrate an example,
	- (a) —best 1, the s_align array is size 1 and only the single best alignment is stored, being the first encountered alignment if multiple alignments of equal score are found.
	- (b) —best 4, the s_align array is size 4. Assume the first 2 alignments score 144 (occupying the first 2 slots on s_align array) and the next 3 alignments score 197 (2 of these alignments will occupy the final 2 slots, where the 3rd alignment will overwrite the first slot holding 144). We will have a situation like: 197, 144, 197, 197. In order to follow the same principle of (a) where the first encountered alignment of the highest score is output, we need to know that this alignment was in slot 3 (not slot 1).
- [affects multiple split databases] moved the declaration + initialization/deletion of int32_t *best_x from outside of "for each index_num loop" to inside the "for each index_part loop". This is required to maintain similar results when using 1 index part (all database indexed as one part) vs. multiple index parts. The difference occurs because of the following situation:

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
 
If min_lis_gv = 2 (best_x[readn] = 2), then ref1 and ref3 will be analyzed in (a) before best_x[readn] = 0 and we stop analysis. However, in (b), if min_lis_gv = 2 outside of "for each index_num loop", only ref1 and ref2 will be analyzed in part 1 at which point best_x[readn] = 0 and sequences in part 2 will not be analyzed. By initializing best_x[readn] = 2 at the start of each index_part, then ref1/ref2 in part1 will be analyzed and ref3/ref4 in part 2, where the correct reference sequence ref3 will be analyzed.

### Fixed

- [affects FASTQ paired reads] fixed the bug regarding --paired_in and --paired_out output, tests added

## [v1.99] - 2014-03-11

### Added 

- SortMeRNA can now perform alignments and output SAM and Blast-like formats (using the SSW Library, see [Zhao M. et al., "SSW Library: An SIMD Smith-Waterman C/C++ Library for Use in Genomic Applications", PLOS ONE, 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138))

### Changed

- Indexing data structures re-written and more optimized for space, requiring considerably less memory than previous versions (integration of the C Minimal Perfect Hashing (CMPH) library, see [http://cmph.sourceforge.net](http://cmph.sourceforge.net))
- Multiple indexes can now be constructed in one command (now indexdb_rna rather than buildtrie)

### Fixed 

- Issues with FASTA/Q output files resolved (thanks to Ali May)

### Removed 

- The $SORTMERNADIR environmental variable no longer used

## [v1.9] - 2013-08-30

### Changed 

- updated merge-paired-reads.sh to work on a cluster (thanks to Nicolas Delhomme)

### Fixed

- fixed a bug for naming output log file (thanks to Shaman Narayanasamy)
- the paths for binaries sortmerna and buildtrie have been corrected to work with `make install` for installation directories other than the default /usr/local

## [v1.8] - 2013-05-13

### Fixed

- fixed a bug to detect last (rRNA) read in fastq files, modified merge/unmerge-paired-reads.sh 

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

- opion -m for specifyinf the amount of memory for loading reads
- local timestamp added to --log statistics file

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

## [v1.0] - 2012-08-15

SortMeRNA v1.0 released
