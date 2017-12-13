/**
* FILE: readstats.hpp
* Created: Nov 17, 2017 Fri
*
* Read file statistics
*/

// standard
#include <chrono>
#include <algorithm> // remove_if
#include <iomanip> // output formatting
#include <locale> // isspace

// 3rd party
#include "zlib.h"
#include "kseq_load.hpp"

// SMR
#include "readstats.hpp"


/**
 * TODO: get rid of kseq, Use stream instead of file. Use Reader for reading.
 *
 * prototype: paralleltraversal.cpp:compute_read_stats
 */
void Readstats::calculate2()
{
#ifdef HAVE_LIBZ
	// Count total number of reads and their combined length
	// (if ZLIB is supported)
	gzFile fp = gzopen(opts.readsfile.c_str(), "r"); // inputreads
	kseq_t *seq = kseq_init(fp);
#else
	// Count total number of reads and their combined length
	// (if ZLIB is not supported)
	FILE* fp = fopen(inputreads, "rb");
	kseq_t *seq = kseq_init(fileno(fp));
#endif
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		full_read_main += seq->seq.l;
		number_total_read++;
		// compute size of all reads to store in memory
		// + 7 (4 possible new lines, 2 symbols > or @ and +, space for comment)
		if (!map_size_set_gv)
			full_file_size += (seq->name.l + seq->comment.l + seq->seq.l + seq->qual.l + 7);
	}
	if (l == -2)
	{
		fprintf(stderr, "  %sERROR%s: Line %d: %s could not read reads file - %s\n\n",
			startColor, "\033[0m", __LINE__, __FILE__, strerror(errno));
		exit(EXIT_FAILURE);
	}
	kseq_destroy(seq);
#ifdef HAVE_LIBZ
	gzclose(fp);
#else
	fclose(fp);
#endif
} // ~Readstats::calculate

void Readstats::calculate()
{
	std::ifstream ifs(opts.readsfile, std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open()) {
		std::cout << "Failed to open Reads file: " << opts.readsfile << "\n";
		exit(EXIT_FAILURE);
	}
	else
	{
		std::string line;
		bool isFastq = true;

		auto t = std::chrono::high_resolution_clock::now();

		for (int count = 0; std::getline(ifs, line); ) // count lines in One record
		{
			// skip empty line
			if (!line.empty())
			{
				// right-trim whitespace in place (removes '\r' too)
				line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
				//line.erase(std::remove_if(begin(line), end(line), [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), end(line));
				// fastq: 0(header), 1(seq), 2(+), 3(quality)
				// fasta: 0(header), 1(seq)
				if (line[0] == FASTA_HEADER_START || line[0] == FASTQ_HEADER_START)
				{
					count = 0; // record start
					isFastq = (line[0] == FASTQ_HEADER_START);
				} // ~if header line
				else {
					++count;
					if (count > 3) {
						std::cout << "Unexpected number of lines: " << count << " in a single Read. Total reads processed: " 
							<< number_total_read <<" Exiting...\n";
						exit(EXIT_FAILURE);
					}
					if (isFastq && (line[0] == '+' || count == 3)) continue; // fastq.quality

					++number_total_read;
					full_read_main += line.length();
				}
			}
		} // ~for getline

		std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
		std::cout << std::setprecision(2) << std::fixed;
		std::cout << "Readstats::calculate done. Elapsed time: " << elapsed.count() << " sec. Reads processed: " << number_total_read << std::endl;
	}
	ifs.close();
} // ~Readstats::calculate

bool Readstats::check_file_format()
{
	bool exit_early = false;
#ifdef HAVE_LIBZ
	// Check file format (if ZLIB supported)
	gzFile fp = gzopen(opts.readsfile.c_str(), "r");
	kseq_t *seq = kseq_init(fp);
#else
	FILE* fp = fopen(opts.readsfile.c_str(), "rb");
	kseq_t *seq = kseq_init(fileno(fp));
#endif
	int l;
	if ((l = kseq_read(seq)) >= 0)
		filesig = seq->last_char;
	else
	{
		fprintf(stderr, "  %sERROR%s: Line %d: %s unrecognized file format or empty file %s\n\n",
			startColor, "\033[0m", __LINE__, __FILE__, opts.readsfile.c_str());
		exit_early = true;
	}
	kseq_destroy(seq);
#ifdef HAVE_LIBZ
	gzclose(fp);
#else
	fclose(fp);
#endif
	return exit_early;
} // ~Readstats::check_file_format

// determine the suffix (fasta, fastq, ...) of aligned strings
void Readstats::calcSuffix()
{
	const std::string suff = opts.readsfile.substr(opts.readsfile.rfind('.') + 1);
	if (suff.length() > 0 && !opts.have_reads_gz)
		suffix.assign(suff);
	else if (filesig == '>')
		suffix.assign("fasta");
	else
		suffix.assign("fastq");
}