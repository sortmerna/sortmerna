#pragma once
 /**
 * @FILE output.hpp
 * Created: Nov 06, 2017 Mon
 * @brief Class for handling SMR output
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2016-20 Clarity Genomics BVBA
 * @copyright 2012-16 Bonsai Bioinformatics Research Group
 * @copyright 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
 *
 * SortMeRNA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SortMeRNA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 * @endparblock
 *
 * @contributors Jenya Kopylova, jenya.kopylov@gmail.com
 *               Laurent Noé, laurent.noe@lifl.fr
 *               Pierre Pericard, pierre.pericard@lifl.fr
 *               Daniel McDonald, wasade@gmail.com
 *               Mikaël Salson, mikael.salson@lifl.fr
 *               Hélène Touzet, helene.touzet@lifl.fr
 *               Rob Knight, robknight@ucsd.edu
 */
#include <fstream>
#include <stdint.h>
#include <string>
#include <vector>

#include "common.hpp"
#include "kvdb.hpp"

// forward
struct Index;
class References;
class Read;
class Refstats;
struct Readstats;
struct Runopts;

/**
 * Summary report (log) data structure
 */
class Summary {
public:
	// vars
	bool is_de_novo_otu;
	bool is_otumapout;
	std::string cmd;
	std::string timestamp;
	std::string pid_str;
	uint64_t total_reads;
	uint64_t total_reads_denovo_clustering;
	uint64_t total_reads_mapped;
	uint64_t total_mapped_sw_id_cov;
	uint32_t min_read_len;
	uint32_t max_read_len;
	uint64_t all_reads_len;
	size_t total_otu;
	std::vector<std::pair<std::string, float>> db_matches;

	// methods
	Summary();
	std::string to_string(Runopts& opts, Refstats& refstats);
};

class Output {
public:
	// output streams
	std::vector<std::ofstream> aligned_os; // fasta/q    20200127 2 files if 'out2', 1 file otherwise
	std::vector<std::ofstream> other_os; // fasta/q non-aligned   20200127 2 files if 'out2', 1 file otherwise
	std::ofstream sam_os; // SAM
	std::ofstream blast_os; // BLAST
	std::ofstream log_os;
	std::ofstream denovo_os;
	std::ofstream biom_os;

	// output file names
	std::vector<std::string> aligned_f; // 20200127 2 files if 'out2'
	std::vector<std::string> other_f; // 20200127 2 files if 'out2'
	std::string sam_f;
	std::string blast_f;
	std::string log_f;
	std::string denovo_otus_f;
	std::string otumap_f;
	std::string biom_f;

	Summary summary;

	Output(Runopts& opts, Readstats& readstats);
	~Output();

	void report_blast(
		Runopts & opts,
		Refstats & refstats,
		References & refs,
		Read & read
	);

	void report_sam(
		Runopts & opts,
		References & refs,
		Read & read
	);

	void writeSamHeader(Runopts & opts);

	void report_fasta(Runopts & opts, std::vector<Read> &reads);
	void report_denovo(Runopts & opts, std::vector<Read> &reads);
	void report_biom();
	void writeLog(Runopts &opts, Refstats &refstats, Readstats &readstats);

	void openfiles(Runopts & opts);
	void closefiles();

private:
	void init(Runopts & opts, Readstats & readstats);
	void write_a_read(std::ofstream& strm, Read& read);

}; // ~class Output


void generateReports(Runopts & opts, Readstats & readstats, Output & output, KeyValueDatabase &kvdb);