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

// forward
struct Index;
class References;
class Read;
class Refstats;
struct Readstats;
struct Runopts;
class KeyValueDatabase;
class Readfeed;

class Output {
public:
	// output streams
	std::vector<std::ofstream> ofs_aligned; // aligned fasta/q 
	std::vector<std::ofstream> ofs_other; // non-aligned fasta/q 
	std::vector<std::ofstream> ofs_blast; // BLAST
	std::ofstream ofs_sam; // SAM
	std::ofstream ofs_denovo;
	std::ofstream ofs_biom;

	// output file names
	std::vector<std::string> f_aligned;
	std::vector<std::string> f_other;
	std::vector<std::string> f_blast;
	std::string f_sam;
	std::string f_denovo_otus;
	std::string f_biom;

	int num_out; // number of output files

	Output(Readfeed& readfeed, Runopts& opts, Readstats& readstats);
	~Output();

	void report_fasta(int id, std::vector<Read>& reads, Runopts& opts);
	void report_blast(int id, Read& read, References& refs, Refstats& refstats, Runopts& opts);
	void report_sam(Runopts & opts, References & refs, Read & read);
	void writeSamHeader(Runopts & opts);
	void report_denovo(Runopts & opts, std::vector<Read> &reads);
	void report_biom();
	void openfiles(Runopts & opts);
	void closefiles();
	void calc_out_type(Runopts& opts);
	void set_num_out();

private:
	int out_type; //  = 0x00
	const int mask_1_file = 0x01;
	const int mask_paired = 0x02;
	const int mask_2_file = 0x04;
	const int mask_other  = 0x08;
	const int mask_paired_in = 0x10;
	const int mask_paired_out = 0x20;
	const int mask_out2 = 0x40;

private:
	void init(Readfeed& readfeed, Runopts& opts, Readstats& readstats);
	void init_fastx(Readfeed& readfeed, Runopts& opts);
	void init_blast(Readfeed& readfeed, Runopts& opts);
	void init_sam(Runopts& opts);
	void init_denovo_otu(Readfeed& readfeed, Runopts& opts);
	void write_a_read(std::ofstream& strm, Read& read);

}; // ~class Output