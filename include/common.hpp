/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is a free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SortMeRNA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 @endparblock

 @contributors Jenya Kopylova   jenya.kopylov@gmail.com
			   Laurent Noé      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mikaël Salson    mikael.salson@lifl.fr
			   Hélène Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

 /** @file common.hpp */

#pragma once

#include <string>
#include <sstream>
#include <iostream> // std::cout

#include <sys/time.h>
#include "config.h"
#if defined(_WIN32)
#  include "windows.h"
#  include "psapi.h"
#else
#  include <fstream>
#endif

const char FASTA_HEADER_START = '>';
const char FASTQ_HEADER_START = '@';
const std::string FWD = "FWD";
const std::string REV = "REV";

enum class BIO_FORMAT : unsigned { FASTQ = 0, FASTA = 1 };
enum class ZIP_FORMAT : unsigned { GZIP = 0, ZLIB = 1, FLAT = 2, XPRESS = 3 };
enum class FEED_TYPE : unsigned { SPLIT_READS = 0, LOCKLESS = 1, MAX = LOCKLESS };
enum class BlastFormat { TABULAR, REGULAR}; // format of the Blast output

/*! @brief Map nucleotides to integers.
Ambiguous letters map to 4.
{A/a,C/c,G/g,T/t,U/u} = {0,1,2,3,3} respectively.
*/
const char nt_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// see old paralleltraversal.cpp::format_rev
const char rc_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4,
4, 4, 4, 4, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

const char nt_map[5] = { 'A', 'C', 'G', 'T', 'N' };

const char complement[5] = { 3, 2, 1, 0, 4 }; // A <-> T, C <-> G, N <-> N

extern timeval t;

/*! @brief Macro for timing */
#define TIME(x) gettimeofday(&t, NULL); x = t.tv_sec + (t.tv_usec/1000000.0);

/*! @brief start color text red */
#if defined(_WIN32)
#  define RED    ""
#  define GREEN  ""
#  define YELLOW ""
#  define BLUE   ""
#  define BOLD   ""
#  define UNDL   ""
#  define COLOFF ""
const char DELIM = ';';
#else
#  define RED    "\033[0;31m"
#  define GREEN  "\033[0;32m"
#  define YELLOW "\033[0;33m"
#  define BLUE   "\033[0;34m"
#  define BOLD   "\033[1m"
#  define UNDL   "\033[4m" // underline
#  define COLOFF "\033[0m" // color off
const char DELIM = ':';
#endif


/*! @brief Maximum length of input reads
    (not limited to this length algorithmically)
*/
#define MAX_READ_LEN 30000

//#define LOCKQUEUE // Lock queue with mutexes
#define CONCURRENTQUEUE // lockless queue
#define STAMP  "[" << __func__ << ":" << __LINE__ << "] "
#define STAMPL "[" << __FILE__ << ":" << __func__ ":" << __LINE__ << "] "


template<typename ...Args>
static inline std::string fold_to_string(Args&&... args) {
    std::stringstream ss;
    (ss << ... << args);
    return ss.str();
}

// https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
#if defined(_WIN32)
static inline size_t get_memory() {
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    return pmc.PrivateUsage; // process Commit memory (== Resource Monitor:Commit)
}
#else
/*
 * amount of memory that have been mapped into the process' address space
 * cat /proc/self/status | grep 'VmRSS:'
 *   'VmRSS:       336 kB'
 */
static inline size_t get_memory() {
    std::ifstream ifs("/proc/self/status", std::ios::binary);
    size_t mem = 0;
    if (ifs.is_open()) {
        for (std::string line; std::getline(ifs, line); ) {
            unsigned pos_start = 0;
            unsigned pos_end = 0;
            if (line.find("VmRSS:") != std::string::npos) {
                for (unsigned i = 0; i < line.length(); ++i) {
                    if (pos_start == 0) {
                        if (std::isdigit(line[i]))
                            pos_start = i;
                    }
                    else {
                        if (pos_end == 0 && !std::isdigit(line[i])) {
                            pos_end = i;
                            break;
                        }
                    }
                }
                if (pos_end > pos_start) {
                    std::string mstr = line.substr(pos_start, pos_end + 1 - pos_start);
                    mem = std::atoi(mstr.data());
                }
            }
        }
    }
    return mem;
}
#endif

#define INFO(...) \
	{\
		std::stringstream ss; \
		ss << STAMP << fold_to_string(__VA_ARGS__) << std::endl; \
		std::cout << ss.str(); \
	}

// no end line
#define INFO_NE(...) \
	{\
		std::stringstream ss; \
		ss << STAMP << fold_to_string(__VA_ARGS__); \
		std::cout << ss.str(); \
	}

// No Stamp, no endl
#define INFO_NS(...) \
	{\
		std::stringstream ss; \
		ss << fold_to_string(__VA_ARGS__); \
		std::cout << ss.str(); \
	}

#define INFO_MEM(...) \
	{\
		std::stringstream ss; \
		ss << STAMP << fold_to_string(__VA_ARGS__) << " Memory KB: " << (get_memory() >> 10) << std::endl; \
		std::cout << ss.str();\
	}

#define WARN(...) \
	{\
		std::stringstream ss; \
		ss << '\n' << STAMP << YELLOW << "WARNING" << COLOFF << ": " << fold_to_string(__VA_ARGS__) << std::endl; \
		std::cout << ss.str();\
	}

#define ERR(...) \
	{\
		std::stringstream ss; \
		ss << '\n' << STAMP << RED << "ERROR" << COLOFF << ": " << fold_to_string(__VA_ARGS__) << std::endl; \
		std::cerr << ss.str();\
	}

#define PRN_MEM(msg) \
	{\
		std::stringstream ss; \
		ss << STAMP << msg << " Memory KB: " << (get_memory() >> 10) << std::endl; \
		std::cout << ss.str();\
	}

#define PRN_MEM_TIME(msg, time) \
    {\
		std::stringstream ss; \
		ss << STAMP << msg << " Memory KB: " << (get_memory() >> 10) << " Elapsed sec: " << time << std::endl; \
		std::cout << ss.str();\
    }
//~EOF