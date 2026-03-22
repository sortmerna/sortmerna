/*
@copyright 2016-2026 Clarity Genomics BVBA
@copyright 2012-2016 Bonsai Bioinformatics Research Group
@copyright 2014-2016 Knight Lab, Department of Pediatrics, UCSD, La Jolla

@parblock
SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA

This is a free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SortMeRNA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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
              biocodz          biocodz@protonmail.com
*/

/**
* writer.cpp   created: Nov 26, 2017 Sun
*/
#include <iomanip>
#include <sstream>
#include <chrono>
#include <iostream>

#include "writer.hpp"


// write read alignment results to disk using e.g. RocksDB
void Writer::write()
{
	{
		std::stringstream ss;
		ss << STAMP << "Writer " << id << " thread " << std::this_thread::get_id() << " started" << std::endl;
		std::cout << ss.str();
	}

	auto t = std::chrono::high_resolution_clock::now();
	int numPopped = 0;
	std::size_t num_aligned = 0; // num reads with 'read.hit = true' i.e. passing E-value threshold
	for (;;) 
	{
		Read read = writeQueue.pop();
		if (read.isEmpty)
		{
			if (writeQueue.getPushers() == 0)
				break; // no more records in the queue and no pushers => stop processing

			if (!read.isValid) 
				continue;
		}
		++numPopped;
		//std::string matchResultsStr = read.matchesToJson();
		std::string readstr = read.toBinString();
		if (!opts.is_dbg_put_kvdb && readstr.size() > 0)
		{
			if (read.is_hit) ++num_aligned;
			kvdb.put(read.id, readstr);
		}
	}
	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;

	{
		std::stringstream ss;
		ss << STAMP << std::setprecision(2) << std::fixed << id << " thread " << std::this_thread::get_id()
			<< " done. Elapsed time: " << elapsed.count() << " s Reads written: " << numPopped 
			<< " Num aligned reads (passing E-value):" << num_aligned << std::endl;
		std::cout << ss.str();
	}
} // Writer::write