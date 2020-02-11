/**
 * @file main.cpp
 * @brief File containing the main function and argument parsing.
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2016-20 Clarity Genomics BVBA
 * @copyright 2012-16 Bonsai Bioinformatics Research Group
 * @copyright 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
 * @copyright 2016-19 Clarity Genomics BVBA
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
 * along with SortMeRNA.  If not, see <http://www.gnu.org/licenses/>.
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

#include <iostream>
#include <filesystem>

#include "options.hpp"
#include "paralleltraversal.hpp"
#include "output.hpp"
#include "readstats.hpp"
#include "cmd.hpp"
#include "kvdb.hpp"
#include "index.hpp"
#include "indexdb.hpp"

namespace fs = std::filesystem;

// forward
void postProcess(Runopts & opts, Readstats & readstats, Output & output, KeyValueDatabase &kvdb); // processor.cpp
void setup_workspace(Runopts & opts);

/**
 * Setup the directory structure necessary for the sortmerna process 
 */
void setup_workspace(Runopts &opts)
{
	fs::create_directories(opts.workdir + "/" + opts.KVDB_DIR);
	fs::create_directory(opts.workdir + "/" + opts.IDX_DIR);
	fs::create_directory(opts.workdir + "/" + opts.OUT_DIR);
}

/*! @fn main()
	@brief main function, parses command line arguments and launches the processing
	@param int argc
	@param char** argv
	@return none
 */
int main(int argc, char** argv)
{
	bool dryrun = false;
	Runopts opts(argc, argv, dryrun);

	std::cout << STAMP << "Running command:\n" << opts.cmdline << std::endl;

	setup_workspace(opts);
	Index index(opts); // reference index DB
	KeyValueDatabase kvdb(opts.kvdbPath);

	if (opts.is_cmd) {
		CmdSession cmd;
		cmd.run(opts);
	}
	else
	{
		Readstats readstats(opts, kvdb);
		Output output(opts, readstats);

		switch (opts.alirep)
		{
		case Runopts::ALIGN_REPORT::align:
			align(opts, readstats, output, index, kvdb);
			break;
		case Runopts::ALIGN_REPORT::postproc:
			postProcess(opts, readstats, output, kvdb);
			break;
		case Runopts::ALIGN_REPORT::report:
			generateReports(opts, readstats, output, kvdb);
			break;
		case Runopts::ALIGN_REPORT::alipost:
			align(opts, readstats, output, index, kvdb);
			postProcess(opts, readstats, output, kvdb);
			break;
		case Runopts::ALIGN_REPORT::all:
			align(opts, readstats, output, index, kvdb);
			postProcess(opts, readstats, output, kvdb);
			generateReports(opts, readstats, output, kvdb);
			break;
		}
	}

	return 0;
}//~main()
