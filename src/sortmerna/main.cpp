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
#include "readfeed.hpp"
#include "processor.hpp"


/*
  main entry of the sortmerna application
 */
int main(int argc, char** argv)
{
	bool dryrun = false;
	Runopts opts(argc, argv, dryrun);

	INFO("Running command:\n", opts.cmdline);

	if (opts.is_cmd) {
		CmdSession cmd;
		cmd.run(opts);
	}
	else
	{
		// init common objects
		KeyValueDatabase kvdb(opts.kvdbdir.string());
		Readfeed readfeed(opts.feed_type, opts.readfiles, opts.num_proc_thread, opts.readb_dir);
		if (opts.feed_type == FEED_TYPE::SPLIT_READS)
			readfeed.split();
		Readstats readstats(readfeed.num_reads_tot, readfeed.length_all, readfeed.min_read_len, readfeed.max_read_len, kvdb, opts);
		Index index(opts); // reference index DB
		Output output(opts, readstats);

		switch (opts.alirep)
		{
		case Runopts::ALIGN_REPORT::align:
			align(readfeed, readstats, index, kvdb, output, opts);
			break;
		case Runopts::ALIGN_REPORT::postproc:
			postProcess(readfeed, readstats, kvdb, output, opts);
			break;
		case Runopts::ALIGN_REPORT::report:
			generateReports(readfeed, readstats, kvdb, output, opts);
			break;
		case Runopts::ALIGN_REPORT::alipost:
			align(readfeed, readstats, index, kvdb, output, opts);
			postProcess(readfeed, readstats, kvdb, output, opts);
			break;
		case Runopts::ALIGN_REPORT::all:
			align(readfeed, readstats, index, kvdb, output, opts);
			postProcess(readfeed, readstats, kvdb, output, opts);
			generateReports(readfeed, readstats, kvdb, output, opts);
			break;
		}
	}
	return 0;
}//~main()
