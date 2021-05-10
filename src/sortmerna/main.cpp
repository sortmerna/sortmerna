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

/*
 * @file main.cpp
 * @brief the main function and argument parsing.
 */

#include <iostream>
#include <filesystem>

#include "options.hpp"
#include "paralleltraversal.hpp"
#include "readstats.hpp"
#include "cmd.hpp"
#include "kvdb.hpp"
#include "index.hpp"
#include "indexdb.hpp"
#include "readfeed.hpp"
#include "processor.hpp"
#include "summary.hpp"
#include "output.hpp"
#include "otumap.h"
#include "refstats.hpp"


/*
*  main entry of the sortmerna application
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
		Index index(opts); // reference index DB
		if (Runopts::ALIGN_REPORT::index_only == opts.alirep) {
			INFO("Only performed indexing as '", OPT_INDEX, "' = 1 was specified");
			return 0;
		}

		// init common objects
		KeyValueDatabase kvdb(opts.kvdbdir.string());
		Readfeed readfeed(opts.feed_type, opts.readfiles, opts.num_proc_thread, opts.readb_dir, opts.is_paired);
		Readstats readstats(readfeed.num_reads_tot, readfeed.length_all, readfeed.min_read_len, readfeed.max_read_len, kvdb, opts);

		switch (opts.alirep)
		{
		case Runopts::ALIGN_REPORT::index_only:
			break;
		case Runopts::ALIGN_REPORT::align:
			align(readfeed, readstats, index, kvdb, opts);
			break;
		case Runopts::ALIGN_REPORT::summary:
			if (opts.is_otu_map || opts.is_denovo) denovo_stats(readfeed, readstats, kvdb, opts);
			if (opts.is_otu_map) fill_otu_map(readfeed, readstats, kvdb, opts);
			writeSummary(readstats, opts);
			break;
		case Runopts::ALIGN_REPORT::report:
			writeReports(readfeed, readstats, kvdb, opts);
			break;
		case Runopts::ALIGN_REPORT::alnsum:
			align(readfeed, readstats, index, kvdb, opts);
			if (opts.is_otu_map || opts.is_denovo) denovo_stats(readfeed, readstats, kvdb, opts);
			if (opts.is_otu_map) fill_otu_map(readfeed, readstats, kvdb, opts);
			writeSummary(readstats, opts);
			break;
		case Runopts::ALIGN_REPORT::all:
			align(readfeed, readstats, index, kvdb, opts);
			// TODO: combine processing otu map and reports to avoid double run through reads and refs (in this case only) 20201126
			if (opts.is_otu_map || opts.is_denovo) denovo_stats(readfeed, readstats, kvdb, opts);
			if (opts.is_otu_map) fill_otu_map(readfeed, readstats, kvdb, opts);
			writeSummary(readstats, opts);
			writeReports(readfeed, readstats, kvdb, opts);
			break;
		}
	}
	return 0;
}//~main()
