/*
@copyright 2016-2026 Clarity Genomics Inc
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
#include "restart.hpp"
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
		if (opts.task == Runopts::TASK::index_only) {
			INFO("Only performed indexing as '", OPT_TASK, "' = 5 was specified");
			return 0;
		}

		// init common objects
		KeyValueDatabase kvdb(opts.kvdbdir.string());
		// Probe kvdb for prior run state. Fresh kvdb -> stamp identity; existing
		// kvdb with matching fingerprints -> auto-resume; mismatched fingerprints
		// throw std::runtime_error so the user sees what changed.
		restart::State rstate = restart::probe_or_init(kvdb, opts);

		Readfeed readfeed(opts.feed_type, opts.readfiles, opts.num_proc_thread, opts.readb_dir, opts.is_paired);
		Readstats readstats(readfeed.num_reads_tot, readfeed.length_all, readfeed.min_read_len, readfeed.max_read_len, kvdb, opts);

		auto need_align = [&]() {
			if (!rstate.is_resume) return true;
			return rstate.phase <= restart::Phase::Align;
		};
		auto need_post_align = [&]() {
			if (!rstate.is_resume) return true;
			if (rstate.phase >= restart::Phase::Report) return false;
			return !restart::is_denovo_done(kvdb);
		};
		auto run_align = [&]() {
			if (need_align()) {
				align(readfeed, readstats, index, kvdb, opts, &rstate);
				restart::set_phase(kvdb, restart::Phase::PostAlign);
			} else {
				INFO("Restart: skipping alignment, phase already at ",
				     restart::phase_to_string(rstate.phase));
			}
		};
		auto run_post_align = [&]() {
			if (!need_post_align()) {
				INFO("Restart: skipping post-alignment (denovo_done present)");
				return;
			}
			restart::set_phase(kvdb, restart::Phase::Denovo);
			if (opts.is_otu_map || opts.is_denovo) denovo_stats(readfeed, readstats, kvdb, opts);
			if (opts.is_otu_map) fill_otu_map(readfeed, readstats, kvdb, opts);
			writeSummary(readstats, opts);
			restart::mark_denovo_done(kvdb);
		};
		auto run_reports = [&]() {
			restart::set_phase(kvdb, restart::Phase::Report);
			writeReports(readfeed, readstats, kvdb, opts);
			restart::set_phase(kvdb, restart::Phase::Done);
		};

		switch (opts.task)
		{
		case Runopts::TASK::index_only:
			break;
		case Runopts::TASK::align:
			run_align();
			break;
		case Runopts::TASK::summary:
			run_post_align();
			break;
		case Runopts::TASK::report:
			run_reports();
			break;
		case Runopts::TASK::align_summary:
			run_align();
			run_post_align();
			break;
		case Runopts::TASK::all:
			run_align();
			// TODO: combine processing otu map and reports to avoid double run through reads and refs (in this case only) 20201126
			run_post_align();
			run_reports();
			break;
		}
	}
	return 0;
}//~main()
