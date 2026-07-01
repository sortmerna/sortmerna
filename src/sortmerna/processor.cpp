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
 * FILE: processor.cpp
 * Created: Nov 26, 2017 Sun
 *
 * performs the alignment
 */

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <thread>
#include <cmath>

#include "processor.hpp"
#include "read.hpp"
#include "readfeed.hpp"
#include "index.hpp"
#include "references.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "options.hpp"
#include "restart.hpp"

// forward
void traverse(Runopts& opts, Index& index, References& refs, Readstats& readstats, Refstats& refstats, Read& read, bool isLastStrand);

namespace {

// Counter-write cadence inside align2. With 100 reads/write per slot, a worst-
// case crash loses ~100 records of work per slot. Tunable; not worth a CLI flag.
constexpr uint32_t COUNTER_FLUSH_INTERVAL = 100;

} // anonymous

/*
 * One alignment worker. Owns the slot(s) for thread `id`:
 *   - non-paired (num_sense==1):           slot = id
 *   - paired, interleaved (1 input file):  slot = id*2; REV calls delegate
 *   - paired, two input files:             slot_fwd = id*2, slot_rev = id*2+1
 *
 * Resume: if the current pass matches rstate->resume_pass, the worker initializes
 * its per-slot ThreadProgress from rstate->resume_progress and then skip-reads
 * records_consumed records from each owned slot before entering the alignment
 * loop. Skip-read uses the same next() that delivers records during alignment,
 * which transparently rebuilds the parser state (last_header, read_count, format
 * detection) inside vstate_in[slot] — no parser-state serialization required.
 *
 * Durability: the per-slot ThreadProgress entry is advanced atomically with the
 * read blob whenever a new hit lands (put_read_with_progress, single WriteBatch).
 * Otherwise it's flushed every COUNTER_FLUSH_INTERVAL records via put_progress.
 * A separate flush thread (started by align()) fsyncs the WAL periodically.
 */
void align2(int id, Readfeed& readfeed, Readstats& readstats,
            Index& index, References& refs, Refstats& refstats,
            KeyValueDatabase& kvdb, Runopts& opts,
            const restart::State* rstate)
{
	const uint16_t pass_i = static_cast<uint16_t>(index.index_num);
	const uint16_t pass_p = static_cast<uint16_t>(index.part);
	const uint32_t ns = readfeed.num_sense;
	const bool is_interleaved = (readfeed.num_orig_files < ns);
	const uint32_t slot_fwd = static_cast<uint32_t>(id) * ns;
	const uint32_t slot_rev = (ns > 1 && !is_interleaved) ? (slot_fwd + 1) : slot_fwd;
	const std::size_t n_indexes = opts.indexfiles.size();

	// Pick up per-slot progress from the resume state, if any. Either way the
	// deltas are sized to the number of indexes so the snapshot/compare logic
	// below can index them unconditionally.
	restart::ThreadProgress prog_fwd;
	restart::ThreadProgress prog_rev;
	prog_fwd.reads_matched_delta.assign(n_indexes, 0);
	prog_rev.reads_matched_delta.assign(n_indexes, 0);
	const bool is_resume_pass = rstate && rstate->is_resume
	    && rstate->resume_pass.first  == static_cast<int>(pass_i)
	    && rstate->resume_pass.second == static_cast<int>(pass_p);
	if (is_resume_pass) {
		const auto& rp = rstate->resume_progress;
		if (slot_fwd < rp.size() && rp[slot_fwd]) {
			prog_fwd = *rp[slot_fwd];
			if (prog_fwd.reads_matched_delta.size() != n_indexes)
				prog_fwd.reads_matched_delta.resize(n_indexes, 0);
		}
		if (ns > 1 && !is_interleaved && slot_rev < rp.size() && rp[slot_rev]) {
			prog_rev = *rp[slot_rev];
			if (prog_rev.reads_matched_delta.size() != n_indexes)
				prog_rev.reads_matched_delta.resize(n_indexes, 0);
		}
	}

	// Skip-read the slots that already produced records in the prior run.
	// next() advances vstate_in[slot] in lockstep with the reader, so the
	// parser state is rebuilt naturally — no special skip path needed.
	if (prog_fwd.records_consumed > 0 || prog_rev.records_consumed > 0) {
		auto t_skip_start = std::chrono::high_resolution_clock::now();
		std::string throwaway;
		uint32_t skipped_fwd = 0;
		for (uint32_t k = 0; k < prog_fwd.records_consumed; ++k) {
			if (!readfeed.next(static_cast<int>(slot_fwd), throwaway)) break;
			++skipped_fwd;
		}
		uint32_t skipped_rev = 0;
		if (ns > 1 && !is_interleaved) {
			for (uint32_t k = 0; k < prog_rev.records_consumed; ++k) {
				if (!readfeed.next(static_cast<int>(slot_rev), throwaway)) break;
				++skipped_rev;
			}
		}
		std::chrono::duration<double> dt = std::chrono::high_resolution_clock::now() - t_skip_start;
		INFO("Processor ", id, " resume: skip-read FWD=", skipped_fwd,
		     " REV=", skipped_rev, " records in ", dt.count(), " sec");
	}

	// Starting idx after skip — mirrors what the prior run's main loop would
	// have used for its next next() call.
	int idx;
	if (ns == 1) {
		idx = static_cast<int>(slot_fwd);
	} else if (is_interleaved) {
		// Single shared reader; idx alternates each iteration.
		idx = static_cast<int>(slot_fwd + (prog_fwd.records_consumed % ns));
	} else {
		// Two readers; FWD is always read first then toggled to REV. So if
		// counts are equal the next read is FWD; if FWD == REV+1 the next is
		// REV. (FWD < REV cannot happen.)
		idx = (prog_fwd.records_consumed > prog_rev.records_consumed)
		      ? static_cast<int>(slot_rev)
		      : static_cast<int>(slot_fwd);
	}

	uint32_t since_flush_fwd = 0;
	uint32_t since_flush_rev = 0;

	unsigned num_all = 0;
	unsigned num_skipped = 0;
	unsigned num_hit = 0;
	std::string readstr;

	auto t_start = std::chrono::high_resolution_clock::now();
	INFO("Processor ", id, " thread ", std::this_thread::get_id(),
	     " started (start_idx=", idx,
	     ", fwd_records=", prog_fwd.records_consumed,
	     ", rev_records=", prog_rev.records_consumed,
	     ", fwd_short=", prog_fwd.num_short_local,
	     ", rev_short=", prog_rev.num_short_local,
	     ", fwd_aligned=", prog_fwd.num_aligned_local,
	     ", rev_aligned=", prog_rev.num_aligned_local, ")");

	// Dual-file paired workers alternate FWD/REV. A resume checkpoint can have
	// fwd_records > rev_records (interrupt between the FWD read and its REV
	// mate), and chunk boundaries can leave one slot with an extra read. If we
	// break on the first EOF the lagging slot is stranded — drain it instead.
	bool drain_lagging = false;

	for (;;) {
		if (!readfeed.next(idx, readstr)) {
			if (ns > 1 && !is_interleaved) {
				const int other = (idx == static_cast<int>(slot_fwd))
				    ? static_cast<int>(slot_rev)
				    : static_cast<int>(slot_fwd);
				if (readfeed.next(other, readstr)) {
					idx = other;
					drain_lagging = true;
				} else {
					break;
				}
			} else {
				break;
			}
		}

		// Update the slot's progress for this iteration.
		const uint32_t sense = static_cast<uint32_t>(idx) % ns; // 0=FWD, 1=REV
		restart::ThreadProgress& prog = (is_interleaved || sense == 0) ? prog_fwd : prog_rev;
		uint32_t& since_flush          = (is_interleaved || sense == 0) ? since_flush_fwd : since_flush_rev;
		const uint32_t target_slot     = is_interleaved ? slot_fwd : (slot_fwd + sense);
		++prog.records_consumed;
		++since_flush;

		bool needs_blob_write = false;
		std::string blob_to_write;
		std::string read_id_for_write;

		{
			Read read(readstr);
			read.init(opts);
			read.is_too_short = read.sequence.size() < refstats.lnwin[index.index_num];

			if (read.is_too_short) {
				read.isValid = false;
				++prog.num_short_local;
				readstats.num_short.fetch_add(1, std::memory_order_relaxed);
			}

			if (read.isValid) {
				read.load_db(kvdb);
				// On resume, the read may carry alignv entries from a partial
				// prior attempt at the current pass. Trim them so traverse()
				// can rebuild them deterministically — re-alignment is then
				// idempotent regardless of whether the prior attempt finished.
				auto& v = read.alignment.alignv;
				bool had_same_pass = false;
				for (const auto& a : v) {
					if (a.index_num == pass_i && a.part == pass_p) { had_same_pass = true; break; }
				}
				if (had_same_pass) {
					v.erase(std::remove_if(v.begin(), v.end(),
						[&](const s_align2& a) {
							return a.index_num == pass_i && a.part == pass_p;
						}), v.end());
					if (v.empty()) { read.is_done = false; read.is_hit = false; }
					read.max_SW_count = 0;
				}
			}

			if (!read.isEmpty && read.isValid && !read.is_done) {
				// Snapshot the read state that alignment.cpp's update of
				// num_aligned / reads_matched_per_db is driven by. After the
				// strands loop we compare and replicate the same delta into
				// this slot's local counters — so on resume we can restore the
				// in-flight pass's contribution that the blob doesn't carry.
				const bool was_hit = read.is_hit;
				int before_min_idx_num = -1;
				if (!read.alignment.alignv.empty()) {
					const uint32_t mi = read.alignment.min_index;
					if (mi < read.alignment.alignv.size())
						before_min_idx_num = static_cast<int>(read.alignment.alignv[mi].index_num);
				}

				int num_strands = 0;
				const bool search_single_strand = opts.is_forward ^ opts.is_reverse;
				num_strands = search_single_strand ? 1 : 2;

				for (int count = 0; count < num_strands && !read.is_done; ++count) {
					if ((search_single_strand && opts.is_reverse) || count == 1) {
						if (!read.reversed) read.revIntStr();
					}
					traverse(opts, index, refs, readstats, refstats, read,
					         search_single_strand || count == 1);
					read.id_win_hits.clear();
				}

				// Mirror what alignment.cpp did to the global counters.
				if (!was_hit && read.is_hit) ++prog.num_aligned_local;
				int after_min_idx_num = -1;
				if (!read.alignment.alignv.empty()) {
					const uint32_t mi = read.alignment.min_index;
					if (mi < read.alignment.alignv.size())
						after_min_idx_num = static_cast<int>(read.alignment.alignv[mi].index_num);
				}
				if (before_min_idx_num != after_min_idx_num) {
					if (before_min_idx_num >= 0
					    && static_cast<std::size_t>(before_min_idx_num) < prog.reads_matched_delta.size())
						--prog.reads_matched_delta[before_min_idx_num];
					if (after_min_idx_num >= 0
					    && static_cast<std::size_t>(after_min_idx_num) < prog.reads_matched_delta.size())
						++prog.reads_matched_delta[after_min_idx_num];
				}

				if (read.is_hit) ++num_hit;
				if (read.is_new_hit) {
					needs_blob_write   = true;
					blob_to_write      = read.toBinString();
					read_id_for_write  = read.id;
				}
				++num_all;
			} else if (read.is_done) {
				++num_skipped;
			}
		} // ~Read destroyed

		// Atomically commit the alignment blob + the slot's progress, so that
		// "thread_done says K" is durable iff "all hits for records 0..K-1
		// produced by this worker are durable" (modulo OS pagecache for non-
		// fsync'd writes — see flush thread in align()).
		if (needs_blob_write) {
			restart::put_read_with_progress(kvdb,
				read_id_for_write, blob_to_write,
				pass_i, pass_p, target_slot, prog);
			since_flush = 0;
		} else if (since_flush >= COUNTER_FLUSH_INTERVAL) {
			restart::put_progress(kvdb, pass_i, pass_p, target_slot, prog);
			since_flush = 0;
		}

		readstr.resize(0);
		if (!drain_lagging && opts.is_paired) idx ^= 1;
	}

	// Final flush of the slot counters so the end-of-pass commit_pass sees
	// the canonical totals and so any later resume of a sibling pass starts
	// from a clean state.
	restart::put_progress(kvdb, pass_i, pass_p, slot_fwd, prog_fwd);
	if (ns > 1 && !is_interleaved) {
		restart::put_progress(kvdb, pass_i, pass_p, slot_rev, prog_rev);
	}

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t_start;
	INFO("Processor ", id, " thread ", std::this_thread::get_id(), " done. Processed ",
	     num_all, " reads. Skipped already processed: ", num_skipped, " reads",
	     " Aligned reads (passing E-value): ", num_hit,
	     " fwd_records=", prog_fwd.records_consumed,
	     " rev_records=", prog_rev.records_consumed,
	     " fwd_aligned=", prog_fwd.num_aligned_local,
	     " rev_aligned=", prog_rev.num_aligned_local,
	     " Runtime sec: ", elapsed.count());
} // ~align2

/*
 * Launches the alignment threads. Called from main.
 *
 * Concurrency: workers run independently and each owns one or two
 * gz_slots/flat_slots plus the corresponding thread_done/{i}/{p}/{slot} key.
 * There is no lock between workers and a writer; in particular the previous
 * shared_mutex + ticker design (which deadlocked under reader-preference
 * shared_mutex and never produced a watermark commit) is gone. The only
 * background thread here flushes the WAL periodically; it never touches
 * any shared sortmerna state.
 */
void align(Readfeed& readfeed, Readstats& readstats, Index& index, KeyValueDatabase& kvdb,
           Runopts& opts, const restart::State* rstate)
{
	INFO("==== Starting alignment ====");
	if (rstate && rstate->is_resume) {
		INFO("Restart: align() resuming. Completed passes: ", rstate->align_done.size(),
		     " latest=(", rstate->latest_pass.first, ",", rstate->latest_pass.second, ")",
		     " resume_pass=(", rstate->resume_pass.first, ",", rstate->resume_pass.second, ")");
	}
	restart::set_phase(kvdb, restart::Phase::Align);
	INFO("Alignment parameters:  is_best: ", opts.is_best,
		"  num_alignments: ", opts.num_alignments,
		"  min_lis: ", opts.min_lis);
	if (opts.num_alignments == 0) {
		INFO("num_alignments is set to: ", opts.num_alignments,
		    ", so all alignments passing E-value threshold will be reported,"
		    " and the option is_best is ignored.");
	}

	unsigned int numCores = std::thread::hardware_concurrency();
	INFO("Number of cores: ", numCores);

	const int numProcThread = opts.num_proc_thread;
	INFO("Using number of Processor threads: ", numProcThread);
	readfeed.init_reading();

	std::vector<std::thread> tpool;
	tpool.reserve(numProcThread);

	Refstats refstats(opts, readstats);
	References refs;

	int loopCount = 0;

	auto is_pass_done = [&](uint16_t i, uint16_t p) -> bool {
		if (!rstate || !rstate->is_resume) return false;
		for (const auto& pr : rstate->align_done) {
			if (pr.first == i && pr.second == p) return true;
		}
		return false;
	};

	auto start_a = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed;

	// Background WAL fsync thread. Workers issue per-read kvdb writes that hit
	// the WAL but are not fsync'd; this thread makes them durable against a
	// kernel crash. A process-only SIGKILL doesn't need this (OS pagecache
	// survives), but a power loss does — and the work is cheap.
	std::atomic<bool> stop_flush{false};
	std::mutex flush_mtx;
	std::condition_variable flush_cv;
	std::thread flush_thread([&]() {
		const auto delay = std::chrono::seconds(opts.flush_delay);
		while (true) {
			std::unique_lock<std::mutex> ul(flush_mtx);
			if (flush_cv.wait_for(ul, delay, [&]{ return stop_flush.load(); }))
				return;
			ul.unlock();
			kvdb.flush_wal();
			INFO("Restart: WAL fsync");
		}
	});

	for (size_t idx_num = 0; idx_num < opts.indexfiles.size(); ++idx_num)
	{
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[idx_num]; ++idx_part)
		{
			if (is_pass_done(static_cast<uint16_t>(idx_num), idx_part)) {
				INFO("Restart: skipping committed pass idx=", idx_num, " part=", idx_part + 1);
				continue;
			}

			INFO("Loading index: ", idx_num, " part: ", idx_part + 1, "/",
			     refstats.num_index_parts[idx_num], " Memory KB: ", (get_memory() >> 10), " ... ");
			auto start_i = std::chrono::high_resolution_clock::now();
			index.load(idx_num, idx_part, opts.indexfiles, refstats);

			// num_short resets at the start of every pass (per-pass semantics
			// in the original design). num_aligned and reads_matched_per_db
			// are cumulative across passes; they were restored from the last
			// commit_pass blob in the Readstats constructor on resume, but
			// that blob does not include the in-flight pass's contribution.
			// Re-apply the per-slot deltas saved by the prior run before
			// workers start, so the final blob reflects the full count.
			readstats.num_short.store(0, std::memory_order_relaxed);
			if (rstate && rstate->is_resume
			    && rstate->resume_pass.first  == static_cast<int>(idx_num)
			    && rstate->resume_pass.second == static_cast<int>(idx_part))
			{
				uint64_t resumed_short   = 0;
				uint64_t resumed_aligned = 0;
				std::vector<int64_t> resumed_matched(readstats.reads_matched_per_db.size(), 0);
				for (const auto& opt : rstate->resume_progress) {
					if (!opt) continue;
					resumed_short   += opt->num_short_local;
					resumed_aligned += opt->num_aligned_local;
					const auto m = std::min(resumed_matched.size(), opt->reads_matched_delta.size());
					for (std::size_t i = 0; i < m; ++i)
						resumed_matched[i] += opt->reads_matched_delta[i];
				}
				readstats.num_short.fetch_add(resumed_short,   std::memory_order_relaxed);
				readstats.num_aligned.fetch_add(resumed_aligned, std::memory_order_relaxed);
				for (std::size_t i = 0; i < resumed_matched.size(); ++i) {
					// reads_matched_per_db entries are plain uint64_t; apply
					// the int64 delta with wrap-safe arithmetic. Negative
					// deltas net out swaps recorded mid-pass.
					readstats.reads_matched_per_db[i] =
						static_cast<uint64_t>(static_cast<int64_t>(readstats.reads_matched_per_db[i])
						                      + resumed_matched[i]);
				}
				INFO("Restart: pass (", idx_num, ",", idx_part + 1,
				     ") resumed num_short=", resumed_short,
				     " resumed num_aligned=", resumed_aligned);
			}

			elapsed = std::chrono::high_resolution_clock::now() - start_i;
			INFO_MEM("done in [", elapsed.count(), "] sec");

			INFO("Loading references ...");
			start_i = std::chrono::high_resolution_clock::now();
			refs.load(idx_num, idx_part, opts, refstats);
			elapsed = std::chrono::high_resolution_clock::now() - start_i;
			INFO_MEM("done in [", elapsed.count(), "] sec.");

			start_i = std::chrono::high_resolution_clock::now();

			for (int i = 0; i < numProcThread; i++) {
				tpool.emplace_back(std::thread(align2, i, std::ref(readfeed),
				                    std::ref(readstats), std::ref(index), std::ref(refs),
				                    std::ref(refstats),  std::ref(kvdb), std::ref(opts),
				                    rstate));
			}
			for (auto& thr : tpool) thr.join();

			// All workers done for this pass. commit_pass writes align_done +
			// readstats blob in one batch and then clears thread_done/{i}/{p}/*.
			restart::commit_pass(kvdb,
				static_cast<uint16_t>(idx_num), idx_part,
				readstats.dbkey, readstats.toBstring());

			++loopCount;

			elapsed = std::chrono::high_resolution_clock::now() - start_i;
			INFO_MEM("done index: ", idx_num, " part: ", idx_part + 1, " in ", elapsed.count(), " sec");

			start_i = std::chrono::high_resolution_clock::now();
			index.unload();
			refs.unload();
			elapsed = std::chrono::high_resolution_clock::now() - start_i;
			INFO_MEM("Index and References unloaded in ", elapsed.count(), " sec.");
			tpool.clear();

			readfeed.rewind_in();
			readfeed.init_vzlib_in();  // noop for INDEXED feed; SPLIT_READS only
		} // ~for(idx_part)
	} // ~for(idx_num)

	// Stop the flush thread before returning. One last fsync afterwards to
	// catch anything written between its last cycle and now.
	{
		std::lock_guard<std::mutex> lk(flush_mtx);
		stop_flush.store(true);
	}
	flush_cv.notify_one();
	flush_thread.join();
	kvdb.flush_wal();

	elapsed = std::chrono::high_resolution_clock::now() - start_a;
	INFO("==== Done alignment in ", elapsed.count(), " sec ====\n");

	readstats.set_is_set_aligned_id_cov();
	readstats.store_to_db(kvdb);
} // ~align

void denovo_stats_run(const uint32_t& id,
	Readfeed& readfeed,
	Readstats& readstats,
	References& refs,
	KeyValueDatabase& kvdb,
	Runopts& opts)
{
	uint64_t countReads = 0;
	uint64_t num_invalid = 0; // empty or invalid reads count
	uint16_t num_reads = opts.is_paired ? 2 : 1;
	std::string readstr;
	std::vector<Read> reads; // two reads if paired, a single read otherwise

	if (opts.dbg_level == 2)
		INFO_MEM("Denovo stats thread ", id, " : ", std::this_thread::get_id(), " started.");

	for (bool isDone = false; !isDone;)
	{
		reads.clear();
		uint32_t idx = id * readfeed.num_sense; // index into split_files array
		for (uint16_t i = 0; i < num_reads; ++i)
		{
			if (readfeed.next(idx, readstr))
			{
				reads.emplace_back(Read(readstr));
				reads[i].init(opts);
				reads[i].load_db(kvdb);
				readstr.resize(0);
				++countReads;
			}
			else {
				isDone = true;
			}
			if (opts.is_paired) idx ^= 1; // switch fwd-rev
		}

		if (!isDone) {
			if (reads.back().isEmpty || !reads.back().isValid) {
				++num_invalid;
				continue;
			}

			for (auto &read: reads) {
				if (read.is03) read.flip34();
				for (auto const& align : read.alignment.alignv) {
					if (align.index_num == refs.num	&& align.part == refs.part)	{
						auto miss_gap_match = read.calc_miss_gap_match(refs, align);
						auto idr = std::floor(std::get<3>(miss_gap_match) * 1000.0 + 0.5) / 1000.0; // round to 3 decimal
						auto covr = std::floor(std::get<4>(miss_gap_match) * 1000.0 + 0.5) / 1000.0;
						auto is_id = idr >= opts.min_id;
						auto is_cov = covr >= opts.min_cov;
						//auto is_id = std::get<3>(miss_gap_match) >= opts.min_id;
						//auto is_cov = std::get<4>(miss_gap_match)>= opts.min_cov;
						if (is_id && is_cov) {
							++read.c_yid_ycov;
							readstats.n_yid_ycov.fetch_add(1, std::memory_order_relaxed);
						}
						else if (is_id) {
							++read.n_yid_ncov;
							readstats.n_yid_ncov.fetch_add(1, std::memory_order_relaxed);
						}
						else if (is_cov) {
							++read.n_nid_ycov;
							readstats.n_nid_ycov.fetch_add(1, std::memory_order_relaxed);
						}
						else {
							++read.n_denovo;
							readstats.num_denovo.fetch_add(1, std::memory_order_relaxed); // neither ID nor COV
						}
					}
				}
				kvdb.put(read.id, read.toBinString()); // store to DB
			} // ~for reads
		} // ~ if !is_done
	} // ~for

	//std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start; // ~20 sec Debug/Win
	INFO_MEM("Denovo stats thread ", id, " : ", std::this_thread::get_id(), " done. Processed reads: ", countReads,
		" Invalid reads: ", num_invalid); // , " denovo count: ", denovo_n
} // ~denovo_stats_run

void denovo_stats(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts)
{
	INFO("==== processing Denovo statistics ====");
	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed;

	int nthreads = 0;
	//if (readfeed.type == FEED_TYPE::SPLIT_READS || readfeed.type == FEED_TYPE::INDEXED_GZ || readfeed.type == FEED_TYPE::INDEXED_FLAT) {
	nthreads = opts.num_proc_thread;
	readfeed.init_reading(); // prepare readfeed
	//}

	std::vector<std::thread> tpool;
	tpool.reserve(nthreads);

	bool indb = readstats.restoreFromDb(kvdb);
	if (indb) {
		INFO("Restored Readstats from DB: ", indb);
	}

	Refstats refstats(opts, readstats);
	References refs;

	// loop through every reference file passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t ref_idx = 0; ref_idx < opts.indexfiles.size(); ++ref_idx)
	{
		// iterate all parts of the index
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[ref_idx]; ++idx_part)
		{
			INFO_NE("loading reference ", ref_idx, " part ", idx_part + 1, "/", refstats.num_index_parts[ref_idx]);
			auto start_i = std::chrono::high_resolution_clock::now();
			refs.load(ref_idx, idx_part, opts, refstats);
			elapsed = std::chrono::high_resolution_clock::now() - start_i;
			INFO_NS(" ... done in sec ", elapsed.count(), "\n");

			start_i = std::chrono::high_resolution_clock::now(); // index processing starts

			// start threads
			//if (opts.feed_type == FEED_TYPE::SPLIT_READS || opts.feed_type == FEED_TYPE::INDEXED_GZ || opts.feed_type == FEED_TYPE::INDEXED_FLAT) {
			for (int i = 0; i < nthreads; ++i) {
				tpool.emplace_back(std::thread(denovo_stats_run, i, std::ref(readfeed),
					std::ref(readstats), std::ref(refs), std::ref(kvdb), std::ref(opts)));
			}
			//}
			// wait for all threads to finish
			for (auto& thr: tpool) {
				thr.join();
			}

			elapsed = std::chrono::high_resolution_clock::now() - start_i; // index processing done
			INFO("done reference ", ref_idx, " part: ", idx_part + 1, " in ", elapsed.count(), " sec");

			start_i = std::chrono::high_resolution_clock::now();
			refs.unload();
			//read_queue.reset();
			elapsed = std::chrono::high_resolution_clock::now() - start_i;
			INFO_MEM("references unloaded in ", elapsed.count(), " sec");
			tpool.clear();
			// rewind for the next index
			readfeed.rewind_in();
			readfeed.init_vzlib_in();
		} // ~for(idx_part)
	} // ~for(ref_idx)

	elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("num_yid_ycov: ", readstats.n_yid_ycov,
		"\n\t\t   num_yid_ncov: ", readstats.n_yid_ncov,
		"\n\t\t   num_nid_ycov: ", readstats.n_nid_ycov,
		"\n\t\t   num_denovo: ", readstats.num_denovo);
	INFO("=== done Denovo stats in ", elapsed.count(), " sec ===\n");
} // ~denovo_stats
