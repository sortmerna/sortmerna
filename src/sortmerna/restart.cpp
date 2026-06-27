/*
@copyright 2016-2026 Clarity Genomics BVBA

SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
Licensed under the GNU LGPL v3 or later. See COPYING.LESSER for details.
*/

#include "restart.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>

#include "common.hpp"
#include "kvdb.hpp"
#include "options.hpp"
#include "version.h"

namespace restart {

namespace {

// FNV-1a 64-bit. We don't need cryptographic strength here: the fingerprint
// guards against accidental mismatches (different opts, replaced inputs).
// Birthday-collision probability over realistic option-set sizes is negligible.
constexpr uint64_t FNV64_OFFSET = 0xcbf29ce484222325ULL;
constexpr uint64_t FNV64_PRIME  = 0x100000001b3ULL;

uint64_t fnv1a(const std::string& s, uint64_t h = FNV64_OFFSET) {
	for (unsigned char c : s) {
		h ^= c;
		h *= FNV64_PRIME;
	}
	return h;
}

std::string to_hex(uint64_t v) {
	char buf[17];
	std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(v));
	return std::string(buf);
}

std::string writer_version() {
	std::ostringstream ss;
	ss << SORTMERNA_MAJOR << '.' << SORTMERNA_MINOR << '.' << SORTMERNA_PATCH;
	return ss.str();
}

std::string make_run_id() {
	// Time-seeded random 128-bit id. Not for security; uniqueness across one
	// user's runs is sufficient.
	auto now = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937_64 rng(static_cast<uint64_t>(now));
	uint64_t a = rng();
	uint64_t b = rng();
	return to_hex(a) + to_hex(b);
}

// File size + mtime fingerprint. Missing file returns "MISSING".
std::string stat_signature(const std::string& path) {
	std::error_code ec;
	auto canon = std::filesystem::weakly_canonical(path, ec);
	std::string p = ec ? path : canon.string();
	struct stat st {};
	if (::stat(p.c_str(), &st) != 0) {
		return p + "|MISSING";
	}
	std::ostringstream ss;
	ss << p << '|' << static_cast<long long>(st.st_size)
	   << '|' << static_cast<long long>(st.st_mtime);
	return ss.str();
}

std::string pair_key(const char* prefix, uint16_t i, uint16_t p) {
	std::ostringstream ss;
	ss << prefix << i << '/' << p;
	return ss.str();
}

// ThreadProgress layout (little-endian, variable length):
//   [ 0.. 3] uint32 records_consumed
//   [ 4..11] uint64 num_short_local
//   [12..19] uint64 num_aligned_local
//   [20..23] uint32 reads_matched_delta_size (N)
//   [24..  ] N * int64 reads_matched_delta entries
constexpr std::size_t TP_FIXED_BYTES = 24;

void push_u32_le(std::string& buf, uint32_t v) {
	for (int i = 0; i < 4; ++i)
		buf.push_back(static_cast<char>((v >> (i * 8)) & 0xff));
}

void push_u64_le(std::string& buf, uint64_t v) {
	for (int i = 0; i < 8; ++i)
		buf.push_back(static_cast<char>((v >> (i * 8)) & 0xff));
}

void push_i64_le(std::string& buf, int64_t v) {
	push_u64_le(buf, static_cast<uint64_t>(v));
}

uint32_t read_u32_le(const char* p) {
	uint32_t v = 0;
	for (int i = 0; i < 4; ++i)
		v |= static_cast<uint32_t>(static_cast<uint8_t>(p[i])) << (i * 8);
	return v;
}

uint64_t read_u64_le(const char* p) {
	uint64_t v = 0;
	for (int i = 0; i < 8; ++i)
		v |= static_cast<uint64_t>(static_cast<uint8_t>(p[i])) << (i * 8);
	return v;
}

int64_t read_i64_le(const char* p) {
	return static_cast<int64_t>(read_u64_le(p));
}

} // anonymous namespace

std::string encode_thread_progress(const ThreadProgress& tp) {
	std::string buf;
	const uint32_t n = static_cast<uint32_t>(tp.reads_matched_delta.size());
	buf.reserve(TP_FIXED_BYTES + 8 * n);
	push_u32_le(buf, tp.records_consumed);
	push_u64_le(buf, tp.num_short_local);
	push_u64_le(buf, tp.num_aligned_local);
	push_u32_le(buf, n);
	for (auto v : tp.reads_matched_delta) push_i64_le(buf, v);
	return buf;
}

std::optional<ThreadProgress> decode_thread_progress(const std::string& s) {
	if (s.size() < TP_FIXED_BYTES) return std::nullopt;
	ThreadProgress tp;
	tp.records_consumed   = read_u32_le(s.data() + 0);
	tp.num_short_local    = read_u64_le(s.data() + 4);
	tp.num_aligned_local  = read_u64_le(s.data() + 12);
	const uint32_t n      = read_u32_le(s.data() + 20);
	if (s.size() != TP_FIXED_BYTES + 8 * static_cast<std::size_t>(n)) return std::nullopt;
	tp.reads_matched_delta.resize(n);
	for (uint32_t i = 0; i < n; ++i) {
		tp.reads_matched_delta[i] = read_i64_le(s.data() + TP_FIXED_BYTES + 8 * i);
	}
	return tp;
}

std::string phase_to_string(Phase p) {
	switch (p) {
		case Phase::Init:      return "init";
		case Phase::Align:     return "align";
		case Phase::PostAlign: return "post_align";
		case Phase::Denovo:    return "denovo";
		case Phase::Report:    return "report";
		case Phase::Done:      return "done";
	}
	return "init";
}

Phase phase_from_string(const std::string& s) {
	if (s == "align")      return Phase::Align;
	if (s == "post_align") return Phase::PostAlign;
	if (s == "denovo")     return Phase::Denovo;
	if (s == "report")     return Phase::Report;
	if (s == "done")       return Phase::Done;
	return Phase::Init;
}

std::string key_align_done(uint16_t i, uint16_t p) {
	return pair_key(PFX_ALIGN_DONE, i, p);
}
std::string key_thread_done_prefix(uint16_t i, uint16_t p) {
	std::ostringstream ss;
	ss << PFX_THREAD_DONE << i << '/' << p << '/';
	return ss.str();
}
std::string key_thread_done(uint16_t i, uint16_t p, uint32_t slot) {
	std::ostringstream ss;
	ss << PFX_THREAD_DONE << i << '/' << p << '/' << slot;
	return ss.str();
}
std::string key_report_done(const std::string& kind) {
	return std::string(PFX_REPORT_DONE) + kind;
}

// Fingerprint of resume-fingerprinted Runopts fields. Path opts, thread counts,
// logging flags, indexing-only options, and the task selector are intentionally
// omitted: changing any of those is safe across a resume.
std::string opts_fingerprint(const Runopts& opts) {
	std::ostringstream ss;
	// Alignment control.
	ss << "is_best="            << opts.is_best            << ';'
	   << "is_best_id_cov="     << opts.is_best_id_cov     << ';'
	   << "is_min_lis="         << opts.is_min_lis         << ';'
	   << "is_num_alignments="  << opts.is_num_alignments  << ';'
	   << "is_full_search="     << opts.is_full_search     << ';'
	   << "is_forward="         << opts.is_forward         << ';'
	   << "is_reverse="         << opts.is_reverse         << ';'
	   << "is_score_split="     << opts.is_score_split     << ';';
	// Pairing / output structure.
	ss << "is_paired="          << opts.is_paired          << ';'
	   << "is_paired_in="       << opts.is_paired_in       << ';'
	   << "is_paired_out="      << opts.is_paired_out      << ';'
	   << "is_out2="            << opts.is_out2            << ';'
	   << "is_sout="            << opts.is_sout            << ';';
	// Report toggles (strict: changing these between runs invalidates resume).
	ss << "is_otu_map="         << opts.is_otu_map         << ';'
	   << "is_denovo="          << opts.is_denovo          << ';'
	   << "is_print_all_reads=" << opts.is_print_all_reads << ';'
	   << "is_sam="             << opts.is_sam             << ';'
	   << "is_SQ="              << opts.is_SQ              << ';'
	   << "is_blast="           << opts.is_blast           << ';'
	   << "is_fastx="           << opts.is_fastx           << ';'
	   << "is_other="           << opts.is_other           << ';'
	   << "zip_out="            << opts.zip_out            << ';';
	// Numeric thresholds and SW parameters.
	ss << "num_alignments=" << opts.num_alignments << ';'
	   << "num_seeds="      << opts.num_seeds      << ';'
	   << "min_lis="        << opts.min_lis        << ';'
	   << "edges="          << opts.edges          << ';'
	   << "match="          << opts.match          << ';'
	   << "mismatch="       << opts.mismatch       << ';'
	   << "gap_open="       << opts.gap_open       << ';'
	   << "gap_ext="        << opts.gap_extension  << ';'
	   << "score_N="        << opts.score_N        << ';'
	   << "feed_type="      << static_cast<int>(opts.feed_type) << ';'
	   << "evalue="         << opts.evalue         << ';'
	   << "min_id="         << opts.min_id         << ';'
	   << "min_cov="        << opts.min_cov        << ';'
	   << "max_read_len="   << opts.max_read_len   << ';'
	   << "seed_win_len="   << opts.seed_win_len   << ';'
	   // Per-slot resume counters are tied to the [thread_id * num_sense + sense]
	   // slot layout, so num_proc_thread MUST match across the original and
	   // resumed runs. Folding it into the fingerprint makes a mismatch abort
	   // with a diagnostic instead of silently mis-aligning records.
	   << "num_proc_thread=" << opts.num_proc_thread << ';';
	// blastops vector (each entry can affect tabular output schema).
	ss << "blastops=[";
	for (const auto& b : opts.blastops) ss << b << ',';
	ss << "];";
	// skiplengths derived from --passes; included so changes to passes invalidate.
	ss << "skiplengths=[";
	for (const auto& v : opts.skiplengths) {
		ss << '(';
		for (auto x : v) ss << x << ',';
		ss << "),";
	}
	ss << "];";

	return to_hex(fnv1a(ss.str()));
}

// File-presence fingerprint of --reads and --ref. We hash a canonicalized
// "path|size|mtime" line per file so any input swap aborts resume.
std::string inputs_fingerprint(const Runopts& opts) {
	std::ostringstream ss;
	ss << "reads:\n";
	for (const auto& r : opts.readfiles) {
		ss << "  " << stat_signature(r) << '\n';
	}
	ss << "refs:\n";
	for (const auto& ref : opts.indexfiles) {
		ss << "  " << stat_signature(ref.first) << '\n';
	}
	return ss.str();
}

State probe_or_init(KeyValueDatabase& kvdb, const Runopts& opts) {
	State st;
	const std::string current_opts_hash = opts_fingerprint(opts);
	const std::string current_inputs    = inputs_fingerprint(opts);

	if (!kvdb.has(KEY_VERSION)) {
		// Fresh kvdb. Stamp identity + initial phase.
		st.is_resume = false;
		st.phase = Phase::Init;
		st.run_id = make_run_id();
		kvdb.put(KEY_VERSION,    writer_version());
		kvdb.put(KEY_RUN_ID,     st.run_id);
		kvdb.put(KEY_OPTS_HASH,  current_opts_hash);
		kvdb.put(KEY_INPUTS,     current_inputs);
		kvdb.put(KEY_PHASE,      phase_to_string(Phase::Init));
		INFO("Restart: fresh kvdb, stamped run_id=", st.run_id, " writer=", writer_version());
		return st;
	}

	// Existing kvdb. Verify fingerprints before allowing resume.
	const std::string stored_version    = kvdb.get(KEY_VERSION);
	const std::string stored_opts_hash  = kvdb.get(KEY_OPTS_HASH);
	const std::string stored_inputs     = kvdb.get(KEY_INPUTS);
	st.run_id                           = kvdb.get(KEY_RUN_ID);

	if (stored_opts_hash != current_opts_hash) {
		ERR("Cannot resume: stored opts_hash ", stored_opts_hash,
		    " differs from current ", current_opts_hash,
		    ". Either re-run with the original options or wipe the kvdb at ",
		    opts.kvdbdir.string());
		throw std::runtime_error("kvdb resume aborted: opts mismatch");
	}
	if (stored_inputs != current_inputs) {
		ERR("Cannot resume: input file fingerprints changed.\n  stored:\n",
		    stored_inputs, "  current:\n", current_inputs,
		    "Either restore the original inputs or wipe the kvdb at ",
		    opts.kvdbdir.string());
		throw std::runtime_error("kvdb resume aborted: inputs mismatch");
	}

	st.is_resume = true;
	st.phase = phase_from_string(kvdb.get(KEY_PHASE));

	// Harvest completed alignment passes.
	kvdb.iter_prefix(PFX_ALIGN_DONE, [&](const std::string& k, const std::string&) {
		// k = "__sm/v1/align_done/{i}/{p}"
		auto tail = k.substr(std::string(PFX_ALIGN_DONE).size());
		auto slash = tail.find('/');
		if (slash == std::string::npos) return true;
		uint16_t i = static_cast<uint16_t>(std::stoul(tail.substr(0, slash)));
		uint16_t p = static_cast<uint16_t>(std::stoul(tail.substr(slash + 1)));
		st.align_done.emplace_back(i, p);
		std::pair<int, int> here{static_cast<int>(i), static_cast<int>(p)};
		if (here > st.latest_pass) st.latest_pass = here;
		return true;
	});

	auto pass_is_done = [&](uint16_t i, uint16_t p) {
		for (const auto& pr : st.align_done) {
			if (pr.first == i && pr.second == p) return true;
		}
		return false;
	};

	// Harvest the in-flight pass via thread_done/*. Any (i, p) that has at
	// least one thread_done entry and no matching align_done is the pass
	// alignment was inside when the prior run stopped. There can only be one
	// such pass: commit_pass clears thread_done/{i}/{p}/* atomically with
	// stamping align_done.
	std::pair<int, int> in_flight = {-1, -1};
	std::vector<std::pair<uint32_t, ThreadProgress>> tp_pairs;
	kvdb.iter_prefix(PFX_THREAD_DONE, [&](const std::string& k, const std::string& v) {
		// k = "__sm/v1/thread_done/{i}/{p}/{slot}"
		auto tail = k.substr(std::string(PFX_THREAD_DONE).size());
		auto s1 = tail.find('/');
		if (s1 == std::string::npos) return true;
		auto s2 = tail.find('/', s1 + 1);
		if (s2 == std::string::npos) return true;
		uint16_t i = static_cast<uint16_t>(std::stoul(tail.substr(0, s1)));
		uint16_t p = static_cast<uint16_t>(std::stoul(tail.substr(s1 + 1, s2 - s1 - 1)));
		uint32_t slot = static_cast<uint32_t>(std::stoul(tail.substr(s2 + 1)));
		if (pass_is_done(i, p)) return true; // stale: pass was committed afterwards
		if (in_flight.first == -1) in_flight = {i, p};
		if (in_flight.first == i && in_flight.second == p) {
			auto decoded = decode_thread_progress(v);
			if (!decoded) {
				WARN("Restart: ignoring malformed thread_done/", i, "/", p, "/", slot,
				     " (", v.size(), " bytes)");
				return true;
			}
			tp_pairs.emplace_back(slot, *decoded);
		}
		return true;
	});
	st.resume_pass = in_flight;

	if (!tp_pairs.empty()) {
		uint32_t max_slot = 0;
		for (const auto& pr : tp_pairs) if (pr.first > max_slot) max_slot = pr.first;
		st.resume_progress.assign(max_slot + 1, std::nullopt);
		for (auto& pr : tp_pairs) st.resume_progress[pr.first] = pr.second;
	}

	uint64_t total_records = 0;
	uint64_t total_short   = 0;
	uint64_t total_aligned = 0;
	for (const auto& opt : st.resume_progress) {
		if (!opt) continue;
		total_records += opt->records_consumed;
		total_short   += opt->num_short_local;
		total_aligned += opt->num_aligned_local;
	}

	INFO("Restart: resuming kvdb run_id=", st.run_id, " writer(stored)=", stored_version,
	     " phase=", phase_to_string(st.phase),
	     " passes_done=", st.align_done.size(),
	     " resume_pass=(", st.resume_pass.first, ",", st.resume_pass.second, ")",
	     " resume_slots=", st.resume_progress.size(),
	     " resume_records=", total_records,
	     " resume_short=", total_short,
	     " resume_aligned=", total_aligned);
	return st;
}

void mark_align_done(KeyValueDatabase& kvdb, uint16_t i, uint16_t p) {
	kvdb.put(key_align_done(i, p), "");
}

void commit_pass(KeyValueDatabase& kvdb, uint16_t i, uint16_t p,
                 const std::string& readstats_dbkey,
                 const std::string& readstats_blob)
{
	// align_done + Readstats blob in one batch — the consistency boundary.
	kvdb.put_batch({
		{key_align_done(i, p), std::string()},
		{readstats_dbkey,       readstats_blob}
	});
	// Then drop the in-flight progress for this pass. Leaving them around is
	// harmless (they would be ignored thanks to align_done), but cleaning up
	// keeps the kvdb compact and probe_or_init faster on the next resume.
	clear_thread_done(kvdb, i, p);
}

void put_read_with_progress(KeyValueDatabase& kvdb,
                            const std::string& read_id,
                            const std::string& read_blob,
                            uint16_t i, uint16_t p, uint32_t slot,
                            const ThreadProgress& tp)
{
	kvdb.put_batch({
		{read_id,                       read_blob},
		{key_thread_done(i, p, slot),   encode_thread_progress(tp)}
	});
}

void put_progress(KeyValueDatabase& kvdb,
                  uint16_t i, uint16_t p, uint32_t slot,
                  const ThreadProgress& tp)
{
	kvdb.put(key_thread_done(i, p, slot), encode_thread_progress(tp));
}

void clear_thread_done(KeyValueDatabase& kvdb, uint16_t i, uint16_t p) {
	kvdb.delete_prefix(key_thread_done_prefix(i, p));
}

void set_phase(KeyValueDatabase& kvdb, Phase p) {
	kvdb.put(KEY_PHASE, phase_to_string(p));
}

void mark_report_done(KeyValueDatabase& kvdb, const std::string& kind) {
	kvdb.put(key_report_done(kind), "");
}

bool is_report_done(KeyValueDatabase& kvdb, const std::string& kind) {
	return kvdb.has(key_report_done(kind));
}

void mark_denovo_done(KeyValueDatabase& kvdb) {
	kvdb.put(KEY_DENOVO_DONE, "");
}

bool is_denovo_done(KeyValueDatabase& kvdb) {
	return kvdb.has(KEY_DENOVO_DONE);
}

} // namespace restart
