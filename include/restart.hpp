/*
@copyright 2016-2026 Clarity Genomics BVBA

SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
Licensed under the GNU LGPL v3 or later. See COPYING.LESSER for details.
*/

/*
 * file: restart.hpp
 *
 * Persists pipeline progress in the kvdb so an interrupted run resumes
 * automatically when sortmerna is re-invoked against the same kvdb.
 *
 * Schema (all keys prefixed with "__sm/v1/"):
 *   version              writer version (string)
 *   run_id               UUID of the run that created this kvdb (string)
 *   opts_hash            FNV-1a 64 hex of resume-fingerprinted Runopts fields
 *   inputs               canonical "{path,size,mtime}" lines for --reads/--ref
 *   phase                init | align | post_align | denovo | report | done
 *   align_done/{i}/{p}   marker: alignment pass (idx_num=i, idx_part=p) committed.
 *                        Written in the same atomic batch as the Readstats blob
 *                        store, so on restart "align_done present" implies
 *                        "Readstats reflects pass (i,p) and all prior passes".
 *   pending_reads/{i}/{p}/{read.id}  empty marker: read was touched during the
 *                                    in-flight pass (i,p) and may need rollback
 *                                    if pass never reaches align_done
 *   denovo_done          marker: denovo_stats phase finished
 *   report_done/{kind}   marker: report kind (fastx|sam|blast|biom|denovo|otu_map) finished
 *
 * Resume policy:
 *   - Probe kvdb on open. If version absent -> fresh kvdb; stamp run identity.
 *   - If present and opts_hash + inputs match current -> resume.
 *   - If present and either fingerprint differs -> abort with a diagnostic.
 *
 * Allow-list: fields varied between runs without invalidating resume are:
 *   path opts (workdir/kvdbdir/idxdir/outdir/readb_dir/aligned_pfx/other_pfx),
 *   thread counts, verbose/dbg_level/is_log/is_pid/is_cmd/is_dbg_put_kvdb,
 *   cmdline, task, queue_size_max, exit_early, is_index_built,
 *   indexing-only fields (max_file_size, interval, max_pos, findex, tmpdir).
 * Everything else affecting alignment / classification / output structure is
 * folded into opts_hash. See restart.cpp::opts_fingerprint.
 */

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

class KeyValueDatabase;
struct Runopts;

namespace restart {

constexpr const char* SCHEMA_VERSION = "1";

// Sentinel/control keys.
constexpr const char* KEY_VERSION     = "__sm/v1/version";
constexpr const char* KEY_RUN_ID      = "__sm/v1/run_id";
constexpr const char* KEY_OPTS_HASH   = "__sm/v1/opts_hash";
constexpr const char* KEY_INPUTS      = "__sm/v1/inputs";
constexpr const char* KEY_PHASE       = "__sm/v1/phase";
constexpr const char* KEY_DENOVO_DONE = "__sm/v1/denovo_done";

// Prefix-namespaced keys (key builders below).
constexpr const char* PFX_ALIGN_DONE    = "__sm/v1/align_done/";
constexpr const char* PFX_PENDING_READS = "__sm/v1/pending_reads/";
constexpr const char* PFX_REPORT_DONE   = "__sm/v1/report_done/";

enum class Phase {
	Init,
	Align,
	PostAlign,
	Denovo,
	Report,
	Done
};

std::string phase_to_string(Phase p);
Phase phase_from_string(const std::string& s);

// Key builders.
std::string key_align_done(uint16_t i, uint16_t p);
std::string key_pending_reads_prefix(uint16_t i, uint16_t p);
std::string key_pending_read(uint16_t i, uint16_t p, const std::string& read_id);
std::string key_report_done(const std::string& kind);

// Snapshot of what probe_or_init found in the kvdb. Drives the resume logic
// in align() and the higher-level main() phase dispatch.
struct State {
	// false: fresh kvdb (no version sentinel); true: existing kvdb whose
	// opts_hash + inputs matched the current run, so we may resume.
	bool is_resume = false;
	Phase phase = Phase::Init;
	// (idx_num, idx_part) of passes already committed (align_done present).
	std::vector<std::pair<uint16_t, uint16_t>> align_done;
	// Latest (i, p) committed via commit_pass; -1/-1 if none. Used as the
	// rewind target for lastIndex/lastPart on rolled-back reads.
	std::pair<int, int> latest_pass = {-1, -1};
	// At most one partial pass to roll back on resume: a pass with
	// pending_reads present but no matching align_done.
	std::pair<int, int> rollback_pass = {-1, -1};
	// Read ids touched during rollback_pass, harvested from pending_reads/*.
	std::vector<std::string> rollback_read_ids;
	std::string run_id;
};

// On open: stamp run identity if fresh, else load existing progress.
// Throws std::runtime_error on opts_hash / inputs mismatch.
State probe_or_init(KeyValueDatabase& kvdb, const Runopts& opts);

// Progress writes.
//
// commit_pass atomically writes the align_done marker and the Readstats blob
// (stored under readstats_dbkey, which is Readstats::dbkey) in a single
// WriteBatch. This is the consistency boundary: "align_done present" implies
// "Readstats blob reflects pass (i,p) and all prior passes".
void commit_pass(KeyValueDatabase& kvdb, uint16_t i, uint16_t p,
                 const std::string& readstats_dbkey,
                 const std::string& readstats_blob);
void mark_align_done(KeyValueDatabase& kvdb, uint16_t i, uint16_t p);
void add_pending_read(KeyValueDatabase& kvdb, uint16_t i, uint16_t p, const std::string& read_id);
// Atomically write a read kvdb entry + its pending-read marker for pass (i,p).
// Used in place of plain kvdb.put(read.id, blob) inside alignment threads so
// the pending marker can never lag the read mutation.
void record_read_with_pending(KeyValueDatabase& kvdb, uint16_t i, uint16_t p,
                              const std::string& read_id,
                              const std::string& read_blob);
void clear_pending_reads(KeyValueDatabase& kvdb, uint16_t i, uint16_t p);
void set_phase(KeyValueDatabase& kvdb, Phase p);
void mark_report_done(KeyValueDatabase& kvdb, const std::string& kind);
bool is_report_done(KeyValueDatabase& kvdb, const std::string& kind);
void mark_denovo_done(KeyValueDatabase& kvdb);
bool is_denovo_done(KeyValueDatabase& kvdb);

// Fingerprint helpers (exposed for testing).
std::string opts_fingerprint(const Runopts& opts);
std::string inputs_fingerprint(const Runopts& opts);

// Common report kinds used for report_done/{kind} markers.
namespace report_kind {
	constexpr const char* FASTX_ALIGNED = "fastx_aligned";
	constexpr const char* FASTX_OTHER   = "fastx_other";
	constexpr const char* SAM           = "sam";
	constexpr const char* BLAST         = "blast";
	constexpr const char* BIOM          = "biom";
	constexpr const char* DENOVO        = "denovo";
	constexpr const char* OTU_MAP       = "otu_map";
}

} // namespace restart
