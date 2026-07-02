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
 *   thread_done/{i}/{p}/{slot}   binary {uint32 records_consumed,
 *                                uint64 num_short_local}. Each alignment worker
 *                                writes this for the slots it owns. Atomically
 *                                batched with the read+alignment blob when a
 *                                new hit lands, otherwise written every N reads.
 *                                On resume each slot is reopened at bytes_start
 *                                and the worker skip-reads records_consumed
 *                                records before starting alignment, which
 *                                naturally rebuilds the parser state (last_header,
 *                                read_count, format detection) without having to
 *                                serialize it.
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
 *   verbose/dbg_level/is_log/is_pid/is_cmd/is_dbg_put_kvdb, cmdline, task,
 *   queue_size_max, exit_early, is_index_built, indexing-only fields
 *   (max_file_size, interval, max_pos, findex, tmpdir).
 * Note: num_proc_thread IS in the fingerprint — the per-slot resume schema
 * (thread_done/{i}/{p}/{slot}) is layout-dependent and silently breaks if the
 * thread count differs across runs. Everything else affecting alignment /
 * classification / output structure is folded into opts_hash. See
 * restart.cpp::opts_fingerprint.
 *
 * Concurrency / durability model (per-thread, lock-free between workers and
 * the durability path):
 *   - Workers update their own thread_done entry.
 *   - A background helper thread periodically calls kvdb.flush_wal(sync=true).
 *     This makes everything written into the WAL up to that moment durable
 *     against a kernel crash. A process-only crash (SIGKILL) preserves the
 *     OS page cache and is already covered.
 *   - Re-alignment is idempotent. When a worker resumes at counter K and the
 *     kvdb already holds partial same-pass alignv entries for read K+1 from
 *     the prior run, align2 trims same-pass entries right after load_db; the
 *     subsequent traverse rebuilds them.
 */

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

class KeyValueDatabase;
struct Runopts;

namespace restart {

// Per-slot durable progress. Captures everything needed to resume the
// in-flight pass without loss:
//   - records_consumed: where the worker restarts (skip-read count).
//   - num_short_local:  this slot's contribution to Readstats::num_short for
//                       the pass (num_short is reset at the start of each pass).
//   - num_aligned_local: this slot's contribution to Readstats::num_aligned
//                        for the pass. num_aligned is cumulative across passes;
//                        without per-slot tracking, resume restores it from the
//                        last commit_pass blob (which doesn't include the
//                        in-flight pass) and we lose all pre-resume hits.
//   - reads_matched_delta: per-index net delta this slot contributed to
//                          Readstats::reads_matched_per_db during this pass.
//                          Length matches opts.indexfiles.size().
// On resume, sums across slots are added to the global Readstats counters AFTER
// the blob is restored and BEFORE workers start.
struct ThreadProgress {
	uint32_t             records_consumed   = 0;
	uint64_t             num_short_local    = 0;
	uint64_t             num_aligned_local  = 0;
	std::vector<int64_t> reads_matched_delta; // size = opts.indexfiles.size()
};

constexpr const char* SCHEMA_VERSION = "1";

// Sentinel/control keys.
constexpr const char* KEY_VERSION     = "__sm/v1/version";
constexpr const char* KEY_RUN_ID      = "__sm/v1/run_id";
constexpr const char* KEY_OPTS_HASH   = "__sm/v1/opts_hash";
constexpr const char* KEY_INPUTS      = "__sm/v1/inputs";
constexpr const char* KEY_PHASE       = "__sm/v1/phase";
constexpr const char* KEY_DENOVO_DONE = "__sm/v1/denovo_done";

// Prefix-namespaced keys (key builders below).
constexpr const char* PFX_ALIGN_DONE   = "__sm/v1/align_done/";
constexpr const char* PFX_THREAD_DONE  = "__sm/v1/thread_done/";
constexpr const char* PFX_REPORT_DONE  = "__sm/v1/report_done/";

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
std::string key_thread_done_prefix(uint16_t i, uint16_t p);
std::string key_thread_done(uint16_t i, uint16_t p, uint32_t slot);
std::string key_report_done(const std::string& kind);

// Encode/decode of the ThreadProgress blob stored under thread_done/*.
// Encoding is fixed-width little-endian: 12 bytes total.
std::string         encode_thread_progress(const ThreadProgress& tp);
std::optional<ThreadProgress> decode_thread_progress(const std::string& s);

// Snapshot of what probe_or_init found in the kvdb. Drives the resume logic
// in align() and the higher-level main() phase dispatch.
struct State {
	// false: fresh kvdb (no version sentinel); true: existing kvdb whose
	// opts_hash + inputs matched the current run, so we may resume.
	bool is_resume = false;
	Phase phase = Phase::Init;
	// (idx_num, idx_part) of passes already committed (align_done present).
	std::vector<std::pair<uint16_t, uint16_t>> align_done;
	// Latest (i, p) committed via commit_pass; -1/-1 if none.
	std::pair<int, int> latest_pass = {-1, -1};
	// At most one partial pass to resume on restart: the pass with thread_done
	// entries that has no matching align_done. Detected via PFX_THREAD_DONE.
	std::pair<int, int> resume_pass = {-1, -1};
	// Per-slot progress for resume_pass, indexed by slot id (i*num_sense + j).
	// std::nullopt means the slot had no entry stored — its counter is 0 and
	// it should be read from bytes_start.
	std::vector<std::optional<ThreadProgress>> resume_progress;
	std::string run_id;
};

// On open: stamp run identity if fresh, else load existing progress.
// Throws std::runtime_error on opts_hash / inputs mismatch.
State probe_or_init(KeyValueDatabase& kvdb, const Runopts& opts);

// Progress writes.
//
// commit_pass atomically writes the align_done marker, the Readstats blob
// (stored under readstats_dbkey, which is Readstats::dbkey), and clears all
// thread_done/{i}/{p}/* keys for the just-finished pass in a single WriteBatch.
// "align_done present" implies "Readstats blob reflects pass (i,p) and all
// prior passes; thread_done entries for (i,p) are stale".
void commit_pass(KeyValueDatabase& kvdb, uint16_t i, uint16_t p,
                 const std::string& readstats_dbkey,
                 const std::string& readstats_blob);

// Atomic per-read commit: store the read blob AND advance the thread_done
// counter for the worker's slot. Used by align2 whenever a new hit causes a
// kvdb write — the counter advances in the same WriteBatch so a crash between
// the two writes cannot leave them mismatched.
void put_read_with_progress(KeyValueDatabase& kvdb,
                            const std::string& read_id,
                            const std::string& read_blob,
                            uint16_t i, uint16_t p, uint32_t slot,
                            const ThreadProgress& tp);

// Periodic counter-only commit: when no read blob is being written, the worker
// still needs to advance its durable counter. Standalone Put, not a batch.
void put_progress(KeyValueDatabase& kvdb,
                  uint16_t i, uint16_t p, uint32_t slot,
                  const ThreadProgress& tp);

void mark_align_done(KeyValueDatabase& kvdb, uint16_t i, uint16_t p);
void clear_thread_done(KeyValueDatabase& kvdb, uint16_t i, uint16_t p);

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
