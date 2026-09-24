/*
 * smr_api.cpp -- SortMeRNA reentrant C API implementation
 */

#include "smr_api.h"
#include "version.h"
#include "options.hpp"
#include "index.hpp"
#include "kvdb.hpp"
#include "readfeed.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "references.hpp"
#include "read.hpp"
#include "processor.hpp"
#include "output.hpp"
#include "summary.hpp"
#include "otumap.h"
#include "restart.hpp"

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <atomic>
#include <mutex>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

/* stringification helpers for version macros */
#define SMR_STRINGIFY2(x) #x
#define SMR_STRINGIFY(x) SMR_STRINGIFY2(x)

/* process-level atomic counter for unique workdir names across all contexts */
static std::atomic<int> g_run_counter{0};

/* process-level mutex serializing smr_run calls. Required because the
 * dup2-based fd suppression is process-wide and not thread-safe. The
 * mutex is uncontended in the common single-threaded case. */
static std::mutex g_run_mutex;

/* Pipeline log messages come from the calling thread and from SortMeRNA's
 * worker threads (they inherit the thread-local callback declared in
 * common.hpp). Every call of a caller's log_callback takes this mutex, so
 * the callback never runs concurrently. */
static std::mutex g_log_mutex;

/* --- Internal context definition (opaque to callers) --- */

struct smr_context {
    smr_config_t config;
    char last_error[1024];
    int last_error_code;
};

static void set_error(smr_context *ctx, int code, const char *fmt, ...) {
    if (!ctx) return;
    ctx->last_error_code = code;
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(ctx->last_error, sizeof(ctx->last_error), fmt, ap);
    va_end(ap);
}

static void ctx_log(const smr_context *ctx, int level, const char *fmt, ...) {
    if (!ctx || !ctx->config.log_callback) return;
    char buf[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    std::lock_guard<std::mutex> lk(g_log_mutex);
    ctx->config.log_callback(level, buf, ctx->config.log_user_data);
}

/* Internal RAII helpers — anonymous namespace gives them TU-local
 * linkage so their type symbols cannot collide with same-named types
 * in any other translation unit. */
namespace {

/* RAII guard for temp workdir cleanup */
class WorkdirGuard {
    std::string path_;
    bool owned_;
public:
    WorkdirGuard(const std::string &path, bool owned)
        : path_(path), owned_(owned) {}
    ~WorkdirGuard() {
        if (owned_) {
            std::error_code ec;
            std::filesystem::remove_all(path_, ec);
        }
    }
    void release() { owned_ = false; }
};

/* RAII guard for stdout/stderr suppression via dup2.
 * Best-effort: if dup/open fails, silently proceeds without suppression.
 * NOTE: dup2 is process-wide — not safe for concurrent smr_run() calls
 * from multiple threads. The macro-based log routing (thread-local) is
 * the primary isolation mechanism; this is a secondary defense against
 * direct std::cout writes that bypass the macros.
 *
 * Nested instances on the same thread are no-ops: only the outermost
 * guard installs/restores the fd redirect. This lets smr_run_seqs call
 * smr_index_load + smr_run_seqs_with_index without double-redirecting. */
class FdRedirectGuard {
    int saved_out_, saved_err_;
    bool outer_;
    static thread_local int depth_;
public:
    FdRedirectGuard() : saved_out_(-1), saved_err_(-1), outer_(depth_ == 0) {
        ++depth_;
        if (!outer_) return;
        fflush(stdout); fflush(stderr);
        saved_out_ = dup(STDOUT_FILENO);
        saved_err_ = dup(STDERR_FILENO);
        int devnull = open("/dev/null", O_WRONLY);
        if (devnull >= 0) {
            if (saved_out_ >= 0) dup2(devnull, STDOUT_FILENO);
            if (saved_err_ >= 0) dup2(devnull, STDERR_FILENO);
            close(devnull);
        } else {
            /* open failed — clean up saved fds to avoid leak */
            if (saved_out_ >= 0) { close(saved_out_); saved_out_ = -1; }
            if (saved_err_ >= 0) { close(saved_err_); saved_err_ = -1; }
        }
    }
    ~FdRedirectGuard() {
        --depth_;
        if (!outer_) return;
        fflush(stdout); fflush(stderr);
        if (saved_out_ >= 0) { dup2(saved_out_, STDOUT_FILENO); close(saved_out_); }
        if (saved_err_ >= 0) { dup2(saved_err_, STDERR_FILENO); close(saved_err_); }
    }
};
thread_local int FdRedirectGuard::depth_ = 0;

void discard_log(int, const char *, void *) {}

struct LogTarget {
    smr_log_fn cb;
    void *ud;
};

void forward_log(int level, const char *msg, void *target) {
    auto *t = static_cast<LogTarget *>(target);
    std::lock_guard<std::mutex> lk(g_log_mutex);
    t->cb(level, msg, t->ud);
}

/* RAII guard for the thread-local log callback: routes pipeline logging to
 * the caller's callback through forward_log for the duration of an API
 * call, and restores the previous routing afterwards. Worker threads started
 * during the call inherit it and are joined before the call returns, so
 * target_ outlives them. With no callback, log output is discarded rather
 * than written to stdout, so the library stays silent also outside the fd
 * redirect (smr_index_load loads the index and references without it). */
class LogRouteGuard {
    LogTarget target_;
    smr_log_fn prev_cb_;
    void *prev_ud_;
public:
    LogRouteGuard(smr_log_fn cb, void *ud)
        : target_{cb, ud}, prev_cb_(smr_tl_log_callback), prev_ud_(smr_tl_log_user_data) {
        smr_tl_log_callback = cb ? forward_log : discard_log;
        smr_tl_log_user_data = cb ? &target_ : nullptr;
    }
    ~LogRouteGuard() {
        smr_tl_log_callback = prev_cb_;
        smr_tl_log_user_data = prev_ud_;
    }
    LogRouteGuard(const LogRouteGuard &) = delete;
    LogRouteGuard &operator=(const LogRouteGuard &) = delete;
};

/* Deletes the kvdb directory on destruction. Declare it before the
 * KeyValueDatabase so the database is closed first. */
class KvdbScratchGuard {
    std::filesystem::path path_;
public:
    explicit KvdbScratchGuard(std::filesystem::path path) : path_(std::move(path)) {}
    ~KvdbScratchGuard() {
        std::error_code ec;
        std::filesystem::remove_all(path_, ec);
    }
};

} // anonymous namespace

/* --- Configuration --- */

void smr_config_init(smr_config_t *cfg) {
    if (!cfg) return;
    memset(cfg, 0, sizeof(*cfg));
    cfg->struct_size = (uint32_t)sizeof(*cfg);

    /* threading */
    cfg->num_threads = 2;          /* Runopts::num_proc_thread */

    /* alignment scoring -- matches Runopts defaults */
    cfg->match    =  2;            /* Runopts::match */
    cfg->mismatch = -3;            /* Runopts::mismatch */
    cfg->gap_open =  5;            /* Runopts::gap_open */
    cfg->gap_ext  =  2;            /* Runopts::gap_extension */
    cfg->score_N  = -3;            /* Runopts default: equals mismatch */
    cfg->evalue   = -1.0;          /* Runopts::evalue (-1 = unset) */

    /* indexing */
    cfg->seed_win_len   = 18;      /* Runopts::seed_win_len */
    cfg->num_alignments =  1;      /* Runopts::num_alignments */

    /* boolean flags */
    cfg->best = 1;                 /* Runopts::is_best */
}

/* --- Context lifecycle --- */

smr_context_t *smr_ctx_create(const smr_config_t *cfg) {
    if (!cfg) return nullptr;
    if (cfg->struct_size != (uint32_t)sizeof(smr_config_t)) return nullptr;

    auto *ctx = static_cast<smr_context_t *>(malloc(sizeof(smr_context_t)));
    if (!ctx) return nullptr;

    memset(ctx, 0, sizeof(*ctx));
    memcpy(&ctx->config, cfg, sizeof(smr_config_t));

    ctx_log(ctx, SMR_LOG_INFO, "context created");
    return ctx;
}

void smr_ctx_destroy(smr_context_t *ctx) {
    if (!ctx) return; /* NULL is always safe -- contract */
    /* Log before any resource teardown so callback fires against a live context */
    ctx_log(ctx, SMR_LOG_INFO, "context destroyed");
    free(ctx);
}

/* --- Error reporting --- */

const char *smr_strerror(int code) {
    switch (code) {
    case SMR_OK:                  return "Success";
    case SMR_ERR_INVALID_CONFIG:  return "Invalid configuration";
    case SMR_ERR_ALLOC:           return "Memory allocation failed";
    case SMR_ERR_IO:              return "I/O error";
    case SMR_ERR_INDEX:           return "Index error";
    case SMR_ERR_ALIGN:           return "Alignment error";
    case SMR_ERR_NOT_IMPLEMENTED: return "Not implemented";
    default:                      return "Unknown error";
    }
}

const char *smr_last_error(const smr_context_t *ctx) {
    if (!ctx) return "";
    return ctx->last_error;
}

int smr_last_error_code(const smr_context_t *ctx) {
    if (!ctx) return SMR_ERR_INVALID_CONFIG;
    return ctx->last_error_code;
}

/* --- Computation helpers --- */

static bool file_exists(const char *path) {
    struct stat st;
    return stat(path, &st) == 0;
}

static bool file_is_empty(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0) return true;
    return st.st_size == 0;
}

/* Pick a unique temporary workdir under the system temp directory (honours
 * TMPDIR). Returns false and sets the context error on failure. */
static bool make_temp_workdir(smr_context *ctx, std::string &workdir) {
    std::error_code ec;
    auto tmp = std::filesystem::temp_directory_path(ec);
    if (ec) {
        set_error(ctx, SMR_ERR_IO, "no temporary directory: %s", ec.message().c_str());
        return false;
    }
    std::ostringstream ss;
    ss << "smr_api_" << getpid() << "_" << g_run_counter.fetch_add(1);
    workdir = (tmp / ss.str()).string();
    return true;
}

/*
 * The kvdb holds the alignment results of one run. The library never carries
 * it over between calls (see smr_config_t.workdir): a leftover kvdb would be
 * resumed by the CLI restart logic or return an earlier run's reads.
 */
static void clear_kvdb(const std::filesystem::path &kvdbdir) {
    std::error_code ec;
    std::filesystem::remove_all(kvdbdir, ec);
    if (ec)
        throw std::runtime_error("failed to remove " + kvdbdir.string() + ": " + ec.message());
}

/*
 * Build an argv vector from smr_config_t + input paths.
 * Always emits all parameters unconditionally to avoid coupling
 * with Runopts internal defaults.
 */
static std::vector<std::string> build_argv(
    const smr_config_t &cfg,
    const char **ref_paths, int32_t num_refs,
    const char **read_paths, int32_t num_reads,
    const std::string &workdir)
{
    std::vector<std::string> args;
    args.push_back("sortmerna");

    for (int32_t i = 0; i < num_refs; i++) {
        args.push_back("--ref");
        args.push_back(ref_paths[i]);
    }
    for (int32_t i = 0; i < num_reads; i++) {
        args.push_back("--reads");
        args.push_back(read_paths[i]);
    }

    args.push_back("--workdir");
    args.push_back(workdir);
    args.push_back("--aligned");
    args.push_back(workdir + "/aligned");
    args.push_back("--other");
    args.push_back(workdir + "/other");

    args.push_back("--threads");
    args.push_back(std::to_string(cfg.num_threads));
    args.push_back("--num_alignments");
    args.push_back(std::to_string(cfg.num_alignments));

    /* always produce BLAST and SAM output for structured result extraction */
    args.push_back("--blast");
    args.push_back("1 cigar qcov");
    args.push_back("--fastx");
    args.push_back("--sam");

    /* always emit all scoring parameters */
    args.push_back("--match");
    args.push_back(std::to_string(cfg.match));
    args.push_back("--mismatch");
    args.push_back(std::to_string(cfg.mismatch));
    args.push_back("--gap_open");
    args.push_back(std::to_string(cfg.gap_open));
    args.push_back("--gap_ext");
    args.push_back(std::to_string(cfg.gap_ext));
    args.push_back("-N");
    args.push_back(std::to_string(cfg.score_N));
    args.push_back("-L");
    args.push_back(std::to_string(cfg.seed_win_len));

    if (cfg.evalue >= 0.0) {
        /* %.17g round-trips; std::to_string would print 1e-10 as 0.000000 */
        char buf[32];
        snprintf(buf, sizeof(buf), "%.17g", cfg.evalue);
        args.push_back("-e");
        args.push_back(buf);
    }

    if (cfg.forward_only) args.push_back("-F");
    if (cfg.reverse_only) args.push_back("-R");
    if (!cfg.best) args.push_back("--no-best");
    if (cfg.full_search) args.push_back("--full_search");
    if (cfg.paired) args.push_back("--paired_in");

    return args;
}

/*
 * Build a CIGAR string from an s_align2's packed cigar vector + read coords.
 * Matches the format produced by report_blast.cpp.
 */
static std::string build_cigar_string(const s_align2 &align, uint32_t readlen) {
    std::string cigar;
    if (align.read_begin1 != 0)
        cigar += std::to_string(align.read_begin1) + "S";
    for (uint32_t c = 0; c < align.cigar.size(); ++c) {
        uint32_t letter = 0xf & align.cigar[c];
        uint32_t length = (0xfffffff0 & align.cigar[c]) >> 4;
        cigar += std::to_string(length);
        if (letter == 0) cigar += "M";
        else if (letter == 1) cigar += "I";
        else cigar += "D";
    }
    auto end_mask = readlen - align.read_end1 - 1;
    if (end_mask > 0)
        cigar += std::to_string(end_mask) + "S";
    return cigar;
}

/*
 * Extract the read identifier (first token of header, without > or @).
 */
static std::string extract_read_id(const std::string &header) {
    size_t start = (header.size() > 0 && (header[0] == '>' || header[0] == '@')) ? 1 : 0;
    auto end = header.find_first_of(" \t", start);
    return header.substr(start, end - start);
}

/*
 * Populate per-read output arrays by iterating through readfeed and loading
 * alignment results from kvdb. Loads references to compute %ID, %COV, E-value
 * (same pattern as the report and denovo_stats phases).
 */
static bool populate_per_read_output(smr_output_t *o,
                                     Readfeed &readfeed,
                                     Readstats &readstats,
                                     KeyValueDatabase &kvdb,
                                     Runopts &opts,
                                     References *preloaded_refs,
                                     Refstats *preloaded_refstats) {
    uint64_t n = o->num_reads;
    if (n == 0) return true;

    o->read_ids   = static_cast<const char**>(calloc(n, sizeof(char*)));
    o->aligned    = static_cast<int32_t*>(calloc(n, sizeof(int32_t)));
    o->ref_index  = static_cast<int32_t*>(calloc(n, sizeof(int32_t)));
    o->e_value    = static_cast<double*>(calloc(n, sizeof(double)));
    o->identity   = static_cast<double*>(calloc(n, sizeof(double)));
    o->coverage   = static_cast<double*>(calloc(n, sizeof(double)));
    o->ref_start  = static_cast<int32_t*>(calloc(n, sizeof(int32_t)));
    o->ref_end    = static_cast<int32_t*>(calloc(n, sizeof(int32_t)));
    o->cigar      = static_cast<const char**>(calloc(n, sizeof(char*)));
    o->ref_name   = static_cast<const char**>(calloc(n, sizeof(char*)));
    o->strand     = static_cast<int32_t*>(calloc(n, sizeof(int32_t)));
    o->score      = static_cast<int32_t*>(calloc(n, sizeof(int32_t)));
    o->edit_distance = static_cast<int32_t*>(calloc(n, sizeof(int32_t)));

    if (!o->read_ids || !o->aligned || !o->ref_index || !o->e_value ||
        !o->identity || !o->coverage || !o->ref_start || !o->ref_end || !o->cigar || !o->ref_name ||
        !o->strand || !o->score || !o->edit_distance)
        return false;

    /* defaults for all reads (unaligned sentinel) */
    for (uint64_t i = 0; i < n; i++) {
        o->ref_index[i]     = -1;
        o->strand[i]        = -1;
        o->score[i]         = -1;
        o->edit_distance[i] = -1;
    }

    /* Use preloaded Refstats/References from the handle if supplied, else
     * construct locally (file-based smr_run path). */
    std::unique_ptr<Refstats> local_refstats_owner;
    Refstats *refstats_ptr = preloaded_refstats;
    if (!refstats_ptr) {
        local_refstats_owner.reset(new Refstats(opts, readstats));
        refstats_ptr = local_refstats_owner.get();
    }
    Refstats &refstats = *refstats_ptr;
    References local_refs;
    References &refs = preloaded_refs ? *preloaded_refs : local_refs;
    unsigned num_parts = readfeed.num_splits > 0 ? readfeed.num_splits : 1;
    unsigned num_sense = readfeed.num_sense;

    /*
     * Visit every read in input order. Both the file feed and the memory feed
     * give each processing thread a contiguous chunk of the input, so input
     * order is part 0 to the end, then part 1, and so on. Paired reads are
     * split over the part's two slots (fwd = p*2, rev = p*2+1), and reading
     * them alternately restores the interleaved order, as the reports do.
     */
    auto for_each_read = [&](auto callback) {
        readfeed.init_reading();
        uint64_t out_idx = 0;
        std::string readstr;
        for (unsigned p = 0; p < num_parts && out_idx < n; ++p) {
            for (bool part_done = false; !part_done && out_idx < n; ) {
                for (unsigned s = 0; s < num_sense && out_idx < n; ++s) {
                    readstr.clear();
                    if (!readfeed.next(static_cast<int>(p * num_sense + s), readstr)) {
                        part_done = true;
                        break;
                    }
                    Read rd(readstr);
                    callback(out_idx, rd);
                    ++out_idx;
                }
            }
        }
    };

    /* Pass 1: collect read IDs and basic alignment info from kvdb (no refs
     * needed). Retain the Read objects for Pass 2 reuse so Pass 2 doesn't
     * re-issue kvdb.get() for the same keys. */
    std::vector<Read> reads_cache;
    reads_cache.reserve(n);
    uint64_t aligned_count = 0;
    for_each_read([&](uint64_t idx, Read &read) {
        read.init(opts);
        read.load_db(kvdb);

        std::string rid = extract_read_id(read.header);
        o->read_ids[idx] = strdup(rid.c_str());

        if (read.is_hit && !read.alignment.alignv.empty()) {
            const auto &align = read.alignment.alignv[0];
            o->aligned[idx]   = 1;
            o->ref_index[idx] = static_cast<int32_t>(align.index_num);
            o->ref_start[idx] = align.ref_begin1 + 1;
            o->ref_end[idx]   = align.ref_end1 + 1;
            o->strand[idx]    = static_cast<int32_t>(align.strand);
            o->score[idx]     = static_cast<int32_t>(align.score1); /* uint16_t, always fits int32_t */

            /* CIGAR string */
            uint32_t readlen = static_cast<uint32_t>(read.sequence.size());

            if (opts.is_library_mode) {
                /* E-value: textbook per-query Karlin-Altschul form,
                 *     E = K * m * n * exp(-λ * S)
                 * where n is THIS read's length and m is the uncorrected DB
                 * size. Both inputs must be batch-invariant, hence
                 * full_ref_raw (from the .stats file) rather than full_ref
                 * (which has an edge-effect correction m' = m - expect_L*numseq
                 * that depends on batch-level Readstats aggregates). Using the
                 * raw m biases library e-values ~1-3% lower than the
                 * edge-corrected form for SILVA-scale DBs (expect_L is ~20-40
                 * nt, numseq*expect_L / m ~ 0.01-0.03). This is on top of the
                 * per-query vs run-aggregate n divergence below.
                 *
                 * Two distinct sources of library-vs-CLI e-value difference:
                 *   (1) n: per-read length (library) vs summed query space (CLI).
                 *       Library e-values smaller by ~factor of batch size.
                 *   (2) m: raw DB length (library) vs edge-corrected (CLI).
                 *       Library e-values smaller by ~1-3% for SILVA-scale DBs.
                 * Both are intentional — (1) is required for batch invariance,
                 * (2) for strict batch invariance of m (since expect_L depends
                 * on the batch-level full_read). Callers filtering on e-value
                 * thresholds should be aware and calibrate accordingly. */
                o->e_value[idx] = static_cast<double>(refstats.gumbel[align.index_num].second)
                    * static_cast<double>(refstats.full_ref_raw[align.index_num])
                    * static_cast<double>(readlen)
                    * std::exp(-refstats.gumbel[align.index_num].first * align.score1);
            } else {
                /* smr_run: the CLI definition, as in report_blast.cpp */
                o->e_value[idx] = (double)refstats.gumbel[align.index_num].second
                    * refstats.full_ref[align.index_num]
                    * refstats.full_read[align.index_num]
                    * std::exp(-refstats.gumbel[align.index_num].first * align.score1);
            }
            o->cigar[idx] = strdup(build_cigar_string(align, readlen).c_str());

            /* coverage from alignment coordinates (no refs needed) */
            if (readlen > 0) {
                o->coverage[idx] = static_cast<double>(align.read_end1 - align.read_begin1 + 1)
                                 / readlen * 100.0;
            }

            ++aligned_count;
        }
        reads_cache.push_back(std::move(read));
    });
    o->num_aligned = aligned_count;

    /* Pass 2: compute %ID / %COV / ref_name for aligned reads (requires
     * loaded refs for sequence comparison). Uses the Read objects retained
     * from Pass 1 — no re-issued kvdb lookups. */
    auto pass2_body = [&](uint64_t ri, Read &rd, uint16_t active_ref_idx, uint16_t active_part) {
        if (o->aligned[ri] != 1) return;
        if (rd.alignment.alignv.empty()) return;
        const auto &al = rd.alignment.alignv[0];
        if (al.index_num != active_ref_idx || al.part != active_part) return;
        if (rd.is03) rd.flip34();
        if (al.strand == rd.reversed) rd.revIntStr();
        auto mgm = rd.calc_miss_gap_match(refs, al);
        o->identity[ri] = std::get<3>(mgm) * 100.0;
        o->coverage[ri] = std::get<4>(mgm) * 100.0;
        uint32_t ed = std::get<0>(mgm) + std::get<1>(mgm);
        o->edit_distance[ri] = ed > (uint32_t)INT32_MAX ? INT32_MAX : static_cast<int32_t>(ed);
        o->ref_name[ri] = strdup(refs.buffer[al.ref_num].id.c_str());
    };

    if (preloaded_refs) {
        /* Library path: refs already loaded for (refs.num, refs.part) by the
         * handle. Iterate reads once and use them directly; no load/unload. */
        for (uint64_t ri = 0; ri < reads_cache.size(); ++ri) {
            pass2_body(ri, reads_cache[ri], refs.num, static_cast<uint16_t>(refs.part));
        }
    } else {
        /* File-based path: load each (ref_idx, part) in turn, run the body,
         * unload. Each read's alignment belongs to exactly one (ref_idx,
         * part); the filter in pass2_body skips the others. */
        for (uint16_t ref_idx = 0; ref_idx < opts.indexfiles.size(); ++ref_idx) {
            for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[ref_idx]; ++idx_part) {
                refs.load(ref_idx, idx_part, opts, refstats);
                for (uint64_t ri = 0; ri < reads_cache.size(); ++ri) {
                    pass2_body(ri, reads_cache[ri], ref_idx, idx_part);
                }
                refs.unload();
            }
        }
    }

    return true;
}

int smr_run(smr_context_t *ctx,
            const char **ref_paths, int32_t num_refs,
            const char **read_paths, int32_t num_reads,
            smr_output_t **out,
            smr_stats_t *stats) {
    if (!ctx) return SMR_ERR_INVALID_CONFIG;

    /* validate inputs */
    if (!ref_paths || num_refs <= 0) {
        set_error(ctx, SMR_ERR_INVALID_CONFIG, "ref_paths is NULL or num_refs <= 0");
        return SMR_ERR_INVALID_CONFIG;
    }
    if (!read_paths || num_reads <= 0) {
        set_error(ctx, SMR_ERR_INVALID_CONFIG, "read_paths is NULL or num_reads <= 0");
        return SMR_ERR_INVALID_CONFIG;
    }

    /* validate file existence */
    for (int32_t i = 0; i < num_refs; i++) {
        if (!ref_paths[i] || !file_exists(ref_paths[i])) {
            set_error(ctx, SMR_ERR_IO, "reference file not found: %s",
                      ref_paths[i] ? ref_paths[i] : "(null)");
            return SMR_ERR_IO;
        }
        if (file_is_empty(ref_paths[i])) {
            set_error(ctx, SMR_ERR_IO, "reference file is empty: %s", ref_paths[i]);
            return SMR_ERR_IO;
        }
    }
    for (int32_t i = 0; i < num_reads; i++) {
        if (!read_paths[i] || !file_exists(read_paths[i])) {
            set_error(ctx, SMR_ERR_IO, "reads file not found: %s",
                      read_paths[i] ? read_paths[i] : "(null)");
            return SMR_ERR_IO;
        }
        if (file_is_empty(read_paths[i])) {
            set_error(ctx, SMR_ERR_IO, "reads file is empty: %s", read_paths[i]);
            return SMR_ERR_IO;
        }
    }

    /* serialize smr_run calls — dup2 fd suppression is process-wide */
    std::lock_guard<std::mutex> run_lock(g_run_mutex);

    auto wall_start = std::chrono::high_resolution_clock::now();

    std::string workdir;
    bool workdir_is_temp = (ctx->config.workdir == nullptr);
    if (!workdir_is_temp) {
        workdir = ctx->config.workdir;
    } else if (!make_temp_workdir(ctx, workdir)) {
        return SMR_ERR_IO;
    }

    /* RAII guards — C++ destroys in reverse declaration order, so:
     *   declared first:  WorkdirGuard  → destroyed last  (cleanup temp dir)
     *   declared second: LogRouteGuard → destroyed second (clear thread-local)
     *   declared third:  FdRedirectGuard → destroyed first (restore fds)
     * This means fds are restored while log routing is still active,
     * so any teardown logging still has somewhere to go. */
    WorkdirGuard wdguard(workdir, workdir_is_temp);
    LogRouteGuard loguard(ctx->config.log_callback, ctx->config.log_user_data);
    FdRedirectGuard fdguard;

    try {
        /* build argv and construct Runopts */
        auto args = build_argv(ctx->config, ref_paths, num_refs, read_paths, num_reads, workdir);
        std::vector<char*> argv_ptrs;
        for (auto &a : args) argv_ptrs.push_back(const_cast<char*>(a.c_str()));
        argv_ptrs.push_back(nullptr);

        bool dryrun = false;
        Runopts opts(static_cast<int>(argv_ptrs.size() - 1), argv_ptrs.data(), dryrun);

        ctx_log(ctx, SMR_LOG_INFO, "pipeline starting: %d refs, %d reads", num_refs, num_reads);

        /* run the alignment pipeline (same as main.cpp) */
        Index index(opts);

        /* fresh kvdb, removed again when this scope ends (after kvdb closes) */
        clear_kvdb(opts.kvdbdir);
        KvdbScratchGuard kvdb_guard(opts.kvdbdir);
        KeyValueDatabase kvdb(opts.kvdbdir.string());
        restart::State rstate = restart::probe_or_init(kvdb, opts);
        Readfeed readfeed(opts.feed_type, opts.readfiles, opts.num_proc_thread, opts.readb_dir, opts.is_paired);
        Readstats readstats(readfeed.num_reads_tot, readfeed.length_all,
                            readfeed.min_read_len, readfeed.max_read_len, kvdb, opts);

        /* align + report */
        align(readfeed, readstats, index, kvdb, opts, &rstate);
        writeSummary(readstats, opts);
        writeReports(readfeed, readstats, kvdb, opts);

        ctx_log(ctx, SMR_LOG_INFO, "pipeline complete");

        /* populate stats from Readstats (authoritative source) */
        if (stats) {
            memset(stats, 0, sizeof(*stats));
            stats->total_reads = readstats.all_reads_count;
            stats->total_aligned = readstats.num_aligned.load();
            stats->total_id_cov_pass = readstats.n_yid_ycov.load();
            stats->total_denovo = readstats.num_denovo.load();
            stats->min_read_len = readstats.min_read_len;
            stats->max_read_len = readstats.max_read_len;
            auto wall_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = wall_end - wall_start;
            stats->wall_time_sec = elapsed.count();
        }

        /* populate output */
        if (out) {
            auto *o = static_cast<smr_output_t *>(calloc(1, sizeof(smr_output_t)));
            if (!o) {
                set_error(ctx, SMR_ERR_ALLOC, "failed to allocate smr_output_t");
                return SMR_ERR_ALLOC;
            }

            o->num_reads = readstats.all_reads_count;

            if (!populate_per_read_output(o, readfeed, readstats, kvdb, opts, nullptr, nullptr)) {
                smr_output_free(o);
                set_error(ctx, SMR_ERR_ALLOC, "failed to allocate per-read output arrays");
                return SMR_ERR_ALLOC;
            }

            *out = o;
        }

        /* RAII guards handle cleanup on return */
        set_error(ctx, SMR_OK, "");
        return SMR_OK;

    } catch (const smr_exit_requested &) {
        set_error(ctx, SMR_ERR_INVALID_CONFIG, "unexpected --help/--version in library context");
        return SMR_ERR_INVALID_CONFIG;
    } catch (const std::exception &e) {
        set_error(ctx, SMR_ERR_ALIGN, "%s", e.what());
        return SMR_ERR_ALIGN;
    } catch (...) {
        set_error(ctx, SMR_ERR_ALIGN, "unknown exception in pipeline");
        return SMR_ERR_ALIGN;
    }
}

/* --- Pre-loaded index (streaming) API --- */

/* Opaque handle. Owns workdir lifecycle, the cached Runopts, the loaded
 * Index+References+Refstats for a single (ref_idx=0, part=0), and an open
 * KeyValueDatabase. The expensive burst-trie + reference-sequence load is
 * paid exactly once per handle in smr_index_load; smr_run_seqs_with_index
 * only constructs per-batch Readfeed + Readstats and dispatches to
 * align_loaded(). Multiple reference files and multi-part indexes are not
 * supported yet. */
struct smr_index {
    smr_context_t *ctx;
    std::vector<std::string> ref_path_storage;
    std::string workdir;
    bool workdir_is_temp;

    std::unique_ptr<Runopts> opts;
    std::unique_ptr<KeyValueDatabase> kvdb;
    std::unique_ptr<Refstats> refstats;
    std::unique_ptr<Index> index;
    References references;   /* not unique_ptr — no default ctor needed */
    bool index_loaded;       /* true once index.load() + references.load() succeeded */
};

static int validate_seqs(smr_context_t *ctx, const smr_seq_t *seqs, int32_t num_seqs,
                         bool paired, bool *has_qual_out) {
    *has_qual_out = false;
    if (!seqs || num_seqs <= 0) {
        set_error(ctx, SMR_ERR_INVALID_CONFIG, "seqs is NULL or num_seqs <= 0");
        return SMR_ERR_INVALID_CONFIG;
    }
    if (paired && (num_seqs % 2 != 0)) {
        set_error(ctx, SMR_ERR_INVALID_CONFIG,
                  "paired mode requires even num_seqs (interleaved fwd/rev); got %d", num_seqs);
        return SMR_ERR_INVALID_CONFIG;
    }
    bool has_qual = (seqs[0].quality != nullptr);
    for (int32_t i = 0; i < num_seqs; i++) {
        if (!seqs[i].id || !seqs[i].sequence) {
            set_error(ctx, SMR_ERR_INVALID_CONFIG, "seq[%d] has NULL id or sequence", i);
            return SMR_ERR_INVALID_CONFIG;
        }
        if (seqs[i].sequence[0] == '\0') {
            set_error(ctx, SMR_ERR_INVALID_CONFIG, "seq[%d] has empty sequence", i);
            return SMR_ERR_INVALID_CONFIG;
        }
        if (seqs[i].quality && strlen(seqs[i].quality) != strlen(seqs[i].sequence)) {
            set_error(ctx, SMR_ERR_INVALID_CONFIG,
                      "seq[%d] quality length != sequence length", i);
            return SMR_ERR_INVALID_CONFIG;
        }
        bool this_has_qual = (seqs[i].quality != nullptr);
        if (this_has_qual != has_qual) {
            set_error(ctx, SMR_ERR_INVALID_CONFIG,
                      "seq[%d]: all sequences must have quality strings or none", i);
            return SMR_ERR_INVALID_CONFIG;
        }
    }
    *has_qual_out = has_qual;
    return SMR_OK;
}

smr_index_t *smr_index_load(smr_context_t *ctx,
                            const char **ref_paths, int32_t num_refs) {
    if (!ctx) return nullptr;

    if (!ref_paths || num_refs <= 0) {
        set_error(ctx, SMR_ERR_INVALID_CONFIG, "ref_paths is NULL or num_refs <= 0");
        return nullptr;
    }
    /* Each reference file is a separate index, and the handle holds only
     * one loaded index. Concatenate the FASTA files instead. */
    if (num_refs > 1) {
        set_error(ctx, SMR_ERR_NOT_IMPLEMENTED,
                  "smr_index_load supports one reference file (got %d); concatenate them", num_refs);
        return nullptr;
    }

    for (int32_t i = 0; i < num_refs; i++) {
        if (!ref_paths[i] || !file_exists(ref_paths[i])) {
            set_error(ctx, SMR_ERR_IO, "reference file not found: %s",
                      ref_paths[i] ? ref_paths[i] : "(null)");
            return nullptr;
        }
        if (file_is_empty(ref_paths[i])) {
            set_error(ctx, SMR_ERR_IO, "reference file is empty: %s", ref_paths[i]);
            return nullptr;
        }
    }

    auto *idx = new (std::nothrow) smr_index_t;
    if (!idx) {
        set_error(ctx, SMR_ERR_ALLOC, "failed to allocate smr_index_t");
        return nullptr;
    }
    idx->ctx = ctx;
    idx->index_loaded = false;
    idx->ref_path_storage.reserve(num_refs);
    for (int32_t i = 0; i < num_refs; i++) {
        idx->ref_path_storage.emplace_back(ref_paths[i]);
    }

    idx->workdir_is_temp = (ctx->config.workdir == nullptr);
    if (!idx->workdir_is_temp) {
        idx->workdir = ctx->config.workdir;
    } else if (!make_temp_workdir(ctx, idx->workdir)) {
        smr_index_free(idx);
        return nullptr;
    }
    std::error_code ec;
    std::filesystem::create_directories(idx->workdir, ec);
    if (ec) {
        set_error(ctx, SMR_ERR_IO, "failed to create workdir %s: %s",
                  idx->workdir.c_str(), ec.message().c_str());
        smr_index_free(idx);
        return nullptr;
    }

    LogRouteGuard loguard(ctx->config.log_callback, ctx->config.log_user_data);

    try {
        /* Write a FASTA placeholder_reads file for Runopts. The real reads
         * come via in-memory Readfeed at each with_index call; this file is
         * never read for sequence data. See header "Input format" note. */
        std::string placeholder_reads = idx->workdir + "/placeholder_reads.fa";
        {
            std::ofstream ofs(placeholder_reads);
            ofs << ">placeholder\nA\n";
        }
        const char *dummy_read_paths[] = { placeholder_reads.c_str() };

        std::vector<const char *> ref_ptrs;
        ref_ptrs.reserve(idx->ref_path_storage.size());
        for (auto &s : idx->ref_path_storage) ref_ptrs.push_back(s.c_str());

        auto args = build_argv(ctx->config, ref_ptrs.data(), (int32_t)ref_ptrs.size(),
                               dummy_read_paths, 1, idx->workdir);
        std::vector<char*> argv_ptrs;
        for (auto &a : args) argv_ptrs.push_back(const_cast<char*>(a.c_str()));
        argv_ptrs.push_back(nullptr);

        /* Process-wide-stateful ops (Runopts parse + dup2 suppression) are
         * serialized; heavy loads below run concurrently across handles. */
        {
            std::lock_guard<std::mutex> run_lock(g_run_mutex);
            FdRedirectGuard fdguard;

            bool dryrun = false;
            idx->opts.reset(new Runopts(static_cast<int>(argv_ptrs.size() - 1),
                                        argv_ptrs.data(), dryrun));
            idx->opts->is_library_mode = true;

            idx->index.reset(new Index(*idx->opts));
            /* fresh kvdb; smr_index_free removes it */
            clear_kvdb(idx->opts->kvdbdir);
            idx->kvdb.reset(new KeyValueDatabase(idx->opts->kvdbdir.string()));

            /* Library-mode Refstats needs no reads: it uses the raw
             * reference size and zeroes the minimal SW score. */
            Readstats dummy_rs(0, 0, 0, 0, *idx->kvdb, *idx->opts);
            idx->refstats.reset(new Refstats(*idx->opts, dummy_rs));
        }
        /* Mutex + fdguard released. index.load/references.load log only
         * through the thread-local LogRouteGuard, so they can run
         * concurrently across handles. */

        /* MVP: require single-part, single-ref index. Multi-part support
         * means either loading all parts into memory (large) or re-loading
         * per batch (defeats the point of the handle). Reject for now. */
        for (size_t ri = 0; ri < idx->opts->indexfiles.size(); ri++) {
            if (idx->refstats->num_index_parts[ri] > 1) {
                set_error(ctx, SMR_ERR_NOT_IMPLEMENTED,
                          "multi-part index (%u parts for ref %zu) not yet supported by smr_index_load",
                          idx->refstats->num_index_parts[ri], ri);
                smr_index_free(idx);
                return nullptr;
            }
        }

        /* Load ref 0 / part 0 into memory — expensive, done once per handle. */
        ctx_log(ctx, SMR_LOG_INFO, "smr_index_load: loading references and index into memory");
        idx->index->load(0, 0, idx->opts->indexfiles, *idx->refstats);
        idx->references.load(0, 0, *idx->opts, *idx->refstats);
        idx->index_loaded = true;
        ctx_log(ctx, SMR_LOG_INFO, "smr_index_load: handle ready (%zu ref(s), single-part)",
                idx->opts->indexfiles.size());

    } catch (const smr_exit_requested &) {
        set_error(ctx, SMR_ERR_INVALID_CONFIG, "unexpected --help/--version in library context");
        smr_index_free(idx);
        return nullptr;
    } catch (const std::exception &e) {
        set_error(ctx, SMR_ERR_INDEX, "%s", e.what());
        smr_index_free(idx);
        return nullptr;
    } catch (...) {
        set_error(ctx, SMR_ERR_INDEX, "unknown exception in smr_index_load");
        smr_index_free(idx);
        return nullptr;
    }

    set_error(ctx, SMR_OK, "");
    return idx;
}

int smr_run_seqs_with_index(smr_index_t *idx,
                            const smr_seq_t *seqs, int32_t num_seqs,
                            smr_output_t **out,
                            smr_stats_t *stats) {
    if (!idx) return SMR_ERR_INVALID_CONFIG;
    smr_context_t *ctx = idx->ctx;
    if (!idx->index_loaded) {
        set_error(ctx, SMR_ERR_INDEX, "handle index not loaded");
        return SMR_ERR_INDEX;
    }

    bool has_qual = false;
    int vrc = validate_seqs(ctx, seqs, num_seqs, ctx->config.paired != 0, &has_qual);
    if (vrc != SMR_OK) return vrc;
    /* FASTA and FASTQ batches are both accepted. The placeholder_reads file
     * at load time satisfies Runopts validation but is never read for
     * sequence data; Readfeed MEMORY mode derives format from the presence
     * of quality strings in the in-memory batch. */

    std::lock_guard<std::mutex> run_lock(g_run_mutex);
    auto wall_start = std::chrono::high_resolution_clock::now();

    LogRouteGuard loguard(ctx->config.log_callback, ctx->config.log_user_data);
    FdRedirectGuard fdguard;

    try {
        /* The kvdb is shared across batches. Library-mode alignment stores
         * every read of the batch under its slot key, overwriting whatever an
         * earlier batch left, so no stale result is read back. */
        Runopts &opts = *idx->opts;

        ctx_log(ctx, SMR_LOG_INFO, "pipeline starting (in-memory, with loaded index): %d seqs",
                num_seqs);

        std::vector<std::string> ids, sequences, quals;
        ids.reserve(num_seqs);
        sequences.reserve(num_seqs);
        if (has_qual) quals.reserve(num_seqs);
        for (int32_t i = 0; i < num_seqs; i++) {
            ids.emplace_back(seqs[i].id);
            sequences.emplace_back(seqs[i].sequence);
            if (has_qual) quals.emplace_back(seqs[i].quality);
        }

        auto basedir = opts.readb_dir;
        Readfeed readfeed(std::move(ids), std::move(sequences), std::move(quals),
                          opts.num_proc_thread, basedir, opts.is_paired);
        Readstats readstats(readfeed.num_reads_tot, readfeed.length_all,
                            readfeed.min_read_len, readfeed.max_read_len, *idx->kvdb, opts);

        align_loaded(readfeed, readstats, *idx->index, idx->references, *idx->refstats,
                     *idx->kvdb, opts);

        /* Library path: writeSummary/writeReports skipped — library callers
         * consume smr_output_t directly from kvdb, not the report files. */

        ctx_log(ctx, SMR_LOG_INFO, "pipeline complete");

        if (stats) {
            memset(stats, 0, sizeof(*stats));
            stats->total_reads = readstats.all_reads_count;
            stats->total_aligned = readstats.num_aligned.load();
            stats->total_id_cov_pass = readstats.n_yid_ycov.load();
            stats->total_denovo = readstats.num_denovo.load();
            stats->min_read_len = readstats.min_read_len;
            stats->max_read_len = readstats.max_read_len;
            auto wall_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = wall_end - wall_start;
            stats->wall_time_sec = elapsed.count();
        }

        if (out) {
            auto *o = static_cast<smr_output_t *>(calloc(1, sizeof(smr_output_t)));
            if (!o) {
                set_error(ctx, SMR_ERR_ALLOC, "failed to allocate smr_output_t");
                return SMR_ERR_ALLOC;
            }

            o->num_reads = readstats.all_reads_count;

            if (!populate_per_read_output(o, readfeed, readstats, *idx->kvdb, opts,
                                          &idx->references, idx->refstats.get())) {
                smr_output_free(o);
                set_error(ctx, SMR_ERR_ALLOC, "failed to allocate per-read output arrays");
                return SMR_ERR_ALLOC;
            }

            *out = o;
        }

        set_error(ctx, SMR_OK, "");
        return SMR_OK;

    } catch (const smr_exit_requested &) {
        set_error(ctx, SMR_ERR_INVALID_CONFIG, "unexpected --help/--version in library context");
        return SMR_ERR_INVALID_CONFIG;
    } catch (const std::exception &e) {
        set_error(ctx, SMR_ERR_ALIGN, "%s", e.what());
        return SMR_ERR_ALIGN;
    } catch (...) {
        set_error(ctx, SMR_ERR_ALIGN, "unknown exception in pipeline");
        return SMR_ERR_ALIGN;
    }
}

void smr_index_free(smr_index_t *idx) {
    if (!idx) return; /* NULL is always safe -- contract */
    /* The Index destructor frees the loaded index, also after a failed
     * smr_index_load; the unique_ptr destructors close the kvdb and release
     * Refstats / Runopts. */
    if (idx->kvdb) {
        idx->kvdb.reset();
        std::error_code ec;
        std::filesystem::remove_all(idx->opts->kvdbdir, ec);
    }
    idx->refstats.reset();
    idx->index.reset();
    idx->opts.reset();
    if (idx->workdir_is_temp && !idx->workdir.empty()) {
        std::error_code ec;
        std::filesystem::remove_all(idx->workdir, ec);
    }
    delete idx;
}

int smr_run_seqs(smr_context_t *ctx,
                 const char **ref_paths, int32_t num_refs,
                 const smr_seq_t *seqs, int32_t num_seqs,
                 smr_output_t **out,
                 smr_stats_t *stats) {
    if (!ctx) return SMR_ERR_INVALID_CONFIG;
    smr_index_t *idx = smr_index_load(ctx, ref_paths, num_refs);
    /* Propagate the specific error code set by smr_index_load via ctx. Safe
     * because smr_index_load's NULL return paths all call set_error(ctx, code)
     * and no other code runs between that and here to overwrite it. */
    if (!idx) return smr_last_error_code(ctx);
    int rc = smr_run_seqs_with_index(idx, seqs, num_seqs, out, stats);
    smr_index_free(idx);
    return rc;
}

void smr_output_free(smr_output_t *out) {
    if (!out) return; /* NULL is always safe -- contract */
    /* free library-owned string arrays */
    if (out->read_ids) {
        for (uint64_t i = 0; i < out->num_reads; i++)
            free(const_cast<char*>(out->read_ids[i]));
        free(const_cast<char**>(out->read_ids));
    }
    if (out->cigar) {
        for (uint64_t i = 0; i < out->num_reads; i++)
            free(const_cast<char*>(out->cigar[i]));
        free(const_cast<char**>(out->cigar));
    }
    if (out->ref_name) {
        for (uint64_t i = 0; i < out->num_reads; i++)
            free(const_cast<char*>(out->ref_name[i]));
        free(const_cast<char**>(out->ref_name));
    }
    free(out->aligned);
    free(out->ref_index);
    free(out->e_value);
    free(out->identity);
    free(out->coverage);
    free(out->ref_start);
    free(out->ref_end);
    free(out->strand);
    free(out->score);
    free(out->edit_distance);
    free(out);
}

/* --- Version --- */

const char *smr_version(void) {
    return SMR_STRINGIFY(SORTMERNA_MAJOR) "."
           SMR_STRINGIFY(SORTMERNA_MINOR) "."
           SMR_STRINGIFY(SORTMERNA_PATCH);
}
