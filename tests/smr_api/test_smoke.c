/*
 * test_smoke.c -- main test runner for smr_api
 * Compiled as C to verify smr_api.h is C-linkable.
 *
 * Single-TU design: all test files are #include'd here.
 */

#include "test_harness.h"
#include "smr_api.h"
#include "version.h"
#include <stddef.h>  /* offsetof */
#include <stdlib.h>  /* mkdtemp, system */
#include <unistd.h>  /* dup, dup2, rmdir */
#include <pthread.h> /* concurrent test */
#include <sys/wait.h> /* waitpid */
#include <math.h>    /* fabs */

/*
 * Compare a produced output file against a golden reference file.
 * For SAM files, skip @-header lines (contain workdir paths that vary).
 * Returns 1 if files match, 0 otherwise.
 */
static int files_match(const char *produced, const char *golden, int skip_sam_headers) {
    FILE *fp = fopen(produced, "r");
    FILE *fg = fopen(golden, "r");
    if (!fp || !fg) {
        if (fp) fclose(fp);
        if (fg) fclose(fg);
        return 0;
    }
    char lp[4096], lg[4096];
    int match = 1;
    while (1) {
        char *rp = fgets(lp, sizeof(lp), fp);
        char *rg = fgets(lg, sizeof(lg), fg);
        /* skip SAM header lines in both files */
        if (skip_sam_headers) {
            while (rp && lp[0] == '@') rp = fgets(lp, sizeof(lp), fp);
            while (rg && lg[0] == '@') rg = fgets(lg, sizeof(lg), fg);
        }
        if (!rp && !rg) break; /* both EOF */
        if (!rp || !rg) { match = 0; break; } /* one EOF early */
        if (strcmp(lp, lg) != 0) { match = 0; break; }
    }
    fclose(fp);
    fclose(fg);
    return match;
}

/* Temporary files and directories go under TMPDIR, else /tmp. */
static const char *tmp_base(void) {
    const char *t = getenv("TMPDIR");
    return (t && t[0]) ? t : "/tmp";
}

static char _tmpdir_buf[512];
static char *make_tmpdir(void) {
    snprintf(_tmpdir_buf, sizeof(_tmpdir_buf), "%s/smr_test_XXXXXX", tmp_base());
    return mkdtemp(_tmpdir_buf);
}

static void rm_rf(const char *dir) {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "rm -rf '%s'", dir);
    system(cmd);
}

/* ---- Smoke tests ---- */

TEST(test_config_init_sets_struct_size) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT((int)cfg.struct_size, (int)sizeof(smr_config_t));
}

TEST(test_config_struct_size_is_first_field) {
    ASSERT_EQ_INT((int)offsetof(smr_config_t, struct_size), 0);
}

TEST(test_config_init_zeroes_pointers) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_NULL(cfg.workdir);
    ASSERT_NULL(cfg.log_callback);
    ASSERT_NULL(cfg.log_user_data);
}

TEST(test_version_matches_version_h) {
    char expected[32];
    snprintf(expected, sizeof(expected), "%d.%d.%d",
             SORTMERNA_MAJOR, SORTMERNA_MINOR, SORTMERNA_PATCH);
    const char *v = smr_version();
    ASSERT_NOT_NULL(v);
    ASSERT_STR_EQ(v, expected);
}

TEST(test_strerror_returns_string) {
    const char *s = smr_strerror(SMR_OK);
    ASSERT_NOT_NULL(s);
}

TEST(test_output_free_null_safe) {
    smr_output_free(NULL);
    ASSERT_TRUE(1);
}

TEST(test_ctx_destroy_null_safe) {
    smr_ctx_destroy(NULL);
    ASSERT_TRUE(1);
}

/* ---- Config defaults ---- */

TEST(test_config_default_threads) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.num_threads, 2);
}

TEST(test_config_default_num_alignments) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.num_alignments, 1);
}

TEST(test_config_default_match) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.match, 2);
}

TEST(test_config_default_mismatch) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.mismatch, -3);
}

TEST(test_config_default_gap_open) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.gap_open, 5);
}

TEST(test_config_default_gap_ext) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.gap_ext, 2);
}

TEST(test_config_default_score_N) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.score_N, -3);
}

TEST(test_config_default_evalue) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    /* unset: the sortmerna default applies; -1.0 is exactly representable */
    ASSERT_TRUE(cfg.evalue == SMR_EVALUE_DEFAULT);
}

TEST(test_config_default_seed_win_len) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.seed_win_len, 18);
}

TEST(test_config_default_best) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.best, 1);
}

TEST(test_config_default_booleans_off) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    ASSERT_EQ_INT(cfg.paired, 0);
    ASSERT_EQ_INT(cfg.forward_only, 0);
    ASSERT_EQ_INT(cfg.reverse_only, 0);
    ASSERT_EQ_INT(cfg.full_search, 0);
}

TEST(test_config_boolean_type_is_int32) {
    ASSERT_EQ_SZ(sizeof(((smr_config_t *)0)->best), sizeof(int32_t));
    ASSERT_EQ_SZ(sizeof(((smr_config_t *)0)->paired), sizeof(int32_t));
    ASSERT_EQ_SZ(sizeof(((smr_config_t *)0)->forward_only), sizeof(int32_t));
}

TEST(test_config_explicit_width_types) {
    ASSERT_EQ_SZ(sizeof(((smr_config_t *)0)->num_threads), sizeof(int32_t));
    ASSERT_EQ_SZ(sizeof(((smr_config_t *)0)->seed_win_len), sizeof(uint32_t));
    ASSERT_EQ_SZ(sizeof(((smr_config_t *)0)->num_alignments), sizeof(uint32_t));
}

/* ---- Context lifecycle and error reporting ---- */

static void test_log_cb(int level, const char *msg, void *user_data) {
    (void)level; (void)msg;
    int *counter = (int *)user_data;
    if (counter) (*counter)++;
}

TEST(test_ctx_create_returns_non_null) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    smr_ctx_destroy(ctx);
}

TEST(test_ctx_create_null_config_returns_null) {
    smr_context_t *ctx = smr_ctx_create(NULL);
    ASSERT_NULL(ctx);
}

TEST(test_ctx_create_bad_struct_size_returns_null) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.struct_size = 0;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NULL(ctx);
}

TEST(test_ctx_create_multiple_independent) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *a = smr_ctx_create(&cfg);
    smr_context_t *b = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(a);
    ASSERT_NOT_NULL(b);
    ASSERT_TRUE(a != b);
    smr_ctx_destroy(a);
    smr_ctx_destroy(b);
}

TEST(test_ctx_last_error_empty_initially) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *err = smr_last_error(ctx);
    ASSERT_NOT_NULL(err);
    ASSERT_EQ_INT((int)err[0], 0); /* empty string */
    smr_ctx_destroy(ctx);
}

TEST(test_strerror_success_msg) {
    const char *s = smr_strerror(SMR_OK);
    ASSERT_STR_EQ(s, "Success");
}

TEST(test_strerror_invalid_config_msg) {
    const char *s = smr_strerror(SMR_ERR_INVALID_CONFIG);
    ASSERT_NOT_NULL(s);
    ASSERT_TRUE(s[0] != '\0');
}

TEST(test_strerror_unknown_code) {
    const char *s = smr_strerror(99999);
    ASSERT_NOT_NULL(s);
    ASSERT_STR_EQ(s, "Unknown error");
}

TEST(test_log_callback_receives_messages) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    int count = 0;
    cfg.log_callback = test_log_cb;
    cfg.log_user_data = &count;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    ASSERT_TRUE(count > 0);
    smr_ctx_destroy(ctx);
}

TEST(test_last_error_set_after_smr_run) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    int rc = smr_run(ctx, NULL, 0, NULL, 0, NULL, NULL);
    ASSERT_TRUE(rc < 0);
    ASSERT_TRUE(smr_last_error_code(ctx) < 0);
    const char *err = smr_last_error(ctx);
    ASSERT_NOT_NULL(err);
    ASSERT_TRUE(err[0] != '\0');
    smr_ctx_destroy(ctx);
}

TEST(test_last_error_code_initially_zero) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    ASSERT_EQ_INT(smr_last_error_code(ctx), SMR_OK);
    smr_ctx_destroy(ctx);
}

TEST(test_last_error_code_null_ctx) {
    ASSERT_EQ_INT(smr_last_error_code(NULL), SMR_ERR_INVALID_CONFIG);
}

/* ---- Independent contexts ---- */

TEST(test_two_contexts_create_destroy) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *a = smr_ctx_create(&cfg);
    smr_context_t *b = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(a);
    ASSERT_NOT_NULL(b);
    /* destroy in reverse order to test independence */
    smr_ctx_destroy(b);
    smr_ctx_destroy(a);
    /* re-create after destroy to test reuse */
    smr_context_t *c = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(c);
    smr_ctx_destroy(c);
}

/* ---- Error handling ---- */

TEST(test_run_null_refs_returns_error) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, NULL, 0, reads, 1, &out, &stats);
    ASSERT_TRUE(rc < 0);
    smr_ctx_destroy(ctx);
}

TEST(test_run_null_reads_returns_error) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, NULL, 0, &out, &stats);
    ASSERT_TRUE(rc < 0);
    smr_ctx_destroy(ctx);
}

TEST(test_run_nonexistent_ref_returns_error) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { "/nonexistent/ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_ERR_IO);
    smr_ctx_destroy(ctx);
}

TEST(test_run_nonexistent_reads_returns_error) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { "/nonexistent/reads.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_ERR_IO);
    smr_ctx_destroy(ctx);
}

TEST(test_run_empty_ref_returns_error) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/empty_file.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_TRUE(rc < 0);
    smr_ctx_destroy(ctx);
}

TEST(test_last_error_descriptive_after_bad_input) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { "/nonexistent/ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    const char *err = smr_last_error(ctx);
    ASSERT_NOT_NULL(err);
    ASSERT_TRUE(err[0] != '\0');
    smr_ctx_destroy(ctx);
}

TEST(test_run_bad_input_does_not_crash) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    /* multiple bad calls in sequence — none should crash */
    smr_run(ctx, NULL, 0, NULL, 0, NULL, NULL);
    const char *refs[] = { "/bad" };
    const char *reads[] = { "/bad" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    smr_run(ctx, refs, 1, NULL, 0, &out, &stats);
    /* if we got here, nothing crashed */
    ASSERT_TRUE(1);
    smr_ctx_destroy(ctx);
}

/* ---- smr_run computation ---- */

TEST(test_run_tiny_aligned_count) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 1);
    ASSERT_EQ_U64(out->num_aligned, 1);
    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_tiny_stats) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    memset(&stats, 0, sizeof(stats));
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_U64(stats.total_reads, 1);
    ASSERT_EQ_U64(stats.total_aligned, 1);
    ASSERT_TRUE(stats.min_read_len > 0);
    ASSERT_TRUE(stats.max_read_len > 0);
    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_small_aligned_count) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 6);
    ASSERT_EQ_U64(out->num_aligned, 4);
    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_output_free_after_run) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    smr_output_free(out);
    /* if we got here, no crash */
    ASSERT_TRUE(1);
    smr_ctx_destroy(ctx);
}

TEST(test_run_multiple_sequential) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out1 = NULL;
    smr_output_t *out2 = NULL;
    smr_stats_t stats;
    int rc1 = smr_run(ctx, refs, 1, reads, 1, &out1, &stats);
    int rc2 = smr_run(ctx, refs, 1, reads, 1, &out2, &stats);
    ASSERT_EQ_INT(rc1, SMR_OK);
    ASSERT_EQ_INT(rc2, SMR_OK);
    ASSERT_EQ_U64(out1->num_aligned, out2->num_aligned);
    smr_output_free(out1);
    smr_output_free(out2);
    smr_ctx_destroy(ctx);
}

/* ---- I/O isolation ---- */

TEST(test_no_stdout_during_run) {
    /* redirect stdout to /dev/null, run pipeline, check nothing leaked */
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;

    /* capture stdout to a temp file */
    fflush(stdout);
    int old_stdout = dup(1);
    FILE *tmp = tmpfile();
    int tmp_fd = fileno(tmp);
    dup2(tmp_fd, 1);

    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);

    /* restore stdout */
    fflush(stdout);
    dup2(old_stdout, 1);
    close(old_stdout);

    /* check run succeeded and captured output is empty */
    ASSERT_EQ_INT(rc, SMR_OK);

    fseek(tmp, 0, SEEK_END);
    long captured_size = ftell(tmp);
    fclose(tmp);

    smr_output_free(out);
    smr_ctx_destroy(ctx);

    ASSERT_EQ_INT((int)captured_size, 0);
}

static void run_log_cb(int level, const char *msg, void *user_data) {
    (void)level; (void)msg;
    int *counter = (int *)user_data;
    if (counter) (*counter)++;
}

TEST(test_log_callback_fires_during_run) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    int count = 0;
    cfg.log_callback = run_log_cb;
    cfg.log_user_data = &count;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    /* at minimum: "context created", "pipeline starting", "pipeline complete" */
    ASSERT_TRUE(count >= 3);
    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_null_log_callback_silent) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.log_callback = NULL;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

/* ---- Golden file comparison ----
 *
 * tests/golden/{tiny,small} hold the reports of the sortmerna CLI run with the
 * same arguments build_argv() passes (SAM @-header lines are not compared).
 *
 * test_read.fasta and set7_arc_bac_16S_database_match.fasta have no trailing
 * newline, and the CLI's read parser currently drops the last line of such a
 * file: AB271211 is read as 1440 of its 1487 nt, and random2 as an empty
 * sequence. The goldens, and the smr_run expectations below, reflect that. */

TEST(test_golden_tiny_blast) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);

    char produced[512];
    snprintf(produced, sizeof(produced), "%s/aligned.blast", wdir);
    ASSERT_TRUE(files_match(produced, SMR_GOLDEN_DIR "/tiny/aligned.blast", 0));

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

TEST(test_golden_tiny_fasta) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);

    char produced[512];
    snprintf(produced, sizeof(produced), "%s/aligned.fa", wdir);
    ASSERT_TRUE(files_match(produced, SMR_GOLDEN_DIR "/tiny/aligned.fa", 0));

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

TEST(test_golden_tiny_sam) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);

    char produced[512];
    snprintf(produced, sizeof(produced), "%s/aligned.sam", wdir);
    ASSERT_TRUE(files_match(produced, SMR_GOLDEN_DIR "/tiny/aligned.sam", 1));

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

TEST(test_golden_small_blast) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);

    char produced[512];
    snprintf(produced, sizeof(produced), "%s/aligned.blast", wdir);
    ASSERT_TRUE(files_match(produced, SMR_GOLDEN_DIR "/small/aligned.blast", 0));

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

TEST(test_golden_small_fasta) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);

    char produced[512];
    snprintf(produced, sizeof(produced), "%s/aligned.fa", wdir);
    ASSERT_TRUE(files_match(produced, SMR_GOLDEN_DIR "/small/aligned.fa", 0));

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

TEST(test_golden_small_sam) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);

    char produced[512];
    snprintf(produced, sizeof(produced), "%s/aligned.sam", wdir);
    ASSERT_TRUE(files_match(produced, SMR_GOLDEN_DIR "/small/aligned.sam", 1));

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

TEST(test_golden_small_other) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    smr_run(ctx, refs, 1, reads, 1, &out, &stats);

    char produced[512];
    snprintf(produced, sizeof(produced), "%s/other.fa", wdir);
    ASSERT_TRUE(files_match(produced, SMR_GOLDEN_DIR "/small/other.fa", 0));

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

/* ---- Concurrent smr_run ---- */

struct thread_arg {
    int rc;
    uint64_t num_aligned;
};

static void *concurrent_worker(void *arg) {
    struct thread_arg *ta = (struct thread_arg *)arg;
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    ta->rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ta->num_aligned = out ? out->num_aligned : 0;
    smr_output_free(out);
    smr_ctx_destroy(ctx);
    return NULL;
}

TEST(test_concurrent_smr_run) {
    pthread_t t1, t2;
    struct thread_arg a1 = {0, 0}, a2 = {0, 0};
    pthread_create(&t1, NULL, concurrent_worker, &a1);
    pthread_create(&t2, NULL, concurrent_worker, &a2);
    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
    ASSERT_EQ_INT(a1.rc, SMR_OK);
    ASSERT_EQ_INT(a2.rc, SMR_OK);
    ASSERT_EQ_U64(a1.num_aligned, 1);
    ASSERT_EQ_U64(a2.num_aligned, 1);
}

/* ---- Multi-threaded alignment ---- */

TEST(test_run_tiny_multithreaded) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 2;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 1);
    ASSERT_EQ_U64(out->num_aligned, 1);
    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

/* ---- Per-read output ---- */

TEST(test_run_tiny_per_read_output) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 1);
    ASSERT_EQ_U64(out->num_aligned, 1);

    /* per-read arrays should be populated */
    ASSERT_NOT_NULL(out->read_ids);
    ASSERT_NOT_NULL(out->aligned);
    ASSERT_NOT_NULL(out->e_value);
    ASSERT_NOT_NULL(out->identity);
    ASSERT_NOT_NULL(out->coverage);
    ASSERT_NOT_NULL(out->ref_start);
    ASSERT_NOT_NULL(out->ref_end);
    ASSERT_NOT_NULL(out->cigar);

    /* golden tiny/aligned.blast (read truncated to 1440 nt, see above) */
    ASSERT_STR_EQ(out->read_ids[0], "AB271211");
    ASSERT_EQ_INT(out->aligned[0], 1);
    ASSERT_TRUE(out->ref_start[0] == 1);
    ASSERT_TRUE(out->ref_end[0] == 1399);
    ASSERT_DOUBLE_NEAR(out->identity[0], 92.3, 0.05);
    ASSERT_DOUBLE_NEAR(out->coverage[0], 97.2, 0.05);
    ASSERT_NOT_NULL(out->cigar[0]);
    ASSERT_STR_EQ(out->cigar[0],
        "41S16I57M2I12M2D4M2I29M1D11M2I3M2D11M1I7M1D13M5D4M3D9M2D3M7D1213M");

    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_small_per_read_output) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 6);
    ASSERT_EQ_U64(out->num_aligned, 4);

    ASSERT_NOT_NULL(out->read_ids);
    ASSERT_NOT_NULL(out->aligned);

    /* first 4 reads are aligned, last 2 are not */
    ASSERT_STR_EQ(out->read_ids[0], "BD.ERD505_1");
    ASSERT_EQ_INT(out->aligned[0], 1);
    ASSERT_STR_EQ(out->read_ids[1], "BD.NBS1076_0");
    ASSERT_EQ_INT(out->aligned[1], 1);
    ASSERT_STR_EQ(out->read_ids[2], "LD.Glosor1_17");
    ASSERT_EQ_INT(out->aligned[2], 1);
    ASSERT_STR_EQ(out->read_ids[3], "BD.ERD510_20");
    ASSERT_EQ_INT(out->aligned[3], 1);
    ASSERT_STR_EQ(out->read_ids[4], "random1");
    ASSERT_EQ_INT(out->aligned[4], 0);
    ASSERT_STR_EQ(out->read_ids[5], "random2");
    ASSERT_EQ_INT(out->aligned[5], 0);

    /* unaligned reads have ref_index == -1 */
    ASSERT_EQ_INT(out->ref_index[4], -1);
    ASSERT_EQ_INT(out->ref_index[5], -1);
    ASSERT_NULL(out->cigar[4]);
    ASSERT_NULL(out->cigar[5]);

    /* spot-check aligned read identity/coverage (golden small/aligned.blast) */
    ASSERT_DOUBLE_NEAR(out->identity[0], 90.7, 0.05);
    ASSERT_DOUBLE_NEAR(out->coverage[0], 54.0, 0.05);
    ASSERT_DOUBLE_NEAR(out->coverage[1], 100.0, 0.05);

    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

/* ---- In-memory input (smr_run_seqs) ---- */

TEST(test_run_seqs_tiny) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };

    /* AB271211 sequence (same as test_read.fasta) — first 100 chars for a shorter test */
    smr_seq_t seqs[1];
    seqs[0].id = "AB271211";
    seqs[0].sequence =
        "TCCAACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTAG"
        "AGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGG"
        "GAATCTTCGGATTCTAGTGGCGGACGGGTGAGTAACGCGTAAGAATCTAACTTCAGGACG"
        "GGGACAACAGTGGGAAACGACTGCTAATACCCGATGTGCCGCGAGGTGAAACCTAATTGG";
    seqs[0].quality = NULL;

    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run_seqs(ctx, refs, 1, seqs, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 1);
    ASSERT_NOT_NULL(out->read_ids);
    ASSERT_STR_EQ(out->read_ids[0], "AB271211");

    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_seqs_null_error) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run_seqs(ctx, refs, 1, NULL, 0, &out, &stats);
    ASSERT_TRUE(rc < 0);
    smr_ctx_destroy(ctx);
}

TEST(test_run_seqs_empty_seq_error) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_seq_t seqs[1];
    seqs[0].id = "test";
    seqs[0].sequence = "";
    seqs[0].quality = NULL;
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run_seqs(ctx, refs, 1, seqs, 1, &out, &stats);
    ASSERT_TRUE(rc < 0);
    smr_ctx_destroy(ctx);
}

TEST(test_run_seqs_quality_length_mismatch) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_seq_t seqs[1];
    seqs[0].id = "test";
    seqs[0].sequence = "ACGT";
    seqs[0].quality = "II";  /* too short */
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run_seqs(ctx, refs, 1, seqs, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_ERR_INVALID_CONFIG);
    smr_ctx_destroy(ctx);
}

TEST(test_run_seqs_paired_odd_count_error) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.paired = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_seq_t seqs[3];
    seqs[0].id = "r1"; seqs[0].sequence = "ACGT"; seqs[0].quality = NULL;
    seqs[1].id = "r2"; seqs[1].sequence = "TGCA"; seqs[1].quality = NULL;
    seqs[2].id = "r3"; seqs[2].sequence = "AAAA"; seqs[2].quality = NULL;
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run_seqs(ctx, refs, 1, seqs, 3, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_ERR_INVALID_CONFIG);
    smr_ctx_destroy(ctx);
}

/* ---- ref_name output ---- */

TEST(test_run_tiny_ref_name) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_NOT_NULL(out->ref_name);
    ASSERT_NOT_NULL(out->ref_name[0]);
    /* golden: AB271211 aligned to Unc49508 */
    ASSERT_STR_EQ(out->ref_name[0], "Unc49508");
    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_small_ref_name) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_NOT_NULL(out->ref_name);
    /* golden: BD.ERD505_1 → EU602318, unaligned reads have NULL ref_name */
    ASSERT_NOT_NULL(out->ref_name[0]);
    ASSERT_STR_EQ(out->ref_name[0], "EU602318");
    ASSERT_NULL(out->ref_name[4]); /* random1 — unaligned */
    ASSERT_NULL(out->ref_name[5]); /* random2 — unaligned */
    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_seqs_paired) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.paired = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };

    /* 2 interleaved pairs (4 sequences total) from the small dataset */
    smr_seq_t seqs[4];
    seqs[0].id = "BD.ERD505_1";
    seqs[0].sequence = "AACGTAGGTGGCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGCAGGCGGAAAAGCAAGTTGGACGTGAAATCTATGGGCTCAACCCATAGCGTG";
    seqs[0].quality = NULL;
    seqs[1].id = "BD.NBS1076_0";
    seqs[1].sequence = "TACGGAGGGTGCAAGCGTTAATCCGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCACGGGCTCAACCTGG";
    seqs[1].quality = NULL;
    seqs[2].id = "LD.Glosor1_17";
    seqs[2].sequence = "TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTG";
    seqs[2].quality = NULL;
    seqs[3].id = "BD.ERD510_20";
    seqs[3].sequence = "TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTG";
    seqs[3].quality = NULL;

    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run_seqs(ctx, refs, 1, seqs, 4, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 4);
    ASSERT_NOT_NULL(out->read_ids);
    /* verify interleaved input order is preserved */
    ASSERT_STR_EQ(out->read_ids[0], "BD.ERD505_1");
    ASSERT_STR_EQ(out->read_ids[1], "BD.NBS1076_0");
    ASSERT_STR_EQ(out->read_ids[2], "LD.Glosor1_17");
    ASSERT_STR_EQ(out->read_ids[3], "BD.ERD510_20");

    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

/* ---- strand / score / edit_distance output ---- */

TEST(test_run_tiny_strand_score_edit) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[]  = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 1);
    ASSERT_EQ_U64(out->num_aligned, 1);

    /* arrays must be non-NULL */
    ASSERT_NOT_NULL(out->strand);
    ASSERT_NOT_NULL(out->score);
    ASSERT_NOT_NULL(out->edit_distance);

    /* golden tiny/aligned.sam: FLAG=0 -> forward, AS:i:2336, NM:i:110 */
    ASSERT_EQ_INT(out->strand[0], 1);
    ASSERT_EQ_INT(out->score[0], 2336);
    ASSERT_EQ_INT(out->edit_distance[0], 110);

    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_small_strand_score_edit) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[]  = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 6);
    ASSERT_EQ_U64(out->num_aligned, 4);

    ASSERT_NOT_NULL(out->strand);
    ASSERT_NOT_NULL(out->score);
    ASSERT_NOT_NULL(out->edit_distance);

    /* aligned reads: strand=1 (forward); AS:i and NM:i of golden small/aligned.sam */
    ASSERT_EQ_INT(out->strand[0], 1);         /* BD.ERD505_1 */
    ASSERT_EQ_INT(out->score[0], 83);
    ASSERT_EQ_INT(out->edit_distance[0], 5);

    ASSERT_EQ_INT(out->strand[1], 1);         /* BD.NBS1076_0 */
    ASSERT_EQ_INT(out->score[1], 132);
    ASSERT_EQ_INT(out->edit_distance[1], 11);

    ASSERT_EQ_INT(out->strand[2], 1);         /* LD.Glosor1_17 */
    ASSERT_EQ_INT(out->score[2], 154);
    ASSERT_EQ_INT(out->edit_distance[2], 9);

    ASSERT_EQ_INT(out->strand[3], 1);         /* BD.ERD510_20 */
    ASSERT_EQ_INT(out->score[3], 154);
    ASSERT_EQ_INT(out->edit_distance[3], 9);

    /* unaligned reads: sentinel -1 for all three */
    ASSERT_EQ_INT(out->strand[4], -1);        /* random1 */
    ASSERT_EQ_INT(out->score[4], -1);
    ASSERT_EQ_INT(out->edit_distance[4], -1);

    ASSERT_EQ_INT(out->strand[5], -1);        /* random2 */
    ASSERT_EQ_INT(out->score[5], -1);
    ASSERT_EQ_INT(out->edit_distance[5], -1);

    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_reverse_strand) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };

    /* reverse-complement of AB271211 — should align on reverse strand */
    smr_seq_t seqs[1];
    seqs[0].id = "AB271211_rc";
    seqs[0].sequence =
        "CCCCAGTCACTAGCCCTGCCTTAGGCATCCCCCTCCTTGCGGTTGAGGTAATGACTTCGGG"
        "CGTGACCAGCTTCCATGGTGTGACGGGCGGTGTGTACAAGGCCCGGGAACGAATTCACCG"
        "CCGTATGCTGACCGGCGATTACTAGCGATTCCTCCTTCATGCAGGCGAGTTGCAGCCTGC"
        "AATCTGAACTGAGGCCGGGTTTGCTGGGATTCGCTGGCTCTCGCAAGTTCGCTGCCCTTT"
        "GTCCCGACCATTGTAGTACGTGTGTCGCCCAAGACGTAAGGGGCATGCTGACTTGACGTC"
        "ATCCCCACCTTCCTCCGGTTTGTCACCGGCAGTCTCCTTAGAGTCCCCAACTTAATGCTGG"
        "CAACTAAGAACGAGGGTTGCGCTCGTTGCGGGACTTAACCCAACATCTCACGACACGAGC"
        "TGACGACAGCCATGCACCACCTGTGTTCGCGCTCCCGAAGGCACCCCCAGCTTTCACCAGG"
        "GTTCGCGACATGTCAAGTCTTGGTAAGGTTCTTCGCGTTGCATCGAATTAAACCACATAC"
        "TCCACCGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTCACACTTGCGTGCGTACTCC"
        "CCAGGCGGGATACTTAACGCGTTAGCTTCGGCACGGCTCGGGTCGATACAAGCCACGCCTA"
        "GTATCCATCGTTTACGGCTAGGACTACAGGGGTATCTAATCCCTTTCGCTCCCCTAGCTTT"
        "CGTCCCTGAGTGTCAGATACAGCCCAGTAGCACGCTTTCGCCACCGATGTTCTTCCCAATC"
        "TCTACGCATTTCACCGCTACACTGGGAATTCCTGCTACCCCTACTGCTCTCTAGTCTGCCA"
        "GTTTCCACCGCCTTTAGGTCGTTAAGCAACCTGATTTGACGGCAGACTTGGCTGACCACCT"
        "GCGGACGCTTTACGCCCAATAATTCCGGATAACGCTTGCCTCCCCCGTATTACCGCGGCTG"
        "CTGGCACGGAGTTAGCCGAGGCTGATTCCTCAAGTACCGTCAGAACTTCTTCCTTGAGAAA"
        "AGAGGTTTACAATCCAAAGACCTTCCTCCCTCACGCGGCGTTGCTCCGTCAGGCTTTCGCC"
        "CATTGCGGAAAATTCCCCACTGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAG"
        "TGTGGCTGCTCATCCTCTCAGACCAGCTACTGATCGTCGCCTTGGTAGGCTCTTACCCCAC"
        "CAACTAGCTAATCAGACGCAAGCTCCTCTTCAGGCCAATTAGGTTTCACCTCGCGGCACAT"
        "CGGGTATTAGCAGTCGTTTCCCACTGTTGTCCCCGTCCTGAAGTTAGATTCTTACGCGTTA"
        "CTCACCCGTCCGCCACTAGAATCCGAAGATTCCCGTTCGACTTGCATGTGTTAGGCACGCC"
        "GCCAGCGTTCATCCTGAGCCAGGATCAAACTCTAATCACTAGTGCGGCCGCCTGCAGGTCG"
        "ACCATATGGGAGAGCTCCCAACGCGTTGGA";
    seqs[0].quality = NULL;

    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run_seqs(ctx, refs, 1, seqs, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 1);
    ASSERT_EQ_U64(out->num_aligned, 1);

    ASSERT_NOT_NULL(out->strand);
    ASSERT_NOT_NULL(out->score);
    ASSERT_NOT_NULL(out->edit_distance);

    /* reverse-complement read: strand must be 0. Score and edit distance
     * are those of the CLI on the same read (FLAG=16, AS:i:2430, NM:i:110);
     * this is the full 1487 nt read, so it scores higher than the truncated
     * test_read.fasta record. */
    ASSERT_EQ_INT(out->strand[0], 0);
    ASSERT_EQ_INT(out->score[0], 2430);
    ASSERT_EQ_INT(out->edit_distance[0], 110);

    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

TEST(test_run_null_out_does_not_crash) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[]  = { SMR_DATA_DIR "/test_ref.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/test_read.fasta" };
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, NULL, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_EQ_U64(stats.total_aligned, 1);
    smr_ctx_destroy(ctx);
}

/* ---- Pre-loaded index (streaming) API ---- */

TEST(test_index_load_rejects_null_ctx) {
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_index_t *idx = smr_index_load(NULL, refs, 1);
    ASSERT_NULL(idx);
}

TEST(test_index_load_rejects_null_refs) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    smr_index_t *idx = smr_index_load(ctx, NULL, 0);
    ASSERT_NULL(idx);
    ASSERT_EQ_INT(smr_last_error_code(ctx), SMR_ERR_INVALID_CONFIG);
    smr_ctx_destroy(ctx);
}

TEST(test_index_load_rejects_zero_num_refs) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_index_t *idx = smr_index_load(ctx, refs, 0);
    ASSERT_NULL(idx);
    ASSERT_EQ_INT(smr_last_error_code(ctx), SMR_ERR_INVALID_CONFIG);
    smr_ctx_destroy(ctx);
}

TEST(test_index_load_rejects_missing_ref_file) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { "/nonexistent/ref.fasta" };
    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NULL(idx);
    ASSERT_EQ_INT(smr_last_error_code(ctx), SMR_ERR_IO);
    smr_ctx_destroy(ctx);
}

TEST(test_index_free_is_null_safe) {
    smr_index_free(NULL);
    ASSERT_TRUE(1);
}

TEST(test_run_seqs_with_index_rejects_null_handle) {
    smr_seq_t seqs[1];
    seqs[0].id = "r"; seqs[0].sequence = "A"; seqs[0].quality = NULL;
    int rc = smr_run_seqs_with_index(NULL, seqs, 1, NULL, NULL);
    ASSERT_EQ_INT(rc, SMR_ERR_INVALID_CONFIG);
}

TEST(test_index_load_and_free_roundtrip) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);
    ASSERT_EQ_INT(smr_last_error_code(ctx), SMR_OK);
    smr_index_free(idx);
    smr_ctx_destroy(ctx);
}

/* Wrapper-equivalence test: the same inputs run through smr_run_seqs and
 * through (smr_index_load + smr_run_seqs_with_index + smr_index_free) must
 * produce byte-identical per-read output and counter-identical stats. */
TEST(test_run_seqs_with_index_matches_legacy_tiny) {
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_seq_t seqs[1];
    seqs[0].id = "AB271211";
    seqs[0].sequence =
        "TCCAACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTAG"
        "AGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGG"
        "GAATCTTCGGATTCTAGTGGCGGACGGGTGAGTAACGCGTAAGAATCTAACTTCAGGACG"
        "GGGACAACAGTGGGAAACGACTGCTAATACCCGATGTGCCGCGAGGTGAAACCTAATTGG";
    seqs[0].quality = NULL;

    /* Legacy path */
    smr_config_t cfg1;
    smr_config_init(&cfg1);
    cfg1.num_threads = 1;
    smr_context_t *ctx1 = smr_ctx_create(&cfg1);
    ASSERT_NOT_NULL(ctx1);
    smr_output_t *out1 = NULL;
    smr_stats_t stats1;
    int rc1 = smr_run_seqs(ctx1, refs, 1, seqs, 1, &out1, &stats1);
    ASSERT_EQ_INT(rc1, SMR_OK);
    ASSERT_NOT_NULL(out1);

    /* New handle path */
    smr_config_t cfg2;
    smr_config_init(&cfg2);
    cfg2.num_threads = 1;
    smr_context_t *ctx2 = smr_ctx_create(&cfg2);
    ASSERT_NOT_NULL(ctx2);
    smr_index_t *idx = smr_index_load(ctx2, refs, 1);
    ASSERT_NOT_NULL(idx);
    smr_output_t *out2 = NULL;
    smr_stats_t stats2;
    int rc2 = smr_run_seqs_with_index(idx, seqs, 1, &out2, &stats2);
    ASSERT_EQ_INT(rc2, SMR_OK);
    ASSERT_NOT_NULL(out2);

    /* Byte-compare every output field */
    ASSERT_EQ_U64(out1->num_reads, out2->num_reads);
    ASSERT_EQ_U64(out1->num_aligned, out2->num_aligned);
    for (uint64_t i = 0; i < out1->num_reads; i++) {
        ASSERT_STR_EQ(out1->read_ids[i], out2->read_ids[i]);
        ASSERT_EQ_INT(out1->aligned[i], out2->aligned[i]);
        ASSERT_EQ_INT(out1->ref_index[i], out2->ref_index[i]);
        ASSERT_EQ_INT(out1->ref_start[i], out2->ref_start[i]);
        ASSERT_EQ_INT(out1->ref_end[i], out2->ref_end[i]);
        ASSERT_EQ_INT(out1->strand[i], out2->strand[i]);
        ASSERT_EQ_INT(out1->score[i], out2->score[i]);
        ASSERT_EQ_INT(out1->edit_distance[i], out2->edit_distance[i]);
        /* Both paths run identical code; values must be bitwise equal. */
        ASSERT_DOUBLE_BITEQ(out1->e_value[i], out2->e_value[i]);
        ASSERT_DOUBLE_BITEQ(out1->identity[i], out2->identity[i]);
        ASSERT_DOUBLE_BITEQ(out1->coverage[i], out2->coverage[i]);
        if (out1->aligned[i]) {
            ASSERT_STR_EQ(out1->cigar[i], out2->cigar[i]);
            ASSERT_STR_EQ(out1->ref_name[i], out2->ref_name[i]);
        }
    }
    ASSERT_EQ_U64(stats1.total_reads, stats2.total_reads);
    ASSERT_EQ_U64(stats1.total_aligned, stats2.total_aligned);
    ASSERT_EQ_U64(stats1.total_id_cov_pass, stats2.total_id_cov_pass);
    ASSERT_EQ_U64(stats1.total_denovo, stats2.total_denovo);
    ASSERT_EQ_INT((int)stats1.min_read_len, (int)stats2.min_read_len);
    ASSERT_EQ_INT((int)stats1.max_read_len, (int)stats2.max_read_len);

    smr_output_free(out1);
    smr_output_free(out2);
    smr_index_free(idx);
    smr_ctx_destroy(ctx1);
    smr_ctx_destroy(ctx2);
}

/* Fixture: 6 reads from set7. First 4 align, last 2 ("random1"/"random2") do not. */
static const smr_seq_t SET7_SEQS[6] = {
    { "BD.ERD505_1",
      "AACGTAGGTGGCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGCAGGCGGAAAAGCAAGTTGGACGTGAAATCTATGGGCTCAACCCATAGCGTG",
      NULL },
    { "BD.NBS1076_0",
      "TACGGAGGGTGCAAGCGTTAATCCGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCACGGGCTCAACCTGG",
      NULL },
    { "LD.Glosor1_17",
      "TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTG",
      NULL },
    { "BD.ERD510_20",
      "TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTG",
      NULL },
    { "random1",
      "GTGTCACGTCAAATTCTCGGCTGGCTCCCTTAGTCGCATTAGTCCATGCAGAACGCGCACAGTTGAGGCAAGGCCGTAAAACACGTATGGATAAGGGGAT",
      NULL },
    { "random2",
      "TCACTTACGATATGCCTGTCTGGGGCCATCTCTAACGTCGGCGATGTTCCCATTCAGCGGCAAGCTCTCGTTCTGCATGGGTCAACTCCCTCACGAAGAA",
      NULL },
};

/* Three sequential batches on ONE handle must produce byte-identical
 * per-read output to three independent smr_run_seqs calls (fresh handle each
 * time). Catches state leakage across calls — Readstats restoration from
 * kvdb, stale kvdb entries, etc.
 *
 * Batch layout is interleaved (not 2+2+2 aligned-then-unaligned) so that each
 * batch mixes aligned and unaligned reads; a state contamination that
 * affected only the aligned-count would otherwise slip through. */
TEST(test_run_seqs_with_index_repeated_matches_independent) {
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };

    /* Indices into SET7_SEQS. reads[0..3] align, [4..5] do not; interleave so
     * every batch sees a mix. */
    const int B0[] = {0, 4}; /* aligned + unaligned */
    const int B1[] = {1, 5}; /* aligned + unaligned */
    const int B2[] = {2, 3}; /* aligned + aligned */
    const int *batches[3] = { B0, B1, B2 };

    smr_seq_t bufs[3][2];
    for (int b = 0; b < 3; b++) {
        bufs[b][0] = SET7_SEQS[batches[b][0]];
        bufs[b][1] = SET7_SEQS[batches[b][1]];
    }

    /* Path A: three independent smr_run_seqs calls (the control). */
    smr_output_t *ctl[3] = { NULL, NULL, NULL };
    smr_stats_t ctl_stats[3];
    for (int b = 0; b < 3; b++) {
        smr_config_t cfg;
        smr_config_init(&cfg);
        cfg.num_threads = 1;
        smr_context_t *ctx = smr_ctx_create(&cfg);
        ASSERT_NOT_NULL(ctx);
        int rc = smr_run_seqs(ctx, refs, 1, bufs[b], 2, &ctl[b], &ctl_stats[b]);
        ASSERT_EQ_INT(rc, SMR_OK);
        ASSERT_NOT_NULL(ctl[b]);
        smr_ctx_destroy(ctx);
    }

    /* Path B: three sequential with_index calls on one handle. */
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);

    smr_output_t *exp[3] = { NULL, NULL, NULL };
    smr_stats_t exp_stats[3];
    for (int b = 0; b < 3; b++) {
        int rc = smr_run_seqs_with_index(idx, bufs[b], 2, &exp[b], &exp_stats[b]);
        ASSERT_EQ_INT(rc, SMR_OK);
        ASSERT_NOT_NULL(exp[b]);
    }

    /* Byte-compare every batch. */
    for (int b = 0; b < 3; b++) {
        ASSERT_EQ_U64(ctl[b]->num_reads, exp[b]->num_reads);
        ASSERT_EQ_U64(ctl[b]->num_aligned, exp[b]->num_aligned);
        for (uint64_t i = 0; i < ctl[b]->num_reads; i++) {
            ASSERT_STR_EQ(ctl[b]->read_ids[i], exp[b]->read_ids[i]);
            ASSERT_EQ_INT(ctl[b]->aligned[i], exp[b]->aligned[i]);
            ASSERT_EQ_INT(ctl[b]->ref_index[i], exp[b]->ref_index[i]);
            ASSERT_EQ_INT(ctl[b]->ref_start[i], exp[b]->ref_start[i]);
            ASSERT_EQ_INT(ctl[b]->ref_end[i], exp[b]->ref_end[i]);
            ASSERT_EQ_INT(ctl[b]->strand[i], exp[b]->strand[i]);
            ASSERT_EQ_INT(ctl[b]->score[i], exp[b]->score[i]);
            ASSERT_EQ_INT(ctl[b]->edit_distance[i], exp[b]->edit_distance[i]);
            ASSERT_DOUBLE_BITEQ(ctl[b]->e_value[i], exp[b]->e_value[i]);
            ASSERT_DOUBLE_BITEQ(ctl[b]->identity[i], exp[b]->identity[i]);
            ASSERT_DOUBLE_BITEQ(ctl[b]->coverage[i], exp[b]->coverage[i]);
            if (ctl[b]->aligned[i]) {
                ASSERT_STR_EQ(ctl[b]->cigar[i], exp[b]->cigar[i]);
                ASSERT_STR_EQ(ctl[b]->ref_name[i], exp[b]->ref_name[i]);
            }
        }
        ASSERT_EQ_U64(ctl_stats[b].total_reads, exp_stats[b].total_reads);
        ASSERT_EQ_U64(ctl_stats[b].total_aligned, exp_stats[b].total_aligned);
        ASSERT_EQ_U64(ctl_stats[b].total_id_cov_pass, exp_stats[b].total_id_cov_pass);
    }

    /* Aggregate-sum cross-check: summed over batches must match what a
     * monolithic 6-read run would produce. Catches symmetric count leakage
     * between batches that per-batch equality would miss. */
    uint64_t sum_reads = 0, sum_aligned = 0;
    for (int b = 0; b < 3; b++) {
        sum_reads += exp_stats[b].total_reads;
        sum_aligned += exp_stats[b].total_aligned;
    }
    ASSERT_EQ_U64(sum_reads, 6);
    ASSERT_EQ_U64(sum_aligned, 4); /* reads 0..3 align, 4..5 don't */

    for (int b = 0; b < 3; b++) {
        smr_output_free(ctl[b]);
        smr_output_free(exp[b]);
    }
    smr_index_free(idx);
    smr_ctx_destroy(ctx);
}

/* Batch-splitting invariance. Running 6 reads as one batch of 6 vs
 * three batches of 2 on the same handle must produce per-read output that is
 * byte-identical, and summed stats that match the monolithic run. */
TEST(test_batch_split_invariance) {
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };

    /* Monolithic 1x6 run. */
    smr_config_t cfg_m;
    smr_config_init(&cfg_m);
    cfg_m.num_threads = 1;
    smr_context_t *ctx_m = smr_ctx_create(&cfg_m);
    ASSERT_NOT_NULL(ctx_m);
    smr_index_t *idx_m = smr_index_load(ctx_m, refs, 1);
    ASSERT_NOT_NULL(idx_m);
    smr_output_t *mono = NULL;
    smr_stats_t mono_stats;
    int rc = smr_run_seqs_with_index(idx_m, SET7_SEQS, 6, &mono, &mono_stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(mono);
    ASSERT_EQ_U64(mono->num_reads, 6);

    /* Split 3x2 run on a separate handle (interleaved to mix aligned+unaligned). */
    const int B0[] = {0, 4};
    const int B1[] = {1, 5};
    const int B2[] = {2, 3};
    const int *batches[3] = { B0, B1, B2 };
    smr_seq_t bufs[3][2];
    for (int b = 0; b < 3; b++) {
        bufs[b][0] = SET7_SEQS[batches[b][0]];
        bufs[b][1] = SET7_SEQS[batches[b][1]];
    }

    smr_config_t cfg_s;
    smr_config_init(&cfg_s);
    cfg_s.num_threads = 1;
    smr_context_t *ctx_s = smr_ctx_create(&cfg_s);
    ASSERT_NOT_NULL(ctx_s);
    smr_index_t *idx_s = smr_index_load(ctx_s, refs, 1);
    ASSERT_NOT_NULL(idx_s);

    smr_output_t *splits[3] = { NULL, NULL, NULL };
    smr_stats_t split_stats[3];
    for (int b = 0; b < 3; b++) {
        rc = smr_run_seqs_with_index(idx_s, bufs[b], 2, &splits[b], &split_stats[b]);
        ASSERT_EQ_INT(rc, SMR_OK);
        ASSERT_NOT_NULL(splits[b]);
    }

    /* Build a concatenated view keyed by read_id, matching mono's input order. */
    for (uint64_t mi = 0; mi < mono->num_reads; mi++) {
        /* find the same read_id in the split output */
        const char *id = mono->read_ids[mi];
        int found = 0;
        for (int b = 0; b < 3 && !found; b++) {
            for (uint64_t i = 0; i < splits[b]->num_reads; i++) {
                if (strcmp(splits[b]->read_ids[i], id) == 0) {
                    ASSERT_EQ_INT(mono->aligned[mi], splits[b]->aligned[i]);
                    ASSERT_EQ_INT(mono->ref_index[mi], splits[b]->ref_index[i]);
                    ASSERT_EQ_INT(mono->ref_start[mi], splits[b]->ref_start[i]);
                    ASSERT_EQ_INT(mono->ref_end[mi], splits[b]->ref_end[i]);
                    ASSERT_EQ_INT(mono->strand[mi], splits[b]->strand[i]);
                    ASSERT_EQ_INT(mono->score[mi], splits[b]->score[i]);
                    ASSERT_EQ_INT(mono->edit_distance[mi], splits[b]->edit_distance[i]);
                    ASSERT_DOUBLE_BITEQ(mono->e_value[mi], splits[b]->e_value[i]);
                    ASSERT_DOUBLE_BITEQ(mono->identity[mi], splits[b]->identity[i]);
                    ASSERT_DOUBLE_BITEQ(mono->coverage[mi], splits[b]->coverage[i]);
                    if (mono->aligned[mi]) {
                        ASSERT_STR_EQ(mono->cigar[mi], splits[b]->cigar[i]);
                        ASSERT_STR_EQ(mono->ref_name[mi], splits[b]->ref_name[i]);
                    }
                    found = 1;
                    break;
                }
            }
        }
        ASSERT_TRUE(found);
    }

    /* Aggregates match. */
    uint64_t sum_reads = 0, sum_aligned = 0;
    for (int b = 0; b < 3; b++) {
        sum_reads += split_stats[b].total_reads;
        sum_aligned += split_stats[b].total_aligned;
    }
    ASSERT_EQ_U64(sum_reads, mono_stats.total_reads);
    ASSERT_EQ_U64(sum_aligned, mono_stats.total_aligned);

    smr_output_free(mono);
    for (int b = 0; b < 3; b++) smr_output_free(splits[b]);
    smr_index_free(idx_m);
    smr_index_free(idx_s);
    smr_ctx_destroy(ctx_m);
    smr_ctx_destroy(ctx_s);
}

/* Log callback that counts occurrences of a given substring in emitted msgs. */
struct substr_counter {
    const char *needle;
    int count;
};
static void substr_count_cb(int level, const char *msg, void *user_data) {
    (void)level;
    struct substr_counter *sc = (struct substr_counter *)user_data;
    if (sc && msg && sc->needle && strstr(msg, sc->needle)) sc->count++;
}

/* An empty reference file must produce SMR_ERR_IO with a descriptive
 * last_error, not a crash or a successfully-loaded empty handle. */
TEST(test_index_load_empty_ref_returns_null) {
    char empty_path[512];
    snprintf(empty_path, sizeof(empty_path), "%s/smr_empty_ref_XXXXXX", tmp_base());
    int fd = mkstemp(empty_path);
    ASSERT_TRUE(fd >= 0);
    close(fd); /* leave the file empty (zero bytes) */

    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { empty_path };
    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NULL(idx);
    ASSERT_EQ_INT(smr_last_error_code(ctx), SMR_ERR_IO);
    const char *msg = smr_last_error(ctx);
    ASSERT_NOT_NULL(msg);
    ASSERT_TRUE(strstr(msg, "empty") != NULL);

    unlink(empty_path);
    smr_ctx_destroy(ctx);
}

/* Two threads alternating smr_run_seqs_with_index on the same
 * handle must be serialized by g_run_mutex and each produce correct output. */
struct handle_thread_arg {
    smr_index_t *idx;
    const smr_seq_t *seqs;
    int num_seqs;
    int rc;
    uint64_t num_aligned;
};
static void *handle_concurrent_worker(void *arg) {
    struct handle_thread_arg *a = (struct handle_thread_arg *)arg;
    for (int i = 0; i < 5; i++) {
        smr_output_t *out = NULL;
        smr_stats_t stats;
        a->rc = smr_run_seqs_with_index(a->idx, a->seqs, a->num_seqs, &out, &stats);
        if (a->rc != SMR_OK) { if (out) smr_output_free(out); return NULL; }
        a->num_aligned += stats.total_aligned;
        smr_output_free(out);
    }
    return NULL;
}
TEST(test_concurrent_batches_on_shared_handle) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);

    struct handle_thread_arg a1 = { idx, &SET7_SEQS[0], 2, 0, 0 };
    struct handle_thread_arg a2 = { idx, &SET7_SEQS[2], 2, 0, 0 };
    pthread_t t1, t2;
    pthread_create(&t1, NULL, handle_concurrent_worker, &a1);
    pthread_create(&t2, NULL, handle_concurrent_worker, &a2);
    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
    ASSERT_EQ_INT(a1.rc, SMR_OK);
    ASSERT_EQ_INT(a2.rc, SMR_OK);
    /* Each worker ran 5 batches; reads 0,1 both align and reads 2,3 both
     * align; so each worker should report 5 * 2 = 10 alignments total. */
    ASSERT_EQ_U64(a1.num_aligned, 10);
    ASSERT_EQ_U64(a2.num_aligned, 10);

    smr_index_free(idx);
    smr_ctx_destroy(ctx);
}

/* Regression: smr_run_seqs_with_index accepts FASTQ batches (non-NULL
 * quality strings). The placeholder reads file written by smr_index_load
 * is FASTA, but the in-memory Readfeed derives the format from the batch's
 * quality strings, so the placeholder's format doesn't matter. */
TEST(test_run_seqs_with_index_accepts_fastq) {
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta" };
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);
    smr_seq_t seqs[1];
    seqs[0].id = "r1";
    seqs[0].sequence = "ACGT";
    seqs[0].quality  = "IIII"; /* non-NULL → FASTQ */
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run_seqs_with_index(idx, seqs, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);
    ASSERT_EQ_U64(out->num_reads, 1);
    smr_output_free(out);
    smr_index_free(idx);
    smr_ctx_destroy(ctx);
}

/* Index + References load exactly once (at smr_index_load time).
 * Subsequent smr_run_seqs_with_index calls on the same handle must NOT
 * re-emit the "Loading references" or "Loading index:" log lines — that's
 * the performance payoff the handle exists to deliver. */
TEST(test_repeated_calls_do_not_reload_refs) {
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };

    struct substr_counter sc = { "loading references and index into memory", 0 };
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.log_callback = substr_count_cb;
    cfg.log_user_data = &sc;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);

    /* Load: should emit "Loading references" (at least once per index-part). */
    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);
    int load_count = sc.count;
    ASSERT_TRUE(load_count >= 1);

    /* Run 3 batches. The counter must NOT increase — refs stay loaded. */
    for (int b = 0; b < 3; b++) {
        smr_output_t *out = NULL;
        smr_stats_t stats;
        int rc = smr_run_seqs_with_index(idx, &SET7_SEQS[b * 2], 2, &out, &stats);
        ASSERT_EQ_INT(rc, SMR_OK);
        smr_output_free(out);
    }
    ASSERT_EQ_INT(sc.count, load_count);

    smr_index_free(idx);
    smr_ctx_destroy(ctx);
}

/* ---- E-values, thresholds, workdir reuse ---- */

/* smr_run e-values use the CLI definition, so they must match the e-values
 * of the aligned.blast report written by the same run. The report prints
 * three significant digits. */
TEST(test_run_evalue_matches_blast_report) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_NOT_NULL(out);

    char path[1024];
    snprintf(path, sizeof(path), "%s/aligned.blast", wdir);
    FILE *f = fopen(path, "r");
    ASSERT_NOT_NULL(f);
    char line[4096];
    int checked = 0;
    while (fgets(line, sizeof(line), f)) {
        /* columns: qseqid sseqid pident length mismatch gapopen qstart qend
         * sstart send evalue ... */
        char qid[256];
        double ev = -1.0;
        if (sscanf(line, "%255s %*s %*s %*s %*s %*s %*s %*s %*s %*s %lf", qid, &ev) != 2)
            continue;
        for (uint64_t i = 0; i < out->num_reads; i++) {
            if (strcmp(out->read_ids[i], qid) != 0) continue;
            ASSERT_EQ_INT(out->aligned[i], 1);
            ASSERT_DOUBLE_NEAR(out->e_value[i], ev, ev * 5e-3);
            checked++;
        }
    }
    fclose(f);
    ASSERT_EQ_INT(checked, 4);

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

/* cfg.evalue reaches the CLI option parser unchanged. The expected counts
 * are those of the CLI with -e 1e-10 and -e 1e-20: BD.ERD505_1 (E=2.44e-13)
 * passes the first threshold but not the second. */
TEST(test_run_evalue_threshold) {
    const double evalues[2] = { 1e-10, 1e-20 };
    const uint64_t expected[2] = { 4, 3 };
    for (int k = 0; k < 2; k++) {
        smr_config_t cfg;
        smr_config_init(&cfg);
        cfg.num_threads = 1;
        cfg.evalue = evalues[k];
        smr_context_t *ctx = smr_ctx_create(&cfg);
        const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
        const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
        smr_output_t *out = NULL;
        smr_stats_t stats;
        int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
        ASSERT_EQ_INT(rc, SMR_OK);
        ASSERT_NOT_NULL(out);
        ASSERT_EQ_U64(out->num_aligned, expected[k]);
        ASSERT_EQ_U64(stats.total_aligned, expected[k]);
        smr_output_free(out);
        smr_ctx_destroy(ctx);
    }
}

/* Compare one read of two outputs field by field (doubles bitwise). */
static int same_read(const smr_output_t *a, uint64_t i, const smr_output_t *b, uint64_t j) {
    if (strcmp(a->read_ids[i], b->read_ids[j]) != 0) return 0;
    if (a->aligned[i] != b->aligned[j] || a->ref_index[i] != b->ref_index[j]) return 0;
    if (a->ref_start[i] != b->ref_start[j] || a->ref_end[i] != b->ref_end[j]) return 0;
    if (a->strand[i] != b->strand[j] || a->score[i] != b->score[j]) return 0;
    if (a->edit_distance[i] != b->edit_distance[j]) return 0;
    if (memcmp(&a->e_value[i], &b->e_value[j], sizeof(double)) != 0) return 0;
    if (memcmp(&a->identity[i], &b->identity[j], sizeof(double)) != 0) return 0;
    if (memcmp(&a->coverage[i], &b->coverage[j], sizeof(double)) != 0) return 0;
    if (a->aligned[i]) {
        if (strcmp(a->cigar[i], b->cigar[j]) != 0) return 0;
        if (strcmp(a->ref_name[i], b->ref_name[j]) != 0) return 0;
    }
    return 1;
}

static int dir_exists(const char *parent, const char *name) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s", parent, name);
    return access(path, F_OK) == 0;
}

/* Two smr_run calls on one workdir: the second loads the index built by the
 * first, both return the same results, and neither leaves a kvdb behind. */
TEST(test_run_reuses_workdir_index) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    struct substr_counter sc = { "Skipping indexing", 0 };
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    cfg.log_callback = substr_count_cb;
    cfg.log_user_data = &sc;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    smr_output_t *out1 = NULL, *out2 = NULL;
    smr_stats_t stats;

    ASSERT_EQ_INT(smr_run(ctx, refs, 1, reads, 1, &out1, &stats), SMR_OK);
    ASSERT_EQ_INT(sc.count, 0);
    ASSERT_TRUE(!dir_exists(wdir, "kvdb"));
    ASSERT_TRUE(dir_exists(wdir, "idx"));

    ASSERT_EQ_INT(smr_run(ctx, refs, 1, reads, 1, &out2, &stats), SMR_OK);
    ASSERT_EQ_INT(sc.count, 1);
    ASSERT_TRUE(!dir_exists(wdir, "kvdb"));

    ASSERT_EQ_U64(out1->num_reads, out2->num_reads);
    ASSERT_EQ_U64(out1->num_aligned, out2->num_aligned);
    for (uint64_t i = 0; i < out1->num_reads; i++)
        ASSERT_TRUE(same_read(out1, i, out2, i));

    smr_output_free(out1);
    smr_output_free(out2);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

/* Same for handles: a second smr_index_load on the workdir reuses the index,
 * and smr_index_free removes the kvdb. */
TEST(test_index_load_reuses_workdir_index) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    struct substr_counter sc = { "Skipping indexing", 0 };
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    cfg.log_callback = substr_count_cb;
    cfg.log_user_data = &sc;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    smr_output_t *out1 = NULL, *out2 = NULL;
    smr_stats_t stats;

    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);
    ASSERT_EQ_INT(sc.count, 0);
    ASSERT_EQ_INT(smr_run_seqs_with_index(idx, SET7_SEQS, 6, &out1, &stats), SMR_OK);
    smr_index_free(idx);
    ASSERT_TRUE(!dir_exists(wdir, "kvdb"));

    idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);
    ASSERT_EQ_INT(sc.count, 1);
    ASSERT_EQ_INT(smr_run_seqs_with_index(idx, SET7_SEQS, 6, &out2, &stats), SMR_OK);
    smr_index_free(idx);
    ASSERT_TRUE(!dir_exists(wdir, "kvdb"));

    for (uint64_t i = 0; i < out1->num_reads; i++)
        ASSERT_TRUE(same_read(out1, i, out2, i));

    smr_output_free(out1);
    smr_output_free(out2);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

/* A kvdb left in the workdir by an earlier process that died before cleanup
 * must not leak into the next run. The child aligns the four set7 reads that
 * hit and exits without smr_index_free, leaving their results in
 * <workdir>/kvdb under the same internal keys that smr_run then uses for
 * random1 and random2, which do not align. */
TEST(test_run_discards_leftover_kvdb) {
    char *wdir = make_tmpdir();
    ASSERT_NOT_NULL(wdir);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };

    fflush(stdout);
    pid_t pid = fork();
    ASSERT_TRUE(pid >= 0);
    if (pid == 0) {
        smr_config_t cfg;
        smr_config_init(&cfg);
        cfg.num_threads = 1;
        cfg.workdir = wdir;
        smr_context_t *ctx = smr_ctx_create(&cfg);
        smr_index_t *idx = smr_index_load(ctx, refs, 1);
        smr_output_t *out = NULL;
        smr_stats_t stats;
        int rc = idx ? smr_run_seqs_with_index(idx, SET7_SEQS, 4, &out, &stats) : -1;
        _exit(rc == SMR_OK && out && out->num_aligned == 4 ? 0 : 1);
    }
    int status = 0;
    ASSERT_EQ_INT(waitpid(pid, &status, 0), pid);
    ASSERT_TRUE(WIFEXITED(status) && WEXITSTATUS(status) == 0);
    ASSERT_TRUE(dir_exists(wdir, "kvdb"));

    char reads_path[1024];
    snprintf(reads_path, sizeof(reads_path), "%s/unaligned.fa", wdir);
    FILE *f = fopen(reads_path, "w");
    ASSERT_NOT_NULL(f);
    for (int i = 4; i < 6; i++)
        fprintf(f, ">%s\n%s\n", SET7_SEQS[i].id, SET7_SEQS[i].sequence);
    fclose(f);

    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    cfg.workdir = wdir;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *reads[] = { reads_path };
    smr_output_t *out = NULL;
    smr_stats_t stats;
    ASSERT_EQ_INT(smr_run(ctx, refs, 1, reads, 1, &out, &stats), SMR_OK);
    ASSERT_EQ_U64(out->num_reads, 2);
    ASSERT_EQ_U64(out->num_aligned, 0);
    ASSERT_EQ_INT(out->aligned[0], 0);
    ASSERT_EQ_INT(out->aligned[1], 0);

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    rm_rf(wdir);
}

/* A handle holds one loaded index; more than one reference file is refused
 * rather than silently aligning against the first only. */
TEST(test_index_load_rejects_multiple_refs) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    const char *refs[] = { SMR_DATA_DIR "/test_ref.fasta",
                           SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    smr_index_t *idx = smr_index_load(ctx, refs, 2);
    ASSERT_NULL(idx);
    ASSERT_EQ_INT(smr_last_error_code(ctx), SMR_ERR_NOT_IMPLEMENTED);
    smr_ctx_destroy(ctx);
}

/* With no log callback, nothing reaches stdout on the handle path either
 * (index and reference loading run outside the fd redirect). */
TEST(test_no_stdout_during_index_load) {
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 1;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    smr_output_t *out = NULL;
    smr_stats_t stats;

    fflush(stdout);
    int old_stdout = dup(1);
    FILE *tmp = tmpfile();
    dup2(fileno(tmp), 1);

    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    int rc = idx ? smr_run_seqs_with_index(idx, SET7_SEQS, 6, &out, &stats) : -1;
    smr_index_free(idx);

    fflush(stdout);
    dup2(old_stdout, 1);
    close(old_stdout);
    fseek(tmp, 0, SEEK_END);
    long captured_size = ftell(tmp);
    fclose(tmp);

    smr_output_free(out);
    smr_ctx_destroy(ctx);
    ASSERT_EQ_INT(rc, SMR_OK);
    ASSERT_EQ_INT((int)captured_size, 0);
}

/* Batch-splitting invariance with several processing threads. The same
 * reads as one batch, as consecutive smaller batches, and with one thread
 * must give identical per-read output in input order. */
static void check_split_invariance(int paired, int threads, int batch) {
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    const int n = 6;
    smr_output_t *mono = NULL, *single = NULL, *parts[6] = { NULL };
    smr_stats_t stats;
    int nparts = n / batch;

    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = threads;
    cfg.paired = paired;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);
    ASSERT_EQ_INT(smr_run_seqs_with_index(idx, SET7_SEQS, n, &mono, &stats), SMR_OK);
    for (int b = 0; b < nparts; b++)
        ASSERT_EQ_INT(smr_run_seqs_with_index(idx, &SET7_SEQS[b * batch], batch,
                                              &parts[b], &stats), SMR_OK);
    smr_index_free(idx);
    smr_ctx_destroy(ctx);

    smr_config_t cfg1 = cfg;
    cfg1.num_threads = 1;
    smr_context_t *ctx1 = smr_ctx_create(&cfg1);
    ASSERT_NOT_NULL(ctx1);
    ASSERT_EQ_INT(smr_run_seqs(ctx1, refs, 1, SET7_SEQS, n, &single, &stats), SMR_OK);
    smr_ctx_destroy(ctx1);

    ASSERT_EQ_U64(mono->num_reads, (uint64_t)n);
    ASSERT_EQ_U64(single->num_reads, (uint64_t)n);
    ASSERT_EQ_U64(mono->num_aligned, 4);
    for (int i = 0; i < n; i++) {
        ASSERT_STR_EQ(mono->read_ids[i], SET7_SEQS[i].id);
        ASSERT_TRUE(same_read(mono, i, single, i));
        ASSERT_TRUE(same_read(mono, i, parts[i / batch], i % batch));
    }

    smr_output_free(mono);
    smr_output_free(single);
    for (int b = 0; b < nparts; b++) smr_output_free(parts[b]);
}

TEST(test_batch_split_invariance_multithreaded) {
    check_split_invariance(0, 3, 2);
    if (_test_current_failed) return;
    check_split_invariance(0, 4, 3);
}

TEST(test_batch_split_invariance_paired_multithreaded) {
    check_split_invariance(1, 2, 2);
    if (_test_current_failed) return;
    check_split_invariance(1, 3, 6);
}

/* Worker threads log through the caller's callback: each alignment worker
 * logs one "started (start_idx=" line per pass. */
TEST(test_log_callback_gets_worker_messages) {
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    struct substr_counter sc = { " started (start_idx=", 0 };
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 3;
    cfg.log_callback = substr_count_cb;
    cfg.log_user_data = &sc;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);

    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);
    smr_output_t *out = NULL;
    smr_stats_t stats;
    ASSERT_EQ_INT(smr_run_seqs_with_index(idx, SET7_SEQS, 6, &out, &stats), SMR_OK);
    ASSERT_EQ_INT(sc.count, 3);
    smr_output_free(out);
    smr_index_free(idx);

    sc.count = 0;
    const char *reads[] = { SMR_DATA_DIR "/set7_arc_bac_16S_database_match.fasta" };
    out = NULL;
    ASSERT_EQ_INT(smr_run(ctx, refs, 1, reads, 1, &out, &stats), SMR_OK);
    ASSERT_EQ_INT(sc.count, 3);
    smr_output_free(out);
    smr_ctx_destroy(ctx);
}

/* The library never runs the callback concurrently, although several
 * worker threads log at the same time. The callback sleeps to widen the
 * window in which an unserialized second call would overlap. */
struct overlap_probe {
    int active;
    int max_active;
    int calls;
};
static void overlap_cb(int level, const char *msg, void *user_data) {
    (void)level; (void)msg;
    struct overlap_probe *p = (struct overlap_probe *)user_data;
    int now = __atomic_add_fetch(&p->active, 1, __ATOMIC_SEQ_CST);
    if (now > __atomic_load_n(&p->max_active, __ATOMIC_SEQ_CST))
        __atomic_store_n(&p->max_active, now, __ATOMIC_SEQ_CST);
    usleep(2000);
    __atomic_add_fetch(&p->calls, 1, __ATOMIC_SEQ_CST);
    __atomic_sub_fetch(&p->active, 1, __ATOMIC_SEQ_CST);
}
TEST(test_log_callback_not_concurrent) {
    const char *refs[] = { SMR_DATA_DIR "/silva-arc-16s-database-id95.fasta" };
    struct overlap_probe probe = { 0, 0, 0 };
    smr_config_t cfg;
    smr_config_init(&cfg);
    cfg.num_threads = 4;
    cfg.log_callback = overlap_cb;
    cfg.log_user_data = &probe;
    smr_context_t *ctx = smr_ctx_create(&cfg);
    ASSERT_NOT_NULL(ctx);
    smr_index_t *idx = smr_index_load(ctx, refs, 1);
    ASSERT_NOT_NULL(idx);
    for (int b = 0; b < 3; b++) {
        smr_output_t *out = NULL;
        smr_stats_t stats;
        ASSERT_EQ_INT(smr_run_seqs_with_index(idx, SET7_SEQS, 6, &out, &stats), SMR_OK);
        smr_output_free(out);
    }
    smr_index_free(idx);
    smr_ctx_destroy(ctx);
    ASSERT_TRUE(probe.calls > 0);
    ASSERT_EQ_INT(probe.max_active, 1);
}

TEST_MAIN_BEGIN()
    RUN_TEST(test_config_init_sets_struct_size);
    RUN_TEST(test_config_struct_size_is_first_field);
    RUN_TEST(test_config_init_zeroes_pointers);
    RUN_TEST(test_version_matches_version_h);
    RUN_TEST(test_strerror_returns_string);
    RUN_TEST(test_output_free_null_safe);
    RUN_TEST(test_ctx_destroy_null_safe);
    RUN_TEST(test_config_default_threads);
    RUN_TEST(test_config_default_num_alignments);
    RUN_TEST(test_config_default_match);
    RUN_TEST(test_config_default_mismatch);
    RUN_TEST(test_config_default_gap_open);
    RUN_TEST(test_config_default_gap_ext);
    RUN_TEST(test_config_default_score_N);
    RUN_TEST(test_config_default_evalue);
    RUN_TEST(test_config_default_seed_win_len);
    RUN_TEST(test_config_default_best);
    RUN_TEST(test_config_default_booleans_off);
    RUN_TEST(test_config_boolean_type_is_int32);
    RUN_TEST(test_config_explicit_width_types);
    RUN_TEST(test_ctx_create_returns_non_null);
    RUN_TEST(test_ctx_create_null_config_returns_null);
    RUN_TEST(test_ctx_create_bad_struct_size_returns_null);
    RUN_TEST(test_ctx_create_multiple_independent);
    RUN_TEST(test_ctx_last_error_empty_initially);
    RUN_TEST(test_strerror_success_msg);
    RUN_TEST(test_strerror_invalid_config_msg);
    RUN_TEST(test_strerror_unknown_code);
    RUN_TEST(test_log_callback_receives_messages);
    RUN_TEST(test_last_error_set_after_smr_run);
    RUN_TEST(test_last_error_code_initially_zero);
    RUN_TEST(test_last_error_code_null_ctx);
    RUN_TEST(test_two_contexts_create_destroy);
    RUN_TEST(test_run_null_refs_returns_error);
    RUN_TEST(test_run_null_reads_returns_error);
    RUN_TEST(test_run_nonexistent_ref_returns_error);
    RUN_TEST(test_run_nonexistent_reads_returns_error);
    RUN_TEST(test_run_empty_ref_returns_error);
    RUN_TEST(test_last_error_descriptive_after_bad_input);
    RUN_TEST(test_run_bad_input_does_not_crash);
    RUN_TEST(test_run_tiny_aligned_count);
    RUN_TEST(test_run_tiny_stats);
    RUN_TEST(test_run_small_aligned_count);
    RUN_TEST(test_run_output_free_after_run);
    RUN_TEST(test_run_multiple_sequential);
    RUN_TEST(test_no_stdout_during_run);
    RUN_TEST(test_log_callback_fires_during_run);
    RUN_TEST(test_null_log_callback_silent);
    /* golden file comparison */
    RUN_TEST(test_golden_tiny_blast);
    RUN_TEST(test_golden_tiny_fasta);
    RUN_TEST(test_golden_tiny_sam);
    RUN_TEST(test_golden_small_blast);
    RUN_TEST(test_golden_small_fasta);
    RUN_TEST(test_golden_small_sam);
    RUN_TEST(test_golden_small_other);
    /* concurrent smr_run */
    RUN_TEST(test_concurrent_smr_run);
    /* multi-threaded alignment */
    RUN_TEST(test_run_tiny_multithreaded);
    /* per-read output */
    RUN_TEST(test_run_tiny_per_read_output);
    RUN_TEST(test_run_small_per_read_output);
    /* in-memory input */
    RUN_TEST(test_run_seqs_tiny);
    RUN_TEST(test_run_seqs_null_error);
    RUN_TEST(test_run_seqs_empty_seq_error);
    RUN_TEST(test_run_seqs_quality_length_mismatch);
    RUN_TEST(test_run_seqs_paired_odd_count_error);
    RUN_TEST(test_run_seqs_paired);
    /* ref_name output */
    RUN_TEST(test_run_tiny_ref_name);
    RUN_TEST(test_run_small_ref_name);
    /* strand / score / edit_distance */
    RUN_TEST(test_run_tiny_strand_score_edit);
    RUN_TEST(test_run_small_strand_score_edit);
    RUN_TEST(test_run_reverse_strand);
    RUN_TEST(test_run_null_out_does_not_crash);
    /* pre-loaded index (streaming) API */
    RUN_TEST(test_index_load_rejects_null_ctx);
    RUN_TEST(test_index_load_rejects_null_refs);
    RUN_TEST(test_index_load_rejects_zero_num_refs);
    RUN_TEST(test_index_load_rejects_missing_ref_file);
    RUN_TEST(test_index_free_is_null_safe);
    RUN_TEST(test_run_seqs_with_index_rejects_null_handle);
    RUN_TEST(test_index_load_and_free_roundtrip);
    RUN_TEST(test_run_seqs_with_index_matches_legacy_tiny);
    RUN_TEST(test_run_seqs_with_index_repeated_matches_independent);
    RUN_TEST(test_batch_split_invariance);
    RUN_TEST(test_repeated_calls_do_not_reload_refs);
    RUN_TEST(test_run_seqs_with_index_accepts_fastq);
    RUN_TEST(test_index_load_empty_ref_returns_null);
    RUN_TEST(test_concurrent_batches_on_shared_handle);
    RUN_TEST(test_index_load_rejects_multiple_refs);
    RUN_TEST(test_no_stdout_during_index_load);
    RUN_TEST(test_batch_split_invariance_multithreaded);
    RUN_TEST(test_batch_split_invariance_paired_multithreaded);
    /* e-values, thresholds, workdir reuse */
    RUN_TEST(test_run_evalue_matches_blast_report);
    RUN_TEST(test_run_evalue_threshold);
    RUN_TEST(test_run_reuses_workdir_index);
    RUN_TEST(test_index_load_reuses_workdir_index);
    RUN_TEST(test_run_discards_leftover_kvdb);
    RUN_TEST(test_log_callback_gets_worker_messages);
    RUN_TEST(test_log_callback_not_concurrent);
TEST_MAIN_END()
