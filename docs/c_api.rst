C Library API
=============

SortMeRNA provides a reentrant C API for embedding rRNA filtering into other
applications without spawning a subprocess. The API is defined in
``include/smr_api.h`` and links as a static library (``libsmr_api.a``).

The library never calls ``exit()``, ``abort()``, or ``assert()``. All output
is routed through a caller-provided log callback; stdout and stderr are
suppressed during pipeline execution.

There are two ways to use it:

- **One-shot**: ``smr_run`` (reads from files) and ``smr_run_seqs`` (reads
  from memory) build or load the index, align one set of reads, and return
  per-read results.
- **Streaming**: ``smr_index_load`` loads the index and references once;
  ``smr_run_seqs_with_index`` then aligns any number of in-memory batches
  against them; ``smr_index_free`` releases them. This is for tools that
  align many batches against the same references (a database table
  function, a server handling requests, a FASTQ reader feeding chunks).

.. contents:: On this page
   :local:
   :depth: 2

Quick start
-----------

.. code-block:: c

   #include "smr_api.h"

   int main(void) {
       /* 1. Configure */
       smr_config_t cfg;
       smr_config_init(&cfg);
       cfg.num_threads = 4;

       /* 2. Create context */
       smr_context_t *ctx = smr_ctx_create(&cfg);
       if (!ctx) return 1;

       /* 3. Run alignment */
       const char *refs[]  = { "silva-bac-16s.fasta" };
       const char *reads[] = { "sample.fastq" };
       smr_output_t *out = NULL;
       smr_stats_t stats;

       int rc = smr_run(ctx, refs, 1, reads, 1, &out, &stats);
       if (rc != SMR_OK) {
           fprintf(stderr, "Error: %s\n", smr_last_error(ctx));
           smr_ctx_destroy(ctx);
           return 1;
       }

       printf("Total reads: %llu\n", (unsigned long long)stats.total_reads);
       printf("Aligned:     %llu\n", (unsigned long long)out->num_aligned);

       /* 4. Inspect per-read results */
       for (uint64_t i = 0; i < out->num_reads; i++) {
           if (out->aligned[i])
               printf("%s: %.1f%% id, e=%.2g, strand=%s, AS=%d, NM=%d\n",
                      out->read_ids[i], out->identity[i], out->e_value[i],
                      out->strand[i] ? "+" : "-",
                      out->score[i], out->edit_distance[i]);
       }

       /* 5. Clean up */
       smr_output_free(out);
       smr_ctx_destroy(ctx);
       return 0;
   }

In-memory input:

.. code-block:: c

   /* Pass sequences directly without writing temp files */
   smr_seq_t seqs[] = {
       { "read1", "ACGTACGTACGT...", NULL },
       { "read2", "TGCATGCATGCA...", NULL },
   };
   int rc = smr_run_seqs(ctx, refs, 1, seqs, 2, &out, &stats);

Streaming:

.. code-block:: c

   cfg.workdir = "/data/smr_work";  /* keep the index between runs */
   smr_context_t *ctx = smr_ctx_create(&cfg);
   const char *refs[] = { "silva-bac-16s.fasta" };

   smr_index_t *idx = smr_index_load(ctx, refs, 1);
   if (!idx) { /* smr_last_error(ctx) */ }

   for (int b = 0; b < num_batches; b++) {
       smr_output_t *out = NULL;
       smr_stats_t stats;
       /* batch_seqs[b]: caller-owned smr_seq_t array of batch_len[b] reads */
       int rc = smr_run_seqs_with_index(idx, batch_seqs[b], batch_len[b], &out, &stats);
       if (rc != SMR_OK) { /* smr_last_error(ctx) */ continue; }
       /* out->aligned[i], out->ref_name[i], out->e_value[i], ... */
       smr_output_free(out);
   }

   smr_index_free(idx);
   smr_ctx_destroy(ctx);

Building
--------

The library is built together with the ``sortmerna`` binary. Configure the
build as described in :doc:`building`, then::

   cmake --build build --target smr_api

With ``-DWITH_TESTS=ON`` the test program ``test_smr_api`` is built too, and
``ctest`` runs it.

``build/src/smr_api/libsmr_api.a`` contains the API and the SortMeRNA core
objects. An application also links ``libbuild_version.a`` and
``3rdparty/alp/libalp.a`` from the build tree, Parasail, RocksDB and zlib,
plus ``-lpthread -ldl -lstdc++ -lm``. Within the SortMeRNA CMake project,
linking the ``smr_api`` target pulls these in. The SortMeRNA CMake files
expect to be the top-level project, so another CMake project builds it with
``ExternalProject_Add`` and links the archives above.


API reference
-------------

Configuration
#############

.. c:function:: void smr_config_init(smr_config_t *cfg)

   Initialize a configuration struct with sensible defaults. Must be called
   before modifying any fields. Sets ``struct_size`` (``uint32_t``) for ABI
   version detection.

   After calling ``smr_config_init``, the struct is ready for immediate use
   without further modification.

.. c:type:: smr_config_t

   Configuration struct. The first field (``struct_size``) must always be set
   by ``smr_config_init`` -- never initialize this struct manually.

   **Threading:**

   ========================  ===========  ===========
   Field                     Type         Default
   ========================  ===========  ===========
   ``num_threads``           ``int32_t``  2
   ========================  ===========  ===========

   **Alignment scoring:**

   ========================  ===========  ===========
   Field                     Type         Default
   ========================  ===========  ===========
   ``match``                 ``int32_t``  2
   ``mismatch``              ``int32_t``  -3
   ``gap_open``              ``int32_t``  5
   ``gap_ext``               ``int32_t``  2
   ``score_N``               ``int32_t``  -3
   ``evalue``                ``double``   ``SMR_EVALUE_DEFAULT`` (-1.0)
   ``seed_win_len``          ``uint32_t`` 18
   ``num_alignments``        ``uint32_t`` 1
   ========================  ===========  ===========

   ``evalue`` is the E-value threshold of ``smr_run``; ``SMR_EVALUE_DEFAULT``
   selects the sortmerna default (1e-5). The in-memory functions apply no
   threshold (see `E-value semantics`_).

   **Boolean flags** (``int32_t``, 0 = off, nonzero = on):

   ========================  ===========
   Field                     Default
   ========================  ===========
   ``best``                  1 (on)
   ``paired``                0
   ``forward_only``          0
   ``reverse_only``          0
   ``full_search``           0
   ========================  ===========

   **Paths:**

   ========================  ===============  ===========
   Field                     Type             Default
   ========================  ===============  ===========
   ``workdir``               ``const char*``  NULL (auto temp dir)
   ========================  ===============  ===========

   When ``workdir`` is NULL, each ``smr_run`` call and each handle uses a
   new directory under the system temp directory (``TMPDIR``), removed
   afterwards. The index is then rebuilt every time.

   Reuse one ``workdir`` across calls to build the reference index only
   once. The index is written to ``<workdir>/idx`` and loaded from there by
   later calls with the same reference file names. As with the sortmerna
   CLI, the index is not rebuilt when a reference file's contents or
   ``seed_win_len`` change; delete ``<workdir>/idx`` in that case.

   ``<workdir>/kvdb`` is scratch space. The library deletes it before and
   after each ``smr_run`` and when a handle is loaded and freed, so results
   never carry over between calls. The same ``workdir`` must not be used by
   two calls or handles at the same time. ``smr_run`` also writes its
   reports (``aligned.*``, ``other.*``) into the workdir.

   **Logging:**

   ========================  =================================  ===========
   Field                     Type                               Default
   ========================  =================================  ===========
   ``log_callback``          ``void(*)(int,const char*,void*)`` NULL
   ``log_user_data``         ``void*``                          NULL
   ========================  =================================  ===========

   When ``log_callback`` is NULL, the library operates silently. When set,
   all informational, warning, and error messages are routed through the
   callback. The ``log_user_data`` pointer is passed through unchanged.

   The callback also receives the messages of SortMeRNA's worker threads,
   so it may run on a thread other than the caller's. The library never
   runs it concurrently: calls are made one at a time. The callback must
   not call SortMeRNA functions.

   Both ``workdir`` and ``log_user_data`` are borrowed pointers -- the pointee
   must outlive the context.


Context lifecycle
#################

.. c:function:: smr_context_t* smr_ctx_create(const smr_config_t *cfg)

   Create a new context from a configuration struct. Returns an opaque pointer
   on success, or NULL on failure (invalid config or allocation failure).

   The config is copied into the context -- the caller may free or reuse the
   ``smr_config_t`` struct immediately after this call.

   :param cfg: Pointer to an initialized config struct. Must not be NULL.
               ``cfg->struct_size`` must equal ``sizeof(smr_config_t)``.
   :returns: Opaque context pointer, or NULL on failure.

.. c:function:: void smr_ctx_destroy(smr_context_t *ctx)

   Destroy a context and free all associated resources. Passing NULL is safe
   (no-op). The context must not be used after this call.


Error reporting
###############

.. c:function:: const char* smr_strerror(int code)

   Return a static string describing an error code category.

   ================================  =============================
   Code                              Message
   ================================  =============================
   ``SMR_OK`` (0)                    ``"Success"``
   ``SMR_ERR_INVALID_CONFIG`` (-1)   ``"Invalid configuration"``
   ``SMR_ERR_ALLOC`` (-2)            ``"Memory allocation failed"``
   ``SMR_ERR_IO`` (-3)               ``"I/O error"``
   ``SMR_ERR_INDEX`` (-4)            ``"Index error"``
   ``SMR_ERR_ALIGN`` (-5)            ``"Alignment error"``
   ``SMR_ERR_NOT_IMPLEMENTED`` (-99) ``"Not implemented"``
   ================================  =============================

.. c:function:: const char* smr_last_error(const smr_context_t *ctx)

   Return a detailed, context-specific error message from the most recent
   failing call. Returns an empty string if no error has occurred or if
   ``ctx`` is NULL.

.. c:function:: int smr_last_error_code(const smr_context_t *ctx)

   Return the numeric error code from the most recent failing call.
   Returns ``SMR_OK`` if no error has occurred. Returns
   ``SMR_ERR_INVALID_CONFIG`` if ``ctx`` is NULL.


Computation
###########

.. c:function:: int smr_run(smr_context_t *ctx, const char **ref_paths, int32_t num_refs, const char **read_paths, int32_t num_reads, smr_output_t **out, smr_stats_t *stats)

   Run the SortMeRNA alignment pipeline, as the ``sortmerna`` binary does.

   :param ctx: Context created by ``smr_ctx_create``.
   :param ref_paths: Array of reference FASTA file paths.
   :param num_refs: Number of reference files (must be > 0).
   :param read_paths: Array of reads file paths (FASTA or FASTQ, plain or gzipped).
   :param num_reads: Number of reads files (must be > 0).
   :param out: Pointer to receive the output struct. May be NULL if only stats
               are needed. The caller must free the output with ``smr_output_free``.
   :param stats: Pointer to receive run statistics. May be NULL.
   :returns: ``SMR_OK`` on success, or a negative error code on failure.

   The function validates all inputs before running: NULL checks, file
   existence, and empty file detection. On error, ``smr_last_error(ctx)``
   provides a descriptive message.

   A context may be reused for multiple sequential ``smr_run`` calls. Each
   call produces an independent output that must be freed separately.
   E-values are those of the CLI and match its ``aligned.blast`` report.

.. c:function:: int smr_run_seqs(smr_context_t *ctx, const char **ref_paths, int32_t num_refs, const smr_seq_t *seqs, int32_t num_seqs, smr_output_t **out, smr_stats_t *stats)

   In-memory input. Accepts sequences directly instead of file paths; no
   read files are written. Equivalent to ``smr_index_load``, one
   ``smr_run_seqs_with_index`` call and ``smr_index_free``, so the
   `Streaming`_ rules apply (one reference file, per-read e-values).

   :param ctx: Context created by ``smr_ctx_create``.
   :param ref_paths: Array of reference FASTA file paths.
   :param num_refs: Number of reference files (must be 1).
   :param seqs: Array of ``smr_seq_t`` input sequences. For paired-end mode
                (``cfg.paired = 1``), sequences must be interleaved:
                ``[fwd0, rev0, fwd1, rev1, ...]``. ``num_seqs`` must be even.
   :param num_seqs: Number of input sequences (must be > 0).
   :param out: Pointer to receive the output struct. May be NULL.
   :param stats: Pointer to receive run statistics. May be NULL.
   :returns: ``SMR_OK`` on success, or a negative error code on failure.

.. c:type:: smr_seq_t

   A single input sequence for the in-memory functions. All pointers are
   caller-owned and must remain valid for the duration of the call.

   ==================  ===============  ==========================================
   Field               Type             Description
   ==================  ===============  ==========================================
   ``id``              ``const char*``  Identifier (without ``>`` or ``@``)
   ``sequence``        ``const char*``  Nucleotide sequence
   ``quality``         ``const char*``  Quality string, or NULL for FASTA
   ==================  ===============  ==========================================

   Within one call, either every sequence has a quality string or none has.
   Identifiers need not be unique.

.. c:type:: smr_output_t

   Alignment results. Library-allocated; call ``smr_output_free`` when done.
   The output is independent of the context and remains valid after
   ``smr_ctx_destroy``.

   There is one entry per input read, in input order; paired input stays
   interleaved. Coordinates (``ref_start``, ``ref_end``) are 1-based,
   matching BLAST/SAM convention.

   ==================  ================  ===========================================
   Field               Type              Description
   ==================  ================  ===========================================
   ``num_reads``       ``uint64_t``      Total number of input reads
   ``num_aligned``     ``uint64_t``      Number of reads with at least one alignment
   ``read_ids``        ``const char**``  Read identifiers (first header token)
   ``aligned``         ``int32_t*``      1 if aligned, 0 otherwise
   ``ref_index``       ``int32_t*``      Reference index, -1 if unaligned
   ``e_value``         ``double*``       E-value of best alignment
   ``identity``        ``double*``       Percent identity (0--100)
   ``coverage``        ``double*``       Query coverage (0--100)
   ``ref_start``       ``int32_t*``      1-based start on reference
   ``ref_end``         ``int32_t*``      1-based end on reference
   ``cigar``           ``const char**``  CIGAR string, NULL if unaligned
   ``ref_name``        ``const char**``  Reference sequence ID, NULL if unaligned
   ``strand``          ``int32_t*``      1 = forward, 0 = reverse-complement; -1 if unaligned
   ``score``           ``int32_t*``      Smith-Waterman alignment score; -1 if unaligned
   ``edit_distance``   ``int32_t*``      Edit distance (mismatches + gaps); -1 if unaligned
   ==================  ================  ===========================================

.. c:type:: smr_stats_t

   Summary statistics from a run. Value-only struct -- no free required.

   =====================  ==============  ==========================================
   Field                  Type            Description
   =====================  ==============  ==========================================
   ``total_reads``        ``uint64_t``    Total reads in input
   ``total_aligned``      ``uint64_t``    Reads passing E-value threshold (any hit
                                          for the in-memory functions)
   ``total_id_cov_pass``  ``uint64_t``    Reads passing both identity and coverage
   ``total_denovo``       ``uint64_t``    De novo reads (failing ID and coverage)
   ``min_read_len``       ``uint32_t``    Shortest read length
   ``max_read_len``       ``uint32_t``    Longest read length
   ``wall_time_sec``      ``double``      Wall clock time in seconds
   =====================  ==============  ==========================================

.. c:function:: void smr_output_free(smr_output_t *out)

   Free an output struct and all library-owned memory within it.
   Passing NULL is safe (no-op).


Streaming
#########

.. c:function:: smr_index_t* smr_index_load(smr_context_t *ctx, const char **ref_paths, int32_t num_refs)

   Load (building it first if needed) the index of one reference file and
   the reference sequences into memory. Returns NULL on failure;
   ``smr_last_error(ctx)`` and ``smr_last_error_code(ctx)`` give the reason.

   ``num_refs`` must be 1: a handle holds one loaded index, and each
   reference file is a separate index. Concatenate several FASTA files into
   one. An index that does not fit in a single part is also refused with
   ``SMR_ERR_NOT_IMPLEMENTED``.

.. c:function:: int smr_run_seqs_with_index(smr_index_t *idx, const smr_seq_t *seqs, int32_t num_seqs, smr_output_t **out, smr_stats_t *stats)

   Align one in-memory batch against a loaded handle. Parameters and output
   are those of ``smr_run_seqs``. Errors and log messages go to the context
   the handle was loaded with. Batches on one handle may mix FASTA and FASTQ.

.. c:function:: void smr_index_free(smr_index_t *idx)

   Release a handle. Passing NULL is safe. The handle must not be used
   afterwards.

The handle keeps a pointer to its context, which must outlive it. Several
handles may exist at once, bound to the same or different contexts, but not
sharing one workdir. ``smr_index_load`` writes a small placeholder reads file
(``placeholder_reads.fa``) to the workdir to satisfy option parsing; it is
never read for sequence data.

E-value semantics
~~~~~~~~~~~~~~~~~

The in-memory functions report every positive Smith-Waterman hit and compute
its e-value with the per-query Karlin-Altschul form
``E = K * m * n * exp(-lambda * S)``, where ``n`` is that read's length and
``m`` is the reference length stored in the index ``.stats`` file. The CLI,
and ``smr_run``, differ in two ways:

1. **n**: the CLI uses the total length of all reads in the run. Library
   e-values are smaller by roughly the number of reads in a batch.
2. **m**: the CLI subtracts an edge-effect term that depends on the reads of
   the run. Library e-values are a further ~1-3% smaller for SILVA-scale
   references.

The CLI also drops hits below a minimum score derived from the e-value
threshold and the run's reads; the library applies no such filter and
ignores ``cfg.evalue``.

All three differences keep the per-read results independent of how reads
are grouped into batches: the same reads submitted as one batch or as many,
with any ``num_threads``, give identical per-read output. Callers filtering
on e-value should calibrate against library output, not CLI output.


Version
#######

.. c:function:: const char* smr_version(void)

   Return the SortMeRNA version string (e.g., ``"7.0.0"``).


Design notes
------------

ABI stability
#############

The ``smr_config_t`` struct uses ``struct_size`` (``uint32_t``) as its first
field for ABI version detection. ``uint32_t`` is used instead of ``size_t`` to
ensure consistent width across 32-bit and 64-bit platforms. ``smr_ctx_create``
validates that the caller's struct size matches the library's, detecting
mismatches when a caller was compiled against a different header version.

All boolean fields use ``int32_t`` (not ``bool``) for consistent sizing across
C and C++ compilers. Integer fields use explicit-width types from
``<stdint.h>``. The header requires only ``<stdint.h>`` (not ``<stddef.h>``).

Opaque context
##############

The ``smr_context_t`` type is an opaque pointer. Internal fields are hidden
from callers, allowing the library to change its internal layout without
breaking ABI. Contexts are created with ``smr_ctx_create`` and destroyed with
``smr_ctx_destroy``.

Error handling
##############

The library uses a two-level error reporting scheme:

1. **Error codes** (return values): negative integers for errors, zero for
   success. Use ``smr_strerror`` for category descriptions.
2. **Detailed messages** (per-context): ``smr_last_error`` returns a
   context-specific message with file paths, parameter values, or exception
   details from the most recent failing call.

The library never calls ``exit()``, ``abort()``, or ``assert()``. Errors,
including those raised in worker threads, are propagated to the caller via
return codes.

I/O isolation
#############

The library suppresses stdout/stderr output using two mechanisms:

1. **Thread-local log routing**: All internal logging macros (INFO, ERR, WARN)
   check a thread-local callback pointer that every API call sets for its
   duration, and that the worker threads it starts inherit. Messages go to
   the caller's ``log_callback`` under a process-wide lock, or are discarded
   when it is NULL. The ``sortmerna`` binary sets no callback and prints to
   stdout/stderr.

2. **fd-level suppression**: As a secondary defense against direct
   ``std::cout`` / ``std::cerr`` writes that bypass the macros, ``smr_run``
   and ``smr_run_seqs_with_index`` temporarily redirect file descriptors 1
   and 2 to ``/dev/null`` via ``dup2``. This is process-wide.

Log levels passed to the callback:

=================  =====
Constant           Value
=================  =====
``SMR_LOG_DEBUG``  0
``SMR_LOG_INFO``   1
``SMR_LOG_WARN``   2
``SMR_LOG_ERROR``  3
=================  =====

Thread safety
#############

``smr_run`` and ``smr_run_seqs_with_index`` are serialized by a process-level
mutex. They are safe to call from multiple threads with independent contexts
or handles -- calls will execute sequentially. ``smr_index_load`` holds the
mutex while it parses options and, if needed, builds the index; loading the
index and references into memory runs concurrently with other calls. Context creation, destruction, and error
queries on *different* contexts are thread-safe and do not acquire the
mutex. Concurrent access to the *same* context from multiple threads is not
supported.

The mutex exists because stdout/stderr suppression uses ``dup2``, which is
process-wide. The log callback routing itself uses thread-locals and is
lock-free.

Memory management
#################

- ``smr_ctx_create`` allocates the context with ``malloc``.
- ``smr_ctx_destroy`` frees it with ``free``. Passing NULL is safe.
- ``smr_run`` and the in-memory functions allocate the output struct with
  ``calloc``.
- ``smr_output_free`` frees the output and all library-owned arrays within it.
  Passing NULL is safe.
- The output is independent of the context and remains valid after
  ``smr_ctx_destroy``.
- Borrowed pointers in ``smr_config_t`` (``workdir``, ``log_user_data``) are
  not copied as strings -- the pointee must outlive the context.
