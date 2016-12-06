# Makefile.am

# bin are installed by default in /usr/local/bin/
bin_PROGRAMS       = indexdb_rna sortmerna

include_HEADERS    = ./alp/njn_uniform.hpp \
                     ./alp/njn_root.hpp \
                     ./alp/njn_vector.hpp \
                     ./alp/njn_approx.hpp \
                     ./alp/njn_matrix.hpp \
                     ./alp/njn_integer.hpp \
                     ./alp/njn_doubletype.hpp \
                     ./alp/njn_function.hpp \
                     ./alp/njn_memutil.hpp \
                     ./alp/njn_dynprogprob.hpp \
                     ./alp/njn_dynprogproblim.hpp \
                     ./alp/njn_dynprogprobproto.hpp \
                     ./alp/njn_ioutil.hpp \
                     ./alp/njn_localmaxstat.hpp \
                     ./alp/njn_localmaxstatmatrix.hpp \
                     ./alp/njn_localmaxstatutil.hpp \
                     ./alp/njn_random.hpp \
                     ./alp/sls_normal_distr_array.hpp \
                     ./alp/sls_alignment_evaluer.hpp \
                     ./alp/sls_alp.hpp \
                     ./alp/sls_alp_data.hpp \
                     ./alp/sls_alp_regression.hpp \
                     ./alp/sls_alp_sim.hpp \
                     ./alp/sls_basic.hpp \
                     ./alp/sls_pvalues.hpp \
                     ./cmph/bdz.h \
                     ./cmph/bdz_ph.h \
                     ./cmph/bdz_structs.h \
                     ./cmph/bdz_structs_ph.h \
                     ./cmph/bitbool.h \
                     ./cmph/bmz.h \
                     ./cmph/bmz8.h \
                     ./cmph/bmz8_structs.h \
                     ./cmph/bmz_structs.h \
                     ./cmph/brz.h \
                     ./cmph/brz_structs.h \
                     ./cmph/buffer_entry.h \
                     ./cmph/buffer_manager.h \
                     ./cmph/chd.h \
                     ./cmph/chd_ph.h \
                     ./cmph/chd_structs.h \
                     ./cmph/chd_structs_ph.h \
                     ./cmph/chm.h \
                     ./cmph/chm_structs.h \
                     ./cmph/cmph.h \
                     ./cmph/cmph_structs.h \
                     ./cmph/cmph_time.h \
                     ./cmph/cmph_types.h \
                     ./cmph/compressed_rank.h \
                     ./cmph/compressed_seq.h \
                     ./cmph/debug.h \
                     ./cmph/fch.h \
                     ./cmph/fch_buckets.h \
                     ./cmph/fch_structs.h \
                     ./cmph/graph.h \
                     ./cmph/hash.h \
                     ./cmph/hash_state.h \
                     ./cmph/jenkins_hash.h \
                     ./cmph/miller_rabin.h \
                     ./cmph/select.h \
                     ./cmph/select_lookup_tables.h \
                     ./cmph/vqueue.h \
                     ./cmph/vstack.h \
                     ./include/config.h \
                     ./include/common.hpp \
                     ./include/indexdb.hpp \
                     ./include/bitvector.hpp \
                     ./include/load_index.hpp \
                     ./include/traverse_bursttrie.hpp \
                     ./include/kseq.h \
                     ./include/kseq_load.hpp \
                     ./include/mmap.hpp \
                     ./include/alignment.hpp \
                     ./include/paralleltraversal.hpp \
                     ./include/ssw.h \
                     ./include/outputformats.hpp

indexdb_rna_SOURCES    = ./cmph/bdz.c \
                         ./cmph/bdz_ph.c \
                         ./cmph/bmz.c \
                         ./cmph/bmz8.c \
                         ./cmph/brz.c \
                         ./cmph/buffer_entry.c \
                         ./cmph/buffer_manager.c \
                         ./cmph/chd.c \
                         ./cmph/chd_ph.c \
                         ./cmph/chm.c \
                         ./cmph/cmph.c \
                         ./cmph/cmph_structs.c \
                         ./cmph/compressed_rank.c \
                         ./cmph/compressed_seq.c \
                         ./cmph/fch.c \
                         ./cmph/fch_buckets.c \
                         ./cmph/graph.c \
                         ./cmph/hash.c \
                         ./cmph/jenkins_hash.c \
                         ./cmph/miller_rabin.c \
                         ./cmph/select.c \
                         ./cmph/vqueue.c \
                         ./cmph/vstack.c \
                         ./src/indexdb.cpp

sortmerna_SOURCES  = ./alp/njn_dynprogprob.cpp \
                     ./alp/njn_dynprogproblim.cpp \
                     ./alp/njn_dynprogprobproto.cpp \
                     ./alp/njn_ioutil.cpp \
                     ./alp/njn_localmaxstat.cpp \
                     ./alp/njn_localmaxstatmatrix.cpp \
                     ./alp/njn_localmaxstatutil.cpp \
                     ./alp/njn_random.cpp \
                     ./alp/sls_alignment_evaluer.cpp \
                     ./alp/sls_alp.cpp \
                     ./alp/sls_alp_data.cpp \
                     ./alp/sls_alp_regression.cpp \
                     ./alp/sls_alp_sim.cpp \
                     ./alp/sls_basic.cpp \
                     ./alp/sls_pvalues.cpp \
                     ./src/bitvector.cpp \
                     ./src/load_index.cpp \
                     ./src/traverse_bursttrie.cpp \
                     ./src/kseq_load.cpp \
                     ./src/mmap.cpp \
                     ./src/alignment.cpp \
                     ./src/paralleltraversal.cpp \
                     ./src/ssw.c \
                     ./src/outputformats.cpp \
                     ./src/main.cpp

indexdir       =  $(prefix)/@PACKAGE@/rRNA_databases
index_DATA     =  rRNA_databases/rfam-5.8s-database-id98.fasta \
                  rRNA_databases/rfam-5s-database-id98.fasta \
                  rRNA_databases/silva-arc-23s-id98.fasta \
                  rRNA_databases/silva-euk-28s-id98.fasta \
                  rRNA_databases/silva-arc-16s-id95.fasta \
                  rRNA_databases/silva-bac-16s-id90.fasta \
                  rRNA_databases/silva-bac-23s-id98.fasta \
                  rRNA_databases/silva-euk-18s-id95.fasta 

EXTRA_DIST = rRNA_databases/rfam-5.8s-database-id98.fasta \
             rRNA_databases/rfam-5s-database-id98.fasta \
             rRNA_databases/silva-arc-23s-id98.fasta \
             rRNA_databases/silva-arc-16s-id95.fasta \
             rRNA_databases/silva-bac-16s-id90.fasta \
             rRNA_databases/silva-bac-23s-id98.fasta \
             rRNA_databases/silva-euk-18s-id95.fasta \
             rRNA_databases/README.txt \
             rRNA_databases/scripts/edit_U_to_T_rna.py \
             rRNA_databases/scripts/generate_databases.sh \
             rRNA_databases/scripts/hmmer_clean.sh \
             rRNA_databases/silva-euk-28s-id98.fasta \
             scripts/merge-paired-reads.sh \
             scripts/unmerge-paired-reads.sh \
             m4/ax_check_zlib.m4 \
             index \
             SortMeRNA-User-Manual-2.1.pdf 

uninstall: uninstall-am

