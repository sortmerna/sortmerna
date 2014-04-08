#include <iostream>
#include <fstream>
#include <string>
#include <istream>
#include <cstdlib>
#include "../include/ssw.h"
#include "../include/common.hpp"
#include "../include/indexdb.hpp"

using namespace std;

void report_blast (ofstream &fileout,
                   s_align* a,
                   char* read_name,
                   char* read_seq,
                   char* read_qual,
                   char* ref_name,
                   char* ref_seq,
                   double evalue,
                   uint32_t readlen,
                   uint32_t bitscore,
                   bool strand,
                   double id,
                   double coverage,
                   uint32_t mismatches,
                   uint32_t gaps);


void report_sam (ofstream &fileout,
                 s_align* a,
                 char* read_name,
                 char* read_seq,
                 char* read_qual,
                 char* ref_name,
                 char* ref_seq,
                 uint32_t readlen,
                 bool strand,
                 uint32_t diff);

void report_fasta (char* acceptedstrings,
                   char* ptr_filetype_or,
                   char* ptr_filetype_ar,
                   char** reads,
                   int32_t strs,
                   vector<bool>& read_hits,
                   uint32_t file_s,
                   char* finalnt
#ifdef chimera
                   ,vector<bool>& chimeric_reads,
                   char* acceptedchimeras_file
#endif
                   );

void report_denovo(char *denovo_otus_file,
              char *ptr_filetype_or,
              char *ptr_filetype_ar,
              char **reads,
              int32_t strs,
              vector<bool>& read_hits_denovo,
              uint32_t file_s,
                   char *finalnt );

void report_biom (char* biomfile);