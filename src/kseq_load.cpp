/*
 * @file kseq_load.cpp
 * @brief Load input reads directly into RAM (FASTA, FASTQ and compressed formats)
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2012-16 Bonsai Bioinformatics Research Group
 * @copyright 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
 *
 * This file is part of SortMeRNA.
 *
 * SortMeRNA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SortMeRNA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with SortMeRNA.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @contributors Jenya Kopylova, jenya.kopylov@gmail.com
 *               Laurent Noé, laurent.noe@lifl.fr
 *               Pierre Pericard, pierre.pericard@lifl.fr
 *               Daniel McDonald, wasade@gmail.com
 *               Mikaël Salson, mikael.salson@lifl.fr
 *               Hélène Touzet, helene.touzet@lifl.fr
 *               Rob Knight, robknight@ucsd.edu
 */

#include "../include/kseq_load.hpp"


/*
 *
 * FUNCTION : load_reads()
 * PURPOSE  : load reads using generic buffer
 * OUTPUT   : double pointer array to sequences in buffer
 * See complete documentation in include/mmap.hpp
 *
 *******************************************************/
char**
load_reads(char* inputreads,
           char*& raw,
           uint64_t number_total_read,
           off_t full_file_size,
           char*& finalnt)
{
  raw = new char[full_file_size]();
  if ( raw == NULL )
  {
    fprintf(stderr,"    %sERROR%s: [Line %d: %s] could not allocate memory for reference sequence buffer\n",
                   startColor,"\033[0m", __LINE__, __FILE__);
    exit(EXIT_FAILURE);
  }
  // 2 pointers per entry (1 for label + 1 for sequence)
  char** reads = new char*[number_total_read*2]();
  if ( reads == NULL )
  {
    fprintf(stderr,"\n  %sERROR%s: [Line %d: %s] cannot allocate memory for reads\n\n",
                   startColor,"\033[0m", __LINE__, __FILE__);
    exit(EXIT_FAILURE);
  }
#ifdef HAVE_LIBZ
  gzFile fp = gzopen(inputreads, "r");
  kseq_t *seq = kseq_init(fp);
#else
  FILE* fp = fopen(inputreads, "r");
  kseq_t *seq = kseq_init(fileno(fp));
#endif
  int l;
  int i = 0;
  char* input_ptr = raw;
  while ((l = kseq_read(seq)) >= 0) {
    // '>' or '@'
    reads[i++] = input_ptr;
    *input_ptr++ = seq->last_char;
    // sequence label
    memcpy(input_ptr, seq->name.s, seq->name.l);
    input_ptr += seq->name.l;
    // comment
    if (seq->comment.l)
    {
      *input_ptr++ = ' ';
      memcpy(input_ptr, seq->comment.s, seq->comment.l);
      input_ptr += seq->comment.l;
    }
    // '\n'
    *input_ptr++ = '\n';
    // sequence
    reads[i++] = input_ptr;
    memcpy(input_ptr, seq->seq.s, seq->seq.l);
    input_ptr += seq->seq.l;
    // '\n'
    *input_ptr++ = '\n';
    // quality (if FASTQ)
    if (seq->qual.l)
    {
      *input_ptr++ = '+';
      // '\n'
      *input_ptr++ = '\n';      
      // quality
      memcpy(input_ptr, seq->qual.s, seq->qual.l);
      input_ptr += seq->qual.l;
      // '\n'
      *input_ptr++ = '\n';      
    }
  }
  // EOF
  *input_ptr = '\0';
  finalnt = input_ptr;
  input_ptr++;
  /* TEST */
  //char* st = raw;
  //while (*st != '\0') {cout << *st++;} cout << endl;
  /* TEST */
  if (l == -2)
  {
    fprintf(stderr,"  %sERROR%s: [Line %d: %s] could not read reads file - %s\n\n",
                   startColor,"\033[0m", __LINE__, __FILE__, strerror(errno));
    exit(EXIT_FAILURE);    
  }
  kseq_destroy(seq);
#ifdef HAVE_LIBZ
  gzclose(fp);
#else
  fclose(fp);
#endif
  return reads;
}//~load_reads()
