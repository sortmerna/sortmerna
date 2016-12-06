/**
 * @file mmap.cpp
 * @brief Load data (reads) using mmap.
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2012-16 Bonsai Bioinformatics Research Group
 * @copyright 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
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
 * along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 * @endparblock
 *
 * @contributors Jenya Kopylova, jenya.kopylov@gmail.com
 *               Laurent Noé, laurent.noe@lifl.fr
 *               Pierre Pericard, pierre.pericard@lifl.fr
 *               Daniel McDonald, wasade@gmail.com
 *               Mikaël Salson, mikael.salson@lifl.fr
 *               Hélène Touzet, helene.touzet@lifl.fr
 *               Rob Knight, robknight@ucsd.edu
 */

 #include "../include/mmap.hpp"


/*
 *
 * FUNCTION : mmap_reads()
 * PURPOSE  : load reads using mmap
 * OUTPUT   : double pointer array to sequences in mmap
 * See complete documentation in include/mmap.hpp
 *
 *******************************************************/
char**
mmap_reads(off_t partial_file_size,
           char* inputreads,
           off_t offset_map,
           char*& raw,
           char filesig,
           uint32_t file_s,
           uint32_t file_sections,
           int32_t &offset_pair_from_top,
           char* split_read_ptr,
           char* split_read,
           uint64_t& strs,
           char*& finalnt,
           uint32_t &reads_offset_f,
           uint32_t &reads_offset_e,
           uint32_t min_lnwin)
{
  // mmap the reads file into memory
  int fd = open(inputreads, O_RDONLY);
  raw = (char*)mmap(0, partial_file_size, PROT_READ, MAP_SHARED, fd, offset_map);
  if ( raw == MAP_FAILED )
  {
    fprintf(stderr,"  %sERROR%s: cannot mmap file: %s\n\n",startColor,"\033[0m",strerror(errno));
    exit(EXIT_FAILURE);
  }
  close(fd);
  // pointer to last character in mmap'd region
  char* end_of_mmap = &raw[partial_file_size-1];
  // range of *start and *end pointers
  {
    // (FASTA) count the number of strings in a file section
    if ( filesig == '>' )
    {
      // beginning of file section
      char *start = &raw[0];
      // end of file section
      char *end = &raw[partial_file_size-1];
#ifdef debug_mmap
      cout << "partial_file_size = " << partial_file_size << endl; //TESTING
#endif
      // compute the reads offset length at top of file section (if the first char is '>', then reads_offset_f = 0)
      while ( (*start != '>') && (start != end_of_mmap) ) { start++; reads_offset_f++; }
      if (start == end_of_mmap) reads_offset_f++;
#ifdef debug_mmap
      cout << "reads_offset_f (only split read) = " << reads_offset_f << endl; //TESTING
#endif
      // compute the reads offset length at bottom of file section (always assume the read is split, except for the last file section)
      if ( file_s != file_sections-1 )
          while ( (*end != '>') && (end != start) ) { end--; reads_offset_e++; }
#ifdef debug_mmap
      cout << "reads_offset_r (only split read) = " << reads_offset_e << endl; //TESTING
#endif
      // only the split read exists in the file section
      if ( partial_file_size-reads_offset_e-2 < 0 )
      {
        // 0 strings can only exist in the last file section
        if ( file_s != file_sections-1 )
        {
          fprintf(stderr,"   %sERROR%s: 0 sequences mapped in the current file section.\n",startColor,"\033[0m");
          exit(EXIT_FAILURE);
        }
        reads_offset_f = 0;
        reads_offset_e = 0;
        // split-read + associated paired-read
        if ( offset_pair_from_top ) strs+=2;
        // only split-read
        else strs++;
      }
      // >0 strings have been counted in the current file section
      else
      {
        // count the number of strings in the file section
        // significance of -2 is to remove the "\n>" characters from search
        // at the end of the mapped file section
        for ( int64_t i = reads_offset_f; i < partial_file_size-reads_offset_e-2; i++ )
          if ( raw[i] == '>' )
            strs++;
        // the paired-read follows the split read at the top of current file section
        if ( offset_pair_from_top )
        {
          // go forwards one character past '>'
          start++;
          reads_offset_f++;
          // increase the read's offset for the top of the file to enclose the paired-read
          while ( *start != '>' ) { start++; reads_offset_f++; }
          strs--;
        }
        // the paired-read comes before the split read at the bottom of current file section
        if ( (strs%2 == 1) && (file_s != file_sections-1) )
        {
          // go backwards one character past '>'
          end--;
          reads_offset_e++;
          // increase the read's offset for the bottom of the file to enclose the paired-read
          while ( *end != '>' ) { end--; reads_offset_e++; }
          strs--;
          // the paired-read will not follow the split read at top of the following file section
          offset_pair_from_top = false;
        }
        // the paired-read follows the split read at top of the following file section
        else offset_pair_from_top = true;         
      }         
      // account for the split read and its paired-read from the previous file section to be searched in current file section
      if ( file_s > 0 ) strs+=2;
#ifdef debug_mmap
      cout << "reads_offset_f (incl. paired-read) = " << reads_offset_f << endl; //TESTING
      cout << "reads_offset_r (incl. paired-read) = " << reads_offset_e << endl; //TESTING
      cout << "strs = " << strs << endl; //TESTING
#endif
      // one pointer for the read tag, second pointer for the read itself
      strs*=2;            
    } //~(FASTA)            
    // (FASTQ) count the number of strings in a file section
    else
    {
      // beginning of file section
      char *start = &raw[0];
      
      // end of file section
      char *end = &raw[partial_file_size-1];
#ifdef debug_mmap
      cout << "partial_file_size = " << partial_file_size << endl; //TESTING
      cout << "offset_pair_from_top (this file section) = " << offset_pair_from_top << endl; //TESTING
#endif
      // compute the reads offset length at top of file section
      int line = 0;
      // compute number of bytes covering split-read at top of file section and the paired-read (if it exists)
      if ( offset_pair_from_top > 0 )
      {
        while ( line < offset_pair_from_top )
        {
          reads_offset_f++;
          if ( *start++ == '\n' ) line++;
        }
      }
#ifdef debug_mmap
      cout << "reads_offset_f (split read + paired-read) = " << reads_offset_f << endl; //TESTING
#endif
      // count the number of lines in the file section
      for ( int64_t i = reads_offset_f; i < partial_file_size; i++ )
        if ( raw[i] == '\n' ) strs++;
#ifdef debug_mmap
      cout << "strs (not counting reads_offset_f) = " << strs << endl; //TESTING
#endif
      // check FASTQ file has multiple of 4 lines
      if ( (strs%4 != 0) && (file_s == file_sections-1) )
      {
        // odd number of lines, but due to one extra new line at bottom of file
        if ( (strs%4 == 1) && (raw[partial_file_size-1] == '\n') && (raw[partial_file_size-2] == '\n') ) strs--;
        else
        {
          fprintf(stderr,"   %sERROR%s: Your FASTQ reads file has an uneven number of lines: %lld\n",startColor,"\033[0m",strs);
          exit(EXIT_FAILURE);
        }
      }
      int offset_pair_from_bottom = 0;
      // count the number of lines belonging to a split read at bottom of file + a possible paired-read,
      // if there are only 8 strs in the file or we are on the last file, this doesn't apply
      if ( (strs > 8) && (file_s != file_sections-1) ) offset_pair_from_bottom = strs%8;
      strs-=offset_pair_from_bottom;
#ifdef debug_mmap
      cout << "offset_pair_from_bottom (this file section) = " << offset_pair_from_bottom << endl; //TESTING
#endif
      // set the number of lines to offset at the top of the next file section
      offset_pair_from_top = 8-offset_pair_from_bottom;
  #ifdef debug_mmap
      cout << "offset_pair_from_top (next file section) = " << offset_pair_from_top << endl; //TESTING
      cout << "strs (incl. reads_offset_f and reads_offset_e) = " << strs << endl; //TESTING
      if (raw[partial_file_size-1] == '\n') cout << "file section ends with a newline\n"; //TESTING
  #endif
      // count one extra newline for the split read at bottom of file                                                                                               
      // in order to facilitate reads_offset_e count                                     
      // conditions: 1. file section cannot be the last one                                                      
      //             2. file section must contain more than 0 reads                                               
      //             3. file section must not end in a new line while containing                                       
      //                an exact number of paired-reads                          
      if ( (file_s != file_sections-1) && (strs != 0) && !((raw[partial_file_size-1] == '\n') && (offset_pair_from_bottom == 0)) ) offset_pair_from_bottom++;
      // compute the reads offset length at bottom of file section
      line = 0;
      if ( offset_pair_from_bottom > 0 )
      {
        while ( line < offset_pair_from_bottom )
        {
          reads_offset_e++;
          if ( *end-- == '\n' ) line++;
        }
        // backtrack to start of read symbol
        reads_offset_e-=2;
      }
#ifdef debug_mmap
      cout << "reads_offset_r (split read + possible paired-read) = " << reads_offset_e << endl; //TESTING
#endif
      // account for the split read if exists
      if ( file_s > 0 ) strs+=8;
#ifdef debug_mmap
      cout << "strs (incl. split read and paired-read) = " << strs << endl; //TESTING
#endif
      // one pointer for the read tag, second pointer for the read itself
      strs/=2;
#ifdef debug_mmap
      cout << "strs/=2 = " << strs << endl; //TESTING
#endif
    }//~FASTQ split read          
  }//~range of *start and *end pointers
  // create a char* array for pointers to each string in mmap
  char** reads = new char*[strs]();
  if ( reads == NULL )
  {
    fprintf(stderr,"\n  %sERROR%s: cannot allocate memory for reads\n\n",startColor,"\033[0m");
    exit(EXIT_FAILURE);
  }
  // record the end of the split read
  if ( file_s > 0 )
  {
    char *start = &raw[0];
    char *end = &raw[reads_offset_f];
#ifdef debug_mmap
    cout << "split read end: \n"; //TESTING
#endif
    while ( start != end )
    {
#ifdef debug_mmap
      cout << (char)*start; //TESTING
#endif
      *split_read_ptr++ = *start++;
    }
#ifdef debug_mmap
    cout << ".STOP." <<endl; //TESTING
#endif
    // signify end of split_read
    *split_read_ptr = '\0';
    // add the split read and associated paired-read to the start of reads to be filtered
    start = split_read;      
#ifdef debug_mmap
    cout << "full split read: \n"; //TESTING
    char *tt = start;
    while ( *tt != '\0' ) cout << (char)*tt++;
    cout << ".STOP." << endl; //TESTING
#endif     
    // split-read or paired-read
    reads[0] = start;
    while ( *start++ != '\n' );
    reads[1] = start;
    int line = 0;
    if ( filesig == '@' )
    {
      // go to second read in fastq format
      while ( line < 3 ) { if ( *start++ == '\n' ) line++; }
    }
    else
    {
      // go to second read in fasta format
      while ( *start != '>') start++;
    }
    // split-read or paired-read
    reads[2] = start;
    while ( *start++ != '\n' );
    reads[3] = start;         
#ifdef debug_mmap
    cout << "*reads[0] = " << (char)*reads[0] << endl; //TESTING
    cout << "*reads[1] = " << (char)*reads[1] << endl; //TESTING
    cout << "*reads[2] = " << (char)*reads[2] << endl; //TESTING
    cout << "*reads[3] = " << (char)*reads[3] << endl; //TESTING
#endif         
  }//~if ( file_s > 0 )
  // a pointer is added to each header and the directly following read,
  // hence the number of pointers for a fasta or fastq file is the same
  char* line = &raw[reads_offset_f];
  finalnt = &raw[partial_file_size-reads_offset_e-1];
  int minlenread = READLEN;
  int readlen = 0;
#ifdef debug_mmap
  cout << "*line = " << (char)*line << endl; //TESTING
  cout << "*(finalnt-1) = " << (char)*(finalnt-1) << endl; //TESTING
  cout << "*finalnt = " << (char)*finalnt << endl; // TESTING
  if ( *finalnt != '\n' )
    cout << "*(finalnt+1) = " << (char)*(finalnt+1) << endl; //TESTING
  cout << "raw[partial_file_size-1] = " << (char)raw[partial_file_size-1] << endl; //TESTING
#endif
  int index = 0;
  if ( file_s > 0 ) index = 4;
  else index = 0;
  // (FASTA)
  if ( filesig == '>' )
  {
    while ( line != finalnt )
    {
      if ( *line == '>' )
      {
        // the read tag
        reads[index++] = line;
        while ( *line++ != '\n' );
        // the read
        reads[index++] = line;
      }
      readlen = 0;       
      while ( (line != finalnt) && (*line != '>') )
      {
        if ( *line != '\n' ) readlen++;
        line++;
      }
      // compute the minimum length read
      // if readlen in less than minlenread then minlenread = readlen,
      // otherwise minlenread = minlenread
      if ( readlen >= min_lnwin ) readlen < minlenread ? minlenread = readlen : minlenread;
    }  
  }//~if ( filesig == '>' )
  // (FASTQ)
  else
  {
    // line counter (4 lines per fastq entry)
    int count = 0;    
#ifdef debug_mmap
    cout << "index (to set up pointers) = " << index << endl; //TESTING
#endif
    while ( line != finalnt )
    {
      if ( count%4 == 0 )
      {
        // the read tag
        reads[index++] = line;
        while ( *line++ != '\n' );
        // the read
        reads[index++] = line;
        // compute the read length
        while ( *line++ != '\n' ) readlen++;
        if ( readlen >= min_lnwin ) readlen < minlenread ? minlenread = readlen : minlenread; 
        readlen = 0;     
        count+=2;           
        // skip the quality ..
      }         
      if ( *line++ == '\n' ) count++;        
    }     
#ifdef debug_mmap
    cout << "number lines indexed = " << count << endl; //TESTING
    cout << "index size = " << index-1 << endl; //TESTING
#endif
  }//~if ( filesig == '@' )   
  // debug_mmap
  if ( minlenread == READLEN )
  {
    fprintf(stderr,"   %sERROR%s: All reads are too short (<%unt) for further analysis.\n\n",startColor,"\033[0m",min_lnwin);
    exit(EXIT_FAILURE);
  }
  return reads;
}


/*
 *
 * FUNCTION : unmmap_reads()
 * PURPOSE  : load reads using mmap
 * OUTPUT   : none
 * See complete documentation in include/mmap.hpp
 *
 *******************************************************/
void unmmap_reads(char*& raw, off_t partial_file_size)
{
  // free the mmap'd file section
  if ( munmap(raw, partial_file_size ) == -1 )
  {
    fprintf(stderr,"  %sERROR%s: Could not munmap file!\n",startColor,"\033[0m");
    exit(EXIT_FAILURE);
  }
}
