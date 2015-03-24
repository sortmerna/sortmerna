/**
 * @file paralleltraversal.cpp
 * @brief File containing functions for index traversal.
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright Copyright (C) 2012-2014 Bonsai Bioinformatics Research Group, LIFL and 
 * INRIA Nord-Europe, France
 * OTU-picking extensions developed in the Knight Lab, BioFrontiers Institute,
 * University of Colorado at Boulder, Boulder, CO
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
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 * @endparblock
 *
 * @authors jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 */

#include "../include/paralleltraversal.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

// see "heuristic 1" below
//#define HEURISTIC1_OFF

/*! @brief Euler's constant */
#define EXP 2.71828182845904523536


/*! @brief The universal Levenshtein automaton for d=1.

    The maximum length of a characteristic bitvector for d=1 is 2d+2=4.
    The number of states of the this automaton is 15, these are numbered
    as 0-14 and correspond to the following symbolic setup:\n
    <table>
     <tr><th> Table index </th><th> Symbolic state </th><th> Significance </th></tr>
     <tr><td> 0 </td><td> {I^0} </td><td> non-accepting state </td></tr>
     <tr><td> 1 </td><td> {(I-1)^1} </td><td> non-accepting state </td></tr>   
     <tr><td> 2 </td><td> {I^1} </td><td> non-accepting state </td></tr>  
     <tr><td> 3 </td><td> {(I-1)^1, I^1}} </td><td> non-accepting state </td></tr>  
     <tr><td> 4 </td><td> {(I+1)^1}  </td><td> non-accepting state </td></tr> 
     <tr><td> 5 </td><td> {(I-1)^1, (I+1)^1} </td><td> non-accepting state </td></tr> 
     <tr><td> 6 </td><td> {I^1, (I+1)^1}  </td><td> non-accepting state </td></tr> 
     <tr><td> 7 </td><td> {(I-1)^1, I^1, (I+1)^1}  </td><td> non-accepting state </td></tr> 
     <tr><td> 8 </td><td> {(M-1)^0} </td><td> accepting state </td></tr> 
     <tr><td> 9 </td><td> {(M-1)^0} </td><td> accepting state </td></tr> 
     <tr><td> 10 </td><td> {M^1} </td><td> accepting state </td></tr> 
     <tr><td> 11 </td><td> {(M-2)^1, M^1} </td><td> accepting state </td></tr>
     <tr><td> 12 </td><td> {(M-1)^1, M^1} </td><td> accepting state </td></tr>
     <tr><td> 13 </td><td> {(M-2)^1, (M-1)^1, M^1} </td><td> accepting state </td></tr>
     <tr><td> 14 </td><td> NULL </td><td> failure state </td></tr>
    </table>
 */
uint32_t table[4][16][14] = {
    {{3, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14},
    {3, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14},
    {7, 14, 14, 14, 4,  4,  4,  4,  14, 14, 14, 14, 14, 14},
    {7, 14, 14, 14, 4,  4,  4,  4,  14, 14, 14, 14, 14, 14},
    {0, 14, 2,  2,  14, 14, 2,  2,  14, 14, 14, 14, 14, 14},
    {0, 14, 2,  2,  14, 14, 2,  2,  14, 14, 14, 14, 14, 14},
    {0, 14, 2,  2,  4,  4,  6,  6,  14, 14, 14, 14, 14, 14},
    {0, 14, 2,  2,  4,  4,  6,  6,  14, 14, 14, 14, 14, 14},
    {3, 1,  14, 1, 14,  1, 14,  1,  14, 14, 14, 14, 14, 14},
    {3, 1,  14, 1, 14,  1, 14,  1,  14, 14, 14, 14, 14, 14},
    {7, 1,  14, 1, 4,   5, 4,   5,  14, 14, 14, 14, 14, 14},
    {7, 1,  14, 1, 4,   5, 4,   5,  14, 14, 14, 14, 14, 14},
    {0, 1,  2,  3, 14,  1, 2,   3,  14, 14, 14, 14, 14, 14},
    {0, 1,  2,  3, 14,  1, 2,   3,  14, 14, 14, 14, 14, 14},
    {0, 1,  2,  3, 4,   5, 6,   7,  14, 14, 14, 14, 14, 14},
    {0, 1,  2,  3, 4,   5, 6,   7,  14, 14, 14, 14, 14, 14}},
    {{3, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14},
    {13, 14, 14, 14, 10, 10, 10, 10, 14, 14, 14, 14, 14, 14},
    {8, 14, 2,  2,  14, 14, 2,  2,  14, 14, 14, 14, 14, 14},
    {8, 14, 2,  2,  10, 10, 12, 12, 14, 14, 14, 14, 14, 14},
    {3, 1, 14,  1,  14, 1,  14, 1,  14, 14, 14, 14, 14, 14},
    {13, 1, 14, 1,  10, 11, 10, 11, 14, 14, 14, 14, 14, 14},
    {8,  1,  2,  3,  14, 1,  2,  3, 14, 14, 14, 14, 14, 14},
    {8, 1, 2, 3, 10, 11, 12, 13, 14, 14, 14, 14, 14, 14}},
    {{12, 14, 14, 14, 14, 14, 14, 14, 12, 14, 14, 14, 14, 14},
    {9, 14, 10, 10, 14, 14, 10, 10, 9, 14, 14, 14, 10, 10},
    {12, 1, 14, 1, 14, 1, 14, 1, 12, 14, 14, 1, 14, 1},
    {9, 1, 10, 12, 14, 1, 10, 12, 9, 14, 14, 1, 10, 12}},
    {{10, 14, 14, 14, 14, 14, 14, 14, 14, 10, 14, 14, 14, 14},
    {10, 10, 14, 10, 14, 10, 14, 10, 14, 10, 14, 14, 10, 14}}};

/*! @brief Map nucleotides to integers.

    Ambiguous letters map to 4. 
    {A/a,C/c,G/g,T/t,U/u} = {0,1,2,3,3} respectively. 
*/
char nt_table[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

/*! @brief Return complement of a nucleotide in
    integer format.

    <table>
     <tr><th>i</th> <th>complement[i]</th></tr>
     <tr><td>0 (A)</td> <td>3 (T)</td></tr>
     <tr><td>1 (C)</td> <td>2 (G)</td></tr>
     <tr><td>2 (G)</td> <td>1 (C)</td></tr>
     <tr><td>3 (T)</td> <td>0 (A)</td></tr>
    </table>
 */
char complement[4] = {3,2,1,0};


/*! @fn smallest()
    @brief Determine the smallest integer.
    @details The mypair data structure holds two integers,
             the first being the position a k-mer occurs 
             on the reference sequence and the second
             being the position a k-mer occurs on the
             query sequence. This function takes two
             mypair data structures and returns the 

    @param const pair<uint32_t,uint32_t> &a
    @param const pair<uint32_t,uint32_t> &b
    @return smallest integer of a and b, or a if a == b
*/
bool smallest ( const mypair &a, const mypair &b )
{
  if ( a.first == b.first ) return ( a.second < b.second );
  return ( a.first < b.first );
}

/*! @fn largest()
    @brief Return the largest integer of two input integers
    @param const mypair &a
    @param const mypair &b
    @return 'a' goes before 'b' if a.first > b.first, otherwise
            'a'
*/
bool largest ( const mypair &a, const mypair &b )
{
  if ( a.first == b.first ) return ( a.second > b.second );
  return ( a.first > b.first );
}

/* @function format_forward()
 * format the forward read into a string on same alphabet without '\n'
 *
 * */
void format_forward(char* read_seq,char* myread,char filesig)
{
  // FASTA
  if ( filesig == '>' )
  {
    while ( (*read_seq != '\0') && (*read_seq != '>') )
    {
      if (*read_seq != '\n') *myread++ = nt_table[(int)*read_seq];
      read_seq++;
    }
    *myread='\n'; // end of read marked by newline
  }
  // FASTQ
  else
  {
    while ( *read_seq != '\n' ) { *myread++ = nt_table[(int)*read_seq++]; }
    *myread='\n'; //end of read marked by newline
  }
}

/* @function format_rev()
 * format the reverse-complement read into a string without '\n'
 *
 * */
void format_rev(char* start_read,char* end_read,char* myread,char filesig)
{
  int8_t rc_table[128] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
  };
    
  // FASTA
  if ( filesig == '>' )
  {
    while ( end_read != start_read )
    {
      if (*end_read != '\n') *myread++ = rc_table[(int)*end_read];
      end_read--;
    }
    *myread++ = rc_table[(int)*end_read];
    *myread='\n';
  }
  // FASTQ
  else
  {
    while ( *end_read != '\n' ) { *myread++ = rc_table[(int)*end_read--]; }
    *myread='\n';
  }
}


/*! @fn traversetrie_align() */
inline void
traversetrie_align ( NodeElement *trie_t,
                    uint32_t lev_t,
                    unsigned char depth,
                    MYBITSET *win_k1_ptr,
                    MYBITSET *win_k1_full,
                    bool &accept_zero_kmer,
                    vector< id_win > &id_hits,
                    uint32_t readn,
                    uint32_t win_num,
                    uint32_t partialwin
                    )
{
  uint16_t lev_t_trie_pivot = lev_t;
  unsigned char value = 0;
    
  // traverse through the node elements in a trie node
  for ( uint32_t node_element = 0; node_element < 4; node_element++ )
  {
    value = trie_t->flag;
        
    // this node element is empty, go to next node element in trie node
    if ( value == 0 )
    {
      lev_t = lev_t_trie_pivot;
      trie_t++;
    }
        
    // this node element points to a trie node or a bucket, continue traversing
    else
    {
      if ( depth < partialwin-2 )
      {
        // send bv to LEV(1)
        lev_t = table[0][(int)*(win_k1_ptr + (depth<<2) + node_element)][(int)(lev_t)];
      }
      else
      {
        lev_t = table[3-partialwin+depth][(int)(*(win_k1_full+node_element) & ((2<<(partialwin-depth))-1))][(int)(lev_t)];
      }
            
            
      // LEV(1) is in a null state, go to next node element
      if ( lev_t == 14 )
      {
        lev_t = lev_t_trie_pivot;
        trie_t++;
      }
      // LEV(1) is not in a null state, continue parallel traversal
      else
      {
        // (1) the node element holds a pointer to another trie node
        if ( value == 1 )
        {
          traversetrie_align ( trie_t->whichnode.trie,
                                        lev_t,
                                        ++depth,
                                        win_k1_ptr,
                                        win_k1_full,
                                        accept_zero_kmer,
                                        id_hits,
                                        readn,
                                        win_num,
                                        partialwin
                                        );
                    
          // go to next window on the read (0-error match found)
          if ( accept_zero_kmer ) return ;
                    
          --depth;
                    
          // go to the next trie node element
          lev_t = lev_t_trie_pivot;
          trie_t++;
        }//~ if (child node a trie node)                
        // (2) the node element points to a bucket
        else
        {
          // this pivot remembers the target levenshtein state of the terminal trie node, from
          //  which exists a bucket; every element in the bucket takes this lev_t state as an
          //  initial input state
          uint32_t lev_t_bucket_pivot = lev_t;
                    
          // number of characters per entry
          uint32_t s = partialwin-depth;
                    
          unsigned char* start_bucket = (unsigned char*)trie_t->whichnode.bucket;
          if ( start_bucket == NULL )
          {
            fprintf(stderr, "  ERROR: pointer start_bucket == NULL (paralleltraversal.cpp)\n");
            exit(EXIT_FAILURE);
          }
          unsigned char* end_bucket = start_bucket + trie_t->size;
          if ( end_bucket == NULL )
          {
            fprintf(stderr, "  ERROR: pointer end_bucket == NULL (paralleltraversal.cpp)\n");
            exit(EXIT_FAILURE);
          }
                    
          // traverse the bucket
          while ( start_bucket != end_bucket )
          {
            uint32_t depth_b = depth;
            lev_t = lev_t_bucket_pivot;
            bool local_accept_kmer = false;
            uint32_t entry_str = *((uint32_t*)start_bucket);
                        
            // for each nt in the string
            for ( uint32_t j = 0; j < s; j++ )
            {
              uint32_t nt = entry_str&3;
                            
              depth_b++;
                            
              // get bitvector for letter
              if ( depth_b < partialwin-2 )
              {
                // send bv to LEV(_k)
                lev_t = table[0][(int)*(win_k1_ptr + (depth_b<<2) + nt)][(int)(lev_t)];
              }
              else
              {
                lev_t = table[3-partialwin+depth_b][(int)(*(win_k1_full + nt) & ((2<<(partialwin-depth_b))-1) )][(int)(lev_t)];
              }
                            
              // if the target lev_t state is a failure state, go to the next bucket element (tail)
              if ( lev_t == 14 ) break;
                            
              // approaching end of tail
              if ( depth_b >= partialwin-2 )
              {
                // 1-error match
                if ( lev_t >= 8 )
                {
                  local_accept_kmer = true;
                }
                // 0-error match
                if ( depth_b == partialwin-1 )
                {
                  if ( lev_t == 9 )
                  {
                    accept_zero_kmer = true;
                                        
                    // turn off heuristic to stop search after finding 0-error match
                    if ( full_search_gv ) accept_zero_kmer = false;
                  }
                }
              }//~last 3 characters in entry
                            
              if ( local_accept_kmer )
              {
                id_win entry = {0,0};
                entry.id = *((uint32_t*)start_bucket + 1);
                entry.win = win_num;
                                
                // empty id_hits array, add 0-error id and exit
                if ( accept_zero_kmer )
                {
                  id_hits.clear();
                  id_hits.push_back(entry);
                                    
                  return;
                }
                                
                // exact match not found, do not include duplicates of 1-error match (for the same window on read)
                if ( !id_hits.empty() )
                {
                  bool found = false;
                  for ( int f = 0; f < id_hits.size(); f++ )
                  {
                    if ( id_hits[f].id == entry.id )
                    {
                      found = true;
                      break;
                    }
                  }
                  if ( found ) break;
                }
                                
                id_hits.push_back(entry);
                                
              }                         
              entry_str>>=2;
            }//~for each 2 bits
                        
            // next entry
            start_bucket+=ENTRYSIZE;
          }//~for each entry
                    
          lev_t = lev_t_trie_pivot;
          trie_t++;
                    
        }//~else the node element points to a bucket
      }//~else LEV(1) is not in a null state, continue parallel traversal
    }//~else this node element points to a trie node or a bucket, continue traversing
  }//~for all trie nodes
    
  return ;
}//~traversetrie_align()




/*
 *
 * @function check files: load the @SQ lines for SAM output & input all header information to SAM output file
 * @return void
 * @version 2.0 February 4, 2014
 *
 *******************************************************************/
void preprocess_data(vector< pair<string,string> >& myfiles,
                     char** argv,
                     int argc,
                     bool yes_SQ,
                     char* acceptedstrings_sam,
                     int32_t _match,
                     int32_t _mismatch,
                     int32_t _gap_open,
                     int32_t _gap_extension,
                     vector<vector<uint32_t> >& skiplengths,
                     vector<uint16_t>& num_index_parts,
                     vector<vector<index_parts_stats> >& index_parts_stats_vec,
                     vector<uint64_t>& full_ref,
                     vector<uint64_t>& full_read,
                     vector<uint32_t>& lnwin,
                     vector<uint32_t>& partialwin,
                     vector<uint32_t>& minimal_score,
                     uint32_t number_total_read,
                     vector<pair<double, double> >& gumbel,
                     vector<uint32_t>& numbvs,
                     vector<uint32_t>& numseq)
{
  ofstream acceptedsam;
  
  if ( samout_gv )
  {
    acceptedsam.open (acceptedstrings_sam);
    if (!acceptedsam.good())
    {
      fprintf(stderr,"  %sERROR%s: could not open SAM output file for writing.\n","\033[0;31m","\033[0m");
      exit(EXIT_FAILURE);
    }
    // @HD header
    else acceptedsam << "@HD\tVN:1.0\tSO:unsorted\n";
  }
    
  // loop through the .seqs index files for each database
  for ( uint16_t index_num = 0; index_num < (uint16_t)myfiles.size(); index_num++ )
  {
      ifstream stats( (char*)(myfiles[index_num].second + ".stats").c_str(), ios::in | ios::binary );
      if ( !stats.good() )
      {
        fprintf(stderr,"\n  %sERROR%s: The index '%s' does not exist.\n","\033[0;31m","\033[0m",(char*)(myfiles[index_num].second + ".stats").c_str());
        fprintf(stderr,"  Make sure you have constructed your index using the command `indexdb'. See `indexdb -h' for help.\n\n");
        exit(EXIT_FAILURE);
      }
      
      // read the file size for file used to build the index
      size_t filesize = 0;
      stats.read(reinterpret_cast<char*>(&filesize), sizeof(size_t));
      
      // read the fasta file name used to build the index
      uint32_t fastafile_len = 0;
      stats.read(reinterpret_cast<char*>(&fastafile_len), sizeof(uint32_t));
      
      char fastafile_name[2000];
      stats.read(reinterpret_cast<char*>(fastafile_name), sizeof(char)*fastafile_len);
      
      // compute reference database file size for this index
      FILE *fastafile = fopen ((char*)(myfiles[index_num].first).c_str(),"r");
      if ( fastafile == NULL )
      {
        fprintf(stderr,"    %sERROR%s: could not open FASTA reference file: %s .\n","\033[0;31m","\033[0m",(char*)(myfiles[index_num].first).c_str());
        exit(EXIT_FAILURE);
      }
      fseek(fastafile,0L,SEEK_END);
      size_t sz = ftell(fastafile);
      fclose(fastafile);
        
      if ( sz != filesize )
      {
        fprintf(stderr,"    %sERROR%s: Based on file size, the FASTA file (%s) passed to --ref <FASTA file, index name>\n","\033[0;31m","\033[0m",(char*)(myfiles[index_num].first).c_str());
        fprintf(stderr,"    does not appear to be the same FASTA file (%s) used to build the index %s.\n",fastafile_name,(char*)(myfiles[index_num].second).c_str());
        fprintf(stderr,"    Check your --ref list of files and corresponding indexes.\n\n");
        exit(EXIT_FAILURE);
      }
      
      // A,C,G,T background frequencies to compute the Gumbel parameters lambda and K
      double background_freq_gv[4] = {0};
      
      // A/C/G/T distribution frequencies
      stats.read(reinterpret_cast<char*>(&background_freq_gv), sizeof(double)*4);
      
      // total length of sequences in the complete database
      stats.read(reinterpret_cast<char*>(&full_ref[index_num]), sizeof(uint64_t));
      
      // sliding window length lnwin & initialize
      stats.read(reinterpret_cast<char*>(&lnwin[index_num]), sizeof(uint32_t));
      
      // total number of reference sequences in one complete reference database
      stats.read(reinterpret_cast<char*>(&numseq[index_num]), sizeof(uint32_t));
      
      partialwin[index_num] = lnwin[index_num]/2;
        
      // number of bitvectors at depth > 0 in [w_1] reverse or [w_2] forward
      numbvs[index_num] = 4*(partialwin[index_num]-3);
      
      // set the window shift for different seed lengths (if not set by user, or one of the lengths is <= 0)
      if ( (skiplengths[index_num][0] == 0) || (skiplengths[index_num][1] == 0) || (skiplengths[index_num][2] == 0) )
      {
        skiplengths[index_num][0] = lnwin[index_num];
        skiplengths[index_num][1] = partialwin[index_num];
        skiplengths[index_num][2] = 3;
      }
      
      // number of index parts
      stats.read(reinterpret_cast<char*>(&num_index_parts[index_num]), sizeof(uint16_t));
      vector<index_parts_stats> hold;
      
      // information on the location and size of sequences used to build each index part
      for ( uint16_t j = 0; j < num_index_parts[index_num]; j++ )
      {
        index_parts_stats stats_hold;
        stats.read(reinterpret_cast<char*>(&stats_hold), sizeof(index_parts_stats));
        hold.push_back(stats_hold);
      }
        
      index_parts_stats_vec.push_back(hold);
      
      // compute Gumbel parameters
      long int rand_ = 182345345;
      //string randout_= "./alp/random_param.txt";
      string randout_= "";
      
      int gapopen_ = _gap_open;
      int gapopen1_ = _gap_open;
      int gapopen2_ = _gap_open;
      
      int gapextend_ = _gap_extension;
      int gapextend1_ = _gap_extension;
      int gapextend2_ = _gap_extension;
      
      int match = _match;
      int mismatch = _mismatch;
      double A_ = background_freq_gv[0];
      double C_ = background_freq_gv[1];
      double G_ = background_freq_gv[2];
      double T_ = background_freq_gv[3];
      
      string scoremat_file_name_ ="";
      string freqs1_file_name_ ="";
      string freqs2_file_name_ ="";
      double max_time_=1;
      double max_mem_=500;
      double eps_lambda_gv_=0.001;
      double eps_K_gv_=0.005;
      string gumbelparout_file_name_ ="";
      bool gapped_ = true;
      bool insertions_after_deletions_=false;
      
      Sls::set_of_parameters gumbel_params;
        
      CGumbelParamsCalc::Params_Run2(
                                     rand_,//randomization number
                                     randout_,//if true, then the program outputs complete randomization information into a file
                                     
                                     gapopen_,//gap opening penalty
                                     gapopen1_,//gap opening penalty for a gap in the sequence #1
                                     gapopen2_,//gap opening penalty for a gap in the sequence #2
                                     
                                     gapextend_,//gap extension penalty
                                     gapextend1_,//gap extension penalty for a gap in the sequence #1
                                     gapextend2_,//gap extension penalty for a gap in the sequence #2
                                     
                                     scoremat_file_name_,//scoring matrix file name
                                     freqs1_file_name_,//probabilities1 file name
                                     freqs2_file_name_,//probabilities1 file name
                                     max_time_,//maximum allowed calculation time in seconds
                                     max_mem_,//maximum allowed memory usage in MB
                                     eps_lambda_gv_,//relative error for lambda_gv calculation
                                     eps_K_gv_,//relative error for K_gv calculation
                                     gumbelparout_file_name_,
                                     gapped_,
                                     insertions_after_deletions_,//if true, then insertions after deletions are allowed
                                     gumbel_params,
                                     match,//NEW - SW score for a match,
                                     mismatch,//NEW - SW score for a mismatch,
                                     A_,//NEW - background frequency for A
                                     C_,//NEW - background frequency for C
                                     G_,//NEW - background frequency for G
                                     T_,//NEW- background frequency for T
                                     (gumbel[index_num].first), //lambda
                                     (gumbel[index_num].second) // K
                                     );
        
        
      // Shannon's entropy for reference sequence nucleotide distribution
      double entropy_H_gv = -(background_freq_gv[0]*(log(background_freq_gv[0])/log(2)) +
                              background_freq_gv[1]*(log(background_freq_gv[1])/log(2)) +
                              background_freq_gv[2]*(log(background_freq_gv[2])/log(2)) +
                              background_freq_gv[3]*(log(background_freq_gv[3])/log(2)));
      
      // Length correction for Smith-Waterman alignment score
      uint64_t expect_L = log((gumbel[index_num].second)*full_read[index_num]*full_ref[index_num])/entropy_H_gv;
      
      // correct the reads & databases sizes for e-value calculation
      if ( full_ref[index_num] > (expect_L*numseq[index_num]) ) full_ref[index_num]-=(expect_L*numseq[index_num]);
      full_read[index_num]-=(expect_L*number_total_read);
      
      // minimum score required to reach E-value
      minimal_score[index_num] = (log(evalue/((double)(gumbel[index_num].second)*full_ref[index_num]*full_read[index_num])))/-(gumbel[index_num].first);
        
      // SAM @SQ data
      if ( samout_gv )
      {
        // number of nucleotide sequences in the reference file
        uint32_t num_sq = 0;
        stats.read(reinterpret_cast<char*>(&num_sq), sizeof(uint32_t));
        
        // loop through each @SQ
        for ( uint32_t j = 0; j < num_sq; j++ )
        {
          // length of the sequence id
          uint32_t len_id = 0;
          stats.read(reinterpret_cast<char*>(&len_id), sizeof(uint32_t));
          
          // the sequence id string
          char s[len_id+1];
          memset(s,0,len_id+1);
          stats.read(reinterpret_cast<char*>(&s), sizeof(char)*len_id);
          
          // the length of the sequence itself
          uint32_t len_seq = 0;
          stats.read(reinterpret_cast<char*>(&len_seq), sizeof(uint32_t));
          
          // @SQ header
          if ( yes_SQ ) acceptedsam << "@SQ\tSN:" << s << "\tLN:" << len_seq << "\n";
        }
      }
        
      stats.close();
    }//~for all indexes
    
    if ( samout_gv )
    {
      // @PG to sam file
      acceptedsam << "@PG\tID:sortmerna\tVN:1.0\tCL:";
      for ( int j = 0; j < argc; j++ ) acceptedsam << argv[j] << " ";
      acceptedsam << endl;
      
      acceptedsam.close();
    }//~if samout_gv
}//~preprocess_data()





/*
 *
 * @function load_ref: load fasta reference sequences into memory
 * for SW alignment
 * @return void
 * @version 1.0 June 12, 2013
 *
 *******************************************************************/
void
load_ref(char* ptr_dbfile,
         char* buffer,
         char** reference_seq,
         uint32_t* reference_seq_len,
         uint64_t seq_part_size,
         uint32_t numseq_part,
         uint64_t start_part,
         bool load_for_search)
{
  FILE *fp = fopen(ptr_dbfile,"r");
  if ( fp == NULL )
  {
    fprintf(stderr,"  %sERROR%s: could not open file %s\n",ptr_dbfile,"\033[0;31m","\033[0m");
    exit(EXIT_FAILURE);
  }
    
  // set the file pointer to the first sequence added to the index for this index file section
  if ( fseek(fp,start_part,SEEK_SET) != 0 )
    {
      fprintf(stderr,"  %sERROR%s: could not locate the sequences used to construct the index (paralleltraversal.cpp).\n","\033[0;31m","\033[0m");
      fprintf(stderr,"  Check that your --ref <FASTA file, index name> correspond correctly for the FASTA file: %s.\n",ptr_dbfile);
    }
    
  // load references sequences into memory, skipping the new lines & spaces in the fasta format
  uint32_t num_seq_read = 0;
  char *s = buffer;
    
  int i = 0;
  int j = 0;
  char c = fgetc(fp);

  if ( load_for_search )
  {
    do
    {
      // the tag
      reference_seq[i++] = s;
      while ( c != '\n' )
      {
        *s++ = c;
        c = fgetc(fp);
      }
      // new line
      *s++ = c;
      
      if ( *s == '\n' )
      {
        fprintf(stderr,"  %sERROR%s: your reference sequences are not in FASTA format "
                       "(there is an extra new line).","\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
      }
      
      // the sequence
      reference_seq[i++] = s;
      c = fgetc(fp);
      do
      {
        if ( c != '\n' && c != ' ' )
        {
          // keep record of ambiguous character for alignment
          *s++ = nt_table[(int)c];
          
          // record the sequence length as we read it
          reference_seq_len[j]++;
        }
        c = fgetc(fp);
      } while ( (c != '>') && (c != EOF) );
      
      *s++ = '\n';
      j++;
      
      num_seq_read++;
        
    } while ( (num_seq_read != numseq_part) && (c != EOF) );
  }
  else
  {
    do
    {
      // the tag
      reference_seq[i++] = s;
      while ( c != '\n' )
      {
        *s++ = c;
        c = fgetc(fp);
      }
      
      // new line
      *s++ = c;
      
      // the sequence
      reference_seq[i++] = s;
      c = fgetc(fp);
      do
      {
        if ( c != '\n' && c != ' ' )
        {
            // keep record of ambiguous character for alignment
            *s++ = nt_table[(int)c];
        }
        c = fgetc(fp);
      } while ( (c != '>') && (c != EOF) );
      
      *s++ = '\n';
      j++;
      
      num_seq_read++;
        
    } while ( (num_seq_read != numseq_part) && (c != EOF) );
  }
    
  fclose(fp);
}//~load_ref


#ifdef see_binary_output
/*
 *
 * @function traversetrie: collect statistics on the mini-burst trie,
 * its size, number of trie nodes vs. buckets
 * @param NodeElement* trie_node
 * @return void
 * @version 1.0 Jan 14, 2013
 *
 *******************************************************************/
void traversetrie_debug( NodeElement* trie_node, uint32_t depth, uint32_t &total_entries, string &kmer_keep, uint32_t partialwin )
{
  char get_char[4] = {'A','C','G','T'};
    
  // traverse through the node elements in a trie node
  for ( int i = 0; i < 4; i++ )
  {
    unsigned char value = trie_node->flag;
        
    // the node element holds a pointer to another trie node
    if ( value == 1 )
    {
      kmer_keep.push_back((char)get_char[i]); //TESTING
      traversetrie_debug( trie_node->whichnode.trie, ++depth, total_entries, kmer_keep, partialwin );
      kmer_keep.pop_back();
      --depth;
    }
        
    // the node element points to a bucket
    else if ( value == 2 )
    {
      kmer_keep.push_back((char)get_char[i]); //TESTING
      
      unsigned char* start_bucket = (unsigned char*)trie_node->whichnode.bucket;
      if ( start_bucket == NULL )
      {
        fprintf(stderr, "  ERROR: pointer start_bucket == NULL (paralleltraversal.cpp)\n");
        exit(EXIT_FAILURE);
      }
      unsigned char* end_bucket = start_bucket + trie_node->size;
      if ( end_bucket == NULL )
      {
        fprintf(stderr, "  ERROR: pointer end_bucket == NULL (paralleltraversal.cpp)\n");
        exit(EXIT_FAILURE);
      }
      
      cout << "size of bucket = " << trie_node->size << endl; //TESTING
      cout << "end_bucket-start_bucket = " << (end_bucket-start_bucket) << endl; //TESTING
            
      // traverse the bucket
      while ( start_bucket != end_bucket )
      {
        uint32_t entry_str = *((uint32_t*)start_bucket);
        uint32_t s = partialwin-depth;
        total_entries++;
        
        // for each nt in the string
        for ( uint32_t j = 0; j < s; j++ )
        {
          kmer_keep.push_back((char)get_char[entry_str&3]); //TESTING
          entry_str>>=2;
        }
        
        cout << kmer_keep << endl; //TESTING
        
        for ( uint32_t j = 0; j < s; j++ ) kmer_keep.pop_back(); //TESTING
        
        // next entry
        start_bucket+=ENTRYSIZE;
      }//~for each entry
      
      kmer_keep.pop_back(); //TESTING
            
    }
        
    // the node element is empty, go to next node element
    else if ( value == 0 )
    {
      ;
    }
        
    trie_node++;
        
  }//~trie nodes
    
  return;
    
}//~taversetrie_debug()
#endif




/*
 *
 * @function load index: read from binary file the L/2-mer look-up
 * tables, the 19-mer position tables and the mini-burst tries
 * @return void
 * @version 1.0 Jan 16, 2013
 *
 *******************************************************************/
void
load_index( char* ptr_dbindex,
           string part_str,
           kmer*& lookup_tbl,
           kmer_origin*& positions_tbl,
           uint32_t& number_elements,
           uint32_t lnwin)
{
    // STEP 1: load the kmer 'count' variables (dbname.kmer.dat)
    string ptr_dbindex_str = ptr_dbindex;
  ifstream inkmer( (char*)(ptr_dbindex_str + ".kmer_" + part_str + ".dat").c_str(), ios::in | ios::binary );
    
  if ( !inkmer.good() )
  {
    fprintf(stderr,"\n  ERROR: The index '%s' does not exist.\n", (char*)(ptr_dbindex_str + ".kmer_" + part_str + ".dat").c_str());
    fprintf(stderr,"  Make sure you have constructed your index using the command `indexdb'. See `indexdb -h' for help.\n\n");
    exit(EXIT_FAILURE);
  }
    
  lookup_tbl = new kmer[(1<<lnwin)]();
  if ( lookup_tbl == NULL )
  {
    fprintf(stderr,"\n  ERROR: failed to allocate memory for look-up table (paralleltraversal.cpp)\n\n");
    exit(EXIT_FAILURE);
  }
    
  uint32_t limit = 1<<lnwin;
    
  for ( uint32_t i = 0; i < limit; i++ ) inkmer.read(reinterpret_cast<char*>(&(lookup_tbl[i].count)), sizeof(uint32_t));
    
  inkmer.close();
    
  // STEP 2: load the burst tries ( bursttrief.dat, bursttrier.dat )
  ifstream btrie( (char*)(ptr_dbindex_str + ".bursttrie_" + part_str + ".dat").c_str(), ios::in | ios::binary );
    
  if ( !btrie.good() )
  {
    fprintf(stderr,"\n  ERROR: The index '%s' does not exist.\n", (char*)(ptr_dbindex_str + ".bursttrie_" + part_str + ".dat").c_str());
    fprintf(stderr,"  Make sure you have constructed your index using the command `indexdb'. See `indexdb -h' for help.\n\n");
    exit(EXIT_FAILURE);
  }
    
  // loop through all 9-mers
  for ( uint32_t i = 0; i < (uint32_t)(1<<lnwin); i++ )
  {
    uint32_t sizeoftries[2] = {0};
        
#ifdef see_binary_output
        cout << "9-mer = " << i; //TESTING
#endif
        
    // ptr to block of memory for two mini-burst tries
    char *dst = NULL;
        
    // the size of both mini-burst tries
    for ( int j = 0; j < 2; j++ )
    {
      btrie.read(reinterpret_cast<char*>(&sizeoftries[j]), sizeof(uint32_t));
            
#ifdef see_binary_output
            if ( j == 0 ) cout << "\tsizeoftrie f = " << sizeoftries[j]; //TESTING
            else cout << "\tsizeoftrie r = " << sizeoftries[j]; //TESTING
#endif
            
    }
        
#ifdef see_binary_output
        cout << "\tlookup_tbl[i].count = " << lookup_tbl[i].count << endl; //TESTING
#endif
        
    // allocate contiguous memory for both mini-burst tries if they exist
    if ( lookup_tbl[i].count != 0 )
    {
      dst = new char[(sizeoftries[0]+sizeoftries[1])]();
      if ( dst == NULL )
      {
        fprintf(stderr,"  ERROR: failed to allocate memory for mini-burst tries (paralleltraversal.cpp)\n");
        exit(EXIT_FAILURE);
      }
            
      // load 2 burst tries per 9-mer
      for ( int j = 0; j < 2; j++ )
      {
        // mini-burst trie exists
        if ( sizeoftries[j] != 0 )
        {
#ifdef see_binary_output
          if ( j == 0 ) cout << "forward burst-trie \n"; //TESTING
          if ( j == 1 ) cout << "reverse burst-trie \n"; //TESTING
#endif            
          // create a root trie node
          NodeElement newnode[4];
                    
          // copy the root trie node into the beginning of burst trie array
          memcpy( dst, &newnode[0], sizeof(NodeElement)*4);
          memset( dst, 0, sizeof(NodeElement)*4);
                    
          if ( j == 0 ) lookup_tbl[i].trie_F = (NodeElement*)dst;
          else lookup_tbl[i].trie_R = (NodeElement*)dst;
                    
          // queue to store the trie nodes as we create them
          deque<NodeElement*> nodes;
          nodes.push_back( (NodeElement*)dst );
          ((NodeElement *&) dst)+=4;
                    
          // queue to store the flags of node elements given in the binary file
          deque<char> flags;
                    
          // read the first trie node
          for ( int i = 0; i < 4; i++ )
          {
            char tmp;
            btrie.read(reinterpret_cast<char*>(&tmp), sizeof(char));
            flags.push_back(tmp);
#ifdef see_binary_output
                        cout << " " << (int)tmp; //TESTING
#endif
          }
                    
          // build the mini-burst trie
          while ( !nodes.empty() )
          {
            // ptr to traverse each trie node
            NodeElement* node = nodes.front();
                        
            // trie node elements
            for ( int i = 0; i < 4; i++ )
            {
              unsigned char flag = flags.front();
                            
              // what does the node element point to
              switch( flag )
              {
                // set values to 0
                case 0:
                {
                  node->flag = 0;
                  node->size = 0;
                  node->whichnode.trie = NULL;
                }
                  break;
                // trie node
                case 1:
                {
                  // read the trie node
                  for ( int i = 0; i < 4; i++ )
                  {
                    char tmp;
                    btrie.read(reinterpret_cast<char*>(&tmp), sizeof(char));
#ifdef see_binary_output
                    cout << " " << (int)tmp; //TESTING
#endif
                    flags.push_back(tmp);
                  }
                                    
                  node->flag = 1;
                  node->size = 0;
                  NodeElement newnode[4];
                  memcpy( (NodeElement*)dst, &newnode[0], sizeof(NodeElement)*4);
                  nodes.push_back( (NodeElement*)dst );
                  node->whichnode.trie = (NodeElement*)dst;
                  ((NodeElement *&) dst)+=4;
                                    
                }
                  break;
                // bucket
                case 2:
                {
                  uint32_t sizeofbucket = 0;
                                    
                  // read the bucket info
                  btrie.read(reinterpret_cast<char*>(&sizeofbucket), sizeof(uint32_t));
                                    
#ifdef see_binary_output
                  cout << "\tsizeofbucket = " << sizeofbucket; //TESTING
#endif                      
                  char* bucket = new char[sizeofbucket]();
                  if ( bucket == NULL )
                  {
                    fprintf(stderr,"\n  %sERROR%s: failed to allocate memory for allocate bucket (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
                    exit(EXIT_FAILURE);
                  }
                                    
                  btrie.read(reinterpret_cast<char*>(bucket), sizeofbucket);
                                    
                  // copy the bucket into the burst trie array
                  memcpy( (void*)dst, (void*)bucket, sizeofbucket);
                                    
                  delete [] bucket;
                  bucket = NULL;
                                    
                  // assign pointers from trie node to the bucket
                  node->flag = flag;
                  node->whichnode.bucket = dst;
                  node->size = sizeofbucket;
                  dst = ((char *)dst)+sizeofbucket;        
                }
                  break;
                // ?
                default:
                {
                  fprintf(stderr, "\n  %sERROR%s: flag is set to %d (load_index)\n","\033[0;31m","\033[0m",flag);
                  exit(EXIT_FAILURE);
                }
                  break;
              }

              flags.pop_front();
              node++;               
            }//~loop through 4 node elements in a trie node
                        
#ifdef see_binary_output
            cout << "\n"; //TESTING
#endif            
            nodes.pop_front();
                        
          }//~while !nodes.empty()
        }//~if mini-burst trie exists
        else
        {
          if ( j == 0 ) lookup_tbl[i].trie_F = NULL;
          else lookup_tbl[i].trie_R = NULL;
        }
      }//~for both mini-burst tries
    }//~if ( sizeoftries != 0 )
    else
    {
      lookup_tbl[i].trie_F = NULL;
      lookup_tbl[i].trie_R = NULL;
    }
  }//~for all 9-mers in the look-up table
    
  btrie.close();
    
  // STEP 3: load the position reference tables (pos.dat)
  ifstream inreff( (char*)(ptr_dbindex_str + ".pos_" + part_str + ".dat").c_str(), ios::in | ios::binary );
    
  if ( !inreff.good() )
  {
    fprintf(stderr,"\n  ERROR: The database name '%s' does not exist.\n\n", (char*)(ptr_dbindex_str + ".pos_" + part_str + ".dat").c_str());
    exit(EXIT_FAILURE);
  }
    
  uint32_t size = 0;
    
  inreff.read(reinterpret_cast<char*>(&number_elements), sizeof(uint32_t));
    
  positions_tbl = new kmer_origin[number_elements]();
  if ( positions_tbl == NULL )
  {
    fprintf(stderr,"  ERROR: could not allocate memory for positions_tbl (main(), paralleltraversal.cpp)\n");
    exit(EXIT_FAILURE);
  }
    
  for ( uint32_t i = 0; i < number_elements; i++ )
  {
    /* the number of positions */
    inreff.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    positions_tbl[i].size = size;
        
    /* the sequence seq_pos array */
    positions_tbl[i].arr = new seq_pos[size]();
        
    if ( positions_tbl[i].arr == NULL )
    {
      fprintf(stderr, "  ERROR: could not allocate memory for positions_tbl (paralleltraversal.cpp)\n");
      exit(EXIT_FAILURE);
    }
        
    inreff.read(reinterpret_cast<char*>(positions_tbl[i].arr), sizeof(seq_pos)*size);
  }
    
  inreff.close();
    
  return ;
    
}//~load_index()




/*
 *
 * FUNCTION   : find_lis()
 * PARAMETERS : deque<mypair > &a, vector<int> &b
 * PURPOSE  : given a list of matching positions on the read, find the longest strictly increasing
 *            subsequence, O(n log k)
 * INPUT  : a list of matching positions on the read which fall within a range of the read's
 *            length on the genome, and a vector<int> &b to store the starting positions of each longest
 *            subsequence
 * OUTPUT : N/A
 *
 **************************************************************************************************************/
void find_lis( deque<pair<uint32_t, uint32_t> > &a, vector<uint32_t> &b, uint32_t readn )
{
  vector<uint32_t> p(a.size());
  int u, v;
    
  if (a.empty()) return;
    
  b.push_back(0);
    
  for (uint32_t i = 1; i < a.size(); i++)
  {
    // If next element a[i] is greater than last element of current longest subsequence a[b.back()], just push it at back of "b" and continue
    if (a[b.back()].second < a[i].second)
    {
      p[i] = b.back();
      b.push_back(i);
      continue;
    }
        
    // Binary search to find the smallest element referenced by b which is just bigger than a[i]
    // Note : Binary search is performed on b (and not a). Size of b is always <=k and hence contributes O(log k) to complexity.
    for (u = 0, v = b.size()-1; u < v;)
    {
      int c = (u + v) / 2;
      if ( a[b[c]].second < a[i].second )
        u=c+1;
      else
        v=c;
    }
        
    // Update b if new value is smaller then previously referenced value
    if (a[i].second < a[b[u]].second)
    {
      if (u > 0) p[i] = b[u-1];
      b[u] = i;
    }
  }
    
    
  for (u = b.size(), v = b.back(); u--; v = p[v]) b[u] = v;
}


/*! @fn paralleltraversal() */
void
paralleltraversal ( char* inputreads,
                   char* ptr_filetype_ar,
                   char* ptr_filetype_or,
                   int32_t match,
                   int32_t mismatch,
                   int32_t gap_open,
                   int32_t gap_extension,
                   int32_t score_N,
                   vector< vector<uint32_t> >& skiplengths,
                   int argc,
                   char **argv,
                   bool yes_SQ,
                   vector< pair<string,string> >& myfiles,
                   bool exit_early)
{
  // the offset from the start of the reads file for mmap
  off_t offset_map = 0;
  // file descriptor to find statistics on the reads file
  int fd = -1;
  // input reads file (fasta or fastq)
  string fname = inputreads;
  // the size of the full reads file (in bytes)
  off_t full_file_size = 0;
  // total number of nucleotides in all reads
  size_t full_read_main = 0;
  // total number of reads
  uint32_t number_total_read = 0;
  // total number of reads mapped passing E-value threshold
  uint32_t total_reads_mapped = 0;
  // total number of reads mapped passing E-value threshold & %id and/or %query coverage thresholds
  uint32_t total_reads_mapped_cov = 0;
  // total number of reads for de novo clustering
  uint32_t total_reads_denovo_clustering = 0;
  // the minimum occurrences of a (L/2)-mer required to allow search for a seed of length L in the burst tries
  uint32_t minoccur = 0;
  // for timing different processes
  double s,f;
  // the comparing character used in parsing the reads file
  char filesig;
  // minimum read length for log statistics
  uint32_t min_read_len = READLEN;
  // maximum read length for log statistics
  uint32_t max_read_len = 0;
  // mean read length for log statistics
  uint32_t mean_read_len = 0;
  // the size of the sliding window on the full file, ~1GB
  off_t partial_file_size = 0;
  // size of the remainder of the file (last window) which
  // is less than pow(2,30) bytes
  off_t last_part_size = 0;
  // number of file sections to mmap
  uint32_t file_sections = 0;
  // index for file_sections
  uint32_t file_s = 0;
    
  // check file for mmap
  if ((fd = open(fname.c_str(), O_RDONLY)) == -1)
  {
    fprintf(stderr,"  %sERROR%s: Could not open the reads file!\n\n",
                   "\033[0;31m","\033[0m");
    exit(EXIT_FAILURE);
  }

  // check which file format to parse: fasta or fastq
  char c;
  int32_t rb = -1;
    
  if ((rb = read(fd, &c, 1)) == -1)
  {
    fprintf(stderr,"  %sERROR%s: Could not read the first character of the "
                   "reads file!\n\n","\033[0;31m","\033[0m");
    exit(EXIT_FAILURE);
  }
    
  // set the appropriate settings for the file format
  if ( c == '>' )
  {
    // fasta format
    filesig = '>';
  }
  else if ( c == '@' )
  {
    // fastq format
    filesig = '@';
  }
  else exit_early = true;
    
  if ( !exit_early )
  {
    eprintf("\n  Computing read file statistics ...");
    TIME(s);
    // find the total length of all the reads for computing the E-value
    char ch;
    FILE *fp = fopen(inputreads,"r");
    if ( fp == NULL )
    {
      fprintf(stderr,"  %sERROR%s: could not open reads file - %s\n\n",
                     "\033[0;31m","\033[0m",strerror(errno));
      exit(EXIT_FAILURE);
    }
    // FASTA
    if ( filesig == '>' )
    {
      while ( (ch = getc(fp)) != EOF )
      {
        // sequence label
        if ( ch == '>' )
        {
          number_total_read++;
          while ( (ch = getc(fp)) != '\n' );
        }
              
        ch = getc(fp);
        
        // nucleotide sequence    
        while ( (ch != EOF) && (ch != '>') )
        {
          if ( (ch != '\n') && (ch != ' ') ) full_read_main++;
          ch = getc(fp);
        }
        ungetc(ch,fp);
      }
    }
    // FASTQ
    else
    {
      int nc = 0;
          
      while ( (ch = getc(fp)) != EOF )
      {
        if ( ch == '\n' ) nc++;
        if ( ((nc-1)%4 == 0) && (ch != '\n')) full_read_main++;
      }
          
      number_total_read=nc/4;
    }

    // find the mean sequence length
    mean_read_len = full_read_main/number_total_read;
      
    fclose(fp);
    
    // check there are an even number of reads for --paired-in
    // and --paired-out options to work
    if ( (number_total_read%2 != 0) && (pairedin_gv || pairedout_gv) )
    {
      fprintf(stderr,"\n    %sWARNING%s: for --paired-in and --paired-out options, the number of reads must be even.\n","\033[0;33m","\033[0m");
      fprintf(stderr,"    There are %d reads in your file.\n",number_total_read);
      fprintf(stderr,"    Reads will still be processed and output, but paired-reads may be split.\n\n");
      pairedin_gv = false;
      pairedout_gv = false;
    }

    // find the size of the total file
    if ((full_file_size = lseek(fd, 0L, SEEK_END)) == -1)
    {
      fprintf(stderr,"  %sERROR%s: Could not seek the reads file!\n\n","\033[0;31m","\033[0m");
      exit(EXIT_FAILURE);
    }
    if (lseek(fd, 0L, SEEK_SET) == -1)
    {
      fprintf(stderr,"  %sERROR%s: Could not seek set the reads file!\n\n","\033[0;31m","\033[0m");
      exit(EXIT_FAILURE);
    }

    partial_file_size = full_file_size;
    last_part_size = full_file_size%map_size_gv;
    
    // if the full_file_size is bigger than m*PAGE_SIZE, mmap
    // the file by 'windows' of size partial_file_size,
    // otherwise keep the full_file_size
    if ( ( file_sections = ceil( (double)full_file_size/(double)(map_size_gv) ) ) > 1 ) partial_file_size = map_size_gv;
    TIME(f);
      
    eprintf(" done [%.2f sec]\n", (f-s));
    eprintf("  size of reads file: %lu bytes\n", (unsigned long int)full_file_size );
    eprintf("  partial section(s) to be executed: %d of size %lu bytes \n", file_sections,(unsigned long int)partial_file_size );
  }//~if (!exit_early)
    
  // output streams for accepted reads (FASTA/FASTQ, SAM and BLAST-like)
  ofstream acceptedreads;
  ofstream acceptedsam;
  ofstream acceptedblast;
    
  // determine the suffix (fasta,fastq..) of accepted strings
  char suffix[20] = "out";
  char *suff = strrchr( inputreads, '.');
  if ( suff != NULL )
    strcpy( suffix, suff+1 );
  else if ( filesig == '>' )
    strcpy( suffix, "fasta");
  else
    strcpy( suffix, "fastq");
  suff = NULL;
    
  char *acceptedstrings = NULL;
  char *acceptedstrings_sam = NULL;
  char *acceptedstrings_blast = NULL;
  char *logoutfile = NULL;
  char *denovo_otus_file = NULL;
  char *acceptedotumap_file = NULL;
  
  // attach pid to output files
  char pidStr[4000];
  if ( pid_gv )
  {
    int32_t pid = getpid();
    sprintf(pidStr,"%d",pid);
  }
    
  // associate the streams with reference sequence file names
  if ( ptr_filetype_ar != NULL )
  {
    if ( fastxout_gv )
    {
      // fasta/fastq output
      acceptedstrings = new char[1000]();
      if ( acceptedstrings == NULL )
      {
        fprintf(stderr,"  %sERROR%s: could not allocate memory for acceptedstrings (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
      }
      strcpy ( acceptedstrings, ptr_filetype_ar );
      if ( pid_gv )
      {
        strcat ( acceptedstrings, "_");
        strcat ( acceptedstrings, pidStr);
      }
      strcat ( acceptedstrings, ".");
      strcat ( acceptedstrings, suffix);
            
      acceptedreads.open ( acceptedstrings );
      acceptedreads.close();
    }
    if ( samout_gv )
    {
      // sam output
      acceptedstrings_sam = new char[1000]();
      if ( acceptedstrings_sam == NULL )
      {
        fprintf(stderr,"  %sERROR%s: could not allocate memory for acceptedstrings_sam (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
      }
      strcpy( acceptedstrings_sam, ptr_filetype_ar );
      if ( pid_gv )
      {
        strcat( acceptedstrings_sam, "_");
        strcat( acceptedstrings_sam, pidStr);
      }
      strcat( acceptedstrings_sam, ".sam");            
      acceptedsam.open ( acceptedstrings_sam );
      acceptedsam.close();
    }
    if ( blastout_gv )
    {
      // blast output
      acceptedstrings_blast = new char[1000]();
      if ( acceptedstrings_blast == NULL )
      {
        fprintf(stderr,"  %sERROR%s: could not allocate memory for acceptedstrings_blast (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
      }
      strcpy( acceptedstrings_blast, ptr_filetype_ar );
      if ( pid_gv )
      {
        strcat( acceptedstrings_blast, "_");
        strcat( acceptedstrings_blast, pidStr);
      }
      strcat( acceptedstrings_blast, ".blast");            
      acceptedblast.open ( acceptedstrings_blast );
      acceptedblast.close();
    }
    if ( logout_gv )
    {
      // statistics file output
      ofstream logstream;
      logoutfile = new char[1000]();
      if ( logoutfile == NULL )
      {
        fprintf(stderr,"  %sERROR%s: could not allocate memory for acceptedstrings_blast (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
      }
      strcpy( logoutfile, ptr_filetype_ar );
      if ( pid_gv )
      {
        strcat( logoutfile, "_");
        strcat( logoutfile, pidStr);
      }
      strcat( logoutfile, ".log");
            
      logstream.open ( logoutfile );
      logstream.close();
    }
    if ( otumapout_gv )
    {
      // OTU map output file
      ofstream otumap;
      acceptedotumap_file = new char[1000]();
      if ( acceptedotumap_file == NULL )
      {
        fprintf(stderr,"  %sERROR%s: could not allocate memory for acceptedotumap (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
      }
      strcpy( acceptedotumap_file, ptr_filetype_ar );
      if ( pid_gv )
      {
        strcat( acceptedotumap_file, "_");
        strcat( acceptedotumap_file, pidStr);
      }
      strcat( acceptedotumap_file, "_otus.txt");        
      otumap.open ( acceptedotumap_file );
      otumap.close();
    }
    if ( de_novo_otu_gv )
    {
      ofstream denovo_otu;
      denovo_otus_file = new char[1000]();
      if ( denovo_otus_file == NULL )
      {
        fprintf(stderr,"  %sERROR%s: could not allocate memory for denovo_otus_file_name (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
      }
      strcpy( denovo_otus_file, ptr_filetype_ar );
      if ( pid_gv )
      {
        strcat( denovo_otus_file, "_");
        strcat( denovo_otus_file, pidStr);
      }
      strcat( denovo_otus_file, "_denovo.");
      strcat ( denovo_otus_file, suffix);
        
      denovo_otu.open ( denovo_otus_file );
      denovo_otu.close();
    }
  }//~if ( ptr_filetype_ar != NULL )
    
  if ( ptr_filetype_or != NULL )
  {
    if ( fastxout_gv )
    {
      // output stream for other reads
      ofstream otherreads;
            
      // add suffix database name to accepted reads file
      if ( pid_gv )
      {
        strcat( ptr_filetype_or, "_");
        strcat( ptr_filetype_or, pidStr);
      }
      strcat ( ptr_filetype_or, "." );
      strcat ( ptr_filetype_or, suffix );
            
      // create the other reads file
      otherreads.open( ptr_filetype_or );
      otherreads.close();
    }
  }

  // empty output files created, exit program
  if ( exit_early )
  {
    fprintf(stderr, "  The input reads file or reference file is empty, "
                    "or the reads file is not in FASTA or FASTQ format, "
                    "no analysis could be made.\n");
    // output parameters to log file
    if ( (ptr_filetype_ar != NULL) && logout_gv )
    {
      FILE* bilan = fopen(logoutfile,"w");
      if ( bilan == NULL )
      {
        fprintf(stderr,"  %sERROR%s: could not open file %s \n","\033[0;31m","\033[0m",logoutfile);
        exit(EXIT_FAILURE);
      }
      fprintf(bilan, "  The input reads file or reference file is empty, "
                     "or the reads file is not in FASTA or FASTQ format, "
                     "no analysis could be made.\n");
      fclose(bilan);
    }
    exit(EXIT_SUCCESS);
  }
    
  int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
    
  {
    int32_t l,k,m;
    // initialize Smith-Waterman scoring matrix for genome sequences
    for (l = k = 0; l < 4; ++l)
    {
      for (m = 0; m < 4; ++m)
        mat[k++] = l == m ? match : mismatch; // weight_match : weight_mismatch (must be negative)
      mat[k++] = score_N; // ambiguous base
    }
    for ( m = 0; m < 5; ++m ) mat[k++] = score_N; // ambiguous base
  }
       
  // the number of parts an index was divided into to fit into specified memory, for each reference database searched
  vector<uint16_t> num_index_parts(myfiles.size(),0);
  vector<uint64_t> full_ref(myfiles.size(),0);
  vector<uint64_t> full_read(myfiles.size(),full_read_main);
  vector<uint32_t> lnwin(myfiles.size(),0);
  vector<uint32_t> partialwin(myfiles.size(),0);
  // Minimal SW score in order to reach threshold E-value
  vector<uint32_t> minimal_score(myfiles.size(),0);
  // array of structs storing information on which sequences from the original FASTA file were added to each index part
  vector<vector<index_parts_stats> > index_parts_stats_vec;
  // Gumbel parameters lambda and K, respectively
  vector<pair<double,double> > gumbel(myfiles.size(),pair<double,double>(-1.0,-1.0));
  // total number of full bitvectors in [w_1] reverse and [w_2] forward
  vector<uint32_t> numbvs(myfiles.size(),0);
  // total number of reads matched for each database in list --ref
  vector<uint32_t> reads_matched_per_db(myfiles.size(),0);
  // number of reference sequences in each index FASTA file
  vector<uint32_t> numseq(myfiles.size(),0);
        
  // set the same skiplengths for all reference files (if the user uses option --passes)
  if ( skiplengths.empty() )
  {
    vector<uint32_t> skiplengths_v(3,0);
    skiplengths.push_back(skiplengths_v);
  }
  for (int32_t i = 0; i < myfiles.size()-1; i++) skiplengths.push_back(skiplengths[0]);
    
  // add header lines to SAM output file
  preprocess_data( myfiles,
                  argv,
                  argc,
                  yes_SQ,
                  acceptedstrings_sam,
                  match,
                  mismatch,
                  gap_open,
                  gap_extension,
                  skiplengths,
                  num_index_parts,
                  index_parts_stats_vec,
                  full_ref,
                  full_read,
                  lnwin,
                  partialwin,
                  minimal_score,
                  number_total_read,
                  gumbel,
                  numbvs,
                  numseq);
    
  // some info on chosen parameters
  eprintf("  Parameters summary:\n");
  eprintf("    Number of seeds = %d\n",seed_hits_gv);
  eprintf("    Edges = %d",edges_gv);
  if (as_percent_gv)
      eprintf(" (as percent)\n");
  else
      eprintf(" (as integer)\n");
  eprintf("    SW match = %d\n",match);
  eprintf("    SW mismatch = %d\n",mismatch);
  eprintf("    SW gap open penalty = %d\n",gap_open);
  eprintf("    SW gap extend penalty = %d\n",gap_extension);
  eprintf("    SW ambiguous nucleotide = %d",score_N);
  if ( score_N > 0 ) eprintf(" %sWarning!%s Positive score set for ambiguous nucleotides.\n","\033[0;33m","\033[0m");
  else eprintf("\n");
  if ( yes_SQ )
      eprintf("    SQ tags are output\n");
  else
      eprintf("    SQ tags are not output\n");
#ifdef _OPENMP
  eprintf("    Number of threads = %d\n",numcpu_gv);
#else
  eprintf("    Number of threads = 1 (OpenMP is not supported with your current C++ compiler).\n");
#endif
  if ( pid_gv )
  {
    eprintf("    Current process pid = %d\n",getpid());
  }
    
  // output parameters to log file
  if ( (ptr_filetype_ar != NULL) && logout_gv )
  {
    FILE* bilan = fopen(logoutfile,"w");
    if ( bilan == NULL )
    {
      fprintf(stderr,"  %sERROR%s: could not open file %s \n","\033[0;31m","\033[0m",logoutfile);
      exit(EXIT_FAILURE);
    }
    
    time_t q = time(0);
    struct tm * now = localtime(&q);
    
    fprintf(bilan," %s\n",asctime(now));
    fprintf(bilan," Command: ");
    for ( int j = 0; j < argc; j++ ) fprintf(bilan,"%s ",argv[j]);
    fprintf(bilan,"\n");
    // some info on chosen parameters
    fprintf(bilan," Process pid = %d\n",getpid());
    fprintf(bilan," Parameters summary:\n");
    for ( uint32_t index_num = 0; index_num < myfiles.size(); index_num++ )
    {
      fprintf(bilan,"    Index: %s\n",(char*)(myfiles[index_num].second).c_str());
      fprintf(bilan,"     Seed length = %d\n",lnwin[index_num]);
      fprintf(bilan,"     Pass 1 = %d, Pass 2 = %d, Pass 3 = %d\n",skiplengths[index_num][0],skiplengths[index_num][1],skiplengths[index_num][2]);
      fprintf(bilan,"     Gumbel lambda = %f\n",(gumbel[index_num].first));
      fprintf(bilan,"     Gumbel K = %f\n",(gumbel[index_num].second));
      fprintf(bilan,"     Minimal SW score based on E-value = %d\n",minimal_score[index_num]);
    }
    fprintf(bilan,"    Number of seeds = %d\n",seed_hits_gv);
    fprintf(bilan,"    Edges = %d",edges_gv);
    if (as_percent_gv)
      fprintf(bilan," (as percent)\n");
    else
      fprintf(bilan," (as integer)\n");
    fprintf(bilan,"    SW match = %d\n",match);
    fprintf(bilan,"    SW mismatch = %d\n",mismatch);
    fprintf(bilan,"    SW gap open penalty = %d\n",gap_open);
    fprintf(bilan,"    SW gap extend penalty = %d\n",gap_extension);
    fprintf(bilan,"    SW ambiguous nucleotide = %d",score_N);
    if ( score_N > 0 ) fprintf(bilan," <-- %sWarning!%s Positive score set for ambiguous nucleotides.\n","\033[0;33m","\033[0m");
    else fprintf(bilan,"\n");
    if ( yes_SQ )
      fprintf(bilan,"    SQ tags are output\n");
    else
      fprintf(bilan,"    SQ tags are not output\n");
#ifdef _OPENMP
    fprintf(bilan,"    Number of threads = %d\n",numcpu_gv);
#else
    fprintf(bilan,"    Number of threads = 1 (OpenMP is not supported with your current C++ compiler).\n");
#endif
    fprintf(bilan,"    Reads file = %s\n\n",inputreads);
        
    fclose(bilan);
  }

  // pointer to the split read
  // (the read which is split between any two file sections)
  char* split_read = NULL;

  // pointer to the position in the split read where to attach
  // the connecting part of the split read (and possibly its pair)
  char* split_read_ptr = NULL;
    
  // the number of lines to offset at the top of the current
  // file section
  int32_t offset_pair_from_top = 0;
    
  // map<reference sequence, vector<list of reads aligned to reference sequence> > otu_map
  map<string,vector<string> > otu_map;
    
  // Loop through all mmap'd read file sections
  while ( file_s < file_sections )
  {
    if ( samout_gv )
    {
      acceptedsam.open(acceptedstrings_sam, ios::app);
      if (!acceptedsam.good())
      {
        fprintf(stderr,"  %sERROR%s: could not open SAM output file for writing: %s.\n","\033[0;31m","\033[0m",acceptedstrings_sam);
        exit(EXIT_FAILURE);
      }
    }
    if ( blastout_gv )
    {
      acceptedblast.open(acceptedstrings_blast, ios::app);
      if (!acceptedblast.good())
      {
        fprintf(stderr,"  %sERROR%s: could not open BLAST output file for writing: %s.\n","\033[0;31m","\033[0m",acceptedstrings_blast);
        exit(EXIT_FAILURE);
      }
    }
        
    eprintf("\n  %sBegin mmap reads section # %d%s:\n","\033[4m",file_s+1,"\033[0m");
        
    // begin file memory map
    TIME(s);
        
    // mmap the reads file into memory
    char* raw = (char*)mmap ( 0, partial_file_size, PROT_READ, MAP_SHARED, fd, offset_map );
    if ( raw == MAP_FAILED )
    {
      close(fd);
      fprintf(stderr,"  %sERROR%s: cannot mmap file: %s\n\n","\033[0;31m","\033[0m",strerror(errno));
      exit(EXIT_FAILURE);
    }
    // pointer to last character in mmap'd region
    char* end_of_mmap = &raw[partial_file_size-1];
    
    int32_t strs = 0;
    
    // the length of the split read in file part i+1 (from beginning of file part)
    int32_t reads_offset_f = 0;
    // the length of the split read in file part i (from end of file part)
    int32_t reads_offset_e = 0;
        
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
            fprintf(stderr,"   %sERROR%s: 0 sequences mapped in the current file section.\n","\033[0;31m","\033[0m");
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
          for ( uint32_t i = reads_offset_f; i < partial_file_size-reads_offset_e-2; i++ ) if ( raw[i] == '>' ) strs++;
          
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
        for ( uint32_t i = reads_offset_f; i < partial_file_size; i++ ) if ( raw[i] == '\n' ) strs++;
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
            fprintf(stderr,"   %sERROR%s: Your FASTQ reads file has an uneven number of lines: %u\n","\033[0;31m","\033[0m",strs);
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
      fprintf(stderr,"\n  %sERROR%s: cannot allocate memory for reads\n\n","\033[0;31m","\033[0m");
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
        
    // a pointer is added to each header and the directly following read, hence the number of pointers for a fasta or fastq file is the same
    char* line = &raw[reads_offset_f];
    char* finalnt = &raw[partial_file_size-reads_offset_e-1];
    int minlenread = 1000000;
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
        if ( readlen >= 20 ) readlen < minlenread ? minlenread = readlen : minlenread;
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
          if ( readlen >= 20 ) readlen < minlenread ? minlenread = readlen : minlenread;
                    
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
    if ( minlenread == 1000000 )
    {
      fprintf(stderr,"   %sERROR%s: All reads are too short (<22nt) for further analysis.\n\n","\033[0;31m","\033[0m");
      exit(EXIT_FAILURE);
    }
        
    TIME(f);
        
    eprintf("  Time to mmap reads and set up pointers [%.2f sec]\n", (f-s) );
                
    // array of bits to represent all reads
    // a bit representing an accepted read is set to 1
    vector<bool> read_hits(strs, false);
    // array of bits to represent all reads
    // a bit representing an accepted read with < %id 
    // and < %coverage is set to 0
    vector<bool> read_hits_denovo(strs, true);
    
    // array of uint16_t to represent all reads, if the read was aligned with a maximum SW score, its number of alignments is incremeted by 1
    uint16_t *read_max_SW_score = new uint16_t[strs];
    memset(read_max_SW_score, 0, sizeof(uint16_t)*strs);
    
    // map accessed by read number, storing a pair <index for smallest SSW score, pointer to array of num_best_hits_gv>
    map<uint32_t, alignment_struct > read_hits_align_info;
    
    // number of alignments to output per read
    int32_t *num_alignments_x = NULL;
        
    // output num_alignments_gv alignments per read
    if ( num_alignments_gv > 0 )
    {
      num_alignments_x = new int32_t[strs];
      for ( int32_t s = 0; s < strs; s++ ) num_alignments_x[s] = num_alignments_gv;
    }
    
    // loop through every index passed to option --ref (ex. SSU 16S and SSU 18S)
    for ( uint16_t index_num = 0; index_num < (uint16_t)myfiles.size(); index_num++)
    {
      // covert part number into a string
      stringstream prt_str;
      uint16_t part = 0;
      prt_str << part;
      string part_str = prt_str.str();
      
      eprintf("\n  Begin analysis of: %s%s%s\n","\033[0;34m",(char*)(myfiles[index_num].first).c_str(),"\033[0m");
      if ( file_s == 0 )
      {
        eprintf("    Seed length = %d\n",lnwin[index_num]);
        eprintf("    Pass 1 = %d, Pass 2 = %d, Pass 3 = %d\n",skiplengths[index_num][0],skiplengths[index_num][1],skiplengths[index_num][2]);
        eprintf("    Gumbel lambda = %f\n",gumbel[index_num].first);
        eprintf("    Gumbel K = %f\n",gumbel[index_num].second);
        eprintf("    Minimal SW score based on E-value = %d\n",minimal_score[index_num]);
      }
      
      // for each partial file of burst trie index (part_0 .. part_x)
      for ( part = 0; part < num_index_parts[index_num]; part++ )
      {
        eprintf("    Loading index part %d/%u ... ",part+1,num_index_parts[index_num] );
                
        TIME(s);

        // number of reference sequences to search alignments for before choosing the best one
        int32_t *best_x = NULL;
        
        // search min_lis_gv reference sequences for alignments
        if ( min_lis_gv > 0 )
        {
          best_x = new int32_t[strs];
          for ( int32_t s = 0; s < strs; s++ )
            best_x[s] = min_lis_gv;
        }
                
        // memory buffer to store the reference sequence database
        char* buffer = NULL;
        // pointer to start of each sequence in buffer
        char** reference_seq = NULL;
        // length of each sequence in buffer
        uint32_t* reference_seq_len = NULL;
        // 9-mer look-up tables
        kmer *lookup_tbl = NULL;
        // 19-mer position look-up tables
        kmer_origin* positions_tbl = NULL;
        // number of elements in the table
        uint32_t number_elements = 0;
                
        uint64_t seq_part_size = index_parts_stats_vec[index_num][part].seq_part_size;
        uint32_t numseq_part = index_parts_stats_vec[index_num][part].numseq_part;
        uint64_t start_part = index_parts_stats_vec[index_num][part].start_part;
                
#pragma omp master
        {
          // 2. load the index part (9-mer lookup table, mini-burst tries and positions table)
          load_index((char*)(myfiles[index_num].second).c_str(), part_str, lookup_tbl, positions_tbl, number_elements, lnwin[index_num] );
                    
          // block of memory to hold all ids + reference sequences
          buffer = new char[(seq_part_size+1)]();
          if ( buffer == NULL )
          {
            fprintf(stderr,"    %sERROR%s: could not allocate memory for reference sequence buffer (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
                    
          // pointer to the start of every sequence in the buffer
          reference_seq = new char*[(numseq_part<<1)]();
          if ( reference_seq == NULL )
          {
            fprintf(stderr,"    %sERROR%s: could not allocate memory for reference_seq (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
                    
          // length of every sequence in the buffer
          reference_seq_len = new uint32_t[numseq_part]();
          if ( reference_seq_len == NULL )
          {
            fprintf(stderr,"    %sERROR%s: could not allocate memory for reference_seq_len (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
                    
          // load the reference sequences for SW alignment
          load_ref((char*)(myfiles[index_num].first).c_str(),buffer,reference_seq,reference_seq_len,seq_part_size,numseq_part,start_part,1);
        }
        TIME(f);
                
        eprintf(" done [%.2f sec]\n",(f-s) );
                
        eprintf("    Begin index search ... ");
                      
        // begin the parallel traversal
        TIME(s);
                
        uint32_t bit_vector_size = (partialwin[index_num]-2)<<2;
        uint32_t offset = (partialwin[index_num]-3)<<2;
                
        // only search the forward xor reverse strand
        int32_t max = 0;
                
        forward_gv = true;
                
        if ( forward_gv ^ reverse_gv ) max = 1;
        // search both strands
        else max = 2;
                
        // search the forward and/or reverse strands
        for ( int32_t strand = 0; strand < max; strand++ )
        {
          // loop through all of the reads in the file 
#pragma omp parallel for num_threads(numcpu_gv) shared(lookup_tbl,positions_tbl,buffer,reference_seq,reference_seq_len,read_hits_align_info,read_hits,read_max_SW_score) schedule(dynamic,256)
          for ( int32_t readn = 1; readn < strs; readn+=2 )
          {
#ifdef debug_align
            cout << "readn = " << readn << endl; //TESTING
#endif                   
            // for reverse reads
            if ( !forward_gv )
            {
              // output the first num_alignments_gv alignments
              if ( num_alignments_gv > 0 )
              {
                // all num_alignments_gv alignments have been output
                if ( num_alignments_x[readn] < 0 ) continue;
              }
              // the maximum scoring alignment has been found, go to next read
	      // (unless all alignments are being output)
              else if ( (num_best_hits_gv > 0) && (min_lis_gv > 0) && (read_max_SW_score[readn] == num_best_hits_gv) ) continue;
            }
                        
            // read on integer alphabet {0,1,2,3}
            char myread[READLEN] = "";
            char* str = reads[readn];
            char* _myread = &myread[0];
            // keep record of what position on the read an ambiguous char exists,
            // the char at this position will be changed to '4' during
            // sequence alignment
            int32_t ambiguous_nt[READLEN] = {0};
            // size of ambiguous_nt
            uint32_t size_ambiguous_nt = 0;
            // flag to count only 1 alignment per read
            bool read_to_count = true;         
            // length of read
            uint32_t readlen = 0;
                        
            // change the read into an integer alphabet -- FASTA
            if ( filesig == '>' )
            {
              // str != finalnt (finalnt marks the end of mmap'd region)
              // str != '>' ('>' marks the start of the next read)
              // str != '\0' ('\0' marks the end of split-read, which
              // must be terminated by '\0')
              while ( (str != finalnt) && (*str != '>') && (*str != '\0'))
              {
                if ( *str != '\n' )
                {
                  *_myread = nt_table[(int)*str];
                  if ( *_myread == 4 )
                  {
                    *_myread = 0;
                    ambiguous_nt[size_ambiguous_nt++] = readlen;
                  }
                  _myread++;
                  readlen++;
                  if ( readlen > READLEN )
                  {
                    fprintf(stderr,"\n  %sERROR%s: at least one of your reads is > %d nt \n",
                                   "\033[0;31m","\033[0m",READLEN);
                    fprintf(stderr,"  Please check your reads or contact the authors.\n");
                    exit(EXIT_FAILURE);
                  }
                }
                str++;
              }
              *_myread = '\n';
            }
            // change the read into an integer alphabet -- FASTQ
            else
            {
              while ( *str != '\n' )
              {
                *_myread = nt_table[(int)*str++];
                if ( *_myread == 4 )
                {
                  // searching the read using the universal levenshtein automaton
                  // can only be done on 0-3 alphabet as burst trie is built 
                  // on 0-3 alphabet
                  *_myread = 0;
                  ambiguous_nt[size_ambiguous_nt++] = readlen;
                }
                _myread++;
                readlen++;
                if ( readlen > READLEN )
                {
                  fprintf(stderr,"\n  %sERROR%s: at least one of your reads is > %d nt \n",
                                 "\033[0;31m","\033[0m",READLEN);
                  fprintf(stderr,"  Please check your reads or contact the authors.\n");
                  exit(EXIT_FAILURE);
                }
              }
              *_myread = '\n';
            }

            // find the minimum sequence length
            readlen < min_read_len ? min_read_len = readlen : min_read_len;

            // find the maximum sequence length
            readlen > max_read_len ? max_read_len = readlen : max_read_len;
                        
            // the read length is too short
            if ( readlen < lnwin[index_num] )
            {
              fprintf(stderr,"\n  %sWARNING%s: At least one of the reads is shorter "
                             "than %u nucleotides, ","\033[0;33m","\033[0m",
                             lnwin[index_num]);
              fprintf(stderr,"by default it will not be searched\n ");
              continue;
            }
                        
            // create the reverse strand
            if ( !forward_gv )
            {
              // faster than xor algorithm
              char myread_rc[READLEN] = "";
              char* revcomp = &myread[readlen-1];
                            
              for ( uint32_t j = 0; j < readlen; j++ )
                myread_rc[j] = complement[(int)*revcomp--];
              myread_rc[readlen] = '\n';
                            
              memcpy(&myread[0],&myread_rc[0],READLEN);
                            
#ifdef debug_align
              cout << "read (before subst. 4's for N's): " << endl;
              char *cp = myread;
              while ( *cp != '\n' ) cout << (int)*cp++;
              cout << endl;
#endif
                            
            }//~if (REVERSE)
                        
            // array of positions of window hits on the reference sequence
            vector< id_win > id_win_hits;
            // number of windows hit to the reference sequence(s)
            uint32_t readhit = 0;
            uint32_t windowshift = skiplengths[index_num][0];
                        
            // keep track of windows which have been already traversed in the burst trie
            vector<bool> read_index_hits(readlen);
            
            // Pass number (possible value 0,1,2)
            uint32_t pass_n = 0;
            
            // the maximum SW score attainable for this read
            uint32_t max_SW_score = readlen*match;             
                        
            // loop for each new Pass to granulate seed search intervals
            bool search = true;
            do
            {
#ifdef debug_align
              cout << "\tpass = " << pass_n << endl; //TESTING
#endif          
              uint32_t numwin = (readlen-lnwin[index_num]+windowshift)/windowshift;
              uint32_t read_index = 0;
                            
              // iterate over windows of the template string
              for ( uint32_t win_num = 0; win_num < numwin; win_num++ )
              {
                // skip position, seed at this position has already been searched for in a previous Pass
                if ( read_index_hits[read_index] ) goto check_score;
                // search position, set search bit to true
                else read_index_hits[read_index].flip();                                
                                
                {
                  // this flag it set to true if a match is found during
                  // subsearch 1(a), to skip subsearch 1(b)
                  bool accept_zero_kmer = false;
                  // ids for k-mers that hit the database
                  vector< id_win > id_hits;
                                    
                  MYBITSET bitwindowsf[bit_vector_size];
                  memset(&bitwindowsf[0],0,bit_vector_size);
                                    
                  // build the first bitvector window
                  init_win_f( &myread[read_index+partialwin[index_num]],
                              // [w_1] forward k = 1
                              // bitwindows[0][0][0]
                              &bitwindowsf[0],
                              // bitwindows[0][1][0]
                              &bitwindowsf[4],
                              numbvs[index_num]);
                                    
                  uint32_t keyf = 0;
                  char *keyf_ptr = &myread[read_index];
                                    
                  // build hash for first half windows (foward and reverse)
                  for ( uint32_t g = 0; g < partialwin[index_num]; g++ )
                    (keyf <<= 2) |= (uint32_t)*keyf_ptr++;
                
                  // do traversal if the exact half window exists in the burst trie
                  if ( (lookup_tbl[keyf].count > minoccur) && (lookup_tbl[keyf].trie_F != NULL) )
                  {
                    /* subsearch (1)(a) d([p_1],[w_1]) = 0 and d([p_2],[w_2]) <= 1;
                     *
                     *  w = |------ [w_1] ------|------ [w_2] ------|
                     *  p = |------ [p_1] ------|------ [p_2] ----| (0/1 deletion in [p_2])
                     *              or
                     *    = |------ [p_1] ------|------ [p_2] ------| (0/1 match/substitution in [p_2])
                     *        or
                     *    = |------ [p_1] ------|------ [p_2] --------| (0/1 insertion in [p_2])
                     *
                     */
                                        
#ifdef debug_align
                    cout << "\tsearch forward mini-burst trie..\n"; //TESTING
#endif
                                        
                    traversetrie_align ( lookup_tbl[keyf].trie_F,
                                         0,
                                         0,
                                         // win2f_k1_ptr
                                         &bitwindowsf[0],
                                         // win2f_k1_full
                                         &bitwindowsf[offset],
                                         accept_zero_kmer,
                                         id_hits,
                                         readn,
                                         read_index,
                                         partialwin[index_num]);
                                        
#ifdef debug_align
                                        cout << "\tdone!\n"; //TESTING
#endif
                                        
                  }//~if exact half window exists in the burst trie
                                    
                  // only search if an exact match has not been found
                  if ( !accept_zero_kmer )
                  {
                    MYBITSET bitwindowsr[bit_vector_size];
                    memset(&bitwindowsr[0],0,bit_vector_size);
                                        
                    // build the first bitvector window
                    init_win_r( &myread[read_index+partialwin[index_num]-1],
                                // [w_1] reverse k = 1
                                // bitwindows[0][0][0]
                                &bitwindowsr[0],
                                // bitwindows[0][1][0]
                                &bitwindowsr[4],
                                numbvs[index_num]);
                                        
                    uint32_t keyr = 0;
                    char *keyr_ptr = &myread[read_index+partialwin[index_num]];
                                        
                    // build hash for first half windows (foward and reverse)
                    for ( uint32_t g = 0; g < partialwin[index_num]; g++ )
                      (keyr <<= 2) |= (uint32_t)*keyr_ptr++;
                                        
                    // continue subsearch (1)(b)
                    if ( (lookup_tbl[keyr].count > minoccur) && (lookup_tbl[keyr].trie_R != NULL) )
                    {
                      /* subsearch (1)(b) d([p_1],[w_1]) = 1 and d([p_2],[w_2]) = 0;
                       *
                       *  w =    |------ [w_1] ------|------ [w_2] -------|
                       *  p =      |------- [p_1] ---|--------- [p_2] ----| (1 deletion in [p_1])
                       *              or
                       *    =    |------ [p_1] ------|------ [p_2] -------| (1 match/substitution in [p_1])
                       *        or
                       *    = |------- [p_1] --------|---- [p_2] ---------| (1 insertion in [p_1])
                       *
                       */
                                            
#ifdef debug_align
                       cout << "\tsearch reverse mini-burst trie..\n"; //TESTING
#endif
                                            
                       traversetrie_align ( lookup_tbl[keyr].trie_R,
                                            0,
                                            0,
                                            /* win1r_k1_ptr */
                                            &bitwindowsr[0],
                                            /* win1r_k1_full */
                                            &bitwindowsr[offset],
                                            accept_zero_kmer,
                                            id_hits,
                                            readn,
                                            read_index,
                                            partialwin[index_num]);
                                            
#ifdef debug_align
                      cout << "\tdone!\n"; //TESTING
#endif                        
                    }//~if exact half window exists in the reverse burst trie                    
                  }//~if (!accept_zero_kmer)
                                             
                  // associate the ids with the read window number
                  if ( !id_hits.empty() )
                  {
                    for ( uint32_t i = 0; i < id_hits.size(); i++ )
                    {
                      id_win_hits.push_back(id_hits[i]);
                    }
                                        
                    readhit++;
                  }                  
                }
                                
                check_score:
                // boolean set to true if SW alignment succeeded between
                // the read and a candidate reference sequence
                bool aligned = false;
                          
                // output read if matched at more than RATIO windows
                if ( win_num == numwin-1 )
                {
                  // flag to check whether read was formatted to 0-4 alphabet
                  // for alignment allowing N's
                  bool read_edited = false;       
#ifdef debug_align
                  cout << "\t\t\treadhit = " << readhit << endl; //TESTING
                  cout << "\t\t\tseed_hits_gv = " << seed_hits_gv << endl; //TESTING
                  cout << "\t\t\tnum_best_hits_gv = " << num_best_hits_gv << endl; //TESTING
#endif                              
                  // STEP 1: the number of matching windows on the read to the
                  // reference database is greater than the threshold,
                  // continue analysis of read
                  if ( readhit >= (uint32_t)seed_hits_gv )
                  {
                    // map<seq, number of occurrences> most_frequent_seq_t
                    map<uint32_t,uint32_t> most_frequent_seq_t;
                    map<uint32_t,uint32_t>::iterator map_it;
                                        
                    uint32_t max_seq = 0;
                    uint32_t max_occur = 0;
                                        
                    // STEP 2: for every reference sequence, compute the number of
                    // window hits belonging to it
                    for ( uint32_t i = 0; i < id_win_hits.size(); i++ )
                    {
                      uint32_t _id = id_win_hits[i].id;
                                                                                        
                      // number of entries in the positions table for this id
                      uint32_t num_hits = positions_tbl[_id].size;
                                            
                      // pointer to the seq_pos array in the positions table for this id
                      seq_pos* positions_tbl_ptr = positions_tbl[_id].arr;
                                            
                      // loop through every position of id
                      for ( uint32_t j = 0; j < num_hits; j++ )
                      {
                        uint32_t seq = positions_tbl_ptr++->seq;
                        
                        // sequence already exists in the map, increment it's value
                        if ( (map_it=most_frequent_seq_t.find(seq)) != most_frequent_seq_t.end() )
                          map_it->second++;
                        // sequence doesn't exist, add it
                        else most_frequent_seq_t[seq] = 1;
                      }
                    }
                                                          
                    // <mypair> = <number of occurrences of a sequence, index of sequence>
                    vector<mypair> most_frequent_seq;
                                        
                    // copy list of occurrences from map to vector for sorting
                    for ( map_it = most_frequent_seq_t.begin(); map_it != most_frequent_seq_t.end(); map_it++ )
                    {
                      // pass candidate reference sequences for further analyses
                      // only if they have enough seed hits
                      if ( map_it->second >= (uint32_t)seed_hits_gv )
                        most_frequent_seq.push_back(mypair(map_it->second,map_it->first));
                    }
                    
                    most_frequent_seq_t.clear();
                                        
                    // sort the highest scoring sequences to the head of the array
                    sort(most_frequent_seq.begin(), most_frequent_seq.end(), largest);

		    /* tmp
		    cout << "\npass = " << pass_n << endl;
		    for ( int r = 0; r < most_frequent_seq.size(); r++ )
		    {
		      cout << most_frequent_seq[r].second << "\t";
		      char* tt = reference_seq[(2*(int)most_frequent_seq[r].second)];
		      while (*tt != ' ' ) cout << (char)*tt++;
		      cout << "\t " << most_frequent_seq[r].first << endl;
		    }
		    */

                    // STEP 3: for each reference sequence candidate
                    // (starting from highest scoring)
                    for ( uint32_t k = 0; k < most_frequent_seq.size(); k++ )
                    {
                      // the maximum scoring alignment has been found,
                      // do not search for anymore alignments
                      if ( (num_best_hits_gv != 0) && (read_max_SW_score[readn] == num_best_hits_gv) ) break;
                                            
                      max_occur = most_frequent_seq[k].first;
                      max_seq = most_frequent_seq[k].second; 
#ifdef debug_align
                      cout << "\t\t\t\tmax_occur = " << max_occur << endl; //TESTING
                      cout << "\t\t\t\tmax_seq = " << max_seq << endl; //TESTING
                      cout << "\t\t\t\tnumseq_part = " << numseq_part << endl; //TESTING
                      cout << "\t\t\t\tindex size for reference_seq = " << (numseq_part<<1) << endl; //TESTING
		      if ( min_lis_gv > 0 )
			cout << "\t\t\t\tbest_x[" << readn << "] = " << best_x[readn] << endl; //TESTING
#endif                             
                      // not enough window hits, try to collect more hits or go to next read
                      if ( max_occur < (uint32_t)seed_hits_gv ) break;
                                            
                      // update number of reference sequences remaining to check
                      if ( (min_lis_gv > 0) && aligned && (k > 0) )
                      {
                        // only decrement best_x if the next ref sequence to check
                        // has a lower seed count than the previous one
                        if ( max_occur < most_frequent_seq[k-1].first )
                        {
                          #pragma omp critical
                          {
                            best_x[readn]--;
                          }

                          if ( best_x[readn] < 1 ) break;
                        }
                      }
                      
                      // check if the maximum number of alignments per read
                      // (--num_alignments INT) have been output
                      if ( num_alignments_gv > 0 )
                      {
                        if ( num_alignments_x[readn] <= 0 ) break;
                      }  
#ifdef debug_align
                      cout << "\t\t\t\tpassed all checks for ref sequence search.. " << endl; //TESTING
#endif                                      
                      // STEP 4: collect all the genome positions belonging to the
                      // reference candidate from the table of positions computed
                      // during indexing
                      vector<mypair> hits_on_genome;
                      uint32_t count = 0;
                                            
                      for ( uint32_t i = 0; i < id_win_hits.size(); i++ )
                      {
                        uint32_t _id = id_win_hits[i].id;
                        uint32_t num_hits = positions_tbl[_id].size;
                        seq_pos* positions_tbl_ptr = positions_tbl[_id].arr;
                                                
                        // loop through every position of id
                        for ( uint32_t j = 0; j < num_hits; j++ )
                        {
                          if ( positions_tbl_ptr->seq == max_seq )
                          {
                            count++;
                            hits_on_genome.push_back(mypair(positions_tbl_ptr->pos,id_win_hits[i].win));
                          }
                          positions_tbl_ptr++;
                        }
                      }
                                                      
                      // sort the positions on genome in ascending order
                      sort(hits_on_genome.begin(), hits_on_genome.end(), smallest);      
                                            
                      // iterate over the set of hits, output windows of
                      // length == read which have at least ratio hits
                      vector<mypair>::iterator it3 = hits_on_genome.begin();
                      deque<mypair> vi_read;
                      deque<mypair>::iterator deq;
                                            
                      // STEP 5: run a sliding window of read's length across
                      // the genome, and search for windows with enough
                      // k-mer hits
                      uint32_t lcs_ref_start = 0;
                      uint32_t lcs_que_start = 0;                
                      uint32_t begin = it3->first;                                      
#ifdef debug_align
                      cout << "\t\t\t\tnumber of seed hits for this ref seq = " << hits_on_genome.size() << endl; //TESTING
#endif                                
                      while ( it3 != hits_on_genome.end() )
                      {
                        uint32_t stop = begin + readlen - lnwin[index_num] + 1;
                        bool push = false;
                                                
                        while ( (it3 != hits_on_genome.end()) && (it3->first <= stop) )
                        {
                          vi_read.push_back(*it3);
                          push = true;
                          it3++;
                        }
                        // heuristic 1: a new window hit was not pushed back, pop queue until new window can be pushed back
                        // this heuristic significantly speeds up the algorithm because we don't perform alignments for
                        // every sub-LIS of a window if an alignment reaching threshold has already been made. It assumes
                        // that every sub-LIS yields the same alignment score, which is true for 99.99% of cases.
#ifndef HEURISTIC1_OFF
                        if ( !push && aligned ) goto pop;
                        else aligned = false;
#endif
#ifdef HEURISTIC1_OFF
                        aligned = false;
#endif                              
                        // enough windows at this position on genome to search for LIS
                        if ( vi_read.size() >= (uint32_t)seed_hits_gv )
                        {
                          vector<uint32_t> list;
                          find_lis(vi_read, list, readn);                                
#ifdef HEURISTIC1_OFF
                          uint32_t list_n = 0;
                                                    
                          do
                          {
#endif                                      
                            // LIS long enough to perform Smith-Waterman alignment
                            if ( list.size() >= (uint32_t)seed_hits_gv )
                            {
#ifdef debug_align
                              cout << "\t\t\t\tLIS = " << list.size() << endl; //TESTING
                              for ( int p = 0; p < list.size(); p++ )
                                  cout << "\t\t\t\tref: " << vi_read[list[p]].first << "\tread: " << vi_read[list[p]].second << endl;   
#endif
#ifdef HEURISTIC1_OFF
                              lcs_ref_start = vi_read[list[list_n]].first;
                              lcs_que_start = vi_read[list[list_n]].second;
#endif
#ifndef HEURISTIC1_OFF
                              lcs_ref_start = vi_read[list[0]].first;
                              lcs_que_start = vi_read[list[0]].second;
#endif                                    
                              // reference string
                              uint32_t head = 0;
                              uint32_t tail = 0;
                              uint32_t align_ref_start = 0;
                              uint32_t align_que_start = 0;
                              uint32_t align_length = 0;
                              uint32_t reflen = reference_seq_len[max_seq];
                              uint32_t edges = 0;
                              if ( as_percent_gv )
                                edges = (((double)edges_gv/100.0)*readlen);
                              else
                                edges = edges_gv;
                              
                              // part of the read hangs off (or matches exactly) the beginning of the reference seq
                              //            ref |-----------------------------------|
                              // que |-------------------|
                              //             LIS |-----|
                              //
                              if ( lcs_ref_start < lcs_que_start )
                              {
                                align_ref_start = 0;
                                align_que_start = lcs_que_start - lcs_ref_start;
                                head = 0;
                                // the read is longer than the reference sequence
                                //            ref |----------------|
                                // que |---------------------...|
                                //                LIS |-----|
                                //
                                if ( reflen < readlen )
                                {
                                  tail = 0;
                                  // beginning from align_ref_start = 0 and align_que_start = X, the read finishes
                                  // before the end of the reference
                                  //            ref |----------------|
                                  // que |------------------------|
                                  //                  LIS |-----|
                                  //                ^
                                  //                align_que_start
                                  if ( align_que_start > (readlen - reflen) )
                                  {
                                    align_length = reflen - (align_que_start - (readlen - reflen));
                                  }
                                  // beginning from align_ref_start = 0 and align_que_start = X, the read finishes
                                  // after the end of the reference
                                  //            ref |----------------|
                                  // que |------------------------------|
                                  //                  LIS |-----|
                                  //                ^
                                  //                align_que_start
                                  else
                                  {
                                    align_length = reflen;
                                  }
                                }
                                else
                                {
                                  tail = reflen - align_ref_start - readlen;
                                  tail > (edges-1) ? tail = edges : tail;
                                  align_length = readlen+head+tail-align_que_start;
                                }
                              }
                              else
                              {
                                align_ref_start = lcs_ref_start - lcs_que_start;
                                align_que_start = 0;
                                align_ref_start > (edges-1) ? head = edges : head;
                                
                                // part of the read hangs off the end of the reference seq
                                // ref |-----------------------------------|
                                //                          que |-------------------|
                                //                            LIS |-----|
                                //
                                if ( align_ref_start+readlen > reflen )
                                {
                                  tail = 0;
                                  align_length = reflen-align_ref_start-head;
                                }
                                // the reference seq fully covers the read
                                // ref |-----------------------------------|
                                //    que |-------------------|
                                //          LIS |-----|
                                //
                                else
                                {
                                  tail = reflen - align_ref_start - readlen;
                                  tail > (edges-1) ? tail = edges : tail;
                                  align_length = readlen+head+tail;
                                }
                              }
                                                                                                                
                              // edit read in numeric alphabet to include ambiguous characters in alignment
                              if ( (size_ambiguous_nt != 0) && !read_edited )
                              {
                                read_edited = true;
                                
                                if ( !forward_gv )
                                {
                                  for ( uint32_t p = 0; p < size_ambiguous_nt; p++ )
                                  {
                                    myread[(readlen-ambiguous_nt[p]-1)] = 4;
                                  }
                                }
                                else
                                {
                                  for ( uint32_t p = 0; p < size_ambiguous_nt; p++ )
                                  {
                                    myread[ambiguous_nt[p]] = 4;
                                  }
                                }
                              }
                                                            
#ifdef debug_align
                              //if ( forward_gv ) //TESTING
                              //{
                              cout << "\n\t\t\t\t\treadlen = " << readlen << endl;
                              cout << "\t\t\t\t\treflen = " << reflen << endl;
                              cout << "\t\t\t\t\tlcs_que_start = " << lcs_que_start << endl;
                              cout << "\t\t\t\t\tlcs_ref_start = " << lcs_ref_start << endl;
                              cout << "\t\t\t\t\talign_que_start = " << align_que_start << endl;
                              cout << "\t\t\t\t\talign_ref_start = " << align_ref_start << endl;
                              cout << "\t\t\t\t\thead = " << head << "; tail = " << tail << endl;
                              cout << "\t\t\t\t\talign_length = " << align_length << endl;
                              cout << "read tag: ";
                              char* tt = reads[readn-1];
                              while (*tt != '\n' ) cout << (char)*tt++;
                              cout << endl;
                              cout << "read (starting from align_que_start): ";
                              tt = myread+align_que_start;
                              while ( *tt != '\n' ) cout << (int)*tt++;
                              cout << endl;
                              
                              
                              cout << "ref tag: ";
                              tt = reference_seq[(2*(int)max_seq)];
                              while (*tt != '\n' ) cout << (char)*tt++;
                              cout << endl;
                              cout << "ref (starting from align_ref_start-head): ";
                              tt = reference_seq[(2*(int)max_seq)+1]+align_ref_start-head;
                              while (*tt != '\n' ) cout << (int)*tt++;
                              cout << endl;
                              cout << "align_length-head-tail = " << (align_length-head-tail) << endl;
                              //}
#endif                                                     
                              // create profile for read
                              s_profile* profile = 0;
                              profile = ssw_init((int8_t*)(myread+align_que_start), (align_length-head-tail), mat, 5, 2);
                              
                              s_align* result = 0;
                              //uint16_t filters = 0;
                              
                              result = ssw_align( profile,
                                                 (int8_t*)reference_seq[(2*(int)max_seq)+1]+align_ref_start-head,align_length,
                                                 gap_open,
                                                 gap_extension,
                                                 2,
                                                 minimal_score[index_num],
                                                 0,
                                                 0);
                              
                              // deallocate memory for profile, no longer needed
                              if ( profile != 0 ) init_destroy(&profile);                                         
#ifdef debug_align
                              cout << "\t\t\t\t\talign done.\n"; //TESTING
#endif
                              // check alignment satisfies all thresholds
                              if ( result != 0 )
                              {
                                if ( result->score1 > minimal_score[index_num] ) aligned = true;
                              }
                              
                              // STEP 8: alignment succeeded, output (--all) or store (--best)
                              // the alignment and go to next alignment
                              if ( aligned )
                              {
#ifdef debug_align
                                cout << "\t\t\t\t\tscore = " << result->score1 << endl;
                                cout << "\t\t\t\t\taligned!\n"; //TESTING
#endif
#pragma omp critical
                                {
                                  // read has not been yet mapped, set bit to true for this read
                                  // (this is the _only_ place where read_hits must be modified)
                                  if ( !read_hits[readn] )
                                  {
                                    read_hits[readn].flip();
                                    total_reads_mapped++;
                                    reads_matched_per_db[index_num]++;
                                  }
                                  
                                  // output (sam or blast-like) or update alignment information
                                  bool strand = false;
                                  if ( forward_gv ) strand = true;
                                  else strand = false;
                                  
                                  // add the offset calculated by the LCS (from the beginning of the sequence)
                                  // to the offset computed by SW alignment
                                  result->ref_begin1 += (align_ref_start-head);
                                  result->ref_end1 += (align_ref_start-head);
                                  result->read_begin1 += align_que_start;
                                  result->read_end1 += align_que_start;
                                  result->readlen = readlen;
                                                                    
                                  // update best alignment
                                  if ( min_lis_gv > -1 )
                                  {
                                    result->index_num = index_num;
                                    result->ref_seq = max_seq;
                                    result->part = part;
                                    result->strand = strand;

                                    map<uint32_t, alignment_struct>::iterator alignment = read_hits_align_info.find(readn);
#ifdef DEBUG_BEST_N
                                    cout << "\nreadn = " << readn << endl;
                                    cout << "max_seq = " << (2*max_seq) << endl;
                                    cout << "result->score1 = " << result->score1 << endl;
                                    cout << "index_num = " << result->index_num << endl;
                                    cout << " part = " << result->part << endl;
                                    cout << " read = ";
                                    char* ttc = reads[readn-1];
                                    while ( *ttc != '\n' ) cout << (char)*ttc++;
                                    cout << endl;
                                    cout << " refseq = ";
                                    char* ttp = reference_seq[(2*result->ref_seq)]+1;
                                    while ( *ttp != '\n' ) cout << (char)*ttp++;
                                    cout << endl;
#endif
                                    // an alignment for this read already exists
                                    if ( alignment != read_hits_align_info.end() )
                                    {
#ifdef DEBUG_BEST_N
                                      cout << "alignment already exists.\n";
#endif
                                      uint32_t smallest_score_index = alignment->second.min_index;
                                      uint32_t highest_score_index = alignment->second.max_index;
				      uint32_t array_size = alignment->second.size;
				      uint32_t array_max_size = alignment->second.max_size;
#ifdef DEBUG_BEST_N
                                      cout << "smallest_score_index = " << smallest_score_index << endl; 
                                      cout << "highest_score_index = " << highest_score_index << endl;
				      cout << "array_size = " << array_size << endl;
				      cout << "array_max_size = " << array_max_size << endl;
				      cout << "num_best_hits_gv = " << num_best_hits_gv << endl;
#endif
                                      // number of alignments stored per read < num_best_hits_gv, 
                                      // add alignment to array without comparison to other members
				      // of array
				      if ( (num_best_hits_gv == 0) || (array_size < num_best_hits_gv) )
				      {
					// number of alignments stored per read == maximum number of 
					// alignments allowed, resize array by another BEST_HITS_INCREMENT slots 
					if ( array_size == array_max_size )
					{
					  uint32_t new_array_max_size = 0;
					  if ( (num_best_hits_gv == 0) || (array_size + BEST_HITS_INCREMENT <= num_best_hits_gv) )
					    new_array_max_size = array_max_size + BEST_HITS_INCREMENT;
					  else
					    new_array_max_size = num_best_hits_gv;
#ifdef DEBUG_BEST_N
					  cout << "\t\tresize array.\n";
                                          cout << "\t\tnew_array_max_size =  " << new_array_max_size << endl;
#endif
					  s_align* bigger_alignment_array = new s_align[new_array_max_size]();
					  if ( bigger_alignment_array == NULL )
					  {
					    fprintf(stderr, "\t  %sERROR%s: could not allocate memory for "
						            "alignment storage (s_align* bigger_alignment_array "
						            "in paralleltraversal.cpp\n", "\033[0;31m", "\033[0m");
					    exit(EXIT_FAILURE);
					  }

					  // copy smaller array to larger memory slot
					  memcpy(bigger_alignment_array, alignment->second.ptr, sizeof(s_align)*array_max_size);
					  // delete smaller array memory
					  delete [] alignment->second.ptr;
					  // set pointer to new larger array
					  alignment->second.ptr = bigger_alignment_array;
					  // set new max_size
					  alignment->second.max_size = new_array_max_size;
					  // update maximum size of array for downstream computation
					  array_max_size = alignment->second.max_size;
					}
#ifdef DEBUG_BEST_N
                                        cout << "\tadd new alignment to empty array slot.\n";
					cout << "\tarray_max_size = " << array_max_size << endl;
#endif
					s_align *smallest_alignment = alignment->second.ptr + array_size;
                                        *smallest_alignment = *result;
					// increment size of array
					alignment->second.size++;
					array_size++;

                                        // all slots have been filled, find slot with smallest
                                        // alignment score and set the smallest_score_index
					// (this is not done when num_best_hits_gv == 0 since
					// we want to output all alignments for some --min_lis)
                                        if ( array_size == num_best_hits_gv )
                                        {
#ifdef DEBUG_BEST_N
                                          cout << "\t\tfind new smallest_score_index of " << num_best_hits_gv << " slots.\n";
#endif
                                          uint32_t smallest_score = 1000000;
                                          s_align *this_alignment = alignment->second.ptr;
                                          for ( int p = 0 ; p < num_best_hits_gv; p++ )
                                          {
                                            if ( this_alignment[p].score1 < smallest_score )
                                            {
                                              smallest_score = this_alignment[p].score1;
                                              smallest_score_index = p;
                                            }
                                          }
#ifdef DEBUG_BEST_N
                                          cout << "\t\tnew smallest_score_index = " << smallest_score_index << endl;
#endif
                                          alignment->second.min_index = smallest_score_index;
                                        }
	                                // update the index position of the first occurrence of the
					// highest alignment score
					if ( result->score1 > (alignment->second.ptr + highest_score_index)->score1 )
					  alignment->second.max_index = array_size - 1;

#ifdef DEBUG_BEST_N
					cout << "\tmax_index = " << alignment->second.max_index << endl;
#endif
                                        // the maximum possible score for this read has been found
                                        if ( result->score1 == max_SW_score ) read_max_SW_score[readn]++;                                  
#ifdef DEBUG_BEST_N
                                        cout << "\tfree result.\n";
#endif
                                        // free result
                                        free(result);
                                        result = NULL;                                          
                                      }//~if (array_size < num_best_hits_gv)

                                      // all num_best_hits_gv slots have been filled,
                                      // replace the alignment having the lowest score
                                      else if ( result->score1 > (alignment->second.ptr + smallest_score_index)->score1 )
                                      {
                                        // get the alignment with the smallest score
                                        s_align *smallest_alignment = (alignment->second.ptr) + smallest_score_index;
#ifdef DEBUG_BEST_N
                                        cout << "\treplace alignment with higher score.\n";
                                        cout << "\tsmallest_alignment->score1 = " << smallest_alignment->score1 << endl;
                                        cout << "\treplace alignment in an existing slot.\n";
#endif
					// update max_index to the position of the first occurrence
					// of the highest scoring alignment
					if ( result->score1 > (alignment->second.ptr + highest_score_index)->score1 )
					  alignment->second.max_index = smallest_score_index;

                                        // decrement number of reads mapped to database
                                        // with lower score
                                        reads_matched_per_db[smallest_alignment->index_num]--;
                                        // increment number of reads mapped to database
                                        // with higher score
                                        reads_matched_per_db[index_num]++;
                                        // free the old cigar
                                        free(smallest_alignment->cigar);
                                        smallest_alignment->cigar = NULL;
                                        *smallest_alignment = *result;
                                        
                                        // find the new smallest_score_index
                                        uint32_t smallest_score = 1000000;
                                        s_align *this_alignment = alignment->second.ptr;
                                        for ( int p = 0 ; p < num_best_hits_gv; p++ )
                                        {
                                          if ( this_alignment[p].score1 < smallest_score )
                                          {
                                            smallest_score = this_alignment[p].score1;
                                            smallest_score_index = p;
                                          }
                                        }                                                               
#ifdef DEBUG_BEST_N
                                        cout << "\nnew smallest_score_index = " << smallest_score_index << endl;
#endif                                 
                                        alignment->second.min_index = smallest_score_index;
                                        
                                        // the maximum possible score for this read has been found
                                        if ( result->score1 == max_SW_score ) read_max_SW_score[readn]++;
                                        
                                        // free result, except the cigar (now new cigar)
                                        free(result);
                                        result = NULL;
                                      }
                                      else
                                      {
                                        // new alignment has a lower score, destroy it
                                        if (result != NULL) align_destroy(&result);
                                      }
                                    }
                                    // an alignment for this read doesn't exist, add the first alignment
                                    else
                                    {
#ifdef DEBUG_BEST_N
                                      cout << "add first alignment.\n";
#endif                          
				      // maximum size of s_align array
				      uint32_t max_size = 0;
				      
                                      // create new instance of alignments
				      if ( (num_best_hits_gv > 0) && (num_best_hits_gv < BEST_HITS_INCREMENT+1) )
					max_size = num_best_hits_gv;
				      else max_size = BEST_HITS_INCREMENT;

				      s_align *new_alignment = new s_align[max_size]();
#ifdef DEBUG_BEST_N
				      cout << "memory allocated for alignment array = "
					   << sizeof(s_align)*max_size << " bytes" << endl;
#endif
				      if ( new_alignment == NULL )
                                      {
                                        fprintf(stderr,"\n  %sERROR%s: could not allocate memory for alignment "
                                                       "storage (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
                                        exit(EXIT_FAILURE);
                                      }
                                      else
                                      {
                                        new_alignment[0] = *result;
                                        alignment_struct tmp (max_size, 1, 0, 0, new_alignment);
                                        read_hits_align_info.insert( pair<uint32_t, alignment_struct > (readn,tmp) );
                                      }
                                                                            
                                      // the maximum possible score for this read has been found
                                      if ( result->score1 == max_SW_score ) read_max_SW_score[readn]++;
                                      
                                      // free result, except the cigar
                                      free(result);
                                      result = NULL;                
#ifdef DEBUG_BEST_N
				      cout << "max_size = " << max_size << endl;
				      cout << "size = " << 1 << endl;
				      cout << "min_index = " << 0 << endl;
				      cout << "max_index = " << 0 << endl;
                                      cout << "added.\n";
#endif
                                    }                                      
                                  }
                                  // output the Nth alignment (set by --num_alignments [INT] parameter)
                                  else if ( num_alignments_gv > -1 )
                                  {
                                    // update number of alignments to output per read
                                    if ( num_alignments_gv > 0 ) num_alignments_x[readn]--;
                                    
                                    // get the edit distance between reference and read (serves for
                                    // SAM output and computing %id and %query coverage)
                                    double id = 0;
                                    char to_char[5] = {'A','C','G','T','N'};
                                    char* ref_seq_ptr = reference_seq[(2*(int)max_seq)+1];
                                    char* read_seq_ptr = myread;
                                    int32_t qb = result->ref_begin1;
                                    int32_t pb = result->read_begin1;
                                    uint32_t mismatches = 0;
                                    uint32_t gaps = 0;
                                    
                                    for (uint32_t c2 = 0; c2 < result->cigarLen; ++c2) 
                                    {
                                      uint32_t letter = 0xf&*(result->cigar + c2);
                                      uint32_t length = (0xfffffff0&*(result->cigar + c2))>>4;
                                      if (letter == 0) 
                                      {
                                        for (int p = 0; p < length; ++p)
                                        {
                                          if ( (char)to_char[(int)*(ref_seq_ptr + qb)] != (char)to_char[(int)*(read_seq_ptr + pb)] ) ++mismatches;
                                          else ++id;
                                          ++qb;
                                          ++pb;
                                        }
                                      } else if (letter == 1) 
                                      {
                                        pb += length;
                                        gaps += length;
                                      } else 
                                      {
                                        qb += length;
                                        gaps += length;
                                      }
                                    }
                                                                        
                                    int32_t align_len = abs(result->read_end1+1 - result->read_begin1);
                                    int32_t total_pos = mismatches+gaps+id;
                                    stringstream ss;
                                    ss.precision(3);
                                    ss << (double)id/total_pos << ' ' << (double)align_len/readlen;
                                    double align_id_round = 0.0;
                                    double align_cov_round = 0.0;
                                    ss >> align_id_round >> align_cov_round;
//#define debug_id_cov
#ifdef debug_id_cov
                                    cout << "read tag: ";
                                    char* tt = reads[readn-1];
                                    while (*tt != '\n' ) cout << (char)*tt++;
                                    cout << endl;
                                    
                                    cout << "ref tag: ";
                                    tt = reference_seq[(2*(int)max_seq)];
                                    while (*tt != '\n' ) cout << (char)*tt++;
                                    cout << endl;
                                    
                                    cout << "align_len = " << align_len << endl;
                                    cout << "id = " << id << endl;
                                    cout << "Score = " << result->score1 << endl;
                                    cout << "%id = " << (double)id/align_len << endl;
                                    cout << "%cov = " << (double)align_len/readlen << endl;
                                    cout << "align_cov = " << (double)align_cov << endl;
                                    cout << "align_id = " << (double)align_id << endl << endl;
                                    cout << "align_id_round = " << (double)align_id_round << endl;
                                    cout << "align_cov_round = " << (double)align_cov_round << endl;
#endif                                                        
                                    // the alignment passed the %id and %query coverage threshold
                                    // output it (SAM, BLAST and FASTA/Q)
                                    if ( ( align_id_round >= align_id ) &&
                                         ( align_cov_round >= align_cov ) &&
                                          read_to_count )
                                    {
                                      total_reads_mapped_cov++;
                                      read_to_count = false;

                                      // do not output read for de novo OTU clustering
                                      // (it passed the %id/coverage thersholds)
                                      if ( de_novo_otu_gv ) read_hits_denovo[readn].flip();
                                    }
                                                                        
                                    // quality for FASTQ
                                    char* read_qual = NULL;
                                    if ( filesig == '@' )
                                    {
                                      // forward
                                      if ( strand )
                                      {
                                        // if second part of split-read or last read in file
                                        if ( ((readn == 3) && (file_s > 0)) || (readn >= (strs-2)) )
                                        {
                                          read_qual = reads[readn];
                                          int8_t numnewlines = 0;
                                          while ( numnewlines<2 ) { if ( *read_qual++ == '\n' ) numnewlines++; }
                                        }
                                        else read_qual = reads[readn+1]-readlen-1;
                                      }
                                      // reverse-complement
                                      else
                                      {
                                        // if second part of split read or last read in file
                                        if ( ((readn == 3) && (file_s > 0)) || (readn >= (strs-2)) )
                                        {
                                          read_qual = reads[readn];
                                          
                                          // last file section
                                          if ( file_s == file_sections-1 )
                                          {
#ifdef debug_mmap
                                            if (readn >= (strs-2)) cout << "get quality for last (reverse) read in last file section\n"; //TESTING
#endif
                                            while ( *read_qual != '\0' ) read_qual++;
                                            if ( *(read_qual-3) == '\n') read_qual--;
                                            read_qual-=2; //account for '\n\0'
                                          }
                                          // file section > 0 and < last file section
                                          else
                                          {
#ifdef debug_mmap
                                            if (readn >= (strs-2)) cout << "get quality for last (reverse) read in (not last) file section\n"; //TESTING
#endif
                                            while ( read_qual != finalnt ) read_qual++;
                                            read_qual--;
                                          }
                                        }
                                        else read_qual = reads[readn+1]-2;
                                      }
                                    }
                                                             
                                    if ( blastout_gv )
                                    {
                                        uint32_t bitscore = (uint32_t)((float)((gumbel[index_num].first)*(result->score1) - log((gumbel[index_num].second)))/(float)log(2));
                                        double evalue_score = (double)(gumbel[index_num].second)*full_ref[index_num]*full_read[index_num]*pow(EXP,(-(gumbel[index_num].first)*result->score1));
                                        
                                        report_blast (acceptedblast, //blast output file
                                                      result, //SW alignment cigar
                                                      reads[readn-1]+1, //read name
                                                      &myread[0], //read sequence
                                                      read_qual, //read quality
                                                      reference_seq[(2*(int)max_seq)]+1, //reference name
                                                      reference_seq[(2*(int)max_seq)+1], //reference sequence
                                                      evalue_score, //e-value score (only for blast)
                                                      readlen,
                                                      bitscore,
                                                      strand, //forward or reverse complement
                                                      (double)id/total_pos, // %id
                                                      (double)align_len/readlen, // %query coverage
                                                      mismatches, //mismatches
                                                      gaps);
                                    }
                                                                        
                                    if ( samout_gv )
                                    {
                                        report_sam (acceptedsam, //sam output file
                                                    result, //SW alignment cigar
                                                    reads[readn-1]+1, //read name
                                                    &myread[0], //read sequence
                                                    read_qual, //read quality
                                                    reference_seq[(2*(int)max_seq)]+1, //reference name
                                                    reference_seq[(2*(int)max_seq)+1], //reference sequence
                                                    readlen,
                                                    strand, //forward or reverse complement
                                                    mismatches+gaps
                                                    );
                                    }
                                                                        
                                                                        
                                    // free alignment info
                                    if ( result != 0 ) align_destroy(&result);
                                  }//~if output all alignments
                                                                    
                                  // continue to next read (do not need to collect more seeds using another pass)
                                  search = false;
                                                                    
                                }//~pragma omp critical
                                
                                // maximum score possible for the read has been reached,
                                // stop searching for further matches
                                if ( (num_best_hits_gv != 0) && (read_max_SW_score[readn] == num_best_hits_gv) ) break;
                                
                                // stop search after the first num_alignments_gv alignments
                                // for this read
                                if ( num_alignments_gv > 0 )
                                {
                                  // go to next read (do not search for further alignments)
                                  if ( num_alignments_x[readn] <= 0 ) break;
                                }                                                         
                              }//~if read aligned
                              // the read did not align, free alignment info
                              else
                              {
#ifdef debug_align
                                  cout << "\t\t\t\tnot aligned \n"; //TESTING
#endif         
                                  if ( result != 0 ) align_destroy(&result);                            
                              }
                              
                            }//~if LIS long enough
                                                        
#ifdef HEURISTIC1_OFF
                          } while ( (it3 == hits_on_genome.end()) && (++list_n < list.size()) );
#endif                                              
                        }//~if enough window hits
                                                
                        pop:
                                                
                        // get the next candidate reference position 
                        if ( !vi_read.empty() ) vi_read.pop_front();
                                                
                        if ( vi_read.empty() )
                        {
                          if ( it3 != hits_on_genome.end() )
                            begin = it3->first;
                          else break;
                        }
                        else
                        {
                          begin = (vi_read.front()).first;
                        }                      
                      }//~for all reference sequence length                   
                    }//~for all of the reference sequence candidates
                  }//if ( readhitf || readhitr > ratio )
                                    
                                    
                  // edit the read back onto 0-3 alphabet to search for seeds 
                  // (seed search cannot proceed on 0-4 alphabet)
                  if ( (size_ambiguous_nt!=0) && read_edited )
                  {
                    if ( !forward_gv )
                    {
                      for ( uint32_t p = 0; p < size_ambiguous_nt; p++ )
                      {
                        myread[(readlen-ambiguous_nt[p])-1] = 0;
                      }
                    }
                    else
                    {
                      for ( uint32_t p = 0; p < size_ambiguous_nt; p++ )
                      {
                        myread[ambiguous_nt[p]] = 0;
                      }
                    }
                  }
                                    
                  // the read was not accepted at current window skip length,
                  // decrease the window skip length
                  if ( search )
                  {
                    // last (3rd) Pass has been made
                    if ( pass_n == 2 ) search = false;
                    else
                    {
                      // the next interval size equals to the current one, skip it
                      while ( (pass_n < 3) && 
                              (skiplengths[index_num][pass_n] == skiplengths[index_num][pass_n+1]) ) ++pass_n;
                      if ( ++pass_n > 2 ) search = false;
                      // set interval skip length for next Pass
                      else windowshift = skiplengths[index_num][pass_n];
                    }
                  }
                                    
                  // do not offset final window on read
                  break;
                }//~( win_num == NUMWIN-1 )
                                
                read_index+=windowshift;               
              }//~for (each window)                
            //~while all three window skip lengths have not been tested, or a match has not been found
            } while ( search );
                        
            // the read didn't align (for --num_alignments [INT] option),
            // output null alignment string
            if ( !read_hits[readn] && !forward_gv && (num_alignments_gv > -1) )
            {
              // do not output read for de novo OTU clustering
              // (it did not pass the E-value threshold)
              if ( de_novo_otu_gv && read_hits_denovo[readn] ) read_hits_denovo[readn].flip();

              // output null alignment string
              if ( print_all_reads_gv )
              {
#pragma omp critical
                {
                  s_align* null_alignment = NULL;
                  if ( blastout_gv && blast_tabular )
                  {
                    report_blast (acceptedblast, // blast output file
                                  null_alignment, // SW alignment cigar
                                  reads[readn-1]+1, //read name
                                  0, // read sequence (in integer format)
                                  0, // read quality
                                  0, // reference name
                                  0, // reference sequence
                                  0, // e-value score
                                  0, // read length (to compute the masked regions)
                                  0, // bitscore
                                  0, // forward or reverse strand
                                  0, // %id
                                  0, // %query coverage
                                  0, // number of mismatches
                                  0); // number of gaps
                  }
                  if ( samout_gv )
                  {
                    report_sam (acceptedsam, // sam output file
                                null_alignment, // SW alignment cigar
                                reads[readn-1]+1, // read name
                                0, // read sequence (in integer format)
                                0, // read quality
                                0, // reference name
                                0, // reference sequence
                                0, // read length (to compute the masked regions)
                                0, // forward or reverse strand
                                0); // edit distance
                  }
                }//~if print_all_reads_gv
              }// allow writing to file 1 thread at a time
            }//~if read didn't align                  
          }//~pragma omp for (each read)
#pragma omp barrier
                    
          // search the reverse strand (default)
          forward_gv = false;
          
        }//for forward and/or reverse strands
        TIME(f);
        eprintf(" done [%.2f sec]\n", (f-s) );
        
        eprintf("    Freeing index ... ");
        
        
        TIME(s);
#pragma omp master
        {
          // free the positions table
          if ( positions_tbl != NULL )
          {
            for ( uint32_t i = 0; i < number_elements; i++ )
            {
                if ( positions_tbl[i].arr != NULL )
                {
                  delete [] positions_tbl[i].arr;
                  positions_tbl[i].arr = NULL;
                }
            }
            delete [] positions_tbl;
            positions_tbl = NULL;
          }
          
          // free reference sequences loaded into memory
          if ( buffer != NULL )
          {
            delete [] buffer;
            buffer = NULL;
          }
          if ( reference_seq != NULL )
          {
            delete [] reference_seq;
            reference_seq = NULL;
          }
          if ( reference_seq_len != NULL )
          {
            delete [] reference_seq_len;
            reference_seq_len = NULL;
          }
                    
          // free 9-mer look-up tables and mini-burst tries
          if ( lookup_tbl != NULL )
          {
            for ( uint32_t i = 0; i < (1<<lnwin[index_num]); i++ )
            {
              if (lookup_tbl[i].trie_F != NULL )
              {
                delete [] lookup_tbl[i].trie_F;
                lookup_tbl[i].trie_F = NULL;
                lookup_tbl[i].trie_R = NULL;
              }
              if ( lookup_tbl[i].trie_R != NULL )
              {
                delete [] lookup_tbl[i].trie_R;
                lookup_tbl[i].trie_R = NULL;
              }
            }
            delete [] lookup_tbl;
            lookup_tbl = NULL;
          }
        }
        TIME(f);
        
        eprintf(" done [%.2f sec]\n", (f-s));
        
        // increment the index part to next file
        prt_str.str("");
        prt_str << (part+1);
        part_str = prt_str.str();

        // clear array holding number of reference sequences to
        // search per read prior to choosing best alignment
        if ( best_x != NULL ) delete [] best_x;
                
      }//~for all parts of the index            
    }//~for all indexes (provided as a list using option --ref)
        
    // clear array holding number of alignments output for each read
    if ( num_alignments_x != NULL ) delete [] num_alignments_x;

    delete [] read_max_SW_score;
    
    eprintf("    Total number of reads mapped (incl. all reads file sections searched): %u\n",total_reads_mapped);
        
    // filter the sequences by %id and %query coverage, output them if --best INT
    if ( min_lis_gv > -1 )
    {
      if ( samout_gv || blastout_gv ) eprintf("    Writing alignments ... ");

#ifdef debug_mmap
      cout << "total index_num = " << myfiles.size() << endl;
      for ( uint32_t index_num = 0; index_num < myfiles.size(); index_num++ )
          cout << "\t parts = " << num_index_parts[index_num] << endl;
#endif            
      TIME(s);
      // loop through all databases
      for ( uint32_t index_num = 0; index_num < myfiles.size(); index_num++ )
      {
#ifdef debug_output
        cout << "index_num = " << index_num << endl;
#endif
        // loop through each section of a database
        for ( uint16_t part = 0; part < num_index_parts[index_num]; part++ )
        {
          uint64_t seq_part_size = index_parts_stats_vec[index_num][part].seq_part_size;
          uint32_t numseq_part = index_parts_stats_vec[index_num][part].numseq_part;
          uint64_t start_part = index_parts_stats_vec[index_num][part].start_part;
#ifdef debug_output
          cout << "part = " << part << endl;
          cout << "seq_part_size = " << seq_part_size << endl;
          cout << "numseq_part = " << numseq_part << endl;
          cout << "start_part = " << start_part << endl;
#endif                    
          // block of memory to hold all ids + reference sequences
          char* buffer = new char[(seq_part_size+1)]();
          if ( buffer == NULL )
          {
            fprintf(stderr,"  %sERROR%s: could not allocate memory for reference sequence buffer (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
              
          // pointer to the start of every sequence in the buffer
          char** reference_seq = new char*[(numseq_part<<1)]();
          if ( reference_seq == NULL )
          {
            fprintf(stderr,"  %sERROR%s: could not allocate memory for reference_seq (paralleltraversal.cpp)\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
                    
          // length of every sequence in the buffer (is not required here)
          uint32_t* reference_seq_len = NULL;
          
          load_ref((char*)(myfiles[index_num].first).c_str(),buffer,reference_seq,reference_seq_len,seq_part_size,numseq_part,start_part,0);
                
          // run through all the reads, output those which aligned
          for ( uint32_t readn = 1; readn < strs; readn+=2 )
          {
            map<uint32_t, alignment_struct >::iterator alignment = read_hits_align_info.find(readn);

            // this read does not have any alignment
            if ( alignment == read_hits_align_info.end() )
            {
              // (output NULL alignment only once)
              if ( (index_num == 0) && (part == 0) )
              {
                // do not output this read for de novo clustering
                // (it did not pass the E-value threshold)
                if ( de_novo_otu_gv && read_hits_denovo[readn] ) read_hits_denovo[readn].flip();

                // output null string for read alignment
                if ( print_all_reads_gv )
                {
                  s_align* null_alignment = NULL;
                  if ( blastout_gv && blast_tabular )
                  {
                    report_blast (acceptedblast, // blast output file
                                  null_alignment, // SW alignment cigar
                                  reads[readn-1]+1, //read name
                                  0, // read sequence (in integer format)
                                  0, // read quality
                                  0, // reference name
                                  0, // reference sequence
                                  0, // e-value score
                                  0, // read length (to compute the masked regions)
                                  0, // bitscore
                                  0, // forward or reverse strand
                                  0, // %id
                                  0, // %query coverage
                                  0, // number of mismatches
                                  0); // number of gaps
                  }                                
                  if ( samout_gv )
                  {
                    report_sam (acceptedsam, // sam output file
                                null_alignment, // SW alignment cigar
                                reads[readn-1]+1, // read name
                                0, // read sequence (in integer format)
                                0, // read quality
                                0, // reference name
                                0, // reference sequence
                                0, // read length (to compute the masked regions)
                                0, // forward or reverse strand
                                0); // edit distance
                  }
                }
              }
              // go to next read
              continue;
            }
#ifdef debug_output
            cout << "readn = " << readn << endl;
#endif
            // get all alignments for this read
            //s_align* ptr_alignment = alignment->second.second;
            s_align* ptr_alignment = alignment->second.ptr;
            if ( ptr_alignment == NULL )
            {
              fprintf(stderr,"  ERROR: s_align* ptr_alignment == NULL, this should not be possible.\n");
              exit(EXIT_FAILURE);
            }
#ifdef debug_output
            cout << "ptr_alignment->index_num = " << (uint16_t)ptr_alignment->index_num << endl; 
            cout << "ptr_alignment->score1 = " << (int32_t)ptr_alignment->score1 << endl;                 
            cout << "ptr_alignment->part = " << (uint16_t)ptr_alignment->part << endl;
#endif
            // OTU-map: index of alignment holding maximum SW score
            int index_max_score = alignment->second.max_index;

            // loop through all of the best alignments for this read
            for ( int p = 0; p < alignment->second.size; p++ )
            {
#ifdef debug_output
              cout << "best_hit = " << p << endl;
#endif
              // continue loop if the reference sequence in this alignment
              // belongs to the database section currently loaded into RAM
              if ( (ptr_alignment->index_num == index_num) && (ptr_alignment->part == part) )
              {
                // format read & get read length
                char myread[READLEN];
                uint32_t readlen = ptr_alignment->readlen;
                
                // format forward read from char to int
                if ( ptr_alignment->strand )
                {
                  format_forward(reads[readn],&myread[0],filesig);
                }
                // format reverse-complement read from char to int
                else
                {
                  char* end_read = NULL;
                  // FASTA
                  if ( filesig == '>' )
                  {
                    // if split-read (for file_s > 0)
                    // note: the split-read contains two reads (a pair, in case reads are paired) that
                    //   are located at another address space than the contiguous memory mapped
                    //   section of reads
                    if ((readn < 4) && (file_s > 0) )
                    {
#ifdef debug_output
                      cout << "process split-read" << endl; //TESTING
#endif
                      end_read = reads[readn];
                      while (*end_read != '\0') end_read++;
                      if ( *end_read == '\0' ) end_read--; //account for '\0'
                      if ( *end_read == '\n' ) end_read--; //account for '\n'
                    }
                    // last read in the partial file section
                    else if ( readn >= (strs-2) )
                    {
                      end_read = reads[readn];
                      
                      // if processing last file section, the final read will end with '\0'
                      if ( file_s == file_sections-1 )
                      {
#ifdef debug_output
                        cout << "process last read in final file section .." << endl; //TESTING
#endif
                        while ( *end_read++ != '\0' );
                        // back-track 3 nucleotides (\n\0_)
                        if ( *(end_read-2) == '\n' ) end_read--;
                        // back-track 2 nucleotides (\0_)
                        end_read-=2;
                      }
                      // if processing a file section > 0 and < last file section, the final read will end with '>' (beginning of split-read)
                      else
                      {
#ifdef debug_output
                        cout << "process last read in file section > 0 and < final file section .." << endl; //TESTING
#endif
                        while ( *end_read++ != '>' );
                        // back-track 3 nucleotides (\n>_)
                        end_read-=3;
                      }
                    }
                    // all other reads
                    else end_read = reads[readn+1]-2;
                  }
                  // FASTQ
                  else end_read = reads[readn]+readlen-1;
                  
                  format_rev(reads[readn],end_read,&myread[0],filesig);
                }//~reverse-complement read
                                           
                // get the edit distance between reference and read
                uint32_t ref_seq = ptr_alignment->ref_seq;
                double id = 0;
                char to_char[5] = {'A','C','G','T','N'};
                char* ref_seq_ptr = reference_seq[(2*ref_seq)+1];
                char* read_seq_ptr = myread;
                int32_t qb = ptr_alignment->ref_begin1;
                int32_t pb = ptr_alignment->read_begin1;
                uint32_t mismatches = 0;
                uint32_t gaps = 0;
#ifdef debug_output
                cout << "ref_seq = " << ref_seq << endl;
                cout << "ptr_alignment->ref_begin1 = " << ptr_alignment->ref_begin1 << endl;
                cout << "ptr_alignment->read_begin1 = " << ptr_alignment->read_begin1 << endl;
#endif                           
                for (uint32_t c2 = 0; c2 < ptr_alignment->cigarLen; ++c2)
                {
                  uint32_t letter = 0xf&*(ptr_alignment->cigar + c2);
                  uint32_t length = (0xfffffff0&*(ptr_alignment->cigar + c2))>>4;
                  if (letter == 0) 
                  {
                    for (int u = 0; u < length; ++u)
                    {
                      if ( (char)to_char[(int)*(ref_seq_ptr + qb)] != (char)to_char[(int)*(read_seq_ptr + pb)] ) ++mismatches;
                      else ++id;
                      ++qb;
                      ++pb;
                    }
                  } 
                  else if (letter == 1) 
                  {
                    pb += length;
                    gaps += length;
                  } 
                  else
                  {
                    qb += length;
                    gaps += length;
                  }
                }
                                
                int32_t align_len = abs(ptr_alignment->read_end1+1 - ptr_alignment->read_begin1);
                int32_t total_pos = mismatches+gaps+id;
                stringstream ss;
                ss.precision(3);
                ss << (double)id/total_pos << ' ' << (double)align_len/readlen;
                double align_id_round = 0.0;
                double align_cov_round = 0.0;
                ss >> align_id_round >> align_cov_round;
                                     
//#define debug_id_cov
#ifdef debug_id_cov
                cout << "read tag: ";
                char* tt = reads[readn-1];
                while (*tt != '\n' ) cout << (char)*tt++;
                cout << endl;                          
                cout << "align_len = " << align_len << endl;
                cout << "total_pos = " << total_pos << endl;
                cout << "id = " << id << endl;
                cout << "Score = " << ptr_alignment->score1 << endl;
                cout << "%id = " << (double)id/total_pos << endl;
                cout << "%cov = " << (double)align_len/readlen << endl;
                cout << "align_cov = " << (double)align_cov << endl;
                cout << "align_id = " << (double)align_id << endl;
                cout << "align_id_round = " << (double)align_id_round << endl;
                cout << "align_cov_round = " << (double)align_cov_round << endl;
#endif                              
                // alignment with the highest SW score passed
                // %id and %coverage thresholds
                if ( ( p == index_max_score ) &&
                     ( align_id_round >= align_id ) && 
                     ( align_cov_round >= align_cov ) ) 
                {
                  // increment number of reads passing identity
                  // and coverage threshold
                  total_reads_mapped_cov++;

                  // do not output read for de novo OTU construction
                  // (it passed the %id/coverage thresholds)
                  if ( de_novo_otu_gv && read_hits_denovo[readn] ) read_hits_denovo[readn].flip();
                  
                  // fill OTU map with highest-scoring alignment for the read
                  if ( otumapout_gv )
                  {
                    // reference sequence identifier for mapped read
                    char ref_seq_arr[4000] = "";
                    char* ref_seq_arr_ptr = ref_seq_arr;
                    char* ref_seq_id_ptr = reference_seq[(2*ref_seq)]+1;
                    while ( (*ref_seq_id_ptr != ' ') && (*ref_seq_id_ptr != '\n') ) *ref_seq_arr_ptr++ = *ref_seq_id_ptr++;
                    string ref_seq_str = ref_seq_arr;
                    
                    // read identifier
                    char read_seq_arr[4000] = "";
                    char* read_seq_arr_ptr = read_seq_arr;
                    char* read_seq_id_ptr = reads[readn-1]+1;
                    while ( (*read_seq_id_ptr != ' ') && (*read_seq_id_ptr != '\n') ) *read_seq_arr_ptr++ = *read_seq_id_ptr++;
                    string read_seq_str = read_seq_arr;
                    
                    otu_map[ref_seq_str].push_back(read_seq_str);
                  }
                }
                                              
                // output alignment to SAM or Blast-like formats
                if ( samout_gv || blastout_gv )
                {
                  // quality for FASTQ
                  char* read_qual = NULL;
                  if ( filesig == '@' )
                  {
                    // forward
                    if ( ptr_alignment->strand )
                    {
                      // if second part of split read or last read in file
                      if ( ((readn == 3) && (file_s > 0)) || (readn >= (strs-2)) )
                      {
#ifdef debug_output
                        if (readn >= (strs-2)) cout << "get quality for last (forward) read in file section\n"; //TESTING
#endif
                        read_qual = reads[readn];
                        int8_t numnewlines = 0;
                        while ( numnewlines<2 ) { if ( *read_qual++ == '\n' ) numnewlines++; }
                      }
                      else read_qual = reads[readn+1]-readlen-1;
                    }
                    // reverse-complement
                    else
                    {
                      // if second part of split read or last read in file
                      if ( ((readn == 3)&&(file_s > 0)) || (readn >= (strs-2)) )
                      {
                        read_qual = reads[readn];
                          
                        // last file section
                        if ( file_s == file_sections-1 )
                        {
#ifdef debug_output
                          if (readn >= (strs-2)) cout << "get quality for last (reverse) read in last file section\n"; //TESTING
#endif
                          while ( *read_qual != '\0' ) read_qual++;
                          if ( *(read_qual-3) == '\n') read_qual--;
                          // account for '\n\0'
                          read_qual-=2;
                        }
                        // file section > 0 and < last file section
                        else
                        {
#ifdef debug_output
                          if (readn >= (strs-2)) cout << "get quality for last (reverse) read in (not last) file section\n"; //TESTING
#endif
                          while ( read_qual != finalnt ) read_qual++;
                          read_qual--;
                        }          
                      }
                      else read_qual = reads[readn+1]-2;
                    }
                  }//~if filesig == '@'
                        
                  if ( blastout_gv )
                  {
                    uint32_t bitscore = (uint32_t)((float)((gumbel[index_num].first)*(ptr_alignment->score1) - log(gumbel[index_num].second))/(float)log(2));
                    double evalue_score = (double)(gumbel[index_num].second)*full_ref[index_num]*full_read[index_num]*pow(EXP,(-(gumbel[index_num].first)*ptr_alignment->score1));
                        
                    report_blast (acceptedblast, //blast output file
                                  ptr_alignment, //SW alignment cigar
                                  reads[readn-1]+1, //read name
                                  myread, //read sequence (in integer format)
                                  read_qual, //read quality
                                  reference_seq[(2*ref_seq)]+1, //reference name
                                  reference_seq[(2*ref_seq)+1], //reference sequence
                                  evalue_score, //e-value score
                                  readlen, //read length (to compute the masked regions)
                                  bitscore,
                                  ptr_alignment->strand,
                                  (double)id/total_pos, // %id
                                  (double)align_len/readlen, // %query coverage
                                  mismatches,
                                  gaps
                                  );
                  }
                        
                  if ( samout_gv )
                  {
                    report_sam (acceptedsam, //sam output file
                                ptr_alignment, //SW alignment cigar
                                reads[readn-1]+1, //read name
                                myread,//read sequence (in integer format)
                                read_qual, //read quality
                                reference_seq[(2*ref_seq)]+1, //reference name
                                reference_seq[(2*ref_seq)+1], //reference sequence
                                readlen, //read length (to compute the masked regions)
                                ptr_alignment->strand,
                                mismatches+gaps); //edit distance
                  }         
                }//~if (samout_gv || blastout_gv)
              }//~if alignment at current database and index part loaded in RAM
	      ptr_alignment++;
            }//~for all best alignments
          }//~for all reads
                        
          // free buffer
          if ( buffer != NULL )
          {
            delete [] buffer;
            buffer = NULL;
          }
          
          // free pointers to reference sequences
          if ( reference_seq != NULL )
          {
            delete [] reference_seq;
            reference_seq = NULL;
          }                  
        }//~for every database section
      }//~for every database

      // free alignment information for all aligned reads
      for ( uint32_t readn = 1; readn < strs; readn+=2 )
      {
        map<uint32_t, alignment_struct >::iterator alignment = read_hits_align_info.find(readn);

        // this read does not have any alignment
        if ( alignment != read_hits_align_info.end() )
        {
          s_align* ptr_alignment = alignment->second.ptr;
          for ( uint32_t p = 0; p < alignment->second.size; p++ )
          {
            free(ptr_alignment->cigar);
            ptr_alignment->cigar = NULL;

            if ( p+1 < num_best_hits_gv )
            {
              ptr_alignment++;
      
              // check whether an alignment exists
              if ( ptr_alignment->cigar == NULL ) break;
            }
          }

          // free memory for all alignments of this read
          delete [] alignment->second.ptr;
          alignment->second.ptr = NULL;
        }
      }
            
      TIME(f);
      if ( samout_gv || blastout_gv ) eprintf(" done [%.2f sec]\n", (f-s) );      
    }// if (min_lis_gv > -1)
        
    if ( align_cov || align_id )
    {
      eprintf("    Total number of reads mapped with");
      if ( align_id > 0 ) eprintf(" >= %.2lf identity,",align_id);
      if ( align_cov > 0 ) eprintf(" >= %.2lf query coverage",align_cov);
      eprintf(" (incl. all reads file sections searched): %u\n",total_reads_mapped_cov);
    }
        
    if ( blastout_gv )
    {
      if ( acceptedblast.is_open() ) acceptedblast.close();
      else
      {
        fprintf(stderr,"  %sERROR%s: file %s was not opened for writing.\n","\033[0;31m","\033[0m",acceptedstrings_blast);
        exit(EXIT_FAILURE);
      }
    }
    if ( samout_gv )
    {
      if ( acceptedsam.is_open() ) acceptedsam.close();
      else
      {
        fprintf(stderr,"  %sERROR%s: file %s was not opened for writing.\n","\033[0;31m","\033[0m",acceptedstrings_sam);
        exit(EXIT_FAILURE);
      }
    }
                
    // output aligned and non-aligned reads to FASTA/FASTQ file
    report_fasta (acceptedstrings,
                  ptr_filetype_or,
                  ptr_filetype_ar,
                  reads,strs,
                  read_hits,
                  file_s,
                  finalnt);
    
    // output aligned and non-aligned reads with < %id and
    // < %coverage to FASTA/FASTQ file for de novo analysis
    if ( de_novo_otu_gv )
    {
      // count number of reads output for de novo clustering
      for ( int d = 1; d < strs; d+=2 )
        if ( read_hits_denovo[d] ) total_reads_denovo_clustering++;

      report_denovo(denovo_otus_file,
                    reads,
                    strs,
                    read_hits_denovo,
                    file_s,
                    finalnt);
    }
        
    read_hits.clear();
    read_hits_denovo.clear();
    
    // free the split_read
    if ( split_read != NULL )
    {
      delete [] split_read;
      split_read = NULL;
      split_read_ptr = NULL;
    }
        
    // record the start of the split_read if it exists
    if ( file_s < file_sections-1 )
    {
      split_read = new char[(READLEN*2)];
      
      if ( split_read == NULL )
      {
        fprintf(stderr, "  %sERROR%s: could not allocate memory for the bridged read\n","\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
      }
      
      // compute the first half of the split_read
      char *start = &raw[partial_file_size]-reads_offset_e-1;
      char *end = &raw[partial_file_size];
      
      split_read_ptr = split_read;
        
#ifdef debug_mmap
      cout << "split read start: "; //TESTING
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
        
    }// (s < file_sections - 1)
    
    // free the mmap'd file section
    if ( munmap(raw, partial_file_size ) == -1 )
    {
      fprintf(stderr,"  %sERROR%s: Could not munmap file!\n","\033[0;31m","\033[0m");
      exit(EXIT_FAILURE);
    }
    
    offset_map+=map_size_gv;
    
    // last section of the full file, count strings until EOF is reached
    if ( ++file_s == file_sections - 1 ) partial_file_size = last_part_size;
    
    delete [] reads;
    reads = NULL;  
  }//~while ( file_s < file_sections )
    
  // output OTU map to file
  if ( otumapout_gv )
  {
    ofstream outfile (acceptedotumap_file,ios::app);
    
    map<string,vector<string> >::iterator otu_map_it;
    
    for ( otu_map_it = otu_map.begin(); otu_map_it != otu_map.end(); otu_map_it++ )
    {
      // output the ref ID
      outfile << otu_map_it->first;
      // output all reads mapping to ref ID
      for ( uint32_t i = 0; i < otu_map_it->second.size(); i++ ) outfile << "\t" << otu_map_it->second[i];
      outfile << "\n";
    }
    
    if ( outfile.is_open() ) outfile.close();
    else
    {
      fprintf(stderr,"  %sERROR%s: file %s was not opened for writing.\n","\033[0;31m","\033[0m",acceptedotumap_file);
      exit(EXIT_FAILURE);
    }
    
    // free memory for OTU mapping file
    if ( acceptedotumap_file != NULL )
    {
      delete [] acceptedotumap_file;
      acceptedotumap_file = NULL;
    }
  }
    
  // close the reads file descriptor
  close(fd);
    
  free(mat);
  mat = NULL;
    
  // create a bilan (log file)
  if ( (ptr_filetype_ar != NULL) && logout_gv )
  {
    FILE* bilan = fopen(logoutfile,"a");
    if ( bilan == NULL )
    {
      fprintf(stderr,"  %sERROR%s: could not open file %s \n","\033[0;31m","\033[0m",logoutfile);
      exit(EXIT_FAILURE);
    }
    
    // output total number of reads
    fprintf(bilan," Results:\n");
    fprintf(bilan,"    Total reads = %u\n", number_total_read);
    if ( de_novo_otu_gv )
    {
      fprintf(bilan,"    Total reads for de novo clustering = %u\n",total_reads_denovo_clustering);
    }
    // output total non-rrna + rrna reads
    fprintf(bilan,"    Total reads passing E-value threshold = %u (%.2f%%)\n",total_reads_mapped,(float)((float)total_reads_mapped/(float)number_total_read)*100);
    fprintf(bilan,"    Total reads failing E-value threshold = %u (%.2f%%)\n",number_total_read-total_reads_mapped,(1-((float)((float)total_reads_mapped/(float)number_total_read)))*100);
    fprintf(bilan,"    Minimum read length = %u\n", min_read_len);
    fprintf(bilan,"    Maximum read length = %u\n", max_read_len);
    fprintf(bilan,"    Mean read length = %u\n", mean_read_len);
    
    fprintf(bilan," By database:\n");    
    // output stats by database
    for ( uint32_t index_num = 0; index_num < myfiles.size(); index_num++ )
    {
      fprintf(bilan,"    %s\t\t%.2f%%\n",(char*)(myfiles[index_num].first).c_str(),(float)((float)reads_matched_per_db[index_num]/(float)number_total_read)*100);
    }
    
    if ( otumapout_gv )
    {
      fprintf(bilan," Total reads passing %%id and %%coverage thresholds = %u\n", total_reads_mapped_cov);
      fprintf(bilan," Total OTUs = %lu\n", otu_map.size());
    }

    time_t q = time(0);
    struct tm * now = localtime(&q);
    
    fprintf(bilan,"\n %s\n",asctime(now));
    
    fclose(bilan);
    
    // free memory of accepted strings
    if ( logoutfile != NULL )
    {
      delete [] logoutfile;
      logoutfile = NULL;
    }
  }
  else if ( otumapout_gv ) otu_map.clear();
  
  // free memory of accepted strings
  if ( acceptedstrings != NULL )
  {
    delete [] acceptedstrings;
    acceptedstrings = NULL;
  }
  
  if ( acceptedstrings_sam != NULL )  
  {
    delete [] acceptedstrings_sam;
    acceptedstrings_sam = NULL;
  }
  
  if ( acceptedstrings_blast != NULL )
  {
    delete [] acceptedstrings_blast;
    acceptedstrings_blast = NULL;
  }

  if ( denovo_otus_file != NULL )
  {
    delete [] denovo_otus_file;
    denovo_otus_file = NULL;
  }
    
  return ;  
}//~paralleltraversal()
