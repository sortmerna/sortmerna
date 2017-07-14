/**
 * @file traverse_bursttrie.cpp
 * @brief Traverse the burst trie searching for candidate sequences.
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

#include "../include/traverse_bursttrie.hpp"

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

 /*! @fn traversetrie_align() */
void
traversetrie_align ( NodeElement *trie_t,
                     uint32_t lev_t,
                     unsigned char depth,
                     MYBITSET *win_k1_ptr,
                     MYBITSET *win_k1_full,
                     bool &accept_zero_kmer,
                     vector< id_win > &id_hits,
                     int64_t readn,
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
          traversetrie_align(trie_t->whichnode.trie,
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
                  for ( uint32_t f = 0; f < id_hits.size(); f++ )
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