/**
 * @file alignment.cpp
 * @brief File containing functions for sequence alignment.
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

#include "../include/alignment.hpp"

bool
smallest ( const mypair &a, const mypair &b )
{
  if ( a.first == b.first ) return ( a.second < b.second );
  return ( a.first < b.first );
}

bool
largest ( const mypair &a, const mypair &b )
{
  if ( a.first == b.first ) return ( a.second > b.second );
  return ( a.first > b.first );
}

/*
 *
 * FUNCTION   : find_lis()
 * see alignment.hpp for documentation
 *************************************/
void find_lis(deque<pair<uint32_t, uint32_t> > &a,
              vector<uint32_t> &b, uint64_t readn)
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

/*
 *
 * FUNCTION   : compute_lis_alignment()
 * see alignment.hpp for documentation
 **************************************/
void
compute_lis_alignment(uint32_t size_ambiguous_nt,
                      uint32_t readhit,
                      vector< id_win >& id_win_hits,
                      kmer_origin* positions_tbl,
                      uint16_t* read_max_SW_score,
                      bool& search,
                      int32_t* best_x,
                      uint64_t readn,
                      int32_t* num_alignments_x,
                      uint32_t readlen,
                      uint32_t lnwin_index_num,
                      uint16_t index_num,
                      uint64_t* reference_seq_len,
                      char* myread,
                      int32_t* ambiguous_nt,
                      int8_t* mat,
                      char** reference_seq,
                      long gap_open,
                      long gap_extension,
                      uint32_t minimal_score_index_num,
                      vector<bool>& read_hits,
                      uint64_t& total_reads_mapped,
                      vector<uint64_t>& reads_matched_per_db,
                      uint16_t part,
                      map<uint64_t, alignment_struct >& read_hits_align_info,
                      uint32_t max_SW_score,
                      bool& read_to_count,
                      uint64_t& total_reads_mapped_cov,
                      vector<bool>& read_hits_denovo,
                      char filesig,
                      uint64_t strs,
                      uint32_t file_s,
                      uint32_t file_sections,
                      char** reads,
                      char* finalnt,
                      double gumbel_lambda_index_num,
                      double gumbel_K_index_num,
                      uint64_t full_ref_index_num,
                      uint64_t full_read_index_num,
                      ofstream& acceptedblast,
                      ofstream& acceptedsam)
{
// flag to check whether read was formatted to 0-4 alphabet
// for alignment allowing N's
bool read_edited = false;
// boolean set to true if SW alignment succeeded between
// the read and a candidate reference sequence
bool aligned = false;
#ifdef debug_align
  cout << "\t\t\treadhit = " << readhit << endl;
  cout << "\t\t\tseed_hits_gv = " << seed_hits_gv << endl;
  cout << "\t\t\tnum_best_hits_gv = " << num_best_hits_gv << endl;
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
        uint32_t stop = begin + readlen - lnwin_index_num + 1;
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
              for ( uint32_t p = 0; p < list.size(); p++ )
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
                                 (int8_t*)reference_seq[(2*(int)max_seq)+1]+align_ref_start-head,
                                 align_length,
                                 gap_open,
                                 gap_extension,
                                 2,
                                 minimal_score_index_num,
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
                if ( result->score1 > minimal_score_index_num ) aligned = true;
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
                    map<uint64_t, alignment_struct>::iterator alignment = read_hits_align_info.find(readn);
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
                      if ( (num_best_hits_gv == 0) || (array_size < (uint32_t)num_best_hits_gv) )
                      {
                        // number of alignments stored per read == maximum number of 
                        // alignments allowed, resize array by another BEST_HITS_INCREMENT slots 
                        if ( array_size == array_max_size )
                        {
                          uint32_t new_array_max_size = 0;
                          if ( (num_best_hits_gv == 0) || (array_size + BEST_HITS_INCREMENT <= (uint32_t)num_best_hits_gv) )
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
                                      "in paralleltraversal.cpp\n", startColor, "\033[0m");
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
                        if ( array_size == (uint32_t)num_best_hits_gv )
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
                                       "storage (paralleltraversal.cpp)\n",startColor,"\033[0m");
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
                        for (uint32_t p = 0; p < length; ++p)
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
                        uint32_t bitscore = (uint32_t)((float)((gumbel_lambda_index_num)*(result->score1) - log((gumbel_K_index_num)))/(float)log(2));
                        double evalue_score = (double)(gumbel_K_index_num)*full_ref_index_num*full_read_index_num*pow(EXP,(-(gumbel_lambda_index_num)*result->score1));
                        
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
  if ( (size_ambiguous_nt != 0) && read_edited )
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

  return ;
 }//~compute_alignment()
