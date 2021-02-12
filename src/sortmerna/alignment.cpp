/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SortMeRNA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 @endparblock

 @contributors Jenya Kopylova   jenya.kopylov@gmail.com
			   Laurent Noé      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mikaël Salson    mikael.salson@lifl.fr
			   Hélène Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

/*
 * @file alignment.cpp
 * @brief functions for sequence alignment.
 */

#include <sstream>
#include <fstream>

#include "alignment.hpp"
#include "read.hpp"
#include "index.hpp"
#include "refstats.hpp"
#include "references.hpp"
#include "readstats.hpp"

#define ASCENDING <
#define DESCENDING >


// forward
s_align2 copyAlignment(s_align* pAlign);
uint32_t inline findMinIndex(std::vector<s_align2>& alignv);
uint32_t inline findMaxIndex(std::vector<s_align2>& alignv);
std::pair<bool,bool> is_id_cov_pass(std::string& read_iseq, s_align2& alignment, References& refs, Runopts& opts);

void find_lis( deque<pair<uint32_t, uint32_t>>& a, vector<uint32_t>& b )
{
	vector<uint32_t> p(a.size());
	std::size_t u, v;

	if (a.empty()) return;

	b.push_back(0);

	for (std::size_t i = 1; i < a.size(); i++)
	{
		// If next element a[i] is greater than last element of current longest subsequence a[b.back()], just push it at back of "b" and continue
		if (a[b.back()].second < a[i].second)
		{
			p[i] = b.back();
			b.push_back(static_cast<uint32_t>(i));
			continue;
		}

		// Binary search to find the smallest element referenced by b which is just bigger than a[i]
		// Note : Binary search is performed on b (and not a). Size of b is always <=k and hence contributes O(log k) to complexity.
		for (u = 0, v = b.size() - 1; u < v;)
		{
			auto c = (u + v) / 2;
			if (a[b[c]].second < a[i].second)
				u = c + 1;
			else
				v = c;
		}

		// Update b if new value is smaller then previously referenced value
		if (a[i].second < a[b[u]].second)
		{
			if (u > 0) p[i] = b[u - 1];
			b[u] = static_cast<uint32_t>(i);
		}
	}

	for (u = b.size(), v = b.back(); u--; v = p[v]) 
		b[u] = static_cast<uint32_t>(v);
} // ~find_lis

void compute_lis_alignment( Read& read, Runopts& opts,
							Index& index, References& refs, 
							Readstats& readstats, Refstats& refstats,
							bool& search, uint32_t max_SW_score	)
{
	// true if SW alignment between the read and a candidate reference meets the threshold
	bool is_aligned = false;

	map<uint32_t, uint32_t> kmer_count_map;
	//    |         |_number of k-mer hits on the reference
	//    |_reference number/position in the ref file

	vector<uint32pair> kmer_count_vec; // use to sort the 'kmer_count_map' (map cannot be sorted)
	map<uint32_t, uint32_t>::iterator map_it;
	uint32_t max_ref = 0; // reference with max kmer occurrences
	uint32_t max_occur = 0; // number of kmer occurrences on the 'max_ref'

	// 1. Find all candidate references by using Read's kmer hits information.
	//    For every reference, compute the number of kmer hits belonging to it
	for (auto const& hit: read.id_win_hits)
	{
		seq_pos* positions_tbl_ptr = index.positions_tbl[hit.id].arr;
		// loop all positions of id
		for (uint32_t j = 0; j < index.positions_tbl[hit.id].size; j++)
		{
			uint32_t seq = positions_tbl_ptr++->seq;
			if ((map_it = kmer_count_map.find(seq)) != kmer_count_map.end())
				map_it->second++; // sequence already in the map, increment its frequency value
			else
				kmer_count_map[seq] = 1; // sequence not in the map, add it
		}
	}

	// copy frequency map to vector for sorting
	// consider only candidate references that have enough seed hits
	for (auto const& freq_pair: kmer_count_map)
	{
		if (freq_pair.second >= (uint32_t)opts.num_seeds)
			kmer_count_vec.push_back(freq_pair);
	}

	kmer_count_map.clear();

	// sort sequences by frequency in descending order
	auto cmp = [](std::pair<uint32_t, uint32_t> e1, std::pair<uint32_t, uint32_t> e2) {
		if (e1.second == e2.second)
			return e1.first ASCENDING e2.first; // order references ascending for equal frequencies (originally - descending)
		return e1.second DESCENDING e2.second; // order frequencies descending
	}; // comparator
	std::sort(kmer_count_vec.begin(), kmer_count_vec.end(), cmp);

	// 2. loop reference candidates, starting from the one with the highest number of kmer hits.
	bool is_search_candidates = true;
	for (uint32_t k = 0; k < kmer_count_vec.size() && is_search_candidates; k++)
	{
		max_ref = kmer_count_vec[k].first;
		max_occur = kmer_count_vec[k].second;
              
		// not enough hits on the reference, try to collect more hits or next read
		if (max_occur < (uint32_t)opts.num_seeds) {
			break;
		}

		// update number of reference sequences remaining to check
		// only decrement read.best if the current ref sequence to check
		// has a lower seed count than the previous one
		if (is_aligned && opts.min_lis > 0 && k > 0 && max_occur < kmer_count_vec[k - 1].second )
		{
			--read.best;
			if (read.best < 1) break;
		}

		// list of matching kmer pairs on a given reference: 
		//  [pair<1st:k-mer ref pos, 2nd:k-mer read pos>] e.g.
		//  [ (493, 0), ..., (674, 18), ... ]
		//      |   |_k-mer position on the read
		//      |_k-mer position on the reference
		vector<uint32pair> hits_on_ref;

		//
		// 3. populate 'hits_on_ref'
		//
		for ( auto const& hit: read.id_win_hits )
		{
			uint32_t num_hits = index.positions_tbl[hit.id].size;
			seq_pos* positions_tbl_ptr = index.positions_tbl[hit.id].arr;
			// loop through every position of id
			for (uint32_t j = 0; j < num_hits; j++)
			{
				if (positions_tbl_ptr->seq == max_ref)
				{
					hits_on_ref.push_back(uint32pair(positions_tbl_ptr->pos, hit.win));
				}
				positions_tbl_ptr++;
			}
		}

		// sort the positions in ascending order
		std::sort(hits_on_ref.begin(), hits_on_ref.end(), [](uint32pair e1, uint32pair e2) {
			if (e1.first == e2.first) 
				return (e1.second ASCENDING e2.second); // order references ascending for equal reference positions
			return (e1.first ASCENDING e2.first);
		}); // smallest

		// iterate over the set of hits, searching for windows of
		// win.len == read.len which have at least ratio hits
		vector<uint32pair>::iterator hits_on_ref_iter = hits_on_ref.begin();
		deque<uint32pair> match_set; // set of matching k-mers fit within the read length: [pair<1st:on ref pos, 2nd:on read pos>]

		// 4. run a sliding window of read's length along the reference, 
		//    searching for windows with enough k-mer hits
		uint32_t lcs_ref_start = 0; // match (LCS) start position on reference
		uint32_t lcs_que_start = 0; // match (LCS) start position on read
		uint32_t begin_ref = hits_on_ref_iter->first; // hit position on reference
		uint32_t begin_read = hits_on_ref_iter->second; // hit position on read
                       
		// TODO: Always does a single iteration because of the line '++hits_on_ref_iter'. 
		//       It has 3 'break' instructions though. Convoluted.
		while (hits_on_ref_iter != hits_on_ref.end() && is_search_candidates)
		{
			// max possible k-mer start position on reference: 
			//   max start position on the reference of a matching k-mer for 
			//   an overlaid read anchored on the reference using a matching k-mer
			// 
			// ref: |--------|k-mer anchor|-------|k-mer|--------------------|k-mer|-----|
			//               ^begin_ref e.g. 20                              ^end_ref_max
			// read:    |----|k-mer anchor|----|k-mer|----------|k-mer|----|---|
			//          |    ^begin_read e.g 10                 ^end_read  |___|_ lnwin
			//          ^read start pos                                    ^   ^read end pos
			//                                                             |_max possible k-mer start position 'end_ref_max'
			// 
			auto end_ref_max = begin_ref - begin_read + read.sequence.length() - refstats.lnwin[index.index_num];
			//auto end_ref_max = begin_ref + read.sequence.length() - refstats.lnwin[index.index_num] + 1; // TODO: original - wrong?
			bool push = false;
			while ( hits_on_ref_iter != hits_on_ref.end() && hits_on_ref_iter->first <= end_ref_max )
			{
				match_set.push_back(*hits_on_ref_iter);
				push = true;
				++hits_on_ref_iter;
			}
			// heuristic 1: a new window hit was not pushed back, pop queue until new window can be pushed back
			// this heuristic significantly speeds up the algorithm because we don't perform alignments for
			// every sub-LIS of a window if an alignment reaching threshold has already been made. It assumes
			// that every sub-LIS yields the same alignment score, which is true for 99.99% of cases.
#ifndef HEURISTIC1_OFF
			if (!push && is_aligned) goto pop;
			else is_aligned = false;
#endif
#ifdef HEURISTIC1_OFF
			aligned = false;
#endif                              
			// enough windows at this position on genome to search for LIS
			if (match_set.size() >= (uint32_t)opts.num_seeds)
			{
				vector<uint32_t> lis_arr; // array of Indices of matches from the match_set comprising the LIS
				find_lis(match_set, lis_arr);
#ifdef HEURISTIC1_OFF
				uint32_t list_n = 0;
				do
				{
#endif                                      
					// LIS long enough to perform Smith-Waterman alignment
					if (lis_arr.size() >= (size_t)opts.min_lis)
					{
#ifdef HEURISTIC1_OFF
						lcs_ref_start = match_set[lis_arr[list_n]].first;
						lcs_que_start = match_set[lis_arr[list_n]].second;
#endif
#ifndef HEURISTIC1_OFF
						lcs_ref_start = match_set[lis_arr[0]].first;
						lcs_que_start = match_set[lis_arr[0]].second;
#endif                                    
						// reference string
						std::size_t head = 0;
						std::size_t tail = 0;
						std::size_t align_ref_start = 0;
						std::size_t align_que_start = 0;
						std::size_t align_length = 0;
						auto reflen = refs.buffer[max_ref].sequence.length();
						uint32_t edges = 0;
						if (opts.is_as_percent)
							edges = (((double)opts.edges / 100.0)*read.sequence.length());
						else
							edges = opts.edges;
						// part of the read hangs off (or matches exactly) the beginning of the reference seq
						//            ref |-----------------------------------|
						// que |-------------------|
						//             LIS |-----|
						//
						if (lcs_ref_start < lcs_que_start)
						{
							align_ref_start = 0;
							align_que_start = lcs_que_start - lcs_ref_start;
							head = 0;
							// the read is longer than the reference sequence
							//            ref |----------------|
							// que |---------------------...|
							//                LIS |-----|
							//
							if (reflen < read.sequence.length())
							{
								tail = 0;
								// beginning from align_ref_start = 0 and align_que_start = X, the read finishes
								// before the end of the reference
								//            ref |----------------|
								// que |------------------------|
								//                  LIS |-----|
								//                ^
								//                align_que_start
								if (align_que_start >(read.sequence.length() - reflen))
								{
									align_length = reflen - (align_que_start - (read.sequence.length() - reflen));
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
								tail = reflen - align_ref_start - read.sequence.length();
								tail > (edges - 1) ? tail = edges : tail;
								align_length = read.sequence.length() + head + tail - align_que_start;
							}
						}
						else
						{
							align_ref_start = lcs_ref_start - lcs_que_start;
							align_que_start = 0;
							align_ref_start > (edges - 1) ? head = edges : head;
							// part of the read hangs off the end of the reference seq
							// ref |-----------------------------------|
							//                          que |-------------------|
							//                            LIS |-----|
							//
							if (align_ref_start + read.sequence.length() > reflen) // readlen
							{
								tail = 0;
								align_length = reflen - align_ref_start - head;
							}
							// the reference seq fully covers the read
							// ref |-----------------------------------|
							//    que |-------------------|
							//          LIS |-----|
							//
							else
							{
								tail = reflen - align_ref_start - read.sequence.length();
								tail > (edges - 1) ? tail = edges : tail;
								align_length = read.sequence.length() + head + tail;
							}
						}

						// put read into 04 encoding before SSW
						if (read.is03) 
							read.flip34();
                       
						// create profile for read
						s_profile* profile = 0;
						profile = ssw_init((int8_t*)(&read.isequence[0] + align_que_start), (align_length - head - tail), &read.scoring_matrix[0], 5, 2);

						s_align* result = 0;

						result = ssw_align(
							profile,
							(int8_t*)refs.buffer[max_ref].sequence.c_str() + align_ref_start - head,
							align_length,
							opts.gap_open,
							opts.gap_extension,
							2,
							refstats.minimal_score[index.index_num], // minimal_score_index_num
							0,
							0
						);

						// deallocate memory for profile, no longer needed
						if (profile != 0) 
							init_destroy(&profile);

						// check alignment passes the threshold
						is_aligned = (result != 0 && result->score1 > refstats.minimal_score[index.index_num]);
						if (is_aligned)
						{
							++read.num_hits;
							if (result->score1 == max_SW_score) 
								++read.max_SW_count; // a max possible score has been found

							// add the offset calculated by the LCS (from the beginning of the sequence)
							// to the offset computed by SW alignment
							result->ref_begin1 += (align_ref_start - head);
							result->ref_end1 += (align_ref_start - head);
							result->read_begin1 += align_que_start;
							result->read_end1 += align_que_start;
							result->readlen = read.sequence.length();
							result->ref_num = max_ref;

							result->index_num = index.index_num;
							result->part = index.part;
							result->strand = !read.reversed; // flag whether the alignment was done on a forward or reverse strand

							s_align2 alignment = copyAlignment(result); // new alignment

							// read has not yet been mapped, set bit to true for this read
							// (this is the Only place where read.is_hit can be modified)
							if (!read.is_hit)
							{
								read.is_hit = true;
								readstats.total_aligned.fetch_add(1, std::memory_order_relaxed);
								++readstats.reads_matched_per_db[index.index_num];
							}

							// calculate %ID and %COV if OTU map or denovo were requested
							if ((opts.is_otu_map || opts.is_denovo) && !(read.is_id && read.is_cov))
							{
								std::pair<bool,bool> is_id_cov = is_id_cov_pass(read.isequence, alignment, refs, opts);
								if (!read.is_id && is_id_cov.first) {
									read.is_id = true;
								}
								if (!read.is_cov && is_id_cov.second) {
									read.is_cov = true;
								}
							}

							// if 'N == 0' or 'Not is_best' or 'is_best And read.alignments.size < N' => 
							//   simply add the new alignment to read.alignments
							if (opts.num_alignments == 0 || !opts.is_best || (opts.is_best && read.alignment.alignv.size() < opts.num_alignments))
							{
								read.alignment.alignv.emplace_back(alignment);
								read.is_new_hit = true; // flag to store in DB
							}
							else if ( opts.is_best 
									&& read.alignment.alignv.size() == opts.num_alignments 
									&& read.alignment.alignv[read.alignment.min_index].score1 < result->score1 )
							{
								if (opts.is_best_id_cov) {
									// TODO: new case to implement 20200703
								}
								else {
									// set min and max pointers - just once, after all the reads' alignments were filled
									if (opts.num_alignments > 1 && read.alignment.max_index == 0 && read.alignment.min_index == 0) {
										read.alignment.min_index = findMinIndex(read.alignment.alignv);
										read.alignment.max_index = findMaxIndex(read.alignment.alignv);
									}

									uint32_t min_score_index = read.alignment.min_index;
									uint32_t max_score_index = read.alignment.max_index;

									// replace the old smallest scored alignment with the new one
									read.alignment.alignv[min_score_index] = alignment;
									read.is_new_hit = true; // flag to store in DB

									// if new_hit > max_hit: the old min_hit_idx becomes the new max_hit_idx
									// only do if num_alignments > 1 i.e. max_idx != min_idx
									if (result->score1 > read.alignment.alignv[max_score_index].score1 && read.alignment.alignv.size() > 1) {
										read.alignment.max_index = min_score_index; // new max index
										read.alignment.min_index = findMinIndex(read.alignment.alignv); // new min index
									}

									// decrement number of reads mapped to database with lower score
									--readstats.reads_matched_per_db[read.alignment.alignv[min_score_index].index_num];
									//                                                           |_old min index
									// increment number of reads mapped to database with higher score
									++readstats.reads_matched_per_db[index.index_num];
								}
							}//~if

							// all alignments have been found - stop searching
							if (opts.num_alignments > 0) {
								if (opts.is_best) {
									if (opts.num_alignments == read.max_SW_count)
										is_search_candidates = false;
								}
								else if (opts.num_alignments == read.alignment.alignv.size())
									is_search_candidates = false;
							}

							// continue to next read (do not need to collect more seeds using another pass)
							search = false;
						}//~if aligned
						
						// free alignment info
						if(result != 0)
						{
							free(result);
							result = 0;
						}
					}//~if LIS long enough                               
#ifdef HEURISTIC1_OFF
				} while ((hits_on_ref_iter == hits_per_ref.end()) && (++list_n < lis_arr.size()));
#endif                                              
			}//~if enough window hits                                                
		pop:
			// get the next candidate reference position 
			if (!match_set.empty())
			{
				match_set.pop_front();
			}

			if (match_set.empty())
			{
				if (hits_on_ref_iter != hits_on_ref.end()) // TODO: seems Always false
				{
					begin_ref = hits_on_ref_iter->first; // TODO: seems never reached
					begin_read = hits_on_ref_iter->second;
				}
				else break;
			}
			else
			{
				begin_ref = (match_set.front()).first;
				begin_read = (match_set.front()).second;
			}
		}//~for all reference sequence length                   
	}//~for all of the reference sequence candidates
} // ~compute_lis_alignment

s_align2 copyAlignment(s_align* pAlign)
{
	s_align2 ret_align;
	for (int i = 0; i < pAlign->cigarLen; i++)
		ret_align.cigar.push_back(*pAlign->cigar++);
	ret_align.index_num = pAlign->index_num;
	ret_align.part = pAlign->part;
	ret_align.readlen = pAlign->readlen;
	ret_align.read_begin1 = pAlign->read_begin1;
	ret_align.read_end1 = pAlign->read_end1;
	ret_align.ref_begin1 = pAlign->ref_begin1;
	ret_align.ref_end1 = pAlign->ref_end1;
	ret_align.ref_num = pAlign->ref_num;
	ret_align.score1 = pAlign->score1;
	ret_align.strand = pAlign->strand;

	return ret_align;
}

/* 
 * find the index of the alignment with the smallest score 
 */
uint32_t inline findMinIndex(std::vector<s_align2>& alignv)
{
	uint32_t min_score = alignv[0].score1;
	uint32_t min_idx = 0;
	for (int i = 0; i < alignv.size(); ++i)
	{
		if (alignv[i].score1 < min_score)
		{
			min_score = alignv[i].score1;
			min_idx = i;
		}
	}
	return min_idx;
}

uint32_t inline findMaxIndex(std::vector<s_align2>& alignv)
{
	uint32_t max_idx = 0;
	uint32_t max_score = alignv[0].score1;
	for (int i = 0; i < alignv.size(); ++i)
	{
		if (alignv[i].score1 > max_score)
		{
			max_score = alignv[i].score1;
			max_idx = i;
		}
	}
	return max_idx;
}

/* 
 * calculate whether the alignment passes ID and COV thresholds 
 *
 * @param read_iseq  read sequence in integer alphabet, see 'read.isequence'
 * @param alignment to check
 * @return pair (is_ID, is_COV)  is_ID: true | false, is_COV: true | false
 *
 */
std::pair<bool,bool> inline is_id_cov_pass(std::string& read_iseq, s_align2& alignment, References& refs, Runopts& opts)
{
	// calculate id, mismatches, gaps for the given alignment
	int id = 0; // count of mismatched characters
	int mismatches = 0; // count of gaps
	int gaps = 0; // count of matched characters

	int32_t ridx = alignment.ref_begin1; // index of the first char in the reference matched part
	int32_t qidx = alignment.read_begin1; // index of the first char in the read matched part

	std::string refseq = refs.buffer[alignment.ref_num].sequence; // reference sequence
	int32_t align_len = abs(alignment.read_end1 + 1 - alignment.read_begin1); // alignment length

	for (uint32_t cidx = 0; cidx < alignment.cigar.size(); ++cidx)
	{
		uint32_t letter = 0xf & alignment.cigar[cidx]; // 4 low bits
		uint32_t length = (0xfffffff0 & alignment.cigar[cidx]) >> 4; // high 28 bits i.e. 32-4=28
		if (letter == 0)
		{
			for (uint32_t u = 0; u < length; ++u)
			{
				if (refseq[ridx] != read_iseq[qidx]) ++mismatches;
				else ++id;
				++ridx;
				++qidx;
			}
		}
		else if (letter == 1)
		{
			qidx += length;
			gaps += length;
		}
		else
		{
			ridx += length;
			gaps += length;
		}
	}

	// round to 3 decimal places
	stringstream ss;
	ss.precision(3);
	ss << (double)id / ((double)mismatches + gaps + id) << ' ' << (double)align_len / read_iseq.length();

	double align_id_round = 0.0;
	double align_cov_round = 0.0;
	ss >> align_id_round >> align_cov_round;

	return { align_id_round >= opts.min_id, align_cov_round >= opts.min_cov };
} // ~is_id_cov_pass
