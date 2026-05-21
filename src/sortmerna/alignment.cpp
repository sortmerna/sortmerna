/*
@copyright 2016-2026 Clarity Genomics Inc
@copyright 2012-2016 Bonsai Bioinformatics Research Group
@copyright 2014-2016 Knight Lab, Department of Pediatrics, UCSD, La Jolla

@parblock
SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA

This is a free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SortMeRNA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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

#include "parasail.h"
#include "alignment.hpp"
#include "read.hpp"
#include "index.hpp"
#include "refstats.hpp"
#include "references.hpp"
#include "readstats.hpp"

#define ASCENDING <
#define DESCENDING >


// forward
uint32_t inline findMinIndex(std::vector<s_align2>& alignv);
uint32_t inline findMaxIndex(std::vector<s_align2>& alignv);
std::pair<bool,bool> is_id_cov_pass(std::string& read_iseq, s_align2& alignment, References& refs, Runopts& opts);

/* Build a parasail scoring matrix for the ACGTN alphabet with SortMeRNA's scoring scheme.
   score_N is applied to all pairs involving N (index 4), which may differ from mismatch. */
static parasail_matrix_t* smr_create_matrix(int match, int mismatch, int score_N)
{
    parasail_matrix_t* mat = parasail_matrix_create("ACGTN", match, mismatch);
    for (int j = 0; j < 5; ++j)
        parasail_matrix_set_value(mat, 4, j, score_N); // N row
    for (int i = 0; i < 4; ++i)
        parasail_matrix_set_value(mat, i, 4, score_N); // N column
    return mat;
}

/* Decode an integer-encoded (0-4 = A,C,G,T,N) sequence slice to character form for parasail. */
static std::string smr_decode_seq(const char* src, int len)
{
    static const char nt_decode[5] = {'A', 'C', 'G', 'T', 'N'};
    std::string s(len, 'N');
    for (int i = 0; i < len; ++i) {
        const auto b = static_cast<uint8_t>(src[i]);
        s[i] = nt_decode[b < 5u ? b : 4u];
    }
    return s;
}

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

	parasail_matrix_t* matrix = smr_create_matrix(opts.match, opts.mismatch, opts.score_N);

	map<uint32_t, uint32_t> refs_kmer_count_map; // map of kmer hits on candidate references
	//    |         |_number of k-mer hits on the reference
	//    |_reference number/position in the ref file

	vector<uint32pair> refs_kmer_count_vec; // use to sort the 'refs_kmer_count_map' (map cannot be sorted)
	map<uint32_t, uint32_t>::iterator map_it;
	uint32_t max_ref = 0; // reference with max kmer occurrences
	uint32_t max_occur = 0; // number of kmer occurrences on the 'max_ref'

	// 1. For each candidate reference compute the number of kmer hits belonging to it
	for (auto const& hit: read.id_win_hits)
	{
		seq_pos* positions_tbl_ptr = index.positions_tbl[hit.id].arr;
		// loop all positions of id
		for (uint32_t j = 0; j < index.positions_tbl[hit.id].size; j++)
		{
			uint32_t seq = positions_tbl_ptr++->seq;
			if ((map_it = refs_kmer_count_map.find(seq)) != refs_kmer_count_map.end())
				map_it->second++; // sequence already in the map, increment its frequency value
			else
				refs_kmer_count_map[seq] = 1; // sequence not in the map, add it
		}
	}

	// copy frequency map to vector for sorting
	// consider only candidate references that have enough seed hits
	for (auto const& freq_pair: refs_kmer_count_map)
	{
		if (freq_pair.second >= (uint32_t)opts.num_seeds)
			refs_kmer_count_vec.push_back(freq_pair);
	}

	refs_kmer_count_map.clear();

	// sort sequences by frequency in descending order
	auto cmp = [](std::pair<uint32_t, uint32_t> e1, std::pair<uint32_t, uint32_t> e2) {
		if (e1.second == e2.second)
			return e1.first ASCENDING e2.first; // order references ascending for equal frequencies (originally - descending)
		return e1.second DESCENDING e2.second; // order frequencies descending
	}; // comparator
	std::sort(refs_kmer_count_vec.begin(), refs_kmer_count_vec.end(), cmp);

	// 2. loop reference candidates, starting from the one with the highest number of kmer hits.
	auto is_search_candidates = true;
	for (uint32_t k = 0; k < refs_kmer_count_vec.size() && is_search_candidates; k++)
	{
		max_ref = refs_kmer_count_vec[k].first;
		max_occur = refs_kmer_count_vec[k].second;
              
		// not enough hits on the reference, try to collect more hits or next read
		if (max_occur < (uint32_t)opts.num_seeds) {
			break; // 'search = true' here
		}

		// update number of reference sequences remaining to check
		// only decrement read.best if the current ref sequence to check
		// has a lower seed count than the previous one
		if (is_aligned && opts.min_lis > 0 && k > 0 && max_occur < refs_kmer_count_vec[k - 1].second )
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
			//auto end_ref_max = begin_ref + read.sequence.length() - refstats.lnwin[index.index_num] + 1; // TODO: original - remove
			auto end_ref_max = begin_ref + read.sequence.length() - begin_read - refstats.lnwin[index.index_num] + 1;
			auto push = false;
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
							edges = static_cast<decltype(edges)>((opts.edges / 100.0) * read.sequence.length());
						else
							edges = static_cast<decltype(edges)>(opts.edges);
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

						// put read into 04 encoding before alignment
						if (read.is03)
							read.flip34();

						// decode integer-encoded slices to character form for parasail
						const int qlen = static_cast<int>(align_length - head - tail);
						const int rlen = static_cast<int>(align_length);
						const std::string qchars = smr_decode_seq(&read.isequence[0] + align_que_start, qlen);
						const std::string rchars = smr_decode_seq(
							refs.buffer[max_ref].sequence.c_str() + align_ref_start - head, rlen);

						parasail_profile_t* pprofile = parasail_profile_create_sat(
							qchars.c_str(), qlen, matrix);

						parasail_result_t* presult = parasail_sw_trace_striped_profile_sat(
							pprofile, rchars.c_str(), rlen, opts.gap_open, opts.gap_extension);

						parasail_profile_free(pprofile);

						const int pscore = parasail_result_get_score(presult);
						is_aligned = (presult != nullptr && pscore > refstats.minimal_score[index.index_num]);
						if (is_aligned)
						{
							if (static_cast<uint32_t>(pscore) == max_SW_score)
								++read.max_SW_count; // a max possible score has been found

							const int32_t ref_offset = static_cast<int32_t>(align_ref_start - head);
							const int32_t que_offset = static_cast<int32_t>(align_que_start);

							s_align2 alignment;
							parasail_cigar_t* pcigar = parasail_result_get_cigar(
								presult, qchars.c_str(), qlen, rchars.c_str(), rlen, matrix);
							if (pcigar) {
								// Normalise parasail CIGAR to BAM encoding expected by callers:
								// parasail uses '='(op=7) for exact match and 'X'(op=8) for mismatch;
								// all callers expect combined M(op=0). Merge adjacent runs.
								// I(op=1) and D(op=2) are already in standard orientation.
								for (int ci = 0; ci < pcigar->len; ++ci) {
									uint32_t op  = pcigar->seq[ci] & 0xfu;
									uint32_t len = pcigar->seq[ci] >> 4u;
									if (op == 7 || op == 8) op = 0; // '='/'X' -> M
									// merge consecutive M blocks from adjacent '='+'X' runs
									if (op == 0 && !alignment.cigar.empty()
											&& (alignment.cigar.back() & 0xfu) == 0)
										alignment.cigar.back() += (len << 4u);
									else
										alignment.cigar.push_back((len << 4u) | op);
								}
								alignment.ref_begin1  = pcigar->beg_ref + ref_offset;
								alignment.read_begin1 = pcigar->beg_query + que_offset;
								parasail_cigar_free(pcigar);
								// Parasail adds leading D ops when traceback exits through the query
								// boundary (i<0): remaining reference-window positions are dumped as D.
								// SW local alignment never genuinely starts with a deletion (negative
								// score). Strip them and advance ref_begin1 past the head-extension bases.
								while (!alignment.cigar.empty() && (alignment.cigar[0] & 0xfu) == 2) {
									alignment.ref_begin1 += alignment.cigar[0] >> 4u;
									alignment.cigar.erase(alignment.cigar.begin());
								}
							}
							alignment.ref_end1  = parasail_result_get_end_ref(presult)   + ref_offset;
							alignment.read_end1 = parasail_result_get_end_query(presult) + que_offset;
							alignment.score1    = static_cast<uint16_t>(pscore);
							alignment.readlen   = static_cast<uint32_t>(read.sequence.length());
							alignment.ref_num   = max_ref;
							alignment.index_num = index.index_num;
							alignment.part      = index.part;
							alignment.strand    = !read.reversed;

							// read has not yet been mapped, set bit to true for this read
							// (this is the Only place where read.is_hit can be modified)
							if (!read.is_hit)
							{
								read.is_hit = true;
								readstats.num_aligned.fetch_add(1, std::memory_order_relaxed);
								++readstats.reads_matched_per_db[index.index_num];
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
									&& read.alignment.alignv[read.alignment.min_index].score1 < pscore )
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
									if (pscore > read.alignment.alignv[max_score_index].score1 && read.alignment.alignv.size() > 1) {
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

							// if all alignments have been found - stop searching
							if (opts.num_alignments > 0) {
								if (opts.is_best) {
									if (opts.num_alignments == read.max_SW_count)
										is_search_candidates = false;
								}
								else if (opts.num_alignments == read.alignment.alignv.size())
									is_search_candidates = false;
							}

							// continue to next read (no need to collect more seeds using another pass)
							search = false;
						}//~if aligned

						parasail_result_free(presult);
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
		}//~for all matching k-mers on a reference
	}//~for all reference candidates

	parasail_matrix_free(matrix);
} // ~compute_lis_alignment

/* 
 * find the index of the alignment with the smallest score 
 */
uint32_t inline findMinIndex(std::vector<s_align2>& alignv)
{
	uint32_t min_score = alignv[0].score1;
	uint32_t min_idx = 0;
	for (unsigned i = 0; i < alignv.size(); ++i)
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
	for (unsigned i = 0; i < alignv.size(); ++i)
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
 * @return pair<is_ID, is_COV>  is_ID: true | false, is_COV: true | false
 *
 */
std::pair<bool,bool> inline is_id_cov_pass(std::string& read_iseq, s_align2& alignment, References& refs, Runopts& opts)
{
	// calculate id, mismatches, gaps for the given alignment
	int id = 0; // count of matched nt
	int mismatches = 0; // count of mismatched nt
	int gaps = 0; // count of gaps

	auto read_i = alignment.ref_begin1; // index of the first char in the reference matched part
	auto query_i = alignment.read_begin1; // index of the first char in the read matched part

	std::string refseq = refs.buffer[alignment.ref_num].sequence; // reference sequence
	int32_t align_len = abs(alignment.read_end1 + 1 - alignment.read_begin1); // alignment length

	for (uint32_t cigar_i = 0; cigar_i < alignment.cigar.size(); ++cigar_i)
	{
		uint32_t letter = 0xf & alignment.cigar[cigar_i]; // 4 low bits
		uint32_t length = (0xfffffff0 & alignment.cigar[cigar_i]) >> 4; // high 28 bits i.e. 32-4=28
		if (letter == 0)
		{
			for (decltype(length) i = 0; i < length; ++i)
			{
				if (refseq[read_i] != read_iseq[query_i]) ++mismatches;
				else ++id;
				++read_i;
				++query_i;
			}
		}
		else if (letter == 1)
		{
			query_i += length;
			gaps += length;
		}
		else
		{
			read_i += length;
			gaps += length;
		}
	}

	auto rid = (double)id / ((double)mismatches + gaps + id); // %ID
	auto rcov = (double)align_len / read_iseq.length(); // %COV
	// round to 3 decimal places
	//stringstream ss;
	//ss.precision(3);
	//ss << (double)id / ((double)mismatches + gaps + id) << ' ' << (double)align_len / read_iseq.length();
	//double align_id_round = 0.0;
	//double align_cov_round = 0.0;
	//ss >> align_id_round >> align_cov_round;

	return { rid >= opts.min_id, rcov >= opts.min_cov };
} // ~is_id_cov_pass
