/**
 * @file alignment.cpp
 * @brief File containing functions for sequence alignment.
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2012-17 Bonsai Bioinformatics Research Group
 * @copyright 2014-17 Knight Lab, Department of Pediatrics, UCSD, La Jolla
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
 * @contributors Jenya Kopylova   jenya.kopylov@gmail.com
 *               Laurent Noé      laurent.noe@lifl.fr
 *               Pierre Pericard  pierre.pericard@lifl.fr
 *               Daniel McDonald  wasade@gmail.com
 *               Mikaël Salson    mikael.salson@lifl.fr
 *               Hélène Touzet    helene.touzet@lifl.fr
 *               Rob Knight       robknight@ucsd.edu
 */

#include <sstream>

#include "alignment.hpp"
#include "read.hpp"
#include "index.hpp"
#include "refstats.hpp"
#include "references.hpp"
#include "readstats.hpp"


bool
smallest(const mypair &a, const mypair &b)
{
	if (a.first == b.first) return (a.second < b.second);
	return (a.first < b.first);
}

bool
largest(const mypair &a, const mypair &b)
{
	if (a.first == b.first) return (a.second > b.second);
	return (a.first > b.first);
}

/*
 *
 * FUNCTION   : find_lis()
 * see alignment.hpp for documentation
 *************************************/
void find_lis(
	deque<pair<uint32_t, uint32_t>> &a,
	vector<uint32_t> &b, uint64_t readn
)
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
		for (u = 0, v = b.size() - 1; u < v;)
		{
			int c = (u + v) / 2;
			if (a[b[c]].second < a[i].second)
				u = c + 1;
			else
				v = c;
		}

		// Update b if new value is smaller then previously referenced value
		if (a[i].second < a[b[u]].second)
		{
			if (u > 0) p[i] = b[u - 1];
			b[u] = i;
		}
	}


	for (u = b.size(), v = b.back(); u--; v = p[v]) b[u] = v;
}


void compute_lis_alignment
	(
		Read & read, Runopts & opts, Index & index, References & refs, Readstats & readstats, Refstats & refstats, Output & output,
		bool & search,
		uint32_t max_SW_score,
		bool& read_to_count
	)
{
	// boolean set to true if SW alignment succeeded between
	// the read and a candidate reference sequence
	bool aligned = false;

	// STEP 1: the number of matching windows on the read to the
	// reference database is greater than the threshold,
	// continue analysis of read
	if (read.readhit >= (uint32_t)opts.seed_hits) // default seed_hits_gv = 2
	{
		// map<seq, number of occurrences> most_frequent_seq_t
		map<uint32_t, uint32_t> most_frequent_seq_t;
		map<uint32_t, uint32_t>::iterator map_it;
		uint32_t max_seq = 0;
		uint32_t max_occur = 0;

		// STEP 2: for every reference sequence, compute the number of
		// window hits belonging to it
		for (uint32_t i = 0; i < read.id_win_hits.size(); i++)
		{
			uint32_t _id = read.id_win_hits[i].id;
			// number of entries in the positions table for this id
			uint32_t num_hits = index.positions_tbl[_id].size;
			// pointer to the seq_pos array in the positions table for this id
			seq_pos* positions_tbl_ptr = index.positions_tbl[_id].arr;
			// loop through every position of id
			for (uint32_t j = 0; j < num_hits; j++)
			{
				uint32_t seq = positions_tbl_ptr++->seq;
				// sequence already exists in the map, increment it's value
				if ((map_it = most_frequent_seq_t.find(seq)) != most_frequent_seq_t.end())
					map_it->second++;
				// sequence doesn't exist, add it
				else most_frequent_seq_t[seq] = 1;
			}
		}

		// <mypair> = <number of occurrences of a sequence, index of sequence>
		vector<mypair> most_frequent_seq;
		// copy list of occurrences from map to vector for sorting
		for (map_it = most_frequent_seq_t.begin(); map_it != most_frequent_seq_t.end(); map_it++)
		{
			// pass candidate reference sequences for further analyses
			// only if they have enough seed hits
			if (map_it->second >= (uint32_t)opts.seed_hits)
				most_frequent_seq.push_back(mypair(map_it->second, map_it->first));
		}
		most_frequent_seq_t.clear();
		// sort the highest scoring sequences to the head of the array
		sort(most_frequent_seq.begin(), most_frequent_seq.end(), largest);

		// STEP 3: for each reference sequence candidate
		// (starting from highest scoring)
		for (uint32_t k = 0; k < most_frequent_seq.size(); k++)
		{
			// the maximum scoring alignment has been found,
			// do not search for anymore alignments
			if ((opts.num_best_hits != 0) && (read.max_SW_score == opts.num_best_hits)) break;
			max_occur = most_frequent_seq[k].first;
			max_seq = most_frequent_seq[k].second;
              
			// not enough window hits, try to collect more hits or go to next read
			if (max_occur < (uint32_t)opts.seed_hits) break;

			// update number of reference sequences remaining to check
			if ((opts.min_lis > 0) && aligned && (k > 0))
			{
				// only decrement best_x if the next ref sequence to check
				// has a lower seed count than the previous one
				if (max_occur < most_frequent_seq[k - 1].first)
				{
					{
						read.best--;
					}
					if (read.best < 1) break;
				}
			}

			// check if the maximum number of alignments per read
			// (--num_alignments INT) have been output
			if (opts.num_alignments > 0)
			{
				if (read.num_alignments <= 0) break;
			}

			// STEP 4: collect all the genome positions belonging to the
			// reference candidate from the table of positions computed
			// during indexing
			vector<mypair> hits_on_genome;
			uint32_t count = 0;
			for (uint32_t i = 0; i < read.id_win_hits.size(); i++)
			{
				uint32_t _id = read.id_win_hits[i].id;
				uint32_t num_hits = index.positions_tbl[_id].size;
				seq_pos* positions_tbl_ptr = index.positions_tbl[_id].arr;
				// loop through every position of id
				for (uint32_t j = 0; j < num_hits; j++)
				{
					if (positions_tbl_ptr->seq == max_seq)
					{
						count++;
						hits_on_genome.push_back(mypair(positions_tbl_ptr->pos, read.id_win_hits[i].win)); // id_win_hits[i]
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
                       
			while (it3 != hits_on_genome.end())
			{
				uint32_t stop = begin + read.sequence.length() - refstats.lnwin[index.index_num] + 1; // lnwin_index_num
				bool push = false;
				while ((it3 != hits_on_genome.end()) && (it3->first <= stop))
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
				if (!push && aligned) goto pop;
				else aligned = false;
#endif
#ifdef HEURISTIC1_OFF
				aligned = false;
#endif                              
				// enough windows at this position on genome to search for LIS
				if (vi_read.size() >= (uint32_t)opts.seed_hits)
				{
					vector<uint32_t> list;
					find_lis(vi_read, list, read.id); // TODO: read.id is not used in this function.
#ifdef HEURISTIC1_OFF
					uint32_t list_n = 0;
					do
					{
#endif                                      
						// LIS long enough to perform Smith-Waterman alignment
						if (list.size() >= (uint32_t)opts.seed_hits)
						{
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
							uint32_t reflen = refs.buffer[max_seq].sequence.length(); // reference_seq_len[max_seq]
							uint32_t edges = 0;
							if (opts.as_percent)
								edges = (((double)opts.edges / 100.0)*read.sequence.length()); // readlen
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
								if (reflen < read.sequence.length()) // readlen
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

							if (read.is04) read.flip34(opts);
                       
							// create profile for read
							s_profile* profile = 0;
							profile = ssw_init((int8_t*)(&read.isequence[0] + align_que_start), (align_length - head - tail), &read.scoring_matrix[0], 5, 2);

							s_align* result = 0;

							result = ssw_align(
								profile,
								(int8_t*)refs.buffer[max_seq].sequence.c_str() + align_ref_start - head, // TODO: review this
								align_length,
								opts.gap_open,
								opts.gap_extension,
								2,
								refstats.minimal_score[index.index_num], // minimal_score_index_num
								0,
								0
							);

							// deallocate memory for profile, no longer needed
							if (profile != 0) init_destroy(&profile);

							// check alignment satisfies all thresholds
							if (result != 0)
							{
								if (result->score1 > refstats.minimal_score[index.index_num]) aligned = true;
							}

							result->index_num = index.index_num;
							result->part = index.part;
							result->strand = !read.reversed; // flag whether the alignment was done on a forward or reverse strand

							// STEP 8: alignment succeeded, output (--all) or store (--best)
							// the alignment and go to next alignment
							if (aligned)
							{
								// read has not been yet mapped, set bit to true for this read
								// (this is the _only_ place where read_hits must be modified)
								if (!read.hit)
								{
									read.hit = true;
									readstats.total_reads_mapped++;
									readstats.reads_matched_per_db[index.index_num]++;
								}

								// add the offset calculated by the LCS (from the beginning of the sequence)
								// to the offset computed by SW alignment
								result->ref_begin1 += (align_ref_start - head);
								result->ref_end1 += (align_ref_start - head);
								result->read_begin1 += align_que_start;
								result->read_end1 += align_que_start;
								result->readlen = read.sequence.length();
								result->ref_seq = max_seq; // TODO: Monitor - moved here from (min_lis_gv > -1)

								// update best alignment
								if (opts.min_lis > -1)
								{
									// an alignment for this read already exists
									//if (alignment != read_hits_align_info.end())
									if (read.hits_align_info.size > 0)
									{
										uint32_t smallest_score_index = read.hits_align_info.min_index;
										uint32_t highest_score_index = read.hits_align_info.max_index;
										uint32_t array_size = read.hits_align_info.size;
										uint32_t array_max_size = read.hits_align_info.max_size;

										// number of alignments stored per read < num_best_hits_gv, 
										// add alignment to array without comparison to other members of array
										if ((opts.num_best_hits == 0) || (array_size < (uint32_t)opts.num_best_hits))
										{
											// all slots have been filled, find slot with smallest
											// alignment score and set the smallest_score_index
											// (this is not done when num_best_hits_gv == 0 since
											// we want to output all alignments for some --min_lis)
											if (array_size == (uint32_t)opts.num_best_hits)
											{
												uint32_t smallest_score = 1000000;
												//s_align *this_alignment = alignment->second.ptr;
												for (int p = 0; p < opts.num_best_hits; p++)
												{
													if (read.hits_align_info.alignv[p].score1  < smallest_score) // this_alignment[p].score1
													{
														smallest_score = read.hits_align_info.alignv[p].score1;
														smallest_score_index = p;
													}
												}
												read.hits_align_info.min_index = smallest_score_index; // alignment->second.min_index
											}

											// update the index position of the first occurrence of the
											// highest alignment score
											if (result->score1 > read.hits_align_info.alignv[highest_score_index].score1)
												read.hits_align_info.max_index = array_size - 1;

											// the maximum possible score for this read has been found
											if (result->score1 == max_SW_score) read.max_SW_score++;

											// free result
											free(result);
											result = NULL;
										}//~if (array_size < num_best_hits_gv)

										// all num_best_hits_gv slots have been filled,
										// replace the alignment having the lowest score
										else if (result->score1 > read.hits_align_info.alignv[smallest_score_index].score1)
										{
											// update max_index to the position of the first occurrence
											// of the highest scoring alignment
											if (result->score1 > read.hits_align_info.alignv[highest_score_index].score1)
												read.hits_align_info.max_index = smallest_score_index;

											// decrement number of reads mapped to database
											// with lower score
											readstats.reads_matched_per_db[read.hits_align_info.alignv[smallest_score_index].index_num]--;

											// increment number of reads mapped to database with higher score
											readstats.reads_matched_per_db[index.index_num]++;

											s_align2 align;
											for (int i = 0; i < result->cigarLen; i++)
												align.cigar.push_back(*result->cigar++);

											align.index_num = result->index_num;
											align.part = result->part;
											align.readlen = result->readlen;
											align.read_begin1 = result->read_begin1;
											align.read_end1 = result->read_end1;
											align.ref_begin1 = result->ref_begin1;
											align.ref_end1 = result->ref_end1;
											align.ref_seq = result->ref_seq;
											align.score1 = result->score1;
											align.strand = result->strand;
											read.hits_align_info.alignv[smallest_score_index] = align; // *result

											// find the new smallest_score_index
											uint32_t smallest_score = 1000000;

											for (int p = 0; p < opts.num_best_hits; p++)
											{
												if (read.hits_align_info.alignv[p].score1 < smallest_score)
												{
													smallest_score = read.hits_align_info.alignv[p].score1;
													smallest_score_index = p;
												}
											}

											read.hits_align_info.min_index = smallest_score_index; // alignment->second.min_index
											// the maximum possible score for this read has been found
											if (result->score1 == max_SW_score) read.max_SW_score++;
											// free result, except the cigar (now new cigar)
											free(result);
											result = NULL;
										}
										else
										{
											// new alignment has a lower score, destroy it
											if (result != NULL) free(result);
										}
									}
									// an alignment for this read doesn't exist, add the first alignment
									else
									{
										// maximum size of s_align array
										uint32_t max_size = 0;
										// create new instance of alignments
										if ((opts.num_best_hits > 0) && (opts.num_best_hits < BEST_HITS_INCREMENT + 1))
											max_size = opts.num_best_hits;
										else 
											max_size = BEST_HITS_INCREMENT;

										s_align2 align;
										for (int i = 0; i < result->cigarLen; i++)
											align.cigar.push_back(*result->cigar++);
										align.index_num = result->index_num;
										align.part = result->part;
										align.readlen = result->readlen;
										align.read_begin1 = result->read_begin1;
										align.read_end1 = result->read_end1;
										align.ref_begin1 = result->ref_begin1;
										align.ref_end1 = result->ref_end1;
										align.ref_seq = result->ref_seq;
										align.score1 = result->score1;
										align.strand = result->strand;
										read.hits_align_info.alignv.push_back(align); // *result

										// the maximum possible score for this read has been found
										if (result->score1 == max_SW_score) read.max_SW_score++;

										// free result, except the cigar
										free(result);
										result = NULL;
									}
								}
								// output the Nth alignment (set by --num_alignments [INT] parameter)
								else if (opts.num_alignments > -1)
								{
									// update number of alignments to output per read
									if (opts.num_alignments > 0) read.num_alignments--; // TODO: why decrement?

									// get the edit distance between reference and read (serves for
									// SAM output and computing %id and %query coverage)
									double id = 0;
									char to_char[5] = { 'A','C','G','T','N' };
									const char* ref_seq_ptr = refs.buffer[max_seq].sequence.data(); // reference_seq  [(2 * (int)max_seq) + 1]
									const char* read_seq_ptr = read.isequence.data(); // myread
									int32_t qb = result->ref_begin1;
									int32_t pb = result->read_begin1;
									uint32_t mismatches = 0;
									uint32_t gaps = 0;

									for (uint32_t c2 = 0; c2 < result->cigarLen; ++c2)
									{
										uint32_t letter = 0xf & *(result->cigar + c2);
										uint32_t length = (0xfffffff0 & *(result->cigar + c2)) >> 4;
										if (letter == 0)
										{
											for (uint32_t p = 0; p < length; ++p)
											{
												if ((char)to_char[(int)*(ref_seq_ptr + qb)] != (char)to_char[(int)*(read_seq_ptr + pb)]) ++mismatches;
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

									int32_t align_len = abs(result->read_end1 + 1 - result->read_begin1);
									int32_t total_pos = mismatches + gaps + id;
									stringstream ss;
									ss.precision(3);
									ss << (double)id / total_pos << ' ' << (double)align_len / read.sequence.length();
									double align_id_round = 0.0;
									double align_cov_round = 0.0;
									ss >> align_id_round >> align_cov_round;

									// the alignment passed the %id and %query coverage threshold
									// output it (SAM, BLAST and FASTA/Q)
									if ((align_id_round >= opts.align_id) && (align_cov_round >= opts.align_cov) && read_to_count)
									{
										readstats.total_reads_mapped_cov++;
										read_to_count = false;

										// do not output read for de novo OTU clustering
										// (it passed the %id/coverage thersholds)
										if (opts.de_novo_otu) read.hit_denovo = !read.hit_denovo; // read_hits_denovo[readn].flip()
									}

									// add alignment information to the read. TODO: check how this affects the old logic
									s_align2 align;
									for (int i = 0; i < result->cigarLen; i++)
										align.cigar.push_back(*result->cigar++);
									align.index_num = result->index_num;
									align.part = result->part;
									align.readlen = result->readlen;
									align.read_begin1 = result->read_begin1;
									align.read_end1 = result->read_end1;
									align.ref_begin1 = result->ref_begin1;
									align.ref_end1 = result->ref_end1;
									align.ref_seq = result->ref_seq;
									align.score1 = result->score1;
									align.strand = result->strand;
									read.hits_align_info.alignv.push_back(align);

									// the maximum possible score for this read has been found
									if (result->score1 == max_SW_score) read.max_SW_score++;

									if (result != 0) free(result); // free alignment info
								}//~if output all alignments

								// continue to next read (do not need to collect more seeds using another pass)
								search = false;

								 // maximum score possible for the read has been reached,
								 // stop searching for further matches
								if ((opts.num_best_hits != 0) && (read.max_SW_score == opts.num_best_hits)) break;

								// stop search after the first num_alignments_gv alignments
								// for this read
								if (opts.num_alignments > 0)
								{
									// go to next read (do not search for further alignments)
									if (read.num_alignments <= 0) break; // num_alignments_x[readn]
								}
							}//~if read aligned
							else // the read did not align
							{
								if (result != 0) free(result); // free alignment info
							}
						}//~if LIS long enough                               
#ifdef HEURISTIC1_OFF
					} while ((it3 == hits_on_genome.end()) && (++list_n < list.size()));
#endif                                              
				}//~if enough window hits                                                
			pop:

				// get the next candidate reference position 
				if (!vi_read.empty()) vi_read.pop_front();

				if (vi_read.empty())
				{
					if (it3 != hits_on_genome.end())
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

	if (read.is04) read.flip34(opts);

} // ~compute_lis_alignment
