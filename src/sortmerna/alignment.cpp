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
uint32_t findMinIndex(Read & read);

//#define DEBUG_ALIGN
#ifdef DEBUG_ALIGN
#define DBG_READ_ID 10017
#define DBG_IDX_PART 0

std::string LOGF = "C:/a01_projects/clarity_genomics/logs/debug.log";
void debug(std::string text)
{
	std::ofstream log;
	log.open(LOGF, std::ios::app);
	log << text;
	log.close();
}

void debug_seq_kmer_freq_vec(int readid, int part, vector<uint32pair> & seq_kmer_freq_vec)
{
	std::stringstream ss;
	ss << "=== seq_kmer_freq_vec ===" << std::endl;
	ss << "Read ID: " << readid << " Part: " << part << std::endl;
	ss << "Reference #, Kmer frequency" << std::endl;
	ss << "---------------------------" << std::endl;
	for (auto pair : seq_kmer_freq_vec)
		ss << pair.first << " " << pair.second << std::endl;
	debug(ss.str());
}

void debug_hits_on_genome(vector<uint32pair> & hits_on_genome)
{
	std::stringstream ss;
	ss << "=== hits_on_genome ===" << std::endl;
	ss << "Position on reference, position on read" << std::endl;
	ss << "---------------------------------------" << std::endl;
	for (uint32pair pair : hits_on_genome)
	{
		ss << pair.first << "  " << pair.second << std::endl;
	}
	debug(ss.str());
}

void debug_hits_on_genome2(std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> hits)
{
	std::stringstream ss;
	ss << "=== hits_on_genome ===" << std::endl;
	ss << "Reference, Position on reference, Position on read" << std::endl;
	ss << "---------------------------------------" << std::endl;
	for (auto triple : hits)
	{
		ss << std::get<0>(triple) << "  " << std::get<1>(triple) << "  " << std::get<2>(triple) << std::endl;
	}
	debug(ss.str());
}

void debug_id_win_hits(Read & read)
{
	std::stringstream ss;
	ss << "=== read.id_win_hits ===" << std::endl;
	ss << "kmer id, kmer position on read" << std::endl;
	for (auto hit : read.id_win_hits)
		ss << hit.id << " " << hit.win << std::endl;
	debug(ss.str());
}
#endif

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
} // ~find_lis

/* 
 * called on each idx * part * read * strand * [1..max opts.skiplengths[index_num].size (3 by default)] 
 */
void compute_lis_alignment
	(
		Read & read, Runopts & opts, Index & index, References & refs, Readstats & readstats, Refstats & refstats,
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
		// frequency map of k-mer occurrences on the references i.e.
		// <reference number : number of the k-mer occurrences>
		map<uint32_t, uint32_t> seq_kmer_freq_map;
		// vector to hold frequency map content for Sorting (map cannot be sorted)
		vector<uint32pair> seq_kmer_freq_vec;
		map<uint32_t, uint32_t>::iterator map_it;
		uint32_t max_seq = 0;
		uint32_t max_occur = 0;

#ifdef DEBUG_ALIGN
		if (read.id == DBG_READ_ID && !read.reversed && refs.part == DBG_IDX_PART)
			debug_id_win_hits(read);
#endif

		// STEP 2: Find all candidate references by using Read's kmers hits information.
		//         For every reference, compute the number of kmers' hits belonging to it
		for (auto hit : read.id_win_hits)
		{
			seq_pos* positions_tbl_ptr = index.positions_tbl[hit.id].arr;
			// loop all positions of id
			for (uint32_t j = 0; j < index.positions_tbl[hit.id].size; j++)
			{
				uint32_t seq = positions_tbl_ptr++->seq;
				if ((map_it = seq_kmer_freq_map.find(seq)) != seq_kmer_freq_map.end())
					map_it->second++; // sequence already in the map, increment its frequency value
				else
					seq_kmer_freq_map[seq] = 1; // sequence not in the map, add it
			}
		}

		// copy frequency map to vector for sorting
		// consider only candidate references that have enough seed hits
		for (auto freq_pair: seq_kmer_freq_map)
		{
			if (freq_pair.second >= (uint32_t)opts.seed_hits)
				seq_kmer_freq_vec.push_back(freq_pair);
		}

		seq_kmer_freq_map.clear();

		// sort sequences by frequency in descending order
		std::sort(seq_kmer_freq_vec.begin(), seq_kmer_freq_vec.end(), 
			[](std::pair<uint32_t, uint32_t> e1, std::pair<uint32_t, uint32_t> e2) {
			if (e1.second == e2.second) 
				return e1.first ASCENDING e2.first; // order references ascending for equal frequencies (originally: descending)
			return e1.second DESCENDING e2.second; // order frequencies descending
		});

#ifdef DEBUG_ALIGN
		if (read.id == DBG_READ_ID && !read.reversed && refs.part == DBG_IDX_PART)
			debug_seq_kmer_freq_vec(read.id, refs.part, seq_kmer_freq_vec);
#endif

		// STEP 3: for each reference sequence candidate,
		//		starting from highest scoring.
		for (uint32_t k = 0; k < seq_kmer_freq_vec.size(); k++)
		{
			// the maximum scoring alignment has been found - stop searching for more alignments
			if (opts.num_best_hits != 0 && read.max_SW_score == opts.num_best_hits) {
#ifdef DEBUG_ALIGN
				if (read.id == DBG_READ_ID && !read.reversed && refs.part == DBG_IDX_PART)
				{
					std::stringstream ss;
					ss << " Breaking OUT num_best_hits: " << opts.num_best_hits << " max_SW_score: " << read.max_SW_score << std::endl;
					debug(ss.str());
				}
#endif
				break;
			}

			max_seq = seq_kmer_freq_vec[k].first;
			max_occur = seq_kmer_freq_vec[k].second;
              
			// not enough window hits, try to collect more hits or next read
			if (max_occur < (uint32_t)opts.seed_hits) {
#ifdef DEBUG_ALIGN
				if (read.id == DBG_READ_ID && !read.reversed && refs.part == DBG_IDX_PART)
				{
					std::stringstream ss;
					ss << " Breaking OUT max_occur: " << max_occur << " seed_hits: " << opts.seed_hits << std::endl;
					debug(ss.str());
				}
#endif
				break;
			}

			// update number of reference sequences remaining to check
			// only decrement read.best if the current ref sequence to check
			// has a lower seed count than the previous one
			if ( opts.min_lis > 0 && aligned && k > 0 && max_occur < seq_kmer_freq_vec[k - 1].second )
			{
				read.best--;
				if (read.best < 1) {
#ifdef DEBUG_ALIGN
					if (read.id == DBG_READ_ID && !read.reversed && refs.part == DBG_IDX_PART)
					{
						std::stringstream ss;
						ss << " Breaking OUT read.best: " << read.best << " opts.min_lis: " << opts.min_lis << std::endl;
						debug(ss.str());
					}
#endif
					break;
				}
			}

			// check if the maximum number of alignments per read
			// (--num_alignments INT) have been output
			if (opts.num_alignments > 0 && read.num_alignments <= 0)
			{
#ifdef DEBUG_ALIGN
				if (read.id == DBG_READ_ID && !read.reversed && refs.part == DBG_IDX_PART)
				{
					std::stringstream ss;
					ss << " Breaking OUT opts.num_alignments: " << opts.num_alignments << " read.num_alignments" << read.num_alignments << std::endl;
					debug(ss.str());
				}
#endif
				break;
			}

			// STEP 4: 
			// collect all the genome positions belonging to the reference candidate
			vector<uint32pair> hits_on_genome;
#ifdef DEBUG_ALIGN
			std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> hits_on_genome2;
			if (read.id == DBG_READ_ID && !read.reversed && refs.part == DBG_IDX_PART)
			{
				std::stringstream ss;
				ss << "=== Analyzing Reference ID: " << max_seq << " Frequency: " << max_occur << " ref_kmer_freq index: " << k << std::endl;
				debug(ss.str());
			}
#endif
			for ( auto hit: read.id_win_hits )
			{
				uint32_t num_hits = index.positions_tbl[hit.id].size;
				seq_pos* positions_tbl_ptr = index.positions_tbl[hit.id].arr;
				// loop through every position of id
				for (uint32_t j = 0; j < num_hits; j++)
				{
					if (positions_tbl_ptr->seq == max_seq)
					{
						hits_on_genome.push_back(uint32pair(positions_tbl_ptr->pos, hit.win));
#ifdef DEBUG_ALIGN
						hits_on_genome2.push_back(std::make_tuple(max_seq, positions_tbl_ptr->pos, hit.win));
#endif
					}
					positions_tbl_ptr++;
				}
			}

			// sort the positions in ascending order
			std::sort(hits_on_genome.begin(), hits_on_genome.end(), [](uint32pair e1, uint32pair e2) {
				if (e1.first == e2.first) 
					return (e1.second ASCENDING e2.second); // order references ascending for equal reference positions
				return (e1.first ASCENDING e2.first);
			}); // smallest

#ifdef DEBUG_ALIGN
			if (read.id == DBG_READ_ID && !read.reversed && refs.part == DBG_IDX_PART)
			{
				//debug_hits_on_genome(hits_on_genome);
				debug_hits_on_genome2(hits_on_genome2);
			}

#endif
			// iterate over the set of hits, output windows of
			// length == read which have at least ratio hits
			vector<uint32pair>::iterator it3 = hits_on_genome.begin();
			deque<uint32pair> vi_read;
			deque<uint32pair>::iterator deq;

			// STEP 5: run a sliding window of read's length across
			// the genome, and search for windows with enough k-mer hits
			uint32_t lcs_ref_start = 0;
			uint32_t lcs_que_start = 0;
			uint32_t begin = it3->first;
                       
			// TODO: remove this 'while'? - Always does a single iteration because of line '++it3'. 
			// It has 3 'break' instructions though. Convoluted.
			while (it3 != hits_on_genome.end())
			{
				uint32_t stop = begin + read.sequence.length() - refstats.lnwin[index.index_num] + 1; // lnwin_index_num
				bool push = false;
				while ( it3 != hits_on_genome.end() && it3->first <= stop )
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

							// put read into 04 encoding before SSW
							if (read.is03) read.flip34();
                       
							// create profile for read
							s_profile* profile = 0;
							profile = ssw_init((int8_t*)(&read.isequence[0] + align_que_start), (align_length - head - tail), &read.scoring_matrix[0], 5, 2);

							s_align* result = 0;

							result = ssw_align(
								profile,
								(int8_t*)refs.buffer[max_seq].sequence.c_str() + align_ref_start - head,
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
							if ( result != 0 && result->score1 > refstats.minimal_score[index.index_num] )
									aligned = true;

							result->index_num = index.index_num;
							result->part = index.part;
							result->strand = !read.reversed; // flag whether the alignment was done on a forward or reverse strand

							// STEP 8: alignment succeeded, output (--all) or store (--best)
							// the alignment and go to next alignment
							if (aligned)
							{
								// read has not been yet mapped, set bit to true for this read
								// (this is the Only place where read_hits must be modified)
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
									if (read.hits_align_info.alignv.size() > 0)
									{
										uint32_t smallest_score_index = read.hits_align_info.min_index;
										uint32_t highest_score_index = read.hits_align_info.max_index;
										uint32_t hits_size = read.hits_align_info.alignv.size();

										// number of alignments stored per read < 'num_best_hits_gv', 
										// add alignment to array without comparison to other members of array
										if (opts.num_best_hits == 0 || hits_size < (uint32_t)opts.num_best_hits)
										{
											// add alignment
											read.hits_align_info.alignv.push_back( copyAlignment(result) );
											++hits_size;

											// read alignments are filled to max size, find the smallest
											// alignment score and set the smallest_score_index
											// (this is not done when num_best_hits_gv == 0 since
											// we want to output all alignments for some --min_lis)
											if (read.hits_align_info.alignv.size() == (uint32_t)opts.num_best_hits)
											{
												read.hits_align_info.min_index = findMinIndex(read);
											}

											// update the index position of the first occurrence of the
											// highest alignment score
											if (result->score1 > read.hits_align_info.alignv[highest_score_index].score1)
												read.hits_align_info.max_index = hits_size - 1;

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

											// replace an old smallest scored alignment with the new one
											read.hits_align_info.alignv[smallest_score_index] = copyAlignment(result);

											read.hits_align_info.min_index = findMinIndex(read);
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

										read.hits_align_info.alignv.push_back( copyAlignment(result) );

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
									// add alignment to the read. TODO: check how this affects the old logic
									read.hits_align_info.alignv.push_back(copyAlignment(result));

									// the maximum possible score for this read has been found
									if (result->score1 == max_SW_score) read.max_SW_score++;

									// update number of alignments to output per read
									if (opts.num_alignments > 0) read.num_alignments--; // TODO: why decrement?

									// get the edit distance between reference and read (serves for
									// SAM output and computing %id and %query coverage)
									uint32_t id = 0;
									uint32_t mismatches = 0;
									uint32_t gaps = 0;
									read.calcMismatchGapId(refs, read.hits_align_info.alignv.size()-1, mismatches, gaps, id);

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
									if ( align_id_round >= opts.align_id && align_cov_round >= opts.align_cov && read_to_count)
									{
										++readstats.total_reads_mapped_cov;
										read_to_count = false;

										// do not output read for de novo OTU clustering
										// it passed the %id/coverage thersholds
										if (opts.de_novo_otu) read.hit_denovo = false;
									}

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
					if (it3 != hits_on_genome.end()) // TODO: seems Always false
						begin = it3->first; // TODO: seems never reached
					else break;
				}
				else
				{
					begin = (vi_read.front()).first;
				}
			}//~for all reference sequence length                   
		}//~for all of the reference sequence candidates
	}//if ( readhitf || readhitr > ratio )
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
	ret_align.ref_seq = pAlign->ref_seq;
	ret_align.score1 = pAlign->score1;
	ret_align.strand = pAlign->strand;

	return ret_align;
}

uint32_t findMinIndex(Read & read)
{
	uint32_t smallest_score = read.hits_align_info.alignv[0].score1;
	uint32_t index = 0;
	for (int i = 0; i < read.hits_align_info.alignv.size(); ++i)
	{
		if (read.hits_align_info.alignv[i].score1 < smallest_score)
		{
			smallest_score = read.hits_align_info.alignv[i].score1;
			index = i;
		}
	}
	return index;
}
