/**
 * FILE: read.cpp
 * Created: Nov 26, 2017 Sun
 */

#include "read.hpp"


 // initialize Smith-Waterman scoring matrix for genome sequences
void Read::initScoringMatrix(Runopts & opts)
{
	int l, k, m;
	int8_t val = 0;
	for (l = k = 0; l < 4; ++l)
	{
		for (m = 0; m < 4; ++m) {
			val = l == m ? opts.match : opts.mismatch; // weight_match : weight_mismatch (must be negative)
			scoring_matrix.push_back(val);
		}
		scoring_matrix.push_back(opts.score_N); // ambiguous base
	}
	for (m = 0; m < 5; ++m) {
		scoring_matrix.push_back(opts.score_N); // ambiguous base
	}
}

void Read::restoreFromDb(KeyValueDatabase & kvdb)
{
	std::string readstr = kvdb.get(std::to_string(id));
}

void Read::unmarshallJson(KeyValueDatabase & kvdb)
{
	printf("Read::unmarshallJson: Not yet Implemented\n");
}