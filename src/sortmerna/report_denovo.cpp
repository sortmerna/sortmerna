/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is a free software: you can redistribute it and/or modify
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

#include "report_denovo.h"
#include "common.hpp"
#include "options.hpp"
#include "read.hpp"
#include "references.hpp"
#include "refstats.hpp"
#include "readfeed.hpp"

ReportDenovo::ReportDenovo(Runopts& opts) : Report(opts), base() {}

ReportDenovo::ReportDenovo(Readfeed& readfeed, Runopts& opts) : ReportDenovo(opts)
{
	init(readfeed, opts);
}

void ReportDenovo::init(Readfeed& readfeed, Runopts& opts)
{
	base.init(opts);
	base.init(readfeed, opts, fv, fsv, opts.aligned_pfx.string() + "_denovo", pid_str);
	openfw(); // open output files for writing
	is_zip = readfeed.orig_files[0].isZip;
	if (is_zip)	init_zip();
}

void ReportDenovo::append(int id, std::vector<Read>& reads, Runopts& opts, bool is_last)
{
	if (opts.is_paired) {
		// caclulate the index of the output file to write to
		for (int i = 0, idx = 0; i < reads.size(); ++i)
		{
			// 1 output file a (aligned reads)
			if (base.num_out == 1) {
				if (opts.is_paired_in || 
					(reads[i].n_denovo > 0 && reads[i].c_yid_ycov == 0
						&& reads[i].n_yid_ncov == 0 && reads[i].n_nid_ycov == 0))
					idx = id;
				else
					continue;
			}
			// 2 output files dp,ds (sout) | df,dr (out2)
			else if (base.num_out == 2) {
				if (opts.is_out2) {
					auto is_both_dn = (reads[0].n_denovo > 0 && reads[0].c_yid_ycov == 0
										&& reads[0].n_yid_ncov == 0 && reads[0].n_nid_ycov == 0) 
									&& (reads[1].n_denovo > 0 && reads[1].c_yid_ycov == 0
										&& reads[1].n_yid_ncov == 0 && reads[1].n_nid_ycov == 0);
					if (opts.is_paired_out && !is_both_dn)
						break;
					else if (opts.is_paired_in ||
								(reads[i].n_denovo > 0 && reads[i].c_yid_ycov == 0
									&& reads[i].n_yid_ncov == 0 && reads[i].n_nid_ycov == 0))
						idx = id * base.num_out + i;
				}
				else if (opts.is_sout) {
					auto is_both_dn = (reads[0].n_denovo > 0 && reads[0].c_yid_ycov == 0
										&& reads[0].n_yid_ncov == 0 && reads[0].n_nid_ycov == 0)
									&& (reads[1].n_denovo > 0 && reads[1].c_yid_ycov == 0
										&& reads[1].n_yid_ncov == 0 && reads[1].n_nid_ycov == 0);
					if (is_both_dn)
						idx = id * base.num_out; // both to 'dp' [0]
					else if (reads[i].n_denovo > 0 && reads[i].c_yid_ycov == 0
								&& reads[i].n_yid_ncov == 0 && reads[i].n_nid_ycov == 0)
						idx = id * base.num_out + 1; // hit to 'ds' [1]
					else
						continue; // ignore non-aligned read
				}
			}
			// 4 output files: dpf, dpr, dsf, dsr
			// being here means both is_out2 & is_sout were set => No paired_in/out
			else if (base.num_out == 4) {
				auto is_both_dn = (reads[0].n_denovo > 0 && reads[0].c_yid_ycov == 0
					&& reads[0].n_yid_ncov == 0 && reads[0].n_nid_ycov == 0)
					&& (reads[1].n_denovo > 0 && reads[1].c_yid_ycov == 0
						&& reads[1].n_yid_ncov == 0 && reads[1].n_nid_ycov == 0);
				if (is_both_dn)
					idx = id * base.num_out + i; // 0 -> apf, 1 -> apr
				else if (reads[i].n_denovo > 0 && reads[i].c_yid_ycov == 0
							&& reads[i].n_yid_ncov == 0 && reads[i].n_nid_ycov == 0)
					idx = id * base.num_out + i + 2; // 0 -> asf, 1 -> asr
				else
					continue; // ignore a non-aligned singleton
			}
			else {
				ERR("min number of output files is 1, max number of output files is 4. The current value is ", base.num_out);
				exit(1);
			}

			if (is_zip)
				base.write_a_read(fsv[idx], reads[i], vstate_out[idx], vzlib_out[idx], is_last);
			else
				base.write_a_read(fsv[idx], reads[i]);
		}
	}//~if paired
	// non-paired
	else
	{
		if (is_zip)
			base.write_a_read(fsv[id], reads[0], vstate_out[id], vzlib_out[id], is_last);
		else
			base.write_a_read(fsv[id], reads[0]);
	}
}

void ReportDenovo::merge(int num_splits)
{
	for (int i = 0; i < base.num_out; ++i) {
		openfw(i);
		for (int j = 1; j < num_splits; ++j) {
			auto idx = i + j * base.num_out;
			openfr(idx);
			fsv[i] << fsv[idx].rdbuf();
			INFO("merged ", fv[idx], " -> ", fv[i]);
			closef(idx);
			std::filesystem::remove(fv[idx]);
			INFO("deleted ", fv[idx]);
		}
		closef(i);
		strip_path_sfx(fv[i]);
	}
}