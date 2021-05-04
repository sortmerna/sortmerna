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

#include "report_fx_other.h"
#include "readfeed.hpp"
#include "options.hpp"
#include "read.hpp"

ReportFxOther::ReportFxOther(Runopts& opts) : Report(opts), base() {}
ReportFxOther::ReportFxOther(Readfeed& readfeed, Runopts& opts) : ReportFxOther(opts) { init(readfeed, opts); }

void ReportFxOther::init(Readfeed& readfeed, Runopts& opts)
{
	base.init(opts);
	base.init(readfeed, opts, fv, fsv, opts.other_pfx.string(), pid_str);
	openfw(opts.dbg_level); // open output files for writing
	is_zip = readfeed.orig_files[0].isZip;
	if (is_zip)	init_zip();
} // ~ReportFxOther::init

void ReportFxOther::append(int id, std::vector<Read>& reads, Runopts& opts, bool is_last)
{
	if (opts.is_paired) {
		if (reads[0].is_hit && reads[1].is_hit)
			return; // both aligned - nothing to write

		// caclulate the index of the output file to write to
		for (int i = 0, idx = 0; i < reads.size(); ++i, idx = 0)
		{
			// 1 output file (non-aligned reads)
			if (base.num_out == 1) {
				if (opts.is_paired_in)
					if (reads[0].is_hit || reads[1].is_hit)
						continue; // both go to aligned
					else
						idx = id; // none hit -> other
				else if (opts.is_paired_out || !reads[i].is_hit)
					idx = id;
				else
					continue;
			}
			// 2 output files ap,as (sout) | af,ar (out2)
			else if (base.num_out == 2) {
				if (opts.is_out2) {
					if (opts.is_paired_in) {
						if (reads[0].is_hit || reads[1].is_hit)
							break; // if not both non-aligned -> aligned
						else
							idx = id * base.num_out + i;
					}
					else if (opts.is_paired_out || !reads[i].is_hit)
						idx = id * base.num_out + i;
					else
						continue; // ignore aligned
				}
				else if (opts.is_sout) {
					if (!reads[0].is_hit && !reads[1].is_hit)
						idx = id * base.num_out; // both to 'op' [0]
					else if (!reads[i].is_hit)
						idx = id * base.num_out + 1; // non-aligned to 'os' [1]
					else
						continue; // ignore an aligned singleton
				}
			}
			// 4 output files: apf, apr, asf, asr
			// being here means both is_out2 & is_sout were set => No paired_in/out
			else if (base.num_out == 4) {
				if (!reads[0].is_hit && !reads[1].is_hit)
					idx = id * base.num_out + i; // 0 -> opf, 1 -> opr
				else if (!reads[i].is_hit)
					idx = id * base.num_out + i + 2; // 0 -> osf, 1 -> osr
				else
					continue; // ignore an aligned singleton
			}
			else {
				ERR("min number of output files is 1, max number of output files is 4. The current value is ", base.num_out);
				exit(1);
			}

			if (is_zip)
				base.write_a_read(fsv[idx], reads[i], vstate_out[idx], vzlib_out[idx], is_last, opts.dbg_level);
			else
				base.write_a_read(fsv[idx], reads[i], opts.dbg_level);
		}
	}//~if paired
	// non-paired
	else
	{
		if (!reads[0].is_hit)
		{
			if (is_zip)
				base.write_a_read(fsv[0], reads[0], vstate_out[0], vzlib_out[0], is_last, opts.dbg_level);
			else
				base.write_a_read(fsv[0], reads[0], opts.dbg_level);
		}
	}
} // ~ReportFxOther::append

ReportFxBase& ReportFxOther::getBase()
{
	return base;
}