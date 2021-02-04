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

/* 
 * file: cmd.hpp
 * created: Jan 17, 2018 Wed
 *
 * Interactive command-line session:
 *    - query records in Key-Value database
 */

#pragma once

#include <string>

// forward
struct Runopts;

enum CMD { EXIT, READ, INDEX };

class CmdSession
{
public:
	CmdSession(){}
	void run(Runopts & opts);
private:
	void cmdRead(Runopts & opts, std::string & cmd);
	void cmdIndex(Runopts & opts, std::string & cmd);
	void cmdTest(Runopts & opts, std::string & cmd);
	void cmd_max_ref_part(Runopts & opts, std::string & cmd); // ref idx=0 part=1
};