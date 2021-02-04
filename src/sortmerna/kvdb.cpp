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
 * FILE: kvdb.cpp
 * Created: Jun 05, 2018
 */
#include "kvdb.hpp"
#include "common.hpp"

#include <iostream>
#include <filesystem>

KeyValueDatabase::KeyValueDatabase(std::string const &kvdbPath) 
{
	// init and open key-value database for read matches
	options.IncreaseParallelism();
#if defined(_WIN32)
	options.compression = rocksdb::kXpressCompression;
#else
	options.compression = rocksdb::kZlibCompression;
#endif
	options.create_if_missing = true;
	rocksdb::Status s = rocksdb::DB::Open(options, kvdbPath, &kvdb);
	assert(s.ok());
}

/* 
 * Remove database files from the given location
 */
int KeyValueDatabase::clear(std::string dbpath)
{
	return 0;
} // ~KeyValueDatabase::clear

void KeyValueDatabase::put(std::string key, std::string val)
{
	rocksdb::Status s = kvdb->Put(rocksdb::WriteOptions(), key, val);
}

std::string KeyValueDatabase::get(std::string key)
{
	std::string val;
	rocksdb::Status s = kvdb->Get(rocksdb::ReadOptions(), key, &val);
	return val;
}