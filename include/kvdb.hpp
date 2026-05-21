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

/**
 * FILE: kvdb.hpp
 * Created: Nov 06, 2017 Mon
 */

#pragma once

#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

class KeyValueDatabase {
public:
	KeyValueDatabase(std::string const &kvdbPath);
	~KeyValueDatabase() { delete kvdb; }

	void put(std::string key, std::string val);
	std::string get(std::string key);
	bool has(const std::string& key);
	void del(const std::string& key);
	// Atomically apply a batch of put operations.
	void put_batch(const std::vector<std::pair<std::string, std::string>>& kvs);
	// Iterate all keys with the given prefix. fn receives (key, value).
	// Return false from fn to stop iteration early.
	void iter_prefix(const std::string& prefix,
	                 const std::function<bool(const std::string&, const std::string&)>& fn);
	// Delete every key with the given prefix.
	void delete_prefix(const std::string& prefix);
	int clear(std::string dbPath);
private:
	rocksdb::DB* kvdb;
	rocksdb::Options options;
};
