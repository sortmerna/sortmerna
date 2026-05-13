/*
@copyright 2016-2026 Clarity Genomics BVBA
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
              biocodz          biocodz@protonmail.com
*/

/*
 * FILE: kvdb.cpp
 * Created: Jun 05, 2018
 */
#include "kvdb.hpp"
#include "common.hpp"

#include <iostream>
#include <filesystem>
#include <memory>

#include "rocksdb/write_batch.h"
#include "rocksdb/iterator.h"

KeyValueDatabase::KeyValueDatabase(std::string const &kvdbPath)
{
	// init and open key-value database for read matches
	options.IncreaseParallelism();
	options.compression = rocksdb::kZlibCompression;
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

bool KeyValueDatabase::has(const std::string& key)
{
	std::string val;
	rocksdb::Status s = kvdb->Get(rocksdb::ReadOptions(), key, &val);
	return s.ok();
}

void KeyValueDatabase::del(const std::string& key)
{
	kvdb->Delete(rocksdb::WriteOptions(), key);
}

void KeyValueDatabase::put_batch(const std::vector<std::pair<std::string, std::string>>& kvs)
{
	rocksdb::WriteBatch batch;
	for (const auto& kv : kvs) {
		batch.Put(kv.first, kv.second);
	}
	kvdb->Write(rocksdb::WriteOptions(), &batch);
}

void KeyValueDatabase::iter_prefix(const std::string& prefix,
                                   const std::function<bool(const std::string&, const std::string&)>& fn)
{
	std::unique_ptr<rocksdb::Iterator> it(kvdb->NewIterator(rocksdb::ReadOptions()));
	for (it->Seek(prefix); it->Valid(); it->Next()) {
		auto key = it->key();
		if (key.size() < prefix.size() || std::memcmp(key.data(), prefix.data(), prefix.size()) != 0)
			break;
		std::string k(key.data(), key.size());
		auto val = it->value();
		std::string v(val.data(), val.size());
		if (!fn(k, v))
			break;
	}
}

void KeyValueDatabase::delete_prefix(const std::string& prefix)
{
	std::unique_ptr<rocksdb::Iterator> it(kvdb->NewIterator(rocksdb::ReadOptions()));
	rocksdb::WriteBatch batch;
	for (it->Seek(prefix); it->Valid(); it->Next()) {
		auto key = it->key();
		if (key.size() < prefix.size() || std::memcmp(key.data(), prefix.data(), prefix.size()) != 0)
			break;
		batch.Delete(key);
	}
	kvdb->Write(rocksdb::WriteOptions(), &batch);
}
