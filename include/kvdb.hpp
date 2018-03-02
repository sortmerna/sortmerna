#pragma once
/**
 * FILE: kvdb.hpp
 * Created: Nov 06, 2017 Mon
 */

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

class KeyValueDatabase {
public:
	KeyValueDatabase(std::string kvdbPath) {
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
	~KeyValueDatabase() { delete kvdb; }

	void put(std::string key, std::string val)
	{
		rocksdb::Status s = kvdb->Put(rocksdb::WriteOptions(), key, val);
	}

	std::string get(std::string key)
	{
		std::string val;
		rocksdb::Status s = kvdb->Get(rocksdb::ReadOptions(), key, &val);
		return val;
	}
private:
	rocksdb::DB* kvdb;
	rocksdb::Options options;
};
