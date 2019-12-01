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
	KeyValueDatabase(std::string const &kvdbPath);
	~KeyValueDatabase() { delete kvdb; }

	void put(std::string key, std::string val);
	std::string get(std::string key);
	int clear(std::string dbPath);
private:
	rocksdb::DB* kvdb;
	rocksdb::Options options;
};
