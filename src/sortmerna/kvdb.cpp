/* 
 * FILE: kvdb.cpp
 * Created: Jun 05, 2018
 * @copyright 2016-19 Clarity Genomics BVBA
 */
#include "kvdb.hpp"

KeyValueDatabase::KeyValueDatabase(std::string const &kvdbPath) {
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