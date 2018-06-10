/* 
 * FILE: kvdb.cpp
 * Created: Jun 06, 2018 Wed
 */
#include <iostream>
#include <cassert>

#include "kvdb.hpp"
#include "options.hpp"

void test_kvdb_clear(KeyValueDatabase & kvdb, std::string & dbpath)
{
	int ret = kvdb.clear(dbpath);
	assert(ret == 0);
}

int main(int argc, char** argv)
{
	std::string dbpath = "C:/a01_projects/clarity_genomics/data/kvdb";
	//Runopts opts(argc, argv, false);
	KeyValueDatabase kvdb(dbpath);

	test_kvdb_clear(kvdb, dbpath);

	return 0;
}