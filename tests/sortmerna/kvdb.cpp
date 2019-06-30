/* 
 * FILE: kvdb.cpp
 * Created: Jun 06, 2018 Wed
 */
#include <iostream>
#include <cassert>

#include "kvdb.hpp"
#include "options.hpp"

void kvdb_clear()
{
	std::string dbpath = "C:/a01_projects/clarity_genomics/data/kvdb";
	KeyValueDatabase kvdb(dbpath);
	int ret = kvdb.clear(dbpath);
	assert(ret == 0);
}