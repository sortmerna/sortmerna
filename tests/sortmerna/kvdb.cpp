/* 
 * FILE: kvdb.cpp
 * Created: Jun 06, 2018 Wed
 */
#include <iostream>
#include <cassert>

#include "kvdb.hpp"
#include "options.hpp"

int main(int argc, char** argv)
{
	Runopts opts(argc, argv, false);
	KeyValueDatabase kvdb(opts.kvdbPath);
	int ret = kvdb.clear(opts.kvdbPath);
	assert(ret == 0);
}