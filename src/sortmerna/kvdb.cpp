/* 
 * FILE: kvdb.cpp
 * Created: Jun 05, 2018
 */
#include <dirent.h>
#include <iostream> // cerr

#include "kvdb.hpp"

/* 
 * Remove database files from the given location
 */
int KeyValueDatabase::clear(std::string dbpath)
{
	DIR *pdir = opendir(dbpath.data());
	struct dirent *next_file;
	std::string fpath;

	if (pdir == NULL)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " Failed to open path " << dbpath << std::endl;
		return 1;
	}

	while ((next_file = readdir(pdir)) != NULL)
	{
		fpath = dbpath + "/" + next_file->d_name;
		if (0 == strcmp(next_file->d_name, ".") || 0 == strcmp(next_file->d_name, ".."))
			continue; // skip '.' and '..'
		remove(fpath.c_str());
	}
	closedir(pdir);
	return 0;
} // ~KeyValueDatabase::clear