/* 
 * FILE: util.cpp
 * Created: Jun 10, 2018 Sun
 */
#include <iostream> // cerr
#include <string>
#include <dirent.h>

int clear_dir(std::string dpath)
{
	DIR *pdir = opendir(dpath.data());
	struct dirent *next_file;
	std::string fpath;

	std::cout << "Cleaning directory: " << dpath << std::endl;

	if (pdir == NULL)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " Failed to open path " << dpath << std::endl;
		return 1;
	}

	while ((next_file = readdir(pdir)) != NULL)
	{
		fpath = dpath + "/" + next_file->d_name;
		if (0 == strcmp(next_file->d_name, ".") || 0 == strcmp(next_file->d_name, ".."))
			continue; // skip '.' and '..'
		int st = remove(fpath.c_str());
		std::cout << "File: " << next_file->d_name << " deleted: " << st << std::endl;
	}
	closedir(pdir);
	return 0;
} /// ~clear_dir