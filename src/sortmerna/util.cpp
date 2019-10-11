/* 
 * FILE: util.cpp
 * Created: Jun 10, 2018 Sun
 * @copyright 2016-19 Clarity Genomics BVBA
 */
#include <iostream> // cerr
#include <fstream>
#include <string>
#include <cstring>
#include <dirent.h>
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>

#include "common.hpp"

// forward
unsigned int check_dir(std::string dpath);
unsigned int list_dir(std::string dpath);
int clear_dir(std::string dpath);
bool dirExists(std::string dpath);
std::string get_user_home();

unsigned int check_dir(std::string dpath)
{
	unsigned int retval = 0;
	auto count = list_dir(dpath);
	if (count > 0) {
		std::cout << STAMP << "Directory " << dpath << " exists and is not empty" << std::endl;
		retval = 1;
	}
	return retval;
}

unsigned int list_dir(std::string dpath)
{
	unsigned int count = 0;
	if (!dirExists(dpath)) return count;

	DIR *pdir = opendir(dpath.data());
	struct dirent *next_file;
	std::string fpath;

	if (pdir == NULL)
	{
		std::cerr << STAMP << "Failed to open path " << dpath << std::endl;
		exit(1);
	}

	while ((next_file = readdir(pdir)) != NULL)
	{
		if (0 == strcmp(next_file->d_name, ".") || 0 == strcmp(next_file->d_name, ".."))
			continue; // skip '.' and '..'
		++count;
	}
	closedir(pdir);

	std::cout << STAMP << "Directory " << dpath << " has " << count << " files" << std::endl;
 	return count;
}

int clear_dir(std::string dpath)
{
	if (!dirExists(dpath)) return 0;

	DIR *pdir = opendir(dpath.data());
	struct dirent *next_file;
	std::string fpath;

	std::cout << "Cleaning directory: " << dpath << std::endl;

	if (pdir == NULL)
	{
		std::cerr << STAMP << "Failed to open path " << dpath << std::endl;
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

/* 
 * @param dpath String - Directory Path
 */
bool dirExists(std::string dpath)
{
	bool exists = false;
	struct stat info;

	if (stat(dpath.data(), &info) != 0)
		std::cout << __FUNCTION__ << ": Path does not exist: " << dpath << std::endl;
	else if (info.st_mode & S_IFDIR)  // S_ISDIR() doesn't exist on my windows 
	{
		std::cout << __FUNCTION__ << ": Path is a directory: " << dpath << std::endl;
		exists = true;
	}
	else
		std::cout << __FUNCTION__ << ": Path is Not a directory: " << dpath << std::endl;
	return exists;
}

std::string get_user_home()
{
	std::string homedir;
#if defined(_WIN32)
	homedir.append(getenv("USERPROFILE"));
	std::replace(homedir.begin(), homedir.end(), '\\', '/');
	//homedir.append(getenv("HOMEDRIVE"));
	//homedir.append(getenv("HOMEPATH"));
#else
	homedir.append(getenv("HOME"));
#endif
	return homedir;
}