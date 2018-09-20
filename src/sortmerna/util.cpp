/* 
 * FILE: util.cpp
 * Created: Jun 10, 2018 Sun
 */
#include <iostream> // cerr
#include <string>
#include <cstring>
#include <dirent.h>

#include <sys/types.h>
#include <sys/stat.h>

// forward
int clear_dir(std::string dpath);
bool dirExists(std::string dpath);

int clear_dir(std::string dpath)
{
	if (!dirExists(dpath)) return 0;

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