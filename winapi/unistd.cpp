#if defined(_WIN32)
#include <Windows.h>
#include "unistd.h"
#include <cstdio>

long sysconf(int opt)
{
	long ret = 0;
	switch (opt)
	{
	case _SC_PAGE_SIZE:
		ret = 4096;
		break;
	case _SC_PHYS_PAGES:
		MEMORYSTATUSEX status;
		status.dwLength = sizeof(status);
		GlobalMemoryStatusEx(&status);
		DWORDLONG total_memsize = (size_t)status.ullTotalPhys;
		ret = total_memsize / 4096;
		break;
//	default:
//		printf("don't know case %s", opt);
//		ret = -1;
//		break;
	}

	return ret;
}
#endif /* _WIN32 */