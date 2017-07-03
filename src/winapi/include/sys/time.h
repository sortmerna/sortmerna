//#pragma once
/**
 * FILE: sys/time.h
 * Created: May 29, 2017
 */
#ifndef _WIN_TIMES_H
#if defined(_WIN32)

#include <windows.h>

struct timezone
{
	int  tz_minuteswest; /* minutes W of Greenwich */
	int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz);

#endif /*_WIN32*/
#endif /*_WIN_TIMES_H*/