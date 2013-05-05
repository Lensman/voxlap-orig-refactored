#pragma once

/** Generic, high precision timer functions */
//================== Fast & accurate TIMER FUNCTIONS begins ==================
#include "porthacks.h"

#ifdef _WIN32
    #include <windows.h>
#else
    #include <sys/time.h>
#endif

static __int64 pertimbase, rdtimbase, nextimstep;
static double perfrq, klockmul, klockadd, klocksub;

static MUST_INLINE uint64_t rdtsc64(void)
{
    #if defined(__GNUC__)
	uint64_t q;
	__asm__ __volatile__ ("rdtsc" : "=A" (q) : : "cc");
	return q;
	#elif defined(_MSC_VER)
    _asm rdtsc
    #endif
}

#if 0
extern void initklock ();
extern void readklock (double *tim);
#if defined(__WATCOMC__)
__int64 rdtsc64 ();
#pragma aux rdtsc64 = "rdtsc" value [edx eax] modify nomemory parm nomemory;

#elif defined(_MSC_VER)

void initklock ()
{
	__int64 q;
	QueryPerformanceFrequency((LARGE_INTEGER *)&q);
	perfrq = (double)q;
	rdtimbase = rdtsc64();
	QueryPerformanceCounter((LARGE_INTEGER *)&pertimbase);
	nextimstep = 4194304; klockmul = 0.000000001; klockadd = 0.0;
}

void readklock (double *tim)
{
	__int64 q = rdtsc64()-rdtimbase;
	if (q > nextimstep)
	{
		__int64 p;
		double d;
		QueryPerformanceCounter((LARGE_INTEGER *)&p);
		d = klockmul; klockmul = ((double)(p-pertimbase))/(((double)q)*perfrq);
		klockadd += (d-klockmul)*((double)q);
		do { nextimstep <<= 1; } while (q > nextimstep);
	}
	(*tim) = ((double)q)*klockmul + klockadd;
}
#elif defined(__GNUC__)
timeval start_time;
timeval end_time;

void initklock(){
    gettimeofday(&start_time, NULL);
}

void readklock (double *tim)
{
    double time_micro_sec;
    gettimeofday(&end_time, NULL);
    // start time and end time in millisec
    double start_ms = (start_time.tv_sec * 1000000.0) + start_time.tv_usec;
    double end_ms = (end_time.tv_sec * 1000000.0) + end_time.tv_usec;

    (*tim) = (double) end_ms - start_ms;
}

#else
#error Fatal : No rdtsc64 function defined.
#endif

#endif //0
//=================== Fast & accurate TIMER FUNCTIONS ends ===================

