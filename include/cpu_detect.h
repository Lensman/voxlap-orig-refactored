#pragma once

#if defined(CPU_DETECT_C) && defined(__cplusplus)
	extern "C" {
#endif

extern inline long testflag (long);
extern inline void cpuid (long, long *);
extern inline long getcputype ();
extern inline void fpuinit (long);
extern inline void code_rwx_unlock ( void *, void *);
extern inline int check_fpu_rdtsc();

#if defined(CPU_DETECT_C) && defined(__cplusplus)
}
#endif