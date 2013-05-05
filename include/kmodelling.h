#pragma once
#include "voxlap5.h"

extern void insslab (long *, const long, long);
extern void delslab (long *, const long, long);

extern void setsphere (const lpoint3d *, long, const long &);
extern void setcube (const long &, const long &, const long &, const long &);

extern void expandbit256 (const void *, void *);
extern void expandstack (const long &, const long &, long *);
extern void expandbitstack (const long &, const long &, int64_t *);

extern void scum (const long &, const long &, const long &, const long &, long *);
extern long *scum2 (const long &, const long &);

extern void scum2finish ();
extern void scumfinish (); 

extern void setcylinder (lpoint3d *, lpoint3d *, long, long, long);
extern void setrect (lpoint3d *, lpoint3d *, long);
extern void setspans (vspans *, long, lpoint3d *, long);
extern void setheightmap (const unsigned char *, long, long, long, long, long, long, long);
