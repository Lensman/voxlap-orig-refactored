#pragma once
#include "voxlap5.h"

	//Falling voxels functions:
extern void checkfloatinbox (long, long, long, long, long, long);
extern void startfalls ();
extern void dofall (long);
extern long meltfall (vx5sprite *, const long, const long);

extern long meltsphere (vx5sprite *, lpoint3d *, long);
extern long meltspans (vx5sprite *, vspans *, long, lpoint3d *);

extern void finishfalls ();

