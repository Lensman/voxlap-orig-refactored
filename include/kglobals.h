#pragma once
#include "voxlap5.h"
#include "porthacks.h"
#define VOXSIZ VSID*VSID*64
#define PI 3.141592653589793
typedef struct { long f, p, x, y; } tiletype;
#define SCISDIST 1.0
#define MAX_PATH 260

extern __int64 flashbrival;
extern long cputype;
extern long fpuasm[2];
extern long frameplace, bytesperline;
extern char nullst; //nullst always NULL string

// File io
extern long backtag, backedup, bacx0, bacy0, bacx1, bacy1;
extern char *sxlbuf;
extern long sxlparspos, sxlparslen;

extern long skypic, nskypic, skybpl, skyysiz, skycurlng, skycurdir;
extern float skylngmul;
extern point2d *skylng;
// Used by kscreen
extern long ylookup[];
extern point3d gipos, gistr, gihei, gifor;
extern float gihx, gihy, gihz, gposxfrac[2], gposyfrac[2], grd;
extern long lastx[], uurendmem[], *uurend;

extern __int64 foglut[], fogcol;
extern long ofogdist;

extern long totclk;
extern long numsprites;
extern char sxlfilnam[];
extern char vxlfilnam[];
extern char skyfilnam[];
extern long sxlmallocsiz, sxlind[];


	//Displayed text message:
extern char message[256];
extern long messagetimeout;

// Located in external assembly object
#ifdef __cplusplus
EXTERN_C {
#endif
EXTERN_C char *sptr[]; 
EXTERN_C long gmipnum; 
EXTERN_C long skyoff, skyxsiz, *skylat;
EXTERN_C long zbufoff;
#ifdef __cplusplus
}
#endif

