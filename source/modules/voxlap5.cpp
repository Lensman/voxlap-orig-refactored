/** Notes
    zbuffer is updated between grouscan and opticast
    pass pixel buffer to grouscan to draw zbuffer in GPU mem
*/
#include <math.h>
#include <stdio.h>

#define VOXLAP5
#include "../../include/voxlap5.h"
#include "../../include/kplib.h"
#include "../../include/ksnippits.h"
#include "../../include/kfonts.h"
#include "../../include/kmodelling.h"
#include "../../include/kfalling.h"

//VOXLAP engine by Ken Silverman (http://advsys.net/ken)

#define USE_INTRINSICS
#define PREC (1024*PREC_TABLE)
#define CMPPREC (256*PREC_TABLE)
#define FPREC (256*PREC_TABLE)
#define USEV5ASM 1
#define SCISDIST 1.0
#define GOLDRAT 0.3819660112501052 //Golden Ratio: 1 - 1/((sqrt(5)+1)/2)
#define ESTNORMRAD 2 //Specially optimized for 2: DON'T CHANGE unless testing!

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <conio.h>
#include <dos.h>
#define MAX_PATH 260
#endif
#include <stdlib.h>

extern char keystatus[256];
extern void readkeyboard ();
extern void breath ();
extern long startdirectdraw (long *, long *, long *, long *);
extern void stopdirectdraw ();
extern void nextpage ();
extern void evilquit (const char *);
extern char *stripdir (char *filnam);

float gihx, gihy, gihz, gposxfrac[2], gposyfrac[2], grd; // was static
point3d gipos, gistr, gihei, gifor;// was static
long lastx[max(MAXYDIM,VSID)];

long ylookup[MAXYDIM+1];
long *vbuf = 0, *vbit = 0, vbiti;


#ifdef __cplusplus
extern "C" {
#endif
__ALIGN(16) char *sptr[(VSID*VSID*4)/3]; // Span pointers
#ifdef __cplusplus
}
#endif
/** voxel buffer format

	\warning : loaddta uses last 2MB of vbuf; vbuf:[VOXSIZ>>2], vbit:[VOXSIZ>>7]
	\warning : loadpng uses last 4MB of vbuf; vbuf:[VOXSIZ>>2], vbit:[VOXSIZ>>7]

    ╔════════════════════╤═════════╤════════╤════════╤════════╗
    ║       vbuf format: │   0:    │   1:   │   2:   │   3:   ║
    ╠════════════════════╪═════════╪════════╪════════╪════════╣
    ║      First header: │ nextptr │   z1   │   z1c  │  dummy ║
    ║           Color 1: │    b    │    g   │    r   │ intens ║
    ║           Color 2: │    b    │    g   │    r   │ intens ║
    ║             ...    │    b    │    g   │    r   │ intens ║
    ║           Color n: │    b    │    g   │    r   │ intens ║
    ║ Additional header: │ nextptr │   z1   │   z1c  │   z0   ║
    ╚════════════════════╧═════════╧════════╧════════╧════════╝
      nextptr: add this # <<2 to index to get to next header (0 if none)
           z1: z floor (top of floor color list)
          z1c: z bottom of floor color list MINUS 1! - needed to calculate
                 slab size with slng() and used as a separator for fcol/ccol
           z0: z ceiling (bottom of ceiling color list)
*/

	//Memory management variables:
#define MAXCSIZ 1028
char tbuf[MAXCSIZ];
long tbuf2[MAXZDIM*3];
long templongbuf[MAXZDIM];

extern long cputype; //bit25=1: SSE, bits30&31=1,1:3DNow!+
char nullst = 0; //nullst always NULL string

#pragma pack(push,1)
	//Rendering variables:
#if (USEZBUFFER == 0)
typedef struct { long col; } castdat;
#else
typedef struct { long col, dist; } castdat;
#endif
__ALIGN(16) typedef struct { castdat *i0, *i1; long z0, z1, cx0, cy0, cx1, cy1; } cftype;
//typedef struct { unsigned short x, y; } uspoint2d;
//typedef struct { long x, y; } lpoint2d;
//typedef struct { float x, y; } point2d;
#pragma pack(pop)

#if (USEV5ASM == 1)
#ifndef __cplusplus
	extern void *cfasm;
	extern castdat skycast;
#else
	extern "C" void *cfasm;
	extern "C" castdat skycast;
#endif
	#define cf ((cftype *)&cfasm)
#else
	cftype cf[256];
#endif

	//Screen related variables:
static long xres, yres;
long bytesperline, frameplace, xres4;
extern long ylookup[];

static lpoint3d glipos;

static point3d gixs, giys, gizs, giadd;

static long gposz, giforzsgn, gstartz0, gstartz1, gixyi[2];
static char *gstartv;

long backtag, backedup = -1, bacx0, bacy0, bacx1, bacy1;
char *bacsptr[262144];

	//Flash variables
#define LOGFLASHVANG 9
static lpoint2d gfc[(1<<LOGFLASHVANG)*8];
static long gfclookup[8] = {4,7,2,5,0,3,6,1}, flashcnt = 0;
int64_t flashbrival;

	//Norm flash variables
#define GSIZ 512  //NOTE: GSIZ should be 1<<x, and must be <= 65536
__ALIGN(16) long bbuf[GSIZ][GSIZ>>5], p2c[32], p2m[32];      //bbuf: 2.0K // was static
static __ALIGN(16) uspoint2d ffx[((GSIZ>>1)+2)*(GSIZ>>1)], *ffxptr; // ffx:16.5K
static long xbsox = -17, xbsoy, xbsof;
static __ALIGN(16) int64_t xbsbuf[25*5+1]; //need few bits before&after for protection

	//Look tables for expandbitstack256:
extern __ALIGN(16) long xbsceil[], xbsflor[]; // was static now in kmodelling



	//Opticast global variables:
	//radar: 320x200 requires  419560*2 bytes (area * 6.56*2)
	//radar: 400x300 requires  751836*2 bytes (area * 6.27*2)
	//radar: 640x480 requires 1917568*2 bytes (area * 6.24*2)
long __ALIGN(16) *radar = 0, *radarmem = 0;

#if (USEZBUFFER == 1)
static long *zbuffermem = 0, zbuffersiz = 0;
#endif

static castdat __ALIGN(16) *angstart[MAXXDIM*4], *gscanptr;
#define CMPRECIPSIZ (MAXXDIM)+32
static float cmprecip[CMPRECIPSIZ], wx0, wy0, wx1, wy1;
static long iwx0, iwy0, iwx1, iwy1;
static point3d gcorn[4];
		 point3d ginor[4]; //Should be static, but... necessary for stupid pingball hack :/
static long  uurendmem[MAXXDIM*2+(LVSID-2)], *uurend; // was static

void mat0(point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *);
void mat1(point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *);
void mat2(point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *, point3d *);

	//Parallaxing sky variables:
long skypic = 0, nskypic = 0, skybpl, skyysiz, skycurlng, skycurdir;
float skylngmul;
point2d *skylng = 0;

extern double logint[];      //was static, now in kmodelling
extern float tempfloatbuf[]; //was static, now in kmodelling
extern long factr[][2];      //was static, now in kmodelling

#ifdef __cplusplus
extern "C" {
#endif

	//Parallaxing sky variables (accessed by assembly code)
long skyoff = 0, skyxsiz, *skylat = 0;

__int64 gi, gcsub[8] =
{
	0xff00ff00ff00ff,0xff00ff00ff00ff,0xff00ff00ff00ff,0xff00ff00ff00ff,
	0xff00ff00ff00ff,0xff00ff00ff00ff,0xff00ff00ff00ff,0xff00ff00ff00ff
};
long gylookup[VSID+40], gmipnum = 0; //256+4+128+4+64+4+...
long gpz[2], gdz[2], gxmip, gxmax, gixy[2], gpixy;
static long gmaxscandist;
tiletype gdd;
//long reax, rebx, recx, redx, resi, redi, rebp, resp, remm[16];
void v5_asm_dep_unlock();
void grouscanasm (long);
long zbufoff;
#if (USEZBUFFER == 1)

#endif
#ifdef __cplusplus
}
#endif
#define gi0 (((long *)&gi)[0])
#define gi1 (((long *)&gi)[1])


	//if (a < 0) return(0); else if (a > b) return(b); else return(a);
inline long lbound0 (const long a, const long b) //b MUST be >= 0
{
	return ((unsigned long)a <= b) 
		? a 
		: ((~(a>>31))&b);
}

	//if (a < b) return(b); else if (a > c) return(c); else return(a);
inline long lbound (const long a, const long b, long c) //c MUST be >= b
{
	c -= b;
	return (unsigned long)(a-b) <= c 
		? a 
		: ((((b-a)>>31)&c) + b);
}

#define LSINSIZ 8 //Must be >= 2!
static point2d usintab[(1<<LSINSIZ)+(1<<(LSINSIZ-2))];
static void ucossininit ()
{
	long i, j;
	double a, ai, s, si, m;

	j = 0; usintab[0].y = 0.0;
	i = (1<<LSINSIZ)-1;
	ai = PI*(-2)/((float)(1<<LSINSIZ)); a = ((float)(-i))*ai;
	ai *= .5; m = sin(ai)*2; s = sin(a); si = cos(a+ai)*m; m = -m*m;
	for(;i>=0;i--)
	{
		usintab[i].y = s; s += si; si += s*m; //MUCH faster than next line :)
		//usintab[i].y = sin(i*PI*2/((float)(1<<LSINSIZ)));
		usintab[i].x = (usintab[j].y-usintab[i].y)/((float)(1<<(32-LSINSIZ)));
		j = i;
	}
	for(i=(1<<(LSINSIZ-2))-1;i>=0;i--) usintab[i+(1<<LSINSIZ)] = usintab[i];
}

	//Calculates cos & sin of 32-bit unsigned long angle in ~15 clock cycles
	//  Accuracy is approximately +/-.0001
static inline void ucossin (unsigned long a, float *cosin)
{
	float f = ((float)(a&((1<<(32-LSINSIZ))-1))); a >>= (32-LSINSIZ);
	cosin[0] = usintab[a+(1<<(LSINSIZ-2))].x*f+usintab[a+(1<<(LSINSIZ-2))].y;
	cosin[1] = usintab[a                 ].x*f+usintab[a                 ].y;
}



static long gkrand = 0;
long colorjit (long i, long jitamount)
{
	gkrand = (gkrand*27584621)+1;
	return((gkrand&jitamount)^i);
}



	//Note: ebx = 512 is no change
	//If PENTIUM III:1.Replace punpcklwd&punpckldq with: pshufw mm1, mm1, 0
	//               2.Use pmulhuw, shift by 8 & mul by 256
	//  :(  Can't mix with floating point
//#pragma aux colormul =
//   "movd mm0, eax"
//   "pxor mm1, mm1"
//   "punpcklbw mm0, mm1"
//   "psllw mm0, 7"
//   "movd mm1, ebx"
//   "punpcklwd mm1, mm1"
//   "punpckldq mm1, mm1"
//   "pmulhw mm0, mm1"
//   "packsswb mm0, mm0"
//   "movd eax, mm0"
//   parm [eax][ebx]
//   modify exact [eax]
//   value [eax]

long colormul (long i, long mulup8)
{
	long r, g, b;

	r = ((((i>>16)&255)*mulup8)>>8); if (r > 255) r = 255;
	g = ((((i>>8 )&255)*mulup8)>>8); if (g > 255) g = 255;
	b = ((((i    )&255)*mulup8)>>8); if (b > 255) b = 255;
	return((i&0xff000000)+(r<<16)+(g<<8)+b);
}

long curcolfunc (lpoint3d * __restrict p) 
{ 
	return(vx5.curcol); 
}

long floorcolfunc (lpoint3d * __restrict p)
{
	char *v;
	for(v=sptr[p->y*VSID+p->x];(p->z>v[2]) && (v[0]);v+=v[0]*4);
	return(*(long *)&v[4]);
}

long jitcolfunc (lpoint3d * __restrict p) { return(colorjit(vx5.curcol,vx5.amount)); }

static const long manycolukup[64] =
{
	  0,  1,  2,  5, 10, 15, 21, 29, 37, 47, 57, 67, 79, 90,103,115,
	127,140,152,165,176,188,198,208,218,226,234,240,245,250,253,254,
	255,254,253,250,245,240,234,226,218,208,198,188,176,165,152,140,
	128,115,103, 90, 79, 67, 57, 47, 37, 29, 21, 15, 10,  5,  2,  1
};
long manycolfunc (lpoint3d * __restrict p)
{
	return((manycolukup[p->x&63]<<16)+(manycolukup[p->y&63]<<8)+manycolukup[p->z&63]+0x80000000);
}

long sphcolfunc (lpoint3d * __restrict p)
{
	long i;
	ftol(sin((p->x+p->y+p->z-vx5.cen)*vx5.daf)*-96,&i);
	return(((i+128)<<24)|(vx5.curcol&0xffffff));
}

#define WOODXSIZ 46
#define WOODYSIZ 24
#define WOODZSIZ 24
static float wx[256], wy[256], wz[256], vx[256], vy[256], vz[256];
long woodcolfunc (lpoint3d * __restrict p)
{
	float col, u, a, f, dx, dy, dz;
	long i, c, xof, yof, tx, ty, xoff;

	if (*(long *)&wx[0] == 0)
	{
		for(i=0;i<256;i++)
		{
			wx[i] = WOODXSIZ * ((float)rand()/32768.0f-.5f) * .5f;
			wy[i] = WOODXSIZ * ((float)rand()/32768.0f-.5f) * .5f;
			wz[i] = WOODXSIZ * ((float)rand()/32768.0f-.5f) * .5f;

				//UNIFORM spherical randomization (see spherand.c)
			dz = 1.0f-(float)rand()/32768.0f*.04f;
			a = (float)rand()/32768.0f*PI*2.0f; fcossin(a,&dx,&dy);
			f = sqrt(1.0f-dz*dz); dx *= f; dy *= f;
				//??z: rings,  ?z?: vertical,  z??: horizontal (nice)
			vx[i] = dz; vy[i] = fabs(dy); vz[i] = dx;
		}
	}

		//(tx&,ty&) = top-left corner of current panel
	ty = p->y - (p->y%WOODYSIZ);
	xoff = ((ty/WOODYSIZ)*(ty/WOODYSIZ)*51721 + (p->z/WOODZSIZ)*357) % WOODXSIZ;
	tx = ((p->x+xoff) - (p->x+xoff)%WOODXSIZ) - xoff;

	xof = p->x - (tx + (WOODXSIZ>>1));
	yof = p->y - (ty + (WOODYSIZ>>1));

	c = ((((tx*429 + 4695) ^ (ty*341 + 4355) ^ 13643) * 2797) & 255);
	dx = xof - wx[c];
	dy = yof - wy[c];
	dz = (p->z%WOODZSIZ) - wz[c];

		//u = distance to center of randomly oriented cylinder
	u = vx[c]*dx + vy[c]*dy + vz[c]*dz;
	u = sqrt(dx*dx + dy*dy + dz*dz - u*u);

		//ring randomness
	u += sin((float)xof*.12 + (float)yof*.15) * .5;
	u *= (sin(u)*.05 + 1);

		//Ring function: smooth saw-tooth wave
	col = sin(u*2)*24;
	col *= pow(1.f-vx[c],.3f);

		//Thin shaded borders
	if ((p->x-tx == 0) || (p->y-ty == 0)) col -= 5;
	if ((p->x-tx == WOODXSIZ-1) || (p->y-ty == WOODYSIZ-1)) col -= 3;

	//f = col+c*.12+72; i = ftolp3(&f);
	  ftol(col+c*.12f+72.0f,&i);

	return(colormul(vx5.curcol,i<<1));
}

long gxsizcache = 0, gysizcache = 0;
long pngcolfunc (lpoint3d * __restrict p)
{
	long x, y, z, u, v;
	float fx, fy, fz, rx, ry, rz;

	if (!vx5.pic) return(vx5.curcol);
	switch(vx5.picmode)
	{
		case 0:
			x = p->x-vx5.pico.x; y = p->y-vx5.pico.y; z = p->z-vx5.pico.z;
			u = (((x&vx5.picu.x) + (y&vx5.picu.y) + (z&vx5.picu.z))^vx5.xoru);
			v = (((x&vx5.picv.x) + (y&vx5.picv.y) + (z&vx5.picv.z))^vx5.xorv);
			break;
		case 1: case 2:
			fx = (float)p->x-vx5.fpico.x;
			fy = (float)p->y-vx5.fpico.y;
			fz = (float)p->z-vx5.fpico.z;
			rx = vx5.fpicu.x*fx + vx5.fpicu.y*fy + vx5.fpicu.z*fz;
			ry = vx5.fpicv.x*fx + vx5.fpicv.y*fy + vx5.fpicv.z*fz;
			rz = vx5.fpicw.x*fx + vx5.fpicw.y*fy + vx5.fpicw.z*fz;
			ftol(atan2(ry,rx)*vx5.xoru/(PI*2),&u);
			if (vx5.picmode == 1) ftol(rz,&v);
			else ftol((atan2(rz,sqrt(rx*rx+ry*ry))/PI+.5)*vx5.ysiz,&v);
			break;
		default: //case 3:
			fx = (float)p->x-vx5.fpico.x;
			fy = (float)p->y-vx5.fpico.y;
			fz = (float)p->z-vx5.fpico.z;
			ftol(vx5.fpicu.x*fx + vx5.fpicu.y*fy + vx5.fpicu.z*fz,&u);
			ftol(vx5.fpicv.x*fx + vx5.fpicv.y*fy + vx5.fpicv.z*fz,&v);
			break;
	}
	if ((unsigned long)(u-gxsizcache) >= (unsigned long)vx5.xsiz){
		if (u < 0) { gxsizcache = u-(u+1)%vx5.xsiz-vx5.xsiz+1; 
		} else { gxsizcache = u-(u%vx5.xsiz); }
	}
	if ((unsigned long)(v-gysizcache) >= (unsigned long)vx5.ysiz){
		if (v < 0) { gysizcache = v-(v+1)%vx5.ysiz-vx5.ysiz+1;
		} else { gysizcache = v-(v%vx5.ysiz);}
	}
	return((vx5.pic[(v-gysizcache)*(vx5.bpl>>2)+(u-gxsizcache)]&0xffffff)|0x80000000);
}

	//Special case for SETSEC & SETCEI bumpmapping (vx5.picmode == 3)
	//no safety checks, returns alpha as signed char in range: (-128 to 127)
long hpngcolfunc (point3d * p)
{
	long u, v;
	float fx, fy, fz;

	fx = p->x-vx5.fpico.x;
	fy = p->y-vx5.fpico.y;
	fz = p->z-vx5.fpico.z;
	ftol(vx5.fpicu.x*fx + vx5.fpicu.y*fy + vx5.fpicu.z*fz,&u);
	ftol(vx5.fpicv.x*fx + vx5.fpicv.y*fy + vx5.fpicv.z*fz,&v);

	if ((unsigned long)(u-gxsizcache) >= (unsigned long)vx5.xsiz)
		if (u < 0) gxsizcache = u-(u+1)%vx5.xsiz-vx5.xsiz+1; else gxsizcache = u-(u%vx5.xsiz);
	if ((unsigned long)(v-gysizcache) >= (unsigned long)vx5.ysiz)
		if (v < 0) gysizcache = v-(v+1)%vx5.ysiz-vx5.ysiz+1; else gysizcache = v-(v%vx5.ysiz);
	return(vx5.pic[(v-gysizcache)*(vx5.bpl>>2)+(u-gxsizcache)]>>24);
}

long slng (const char *s)
{
	const char *v;

	for(v=s;v[0];v+=v[0]*4);
	return((long)v-(long)s+(v[2]-v[1]+1)*4+4);
}

void voxdealloc (const char *v)
{
	long i, j;
	i = (((long)v-(long)vbuf)>>2); j = (slng(v)>>2)+i;
#if 0
	while (i < j) { vbit[i>>5] &= ~(1<<i); i++; }
#else
	if (!((j^i)&~31))
		vbit[i>>5] &= ~(p2m[j&31]^p2m[i&31]);
	else
	{
		vbit[i>>5] &=   p2m[i&31];  i >>= 5;
		vbit[j>>5] &= (~p2m[j&31]); j >>= 5;
		for(j--;j>i;j--) vbit[j] = 0;
	}
#endif
}

	//Note: danum MUST be a multiple of 4!
char *voxalloc (long danum)
{
	long i, badcnt, p0, p1, vend;

	badcnt = 0; danum >>= 2; vend = (VOXSIZ>>2)-danum;
	do
	{
		for(;vbiti<vend;vbiti+=danum)
		{
			if (vbit[vbiti>>5]&(1<<vbiti)) continue;
			for(p0=vbiti;(!(vbit[(p0-1)>>5]&(1<<(p0-1))));p0--);
			for(p1=p0+danum-1;p1>vbiti;p1--)
				if (vbit[p1>>5]&(1<<p1)) goto allocnothere;

			vbiti = p0+danum;
			for(i=p0;i<vbiti;i++) vbit[i>>5] |= (1<<i);
			return((char *)(&vbuf[p0]));
allocnothere:;
		}
		vbiti = 0; badcnt++;
	} while (badcnt < 2);
	evilquit("voxalloc: vbuf full"); return(0);
}

inline long isvoxelsolid (const long& x, const long& y, const long& z)
{
	if ((unsigned long)(x|y) >= VSID) return(0);
	char const * __restrict v = sptr[y*VSID+x];
	while (1)
	{
		if (z < v[1]) return(0);
		if (!v[0]) return(1);
		v += v[0]*4;
		if (z < v[3]) return(1);
	}
}

	//Returns 1 if any voxels in range (x,y,z0) to (x,y,z1-1) are solid, else 0
long anyvoxelsolid (const long& x, const long& y, const long& z0, const long& z1)
{
	if ((unsigned long)(x|y) >= VSID) return(0);
		//         v1.....v3   v1.....v3    v1.......................>
		//                z0.........z1
	char const * __restrict v;
	v = sptr[y*VSID+x];
	while (1)
	{
		if (z1 <= v[1]) return(0);
		if (!v[0]) return(1);
		v += v[0]*4;
		if (z0 < v[3]) return(1);
	}
}

	//Returns 1 if any voxels in range (x,y,z0) to (x,y,z1-1) are empty, else 0
long anyvoxelempty (const long& x, const long& y, const long& z0, const long& z1)
{
	//         v1.....v3   v1.....v3    v1.......................>
	//                z0.........z1
	if ((unsigned long)(x|y) >= VSID) return(1);
	char const *v;

	v = sptr[y*VSID+x];
	while (1)
	{
		if (z0 < v[1]) return(1);
		if (!v[0]) return(0);
		v += v[0]*4;
		if (z1 <= v[3]) return(0);
	}
}

	//Returns z of first solid voxel under (x,y,z). Returns z if in solid.
long getfloorz (const long& x, const long& y, const long& z)
{
	if ((unsigned long)(x|y) >= VSID) return(z);
	char const * v;

	v = sptr[y*VSID+x];
	while (1)
	{
		if (z <= v[1]) return(v[1]);
		if (!v[0]) break;
		v += v[0]*4;
		if (z < v[3]) break;
	}
	return(z);
}

	//Returns:
	//   0: air
	//   1: unexposed solid
	//else: address to color in vbuf (this can never be 0 or 1)
long getcube (const long& x, const long& y, const long& z)
{
	if ((unsigned long)(x|y) >= VSID) return(0);
	long ceilnum;
	register char const *v;

	v = sptr[y*VSID+x];
	while (1)
	{
		if (z <= v[2])
		{
			if (z < v[1]) return(0);
			return((long)&v[(z-v[1])*4+4]);
		}
		ceilnum = v[2]-v[1]-v[0]+2;

		if (!v[0]) return(1);
		v += v[0]*4;

		if (z < v[3])
		{
			if (z-v[3] < ceilnum) return(1);
			return((long)&v[(z-v[3])*4]);
		}
	}
}
void gline (long leng, float x0, float y0, float x1, float y1)
{
	uint64_t q;
	float f, f1, f2, vd0, vd1, vz0, vx1, vy1, vz1;
	long j;
	cftype *c;
#if (USEV5ASM == 0)
	long gx, ogx, gy, ixy, col, dax, day;
	cftype *c2, *ce;
	char *v;
#endif

	vd0 = x0*gistr.x + y0*gihei.x + gcorn[0].x;
	vd1 = x0*gistr.y + y0*gihei.y + gcorn[0].y;
	vz0 = x0*gistr.z + y0*gihei.z + gcorn[0].z;
	vx1 = x1*gistr.x + y1*gihei.x + gcorn[0].x;
	vy1 = x1*gistr.y + y1*gihei.y + gcorn[0].y;
	vz1 = x1*gistr.z + y1*gihei.z + gcorn[0].z;

	f = sqrt(vx1*vx1 + vy1*vy1);
	f1 = f / vx1;
	f2 = f / vy1;
	if (fabs(vx1) > fabs(vy1)) { 
		vd0 = vd0*f1; 
	} else { 
		vd0 = vd1*f2; }
	
	if (*(long *)&vd0 < 0) vd0 = 0; //vd0 MUST NOT be negative: bad for asm
	vd1 = f;
	ftol(fabs(f1)*PREC,&gdz[0]);
	ftol(fabs(f2)*PREC,&gdz[1]);

	gixy[0] = (((*(signed long *)&vx1)>>31)<<3)+4; //=sgn(vx1)*4
	gixy[1] = gixyi[(*(unsigned long *)&vy1)>>31]; //=sgn(vy1)*4*VSID

	if (gdz[0] <= 0) { 	
		ftol( gposxfrac[(*(unsigned long *)&vx1)>>31]*fabs(f1)*PREC,&gpz[0]); 
		if (gpz[0] <= 0) { gpz[0] = 0x7fffffff; } 
		gdz[0] = 0x7fffffff-gpz[0]; //Hack for divide overflow
	} else { 
		ftol( gposxfrac[(*(unsigned long *)&vx1)>>31]*(float)gdz[0],&gpz[0]); 
	}
	
	if (gdz[1] <= 0) { 
		ftol( gposyfrac[(*(unsigned long *)&vy1)>>31]*fabs(f2)*PREC,&gpz[1]); 
		if (gpz[1] <= 0) gpz[1] = 0x7fffffff; 
		gdz[1] = 0x7fffffff-gpz[1]; //Hack for divide overflow
	} else { 
		ftol( gposyfrac[(*(unsigned long *)&vy1)>>31]*(float)gdz[1], &gpz[1] ); 
	}

	c = &cf[128];
	c->i0 = gscanptr; c->i1 = &gscanptr[leng];
	c->z0 = gstartz0; c->z1 = gstartz1;
	if (giforzsgn < 0)
	{
		ftol((vd1-vd0)*cmprecip[leng],&gi0); ftol(vd0*CMPPREC,&c->cx0);
		ftol((vz1-vz0)*cmprecip[leng],&gi1); ftol(vz0*CMPPREC,&c->cy0);
	}
	else
	{
		ftol((vd0-vd1)*cmprecip[leng],&gi0); ftol(vd1*CMPPREC,&c->cx0);
		ftol((vz0-vz1)*cmprecip[leng],&gi1); ftol(vz1*CMPPREC,&c->cy0);
	}
	c->cx1 = leng*gi0 + c->cx0;
	c->cy1 = leng*gi1 + c->cy0;

	gxmax = gmaxscandist;

		//Hack for early-out case when looking up towards sky
#if 0  //DOESN'T WORK WITH LOWER MIPS!
	if (c->cy1 < 0)
		if (gposz > 0)
		{
			if (dmulrethigh(-gposz,c->cx1,c->cy1,gxmax) >= 0)
			{
				j = scale(-gposz,c->cx1,c->cy1)+PREC; //+PREC for good luck
				if ((unsigned long)j < (unsigned long)gxmax) gxmax = j;
			}
		} else gxmax = 0;
#endif

		//Clip borders safely (MUST use integers!) - don't wrap around
#if ((USEZBUFFER == 1) && (USEV5ASM != 0))
	skycast.dist = gxmax;
#endif
	
	if (gixy[0] < 0) { j = glipos.x; } else{ j = VSID-1-(glipos.x); }
	q = mul64(gdz[0],j); q += (uint64_t)gpz[0];
	
	if (q < (uint64_t)gxmax)
	{
		gxmax = (long)q;
#if ((USEZBUFFER == 1) && (USEV5ASM != 0))
		skycast.dist = 0x7fffffff;
#endif
	}
	if (gixy[1] < 0) { j = glipos.y; } else {j = VSID-1-(glipos.y);}
	q = mul64(gdz[1],j); q += (uint64_t)gpz[1];

	if (q < (uint64_t)gxmax)
	{
		gxmax = (long)q;
#if ((USEZBUFFER == 1) && (USEV5ASM != 0))
		skycast.dist = 0x7fffffff;
#endif
	}

	if (vx5.sideshademode)
	{
		gcsub[0] = gcsub[(((unsigned long)gixy[0])>>31)+4];
		gcsub[1] = gcsub[(((unsigned long)gixy[1])>>31)+6];
	}

#if (USEV5ASM == 1)
	if (nskypic)
	{
		if (skycurlng < 0)
		{
			ftol((atan2(vy1,vx1)+PI)*skylngmul-.5,&skycurlng);
			if ((unsigned long)skycurlng >= skyysiz)
				skycurlng = ((skyysiz-1)&(j>>31));
		}
		else if (skycurdir < 0)
		{
			j = skycurlng+1; if (j >= skyysiz) j = 0;
			while (skylng[j].x*vy1 > skylng[j].y*vx1)
				{ skycurlng = j++; if (j >= skyysiz) j = 0; }
		}
		else
		{
			while (skylng[skycurlng].x*vy1 < skylng[skycurlng].y*vx1)
				{ skycurlng--; if (skycurlng < 0) skycurlng = skyysiz-1; }
		}
		skyoff = skycurlng*skybpl + nskypic;
	}

	//resp = 0;
	grouscanasm((long)gstartv);
	//if (resp)
	//{
	//   static char tempbuf[2048], tempbuf2[256];
	//   sprintf(tempbuf,"eax:%08x\tmm0:%08x%08x\nebx:%08x\tmm1:%08x%08x\necx:%08x\tmm2:%08x%08x\nedx:%08x\tmm3:%08x%08x\nesi:%08x\tmm4:%08x%08x\nedi:%08x\tmm5:%08x%08x\nebp:%08x\tmm6:%08x%08x\nesp:%08x\tmm7:%08x%08x\n",
	//      reax,remm[ 1],remm[ 0], rebx,remm[ 3],remm[ 2],
	//      recx,remm[ 5],remm[ 4], redx,remm[ 7],remm[ 6],
	//      resi,remm[ 9],remm[ 8], redi,remm[11],remm[10],
	//      rebp,remm[13],remm[12], resp,remm[15],remm[14]);
	//
	//   for(j=0;j<3;j++)
	//   {
	//      sprintf(tempbuf2,"%d i0:%d i1:%d z0:%ld z1:%ld cx0:%08x cy0:%08x cx1:%08x cy1:%08x\n",
	//         j,(long)cf[j].i0-(long)gscanptr,(long)cf[j].i1-(long)gscanptr,cf[j].z0,cf[j].z1,cf[j].cx0,cf[j].cy0,cf[j].cx1,cf[j].cy1);
	//      strcat(tempbuf,tempbuf2);
	//   }
	//   evilquit(tempbuf);
	//}
#else
//------------------------------------------------------------------------
	ce = c; v = gstartv;
	j = (((unsigned long)(gpz[1]-gpz[0]))>>31);
	gx = gpz[j];
	ixy = gpixy;
	if (v == (char *)*(long *)gpixy) goto drawflor; goto drawceil;

	while (1)
	{

drawfwall:;
		if (v[1] != c->z1)
		{
			if (v[1] > c->z1) c->z1 = v[1];
			else { do
			{
				c->z1--; col = *(long *)&v[(c->z1-v[1])*4+4];
				while (dmulrethigh(gylookup[c->z1],c->cx1,c->cy1,ogx) < 0)
				{
					c->i1->col = col; c->i1--; if (c->i0 > c->i1) goto deletez;
					c->cx1 -= gi0; c->cy1 -= gi1;
				}
			} while (v[1] != c->z1); }
		}

		if (v == (char *)*(long *)ixy) goto drawflor;

//drawcwall:;
		if (v[3] != c->z0)
		{
			if (v[3] < c->z0) c->z0 = v[3];
			else { do
			{
				c->z0++; col = *(long *)&v[(c->z0-v[3])*4-4];
				while (dmulrethigh(gylookup[c->z0],c->cx0,c->cy0,ogx) >= 0)
				{
					c->i0->col = col; c->i0++; if (c->i0 > c->i1) goto deletez;
					c->cx0 += gi0; c->cy0 += gi1;
				}
			} while (v[3] != c->z0); }
		}

drawceil:;
		while (dmulrethigh(gylookup[c->z0],c->cx0,c->cy0,gx) >= 0)
		{
			c->i0->col = (*(long *)&v[-4]); c->i0++; if (c->i0 > c->i1) goto deletez;
			c->cx0 += gi0; c->cy0 += gi1;
		}

drawflor:;
		while (dmulrethigh(gylookup[c->z1],c->cx1,c->cy1,gx) < 0)
		{
			c->i1->col = *(long *)&v[4]; c->i1--; if (c->i0 > c->i1) goto deletez;
			c->cx1 -= gi0; c->cy1 -= gi1;
		}

afterdelete:;
		c--;
		if (c < &cf[128])
		{
			ixy += gixy[j];
			gpz[j] += gdz[j];
			j = (((unsigned long)(gpz[1]-gpz[0]))>>31);
			ogx = gx; gx = gpz[j];

			if (gx > gxmax) break;
			v = (char *)*(long *)ixy; c = ce;
		}
			//Find highest intersecting vbuf slab
		while (1)
		{
			if (!v[0]) goto drawfwall;
			if (dmulrethigh(gylookup[v[2]+1],c->cx0,c->cy0,ogx) >= 0) break;
			v += v[0]*4;
		}
			//If next slab ALSO intersects, split cf!
		gy = gylookup[v[v[0]*4+3]];
		if (dmulrethigh(gy,c->cx1,c->cy1,ogx) < 0)
		{
			col = (long)c->i1; dax = c->cx1; day = c->cy1;
			while (dmulrethigh(gylookup[v[2]+1],dax,day,ogx) < 0)
				{ col -= sizeof(castdat); dax -= gi0; day -= gi1; }
			ce++; if (ce >= &cf[192]) return; //Give it max=64 entries like ASM
			for(c2=ce;c2>c;c2--) c2[0] = c2[-1];
			c[1].i1 = (castdat *)col; c->i0 = ((castdat *)col)+1;
			c[1].cx1 = dax; c->cx0 = dax+gi0;
			c[1].cy1 = day; c->cy0 = day+gi1;
			c[1].z1 = c->z0 = v[v[0]*4+3];
			c++;
		}
	}
//------------------------------------------------------------------------

	for(c=ce;c>=&cf[128];c--)
		while (c->i0 <= c->i1) { c->i0->col = 0; c->i0++; }
	return;

deletez:;
	ce--; if (ce < &cf[128]) return;
	for(c2=c;c2<=ce;c2++) c2[0] = c2[1];
	goto afterdelete;
#endif
}

static inline void addusb (char *a, long b)
{
	(*a) += b; if ((*a) < b) (*a) = 255;
}
/** (cone diameter vs. % 3D angular area) or: (a vs. 2/(1-cos(a*.5*PI/180)))
		╔═════════════╤═══════════╤═══════════╤══════════╤═══════════╤═════════╗
		║  0: inf     │ 25: 84.37 │ 50: 21.35 │ 75: 9.68 │ 100: 5.60 │ 180: 2  ║
		║  5: 2101.33 │ 30: 58.70 │ 55: 17.70 │ 80: 8.55 │ 105: 5.11 │ 360: 1  ║
		║ 10:  525.58 │ 35: 43.21 │ 60: 14.93 │ 85: 7.61 │ 110: 4.69 ╠═════════╝
		║ 15:  233.78 │ 40: 33.16 │ 65: 12.77 │ 90: 6.83 │ 115: 4.32 ║
		║ 20:  131.65 │ 45: 26.27 │ 70: 11.06 │ 95: 6.17 │ 120: 4    ║
		╚═════════════╧═══════════╧═══════════╧══════════╧═══════════╝
 */
void setflash (float px, float py, float pz, long flashradius, long numang, long intens)
{
	uint64_t q;
	float vx, vy;
	long i, j, gx, ogx, ixy, col, angoff;
	long ipx, ipy, ipz, sz0, sz1;
	cftype *c, *c2, *ce;
	char *v, *vs;

	ipx = (long)px; ipy = (long)py; ipz = (long)pz;
	vx5.minx = ipx-flashradius; vx5.maxx = ipx+flashradius+1;
	vx5.miny = ipy-flashradius; vx5.maxy = ipy+flashradius+1;
	vx5.minz = ipz-flashradius; vx5.maxz = ipz+flashradius+1;

	if (flashradius > 2047) flashradius = 2047;
	flashradius *= FPREC;

	flashbrival = (intens<<24);
	angoff = ((0x52741630>>(flashcnt<<2))&15); flashcnt++;

	gposxfrac[1] = px - (float)(ipx); gposxfrac[0] = 1 - gposxfrac[1];
	gposyfrac[1] = py - (float)(ipy); gposyfrac[0] = 1 - gposyfrac[1];
	gpixy = (long)&sptr[ipy*VSID + ipx];
	ftol(pz*FPREC-.5f,&gposz);
	for(gylookup[0]=-gposz,i=1;i<260;i++) gylookup[i] = gylookup[i-1]+FPREC;

	vs = (char *)*(long *)gpixy;
	if (ipz >= vs[1])
	{
		do
		{
			if (!vs[0]) return;
			vs += vs[0]*4;
		} while (ipz >= vs[1]);
		if (ipz < vs[3]) return;
		sz0 = vs[3];
	} else sz0 = 0;
	sz1 = vs[1];

	for(i=0;i<numang;i++)
	{
		clearMMX();

		fcossin(((float)i+(float)angoff*.125f)*PI*2.0f/(float)numang,&vx,&vy);

		ftol(FPREC/fabs(vx),&gdz[0]);
		ftol(FPREC/fabs(vy),&gdz[1]);

		gixy[0] = (((*(signed long *)&vx)>>31) & (     -8)) +      4;
		gixy[1] = (((*(signed long *)&vy)>>31) & (VSID*-8)) + VSID*4;
		if (gdz[0] < 0) { gpz[0] = 0x7fffffff; gdz[0] = 0; } //Hack for divide overflow
		else ftol(gposxfrac[(*(unsigned long *)&vx)>>31]*(float)gdz[0],&gpz[0]);
		if (gdz[1] < 0) { gpz[1] = 0x7fffffff; gdz[1] = 0; } //Hack for divide overflow
		else ftol(gposyfrac[(*(unsigned long *)&vy)>>31]*(float)gdz[1],&gpz[1]);

		c = ce = &cf[128];
		v = vs; c->z0 = sz0; c->z1 = sz1;
			//Note!  These substitions are used in flashscan:
			//   c->i0 in flashscan is now: c->cx0
			//   c->i1 in flashscan is now: c->cx1
		c->cx0 = (((i+flashcnt+rand())&7)<<LOGFLASHVANG);
		c->cx1 = c->cx0+(1<<LOGFLASHVANG)-1;

		gxmax = flashradius;

			//Clip borders safely (MUST use integers!) - don't wrap around
		
		if (gixy[0] < 0) j = ipx; else j = VSID-1-ipx;
		q = mul64(gdz[0],j); q += (uint64_t)gpz[0];
		if (q < (uint64_t)gxmax) gxmax = (long)q;
		if (gixy[1] < 0) j = ipy; else j = VSID-1-ipy;
		q = mul64(gdz[1],j); q += (uint64_t)gpz[1];
		if (q < (uint64_t)gxmax) gxmax = (long)q;

	//------------------------------------------------------------------------
		j = (((unsigned long)(gpz[1]-gpz[0]))>>31);
		gx = gpz[j];
		ixy = gpixy;
		if (v == (char *)*(long *)gpixy) goto fdrawflor; goto fdrawceil;

		while (1)
		{

fdrawfwall:;
			if (v[1] != c->z1)
			{
				if (v[1] > c->z1) c->z1 = v[1];
				else { do
				{
					c->z1--; col = (long)&v[(c->z1-v[1])*4+4];
					while (dmulrethigh(gylookup[c->z1],gfc[c->cx1].x,gfc[c->cx1].y,ogx) < 0)
					{
						mmxcoloradd((long *)col); c->cx1--;
						if (c->cx0 > c->cx1) goto fdeletez;
					}
				} while (v[1] != c->z1); }
			}

			if (v == (char *)*(long *)ixy) goto fdrawflor;

//fdrawcwall:;
			if (v[3] != c->z0)
			{
				if (v[3] < c->z0) c->z0 = v[3];
				else { do
				{
					c->z0++; col = (long)&v[(c->z0-v[3])*4-4];
					while (dmulrethigh(gylookup[c->z0],gfc[c->cx0].x,gfc[c->cx0].y,ogx) >= 0)
					{
						mmxcoloradd((long *)col); c->cx0++;
						if (c->cx0 > c->cx1) goto fdeletez;
					}
				} while (v[3] != c->z0); }
			}

fdrawceil:;
			while (dmulrethigh(gylookup[c->z0],gfc[c->cx0].x,gfc[c->cx0].y,gx) >= 0)
			{
				mmxcoloradd((long *)&v[-4]); c->cx0++;
				if (c->cx0 > c->cx1) goto fdeletez;
			}

fdrawflor:;
			while (dmulrethigh(gylookup[c->z1],gfc[c->cx1].x,gfc[c->cx1].y,gx) < 0)
			{
				mmxcoloradd((long *)&v[4]); c->cx1--;
				if (c->cx0 > c->cx1) goto fdeletez;
			}

fafterdelete:;
			c--;
			if (c < &cf[128])
			{
				ixy += gixy[j];
				gpz[j] += gdz[j];
				j = (((unsigned long)(gpz[1]-gpz[0]))>>31);
				ogx = gx; gx = gpz[j];

				if (gx > gxmax) break;
				v = (char *)*(long *)ixy; c = ce;
			}
				//Find highest intersecting vbuf slab
			while (1)
			{
				if (!v[0]) goto fdrawfwall;
				if (dmulrethigh(gylookup[v[2]+1],gfc[c->cx0].x,gfc[c->cx0].y,ogx) >= 0) break;
				v += v[0]*4;
			}
				//If next slab ALSO intersects, split cf!
			if (dmulrethigh(gylookup[v[v[0]*4+3]],gfc[c->cx1].x,gfc[c->cx1].y,ogx) < 0)
			{
				col = c->cx1;
				while (dmulrethigh(gylookup[v[2]+1],gfc[col].x,gfc[col].y,ogx) < 0)
					col--;
				ce++; if (ce >= &cf[192]) break; //Give it max=64 entries like ASM
				for(c2=ce;c2>c;c2--) c2[0] = c2[-1];
				c[1].cx1 = col; c->cx0 = col+1;
				c[1].z1 = c->z0 = v[v[0]*4+3];
				c++;
			}
		}
fcontinue:;
	}

	clearMMX();
	updatebbox(vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,0);
	return;

fdeletez:;
	ce--; if (ce < &cf[128]) goto fcontinue;
	for(c2=c;c2<=ce;c2++) c2[0] = c2[1];
	goto fafterdelete;
}

#if (ESTNORMRAD == 2)
static const signed char bitnum[32] =
{
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5
};

static const long __ALIGN(16) bitsnum[32] =
{
	0        ,1-(2<<16),1-(1<<16),2-(3<<16),
	1        ,2-(2<<16),2-(1<<16),3-(3<<16),
	1+(1<<16),2-(1<<16),2        ,3-(2<<16),
	2+(1<<16),3-(1<<16),3        ,4-(2<<16),
	1+(2<<16),2        ,2+(1<<16),3-(1<<16),
	2+(2<<16),3        ,3+(1<<16),4-(1<<16),
	2+(3<<16),3+(1<<16),3+(2<<16),4,
	3+(3<<16),4+(1<<16),4+(2<<16),5
};
static float __ALIGN(16) fsqrecip[5860]; //75*75 + 15*15 + 3*3 = 5859 is max value (5*5*5 box)
#endif

void estnorm (const long x, const long y, long z, point3d *fp)
{
	__ALIGN(16) lpoint3d n;
	__ALIGN(16) long *lptr, xx, yy, zz, b[5], i, j, k;
	float f;

	n.x = 0; n.y = 0; n.z = 0;

#if (ESTNORMRAD == 2)
	if (labs(x-xbsox) + labs(y-xbsoy) > 1)
	{
			//x,y not close enough to cache: calls expandbitstack 25 times :(
		xbsox = x; xbsoy = y; xbsof = 24*5;
		lptr = (long *)(&xbsbuf[24*5+1]);
		for(yy=-2;yy<=2;yy++)
			for(xx=-2;xx<=2;xx++,lptr-=10)
				expandbitstack(x+xx,y+yy,(int64_t *)lptr);
	}
	else if (x != xbsox)
	{
			//shift xbsbuf cache left/right: calls expandbitstack 5 times :)
		if (x < xbsox) { xx = -2; xbsof -= 24*5; lptr = (long *)(&xbsbuf[xbsof+1]); }
					 else { xx = 2; lptr = (long *)(&xbsbuf[xbsof-5*5+1]); xbsof -= 1*5; }
		xbsox = x; if (xbsof < 0) xbsof += 25*5;
		for(yy=-2;yy<=2;yy++)
		{
			if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			expandbitstack(x+xx,y+yy,(int64_t *)lptr);
			lptr -= 5*10;
		}
	}
	else if (y != xbsoy)
	{
			//shift xbsbuf cache up/down: calls expandbitstack 5 times :)
		if (y < xbsoy) { yy = -2; xbsof -= 20*5; lptr = (long *)(&xbsbuf[xbsof+1]); }
					 else { yy = 2; lptr = (long *)(&xbsbuf[xbsof+1]); xbsof -= 5*5; }
		xbsoy = y; if (xbsof < 0) xbsof += 25*5;
		for(xx=-2;xx<=2;xx++)
		{
			if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			expandbitstack(x+xx,y+yy,(int64_t *)lptr);
			lptr -= 1*10;
		}
	}

	z -= 2;
	if ((z&31) <= 27) //2 <= (z&31) <= 29
		{ lptr = (long *)((long)(&xbsbuf[xbsof+1]) + ((z&~31)>>3)); z &= 31; }
	else
		{ lptr = (long *)((long)(&xbsbuf[xbsof+1]) + (z>>3)); z &= 7; }

	for(yy=-2;yy<=2;yy++)
	{
		if (lptr >= (long *)&xbsbuf[1+10*5])
		{
			b[0] = ((lptr[  0]>>z)&31); b[1] = ((lptr[-10]>>z)&31);
			b[2] = ((lptr[-20]>>z)&31); b[3] = ((lptr[-30]>>z)&31);
			b[4] = ((lptr[-40]>>z)&31); lptr -= 50;
		}
		else
		{
			b[0] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			b[1] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			b[2] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			b[3] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			b[4] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
		}

			//Make filter spherical
		//if (yy&1) { b[0] &= 0xe; b[4] &= 0xe; }
		//else if (yy) { b[0] &= 0x4; b[1] &= 0xe; b[3] &= 0xe; b[4] &= 0x4; }

		n.x += ((bitnum[b[4]]-bitnum[b[0]])<<1)+bitnum[b[3]]-bitnum[b[1]];
		j = bitsnum[b[0]]+bitsnum[b[1]]+bitsnum[b[2]]+bitsnum[b[3]]+bitsnum[b[4]];
		n.z += j; n.y += (*(signed short *)&j)*yy;
	}
	n.z >>= 16;
#else
	for(yy=-ESTNORMRAD;yy<=ESTNORMRAD;yy++)
		for(xx=-ESTNORMRAD;xx<=ESTNORMRAD;xx++)
			for(zz=-ESTNORMRAD;zz<=ESTNORMRAD;zz++)
				if (isvoxelsolid(x+xx,y+yy,z+zz))
					{ n.x += xx; n.y += yy; n.z += zz; }
#endif

#if 1
	f = fsqrecip[n.x*n.x + n.y*n.y + n.z*n.z];
	fp->x = ((float)n.x)*f; fp->y = ((float)n.y)*f; fp->z = ((float)n.z)*f;
#else
	//f = 1.0 / sqrt((double)(n.x*n.x + n.y*n.y + n.z*n.z));
		//fp->x = f*(float)n.x; fp->y = f*(float)n.y; fp->z = f*(float)n.z;
	zz = n.x*n.x + n.y*n.y + n.z*n.z;
	if (cputype&(1<<25))
	{

		#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM) //gcc inline asm
		__asm__ __volatile__
		(
			"cvtsi2ss	zz,%xmm0\n"
			"rsqrtss	%xmm0,%xmm0\n"
			//"movss	%xmm0, f\n"

				//fp->x = f*(float)n.x; fp->y = f*(float)n.y; fp->z = f*(float)n.z;\n"
			"cvtsi2ss	%xmm1,n.z\n"
			"shufps	0,%xmm0,%xmm0\n"
			"mov	fp,%eax\n"
			"movlhps	%xmm1,%xmm1\n"
			"cvtpi2ps	n,%xmm1\n"
			"mulps	%xmm1,%xmm0\n"
			"movlps	%xmm0,(%eax)\n"
			"movhlps	%xmm0,%xmm0\n"
			"movss	%xmm0,8(%eax)\n"
		);
		#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
		_asm
		{
			cvtsi2ss	xmm0, zz
			rsqrtss	xmm0, xmm0
			//movss	f, xmm0

				//fp->x = f*(float)n.x; fp->y = f*(float)n.y; fp->z = f*(float)n.z;
			cvtsi2ss	xmm1, n.z
			shufps	xmm0, xmm0, 0
			mov	eax, fp
			movlhps	xmm1, xmm1
			cvtpi2ps	xmm1, n
			mulps	xmm0, xmm1
			movlps	[eax], xmm0
			movhlps	xmm0, xmm0
			movss	[eax+8], xmm0
		}
		#else // C Default
		#error No C Default yet defined
		#endif
	}
	else
	{
		#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM) //gcc inline asm
		__asm__ __volatile__
		(
			"pi2fd	zz,%mm0\n"       //mm0:     0          zz
			"pfrsqrt	%mm0,%mm0\n" //mm0: 1/sqrt(zz) 1/sqrt(zz)
			"pi2fd	n.x,%mm1\n"      //mm1:     0         n.x
			"pi2fd	n.y,%mm2\n"      //mm2:     0         n.y
			"punpckldq	%mm2,%mm1\n" //mm1:    n.y        n.x
			"pi2fd	n.z,%mm2\n"      //mm2:     0         n.z
			"pfmul	%mm0,%mm1\n"     //mm1:n.y/sqrt(zz) n.x/sqrt(zz)
			"pfmul	%mm0,%mm2\n"     //mm2:     0       n.z/sqrt(zz)
			"mov	fp,%eax\n"
			"movq	%mm1,(%eax)\n"
			"movl	%mm2,8(%eax)\n"
			"femms\n"
		);
		#elif defined(_MSC_VER) && && defined(__i386__) && !defined(NOASM)
		_asm
		{
			pi2fd	mm0, zz      //mm0:      0          zz
			pfrsqrt	mm0, mm0     //mm0:  1/sqrt(zz) 1/sqrt(zz)
			pi2fd	mm1, n.x     //mm1:      0         n.x
			pi2fd	mm2, n.y     //mm2:      0         n.y
			punpckldq	mm1, mm2 //mm1:     n.y        n.x
			pi2fd	mm2, n.z     //mm2:      0         n.z
			pfmul	mm1, mm0     //mm1: n.y/sqrt(zz) n.x/sqrt(zz)
			pfmul	mm2, mm0     //mm2:      0       n.z/sqrt(zz)
			mov	eax, fp
			movq	[eax], mm1
			movd	[eax+8], mm2
			femms
		}
		#else // C Default
		#error No C Default yet defined
		#endif
	}
#endif
}

static long vspan (long x, long y0, long y1)
{
	long y, yy, *bbufx;

	y = (y0>>5); bbufx = &bbuf[x][0];
	if ((y1>>5) == y)
	{
		yy = bbufx[y]; bbufx[y] &= ~(p2m[y1&31]^p2m[y0&31]);
		return(bbufx[y] ^ yy);
	}

	if (!(bbufx[y]&(~p2m[y0&31])))
		if (!(bbufx[y1>>5]&p2m[y1&31]))
		{
			for(yy=(y1>>5)-1;yy>y;yy--)
				if (bbufx[yy]) goto vspan_skip;
			return(0);
		}
vspan_skip:;
	bbufx[y] &= p2m[y0&31];
	bbufx[y1>>5] &= (~p2m[y1&31]);
	for(yy=(y1>>5)-1;yy>y;yy--) bbufx[yy] = 0;
	return(1);
}

static long docube (long x, long y, long z)
{
	long x0, y0, x1, y1, g;

	ffxptr = &ffx[(z+1)*z-1];
	x0 = (long)ffxptr[x].x; x1 = (long)ffxptr[x].y;
	y0 = (long)ffxptr[y].x; y1 = (long)ffxptr[y].y;
	for(g=0;x0<x1;x0++) g |= vspan(x0,y0,y1);
	return(g);
}

void setnormflash (float px, float py, float pz, long flashradius, long intens)
{
	point3d fp;
	float f, fintens;
	long i, j, k, l, m, x, y, z, xx, yy, xi, yi, xe, ye, ipx, ipy, ipz;
	long ceilnum, sq;
	char *v;

	ipx = (long)px; ipy = (long)py; ipz = (long)pz;
	vx5.minx = ipx-flashradius+1; vx5.maxx = ipx+flashradius;
	vx5.miny = ipy-flashradius+1; vx5.maxy = ipy+flashradius;
	vx5.minz = ipz-flashradius+1; vx5.maxz = ipz+flashradius;

	if (isvoxelsolid(ipx,ipy,ipz)) return;

	fintens = intens;
	if (flashradius > (GSIZ>>1)) flashradius = (GSIZ>>1);

	xbsox = -17;

		/**      ╔═7═╗
		        11   8
		    ╔═11═╪═4═╪═8═╤═7═╗
		    3    0   1   2   3 
		    ╚═10═╪═5═╪═9═╧═6═╝
		         10  9
		         ╚═6═╝
        */
		//Do left&right faces of the cube
	for(j=1;j>=0;j--)
	{
		clearbuf((void *)bbuf,GSIZ*(GSIZ>>5),0xffffffff);
		for(y=1;y<flashradius;y++)
		{
			if (j) yy = ipy-y; else yy = ipy+y;
			for(xi=1,xe=y+1;xi>=-1;xi-=2,xe=-xe)
				for(x=(xi>>1);x!=xe;x+=xi)
				{
					xx = ipx+x;
					if ((unsigned long)(xx|yy) >= VSID) continue;
					v = sptr[yy*VSID+xx]; i = 0; sq = x*x+y*y;
					while (1)
					{
						for(z=v[1];z<=v[2];z++)
						{
							if (z-ipz < 0) { tbuf2[i] = z-ipz; tbuf2[i+1] = (long)&v[(z-v[1])*4+4]; i += 2; }
							else
							{
								//if (z-ipz < -y) continue; //TEMP HACK!!!
								if (z-ipz > y) goto normflash_exwhile1;
								if (!docube(x,z-ipz,y)) continue;
								estnorm(xx,yy,z,&fp); if (j) fp.y = -fp.y;
								f = fp.x*x + fp.y*y + fp.z*(z-ipz);
								if (*(long *)&f > 0) addusb(&v[(z-v[1])*4+7],f*fintens/((z-ipz)*(z-ipz)+sq));
							}
						}
						if (!v[0]) break;
						ceilnum = v[2]-v[1]-v[0]+2; v += v[0]*4;
						for(z=v[3]+ceilnum;z<v[3];z++)
						{
							if (z < ipz) { tbuf2[i] = z-ipz; tbuf2[i+1] = (long)&v[(z-v[3])*4]; i += 2; }
							else
							{
								//if (z-ipz < -y) continue; //TEMP HACK!!!
								if (z-ipz > y) goto normflash_exwhile1;
								if (!docube(x,z-ipz,y)) continue;
								estnorm(xx,yy,z,&fp); if (j) fp.y = -fp.y;
								f = fp.x*x + fp.y*y + fp.z*(z-ipz);
								if (*(long *)&f > 0) addusb(&v[(z-v[3])*4+3],f*fintens/((z-ipz)*(z-ipz)+sq));
							}
						}
					}
normflash_exwhile1:;
					while (i > 0)
					{
						i -= 2; if (tbuf2[i] < -y) break;
						if (!docube(x,tbuf2[i],y)) continue;
						estnorm(xx,yy,tbuf2[i]+ipz,&fp); if (j) fp.y = -fp.y;
						f = fp.x*x + fp.y*y + fp.z*tbuf2[i];
						if (*(long *)&f > 0) addusb(&((char *)tbuf2[i+1])[3],f*fintens/(tbuf2[i]*tbuf2[i]+sq));
					}
				}
		}
	}

		//Do up&down faces of the cube
	for(j=1;j>=0;j--)
	{
		clearbuf((void *)bbuf,GSIZ*(GSIZ>>5),0xffffffff);
		for(y=1;y<flashradius;y++)
		{
			if (j) xx = ipx-y; else xx = ipx+y;
			for(xi=1,xe=y+1;xi>=-1;xi-=2,xe=-xe)
				for(x=(xi>>1);x!=xe;x+=xi)
				{
					yy = ipy+x;
					if ((unsigned long)(xx|yy) >= VSID) continue;
					v = sptr[yy*VSID+xx]; i = 0; sq = x*x+y*y; m = x+xi-xe;
					while (1)
					{
						for(z=v[1];z<=v[2];z++)
						{
							if (z-ipz < 0) { tbuf2[i] = z-ipz; tbuf2[i+1] = (long)&v[(z-v[1])*4+4]; i += 2; }
							else
							{
								//if (z-ipz < -y) continue; //TEMP HACK!!!
								if (z-ipz > y) goto normflash_exwhile2;
								if ((!docube(x,z-ipz,y)) || (!m)) continue;
								estnorm(xx,yy,z,&fp); if (j) fp.x = -fp.x;
								f = fp.x*y + fp.y*x + fp.z*(z-ipz);
								if (*(long *)&f > 0) addusb(&v[(z-v[1])*4+7],f*fintens/((z-ipz)*(z-ipz)+sq));
							}
						}
						if (!v[0]) break;
						ceilnum = v[2]-v[1]-v[0]+2; v += v[0]*4;
						for(z=v[3]+ceilnum;z<v[3];z++)
						{
							if (z < ipz) { tbuf2[i] = z-ipz; tbuf2[i+1] = (long)&v[(z-v[3])*4]; i += 2; }
							else
							{
								//if (z-ipz < -y) continue; //TEMP HACK!!!
								if (z-ipz > y) goto normflash_exwhile2;
								if ((!docube(x,z-ipz,y)) || (!m)) continue;
								estnorm(xx,yy,z,&fp); if (j) fp.x = -fp.x;
								f = fp.x*y + fp.y*x + fp.z*(z-ipz);
								if (*(long *)&f > 0) addusb(&v[(z-v[3])*4+3],f*fintens/((z-ipz)*(z-ipz)+sq));
							}
						}
					}
normflash_exwhile2:;
					while (i > 0)
					{
						i -= 2; if (tbuf2[i] < -y) break;
						if ((!docube(x,tbuf2[i],y)) || (!m)) continue;
						estnorm(xx,yy,tbuf2[i]+ipz,&fp); if (j) fp.x = -fp.x;
						f = fp.x*y + fp.y*x + fp.z*tbuf2[i];
						if (*(long *)&f > 0) addusb(&((char *)tbuf2[i+1])[3],f*fintens/(tbuf2[i]*tbuf2[i]+sq));
					}
				}
		}
	}

		//Do the bottom face of the cube
	clearbuf((void *)bbuf,GSIZ*(GSIZ>>5),0xffffffff);
	for(yi=1,ye=flashradius+1;yi>=-1;yi-=2,ye=-ye)
		for(y=(yi>>1);y!=ye;y+=yi)
			for(xi=1,xe=flashradius+1;xi>=-1;xi-=2,xe=-xe)
				for(x=(xi>>1);x!=xe;x+=xi)
				{
					xx = ipx+x; yy = ipy+y;
					if ((unsigned long)(xx|yy) >= VSID) goto normflash_exwhile3;
					k = max(labs(x),labs(y));

					v = sptr[yy*VSID+xx]; sq = x*x+y*y;
					while (1)
					{
						for(z=v[1];z<=v[2];z++)
						{
							if (z-ipz < k) continue;
							if (z-ipz >= flashradius) goto normflash_exwhile3;
							if ((!docube(x,y,z-ipz)) || (z-ipz == k)) continue;
							estnorm(xx,yy,z,&fp);
							f = fp.x*x + fp.y*y + fp.z*(z-ipz);
							if (*(long *)&f > 0) addusb(&v[(z-v[1])*4+7],f*fintens/((z-ipz)*(z-ipz)+sq));
						}
						if (!v[0]) break;
						ceilnum = v[2]-v[1]-v[0]+2; v += v[0]*4;
						for(z=v[3]+ceilnum;z<v[3];z++)
						{
							if (z-ipz < k) continue;
							if (z-ipz >= flashradius) goto normflash_exwhile3;
							if ((!docube(x,y,z-ipz)) || (z-ipz <= k)) continue;
							estnorm(xx,yy,z,&fp);
							f = fp.x*x + fp.y*y + fp.z*(z-ipz);
							if (*(long *)&f > 0) addusb(&v[(z-v[3])*4+3],f*fintens/((z-ipz)*(z-ipz)+sq));
						}
					}
normflash_exwhile3:;
				}


		//Do the top face of the cube
	clearbuf((void *)bbuf,GSIZ*(GSIZ>>5),0xffffffff);
	for(yi=1,ye=flashradius+1;yi>=-1;yi-=2,ye=-ye)
		for(y=(yi>>1);y!=ye;y+=yi)
			for(xi=1,xe=flashradius+1;xi>=-1;xi-=2,xe=-xe)
				for(x=(xi>>1);x!=xe;x+=xi)
				{
					xx = ipx+x; yy = ipy+y;
					if ((unsigned long)(xx|yy) >= VSID) goto normflash_exwhile4;
					k = max(labs(x),labs(y)); m = ((x+xi != xe) && (y+yi != ye));

					v = sptr[yy*VSID+xx]; i = 0; sq = x*x+y*y;
					while (1)
					{
						for(z=v[1];z<=v[2];z++)
						{
							if (ipz-z >= flashradius) continue;
							if (ipz-z < k) goto normflash_exwhile4;
							tbuf2[i] = ipz-z; tbuf2[i+1] = (long)&v[(z-v[1])*4+4]; i += 2;
						}
						if (!v[0]) break;
						ceilnum = v[2]-v[1]-v[0]+2; v += v[0]*4;
						for(z=v[3]+ceilnum;z<v[3];z++)
						{
							if (ipz-z >= flashradius) continue;
							if (ipz-z < k) goto normflash_exwhile4;
							tbuf2[i] = ipz-z; tbuf2[i+1] = (long)&v[(z-v[3])*4]; i += 2;
						}
					}
normflash_exwhile4:;
					while (i > 0)
					{
						i -= 2;
						if ((!docube(x,y,tbuf2[i])) || (tbuf2[i] <= k)) continue;
						estnorm(xx,yy,ipz-tbuf2[i],&fp);
						f = fp.x*x + fp.y*y - fp.z*tbuf2[i];
						if (*(long *)&f > 0) addusb(&((char *)tbuf2[i+1])[3],f*fintens/(tbuf2[i]*tbuf2[i]+sq));
					}
				}
	
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,0 );
}

void hline (float x0, float y0, float x1, float y1, long *ix0, long *ix1)
{
	float dyx;

	dyx = (y1-y0) * grd; //grd = 1/(x1-x0)

		  if (y0 < wy0) ftol((wy0-y0)/dyx+x0,ix0);
	else if (y0 > wy1) ftol((wy1-y0)/dyx+x0,ix0);
	else ftol(x0,ix0);
		  if (y1 < wy0) ftol((wy0-y0)/dyx+x0,ix1);
	else if (y1 > wy1) ftol((wy1-y0)/dyx+x0,ix1);
	else ftol(x1,ix1);
	if ((*ix0) < iwx0) (*ix0) = iwx0;
	if ((*ix0) > iwx1) (*ix0) = iwx1; //(*ix1) = min(max(*ix1,wx0),wx1);
	gline(labs((*ix1)-(*ix0)),(float)(*ix0),((*ix0)-x1)*dyx + y1,
									  (float)(*ix1),((*ix1)-x1)*dyx + y1);
}

void vline (float x0, float y0, float x1, float y1, long *iy0, long *iy1)
{
	float dxy;

	dxy = (x1-x0) * grd; //grd = 1/(y1-y0)

		  if (x0 < wx0) ftol((wx0-x0)/dxy+y0,iy0);
	else if (x0 > wx1) ftol((wx1-x0)/dxy+y0,iy0);
	else ftol(y0,iy0);
		  if (x1 < wx0) ftol((wx0-x0)/dxy+y0,iy1);
	else if (x1 > wx1) ftol((wx1-x0)/dxy+y0,iy1);
	else ftol(y1,iy1);
	if ((*iy0) < iwy0) (*iy0) = iwy0;
	if ((*iy0) > iwy1) (*iy0) = iwy1;
	gline(labs((*iy1)-(*iy0)),((*iy0)-y1)*dxy + x1,(float)(*iy0),
									  ((*iy1)-y1)*dxy + x1,(float)(*iy1));
}

static float __ALIGN(16) optistrx, optistry, optiheix, optiheiy, optiaddx, optiaddy;

__int64 foglut[2048], fogcol;
long ofogdist = 700;

#ifdef _MSC_VER

#ifdef __cplusplus
extern "C" {
#endif
extern __ALIGN(16) void *opti4asm;
#define opti4 ((point4d *)&opti4asm)
#ifdef __cplusplus
}
#endif

void (*hrend)(long,long,long,long,long,long);
void (*vrend)(long,long,long,long,long);

inline float f_rsqrt( const float number )
{
        long i;
        float x2, y;
        const float threehalfs = 1.5F;
 
        x2 = number * 0.5F;
        y  = number;
        i  = * ( long * ) &y;
        i  = 0x5f3759df - ( i >> 1 );
        y  = * ( float * ) &i;
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
        // 2nd iteration, this can be removed if you don't need accuracy
      	y  = y * ( threehalfs - ( x2 * y * y ) );   
        return y;
}

inline void hrendz (long sx, long sy, long p1, long plc, long incr, long j)
{
	long p0, i; 
	castdat * c0;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;

	__m128 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;

	// this big chunk of scalar code looks like it needs to be vectorized
	// eg something like this;
	/*
		= movaps[ sx,       sy,       sx',      sy'      ] 
		* mul   [ optistrx, optiheix, optistry, optiheiy ] 
		+ addps [ optiaddx, 0,        optiaddy, 0        ] 
		+ addps [ sy,       0,        sy',      0        ] 
	*/

	xmm0 = _mm_cvtsi32_ss( xmm0, sx & 0xfffffffc );
	xmm4 = _mm_cvtsi32_ss( xmm4, sy );


	xmm1 = _mm_move_ss( xmm1 , xmm0 );
	xmm5 = _mm_move_ss( xmm5 , xmm4 );
	
	xmm0 = _mm_mul_ss( xmm0, *(__m128*)&optistrx );
	xmm1 = _mm_mul_ss( xmm1, *(__m128*)&optistry );
	xmm4 = _mm_mul_ss( xmm4, *(__m128*)&optiheix );
	xmm5 = _mm_mul_ss( xmm5, *(__m128*)&optiheiy );
	xmm0 = _mm_add_ss( xmm0, *(__m128*)&optiaddx );
	xmm1 = _mm_add_ss( xmm1, *(__m128*)&optiaddy );
	xmm0 = _mm_add_ss( xmm0, xmm4 );
	xmm1 = _mm_add_ss( xmm1, xmm5 );

	xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0 );         //[sx, sx, sx, sx]
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0 );	        //[sy, sy, sy, sy]
	xmm0 = _mm_add_ps( xmm0, *(__m128*)&opti4[0] ); //[sx, sx + optistrx, sx + optistrx*2, sx + optistrx*3 ] 
	xmm1 = _mm_add_ps( xmm1, *(__m128*)&opti4[1] ); //[sy, sy + optistry, sy + optistry*2, sy + optistry*3 ] 
	xmm2 = _mm_load_ps( (const float *)&opti4[2] ); //[optistrx*4.0f ... ]
	xmm3 = _mm_load_ps( (const float *)&opti4[3] ); //[optistry*4.0f ... ]
	xmm2 = _mm_add_ps( xmm2, xmm0 );                //[sx, sx, sx, sx]
	xmm3 = _mm_add_ps( xmm3, xmm1 );

	xmm0 = _mm_mul_ps( xmm0, xmm0 );
	xmm1 = _mm_mul_ps( xmm1, xmm1 );
	xmm2 = _mm_mul_ps( xmm2, xmm2 );
	xmm3 = _mm_mul_ps( xmm3, xmm3 );
	xmm0 = _mm_add_ps( xmm0, xmm1 );
	xmm1 = _mm_load_ps( (const float *)&opti4[4] ); //[ (optistrx*optistrx + optistry*optistry)*32.0f ... ]  
	xmm2 = _mm_add_ps( xmm2, xmm3 );
	xmm2 = _mm_sub_ps( xmm2, xmm0 );
	i = zbufoff;

	// --------------------------------------------------
	// Check for pixel alignment of 4, and render blocks of pixels 
	// based on that alignment.
	// --------------------------------------------------
	if ( ( *(long *)&p0 & 15 ) != 0 ) {
		
		if ( ( *(long*)&p0 & 8 ) != 0 ) {
			xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x4e  ); // rotate right by 2
		}
		if ( ( *(long*)&p0 & 4 ) != 0 ) {
			xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x39  ); // rotate right by 1
		}

		//;Do first 0-3 pixels to align unrolled loop of 4
		do {
			__m128i pixel = _mm_loadl_epi64( (const __m128i*)(&angstart[plc>>16][j]) );
			*(int *)p0 = _mm_cvtsi128_si32(pixel);

			xmm3 = _mm_rsqrt_ss(xmm0);
			pixel = _mm_shuffle_epi32(pixel, _MM_SHUFFLE(0,2,3,1) );
			xmm7 = _mm_cvtepi32_ps( pixel );
			xmm7 = _mm_mul_ss( xmm7, xmm3 );
			xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x39  );
			_mm_store_ss( (float*)(p0+i), xmm7 );
			plc += incr; p0 += 4;
			if (p0 == p1) {	return;	}
		} while ( ( *(long *)&p0 & 15 ) != 0 ); 
		
		#if USEZBUFFER ==1
		xmm0 = _mm_add_ps( xmm0, xmm2 );
		xmm2 = _mm_add_ps( xmm2, xmm1 );
		#endif
	}
	// --------------------------------------------------
	// Do 4 pixels at a time 
    // --------------------------------------------------
	if (p0+16 < p1) {
		do {

			// Could unwrap further on 64 bit systems with another 8 simd regs to play with
			// This could be converted to PINSW / PEXTRW method 

			// This block should end up as movq's
			__m128i col0i = _mm_loadl_epi64( (const __m128i*)&angstart[plc>>16][j] );
			__m128i col1i = _mm_loadl_epi64( (const __m128i*)&angstart[(plc+(incr))>>16][j] );
			__m128i col2i = _mm_loadl_epi64( (const __m128i*)&angstart[(plc+(incr*2))>>16][j] );						
			__m128i col3i = _mm_loadl_epi64( (const __m128i*)&angstart[(plc+(incr*3))>>16][j] );
			
			// Combine 4 128bit registers into 2 128bit registers
			__m128i comb1 = _mm_unpacklo_epi32( col0i, col1i );
			__m128i comb2 = _mm_unpacklo_epi32( col2i, col3i );
			__m128i pixel = _mm_unpacklo_epi32( comb1, comb2 );

			// Store color for 4 pixels
			_mm_stream_si128( (__m128i*)(p0), _mm_unpacklo_epi32( comb1, comb2 ) );		

			__m128 zpixel = _mm_cvtepi32_ps( _mm_unpackhi_epi32( comb1, comb2 ) );
			xmm3 = _mm_rsqrt_ps(xmm0); 
			xmm0 = _mm_add_ps(xmm0, xmm2);        
			xmm2 = _mm_add_ps(xmm2, xmm1);

			zpixel = _mm_mul_ps( zpixel, xmm3 );

			// Store 4 zbuffer pixels
			_mm_stream_ps( (float*)(p0+i), zpixel  );

			p0+=16;// Move raster
			plc+=incr*4;

		} while( p0+16 < p1 );
	}
    // --------------------------------------------------
	// Render remaining individual pixels
    // --------------------------------------------------
	for ( ;p0<p1;p0+=4,plc+=incr ) {
		c0 = &angstart[plc>>16][j];
		__m128i castdat = _mm_loadl_epi64( (const __m128i*)(c0) );
		*(int *)p0 = _mm_cvtsi128_si32(castdat);
		xmm3 = _mm_rsqrt_ss(xmm0); 
		castdat = _mm_shuffle_epi32(castdat, _MM_SHUFFLE(0,2,3,1) );
		xmm7 = _mm_cvtepi32_ps( castdat );         
		xmm7 = _mm_mul_ss( xmm7, xmm3 );              
		xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x39  );   // rotate right by 1
		_mm_store_ss( (float*)(p0+i), xmm7 );       
	} 
}
/*
void vrendz (long sx, long sy, long p1, long iplc, long iinc)
{
	float dirx, diry; long i, p0;
	castdat * c0;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;
	dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx;
	diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;
	i = zbufoff;
	while (p0 < p1)
	{
		c0 = &angstart[uurend[sx]>>16][iplc];
		*(long *)p0 = c0->col;
		*(float *)(p0+i) = (float)c0->dist*f_rsqrt(dirx*dirx+diry*diry);
		dirx += optistrx; diry += optistry; 
		uurend[sx] += uurend[sx+MAXXDIM]; 
		p0 += 4; iplc += iinc; sx++;
	}
}
*/
inline void vrendz (long sx, long sy, long p1, long iplc, long iinc)
{
	long p0, i; 
	castdat * c0;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;

	__m128 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;


	xmm0 = _mm_cvtsi32_ss( xmm0, sx & 0xfffffffc );
	xmm4 = _mm_cvtsi32_ss( xmm4, sy );


	xmm1 = _mm_move_ss( xmm1 , xmm0 );
	xmm5 = _mm_move_ss( xmm5 , xmm4 );
	
	xmm0 = _mm_mul_ss( xmm0, *(__m128*)&optistrx );
	xmm1 = _mm_mul_ss( xmm1, *(__m128*)&optistry );
	xmm4 = _mm_mul_ss( xmm4, *(__m128*)&optiheix );
	xmm5 = _mm_mul_ss( xmm5, *(__m128*)&optiheiy );
	xmm0 = _mm_add_ss( xmm0, *(__m128*)&optiaddx );
	xmm1 = _mm_add_ss( xmm1, *(__m128*)&optiaddy );
	xmm0 = _mm_add_ss( xmm0, xmm4 );
	xmm1 = _mm_add_ss( xmm1, xmm5 );

	xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0 );         //[sx, sx, sx, sx]
	xmm1 = _mm_shuffle_ps( xmm1, xmm1, 0 );	        //[sy, sy, sy, sy]
	xmm0 = _mm_add_ps( xmm0, *(__m128*)&opti4[0] ); //[sx, sx + optistrx, sx + optistrx*2, sx + optistrx*3 ] 
	xmm1 = _mm_add_ps( xmm1, *(__m128*)&opti4[1] ); //[sy, sy + optistry, sy + optistry*2, sy + optistry*3 ] 
	xmm2 = _mm_load_ps( (const float *)&opti4[2] ); //[optistrx*4.0f ... ]
	xmm3 = _mm_load_ps( (const float *)&opti4[3] ); //[optistry*4.0f ... ]
	xmm2 = _mm_add_ps( xmm2, xmm0 );                //[sx, sx, sx, sx]
	xmm3 = _mm_add_ps( xmm3, xmm1 );

	xmm0 = _mm_mul_ps( xmm0, xmm0 );
	xmm1 = _mm_mul_ps( xmm1, xmm1 );
	xmm2 = _mm_mul_ps( xmm2, xmm2 );
	xmm3 = _mm_mul_ps( xmm3, xmm3 );
	xmm0 = _mm_add_ps( xmm0, xmm1 );
	xmm1 = _mm_load_ps( (const float *)&opti4[4] ); //[ (optistrx*optistrx + optistry*optistry)*32.0f ... ]  
	xmm2 = _mm_add_ps( xmm2, xmm3 );
	xmm2 = _mm_sub_ps( xmm2, xmm0 );
	i = zbufoff;

	// --------------------------------------------------
	// Check for pixel alignment of 4, and render blocks of pixels 
	// based on that alignment.
	// --------------------------------------------------
	if ( ( *(long *)&p0 & 15 ) != 0 ) {
		
		if ( ( *(long*)&p0 & 8 ) != 0 ) {
			xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x4e  ); // rotate right by 2
		}
		if ( ( *(long*)&p0 & 4 ) != 0 ) {
			xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x39  ); // rotate right by 1
		}

		//;Do first 0-3 pixels to align unrolled loop of 4
		do {
			__m128i pixel = _mm_loadl_epi64( (const __m128i*)( &angstart[uurend[sx]>>16][iplc] ));
			*(int *)p0 = _mm_cvtsi128_si32(pixel);

			xmm3 = _mm_rsqrt_ss(xmm0);
			pixel = _mm_shuffle_epi32(pixel, _MM_SHUFFLE(0,2,3,1) );
			xmm7 = _mm_cvtepi32_ps( pixel );
			xmm7 = _mm_mul_ss( xmm7, xmm3 );
			xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x39  );
			_mm_store_ss( (float*)(p0+i), xmm7 );
			uurend[sx] += uurend[sx+MAXXDIM]; 
			iplc += iinc; p0 += 4; sx++;
			if (p0 == p1) {	return;	}
		} while ( ( *(long *)&p0 & 15 ) != 0 ); 

		xmm0 = _mm_add_ps( xmm0, xmm2 );
		xmm2 = _mm_add_ps( xmm2, xmm1 );
	}
	// --------------------------------------------------
	// Do 4 pixels at a time 
    // --------------------------------------------------

	while( p0+16 < p1 ) 
	{

		// This block should end up as movq's
		__m128i sx_offsets = _mm_load_si128( (const __m128i *)&uurend[sx] );
		/* sse4
		int idx0 = _mm_extract_epi32( sx_offsets, 0 );
		int idx1 = _mm_extract_epi32( sx_offsets, 2 );
		int idx2 = _mm_extract_epi32( sx_offsets, 4 );
		int idx3 = _mm_extract_epi32( sx_offsets, 6 );
		
		__m128i col0i = _mm_loadl_epi64( (const __m128i*)&angstart[idx0>>16][iplc] );
		__m128i col1i = _mm_loadl_epi64( (const __m128i*)&angstart[idx1>>16][iplc+iinc] );
		__m128i col2i = _mm_loadl_epi64( (const __m128i*)&angstart[idx2>>16][iplc+iinc*2] );						
		__m128i col3i = _mm_loadl_epi64( (const __m128i*)&angstart[idx3>>16][iplc+iinc*3] );
		*/
		
		__m128i col0i = _mm_loadl_epi64( (const __m128i*)&angstart[uurend[sx]>>16][iplc] );
		__m128i col1i = _mm_loadl_epi64( (const __m128i*)&angstart[uurend[sx+1]>>16][iplc+iinc] );
		__m128i col2i = _mm_loadl_epi64( (const __m128i*)&angstart[uurend[sx+2]>>16][iplc+iinc*2] );						
		__m128i col3i = _mm_loadl_epi64( (const __m128i*)&angstart[uurend[sx+3]>>16][iplc+iinc*3] );
		

		// Combine 4 128bit registers into 2 128bit registers
		__m128i comb1 = _mm_unpacklo_epi32( col0i, col1i );
		__m128i comb2 = _mm_unpacklo_epi32( col2i, col3i );
		__m128i pixel = _mm_unpacklo_epi32( comb1, comb2 );

		// Store color for 4 pixels
		_mm_stream_si128( (__m128i*)(p0), _mm_unpacklo_epi32( comb1, comb2 ) );		

		__m128 zpixel = _mm_cvtepi32_ps( _mm_unpackhi_epi32( comb1, comb2 ) );
		xmm3 = _mm_rsqrt_ps(xmm0); 
		xmm0 = _mm_add_ps(xmm0, xmm2);        
		xmm2 = _mm_add_ps(xmm2, xmm1);

		zpixel = _mm_mul_ps( zpixel, xmm3 );

		// Store 4 zbuffer pixels
		_mm_stream_ps( (float*)(p0+i), zpixel  );
		__m128i sx_offsets_to_add = _mm_loadu_si128( (const __m128i *)&uurend[sx+MAXXDIM] );
		sx_offsets = _mm_add_epi32(sx_offsets, sx_offsets_to_add );
		_mm_store_si128 ( (__m128i *)&uurend[sx], sx_offsets);
		iplc += iinc*4; p0 += 16; sx+=4;

	} 

    // --------------------------------------------------
	// Render remaining individual pixels
    // --------------------------------------------------
	while ( p0!=p1 ) {
		c0 = &angstart[uurend[sx]>>16][iplc];
		
		__m128i castdat = _mm_loadl_epi64( (const __m128i*)(c0) );
		*(int *)p0 = _mm_cvtsi128_si32(castdat);
		xmm3 = _mm_rsqrt_ss(xmm0); 
		castdat = _mm_shuffle_epi32(castdat, _MM_SHUFFLE(0,2,3,1) );
		xmm7 = _mm_cvtepi32_ps( castdat );         
		xmm7 = _mm_mul_ss( xmm7, xmm3 );              
		xmm0 = _mm_shuffle_ps( xmm0, xmm0, 0x39  );   // rotate right by 1
		_mm_store_ss( (float*)(p0+i), xmm7 );    
		uurend[sx] += uurend[sx+MAXXDIM];   
		p0+=4; iplc+=iinc; sx++;
	} 
}


#if (USEZBUFFER != 1)
void hrendnoz (long sx, long sy, long p1, long plc, long incr, long j)
{
	sy = ylookup[sy]+frameplace; p1 = sy+(p1<<2); sy += (sx<<2);
	do
	{
		*(long *)sy = angstart[plc>>16][j].col;
		plc += incr; sy += 4;
	} while (sy != p1);
}

void vrendnoz (long sx, long sy, long p1, long iplc, long iinc)
{
	sy = ylookup[sy]+(sx<<2)+frameplace;	
	for(;sx<p1;sx++)
	{
		*(long *)sy = angstart[uurend[sx]>>16][iplc].col;
		uurend[sx] += uurend[sx+MAXXDIM]; sy += 4; iplc += iinc;
	}
}

#else

#if 0
	//Example C code


	//Example C code
void hrendzfog (long sx, long sy, long p1, long plc, long incr, long j)
{
	long p0, i, k, l; float dirx, diry;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;
	dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx;
	diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;
	i = zbufoff;
	do
	{
		k = angstart[plc>>16][j].col;
		l = angstart[plc>>16][j].dist;
		l = (foglut[l>>20]&32767);
		*(long *)p0 = ((((( vx5.fogcol     &255)-( k     &255))*l)>>15)    ) +
						  ((((((vx5.fogcol>> 8)&255)-((k>> 8)&255))*l)>>15)<< 8) +
						  ((((((vx5.fogcol>>16)&255)-((k>>16)&255))*l)>>15)<<16)+k;
		*(float *)(p0+i) = (float)angstart[plc>>16][j].dist/sqrt(dirx*dirx+diry*diry);
		dirx += optistrx; diry += optistry; plc += incr; p0 += 4;
	} while (p0 != p1);
}

	//Example C code
void vrendzfog (long sx, long sy, long p1, long iplc, long iinc)
{
	float dirx, diry; long i, k, l, p0;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;
	dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx;
	diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;
	i = zbufoff;
	while (p0 < p1)
	{
		k = angstart[uurend[sx]>>16][iplc].col;
		l = angstart[uurend[sx]>>16][iplc].dist;
		l = (foglut[l>>20]&32767);
		*(long *)p0 = ((((( vx5.fogcol     &255)-( k     &255))*l)>>15)    ) +
						  ((((((vx5.fogcol>> 8)&255)-((k>> 8)&255))*l)>>15)<< 8) +
						  ((((((vx5.fogcol>>16)&255)-((k>>16)&255))*l)>>15)<<16)+k;
		*(float *)(p0+i) = (float)angstart[uurend[sx]>>16][iplc].dist/sqrt(dirx*dirx+diry*diry);
		dirx += optistrx; diry += optistry; uurend[sx] += uurend[sx+MAXXDIM]; p0 += 4; iplc += iinc; sx++;
	}
}

#endif


void hrendzfog (long sx, long sy, long p1, long plc, long incr, long j)
{
	long p0, i, k, l; float dirx, diry;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;
	dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx;
	diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;

	i = zbufoff;
	do
	{
		k = angstart[plc>>16][j].col;
		l = angstart[plc>>16][j].dist;
		l = (foglut[l>>20]&32767);
		*(long *)p0 = ((((( vx5.fogcol     &255)-( k     &255))*l)>>15)    ) +
						  ((((((vx5.fogcol>> 8)&255)-((k>> 8)&255))*l)>>15)<< 8) +
						  ((((((vx5.fogcol>>16)&255)-((k>>16)&255))*l)>>15)<<16)+k;
		*(float *)(p0+i) = (float)angstart[plc>>16][j].dist/sqrt(dirx*dirx+diry*diry);
		dirx += optistrx; diry += optistry; plc += incr; p0 += 4;
	} while (p0 != p1);
}

	//Example C code
void vrendzfog (long sx, long sy, long p1, long iplc, long iinc)
{
	float dirx, diry; long i, k, l, p0;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;
	dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx;
	diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;
	i = zbufoff;
	while (p0 < p1)
	{
		k = angstart[uurend[sx]>>16][iplc].col;
		l = angstart[uurend[sx]>>16][iplc].dist;
		l = (foglut[l>>20]&32767);
		*(long *)p0 = ((((( vx5.fogcol     &255)-( k     &255))*l)>>15)    ) +
						  ((((((vx5.fogcol>> 8)&255)-((k>> 8)&255))*l)>>15)<< 8) +
						  ((((((vx5.fogcol>>16)&255)-((k>>16)&255))*l)>>15)<<16)+k;
		*(float *)(p0+i) = (float)angstart[uurend[sx]>>16][iplc].dist/sqrt(dirx*dirx+diry*diry);
		dirx += optistrx; diry += optistry; uurend[sx] += uurend[sx+MAXXDIM]; p0 += 4; iplc += iinc; sx++;
	}
}
/** Similar to a 3DDDA */
void hrendzsse (long sx, long sy, long p1, long plc, long incr, long j)
{
	_asm
	{
		push esi
		push edi
		prefetchnta optiheix	
beghasm_p3:
		prefetchnta frameplace
		prefetchnta opti4asm[0*16]
		mov eax, sx
		mov ecx, sy
		mov esi, p1

		mov edx, ylookup[ecx*4]		; " edx = ylookup[sy]"
		add edx, frameplace         ; " edx += frameplace"
		
		lea edi, [edx+eax*4]	    ; " edi = ylookup[sy]+(sx<<2)+frameplace; "
		lea esi, [edx+esi*4]        ; " esi = ylookup[sy]+(p1<<2)+frameplace; "

		and eax, 0xfffffffc         ; " sx & 0xfffffffc"
		cvtsi2ss xmm0, eax          ; " [ 0, sx ]         ==  = movaps[ sx,       sy,       sx',      sy'      ] "
		cvtsi2ss xmm4, ecx          ; " [ 0, sy ]         ==  * mul   [ optistrx, optiheix, optistry, optiheiy ] "
		movss xmm1, xmm0            ; " [ 0, sx']         ==  + addps [ optiaddx, 0,        optiaddy, 0        ] "
		movss xmm5, xmm4            ; " [ 0, sy']         ==  + addps [ sy,       0,        sy',      0        ] "
		mulss xmm0, optistrx        ; " [ 0, sx ] * optistrx "
		mulss xmm1, optistry        ; " [ 0, sx'] * optistry "
		mulss xmm4, optiheix        ; " [ 0, sy ] * optiheix "
		mulss xmm5, optiheiy        ; " [ 0, sy'] * optiheiy "
		addss xmm0, optiaddx        ; " [ 0, sx ] + optiaddx "
		addss xmm1, optiaddy        ; " [ 0, sx'] + optiaddy "
		addss xmm0, xmm4            ; " [ 0, sx ] + [ 0, sy ]  == dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx; "
		addss xmm1, xmm5            ; " [ 0, sx'] + [ 0, sy']  == diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;"

		mov ecx, zbufoff
		mov edx, j
		movd mm6, plc
		movd mm7, incr

		shufps xmm0, xmm0, 0        ; " [ 0,0,0,sx ] = [ sx,sx,sx,sx ] "
		shufps xmm1, xmm1, 0        ; " [ 0,0,0,sy ] = [ sy,sy,sy,sy ] "

		addps xmm0, opti4asm[0*16]  ; " [ sx,sx,sx,sx ] + [ 0, optistrx, optistrx*2, optistrx*3 ] "
		addps xmm1, opti4asm[1*16]  ; " [ sy,sy,sy,sy ] + [ 0, optistry, optistry*2, optistry*3 ] " 
		movaps xmm2, opti4asm[2*16] ; " [ optistrx*4.0f, optistrx*4.0f, optistrx*4.0f, optistrx*4.0f ] "
		movaps xmm3, opti4asm[3*16] ; " [ optistry*4.0f, optistry*4.0f, optistry*4.0f, optistry*4.0f ] "
		
		;" Overview of algorithm "
		;" xmm0 =  xmm0      ^2 +  xmm1      ^2        (p)"
		;" xmm2 = (xmm0+xmm2)^2 + (xmm1+xmm3)^2 - xmm0 (v)"
		;" xmm1 = ...                                  (a)"
		;" This block converts inner loop...              "
		;" from: 1 / sqrt(x*x + y*y), x += xi, y += yi;   "
		;"   to: 1 / sqrt(p), p += v, v += a;             "

		addps xmm2, xmm0            ; " sxT = [ sx'',sx'',sx'',sx'' ] + [ optistrx*4.0f, optistrx*4.0f, optistrx*4.0f, optistrx*4.0f ] "
		addps xmm3, xmm1            ; " syT = [ sy'',sy'',sy'',sy'' ] + [ optistry*4.0f, optistry*4.0f, optistry*4.0f, optistry*4.0f ] "
		mulps xmm0, xmm0            ; " sx'' = [ sx'',sx'',sx'',sx'' ] * [ sx'',sx'',sx'',sx'' ] "
		mulps xmm1, xmm1            ; " sy'' = [ sy'',sy'',sy'',sy'' ] * [ sy'',sy'',sy'',sy'' ] "
		mulps xmm2, xmm2            ; " [ sxT ] * [ sxT ] "
		mulps xmm3, xmm3            ; " [ syT ] * [ syT ] "
		addps xmm0, xmm1            ; " sxsy = [ sy''... ] + [sx''...] "
		movaps xmm1, opti4asm[4*16] ; " (optistrx*optistrx + optistry*optistry)*32.0f  "
		addps xmm2, xmm3            ; " [ sxsyT ... ]  = [ sxT...] + [syT...]"
		subps xmm2, xmm0            ; " [ sxsyT ... ] -= [ sxsy ... ] "

			;Do first 0-3 pixels to align unrolled loop of 4
		test edi, 15
		jz short skip1ha

		test edi, 8
		jz short skipshufa
		shufps xmm0, xmm0, 0x4e ;rotate right by 2
skipshufa:
		test edi, 4
		jz short skipshufb
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
skipshufb:

beg1ha:
		pextrw eax, mm6, 1
		paddd mm6, mm7
		mov eax, angstart[eax*4]
		movd mm0, [eax+edx*8]
		movd [edi], mm0
		cvtsi2ss xmm7, [eax+edx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7                             ;zbufoff
		add edi, 4
		cmp edi, esi
		jz short endh
		test edi, 15
		jnz short beg1ha

		addps xmm0, xmm2
		addps xmm2, xmm1
skip1ha:
		lea eax, [edi+16]      ;these 3 lines re-ordered
		cmp eax, esi
		ja short skip4h

		movq mm0, mm6          ;mm0: 0,plc
		paddd mm0, mm7         ;mm0: 0,plc+inc
		punpckldq mm7, mm7     ;mm7: inc,inc
		paddd mm7, mm7         ;mm7: inc+inc,inc+inc
		punpckldq mm6, mm0     ;mm6: plc+inc,plc
		sub esi, 16

		 ;eax: temp   ³ mm0:  z0 argb0   argb1 argb0 ³ xmm0: plc3 plc2 plc1 plc0
		 ;ebx:  -     ³ mm1:  z1 argb1               ³ xmm1: acc3 acc2 acc1 acc0
		 ;ecx:zbufoff ³ mm2:  z2 argb2   argb3 argb2 ³ xmm2: inc3 inc2 inc1 inc0
		 ;edx:  j     ³ mm3:  z3 argb3               ³ xmm3:  r3   r2   r1   r0
		 ;esi:  -     ³ mm4:              z1    z0   ³ xmm4:            z3   z2
		 ;edi:scroff  ³ mm5:              z3    z2   ³ xmm5:
		 ;ebp:  -     ³ mm6: plc1 plc0               ³ xmm6:
beg4h:   ;esp:  -     ³ mm7: inc1 inc0               ³ xmm7:  z3   z2   z1   z0
		; "Lots of cache misses here"
		pextrw eax, mm6, 1           ; "plc1 Extract the 1st word from mm6 and into eax."
		mov eax, angstart[eax*4]
		movq mm0, [eax+edx*8]        ; "plc0 uurend[sx]>>16"
		pextrw eax, mm6, 3
		mov eax, angstart[eax*4]
		movq mm1, [eax+edx*8]
		paddd mm6, mm7               ; inc
		pextrw eax, mm6, 1           ; plc1
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]        ; z2
		pextrw eax, mm6, 3           ; plc0
		mov eax, angstart[eax*4]     ;
		movq mm3, [eax+edx*8]        ; z3
		paddd mm6, mm7               ; inc

		movq mm4, mm0                ; "save copy for zbuffer"
		movq mm5, mm2 
		punpckldq mm0, mm1
		punpckldq mm2, mm3
		movntq [edi], mm0        ; "Store"
		movntq [edi+8], mm2      ; "Store"

		punpckhdq mm4, mm1       ; "Interleave high doublewords of mm and mm/m64 into mm."
		punpckhdq mm5, mm3       ;
		cvtpi2ps xmm7, mm4       ; "Convert two signed doubleword integers to two SP-FP values"
		
		rsqrtps xmm3, xmm0       ; "packed approximations of the reciprocals"
		cvtpi2ps xmm4, mm5

		addps xmm0, xmm2         ; "Add packed single-precision"
		addps xmm2, xmm1

		movlhps xmm7, xmm4       ; "Move Low to High Packed Single-FP."
		mulps xmm7, xmm3

		movntps [edi+ecx], xmm7  ; "Move Aligned Four Packed Single-FP Non Temporal : zbufoff"


		add edi, 16
		cmp edi, esi
		
		jbe short beg4h
		add esi, 16
		cmp edi, esi
		jae endh

		psrad mm7, 1             ; "Restore mm7 from incr*2 to just incr for single loop"
skip4h:
beg1h:
		pextrw eax, mm6, 1
		paddd mm6, mm7
		mov eax, angstart[eax*4]
		movd mm0, [eax+edx*8]
		movd [edi], mm0
		cvtsi2ss xmm7, [eax+edx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39       ; rotate right by 1
		movss [edi+ecx], xmm7         ; store value in zbuff
		add edi, 4
        cmp edi, esi                  ; p0 < p1
        jb short beg1h
endh:   pop edi
		pop esi
	}
}

void hrendzfogsse (long sx, long sy, long p1, long plc, long incr, long j)
{
	static int64_t mm7bak;
	_asm
	{
		push esi
		push edi
beghasm_p3:
		mov eax, sx
		mov ecx, sy
		mov esi, p1
		mov edx, ylookup[ecx*4]
		add edx, frameplace
		lea edi, [edx+eax*4]
		lea esi, [edx+esi*4]

		and eax, 0xfffffffc
		cvtsi2ss xmm0, eax
		cvtsi2ss xmm4, ecx
		movss xmm1, xmm0
		movss xmm5, xmm4
		mulss xmm0, optistrx
		mulss xmm1, optistry
		mulss xmm4, optiheix
		mulss xmm5, optiheiy
		addss xmm0, optiaddx
		addss xmm1, optiaddy
		addss xmm0, xmm4
		addss xmm1, xmm5

		mov ecx, zbufoff
		mov edx, j
		movd mm6, plc
		movd mm7, incr

		shufps xmm0, xmm0, 0
		shufps xmm1, xmm1, 0
		movaps xmm2, opti4asm[2*16]
		movaps xmm3, opti4asm[3*16]
		addps xmm0, opti4asm[0*16]
		addps xmm1, opti4asm[1*16]
			;xmm0 =  xmm0      ^2 +  xmm1      ^2        (p)
			;xmm2 = (xmm0+xmm2)^2 + (xmm1+xmm3)^2 - xmm0 (v)
			;xmm1 = ...                                  (a)
		addps xmm2, xmm0  ;This block converts inner loop...
		addps xmm3, xmm1  ;from: 1 / sqrt(x*x + y*y), x += xi, y += yi;
		mulps xmm0, xmm0  ;  to: 1 / sqrt(p), p += v, v += a;
		mulps xmm1, xmm1
		mulps xmm2, xmm2
		mulps xmm3, xmm3
		addps xmm0, xmm1
		movaps xmm1, opti4asm[4*16]
		addps xmm2, xmm3
		subps xmm2, xmm0

			;Do first 0-3 pixels to align unrolled loop of 4
		test edi, 15
		jz short skip1ha

		test edi, 8
		jz short skipshufa
		shufps xmm0, xmm0, 0x4e ;rotate right by 2
skipshufa:
		test edi, 4
		jz short skipshufb
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
skipshufb:

beg1ha:
		pextrw eax, mm6, 1
		paddd mm6, mm7
		mov eax, angstart[eax*4]

			;Z
		cvtsi2ss xmm7, [eax+edx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7

			;Col
		punpcklbw mm0, [eax+edx*8]
		psrlw mm0, 8
		movq mm1, fogcol
		psubw mm1, mm0
		paddw mm1, mm1
		mov eax, [eax+edx*8+4]
		shr eax, 16+4
		pmulhw mm1, foglut[eax*8]
		paddw mm0, mm1
		packuswb mm0, mm1
		movd [edi], mm0

		add edi, 4
		cmp edi, esi
		jz short endh
		test edi, 15
		jnz short beg1ha

		addps xmm0, xmm2
		addps xmm2, xmm1
skip1ha:
		lea eax, [edi+16]      ;these 3 lines re-ordered
		cmp eax, esi
		ja short skip4h

		movq mm0, mm6          ;mm0: 0,plc
		paddd mm0, mm7         ;mm0: 0,plc+inc
		punpckldq mm7, mm7     ;mm7: inc,inc
		punpckldq mm6, mm0     ;mm6: plc+inc,plc
		paddd mm7, mm7         ;mm7: inc+inc,inc+inc

		sub esi, 16

		 ;eax: temp   ³ mm0:  z0 argb0   argb1 argb0 ³ xmm0: plc3 plc2 plc1 plc0
		 ;ebx:  -     ³ mm1:  z1 argb1               ³ xmm1: acc3 acc2 acc1 acc0
		 ;ecx:zbufoff ³ mm2:  z2 argb2   argb3 argb2 ³ xmm2: inc3 inc2 inc1 inc0
		 ;edx:  j     ³ mm3:  z3 argb3               ³ xmm3:  r3   r2   r1   r0
		 ;esi:  -     ³ mm4:              z1    z0   ³ xmm4:            z3   z2
		 ;edi:scroff  ³ mm5:              z3    z2   ³ xmm5:
		 ;ebp:  -     ³ mm6: plc1 plc0               ³ xmm6:
		 ;esp:  -     ³ mm7: inc1 inc0               ³ xmm7:  z3   z2   z1   z0

		movq mm7bak, mm7
beg4h:pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm4, [eax+edx*8]
		pextrw eax, mm6, 3
		mov eax, angstart[eax*4]
		movq mm1, [eax+edx*8]
		paddd mm6, mm7bak
		pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm5, [eax+edx*8]
		pextrw eax, mm6, 3
		mov eax, angstart[eax*4]
		movq mm3, [eax+edx*8]
		paddd mm6, mm7bak

		movq mm0, mm4
		movq mm2, mm5

			;Do Z
		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

			;Do colors
			;mm4:dist1 dist0
			;mm5:dist3 dist2
		pxor mm7, mm7
		punpcklbw mm0, mm7
		punpcklbw mm1, mm7
		punpcklbw mm2, mm7
		punpcklbw mm3, mm7

		movq mm7, fogcol
		psubw mm7, mm0
		paddw mm7, mm7
		pextrw eax, mm4, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm0, mm7

		movq mm7, fogcol
		psubw mm7, mm1
		paddw mm7, mm7
		pextrw eax, mm4, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm1, mm7

		movq mm7, fogcol
		psubw mm7, mm2
		paddw mm7, mm7
		pextrw eax, mm5, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm2, mm7

		movq mm7, fogcol
		psubw mm7, mm3
		paddw mm7, mm7
		pextrw eax, mm5, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm3, mm7

		packuswb mm0, mm1
		packuswb mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		add edi, 16
		cmp edi, esi
		jbe short beg4h
		add esi, 16
		cmp edi, esi
		jae endh

		movq mm7, mm7bak
		psrad mm7, 1    ;Restore mm7 from incr*2 to just incr for single loop
skip4h:
beg1h:
		pextrw eax, mm6, 1
		paddd mm6, mm7
		mov eax, angstart[eax*4]

			;Z
		cvtsi2ss xmm7, [eax+edx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7

			;Col
		punpcklbw mm0, [eax+edx*8]
		psrlw mm0, 8
		movq mm1, fogcol
		psubw mm1, mm0
		paddw mm1, mm1
		mov eax, [eax+edx*8+4]
		shr eax, 16+4
		pmulhw mm1, foglut[eax*8]
		paddw mm0, mm1
		packuswb mm0, mm1
		movd [edi], mm0

		add edi, 4
		cmp edi, esi
		jb short beg1h
endh: pop edi
		pop esi
	}
}

void hrendz3dn (long sx, long sy, long p1, long plc, long incr, long j)
{
	_asm
	{
		push esi
		push edi
		mov eax, sy
		mov eax, ylookup[eax*4]
		add eax, frameplace
		mov esi, p1
		lea esi, [eax+esi*4]    ;esi = p1
		mov edi, sx
		lea edi, [eax+edi*4]    ;edi = p0

		movd mm0, sx
		punpckldq mm0, sy
		pi2fd mm0, mm0          ;mm0: (float)sy (float)sx
		pshufw mm2, mm0, 0xee   ;mm2: (float)sy (float)sy
		punpckldq mm0, mm0      ;mm0: (float)sx (float)sx
		movd mm1, optistrx
		punpckldq mm1, optistry
		pfmul mm0, mm1          ;mm0: (float)sx*optistry (float)sx*optistrx
		movd mm3, optiheix
		punpckldq mm3, optiheiy
		pfmul mm2, mm3          ;mm2: (float)sy*optiheiy (float)sy*optiheix
		pfadd mm0, mm2
		movd mm3, optiaddx
		punpckldq mm3, optiaddy ;mm3: optiaddy optiaddx
		pfadd mm0, mm3          ;mm0: diry diry

		movd mm6, plc
		movd mm7, incr
		mov ecx, zbufoff
		mov edx, j

beg:  pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]   ;mm2:      dist       col
		pshufw mm3, mm2, 0xee   ;mm3:         ?      dist
		pi2fd mm3, mm3          ;mm3:         ?   (f)dist
		movq mm4, mm0           ;mm4:      diry      dirx
		pfmul mm4, mm4          ;mm4:    diry^2    dirx^2
		pfadd mm0, mm1          ;mm0: dirx+optx diry+opty (unrelated)
		pfacc mm4, mm4          ;mm4: (x^2+y^2)   x^2+y^2
		pfrsqrt mm4, mm4        ;mm4: 1/sqrt(*) 1/sqrt(*)
		pfmul mm3, mm4          ;mm3:         0    zvalue
		paddd mm6, mm7          ;mm6:            plc+incr (unrelated)
		movd [edi], mm2
		movd [edi+ecx], mm3
		add edi, 4
		cmp edi, esi
		jb short beg
		pop edi
		pop esi
	}
}

void hrendzfog3dn (long sx, long sy, long p1, long plc, long incr, long j)
{
	_asm
	{
		push esi
		push edi
		mov eax, sy
		mov eax, ylookup[eax*4]
		add eax, frameplace
		mov esi, p1
		lea esi, [eax+esi*4]    ;esi = p1
		mov edi, sx
		lea edi, [eax+edi*4]    ;edi = p0

		movd mm0, sx
		punpckldq mm0, sy
		pi2fd mm0, mm0          ;mm0: (float)sy (float)sx
		pshufw mm2, mm0, 0xee   ;mm2: (float)sy (float)sy
		punpckldq mm0, mm0      ;mm0: (float)sx (float)sx
		movd mm1, optistrx
		punpckldq mm1, optistry
		pfmul mm0, mm1          ;mm0: (float)sx*optistry (float)sx*optistrx
		movd mm3, optiheix
		punpckldq mm3, optiheiy
		pfmul mm2, mm3          ;mm2: (float)sy*optiheiy (float)sy*optiheix
		pfadd mm0, mm2
		movd mm3, optiaddx
		punpckldq mm3, optiaddy ;mm3: optiaddy optiaddx
		pfadd mm0, mm3          ;mm0: diry diry

		pxor mm5, mm5

		movd mm6, plc
		movd mm7, incr
		mov ecx, zbufoff
		mov edx, j

beg:  pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]   ;mm2:      dist       col
		pshufw mm3, mm2, 0xee   ;mm3:         ?      dist
		pi2fd mm3, mm3          ;mm3:         ?   (f)dist
		movq mm4, mm0           ;mm4:      diry      dirx
		pfmul mm4, mm4          ;mm4:    diry^2    dirx^2
		pfadd mm0, mm1          ;mm0: dirx+optx diry+opty (unrelated)
		pfacc mm4, mm4          ;mm4: (x^2+y^2)   x^2+y^2
		pfrsqrt mm4, mm4        ;mm4: 1/sqrt(*) 1/sqrt(*)
		pfmul mm3, mm4          ;mm3:         0    zvalue
		paddd mm6, mm7          ;mm6:            plc+incr (unrelated)

			;Extra calculations for fog
		pextrw eax, mm2, 3
		punpcklbw mm2, mm5
		movq mm4, fogcol
		psubw mm4, mm2
		paddw mm4, mm4
		shr eax, 4
		pmulhw mm4, foglut[eax*8]
		paddw mm2, mm4
		packuswb mm2, mm4

		movd [edi], mm2
		movd [edi+ecx], mm3
		add edi, 4
		cmp edi, esi
		jb short beg
		pop edi
		pop esi
	}
}

void vrendzsse (long sx, long sy, long p1, long iplc, long iinc)
{
	_asm
	{
		push ebx
		push esi
		push edi
begvasm_p3:
		mov esi, sx
		mov eax, sy
		mov edx, p1
		mov ecx, ylookup[eax*4]
		add ecx, frameplace
		lea edx, [ecx+edx*4]
		lea edi, [ecx+esi*4]

		mov ecx, esi
		and ecx, 0xfffffffc
		cvtsi2ss xmm0, ecx
		cvtsi2ss xmm4, eax
		movss xmm1, xmm0
		movss xmm5, xmm4
		mulss xmm0, optistrx
		mulss xmm1, optistry
		mulss xmm4, optiheix
		mulss xmm5, optiheiy
		addss xmm0, optiaddx
		addss xmm1, optiaddy
		addss xmm0, xmm4
		addss xmm1, xmm5

		shufps xmm0, xmm0, 0
		shufps xmm1, xmm1, 0
		movaps xmm2, opti4asm[2*16]
		movaps xmm3, opti4asm[3*16]
		addps xmm0, opti4asm[0*16]
		addps xmm1, opti4asm[1*16]
			;xmm0 =  xmm0      ^2 +  xmm1      ^2        (p)
			;xmm2 = (xmm0+xmm2)^2 + (xmm1+xmm3)^2 - xmm0 (v)
			;xmm1 = ...                                  (a)
		addps xmm2, xmm0  ;This block converts inner loop...
		addps xmm3, xmm1  ;from: 1 / sqrt(x*x + y*y), x += xi, y += yi;
		mulps xmm0, xmm0  ;  to: 1 / sqrt(p), p += v, v += a;
		mulps xmm1, xmm1
		mulps xmm2, xmm2
		mulps xmm3, xmm3
		addps xmm0, xmm1
		movaps xmm1, opti4asm[4*16]
		addps xmm2, xmm3
		subps xmm2, xmm0

		mov p1, edx
		mov ecx, zbufoff
		shl esi, 2
		add esi, uurend
		mov ebx, iplc

		cmp edi, edx
		jae short endv

			;Do first 0-3 pixels to align unrolled loop of 4
		test edi, 15
		jz short skip1va

		test edi, 8
		jz short skipshufc
		shufps xmm0, xmm0, 0x4e ;rotate right by 2
skipshufc:
		test edi, 4
		jz short skipshufd
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
skipshufd:

beg1va:
		mov edx, [esi]
		mov eax, [esi+MAXXDIM*4]
		add eax, edx
		sar edx, 16
		mov edx, angstart[edx*4]
		mov [esi], eax
		mov eax, [edx+ebx*8]
		mov [edi], eax
		cvtsi2ss xmm7, [edx+ebx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7
		add ebx, iinc
		add esi, 4
		add edi, 4
		cmp edi, p1
		jz short endv
		test edi, 15
		jnz short beg1va

		addps xmm0, xmm2
		addps xmm2, xmm1
skip1va:
		lea edx, [edi+16]
		cmp edx, p1
		ja short prebeg1v

		cmp iinc, 0
		jl short beg4vn

beg4vp:
		movq mm6, [esi]
		movq mm7, [esi+8]
		pextrw eax, mm6, 1
		pextrw edx, mm6, 3
		paddd mm6, [esi+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm0, [eax+ebx*8]
		movq mm1, [edx+ebx*8+8]
		pextrw eax, mm7, 1
		pextrw edx, mm7, 3
		paddd mm7, [esi+8+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm2, [eax+ebx*8+16]
		movq mm3, [edx+ebx*8+24]
		add ebx, 4

		movq mm4, mm0
		movq mm5, mm2
		punpckldq mm0, mm1
		punpckldq mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

		movq [esi], mm6
		movq [esi+8], mm7

		add esi, 16
		add edi, 16
		lea edx, [edi+16]
		cmp edx, p1
		jbe short beg4vp
		cmp edi, p1
		jae short endv
		jmp short prebeg1v

beg4vn:
		movq mm6, [esi]
		movq mm7, [esi+8]
		pextrw eax, mm6, 1
		pextrw edx, mm6, 3
		paddd mm6, [esi+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm0, [eax+ebx*8]
		movq mm1, [edx+ebx*8-8]
		pextrw eax, mm7, 1
		pextrw edx, mm7, 3
		paddd mm7, [esi+8+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm2, [eax+ebx*8-16]
		movq mm3, [edx+ebx*8-24]
		sub ebx, 4

		movq mm4, mm0
		movq mm5, mm2
		punpckldq mm0, mm1
		punpckldq mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

		movq [esi], mm6
		movq [esi+8], mm7

		add esi, 16
		add edi, 16
		lea edx, [edi+16]
		cmp edx, p1
		jbe short beg4vn
		cmp edi, p1
		jae short endv

prebeg1v:
beg1v:
		mov edx, [esi]
		mov eax, [esi+MAXXDIM*4]
		add eax, edx
		sar edx, 16
		mov edx, angstart[edx*4]
		mov [esi], eax
		mov eax, [edx+ebx*8]
		mov [edi], eax
		cvtsi2ss xmm7, [edx+ebx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7
		add ebx, iinc
		add esi, 4
		add edi, 4
		cmp edi, p1
		jne short beg1v
endv: pop edi
		pop esi
		pop ebx
	}
}

void vrendzfogsse (long sx, long sy, long p1, long iplc, long iinc)
{
	_asm
	{
		push ebx
		push esi
		push edi
begvasm_p3:
		mov esi, sx
		mov eax, sy
		mov edx, p1
		mov ecx, ylookup[eax*4]
		add ecx, frameplace
		lea edx, [ecx+edx*4]
		lea edi, [ecx+esi*4]

		mov ecx, esi
		and ecx, 0xfffffffc
		cvtsi2ss xmm0, ecx
		cvtsi2ss xmm4, eax
		movss xmm1, xmm0
		movss xmm5, xmm4
		mulss xmm0, optistrx
		mulss xmm1, optistry
		mulss xmm4, optiheix
		mulss xmm5, optiheiy
		addss xmm0, optiaddx
		addss xmm1, optiaddy
		addss xmm0, xmm4
		addss xmm1, xmm5

		shufps xmm0, xmm0, 0
		shufps xmm1, xmm1, 0
		movaps xmm2, opti4asm[2*16]
		movaps xmm3, opti4asm[3*16]
		addps xmm0, opti4asm[0*16]
		addps xmm1, opti4asm[1*16]
			;xmm0 =  xmm0      ^2 +  xmm1      ^2        (p)
			;xmm2 = (xmm0+xmm2)^2 + (xmm1+xmm3)^2 - xmm0 (v)
			;xmm1 = ...                                  (a)
		addps xmm2, xmm0  ;This block converts inner loop...
		addps xmm3, xmm1  ;from: 1 / sqrt(x*x + y*y), x += xi, y += yi;
		mulps xmm0, xmm0  ;  to: 1 / sqrt(p), p += v, v += a;
		mulps xmm1, xmm1
		mulps xmm2, xmm2
		mulps xmm3, xmm3
		addps xmm0, xmm1
		movaps xmm1, opti4asm[4*16]
		addps xmm2, xmm3
		subps xmm2, xmm0

		mov p1, edx
		mov ecx, zbufoff
		shl esi, 2
		add esi, uurend
		mov ebx, iplc

		cmp edi, edx
		jae short endv

			;Do first 0-3 pixels to align unrolled loop of 4
		test edi, 15
		jz short skip1va

		test edi, 8
		jz short skipshufc
		shufps xmm0, xmm0, 0x4e ;rotate right by 2
skipshufc:
		test edi, 4
		jz short skipshufd
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
skipshufd:

beg1va:
		mov edx, [esi]
		mov eax, [esi+MAXXDIM*4]
		add eax, edx
		sar edx, 16
		mov edx, angstart[edx*4]
		mov [esi], eax

			;Z
		cvtsi2ss xmm7, [edx+ebx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7

			;Col
		punpcklbw mm0, [edx+ebx*8]
		psrlw mm0, 8
		movq mm1, fogcol
		psubw mm1, mm0
		paddw mm1, mm1
		mov eax, [edx+ebx*8+4]
		shr eax, 16+4
		pmulhw mm1, foglut[eax*8]
		paddw mm0, mm1
		packuswb mm0, mm1
		movd [edi], mm0

		add ebx, iinc
		add esi, 4
		add edi, 4
		cmp edi, p1
		jz short endv
		test edi, 15
		jnz short beg1va

		addps xmm0, xmm2
		addps xmm2, xmm1
skip1va:
		lea edx, [edi+16]
		cmp edx, p1
		ja short prebeg1v

		cmp iinc, 0
		jl short beg4vn

beg4vp:
		movq mm6, [esi]
		movq mm7, [esi+8]
		pextrw eax, mm6, 1
		pextrw edx, mm6, 3
		paddd mm6, [esi+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm4, [eax+ebx*8]
		movq mm1, [edx+ebx*8+8]
		pextrw eax, mm7, 1
		pextrw edx, mm7, 3
		paddd mm7, [esi+8+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm5, [eax+ebx*8+16]
		movq mm3, [edx+ebx*8+24]
		add ebx, 4

			;Do Z
		movq mm0, mm4
		movq mm2, mm5
		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

		movq [esi], mm6
		movq [esi+8], mm7

			;Do color
		pxor mm7, mm7
		punpcklbw mm0, mm7
		punpcklbw mm1, mm7
		punpcklbw mm2, mm7
		punpcklbw mm3, mm7

		movq mm7, fogcol
		psubw mm7, mm0
		paddw mm7, mm7
		pextrw eax, mm4, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm0, mm7

		movq mm7, fogcol
		psubw mm7, mm1
		paddw mm7, mm7
		pextrw eax, mm4, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm1, mm7

		movq mm7, fogcol
		psubw mm7, mm2
		paddw mm7, mm7
		pextrw eax, mm5, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm2, mm7

		movq mm7, fogcol
		psubw mm7, mm3
		paddw mm7, mm7
		pextrw eax, mm5, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm3, mm7

		packuswb mm0, mm1
		packuswb mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		add esi, 16
		add edi, 16
		lea edx, [edi+16]
		cmp edx, p1
		jbe short beg4vp
		cmp edi, p1
		jae short endv
		jmp short prebeg1v

beg4vn:
		movq mm6, [esi]
		movq mm7, [esi+8]
		pextrw eax, mm6, 1
		pextrw edx, mm6, 3
		paddd mm6, [esi+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm4, [eax+ebx*8]
		movq mm1, [edx+ebx*8-8]
		pextrw eax, mm7, 1
		pextrw edx, mm7, 3
		paddd mm7, [esi+8+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm5, [eax+ebx*8-16]
		movq mm3, [edx+ebx*8-24]
		sub ebx, 4

			;Do Z
		movq mm0, mm4
		movq mm2, mm5
		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

		movq [esi], mm6
		movq [esi+8], mm7

			;Do color
		pxor mm7, mm7
		punpcklbw mm0, mm7
		punpcklbw mm1, mm7
		punpcklbw mm2, mm7
		punpcklbw mm3, mm7

		movq mm7, fogcol
		psubw mm7, mm0
		paddw mm7, mm7
		pextrw eax, mm4, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm0, mm7

		movq mm7, fogcol
		psubw mm7, mm1
		paddw mm7, mm7
		pextrw eax, mm4, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm1, mm7

		movq mm7, fogcol
		psubw mm7, mm2
		paddw mm7, mm7
		pextrw eax, mm5, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm2, mm7

		movq mm7, fogcol
		psubw mm7, mm3
		paddw mm7, mm7
		pextrw eax, mm5, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm3, mm7

		packuswb mm0, mm1
		packuswb mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		add esi, 16
		add edi, 16
		lea edx, [edi+16]
		cmp edx, p1
		jbe short beg4vn
		cmp edi, p1
		jae short endv

prebeg1v:
beg1v:
		mov edx, [esi]
		mov eax, [esi+MAXXDIM*4]
		add eax, edx
		sar edx, 16
		mov edx, angstart[edx*4]
		mov [esi], eax

			;Z
		cvtsi2ss xmm7, [edx+ebx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7

			;Col
		punpcklbw mm0, [edx+ebx*8]
		psrlw mm0, 8
		movq mm1, fogcol
		psubw mm1, mm0
		paddw mm1, mm1
		mov eax, [edx+ebx*8+4]
		shr eax, 16+4
		pmulhw mm1, foglut[eax*8]
		paddw mm0, mm1
		packuswb mm0, mm1
		movd [edi], mm0

		add ebx, iinc
		add esi, 4
		add edi, 4
		cmp edi, p1
		jne short beg1v
endv: pop edi
		pop esi
		pop ebx
	}
}

void vrendz3dn (long sx, long sy, long p1, long iplc, long iinc)
{
	_asm
	{
		push ebx
		push esi
		push edi
		mov esi, p1
		mov edi, sx
		cmp edi, esi
		jae short endv
		mov eax, sy
		mov eax, ylookup[eax*4]
		add eax, frameplace
		lea esi, [eax+esi*4]    ;esi = p1
		lea edi, [eax+edi*4]    ;edi = p0

		movd mm0, sx
		punpckldq mm0, sy
		pi2fd mm0, mm0          ;mm0: (float)sy (float)sx
		pshufw mm2, mm0, 0xee   ;mm2: (float)sy (float)sy
		punpckldq mm0, mm0      ;mm0: (float)sx (float)sx
		movd mm1, optistrx
		punpckldq mm1, optistry
		pfmul mm0, mm1          ;mm0: (float)sx*optistry (float)sx*optistrx
		movd mm3, optiheix
		punpckldq mm3, optiheiy
		pfmul mm2, mm3          ;mm2: (float)sy*optiheiy (float)sy*optiheix
		pfadd mm0, mm2
		movd mm3, optiaddx
		punpckldq mm3, optiaddy ;mm3: optiaddy optiaddx
		pfadd mm0, mm3          ;mm0: diry diry

		mov ecx, zbufoff
		mov edx, iplc
		mov ebx, sx
		mov eax, uurend
		lea ebx, [eax+ebx*4]

begv_3dn:
		movd mm5, [ebx]
		pextrw eax, mm5, 1
		paddd mm5, [ebx+MAXXDIM*4]
		movd [ebx], mm5
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]   ;mm2:      dist       col
		pshufw mm3, mm2, 0xee   ;mm3:         ?      dist
		pi2fd mm3, mm3          ;mm3:         ?   (f)dist
		movq mm4, mm0           ;mm4:      diry      dirx
		pfmul mm4, mm4          ;mm4:    diry^2    dirx^2
		pfadd mm0, mm1          ;mm0: dirx+optx diry+opty (unrelated)
		pfacc mm4, mm4          ;mm4: (x^2+y^2)   x^2+y^2
		pfrsqrt mm4, mm4        ;mm4: 1/sqrt(*) 1/sqrt(*)
		pfmul mm3, mm4          ;mm3:         0    zvalue
		movd [edi], mm2
		movd [edi+ecx], mm3
		add edx, iinc
		add ebx, 4
		add edi, 4
		cmp edi, esi
		jb short begv_3dn
endv: pop edi
		pop esi
		pop ebx
	}
}

void vrendzfog3dn (long sx, long sy, long p1, long iplc, long iinc)
{
	_asm
	{
		push ebx
		push esi
		push edi
		mov esi, p1
		mov edi, sx
		cmp edi, esi
		jae short endv
		mov eax, sy
		mov eax, ylookup[eax*4]
		add eax, frameplace
		lea esi, [eax+esi*4]    ;esi = p1
		lea edi, [eax+edi*4]    ;edi = p0

		movd mm0, sx
		punpckldq mm0, sy
		pi2fd mm0, mm0          ;mm0: (float)sy (float)sx
		pshufw mm2, mm0, 0xee   ;mm2: (float)sy (float)sy
		punpckldq mm0, mm0      ;mm0: (float)sx (float)sx
		movd mm1, optistrx
		punpckldq mm1, optistry
		pfmul mm0, mm1          ;mm0: (float)sx*optistry (float)sx*optistrx
		movd mm3, optiheix
		punpckldq mm3, optiheiy
		pfmul mm2, mm3          ;mm2: (float)sy*optiheiy (float)sy*optiheix
		pfadd mm0, mm2
		movd mm3, optiaddx
		punpckldq mm3, optiaddy ;mm3: optiaddy optiaddx
		pfadd mm0, mm3          ;mm0: diry diry

		pxor mm6, mm6

		mov ecx, zbufoff
		mov edx, iplc
		mov ebx, sx
		mov eax, uurend
		lea ebx, [eax+ebx*4]

begv_3dn:
		movd mm5, [ebx]
		pextrw eax, mm5, 1
		paddd mm5, [ebx+MAXXDIM*4]
		movd [ebx], mm5
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]   ;mm2:      dist       col
		pshufw mm3, mm2, 0xee   ;mm3:         ?      dist
		pi2fd mm3, mm3          ;mm3:         ?   (f)dist
		movq mm4, mm0           ;mm4:      diry      dirx
		pfmul mm4, mm4          ;mm4:    diry^2    dirx^2
		pfadd mm0, mm1          ;mm0: dirx+optx diry+opty (unrelated)
		pfacc mm4, mm4          ;mm4: (x^2+y^2)   x^2+y^2
		pfrsqrt mm4, mm4        ;mm4: 1/sqrt(*) 1/sqrt(*)
		pfmul mm3, mm4          ;mm3:         0    zvalue

			;Extra calculations for fog
		pextrw eax, mm2, 3
		punpcklbw mm2, mm6
		movq mm4, fogcol
		psubw mm4, mm2
		paddw mm4, mm4
		shr eax, 4
		pmulhw mm4, foglut[eax*8]
		paddw mm2, mm4
		packuswb mm2, mm4

		movd [edi], mm2
		movd [edi+ecx], mm3
		add edx, iinc
		add ebx, 4
		add edi, 4
		cmp edi, esi
		jb short begv_3dn
endv: pop edi
		pop esi
		pop ebx
	}
}

#endif
/** Set global camera position for future voxlap5 engine calls.
 *  Functions that depend on this include: opticast, drawsprite, spherefill, etc...
 *
 *  The 5th & 6th parameters define the center of the screen projection. This
 *  is the point on the screen that intersects the <ipos + ifor*t> vector.
 *
 *  The last parameter is the focal length - use it to control zoom. If you
 *  want a 90 degree field of view (left to right side of screen), then
 *  set it to half of the screen's width: (xdim*.5).
 *
 *  @param ipo camera position
 *  @param ist camera's unit RIGHT vector
 *  @param ihe camera's unit DOWN vector
 *  @param ifo camera's unit FORWARD vector
 *  @param dahx x-dimension of viewing window
 *  @param dahy y-dimension of viewing window
 *  @param dahz z-dimension of viewing window (should = dahx for 90 degree FOV)
 */
void setcamera (dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo,
					 float dahx, float dahy, float dahz)
{
	long i, j;

	gipos.x = ipo->x; gipos.y = ipo->y; gipos.z = ipo->z;
	gistr.x = ist->x; gistr.y = ist->y; gistr.z = ist->z;
	gihei.x = ihe->x; gihei.y = ihe->y; gihei.z = ihe->z;
	gifor.x = ifo->x; gifor.y = ifo->y; gifor.z = ifo->z;
	gihx = dahx; gihy = dahy; gihz = dahz;

	gixs.x = gistr.x; gixs.y = gihei.x; gixs.z = gifor.x;
	giys.x = gistr.y; giys.y = gihei.y; giys.z = gifor.y;
	gizs.x = gistr.z; gizs.y = gihei.z; gizs.z = gifor.z;
	giadd.x = -(gipos.x*gistr.x + gipos.y*gistr.y + gipos.z*gistr.z);
	giadd.y = -(gipos.x*gihei.x + gipos.y*gihei.y + gipos.z*gihei.z);
	giadd.z = -(gipos.x*gifor.x + gipos.y*gifor.y + gipos.z*gifor.z);

	gcorn[0].x = -gihx*gistr.x - gihy*gihei.x + gihz*gifor.x;
	gcorn[0].y = -gihx*gistr.y - gihy*gihei.y + gihz*gifor.y;
	gcorn[0].z = -gihx*gistr.z - gihy*gihei.z + gihz*gifor.z;
	gcorn[1].x = xres*gistr.x+gcorn[0].x;
	gcorn[1].y = xres*gistr.y+gcorn[0].y;
	gcorn[1].z = xres*gistr.z+gcorn[0].z;
	gcorn[2].x = yres*gihei.x+gcorn[1].x;
	gcorn[2].y = yres*gihei.y+gcorn[1].y;
	gcorn[2].z = yres*gihei.z+gcorn[1].z;
	gcorn[3].x = yres*gihei.x+gcorn[0].x;
	gcorn[3].y = yres*gihei.y+gcorn[0].y;
	gcorn[3].z = yres*gihei.z+gcorn[0].z;
	for(j=0,i=3;j<4;i=j++)
	{
		ginor[i].x = gcorn[i].y*gcorn[j].z - gcorn[i].z*gcorn[j].y;
		ginor[i].y = gcorn[i].z*gcorn[j].x - gcorn[i].x*gcorn[j].z;
		ginor[i].z = gcorn[i].x*gcorn[j].y - gcorn[i].y*gcorn[j].x;
	}
}

/** Render VXL screen (this is where it all happens!)
 *
 *  Make sure you have .VXL loaded in memory by using one of the loadnul(),
 *  loadvxl(), loadbsp(), loaddta() functions.
 *  Also make sure to call setcamera() and setvoxframebuffer() before this.
 */
void opticast ()
{
	float f, ff, cx, cy, fx, fy, gx, gy, x0, y0, x1, y1, x2, y2, x3, y3;
	long i, j, sx, sy, p0, p1, cx16, cy16, kadd, kmul, u, u1, ui;

	if (gifor.z < 0) giforzsgn = -1; else giforzsgn = 1; //giforzsgn = (gifor.z < 0);

	gixyi[0] = (VSID<<2); gixyi[1] = -gixyi[0];
	glipos.x = ((long)gipos.x);
	glipos.y = ((long)gipos.y);
	glipos.z = ((long)gipos.z);
	gpixy = (long)&sptr[glipos.y*VSID + glipos.x];
	ftol(gipos.z*PREC-.5f,&gposz);
	gposxfrac[1] = gipos.x - (float)glipos.x; gposxfrac[0] = 1-gposxfrac[1];
	gposyfrac[1] = gipos.y - (float)glipos.y; gposyfrac[0] = 1-gposyfrac[1];
#if (USEV5ASM == 1)
	for(j=u=0;j<gmipnum;j++,u+=i)
		for(i=0;i<(256>>j)+4;i++)
			gylookup[i+u] = ((((gposz>>j)-i*PREC)>>(16-j))&0x0000ffff);
	gxmip = max(vx5.mipscandist,4)*PREC;
#else
	for(i=0;i<256+4;i++) gylookup[i] = (i*PREC-gposz);
#endif
	gmaxscandist = min(max(vx5.maxscandist,1),2047)*PREC;

#if (USEZBUFFER != 1)
	hrend = hrendnoz; vrend = vrendnoz;
#else
	if (ofogdist < 0)
	{
#if !defined( USE_INTRINSICS )
		if (cputype&(1<<25)) { hrend = hrendzsse; vrend = vrendzsse;}
							 else { hrend = hrendz3dn; vrend = vrendz3dn; }
#else
		hrend = hrendz; vrend = vrendz;	
#endif	
	}
	else
	{
		if (cputype&(1<<25)) { hrend = hrendzfogsse; vrend = vrendzfogsse; }
							 else { hrend = hrendzfog3dn; vrend = vrendzfog3dn; }

	}
#endif
	if (ofogdist < 0) nskypic = skypic;
				  else { nskypic = skyoff = 0; } //Optimization hack: draw sky as pure black when using fog

	gstartv = (char *)*(long *)gpixy;
	if (glipos.z >= gstartv[1])
	{
		do
		{
			if (!gstartv[0]) return;
			gstartv += gstartv[0]*4;
		} while (glipos.z >= gstartv[1]);
		if (glipos.z < gstartv[3]) return;
		gstartz0 = gstartv[3];
	} else gstartz0 = 0;
	gstartz1 = gstartv[1];

	if (gifor.z == 0) f = 32000; else f = gihz/gifor.z;
	f = min(max(f,-32000),32000);
	cx = gistr.z*f + gihx;
	cy = gihei.z*f + gihy;

	wx0 = (float)(-(vx5.anginc)); wx1 = (float)(xres-1+(vx5.anginc));
	wy0 = (float)(-(vx5.anginc)); wy1 = (float)(yres-1+(vx5.anginc));
	ftol(wx0,&iwx0); ftol(wx1,&iwx1);
	ftol(wy0,&iwy0); ftol(wy1,&iwy1);

	fx = wx0-cx; fy = wy0-cy; gx = wx1-cx; gy = wy1-cy;
	x0 = x3 = wx0; y0 = y1 = wy0; x1 = x2 = wx1; y2 = y3 = wy1;
	if (fy < 0)
	{
		if (fx < 0) { f = sqrt(fx*fy); x0 = cx-f; y0 = cy-f; }
		if (gx > 0) { f = sqrt(-gx*fy); x1 = cx+f; y1 = cy-f; }
	}
	if (gy > 0)
	{
		if (gx > 0) { f = sqrt(gx*gy); x2 = cx+f; y2 = cy+f; }
		if (fx < 0) { f = sqrt(-fx*gy); x3 = cx-f; y3 = cy+f; }
	}
	if (x0 > x1) { if (fx < 0) y0 = fx/gx*fy + cy; else y1 = gx/fx*fy + cy; }
	if (y1 > y2) { if (fy < 0) x1 = fy/gy*gx + cx; else x2 = gy/fy*gx + cx; }
	if (x2 < x3) { if (fx < 0) y3 = fx/gx*gy + cy; else y2 = gx/fx*gy + cy; }
	if (y3 < y0) { if (fy < 0) x0 = fy/gy*fx + cx; else x3 = gy/fy*fx + cx; }
		//This makes precision errors cause pixels to overwrite rather than omit
	x0 -= .01; x1 += .01;
	y1 -= .01; y2 += .01;
	x3 -= .01; x2 += .01;
	y0 -= .01; y3 += .01;

	f = (float)PREC / gihz;
	optistrx = gistr.x*f; optiheix = gihei.x*f; optiaddx = gcorn[0].x*f;
	optistry = gistr.y*f; optiheiy = gihei.y*f; optiaddy = gcorn[0].y*f;
#ifdef _MSC_VER
	opti4[0].y = optistrx; opti4[0].z = optistrx*2; opti4[0].z2 = optistrx*3;
	opti4[1].y = optistry; opti4[1].z = optistry*2; opti4[1].z2 = optistry*3;
	opti4[2].x = opti4[2].y = opti4[2].z = opti4[2].z2 = optistrx*4.0f;
	opti4[3].x = opti4[3].y = opti4[3].z = opti4[3].z2 = optistry*4.0f;
	opti4[4].x = opti4[4].y = opti4[4].z = opti4[4].z2 = (optistrx*optistrx + optistry*optistry)*32.0f; //NEW ALGO!
#endif

	ftol(cx*65536,&cx16);
	ftol(cy*65536,&cy16);

	ftol((x1-x0)/vx5.anginc,&j);
	if ((fy < 0) && (j > 0)) //(cx,cy),(x0,wy0),(x1,wy0)
	{
		ff = (x1-x0) / (float)j; grd = 1.0f / (wy0-cy);
		gscanptr = (castdat *)radar; skycurlng = -1; skycurdir = -giforzsgn;
		for(i=0,f=x0+ff*.5f;i<j;f+=ff,i++)
		{
			vline(cx,cy,f,wy0,&p0,&p1);
			if (giforzsgn < 0) angstart[i] = gscanptr+p0; else angstart[i] = gscanptr-p1;
			gscanptr += labs(p1-p0)+1;
		}

		j <<= 16; f = (float)j / ((x1-x0)*grd); ftol((cx-x0)*grd*f,&kadd);
		ftol(cx-.5f,&p1); p0 = lbound0(p1+1,xres); p1 = lbound0(p1,xres);
		ftol(cy-0.50005f,&sy); if (sy >= yres) sy = yres-1;
		ff = (fabs((float)p1-cx)+1)*f/2147483647.0 + cy; //Anti-crash hack
		while ((ff < sy) && (sy >= 0)) sy--;
		if (sy >= 0)
		{
			ftol(f,&kmul);
			for(;sy>=0;sy--) if (isshldiv16safe(kmul,(sy<<16)-cy16)) break; //Anti-crash hack
			if (giforzsgn < 0) i = -sy; else i = sy;
			for(;sy>=0;sy--,i-=giforzsgn)
			{
				ui = shldiv16(kmul,(sy<<16)-cy16);
				u = mulshr16((p0<<16)-cx16,ui)+kadd;
				while ((p0 > 0) && (u >= ui)) { u -= ui; p0--; }
				u1 = (p1-p0)*ui + u;
				while ((p1 < xres) && (u1 < j)) { u1 += ui; p1++; }
				if (p0 < p1) hrend(p0,sy,p1,u,ui,i);
			}
			clearMMX();
		}
	}

	ftol((y2-y1)/vx5.anginc,&j);
	if ((gx > 0) && (j > 0)) //(cx,cy),(wx1,y1),(wx1,y2)
	{
		ff = (y2-y1) / (float)j; grd = 1.0f / (wx1-cx);
		gscanptr = (castdat *)radar; skycurlng = -1; skycurdir = -giforzsgn;
		for(i=0,f=y1+ff*.5f;i<j;f+=ff,i++)
		{
			hline(cx,cy,wx1,f,&p0,&p1);
			if (giforzsgn < 0) angstart[i] = gscanptr-p0; else angstart[i] = gscanptr+p1;
			gscanptr += labs(p1-p0)+1;
		}

		j <<= 16; f = (float)j / ((y2-y1)*grd); ftol((cy-y1)*grd*f,&kadd);
		ftol(cy-.5f,&p1); p0 = lbound0(p1+1,yres); p1 = lbound0(p1,yres);
		ftol(cx+0.50005f,&sx); if (sx < 0) sx = 0;
		ff = (fabs((float)p1-cy)+1)*f/2147483647.0 + cx; //Anti-crash hack
		while ((ff > sx) && (sx < xres)) sx++;
		if (sx < xres)
		{
			ftol(f,&kmul);
			for(;sx<xres;sx++) if (isshldiv16safe(kmul,(sx<<16)-cx16)) break; //Anti-crash hack
			for(;sx<xres;sx++)
			{
				ui = shldiv16(kmul,(sx<<16)-cx16);
				u = mulshr16((p0<<16)-cy16,ui)+kadd;
				while ((p0 > 0) && (u >= ui)) { u -= ui; lastx[--p0] = sx; }
				uurend[sx] = u; uurend[sx+MAXXDIM] = ui; u += (p1-p0)*ui;
				while ((p1 < yres) && (u < j)) { u += ui; lastx[p1++] = sx; }
			}
			if (giforzsgn < 0)
				  { for(sy=p0;sy<p1;sy++) vrend(lastx[sy],sy,xres,lastx[sy],1); }
			else { for(sy=p0;sy<p1;sy++) vrend(lastx[sy],sy,xres,-lastx[sy],-1); }
			clearMMX();
		}
	}

	ftol((x2-x3)/vx5.anginc,&j);
	if ((gy > 0) && (j > 0)) //(cx,cy),(x2,wy1),(x3,wy1)
	{
		ff = (x2-x3) / (float)j; grd = 1.0f / (wy1-cy);
		gscanptr = (castdat *)radar; skycurlng = -1; skycurdir = giforzsgn;
		for(i=0,f=x3+ff*.5f;i<j;f+=ff,i++)
		{
			vline(cx,cy,f,wy1,&p0,&p1);
			if (giforzsgn < 0) angstart[i] = gscanptr-p0; else angstart[i] = gscanptr+p1;
			gscanptr += labs(p1-p0)+1;
		}

		j <<= 16; f = (float)j / ((x2-x3)*grd); ftol((cx-x3)*grd*f,&kadd);
		ftol(cx-.5f,&p1); p0 = lbound0(p1+1,xres); p1 = lbound0(p1,xres);
		ftol(cy+0.50005f,&sy); if (sy < 0) sy = 0;
		ff = (fabs((float)p1-cx)+1)*f/2147483647.0 + cy; //Anti-crash hack
		while ((ff > sy) && (sy < yres)) sy++;
		if (sy < yres)
		{
			ftol(f,&kmul);
			for(;sy<yres;sy++) if (isshldiv16safe(kmul,(sy<<16)-cy16)) break; //Anti-crash hack
			if (giforzsgn < 0) i = sy; else i = -sy;
			for(;sy<yres;sy++,i-=giforzsgn)
			{
				ui = shldiv16(kmul,(sy<<16)-cy16);
				u = mulshr16((p0<<16)-cx16,ui)+kadd;
				while ((p0 > 0) && (u >= ui)) { u -= ui; p0--; }
				u1 = (p1-p0)*ui + u;
				while ((p1 < xres) && (u1 < j)) { u1 += ui; p1++; }
				if (p0 < p1) hrend(p0,sy,p1,u,ui,i);
			}
			clearMMX();
		}
	}

	ftol((y3-y0)/vx5.anginc,&j);
	if ((fx < 0) && (j > 0)) //(cx,cy),(wx0,y3),(wx0,y0)
	{
		ff = (y3-y0) / (float)j; grd = 1.0f / (wx0-cx);
		gscanptr = (castdat *)radar; skycurlng = -1; skycurdir = giforzsgn;
		for(i=0,f=y0+ff*.5f;i<j;f+=ff,i++)
		{
			hline(cx,cy,wx0,f,&p0,&p1);
			if (giforzsgn < 0) angstart[i] = gscanptr+p0; else angstart[i] = gscanptr-p1;
			gscanptr += labs(p1-p0)+1;
		}

		j <<= 16; f = (float)j / ((y3-y0)*grd); ftol((cy-y0)*grd*f,&kadd);
		ftol(cy-.5f,&p1); p0 = lbound0(p1+1,yres); p1 = lbound0(p1,yres);
		ftol(cx-0.50005f,&sx); if (sx >= xres) sx = xres-1;
		ff = (fabs((float)p1-cy)+1)*f/2147483647.0 + cx; //Anti-crash hack
		while ((ff < sx) && (sx >= 0)) sx--;
		if (sx >= 0)
		{
			ftol(f,&kmul);
			for(;sx>=0;sx--) if (isshldiv16safe(kmul,(sx<<16)-cx16)) break; //Anti-crash hack
			for(;sx>=0;sx--)
			{
				ui = shldiv16(kmul,(sx<<16)-cx16);
				u = mulshr16((p0<<16)-cy16,ui)+kadd;
				while ((p0 > 0) && (u >= ui)) { u -= ui; lastx[--p0] = sx; }
				uurend[sx] = u; uurend[sx+MAXXDIM] = ui; u += (p1-p0)*ui;
				while ((p1 < yres) && (u < j)) { u += ui; lastx[p1++] = sx; }
			}
			for(sy=p0;sy<p1;sy++) vrend(0,sy,lastx[sy]+1,0,giforzsgn);
			clearMMX();
		}
	}
}

/**
	//0: asm temp for current x
	//1: asm temp for current y
	//2: bottom (28)
	//3: top    ( 0)
	//4: left   ( 8)
	//5: right  (24)
	//6: up     (12)
	//7: down   (12)
	//setsideshades(0,0,0,0,0,0);
	//setsideshades(0,28,8,24,12,12);
	*/
/** Shade offset for each face of the cube: useful for editing
 *
 * @param sto top face (z minimum) shade offset
 * @param sbo bottom face (z maximum) shade offset
 * @param sle left face (x minimum) shade offset
 * @param sri right face (x maximum) shade offset
 * @param sup top face (y minimum) shade offset
 * @param sdo bottom face (y maximum) shade offset
 */
void setsideshades (char sto, char sbo, char sle, char sri, char sup, char sdo)
{
	((char *)&gcsub[2])[7] = sbo; ((char *)&gcsub[3])[7] = sto;
	((char *)&gcsub[4])[7] = sle; ((char *)&gcsub[5])[7] = sri;
	((char *)&gcsub[6])[7] = sup; ((char *)&gcsub[7])[7] = sdo;
	if (!(sto|sbo|sle|sri|sup|sdo))
	{
		vx5.sideshademode = 0;
		((char *)&gcsub[0])[7] = ((char *)&gcsub[1])[7] = 0x00;
	}
	else vx5.sideshademode = 1;
}


//------------------- editing backup / restore begins ------------------------

void voxdontrestore ()
{
	long i;

	if (backedup == 1)
	{
		for(i=(bacx1-bacx0)*(bacy1-bacy0)-1;i>=0;i--) voxdealloc(bacsptr[i]);
	}
	backedup = -1;
}

void voxrestore ()
{
	long i, j, x, y;
	char *v, *daptr;

	if (backedup == 1)
	{
		i = 0;
		for(y=bacy0;y<bacy1;y++)
		{
			j = y*VSID;
			for(x=bacx0;x<bacx1;x++)
			{
				v = sptr[j+x]; vx5.globalmass += v[1];
				while (v[0]) { v += v[0]*4; vx5.globalmass += v[1]-v[3]; }

				voxdealloc(sptr[j+x]);
				sptr[j+x] = bacsptr[i]; i++;

				v = sptr[j+x]; vx5.globalmass -= v[1];
				while (v[0]) { v += v[0]*4; vx5.globalmass += v[3]-v[1]; }
			}
		}
		if (vx5.vxlmipuse > 1) genmipvxl(bacx0,bacy0,bacx1,bacy1);
	}
	else if (backedup == 2)
	{
		daptr = (char *)bacsptr;
		for(y=bacy0;y<bacy1;y++)
		{
			j = y*VSID;
			for(x=bacx0;x<bacx1;x++)
			{
				for(v=sptr[j+x];v[0];v+=v[0]*4)
					for(i=1;i<v[0];i++) v[(i<<2)+3] = *daptr++;
				for(i=1;i<=v[2]-v[1]+1;i++) v[(i<<2)+3] = *daptr++;
			}
		}
		if (vx5.vxlmipuse > 1) genmipvxl(bacx0,bacy0,bacx1,bacy1);
	}
	backedup = -1;
}

void voxbackup (long x0, long y0, long x1, long y1, long tag)
{
	long i, j, n, x, y;
	char *v, *daptr;

	voxdontrestore();

	x0 = max(x0-2,0); y0 = max(y0-2,0);
	x1 = min(x1+2,VSID); y1 = min(y1+2,VSID);
	if ((x1-x0)*(y1-y0) > 262144) return;

	bacx0 = x0; bacy0 = y0; bacx1 = x1; bacy1 = y1; backtag = tag;

	if (tag&0x10000)
	{
		backedup = 1;
		i = 0;
		for(y=bacy0;y<bacy1;y++)
		{
			j = y*VSID;
			for(x=bacx0;x<bacx1;x++)
			{
				bacsptr[i] = v = sptr[j+x]; i++;
				n = slng(v); sptr[j+x] = voxalloc(n);

				copybuf((void *)v,(void *)sptr[j+x],n>>2);
			}
		}
	}
	else if (tag&0x20000)
	{
		backedup = 2;
			//WARNING!!! Right now, function will crash if saving more than
			//   1<<20 colors :( This needs to be addressed!!!
		daptr = (char *)bacsptr;
		for(y=bacy0;y<bacy1;y++)
		{
			j = y*VSID;
			for(x=bacx0;x<bacx1;x++)
			{
				for(v=sptr[j+x];v[0];v+=v[0]*4)
					for(i=1;i<v[0];i++) *daptr++ = v[(i<<2)+3];
				for(i=1;i<=v[2]-v[1]+1;i++) *daptr++ = v[(i<<2)+3];
			}
		}
	}
	else backedup = 0;
}

//-------------------- editing backup / restore ends -------------------------

static int64_t qmulmip[8] =
{
	0x7fff7fff7fff7fff,0x4000400040004000,0x2aaa2aaa2aaa2aaa,0x2000200020002000,
	0x1999199919991999,0x1555155515551555,0x1249124912491249,0x1000100010001000
};
static long mixc[MAXZDIM>>1][8]; //4K
static long mixn[MAXZDIM>>1];    //0.5K
void genmipvxl (long x0, long y0, long x1, long y1)
{
	long i, n, oldn, x, y, z, xsiz, ysiz, zsiz, oxsiz, oysiz;
	long cz, oz, nz, zz, besti, cstat, curz[4], curzn[4][4], mipnum, mipmax;
	char *v[4], *tv, **sr, **sw, **ssr, **ssw;

	if ((!(x0|y0)) && (x1 == VSID) && (y1 == VSID)) mipmax = vx5.vxlmipuse;
															 else mipmax = gmipnum;
	if (mipmax <= 0) return;
	mipnum = 1;

	vx5.colfunc = curcolfunc;
	xsiz = VSID; ysiz = VSID; zsiz = MAXZDIM;
	ssr = sptr; ssw = sptr+xsiz*ysiz;
	while ((xsiz > 1) && (ysiz > 1) && (zsiz > 1) && (mipnum < mipmax))
	{
		oxsiz = xsiz; xsiz >>= 1;
		oysiz = ysiz; ysiz >>= 1;
						  zsiz >>= 1;

		x0--; if (x0 < 0) x0 = 0;
		y0--; if (y0 < 0) y0 = 0;
		x1++; if (x1 > VSID) x1 = VSID;
		y1++; if (y1 > VSID) y1 = VSID;

		x0 >>= 1; x1 = ((x1+1)>>1);
		y0 >>= 1; y1 = ((y1+1)>>1);
		for(y=y0;y<y1;y++)
		{
			sr = ssr+oxsiz*(y<<1)+(x0<<1);
			sw = ssw+xsiz*y+x0;
			for(x=x0;x<x1;x++)
			{
					//ÚÄÄÄÂÄÄÄÂÄÄÄÂÄÄÄ¿
					//³npt³z1 ³z1c³dum³
					//³ b ³ g ³ r ³ i ³
					//³ b ³ g ³ r ³ i ³
					//³npt³z1 ³z1c³z0 ³
					//³ b ³ g ³ r ³ i ³
					//ÀÄÄÄÁÄÄÄÁÄÄÄÁÄÄÄÙ
				v[0] = sr[      0];
				v[1] = sr[      1];
				v[2] = sr[oysiz  ];
				v[3] = sr[oysiz+1];
				for(i=3;i>=0;i--)
				{
					curz[i] = curzn[i][0] = (long)v[i][1];
					curzn[i][1] = ((long)v[i][2])+1;

					tv = v[i];
					while (1)
					{
						oz = (long)tv[1];
						for(z=oz;z<=((long)tv[2]);z++)
						{
							nz = (z>>1);
							mixc[nz][mixn[nz]++] = *(long *)(&tv[((z-oz)<<2)+4]);
						}
						z = (z-oz) - (((long)tv[0])-1);
						if (!tv[0]) break;
						tv += (((long)tv[0])<<2);
						oz = (long)tv[3];
						for(;z<0;z++)
						{
							nz = ((z+oz)>>1);
							mixc[nz][mixn[nz]++] = *(long *)(&tv[z<<2]);
						}
					}
				}
				cstat = 0; oldn = 0; n = 4; tbuf[3] = 0; z = 0x80000000;
				while (1)
				{
					oz = z;

						//z,besti = min,argmin(curz[0],curz[1],curz[2],curz[3])
					besti = (((unsigned long)(curz[1]-curz[    0]))>>31);
						 i = (((unsigned long)(curz[3]-curz[    2]))>>31)+2;
					besti +=(((( signed long)(curz[i]-curz[besti]))>>31)&(i-besti));
					z = curz[besti]; if (z >= MAXZDIM) break;

					if ((!cstat) && ((z>>1) >= ((oz+1)>>1)))
					{
						if (oz >= 0)
						{
							tbuf[oldn] = ((n-oldn)>>2);
							tbuf[oldn+2]--;
							tbuf[n+3] = ((oz+1)>>1);
							oldn = n; n += 4;
						}
						tbuf[oldn] = 0;
						tbuf[oldn+1] = tbuf[oldn+2] = (z>>1); cz = -1;
					}
					if (cstat&0x1111)
					{
						if (((((long)tbuf[oldn+2])<<1)+1 >= oz) && (cz < 0))
						{
							while ((((long)tbuf[oldn+2])<<1) < z)
							{
								zz = (long)tbuf[oldn+2];

								_asm //*(long *)&tbuf[n] = mixc[zz][rand()%mixn[zz]];
								{    //mixn[zz] = 0;
									mov eax, zz
									mov ecx, mixn[eax*4]
									mov mixn[eax*4], 0
									shl eax, 5
									pxor mm0, mm0
									movq mm2, qmulmip[ecx*8-8]
									pcmpeqb mm6, mm6
									movq mm7, mm0
					 vxlmipbeg0:movd mm1, mixc[eax+ecx*4-4]
									punpcklbw mm1, mm7
									paddw mm0, mm1
									dec ecx
									jnz short vxlmipbeg0
									paddw mm0, mm0
									psubw mm0, mm6 ;rounding bias
									pmulhw mm0, mm2
									packuswb mm0, mm0
									mov eax, n
									movd tbuf[eax], mm0
								}

								tbuf[oldn+2]++; n += 4;
							}
						}
						else
						{
							if (cz < 0) cz = (oz>>1);
							else if ((cz<<1)+1 < oz)
							{
									//Insert fake slab
								tbuf[oldn] = ((n-oldn)>>2);
								tbuf[oldn+2]--;
								tbuf[n] = 0;
								tbuf[n+1] = tbuf[n+2] = tbuf[n+3] = cz;
								oldn = n; n += 4;
								cz = (oz>>1);
							}
							while ((cz<<1) < z)
							{
								_asm //*(long *)&tbuf[n] = mixc[cz][rand()%mixn[cz]];
								{    //mixn[cz] = 0;
									mov eax, cz
									mov ecx, mixn[eax*4]
									mov mixn[eax*4], 0
									shl eax, 5
									pxor mm0, mm0
									movq mm2, qmulmip[ecx*8-8]
									pcmpeqb mm6, mm6
									movq mm7, mm0
					 vxlmipbeg1:movd mm1, mixc[eax+ecx*4-4]
									punpcklbw mm1, mm7
									paddw mm0, mm1
									dec ecx
									jnz short vxlmipbeg1
									paddw mm0, mm0
									psubw mm0, mm6 ;rounding bias
									pmulhw mm0, mm2
									packuswb mm0, mm0
									mov eax, n
									movd tbuf[eax], mm0
								}

								cz++; n += 4;
							}
						}
					}

					i = (besti<<2);
					cstat = (((1<<i)+cstat)&0x3333); //--33--22--11--00
					switch ((cstat>>i)&3)
					{
						case 0: curz[besti] = curzn[besti][0]; break;
						case 1: curz[besti] = curzn[besti][1]; break;
						case 2:
							if (!(v[besti][0])) { curz[besti] = MAXZDIM; }
							else
							{
								tv = v[besti]; i = (((long)tv[2])-((long)tv[1])+1)-(((long)tv[0])-1);
								tv += (((long)tv[0])<<2);
								curz[besti] = ((long)(tv[3])) + i;
								curzn[besti][3] = (long)(tv[3]);
								curzn[besti][0] = (long)(tv[1]);
								curzn[besti][1] = ((long)tv[2])+1;
								v[besti] = tv;
							}
							break;
						case 3: curz[besti] = curzn[besti][3]; break;
						default: _gtfo(); //tells MSVC default can't be reached
					}
				}
				tbuf[oldn+2]--;
				if (cz >= 0)
				{
					tbuf[oldn] = ((n-oldn)>>2);
					tbuf[n] = 0;
					tbuf[n+1] = tbuf[n+3] = cz;
					tbuf[n+2] = cz-1;
					n += 4;
				}

					//De-allocate column (x,y) if it exists
				if (sw[0]) voxdealloc(sw[0]);

					//Allocate & copy to new column (x,y)
				sw[0] = voxalloc(n);
				copybuf((void *)tbuf,(void *)sw[0],n>>2);
				sw++; sr += 2;
			}
			sr += ysiz*2;
		}
		ssr = ssw; ssw += xsiz*ysiz;
		mipnum++; if (mipnum > gmipnum) gmipnum = mipnum;
	}

		//Remove extra mips (bbox must be 0,0,VSID,VSID to get inside this)
	while ((xsiz > 1) && (ysiz > 1) && (zsiz > 1) && (mipnum < gmipnum))
	{
		xsiz >>= 1; ysiz >>= 1; zsiz >>= 1;
		for(i=xsiz*ysiz;i>0;i--)
		{
			if (ssw[0]) voxdealloc(ssw[0]); //De-allocate column if it exists
			ssw++;
		}
		gmipnum--;
	}

	clearMMX();

#if 0 //TEMP HACK!!!
	{
	FILE *fil;
	dpoint3d dp;
	if (!(fil = fopen("temp512.vxl","wb"))) return;
	i = 0x09072000; fwrite(&i,4,1,fil);  //Version
	i = (VSID>>1); fwrite(&i,4,1,fil);
	i = (VSID>>1); fwrite(&i,4,1,fil);
	dp.x = (double)i*.5; dp.y = (double)i*.5; dp.z = (double)i*.5;
	fwrite(&dp,24,1,fil);
	dp.x = 1.0; dp.y = 0.0; dp.z = 0.0; fwrite(&dp,24,1,fil);
	dp.x = 0.0; dp.y = 0.0; dp.z = 1.0; fwrite(&dp,24,1,fil);
	dp.x = 0.0; dp.y =-1.0; dp.z = 0.0; fwrite(&dp,24,1,fil);
	for(i=0;i<(VSID>>1)*(VSID>>1);i++)
		fwrite((void *)sptr[i+VSID*VSID],slng(sptr[i+VSID*VSID]),1,fil);
	fclose(fil);
	}
	gmipnum = 1;
#endif

}

static long min0[VSID], max0[VSID]; //MAXY
static long min1[VSID], max1[VSID]; //MAXX
static long min2[VSID], max2[VSID]; //MAXY

static void canseerange (point3d *p0, point3d *p1)
{
	lpoint3d a, c, d, p, i;
	point3d f, g;
	long cnt, j;

	ftol(p0->x-.5,&a.x); ftol(p0->y-.5,&a.y); ftol(p0->z-.5,&a.z);
	ftol(p1->x-.5,&c.x); ftol(p1->y-.5,&c.y); ftol(p1->z-.5,&c.z);
	cnt = 0;

		  if (c.x <  a.x) { d.x = -1; f.x = p0->x-a.x;   g.x = (p0->x-p1->x)*1024; cnt += a.x-c.x; }
	else if (c.x != a.x) { d.x =  1; f.x = a.x+1-p0->x; g.x = (p1->x-p0->x)*1024; cnt += c.x-a.x; }
	else f.x = g.x = 0;
		  if (c.y <  a.y) { d.y = -1; f.y = p0->y-a.y;   g.y = (p0->y-p1->y)*1024; cnt += a.y-c.y; }
	else if (c.y != a.y) { d.y =  1; f.y = a.y+1-p0->y; g.y = (p1->y-p0->y)*1024; cnt += c.y-a.y; }
	else f.y = g.y = 0;
		  if (c.z <  a.z) { d.z = -1; f.z = p0->z-a.z;   g.z = (p0->z-p1->z)*1024; cnt += a.z-c.z; }
	else if (c.z != a.z) { d.z =  1; f.z = a.z+1-p0->z; g.z = (p1->z-p0->z)*1024; cnt += c.z-a.z; }
	else f.z = g.z = 0;

	ftol(f.x*g.z - f.z*g.x,&p.x); ftol(g.x,&i.x);
	ftol(f.y*g.z - f.z*g.y,&p.y); ftol(g.y,&i.y);
	ftol(f.y*g.x - f.x*g.y,&p.z); ftol(g.z,&i.z);
	for(;cnt;cnt--)
	{
			//use a.x, a.y, a.z
		if (a.x < min0[a.y]) min0[a.y] = a.x;
		if (a.x > max0[a.y]) max0[a.y] = a.x;
		if (a.z < min1[a.x]) min1[a.x] = a.z;
		if (a.z > max1[a.x]) max1[a.x] = a.z;
		if (a.z < min2[a.y]) min2[a.y] = a.z;
		if (a.z > max2[a.y]) max2[a.y] = a.z;

		if (((p.x|p.y) >= 0) && (a.z != c.z)) { a.z += d.z; p.x -= i.x; p.y -= i.y; }
		else if ((p.z >= 0) && (a.x != c.x))  { a.x += d.x; p.x += i.z; p.z -= i.y; }
		else                                  { a.y += d.y; p.y += i.z; p.z += i.x; }
	}
}

void settri (point3d *p0, point3d *p1, point3d *p2, long bakit)
{
	point3d n;
	float f, x0, y0, z0, x1, y1, z1, rx, ry, k0, k1;
	long i, x, y, z, iz0, iz1, minx, maxx, miny, maxy;

	if (p0->x < p1->x) { x0 = p0->x; x1 = p1->x; } else { x0 = p1->x; x1 = p0->x; }
	if (p2->x < x0) x0 = p2->x;
	if (p2->x > x1) x1 = p2->x;
	if (p0->y < p1->y) { y0 = p0->y; y1 = p1->y; } else { y0 = p1->y; y1 = p0->y; }
	if (p2->y < y0) y0 = p2->y;
	if (p2->y > y1) y1 = p2->y;
	if (p0->z < p1->z) { z0 = p0->z; z1 = p1->z; } else { z0 = p1->z; z1 = p0->z; }
	if (p2->z < z0) z0 = p2->z;
	if (p2->z > z1) z1 = p2->z;

	ftol(x0-.5,&minx); ftol(y0-.5,&miny);
	ftol(x1-.5,&maxx); ftol(y1-.5,&maxy);
	vx5.minx = minx; vx5.maxx = maxx+1;
	vx5.miny = miny; vx5.maxy = maxy+1;
	ftol(z0-.5,&vx5.minz); ftol(z1+.5,&vx5.maxz);
	if (bakit) voxbackup(minx,miny,maxx+1,maxy+1,bakit);

	for(i=miny;i<=maxy;i++) { min0[i] = 0x7fffffff; max0[i] = 0x80000000; }
	for(i=minx;i<=maxx;i++) { min1[i] = 0x7fffffff; max1[i] = 0x80000000; }
	for(i=miny;i<=maxy;i++) { min2[i] = 0x7fffffff; max2[i] = 0x80000000; }

	canseerange(p0,p1);
	canseerange(p1,p2);
	canseerange(p2,p0);

	n.x = (p1->z-p0->z)*(p2->y-p1->y) - (p1->y-p0->y) * (p2->z-p1->z);
	n.y = (p1->x-p0->x)*(p2->z-p1->z) - (p1->z-p0->z) * (p2->x-p1->x);
	n.z = (p1->y-p0->y)*(p2->x-p1->x) - (p1->x-p0->x) * (p2->y-p1->y);
	f = 1.0 / sqrt(n.x*n.x + n.y*n.y + n.z*n.z); if (n.z < 0) f = -f;
	n.x *= f; n.y *= f; n.z *= f;

	if (n.z > .01)
	{
		f = -1.0 / n.z; rx = n.x*f; ry = n.y*f;
		k0 = ((n.x>=0)-p0->x)*rx + ((n.y>=0)-p0->y)*ry - ((n.z>=0)-p0->z) + .5;
		k1 = ((n.x< 0)-p0->x)*rx + ((n.y< 0)-p0->y)*ry - ((n.z< 0)-p0->z) - .5;
	}
	else { rx = 0; ry = 0; k0 = -2147000000.0; k1 = 2147000000.0; }

	for(y=miny;y<=maxy;y++)
		for(x=min0[y];x<=max0[y];x++)
		{
			f = (float)x*rx + (float)y*ry; ftol(f+k0,&iz0); ftol(f+k1,&iz1);
			if (iz0 < min1[x]) iz0 = min1[x];
			if (iz1 > max1[x]) iz1 = max1[x];
			if (iz0 < min2[y]) iz0 = min2[y];
			if (iz1 > max2[y]) iz1 = max2[y];

				//set: (x,y,iz0) to (x,y,iz1) (inclusive)
			insslab(scum2(x,y),iz0,iz1+1);
	}
	scum2finish();
	updatebbox(vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,0);
}

	//Known problems:
	//1. Need to test faces for intersections on p1<->p2 line (not just edges)
	//2. Doesn't guarantee that hit point/line is purely air (but very close)
	//3. Piescan is more useful for parts of rope code :/
static long tripind[24] = {0,4,1,5,2,6,3,7,0,2,1,3,4,6,5,7,0,1,2,3,4,5,6,7};
long triscan (point3d *p0, point3d *p1, point3d *p2, point3d *hit, lpoint3d *lhit)
{
	point3d n, d[8], cp2;
	float f, g, x0, x1, y0, y1, rx, ry, k0, k1, fx, fy, fz, pval[8];
	long i, j, k, x, y, z, iz0, iz1, minx, maxx, miny, maxy, didhit;

	didhit = 0;

	if (p0->x < p1->x) { x0 = p0->x; x1 = p1->x; } else { x0 = p1->x; x1 = p0->x; }
	if (p2->x < x0) x0 = p2->x;
	if (p2->x > x1) x1 = p2->x;
	if (p0->y < p1->y) { y0 = p0->y; y1 = p1->y; } else { y0 = p1->y; y1 = p0->y; }
	if (p2->y < y0) y0 = p2->y;
	if (p2->y > y1) y1 = p2->y;
	ftol(x0-.5,&minx); ftol(y0-.5,&miny);
	ftol(x1-.5,&maxx); ftol(y1-.5,&maxy);
	for(i=miny;i<=maxy;i++) { min0[i] = 0x7fffffff; max0[i] = 0x80000000; }
	for(i=minx;i<=maxx;i++) { min1[i] = 0x7fffffff; max1[i] = 0x80000000; }
	for(i=miny;i<=maxy;i++) { min2[i] = 0x7fffffff; max2[i] = 0x80000000; }

	canseerange(p0,p1);
	canseerange(p1,p2);
	canseerange(p2,p0);

	n.x = (p1->z-p0->z)*(p2->y-p1->y) - (p1->y-p0->y) * (p2->z-p1->z);
	n.y = (p1->x-p0->x)*(p2->z-p1->z) - (p1->z-p0->z) * (p2->x-p1->x);
	n.z = (p1->y-p0->y)*(p2->x-p1->x) - (p1->x-p0->x) * (p2->y-p1->y);
	f = 1.0 / sqrt(n.x*n.x + n.y*n.y + n.z*n.z); if (n.z < 0) f = -f;
	n.x *= f; n.y *= f; n.z *= f;

	if (n.z > .01)
	{
		f = -1.0 / n.z; rx = n.x*f; ry = n.y*f;
		k0 = ((n.x>=0)-p0->x)*rx + ((n.y>=0)-p0->y)*ry - ((n.z>=0)-p0->z) + .5;
		k1 = ((n.x< 0)-p0->x)*rx + ((n.y< 0)-p0->y)*ry - ((n.z< 0)-p0->z) - .5;
	}
	else { rx = 0; ry = 0; k0 = -2147000000.0; k1 = 2147000000.0; }

	cp2.x = p2->x; cp2.y = p2->y; cp2.z = p2->z;

	for(y=miny;y<=maxy;y++)
		for(x=min0[y];x<=max0[y];x++)
		{
			f = (float)x*rx + (float)y*ry; ftol(f+k0,&iz0); ftol(f+k1,&iz1);
			if (iz0 < min1[x]) iz0 = min1[x];
			if (iz1 > max1[x]) iz1 = max1[x];
			if (iz0 < min2[y]) iz0 = min2[y];
			if (iz1 > max2[y]) iz1 = max2[y];
			for(z=iz0;z<=iz1;z++)
			{
				if (!isvoxelsolid(x,y,z)) continue;

				for(i=0;i<8;i++)
				{
					d[i].x = (float)(( i    &1)+x);
					d[i].y = (float)(((i>>1)&1)+y);
					d[i].z = (float)(((i>>2)&1)+z);
					pval[i] = (d[i].x-p0->x)*n.x + (d[i].y-p0->y)*n.y + (d[i].z-p0->z)*n.z;
				}
				for(i=0;i<24;i+=2)
				{
					j = tripind[i+0];
					k = tripind[i+1];
					if (((*(long *)&pval[j])^(*(long *)&pval[k])) < 0)
					{
						f = pval[j]/(pval[j]-pval[k]);
						fx = (d[k].x-d[j].x)*f + d[j].x;
						fy = (d[k].y-d[j].y)*f + d[j].y;
						fz = (d[k].z-d[j].z)*f + d[j].z;

							//         (p0->x,p0->y,p0->z)
							//             _|     |_
							//           _|     .   |_
							//         _|  (fx,fy,fz) |_
							//       _|                 |_
							//(p1->x,p1->y,p1->z)-.----(cp2.x,cp2.y,cp2.z)

						if ((fabs(n.z) > fabs(n.x)) && (fabs(n.z) > fabs(n.y)))
						{ //x,y
						  // ix = p1->x + (cp2.x-p1->x)*t;
						  // iy = p1->y + (cp2.y-p1->y)*t;
						  //(iz = p1->z + (cp2.z-p1->z)*t;)
						  // ix = p0->x + (fx-p0->x)*u;
						  // iy = p0->y + (fy-p0->y)*u;
						  // (p1->x-cp2.x)*t + (fx-p0->x)*u = p1->x-p0->x;
						  // (p1->y-cp2.y)*t + (fy-p0->y)*u = p1->y-p0->y;

							f = (p1->x-cp2.x)*(fy-p0->y) - (p1->y-cp2.y)*(fx-p0->x);
							if ((*(long *)&f) == 0) continue;
							f = 1.0 / f;
							g = ((p1->x-cp2.x)*(p1->y-p0->y) - (p1->y-cp2.y)*(p1->x-p0->x))*f;
							//NOTE: The following trick assumes g not * or / by f!
							//if (((*(long *)&g)-(*(long *)&f))^(*(long *)&f)) >= 0) continue;
							if ((*(long *)&g) < 0x3f800000) continue;
							g = ((p1->x-p0->x)*(fy-p0->y) - (p1->y-p0->y)*(fx-p0->x))*f;
						}
						else if (fabs(n.y) > fabs(n.x))
						{ //x,z
							f = (p1->x-cp2.x)*(fz-p0->z) - (p1->z-cp2.z)*(fx-p0->x);
							if ((*(long *)&f) == 0) continue;
							f = 1.0 / f;
							g = ((p1->x-cp2.x)*(p1->z-p0->z) - (p1->z-cp2.z)*(p1->x-p0->x))*f;
							if ((*(long *)&g) < 0x3f800000) continue;
							g = ((p1->x-p0->x)*(fz-p0->z) - (p1->z-p0->z)*(fx-p0->x))*f;
						}
						else
						{ //y,z
							f = (p1->y-cp2.y)*(fz-p0->z) - (p1->z-cp2.z)*(fy-p0->y);
							if ((*(long *)&f) == 0) continue;
							f = 1.0 / f;
							g = ((p1->y-cp2.y)*(p1->z-p0->z) - (p1->z-cp2.z)*(p1->y-p0->y))*f;
							if ((*(long *)&g) < 0x3f800000) continue;
							g = ((p1->y-p0->y)*(fz-p0->z) - (p1->z-p0->z)*(fy-p0->y))*f;
						}
						if ((*(unsigned long *)&g) >= 0x3f800000) continue;
						(hit->x) = fx; (hit->y) = fy; (hit->z) = fz;
						(lhit->x) = x; (lhit->y) = y; (lhit->z) = z; didhit = 1;
						(cp2.x) = (cp2.x-p1->x)*g + p1->x;
						(cp2.y) = (cp2.y-p1->y)*g + p1->y;
						(cp2.z) = (cp2.z-p1->z)*g + p1->z;
					}
				}
			}
		}
	return(didhit);
}



#define LPATBUFSIZ 14
static lpoint2d *patbuf;
#define LPATHASHSIZ 12
static lpoint3d *pathashdat;
static long *pathashead, pathashcnt, pathashmax;

static void initpathash ()
{
	patbuf = (lpoint2d *)radar;
	pathashead = (long *)(((long)patbuf)+(1<<LPATBUFSIZ)*sizeof(lpoint2d));
	pathashdat = (lpoint3d *)(((long)pathashead)+((1<<LPATHASHSIZ)*4));
	pathashmax = ((max((MAXXDIM*MAXYDIM*27)>>1,(VSID+4)*3*256*4)-((1<<LPATBUFSIZ)*sizeof(lpoint2d))-(1<<LPATHASHSIZ)*4)/12);
	memset(pathashead,-1,(1<<LPATHASHSIZ)*4);
	pathashcnt = 0;
}

static long readpathash (long i)
{
	long j = (((i>>LPATHASHSIZ)-i) & ((1<<LPATHASHSIZ)-1));
	for(j=pathashead[j];j>=0;j=pathashdat[j].x)
		if (pathashdat[j].y == i) return(pathashdat[j].z);
	return(-1);
}

static void writepathash (long i, long v)
{
	long k, j = (((i>>LPATHASHSIZ)-i) & ((1<<LPATHASHSIZ)-1));
	for(k=pathashead[j];k>=0;k=pathashdat[k].x)
		if (pathashdat[k].y == i) { pathashdat[k].z = v; return; }
	pathashdat[pathashcnt].x = pathashead[j]; pathashead[j] = pathashcnt;
	pathashdat[pathashcnt].y = i;
	pathashdat[pathashcnt].z = v;
	pathashcnt++;
}

static signed char cdir[26*4] = //sqrt(2) =~ 58/41, sqrt(3) =~ 71/41;
{
	-1, 0, 0,41,  1, 0, 0,41,  0,-1, 0,41,  0, 1, 0,41,  0, 0,-1,41,  0, 0, 1,41,
	-1,-1, 0,58, -1, 1, 0,58, -1, 0,-1,58, -1, 0, 1,58,  0,-1,-1,58,  0,-1, 1,58,
	 1,-1, 0,58,  1, 1, 0,58,  1, 0,-1,58,  1, 0, 1,58,  0, 1,-1,58,  0, 1, 1,58,
	-1,-1,-1,71, -1,-1, 1,71, -1, 1,-1,71, -1, 1, 1,71,
	 1,-1,-1,71,  1,-1, 1,71,  1, 1,-1,71,  1, 1, 1,71,
};

long findpath (long *pathpos, long pathmax, lpoint3d *p1, lpoint3d *p0)
{
	long i, j, k, x, y, z, c, nc, xx, yy, zz, bufr, bufw, pcnt;

	if (!(getcube(p0->x,p0->y,p0->z)&~1))
	{
		for(i=5;i>=0;i--)
		{
			x = p0->x+(long)cdir[i*4]; y = p0->y+(long)cdir[i*4+1]; z = p0->z+(long)cdir[i*4+2];
			if (getcube(x,y,z)&~1) { p0->x = x; p0->y = y; p0->z = z; break; }
		}
		if (i < 0) return(0);
	}
	if (!(getcube(p1->x,p1->y,p1->z)&~1))
	{
		for(i=5;i>=0;i--)
		{
			x = p1->x+(long)cdir[i*4]; y = p1->y+(long)cdir[i*4+1]; z = p1->z+(long)cdir[i*4+2];
			if (getcube(x,y,z)&~1) { p1->x = x; p1->y = y; p1->z = z; break; }
		}
		if (i < 0) return(0);
	}

	initpathash();
	j = (p0->x*VSID + p0->y)*MAXZDIM+p0->z;
	patbuf[0].x = j; patbuf[0].y = 0; bufr = 0; bufw = 1;
	writepathash(j,0);
	do
	{
		j = patbuf[bufr&((1<<LPATBUFSIZ)-1)].x;
		x = j/(VSID*MAXZDIM); y = ((j/MAXZDIM)&(VSID-1)); z = (j&(MAXZDIM-1));
		c = patbuf[bufr&((1<<LPATBUFSIZ)-1)].y; bufr++;
		for(i=0;i<26;i++)
		{
			xx = x+(long)cdir[i*4]; yy = y+(long)cdir[i*4+1]; zz = z+(long)cdir[i*4+2];
			j = (xx*VSID + yy)*MAXZDIM+zz;

			//nc = c+(long)cdir[i*4+3]; //More accurate but lowers max distance a lot!
			//if (((k = getcube(xx,yy,zz))&~1) && ((unsigned long)nc < (unsigned long)readpathash(j)))

			if (((k = getcube(xx,yy,zz))&~1) && (readpathash(j) < 0))
			{
				nc = c+(long)cdir[i*4+3];
				if ((xx == p1->x) && (yy == p1->y) && (zz == p1->z)) { c = nc; goto pathfound; }
				writepathash(j,nc);
				if (pathashcnt >= pathashmax) return(0);
				patbuf[bufw&((1<<LPATBUFSIZ)-1)].x = (xx*VSID + yy)*MAXZDIM+zz;
				patbuf[bufw&((1<<LPATBUFSIZ)-1)].y = nc; bufw++;
			}
		}
	} while (bufr != bufw);

pathfound:
	if (pathmax <= 0) return(0);
	pathpos[0] = (p1->x*VSID + p1->y)*MAXZDIM+p1->z; pcnt = 1;
	x = p1->x; y = p1->y; z = p1->z;
	do
	{
		for(i=0;i<26;i++)
		{
			xx = x+(long)cdir[i*4]; yy = y+(long)cdir[i*4+1]; zz = z+(long)cdir[i*4+2];
			nc = c-(long)cdir[i*4+3];
			if (readpathash((xx*VSID + yy)*MAXZDIM+zz) == nc)
			{
				if (pcnt >= pathmax) return(0);
				pathpos[pcnt] = (xx*VSID + yy)*MAXZDIM+zz; pcnt++;
				x = xx; y = yy; z = zz; c = nc; break;
			}
		}
	} while (i < 26);
	if (pcnt >= pathmax) return(0);
	pathpos[pcnt] = (p0->x*VSID + p0->y)*MAXZDIM+p0->z;
	return(pcnt+1);
}

//---------------------------------------------------------------------

static unsigned short xyoffs[256][256+1];
void setkvx (const char *filename, long ox, long oy, long oz, long rot, long bakit)
{
	long i, j, x, y, z, xsiz, ysiz, zsiz, longpal[256], zleng, oldz, vis;
	long d[3], k[9], x0, y0, z0, x1, y1, z1;
	char ch, typ;
	FILE *fp;

	typ = filename[strlen(filename)-3]; if (typ == 'k') typ = 'K';

	if (!(fp = fopen(filename,"rb"))) return;

	fseek(fp,-768,SEEK_END);
	for(i=0;i<255;i++)
	{
		longpal[i]  = (((long)fgetc(fp))<<18);
		longpal[i] += (((long)fgetc(fp))<<10);
		longpal[i] += (((long)fgetc(fp))<< 2) + 0x80000000;
	}
	longpal[255] = 0x7ffffffd;

	if (typ == 'K') //Load .KVX file
	{
		fseek(fp,4,SEEK_SET);
		fread(&xsiz,4,1,fp);
		fread(&ysiz,4,1,fp);
		fread(&zsiz,4,1,fp);
		fseek(fp,((xsiz+1)<<2)+28,SEEK_SET);
		for(i=0;i<xsiz;i++) fread(&xyoffs[i][0],(ysiz+1)<<1,1,fp);
	}
	else           //Load .VOX file
	{
		fseek(fp,0,SEEK_SET);
		fread(&xsiz,4,1,fp);
		fread(&ysiz,4,1,fp);
		fread(&zsiz,4,1,fp);
	}

		//rot: low 3 bits for axis negating, high 6 states for axis swapping
		//k[0], k[3], k[6] are indeces
		//k[1], k[4], k[7] are xors
		//k[2], k[5], k[8] are adds
	switch (rot&~7)
	{
		case  0: k[0] = 0; k[3] = 1; k[6] = 2; break; //can use scum!
		case  8: k[0] = 1; k[3] = 0; k[6] = 2; break; //can use scum!
		case 16: k[0] = 0; k[3] = 2; k[6] = 1; break;
		case 24: k[0] = 2; k[3] = 0; k[6] = 1; break;
		case 32: k[0] = 1; k[3] = 2; k[6] = 0; break;
		case 40: k[0] = 2; k[3] = 1; k[6] = 0; break;
		default: _gtfo(); //tells MSVC default can't be reached
	}
	k[1] = ((rot<<31)>>31);
	k[4] = ((rot<<30)>>31);
	k[7] = ((rot<<29)>>31);

	d[0] = xsiz; d[1] = ysiz; d[2] = zsiz;
	k[2] = ox-((d[k[0]]>>1)^k[1]);
	k[5] = oy-((d[k[3]]>>1)^k[4]);
	k[8] = oz-((d[k[6]]>>1)^k[7]); k[8] -= (d[k[6]]>>1);

	d[0] = d[1] = d[2] = 0;
	x0 = x1 = (d[k[0]]^k[1])+k[2];
	y0 = y1 = (d[k[3]]^k[4])+k[5];
	z0 = z1 = (d[k[6]]^k[7])+k[8];
	d[0] = xsiz; d[1] = ysiz; d[2] = zsiz;
	x0 = min(x0,(d[k[0]]^k[1])+k[2]); x1 = max(x1,(d[k[0]]^k[1])+k[2]);
	y0 = min(y0,(d[k[3]]^k[4])+k[5]); y1 = max(y1,(d[k[3]]^k[4])+k[5]);
	z0 = min(z0,(d[k[6]]^k[7])+k[8]); z1 = max(z1,(d[k[6]]^k[7])+k[8]);
	if (x0 < 1) { i = 1-x0; x0 += i; x1 += i; k[2] += i; }
	if (y0 < 1) { i = 1-y0; y0 += i; y1 += i; k[5] += i; }
	if (z0 < 0) { i = 0-z0; z0 += i; z1 += i; k[8] += i; }
	if (x1 > VSID-2)    { i = VSID-2-x1; x0 += i; x1 += i; k[2] += i; }
	if (y1 > VSID-2)    { i = VSID-2-y1; y0 += i; y1 += i; k[5] += i; }
	if (z1 > MAXZDIM-1) { i = MAXZDIM-1-z1; z0 += i; z1 += i; k[8] += i; }

	vx5.minx = x0; vx5.maxx = x1+1;
	vx5.miny = y0; vx5.maxy = y1+1;
	vx5.minz = z0; vx5.maxz = z1+1;
	if (bakit) voxbackup(x0,y0,x1+1,y1+1,bakit);

	j = (!(k[3]|(rot&3))); //if (j) { can use scum/scumfinish! }

	for(x=0;x<xsiz;x++)
	{
		d[0] = x;
		for(y=0;y<ysiz;y++)
		{
			d[1] = y;
			if (k[6] == 2) //can use scum!
			{
				clearbuf((void *)&templongbuf[z0],z1-z0+1,-3);
				if (typ == 'K')
				{
					oldz = -1;
					i = xyoffs[d[0]][d[1]+1] - xyoffs[d[0]][d[1]]; if (!i) continue;
					while (i > 0)
					{
						z = fgetc(fp); zleng = fgetc(fp); i -= (zleng+3);
						vis = fgetc(fp);

						if ((oldz >= 0) && (!(vis&16)))
							for(;oldz<z;oldz++)
								templongbuf[(oldz^k[7])+k[8]] = vx5.curcol;

						for(;zleng>0;zleng--,z++)
							templongbuf[(z^k[7])+k[8]] = longpal[fgetc(fp)];
						oldz = z;
					}
				}
				else
				{
					for(z=0;z<zsiz;z++)
						templongbuf[(z^k[7])+k[8]] = longpal[fgetc(fp)];
				}

				scum((d[k[0]]^k[1])+k[2],(d[k[3]]^k[4])+k[5],z0,z1+1,templongbuf);
				if (!j) scumfinish();
			}
			else
			{
				if (typ == 'K')
				{
					oldz = -1;
					i = xyoffs[d[0]][d[1]+1] - xyoffs[d[0]][d[1]]; if (!i) continue;
					while (i > 0)
					{
						z = fgetc(fp); zleng = fgetc(fp); i -= (zleng+3);
						vis = fgetc(fp);

						if ((oldz >= 0) && (!(vis&16)))
							for(;oldz<z;oldz++)
							{
								d[2] = oldz;
								setcube((d[k[0]]^k[1])+k[2],(d[k[3]]^k[4])+k[5],(d[k[6]]^k[7])+k[8],vx5.curcol);
							}

						for(;zleng>0;zleng--,z++)
						{
							ch = fgetc(fp);
							d[2] = z;
							setcube((d[k[0]]^k[1])+k[2],(d[k[3]]^k[4])+k[5],(d[k[6]]^k[7])+k[8],longpal[ch]);
						}
						oldz = z;
					}
				}
				else
				{
					for(z=0;z<zsiz;z++)
					{
						ch = fgetc(fp);
						if (ch != 255)
						{
							d[2] = z;
							setcube((d[k[0]]^k[1])+k[2],(d[k[3]]^k[4])+k[5],(d[k[6]]^k[7])+k[8],longpal[ch]);
						}
					}
				}
			}
		}
	}
	if (j) scumfinish();
	updatebbox(vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,0);

	fclose(fp);
}

//------------------------- SXL parsing code begins --------------------------


//--------------------------  Name hash code begins --------------------------

	//khashbuf format: (used by getkv6/getkfa to avoid duplicate loads)
	//[long index to next hash or -1][pointer to struct][char type]string[\0]
	//[long index to next hash or -1][pointer to struct][chat type]string[\0]
	//...
	//type:0 = kv6data
	//type:1 = kfatype
#define KHASHINITSIZE 8192
char *khashbuf = 0; // was static
long khashead[256], khashpos = 0, khashsiz = 0;

//-------------------------- KV6 sprite code begins --------------------------

//EQUIVEC code begins -----------------------------------------------------
point3d univec[256];
__declspec(align(8)) short iunivec[256][4];

typedef struct
{
	float fibx[45], fiby[45];
	float azval[20], zmulk, zaddk;
	long fib[47], aztop, npoints;
} equivectyp;
static equivectyp equivec;


void equiind2vec (long i, float *x, float *y, float *z)
{
	float r;
	(*z) = (float)i*equivec.zmulk + equivec.zaddk; r = sqrt(1.f - (*z)*(*z));
	fcossin((float)i*(GOLDRAT*PI*2),x,y); (*x) *= r; (*y) *= r;
}

	//Very fast; good quality
long equivec2indmem (float x, float y, float z)
{
	long b, i, j, k, bestc;
	float xy, zz, md, d;

	xy = atan2(y,x); //atan2 is 150 clock cycles!
	j = ((*(long *)&z)&0x7fffffff);
	bestc = equivec.aztop;
	do
	{
		if (j < *(long *)&equivec.azval[bestc]) break;
		bestc--;
	} while (bestc);

	zz = z + 1.f;
	ftol(equivec.fibx[bestc]*xy + equivec.fiby[bestc]*zz - .5,&i);
	bestc++;
	ftol(equivec.fibx[bestc]*xy + equivec.fiby[bestc]*zz - .5,&j);

	k = dmulshr0(equivec.fib[bestc+2],i,equivec.fib[bestc+1],j);
	if ((unsigned long)k < equivec.npoints)
	{
		md = univec[k].x*x + univec[k].y*y + univec[k].z*z;
		j = k;
	} else md = -2.f;
	b = bestc+3;
	do
	{
		i = equivec.fib[b] + k;
		if ((unsigned long)i < equivec.npoints)
		{
			d = univec[i].x*x + univec[i].y*y + univec[i].z*z;
			if (*(long *)&d > *(long *)&md) { md = d; j = i; }
		}
		b--;
	} while (b != bestc);
	return(j);
}

void equivecinit (long n)
{
	float t0, t1;
	long z;

		//Init constants for ind2vec
	equivec.npoints = n;
	equivec.zmulk = 2 / (float)n; equivec.zaddk = equivec.zmulk*.5 - 1.0;

		//equimemset
	for(z=n-1;z>=0;z--)
		equiind2vec(z,&univec[z].x,&univec[z].y,&univec[z].z);
	if (n&1) //Hack for when n=255 and want a <0,0,0> vector
		{ univec[n].x = univec[n].y = univec[n].z = 0; }

		//Init Fibonacci table
	equivec.fib[0] = 0; equivec.fib[1] = 1;
	for(z=2;z<47;z++) equivec.fib[z] = equivec.fib[z-2]+equivec.fib[z-1];

		//Init fibx/y LUT
	t0 = .5 / PI; t1 = (float)n * -.5;
	for(z=0;z<45;z++)
	{
		t0 = -t0; equivec.fibx[z] = (float)equivec.fib[z+2]*t0;
		t1 = -t1; equivec.fiby[z] = ((float)equivec.fib[z+2]*GOLDRAT - (float)equivec.fib[z])*t1;
	}

	t0 = 1 / ((float)n * PI);
	for(equivec.aztop=0;equivec.aztop<20;equivec.aztop++)
	{
		t1 = 1 - (float)equivec.fib[(equivec.aztop<<1)+6]*t0; if (t1 < 0) break;
		equivec.azval[equivec.aztop+1] = sqrt(t1);
	}
}

//EQUIVEC code ends -------------------------------------------------------

/** Bitmask for multiplying npix at 9 mip levels */
static const long umulmip[9] =
{
	(long)0             // 0b00000000000000000000000000000000
	,(long)4294967295   // 0b11111111111111111111111111111111
	,(long)2147483648   // 0b10000000000000000000000000000000
	,(long)1431655765   // 0b01010101010101010101010101010101
	,(long)1073741824   // 0b01000000000000000000000000000000
	,(long)858993459    // 0b00110011001100110011001100110011
	,(long)715827882    // 0b00101010101010101010101010101010
	,(long)613566756    // 0b00100100100100100100100100100100
	,(long)536870912    // 0b00100000000000000000000000000000
};

/** Generate 1 more mip-level for a .KV6 sprite.
 *  This function generates a
 *  lower MIP level only if kv6->lowermip is NULL, and kv6->xsiz,
 *  kv6->ysiz, and kv6->zsiz are all >= 3. When these conditions are
 *  true, it will generate a new .KV6 sprite with half the resolution in
 *  all 3 dimensions. It will set kv6->lowermip so it points to the newly
 *  generated .KV6 object. You can use freekv6() to de-allocate all levels
 *  of the .KV6 object.
 *
 *  To generate all mip levels use this pseudo-code:
 *  for(kv6data *tempkv6=mykv6;tempkv6=genmipkv6(tempkv6););
 *
 *  @param kv6 pointer to current MIP-level
 *  @return pointer to newly generated half-size MIP-level
 */
kv6data *genmipkv6 (kv6data *kv6)
{
	kv6data *nkv6;
	kv6voxtype *v0[2], *vs[4], *ve[4], *voxptr;
	unsigned short *xyptr, *xyi2, *sxyi2;
	long i, j, x, y, z, xs, ys, zs, xysiz, n, oxn, oxyn, *xptr;
	long xx, yy, zz, r, g, b, vis, npix, sxyi2i, darand = 0;
	char vecbuf[8];

	if ((!kv6) || (kv6->lowermip)) return(0);

	xs = ((kv6->xsiz+1)>>1); ys = ((kv6->ysiz+1)>>1); zs = ((kv6->zsiz+1)>>1);
	if ((xs < 2) || (ys < 2) || (zs < 2)) return(0);
	xysiz = ((((xs*ys)<<1)+3)&~3);
	i = sizeof(kv6data) + (xs<<2) + xysiz + kv6->numvoxs*sizeof(kv6voxtype);
	nkv6 = (kv6data *)malloc(i);
	if (!nkv6) return(0);

	kv6->lowermip = nkv6;
	nkv6->xsiz = xs;
	nkv6->ysiz = ys;
	nkv6->zsiz = zs;
	nkv6->xpiv = kv6->xpiv*.5;
	nkv6->ypiv = kv6->ypiv*.5;
	nkv6->zpiv = kv6->zpiv*.5;
	nkv6->namoff = 0;
	nkv6->lowermip = 0;

	xptr = (long *)(((long)nkv6) + sizeof(kv6data));
	xyptr = (unsigned short *)(((long)xptr) + (xs<<2));
	voxptr = (kv6voxtype *)(((long)xyptr) + xysiz);
	n = 0;

	v0[0] = kv6->vox; sxyi2 = kv6->ylen; sxyi2i = (kv6->ysiz<<1);
	for(x=0;x<xs;x++)
	{
		v0[1] = v0[0]+kv6->xlen[x<<1];

			//vs: start pointer of each of the 4 columns
			//ve: end pointer of each of the 4 columns
		vs[0] = v0[0]; vs[2] = v0[1];

		xyi2 = sxyi2; sxyi2 += sxyi2i;

		oxn = n;
		for(y=0;y<ys;y++)
		{
			oxyn = n;

			ve[0] = vs[1] = vs[0]+xyi2[0];
			if ((x<<1)+1 < kv6->xsiz) { ve[2] = vs[3] = vs[2]+xyi2[kv6->ysiz]; }
			if ((y<<1)+1 < kv6->ysiz)
			{
				ve[1] = vs[1]+xyi2[1];
				if ((x<<1)+1 < kv6->xsiz) ve[3] = vs[3]+xyi2[kv6->ysiz+1];
			}
			xyi2 += 2;

			while (1)
			{
				z = 0x7fffffff;
				for(i=3;i>=0;i--)
					if ((vs[i] < ve[i]) && (vs[i]->z < z)) z = vs[i]->z;
				if (z == 0x7fffffff) break;

				z |= 1;

				r = g = b = vis = npix = 0;
				for(i=3;i>=0;i--)
					for(zz=z-1;zz<=z;zz++)
					{
						if ((vs[i] >= ve[i]) || (vs[i]->z > zz)) continue;
						r += (vs[i]->col&0xff00ff); //MMX-style trick!
						g += (vs[i]->col&  0xff00);
						//b += (vs[i]->col&    0xff);
						vis |= vs[i]->vis;
						vecbuf[npix] = vs[i]->dir;
						npix++; vs[i]++;
					}

				if (npix)
				{
					if (n >= kv6->numvoxs) return(0); //Don't let it crash!

					i = umulmip[npix]; j = (npix>>1);
					voxptr[n].col = (umulshr32(r+(j<<16),i)&0xff0000) +
										 (umulshr32(g+(j<< 8),i)&  0xff00) +
										 (umulshr32((r&0xfff)+ j     ,i));
					voxptr[n].z = (z>>1);
					voxptr[n].vis = vis;
					voxptr[n].dir = vecbuf[umulshr32(darand,npix)]; darand += i;
					n++;
				}
			}
			xyptr[0] = n-oxyn; xyptr++;
			vs[0] = ve[1]; vs[2] = ve[3];
		}
		xptr[x] = n-oxn;
		if ((x<<1)+1 >= kv6->xsiz) break; //Avoid read page fault
		v0[0] = v0[1]+kv6->xlen[(x<<1)+1];
	}

	nkv6->leng = sizeof(kv6data) + (xs<<2) + xysiz + n*sizeof(kv6voxtype);
	nkv6 = (kv6data *)realloc(nkv6,nkv6->leng); if (!nkv6) return(0);
	nkv6->xlen = (unsigned long *)(((long)nkv6) + sizeof(kv6data));
	nkv6->ylen = (unsigned short *)(((long)nkv6->xlen) + (xs<<2));
	nkv6->vox = (kv6voxtype *)(((long)nkv6->ylen) + xysiz);
	nkv6->numvoxs = n;
	kv6->lowermip = nkv6;
	return(nkv6);
}

#ifdef __cplusplus
extern "C" {
#endif

extern void *caddasm;
#define cadd4 ((point4d *)&caddasm)
extern void *ztabasm;
#define ztab4 ((point4d *)&ztabasm)
extern short qsum0[4], qsum1[4], qbplbpp[4];
extern long kv6frameplace, kv6bytesperline;
extern float scisdist;
extern int64_t kv6colmul[256], kv6coladd[256];

char ptfaces16[43][8] =
{
	0, 0, 0,  0,  0, 0, 0,0,  4, 0,32,96, 64, 0,32,0,  4,16,80,112,48, 16,80,0,  0,0,0,0,0,0,0,0,
	4,64,96,112, 80,64,96,0,  6, 0,32,96,112,80,64,0,  6,16,80, 64,96,112,48,0,  0,0,0,0,0,0,0,0,
	4, 0,16, 48, 32, 0,16,0,  6, 0,16,48, 32,96,64,0,  6, 0,16, 80,112,48,32,0,  0,0,0,0,0,0,0,0,
	0, 0, 0,  0,  0, 0, 0,0,  0, 0, 0, 0,  0, 0, 0,0,  0, 0, 0,  0,  0, 0, 0,0,  0,0,0,0,0,0,0,0,
	4, 0,64, 80, 16, 0,64,0,  6, 0,32,96, 64,80,16,0,  6, 0,64, 80,112,48,16,0,  0,0,0,0,0,0,0,0,
	6, 0,64, 96,112,80,16,0,  6, 0,32,96,112,80,16,0,  6, 0,64, 96,112,48,16,0,  0,0,0,0,0,0,0,0,
	6, 0,64, 80, 16,48,32,0,  6,16,48,32, 96,64,80,0,  6, 0,64, 80,112,48,32,0,  0,0,0,0,0,0,0,0,
	0, 0, 0,  0,  0, 0, 0,0,  0, 0, 0, 0,  0, 0, 0,0,  0, 0, 0,  0,  0, 0, 0,0,  0,0,0,0,0,0,0,0,
	4,32,48,112, 96,32,48,0,  6, 0,32,48,112,96,64,0,  6,16,80,112, 96,32,48,0,  0,0,0,0,0,0,0,0,
	6,32,48,112, 80,64,96,0,  6, 0,32,48,112,80,64,0,  6,16,80, 64, 96,32,48,0,  0,0,0,0,0,0,0,0,
	6, 0,16, 48,112,96,32,0,  6, 0,16,48,112,96,64,0,  6, 0,16, 80,112,96,32,0,
};

void drawboundcubesseinit ();
void drawboundcubesse (kv6voxtype *, long);
void drawboundcube3dninit ();
void drawboundcube3dn (kv6voxtype *, long);

#ifdef __cplusplus
}
#endif

//static void initboundcubescr (long dafram, long dabpl, long x, long y, long dabpp)
//{
//   qsum1[3] = qsum1[1] = 0x7fff-y; qsum1[2] = qsum1[0] = 0x7fff-x;
//   qbplbpp[1] = dabpl; qbplbpp[0] = ((dabpp+7)>>3);
//   kv6frameplace = dafram; kv6bytesperline = dabpl;
//}

static __declspec(align(8)) short lightlist[MAXLIGHTS+1][4];
static int64_t all32767 = 0x7fff7fff7fff7fff;

static void updatereflects (vx5sprite *spr)
{
	__int64 fogmul;
	point3d tp;
	float f, g, h, fx, fy, fz;
	long i, j;

#if 0
	KV6 lighting calculations for: fog, white, black, intens(normal dot product), black currently not supported!

	long vx5.kv6black = 0x000000, vx5.kv6white = 0x808080;
	long nw.r = vx5.kv6white.r-vx5.kv6black.r, nb.r = vx5.kv6black.r*2;
	long nw.g = vx5.kv6white.g-vx5.kv6black.g, nb.g = vx5.kv6black.g*2;
	long nw.b = vx5.kv6white.b-vx5.kv6black.b, nb.b = vx5.kv6black.b*2;
	col.r = mulshr7(col.r,nw.r)+nb.r; col.r = mulshr7(col.r,intens); col.r += mulshr15(fogcol.r-col.r,fogmul);
	col.g = mulshr7(col.g,nw.g)+nb.g; col.g = mulshr7(col.g,intens); col.g += mulshr15(fogcol.g-col.g,fogmul);
	col.b = mulshr7(col.b,nw.b)+nb.b; col.b = mulshr7(col.b,intens); col.b += mulshr15(fogcol.b-col.b,fogmul);

	col.r = ((col.r*intens*nw.r*(32767-fogmul))>>29) +
			  ((      intens*nb.r*(32767-fogmul))>>22) + ((fogcol.r*fogmul)>>15);
	col.g = ((col.g*intens*nw.g*(32767-fogmul))>>29) +
			  ((      intens*nb.g*(32767-fogmul))>>22) + ((fogcol.g*fogmul)>>15);
	col.b = ((col.b*intens*nw.b*(32767-fogmul))>>29) +
			  ((      intens*nb.b*(32767-fogmul))>>22) + ((fogcol.b*fogmul)>>15);
#endif

		//Use cylindrical x-y distance for fog
	if (ofogdist >= 0)
	{

		ftol(sqrt((spr->p.x-gipos.x)*(spr->p.x-gipos.x) + (spr->p.y-gipos.y)*(spr->p.y-gipos.y)),&i);
		if (i > 2047) i = 2047;
		fogmul = foglut[i];

#if 0
		i = (long)(*(short *)&fogmul);
		((short *)kv6coladd)[0] = (short)((((long)(((short *)&fogcol)[0]))*i)>>1);
		((short *)kv6coladd)[1] = (short)((((long)(((short *)&fogcol)[1]))*i)>>1);
		((short *)kv6coladd)[2] = (short)((((long)(((short *)&fogcol)[2]))*i)>>1);
#else
		_asm
		{
			movq mm0, fogcol
			paddd mm0, mm0
			pmulhuw mm0, fogmul
			movq kv6coladd[0], mm0
			emms
		}
#endif
	} else { fogmul = 0I64; kv6coladd[0] = 0I64; }

	if (spr->flags&1)
	{
		//((short *)kv6colmul)[0] = 0x010001000100;                                        //Do (nothing)
		  ((short *)kv6colmul)[0] = (short)((vx5.kv6col&255)<<1);                          //Do     white
		//((short *)kv6colmul)[0] = (short)((32767-(short)fogmul)>>7);                     //Do fog
		//((short *)kv6colmul)[0] = (short)(((32767-(short)fogmul)*(vx5.kv6col&255))>>14); //Do fog&white

		((short *)kv6colmul)[1] = ((short *)kv6colmul)[0];
		((short *)kv6colmul)[2] = ((short *)kv6colmul)[0];
		((short *)kv6colmul)[3] = 0;
		for(i=1;i<256;i++) kv6colmul[i] = kv6colmul[0];
		return;
	}

	if (vx5.lightmode < 2)
	{
		fx = 1.0; fy = 1.0; fz = 1.0;
		tp.x = spr->s.x*fx + spr->s.y*fy + spr->s.z*fz;
		tp.y = spr->h.x*fx + spr->h.y*fy + spr->h.z*fz;
		tp.z = spr->f.x*fx + spr->f.y*fy + spr->f.z*fz;

		// 64 / sqrt( tp dot tp ) == 64 * f_rsqrt( tp dot tp )
		f = 64.0 * f_rsqrt( tp.x*tp.x + tp.y*tp.y + tp.z*tp.z );

			//for(i=255;i>=0;i--)
			//{
			//   ftol(univec[i].x*tp.x + univec[i].y*tp.y + univec[i].z*tp.z,&j);
			//   j = (lbound0(j+128,255)<<8);
			//   ((unsigned short *)(&kv6colmul[i]))[0] = j;
			//   ((unsigned short *)(&kv6colmul[i]))[1] = j;
			//   ((unsigned short *)(&kv6colmul[i]))[2] = j;
			//}
		g = ((float)((((long)fogmul)&32767)^32767))*(16.f*8.f/65536.f);
		if (!(((vx5.kv6col&0xffff)<<8)^(vx5.kv6col&0xffff00))) //Cool way to check if R==G==B :)
		{
			g *= ((float)(vx5.kv6col&255))/256.f;
				//This case saves 1 MMX multiply per iteration
			f *= g;
			lightlist[0][0] = (short)(tp.x*f);
			lightlist[0][1] = (short)(tp.y*f);
			lightlist[0][2] = (short)(tp.z*f);
			lightlist[0][3] = (short)(g*128.f);
			#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
			__asm__ __volatile__
			(
			".Lnolighta:\n"
				"movq	%c[uv](%[c]), %[y0]\n"
				"movq	%c[uv]-8(%[c]), %[y1]\n"
				"pmaddwd	%[y6], %[y0]\n"      //mm0: [tp.a*iunivec.a + tp.z*iunivec.z][tp.y*iunivec.y + tp.x*iunivec.x]
				"pmaddwd	%[y6], %[y1]\n"
				"pshufw	$0x4e, %[y0], %[y2]\n"   //Before: mm0: [ 0 ][ a ][   ][   ][ 0 ][ b ][   ][   ]
				"pshufw	$0x4e, %[y1], %[y3]\n"
				"paddd	%[y2], %[y0]\n"
				"paddd	%[y3], %[y1]\n"
				"pshufw $0x55, %[y0], %[y0]\n"
				"pshufw	$0x55, %[y1], %[y1]\n"   //After:  mm0: [   ][   ][   ][a+b][   ][a+b][   ][a+b]
				"movq	%[y0], %c[kvcm](%[c])\n"
				"movq	%[y1], %c[kvcm]-8(%[c])\n"
				"sub	$2*8, %[c]\n"
				"jnc    .Lnolighta\n"
				: [y0] "=y" (reg0), [y1] "=y" (reg1),
				  [y2] "=y" (reg2), [y3] "=y" (reg3),
				  [y6] "=y" (reg6)
				: [c]  "r" (255*8), "4" (*(int64_t *)lightlist),
				  [uv] "p" (iunivec), [kvcm] "p" (kv6colmul)
				:
			);
			#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
			_asm
			{
				movq	mm6, lightlist[0]
				mov	ecx, 255*8
			nolighta:
				movq	mm0, iunivec[ecx]
				movq	mm1, iunivec[ecx-8]
				pmaddwd	mm0, mm6 //mm0: [tp.a*iunivec.a + tp.z*iunivec.z][tp.y*iunivec.y + tp.x*iunivec.x]
				pmaddwd	mm1, mm6
				pshufw	mm2, mm0, 0x4e  //Before: mm0: [ 0 ][ a ][   ][   ][ 0 ][ b ][   ][   ]
				pshufw	mm3, mm1, 0x4e
				paddd	mm0, mm2
				paddd	mm1, mm3
				pshufw	mm0, mm0, 0x55
				pshufw	mm1, mm1, 0x55  //After:  mm0: [   ][   ][   ][a+b][   ][a+b][   ][a+b]
				movq	kv6colmul[ecx], mm0
				movq	kv6colmul[ecx-8], mm1
				sub	ecx, 2*8
				jnc	short nolighta
			}
			#endif 
		}
		else
		{
			f *= g;
			lightlist[0][0] = (short)(tp.x*f);
			lightlist[0][1] = (short)(tp.y*f);
			lightlist[0][2] = (short)(tp.z*f);
			lightlist[0][3] = (short)(g*128.f);
			#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
			__asm__ __volatile__
			(
				"punpcklbw	%[vxpart], %[y5]\n"
			".Lnolightb:\n"
				"movq	%c[uv](%[c]), %[y0]\n"
				"movq	%c[uv]-8(%[c]), %[y1]\n"
				"pmaddwd	%[y6], %[y0]\n"      //mm0: [tp.a*iunivec.a + tp.z*iunivec.z][tp.y*iunivec.y + tp.x*iunivec.x]
				"pmaddwd	%[y6], %[y1]\n"
				"pshufw	$0x4e, %[y0], %[y2]\n"   //Before: mm0: [ 0 ][ a ][   ][   ][ 0 ][ b ][   ][   ]
				"pshufw	$0x4e, %[y1], %[y3]\n"
				"paddd	%[y2], %[y0]\n"
				"paddd	%[y3], %[y1]\n"
				"pshufw $0x55, %[y0], %[y0]\n"
				"pshufw	$0x55, %[y1], %[y1]\n"   //After:  mm0: [   ][   ][   ][a+b][   ][a+b][   ][a+b]
				"pmulhuw	%[y5], %[y0]\n"
				"pmulhuw	%[y5], %[y1]\n"
				"movq	%[y0], %c[kvcm](%[c])\n"
				"movq	%[y1], %c[kvcm]-8(%[c])\n"
				"sub	$2*8, %[c]\n"
				"jnc .Lnolightb\n"
				: [y0] "+y" (reg0), [y1] "+y" (reg1),
				  [y2] "+y" (reg2), [y3] "+y" (reg3),
				  [y5] "=y" (reg5), [y6] "=y" (reg6)
				: "5" (*(int64_t *)lightlist),
				  [c]  "r" (255*8), [vxpart] "m" (vx5.kv6col),
				  [uv] "p" (iunivec), [kvcm] "p" (kv6colmul)
			);
            #elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
			_asm
			{
				punpcklbw	mm5, vx5.kv6col
				movq	mm6, lightlist[0]
				mov	ecx, 255*8
			nolightb:
				movq	mm0, iunivec[ecx]
				movq	mm1, iunivec[ecx-8]
				pmaddwd	mm0, mm6 //mm0: [tp.a*iunivec.a + tp.z*iunivec.z][tp.y*iunivec.y + tp.x*iunivec.x]
				pmaddwd	mm1, mm6
				pshufw	mm2, mm0, 0x4e //Before: mm0: [ 0 ][ a ][   ][   ][ 0 ][ b ][   ][   ]
				pshufw	mm3, mm1, 0x4e
				paddd	mm0, mm2
				paddd	mm1, mm3
				pshufw	mm0, mm0, 0x55
				pshufw	mm1, mm1, 0x55 //After:  mm0: [   ][   ][   ][a+b][   ][a+b][   ][a+b]
				pmulhuw	mm0, mm5
				pmulhuw	mm1, mm5
				movq	kv6colmul[ecx], mm0
				movq	kv6colmul[ecx-8], mm1
				sub	ecx, 2*8
				jnc short nolightb
			}
			#endif
		}
		//NOTE: emms not necessary!
	}
	else
	{
		point3d sprs, sprh, sprf;
		float ff, gg, hh;
		long k, lightcnt;

			//WARNING: this only works properly for orthonormal matrices!
		f = f_rsqrt(spr->s.x*spr->s.x + spr->s.y*spr->s.y + spr->s.z*spr->s.z);
		sprs.x = spr->s.x*f; sprs.y = spr->s.y*f; sprs.z = spr->s.z*f;
		f = f_rsqrt(spr->h.x*spr->h.x + spr->h.y*spr->h.y + spr->h.z*spr->h.z);
		sprh.x = spr->h.x*f; sprh.y = spr->h.y*f; sprh.z = spr->h.z*f;
		f = f_rsqrt(spr->f.x*spr->f.x + spr->f.y*spr->f.y + spr->f.z*spr->f.z);
		sprf.x = spr->f.x*f; sprf.y = spr->f.y*f; sprf.z = spr->f.z*f;

		hh = ((float)((((long)fogmul)&32767)^32767))/65536.f * 2.f;


			//Find which lights are close enough to affect sprite.
		lightcnt = 0;
		for(i=vx5.numlights-1;i>=0;i--)
		{
			fx = vx5.lightsrc[i].p.x-(spr->p.x);
			fy = vx5.lightsrc[i].p.y-(spr->p.y);
			fz = vx5.lightsrc[i].p.z-(spr->p.z);
			gg = fx*fx + fy*fy + fz*fz; ff = vx5.lightsrc[i].r2;
			if (*(long *)&gg < *(long *)&ff)
			{
				f = sqrt(ff); g = sqrt(gg);
				//h = (16.0/(sqrt(gg)*gg) - 16.0/(sqrt(ff)*ff))*vx5.lightsrc[i].sc;
				h = (f*ff - g*gg)/(f*ff*g*gg) * vx5.lightsrc[i].sc*16.0;
				if (g*h > 4096.0) h = 4096.0/g; //Max saturation clipping
				h *= hh;
				lightlist[lightcnt][0] = (short)((fx*sprs.x + fy*sprs.y + fz*sprs.z)*h);
				lightlist[lightcnt][1] = (short)((fx*sprh.x + fy*sprh.y + fz*sprh.z)*h);
				lightlist[lightcnt][2] = (short)((fx*sprf.x + fy*sprf.y + fz*sprf.z)*h);
				lightlist[lightcnt][3] = 0;
				lightcnt++;
			}
		}

		fx = 0.0; fy = 0.5; fz = 1.0;

		#if !defined(NOASM)
		hh *= 16*16.f*8.f/2.f;
		lightlist[lightcnt][0] = (short)((sprs.x*fx + sprs.y*fy + sprs.z*fz)*hh);
		lightlist[lightcnt][1] = (short)((sprh.x*fx + sprh.y*fy + sprh.z*fz)*hh);
		lightlist[lightcnt][2] = (short)((sprf.x*fx + sprf.y*fy + sprf.z*fz)*hh);
		lightlist[lightcnt][3] = (short)(hh*(48/16.0));
        #if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
		__asm__ __volatile__
		(
			"punpcklbw	%[vxpart], %[y5]\n"
			"pxor	%[y6], %[y6]\n"
			"shl	$3, %[d]\n"
		".Lbeglig:\n"
			"movq	%c[uv](%[c]), %[y3]\n" //mm3: 256 u[i].z*256 u[i].y*256 u[i].x*256
			"mov	%[d], %[a]\n"
			"movq	%c[ll](%[d]), %[y0]\n" //mm0: 48*256,0 tp.z*256 tp.y*256 tp.x*256
			"pmaddwd	%[y3], %[y0]\n"
			"pshufw	$0x4e, %[y0], %[y2]\n"
			"paddd	%[y2], %[y0]\n"
			"sub	$8, %[a]\n"
			"js	.Lendlig\n"
		".Lbeglig2:\n"
			"movq	%c[ll](%[a]), %[y1]\n" //mm1: 0 tp.z*256 tp.y*256 tp.x*256
			"pmaddwd	%[y3], %[y1]\n"
			"pshufw	$0x4e, %[y1], %[y2]\n"
			"paddd	%[y2], %[y1]\n"
			"pminsw	%[y6], %[y1]\n"        //16-bits is ugly, but ok here
			"psubd	%[y1], %[y0]\n"
			"sub	$8, %[a]\n"
			"jns	.Lbeglig2\n"           //mm0: 00 II ii ii 00 II ii ii
		".Lendlig:\n"
			"pshufw	$0x55, %[y0], %[y0]\n"  //mm0: 00 II 00 II 00 II 00 II
			"pmulhuw	%[y5], %[y0]\n"
			"movq	%[y0], %c[kvcm](%[c])\n"
			"sub	$8, %[c]\n"
			"jnc	.Lbeglig\n"
			:
			: [y0] "y" (reg0), [y1] "y" (reg1),
			  [y2] "y" (reg2), [y3] "y" (reg3),
			  [y5] "y" (reg5), [y6] "y" (reg6),
			  [c]  "r" (255*8), [d] "r" (lightcnt), [a] "r" (0),
			  [vxpart] "m" (vx5.kv6col),
			  [uv] "p" (iunivec), [kvcm] "p" (kv6colmul),
			  [ll] "p" (lightlist)
			:
		);
        #elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
		_asm
		{
			punpcklbw	mm5, vx5.kv6col
			pxor	mm6, mm6
			mov	edx, lightcnt
			shl	edx, 3
			mov	ecx, 255*8
		beglig:
			movq	mm3, iunivec[ecx]   //mm3: 256 u[i].z*256 u[i].y*256 u[i].x*256
			mov	eax, edx
			movq	mm0, lightlist[edx] //mm0: 48*256,0 tp.z*256 tp.y*256 tp.x*256
			pmaddwd	mm0, mm3
			pshufw	mm2, mm0, 0x4e
			paddd	mm0, mm2
			sub	eax, 8
			js	short endlig
		beglig2:
			movq	mm1, lightlist[eax] //mm1: 0 tp.z*256 tp.y*256 tp.x*256
			pmaddwd	mm1, mm3
			pshufw	mm2, mm1, 0x4e
			paddd	mm1, mm2
			pminsw	mm1, mm6         //16-bits is ugly, but ok here
			psubd	mm0, mm1
			sub	eax, 8
			jns	short beglig2        //mm0: 00 II ii ii 00 II ii ii
		endlig:
			pshufw	mm0, mm0, 0x55   //mm0: 00 II 00 II 00 II 00 II
			pmulhuw	mm0, mm5
			movq	kv6colmul[ecx], mm0
			sub	ecx, 8
			jnc	short beglig
		}
		#endif
		//NOTE: emms not necessary!

		#else // C Default
		tp.x = (sprs.x*fx + sprs.y*fy + sprs.z*fz)*16.0;
		tp.y = (sprh.x*fx + sprh.y*fy + sprh.z*fz)*16.0;
		tp.z = (sprf.x*fx + sprf.y*fy + sprf.z*fz)*16.0;
		for(i=255;i>=0;i--)
		{
			f = tp.x*univec[i].x + tp.y*univec[i].y + tp.z*univec[i].z + 48;
			for(k=lightcnt-1;k>=0;k--)
			{
				h = lightlist[k][0]*univec[i].x + lightlist[k][1]*univec[i].y + lightlist[k][2]*univec[i].z;
				if (*(long *)&h < 0) f -= h;
			}
			if (f > 255) f = 255;
			ftol(f,&j); j <<= 8;
			((unsigned short *)(&kv6colmul[i]))[0] = j;
			((unsigned short *)(&kv6colmul[i]))[1] = j;
			((unsigned short *)(&kv6colmul[i]))[2] = j;
		}
		#endif
	}
}


#define DRAWBOUNDCUBELINE(const) \
	for(;v0<=v1 && v0->z<inz;v0++) drawboundcubesse(v0,const+0x20);\
	for(;v0<=v1 && v1->z>inz;v1--) drawboundcubesse(v1,const+0x10);\
						  if (v0 == v1) drawboundcubesse(v1,const+0x00);

#define DRAWBOUNDCUBELINE_3DN(const) \
	for(;v0<=v1 && v0->z<inz;v0++) drawboundcube3dn(v0,const+0x20);\
	for(;v0<=v1 && v1->z>inz;v1--) drawboundcube3dn(v1,const+0x10);\
						  if (v0 == v1) drawboundcube3dn(v1,const+0x00);

	//Code taken from renderboundcube of SLAB6D (Pentium III version :)
#define MAXZSIZ 1024
static void kv6draw (vx5sprite *spr)
{
	point4d *r0, *r1, *r2;
	kv6voxtype *xv, *yv, *v0, *v1;
	kv6data *kv;
	point3d ts, th, tf;
	point3d npos, nstr, nhei, nfor, tp, tp2;
	float f;
	long x, y, z, inx, iny, inz, nxplanemin, nxplanemax;
	unsigned short *ylenptr;

	kv = spr->voxnum; if (!kv) return;

	z = 0; //Quick & dirty estimation of distance
	ftol((spr->p.x-gipos.x)*gifor.x + (spr->p.y-gipos.y)*gifor.y + (spr->p.z-gipos.z)*gifor.z,&y);
	while ((kv->lowermip) && (y >= vx5.kv6mipfactor)) { kv = kv->lowermip; z++; y >>= 1; }
	if (!z)
	{
		nxplanemin = vx5.xplanemin; nxplanemax = vx5.xplanemax;
		ts = spr->s; th = spr->h; tf = spr->f;
	}
	else
	{
		nxplanemin = (vx5.xplanemin>>z);
		nxplanemax = (vx5.xplanemax>>z); f = (float)(1<<z);
		ts.x = spr->s.x*f; ts.y = spr->s.y*f; ts.z = spr->s.z*f;
		th.x = spr->h.x*f; th.y = spr->h.y*f; th.z = spr->h.z*f;
		tf.x = spr->f.x*f; tf.y = spr->f.y*f; tf.z = spr->f.z*f;
	}

		//View frustrum culling (72*,63+,12fabs,4cmp)
	tp2.x = ((float)kv->xsiz)*.5; tp.x = tp2.x - kv->xpiv;
	tp2.y = ((float)kv->ysiz)*.5; tp.y = tp2.y - kv->ypiv;
	tp2.z = ((float)kv->zsiz)*.5; tp.z = tp2.z - kv->zpiv;
	npos.x = tp.x*ts.x + tp.y*th.x + tp.z*tf.x + (spr->p.x-gipos.x);
	npos.y = tp.x*ts.y + tp.y*th.y + tp.z*tf.y + (spr->p.y-gipos.y);
	npos.z = tp.x*ts.z + tp.y*th.z + tp.z*tf.z + (spr->p.z-gipos.z);
	nstr.x = ts.x*tp2.x; nstr.y = ts.y*tp2.x; nstr.z = ts.z*tp2.x;
	nhei.x = th.x*tp2.y; nhei.y = th.y*tp2.y; nhei.z = th.z*tp2.y;
	nfor.x = tf.x*tp2.z; nfor.y = tf.y*tp2.z; nfor.z = tf.z*tp2.z;
	for(z=3;z>=0;z--) //72*,63+
	{
			//movaps xmm0, nx4      mulaps xmm0, ginor[0].x(dup 4)
			//movaps xmm1, ny4      mulaps xmm1, ginor[0].y(dup 4)
			//movaps xmm2, nz4      mulaps xmm2, ginor[0].z(dup 4)
			//addps xmm0, xmm1      addps xmm0, xmm2
			//andps xmm0, [0x7fffffff7fffffff7fffffffffffffff]
			//movhlps xmm1, xmm0    addps xmm0, xmm1
			//shufps xmm1, xmm0, 1  addss xmm0, xmm1
			//ucomiss xmm0, [0x0]   jnz retfunc
		if (fabs(nstr.x*ginor[z].x + nstr.y*ginor[z].y + nstr.z*ginor[z].z) +
			 fabs(nhei.x*ginor[z].x + nhei.y*ginor[z].y + nhei.z*ginor[z].z) +
			 fabs(nfor.x*ginor[z].x + nfor.y*ginor[z].y + nfor.z*ginor[z].z) +
					npos.x*ginor[z].x + npos.y*ginor[z].y + npos.z*ginor[z].z < 0) return;
	}
#if 0   //There are bugs when some vertices are behind ifor plane
	x = xres; y = xres; inx = 0; iny = 0;
	for(z=7;z>=0;z--) //This is useful for debugging.
	{
		tp.x = -kv->xpiv; if (z&1) tp.x += (float)kv->xsiz;
		tp.y = -kv->ypiv; if (z&2) tp.y += (float)kv->ysiz;
		tp.z = -kv->zpiv; if (z&4) tp.z += (float)kv->zsiz;
		drawspherefill(tp.x*ts.x+tp.y*th.x+tp.z*tf.x+spr->p.x,
							tp.x*ts.y+tp.y*th.y+tp.z*tf.y+spr->p.y,
							tp.x*ts.z+tp.y*th.z+tp.z*tf.z+spr->p.z,.5,0xf08040);
		tp2.x = tp.x*ts.x+tp.y*th.x+tp.z*tf.x+spr->p.x-gipos.x;
		tp2.y = tp.x*ts.y+tp.y*th.y+tp.z*tf.y+spr->p.y-gipos.y;
		tp2.z = tp.x*ts.z+tp.y*th.z+tp.z*tf.z+spr->p.z-gipos.z;
		tp.z = tp2.x*gifor.x + tp2.y*gifor.y + tp2.z*gifor.z; if (tp.z < 2) continue;
		tp.x = tp2.x*gistr.x + tp2.y*gistr.y + tp2.z*gistr.z;
		tp.y = tp2.x*gihei.x + tp2.y*gihei.y + tp2.z*gihei.z;
		ftol(tp.x*gihz/tp.z + gihx,&inz);
		if (inz < x) x = inz;
		if (inz > inx) inx = inz;
		ftol(tp.y*gihz/tp.z+gihy,&inz);
		if (inz < y) y = inz;
		if (inz > iny) iny = inz;
	}
	if (x < 0) x = 0;
	if (inx > xres) inx = xres;
	if (y < 0) y = 0;
	if (iny > yres) iny = yres;
	if (x < inx)
		for(inz=y;inz<iny;inz++)
			clearbuf((void *)(ylookup[inz]+(x<<2)+frameplace),inx-x,0L);
#endif

#if (USEZBUFFER == 0)
	lpoint3d lp;
	if (!cansee(&gipos,&spr->p,&lp)) return; //Very crappy Z-buffer!
#endif

	r0 = &ztab4[MAXZSIZ]; r1 = &ztab4[MAXZSIZ+1]; r2 = &ztab4[MAXZSIZ+2];

		//Rotate sprite from world to screen coordinates:
	mat2(&gixs,&giys,&gizs,&giadd, &ts,&th,&tf,&spr->p, &nstr,&nhei,&nfor,&npos);
	npos.x -= (kv->xpiv*nstr.x + kv->ypiv*nhei.x + kv->zpiv*nfor.x);
	npos.y -= (kv->xpiv*nstr.y + kv->ypiv*nhei.y + kv->zpiv*nfor.y);
	npos.z -= (kv->xpiv*nstr.z + kv->ypiv*nhei.z + kv->zpiv*nfor.z);

		//Find split point by using Cramer's rule
		//Must use Cramer's rule for non-orthonormal input matrices
	tp.x = nhei.y*nfor.z - nfor.y*nhei.z;
	tp.y = nfor.y*nstr.z - nstr.y*nfor.z;
	tp.z = nstr.y*nhei.z - nhei.y*nstr.z;
	f = nstr.x*tp.x + nhei.x*tp.y + nfor.x*tp.z;
	if (f != 0)
	{
		f = -1.0f / f;
		tp2.x = npos.y*nfor.z - nfor.y*npos.z;
		tp2.y = nhei.y*npos.z - npos.y*nhei.z;
		tp2.z = npos.y*nstr.z - nstr.y*npos.z;
		inx = (npos.x*tp.x - nhei.x*tp2.x - nfor.x*tp2.y)*f;
		iny = (npos.x*tp.y + nstr.x*tp2.x - nfor.x*tp2.z)*f;
		inz = (npos.x*tp.z + nstr.x*tp2.y + nhei.x*tp2.z)*f;
	}
	else { inx = iny = inz = -1; }
	inx = lbound(inx,-1,kv->xsiz);
	iny = lbound(iny,-1,kv->ysiz);
	inz = lbound(inz,-1,kv->zsiz);

	f = nhei.x; nhei.x = nfor.x; nfor.x = -f;
	f = nhei.y; nhei.y = nfor.y; nfor.y = -f;
	f = nhei.z; nhei.z = nfor.z; nfor.z = -f;

	if (kv->zsiz >= MAXZSIZ) return; //HACK TO PREVENT CRASHES FOR NOW... FIX!
	qsum0[2] = qsum0[0] = 0x7fff-(xres-(long)gihx);
	qsum0[3] = qsum0[1] = 0x7fff-(yres-(long)gihy);

		//r1->x = nstr.z; r1->y = nhei.z; r1->z = nfor.z;
		//minps(r1,r1,&ztab4[0]); //&ztab4[0] always 0
		//scisdist = -(r1->x + r1->y + r1->z);
	scisdist = 0;
	if (*(long *)&nstr.z < 0) scisdist -= nstr.z;
	if (*(long *)&nhei.z < 0) scisdist -= nhei.z;
	if (*(long *)&nfor.z < 0) scisdist -= nfor.z;

	cadd4[1].x = nstr.x*gihz; cadd4[1].y = nstr.y*gihz; cadd4[1].z = cadd4[1].z2 = nstr.z;
	cadd4[2].x = nhei.x*gihz; cadd4[2].y = nhei.y*gihz; cadd4[2].z = cadd4[2].z2 = nhei.z;
	cadd4[4].x = nfor.x*gihz; cadd4[4].y = nfor.y*gihz; cadd4[4].z = cadd4[4].z2 = nfor.z;
		  r1->x = npos.x*gihz;      r1->y = npos.y*gihz;      r1->z =      r1->z2 = npos.z;

	updatereflects(spr);
	//No more 8087 code after here!!! ----------------------------------------

	if (cputype&(1<<25))
	{
		addps(&cadd4[3],&cadd4[1],&cadd4[2]);
		addps(&cadd4[5],&cadd4[1],&cadd4[4]);
		addps(&cadd4[6],&cadd4[2],&cadd4[4]);
		addps(&cadd4[7],&cadd4[3],&cadd4[4]);

		for(z=1;z<kv->zsiz;z++) addps(&ztab4[z],&ztab4[z-1],&cadd4[2]);
		intss(r2,-kv->ysiz); mulps(r2,r2,&cadd4[4]);

		subps(r1,r1,&cadd4[4]); //ANNOYING HACK!!!

		_asm
		{
			movq mm6, qsum0
			movq mm7, qsum1
		}

		xv = kv->vox; ylenptr = kv->ylen;
		for(x=0;x<inx;x++,ylenptr+=kv->ysiz)
		{
			if ((x < nxplanemin) || (x >= nxplanemax))
				{ xv += kv->xlen[x]; addps(r1,r1,&cadd4[1]); continue; }
			yv = xv+kv->xlen[x]; movps(r0,r1);
			for(y=0;y<iny;y++)
			{
				v0 = xv; xv += ylenptr[y]; v1 = xv-1;
				DRAWBOUNDCUBELINE(0xa)
				subps(r0,r0,&cadd4[4]);
			}
			xv = yv;
			addps(r0,r1,r2);
			addps(r1,r1,&cadd4[1]);
			for(y=kv->ysiz-1;y>iny;y--)
			{
				addps(r0,r0,&cadd4[4]);
				v1 = yv-1; yv -= ylenptr[y]; v0 = yv;
				DRAWBOUNDCUBELINE(0x6)
			}
			if ((unsigned long)iny < (unsigned long)kv->ysiz)
			{
				addps(r0,r0,&cadd4[4]);
				v1 = yv-1; yv -= ylenptr[y]; v0 = yv;
				DRAWBOUNDCUBELINE(0x2)
			}
		}
		xv = &kv->vox[kv->numvoxs]; ylenptr = &kv->ylen[(kv->xsiz-1)*kv->ysiz];
		intss(r0,kv->xsiz-x); mulps(r0,r0,&cadd4[1]); addps(r1,r1,r0);
		for(x=kv->xsiz-1;x>inx;x--,ylenptr-=kv->ysiz)
		{
			if ((x < nxplanemin) || (x >= nxplanemax))
				{ xv -= kv->xlen[x]; subps(r1,r1,&cadd4[1]); continue; }
			yv = xv-kv->xlen[x];
			subps(r1,r1,&cadd4[1]);
			addps(r0,r1,r2);
			for(y=kv->ysiz-1;y>iny;y--)
			{
				addps(r0,r0,&cadd4[4]);
				v1 = xv-1; xv -= ylenptr[y]; v0 = xv;
				DRAWBOUNDCUBELINE(0x5)
			}
			xv = yv; movps(r0,r1);
			for(y=0;y<iny;y++)
			{
				v0 = yv; yv += ylenptr[y]; v1 = yv-1;
				DRAWBOUNDCUBELINE(0x9)
				subps(r0,r0,&cadd4[4]);
			}
			if ((unsigned long)iny < (unsigned long)kv->ysiz)
			{
				v0 = yv; yv += ylenptr[y]; v1 = yv-1;
				DRAWBOUNDCUBELINE(0x1)
			}
		}
		if ((unsigned long)inx < (unsigned long)kv->xsiz)
		{
			if ((x < nxplanemin) || (x >= nxplanemax)) { { clearMMX(); } return; }
			yv = xv-kv->xlen[x];
			subps(r1,r1,&cadd4[1]);
			addps(r0,r1,r2);
			for(y=kv->ysiz-1;y>iny;y--)
			{
				addps(r0,r0,&cadd4[4]);
				v1 = xv-1; xv -= ylenptr[y]; v0 = xv;
				DRAWBOUNDCUBELINE(0x4)
			}
			xv = yv; movps(r0,r1);
			for(y=0;y<iny;y++)
			{
				v0 = yv; yv += ylenptr[y]; v1 = yv-1;
				DRAWBOUNDCUBELINE(0x8)
				subps(r0,r0,&cadd4[4]);
			}
			if ((unsigned long)iny < (unsigned long)kv->ysiz)
			{
				v0 = yv; yv += ylenptr[y]; v1 = yv-1;
				DRAWBOUNDCUBELINE(0x0)
			}
		}
	}
	else
	{
		addps_3dn(&cadd4[3],&cadd4[1],&cadd4[2]);
		addps_3dn(&cadd4[5],&cadd4[1],&cadd4[4]);
		addps_3dn(&cadd4[6],&cadd4[2],&cadd4[4]);
		addps_3dn(&cadd4[7],&cadd4[3],&cadd4[4]);

		for(z=1;z<kv->zsiz;z++) addps_3dn(&ztab4[z],&ztab4[z-1],&cadd4[2]);
		intss_3dn(r2,-kv->ysiz); mulps_3dn(r2,r2,&cadd4[4]);

		subps_3dn(r1,r1,&cadd4[4]); //ANNOYING HACK!!!

		_asm
		{
			movq mm6, qsum0
			movq mm7, qsum1
		}

		xv = kv->vox; ylenptr = kv->ylen;
		for(x=0;x<inx;x++,ylenptr+=kv->ysiz)
		{
			if ((x < nxplanemin) || (x >= nxplanemax))
				{ xv += kv->xlen[x]; addps_3dn(r1,r1,&cadd4[1]); continue; }
			yv = xv+kv->xlen[x]; movps_3dn(r0,r1);
			for(y=0;y<iny;y++)
			{
				v0 = xv; xv += ylenptr[y]; v1 = xv-1;
				DRAWBOUNDCUBELINE_3DN(0xa)
				subps_3dn(r0,r0,&cadd4[4]);
			}
			xv = yv;
			addps_3dn(r0,r1,r2);
			addps_3dn(r1,r1,&cadd4[1]);
			for(y=kv->ysiz-1;y>iny;y--)
			{
				addps_3dn(r0,r0,&cadd4[4]);
				v1 = yv-1; yv -= ylenptr[y]; v0 = yv;
				DRAWBOUNDCUBELINE_3DN(0x6)
			}
			if ((unsigned long)iny < (unsigned long)kv->ysiz)
			{
				addps_3dn(r0,r0,&cadd4[4]);
				v1 = yv-1; yv -= ylenptr[y]; v0 = yv;
				DRAWBOUNDCUBELINE_3DN(0x2)
			}
		}
		xv = &kv->vox[kv->numvoxs]; ylenptr = &kv->ylen[(kv->xsiz-1)*kv->ysiz];
		intss_3dn(r0,kv->xsiz-x); mulps_3dn(r0,r0,&cadd4[1]); addps_3dn(r1,r1,r0);
		for(x=kv->xsiz-1;x>inx;x--,ylenptr-=kv->ysiz)
		{
			if ((x < nxplanemin) || (x >= nxplanemax))
				{ xv -= kv->xlen[x]; subps_3dn(r1,r1,&cadd4[1]); continue; }
			yv = xv-kv->xlen[x];
			subps_3dn(r1,r1,&cadd4[1]);
			addps_3dn(r0,r1,r2);
			for(y=kv->ysiz-1;y>iny;y--)
			{
				addps_3dn(r0,r0,&cadd4[4]);
				v1 = xv-1; xv -= ylenptr[y]; v0 = xv;
				DRAWBOUNDCUBELINE_3DN(0x5)
			}
			xv = yv; movps_3dn(r0,r1);
			for(y=0;y<iny;y++)
			{
				v0 = yv; yv += ylenptr[y]; v1 = yv-1;
				DRAWBOUNDCUBELINE_3DN(0x9)
				subps_3dn(r0,r0,&cadd4[4]);
			}
			if ((unsigned long)iny < (unsigned long)kv->ysiz)
			{
				v0 = yv; yv += ylenptr[y]; v1 = yv-1;
				DRAWBOUNDCUBELINE_3DN(0x1)
			}
		}
		if ((unsigned long)inx < (unsigned long)kv->xsiz)
		{
			if ((x < nxplanemin) || (x >= nxplanemax)) { { clearMMX(); } return; }
			yv = xv-kv->xlen[x];
			subps_3dn(r1,r1,&cadd4[1]);
			addps_3dn(r0,r1,r2);
			for(y=kv->ysiz-1;y>iny;y--)
			{
				addps_3dn(r0,r0,&cadd4[4]);
				v1 = xv-1; xv -= ylenptr[y]; v0 = xv;
				DRAWBOUNDCUBELINE_3DN(0x4)
			}
			xv = yv; movps_3dn(r0,r1);
			for(y=0;y<iny;y++)
			{
				v0 = yv; yv += ylenptr[y]; v1 = yv-1;
				DRAWBOUNDCUBELINE_3DN(0x8)
				subps_3dn(r0,r0,&cadd4[4]);
			}
			if ((unsigned long)iny < (unsigned long)kv->ysiz)
			{
				v0 = yv; yv += ylenptr[y]; v1 = yv-1;
				DRAWBOUNDCUBELINE_3DN(0x0)
			}
		}
	}
	clearMMX();
}

#endif

//-------------------------- KFA sprite code begins --------------------------

static kv6voxtype *getvptr (kv6data *kv, long x, long y)
{
	kv6voxtype *v;
	long i, j;

	v = kv->vox;
	if ((x<<1) < kv->xsiz) { for(i=0         ;i< x;i++) v += kv->xlen[i]; }
	else { v += kv->numvoxs; for(i=kv->xsiz-1;i>=x;i--) v -= kv->xlen[i]; }
	j = x*kv->ysiz;
	if ((y<<1) < kv->ysiz) { for(i=0         ;i< y;i++) v += kv->ylen[j+i]; }
	else { v += kv->xlen[x]; for(i=kv->ysiz-1;i>=y;i--) v -= kv->ylen[j+i]; }
	return(v);
}

#define VFIFSIZ 16384 //SHOULDN'T BE STATIC ALLOCATION!!!
static long vfifo[VFIFSIZ];

void floodsucksprite (vx5sprite *spr, kv6data *kv, long ox, long oy,
									  kv6voxtype *v0, kv6voxtype *v1)
{
	kv6voxtype *v, *ve, *ov, *v2, *v3;
	kv6data *kv6;
	long i, j, x, y, z, x0, y0, z0, x1, y1, z1, n, vfif0, vfif1;

	x0 = x1 = ox; y0 = y1 = oy; z0 = v0->z; z1 = v1->z;

	n = (((long)v1)-((long)v0))/sizeof(kv6voxtype)+1;
	v1->vis &= ~64;

	vfifo[0] = ox; vfifo[1] = oy;
	vfifo[2] = (long)v0; vfifo[3] = (long)v1;
	vfif0 = 0; vfif1 = 4;

	while (vfif0 < vfif1)
	{
		i = (vfif0&(VFIFSIZ-1)); vfif0 += 4;
		ox = vfifo[i]; oy = vfifo[i+1];
		v0 = (kv6voxtype *)vfifo[i+2]; v1 = (kv6voxtype *)vfifo[i+3];

		if (ox < x0) x0 = ox;
		if (ox > x1) x1 = ox;
		if (oy < y0) y0 = oy;
		if (oy > y1) y1 = oy;
		if (v0->z < z0) z0 = v0->z;
		if (v1->z > z1) z1 = v1->z;
		for(v=v0;v<=v1;v++) v->vis |= 128; //Mark as part of current piece

		for(j=0;j<4;j++)
		{
			switch(j)
			{
				case 0: x = ox-1; y = oy; break;
				case 1: x = ox+1; y = oy; break;
				case 2: x = ox; y = oy-1; break;
				case 3: x = ox; y = oy+1; break;
				default: _gtfo(); //tells MSVC default can't be reached
			}
			if ((unsigned long)x >= kv->xsiz) continue;
			if ((unsigned long)y >= kv->ysiz) continue;

			v = getvptr(kv,x,y);
			for(ve=&v[kv->ylen[x*kv->ysiz+y]];v<ve;v++)
			{
				if (v->vis&16) ov = v;
				if (((v->vis&(64+32)) == 64+32) && (v0->z <= v->z) && (v1->z >= ov->z))
				{
					i = (vfif1&(VFIFSIZ-1)); vfif1 += 4;
					if (vfif1-vfif0 >= VFIFSIZ) //FIFO Overflow... make entire object 1 piece :/
					{
						for(i=kv->numvoxs-1;i>=0;i--)
						{
							if ((kv->vox[i].vis&(64+32)) == 64+32) { v1 = &kv->vox[i]; v1->vis &= ~64; }
							if (kv->vox[i].vis&16) for(v=&kv->vox[i];v<=v1;v++) kv->vox[i].vis |= 128;
						}
						x0 = y0 = z0 = 0; x1 = kv->xsiz; y1 = kv->ysiz; z1 = kv->zsiz; n = kv->numvoxs;
						goto floodsuckend;
					}
					vfifo[i] = x; vfifo[i+1] = y;
					vfifo[i+2] = (long)ov; vfifo[i+3] = (long)v;
					n += (((long)v)-((long)ov))/sizeof(kv6voxtype)+1;
					v->vis &= ~64;
				}
			}
		}
	}
	x1++; y1++; z1++;
floodsuckend:;

	i = sizeof(kv6data) + n*sizeof(kv6voxtype) + (x1-x0)*4 + (x1-x0)*(y1-y0)*2;
	if (!(kv6 = (kv6data *)malloc(i))) return;
	kv6->leng = i;
	kv6->xsiz = x1-x0;
	kv6->ysiz = y1-y0;
	kv6->zsiz = z1-z0;
	kv6->xpiv = 0; //Set limb pivots to 0 - don't need it twice!
	kv6->ypiv = 0;
	kv6->zpiv = 0;
	kv6->numvoxs = n;
	kv6->namoff = 0;
	kv6->lowermip = 0;
	kv6->vox = (kv6voxtype *)(((long)kv6)+sizeof(kv6data));
	kv6->xlen = (unsigned long *)(((long)kv6->vox)+n*sizeof(kv6voxtype));
	kv6->ylen = (unsigned short *)(((long)kv6->xlen)+(x1-x0)*4);

		//Extract sub-KV6 to newly allocated kv6data
	v3 = kv6->vox; n = 0;
	for(x=0,v=kv->vox;x<x0;x++) v += kv->xlen[x];
	for(;x<x1;x++)
	{
		v2 = v; ox = n;
		for(y=0;y<y0;y++) v += kv->ylen[x*kv->ysiz+y];
		for(;y<y1;y++)
		{
			oy = n;
			for(ve=&v[kv->ylen[x*kv->ysiz+y]];v<ve;v++)
				if (v->vis&128)
				{
					v->vis &= ~128;
					(*v3) = (*v);
					v3->z -= z0;
					v3++; n++;
				}
			kv6->ylen[(x-x0)*(y1-y0)+(y-y0)] = n-oy;
		}
		kv6->xlen[x-x0] = n-ox;
		v = v2+kv->xlen[x];
	}

	spr->p.x = x0-kv->xpiv;
	spr->p.y = y0-kv->ypiv;
	spr->p.z = z0-kv->zpiv;
	spr->s.x = 1; spr->s.y = 0; spr->s.z = 0;
	spr->h.x = 0; spr->h.y = 1; spr->h.z = 0;
	spr->f.x = 0; spr->f.y = 0; spr->f.z = 1;
	spr->voxnum = kv6;
	spr->flags = 0;
}

void kfasorthinge (hingetype *h, long nh, long *hsort)
{
	long i, j, n;

		//First pass: stick hinges with parent=-1 at end
	n = nh; j = 0;
	for(i=n-1;i>=0;i--)
	{
		if (h[i].parent < 0) hsort[--n] = i;
							 else hsort[j++] = i;
	}
		//Finish accumulation (n*log(n) if tree is perfectly balanced)
	while (n > 0)
	{
		i--; if (i < 0) i = n-1;
		j = hsort[i];
		if (h[h[j].parent].parent < 0)
		{
			h[j].parent = -2-h[j].parent; n--;
			hsort[i] = hsort[n]; hsort[n] = j;
		}
	}
		//Restore parents to original values
	for(i=nh-1;i>=0;i--) h[i].parent = -2-h[i].parent;
}

/** Given vector a, returns b&c that makes (a,b,c) orthonormal */
void genperp (point3d *a, point3d *b, point3d *c)
{
	float t;

	if ((a->x == 0) && (a->y == 0) && (a->z == 0))
		{ b->x = 0; b->y = 0; b->z = 0; return; }
	if ((fabs(a->x) < fabs(a->y)) && (fabs(a->x) < fabs(a->z)))
		{ t = 1.0 / sqrt(a->y*a->y + a->z*a->z); b->x = 0; b->y = a->z*t; b->z = -a->y*t; }
	else if (fabs(a->y) < fabs(a->z))
		{ t = 1.0 / sqrt(a->x*a->x + a->z*a->z); b->x = -a->z*t; b->y = 0; b->z = a->x*t; }
	else
		{ t = 1.0 / sqrt(a->x*a->x + a->y*a->y); b->x = a->y*t; b->y = -a->x*t; b->z = 0; }
	c->x = a->y*b->z - a->z*b->y;
	c->y = a->z*b->x - a->x*b->z;
	c->z = a->x*b->y - a->y*b->x;
}

	//A * B = C, find A   36*, 27ñ
	//[asx ahx agx aox][bsx bhx bgx box]   [csx chx cgx cox]
	//[asy ahy agy aoy][bsy bhy bgy boy] = [csy chy cgy coy]
	//[asz ahz agz aoz][bsz bhz bgz boz]   [csz chz cgz coz]
	//[  0   0   0   1][  0   0   0   1]   [  0   0   0   1]
void mat0 (point3d *as, point3d *ah, point3d *ag, point3d *ao,
			  point3d *bs, point3d *bh, point3d *bg, point3d *bo,
			  point3d *cs, point3d *ch, point3d *cg, point3d *co)
{
	point3d ts, th, tg, to;
	ts.x = bs->x*cs->x + bh->x*ch->x + bg->x*cg->x;
	ts.y = bs->x*cs->y + bh->x*ch->y + bg->x*cg->y;
	ts.z = bs->x*cs->z + bh->x*ch->z + bg->x*cg->z;
	th.x = bs->y*cs->x + bh->y*ch->x + bg->y*cg->x;
	th.y = bs->y*cs->y + bh->y*ch->y + bg->y*cg->y;
	th.z = bs->y*cs->z + bh->y*ch->z + bg->y*cg->z;
	tg.x = bs->z*cs->x + bh->z*ch->x + bg->z*cg->x;
	tg.y = bs->z*cs->y + bh->z*ch->y + bg->z*cg->y;
	tg.z = bs->z*cs->z + bh->z*ch->z + bg->z*cg->z;
	to.x = co->x - bo->x*ts.x - bo->y*th.x - bo->z*tg.x;
	to.y = co->y - bo->x*ts.y - bo->y*th.y - bo->z*tg.y;
	to.z = co->z - bo->x*ts.z - bo->y*th.z - bo->z*tg.z;
	(*as) = ts; (*ah) = th; (*ag) = tg; (*ao) = to;
}

	//A * B = C, find B   36*, 27ñ
	//[asx ahx agx aox][bsx bhx bgx box]   [csx chx cgx cox]
	//[asy ahy agy aoy][bsy bhy bgy boy] = [csy chy cgy coy]
	//[asz ahz agz aoz][bsz bhz bgz boz]   [csz chz cgz coz]
	//[  0   0   0   1][  0   0   0   1]   [  0   0   0   1]
void mat1 (point3d *as, point3d *ah, point3d *ag, point3d *ao,
			  point3d *bs, point3d *bh, point3d *bg, point3d *bo,
			  point3d *cs, point3d *ch, point3d *cg, point3d *co)
{
	point3d ts, th, tg, to;
	float x = co->x-ao->x, y = co->y-ao->y, z = co->z-ao->z;
	ts.x = cs->x*as->x + cs->y*as->y + cs->z*as->z;
	ts.y = cs->x*ah->x + cs->y*ah->y + cs->z*ah->z;
	ts.z = cs->x*ag->x + cs->y*ag->y + cs->z*ag->z;
	th.x = ch->x*as->x + ch->y*as->y + ch->z*as->z;
	th.y = ch->x*ah->x + ch->y*ah->y + ch->z*ah->z;
	th.z = ch->x*ag->x + ch->y*ag->y + ch->z*ag->z;
	tg.x = cg->x*as->x + cg->y*as->y + cg->z*as->z;
	tg.y = cg->x*ah->x + cg->y*ah->y + cg->z*ah->z;
	tg.z = cg->x*ag->x + cg->y*ag->y + cg->z*ag->z;
	to.x = as->x*x + as->y*y + as->z*z;
	to.y = ah->x*x + ah->y*y + ah->z*z;
	to.z = ag->x*x + ag->y*y + ag->z*z;
	(*bs) = ts; (*bh) = th; (*bg) = tg; (*bo) = to;
}

	//A * B = C, find C   36*, 27ñ
	//[asx ahx afx aox][bsx bhx bfx box]   [csx chx cfx cox]
	//[asy ahy afy aoy][bsy bhy bfy boy] = [csy chy cfy coy]
	//[asz ahz afz aoz][bsz bhz bfz boz]   [csz chz cfz coz]
	//[  0   0   0   1][  0   0   0   1]   [  0   0   0   1]

// This might be faster using point4d's instead

void mat2 (point3d *a_s, point3d *a_h, point3d *a_f, point3d *a_o,
		   point3d *b_s, point3d *b_h, point3d *b_f, point3d *b_o,
		   point3d *c_s, point3d *c_h, point3d *c_f, point3d *c_o)
{
	#if !defined(_MSVC_VER) && defined(__i386__) && !defined(NOASM)
	if (cputype&(1<<25))
	{
		_asm
		{
			mov eax, b_s
			mov edx, b_h
			movups xmm0, [eax]      ;"xmm0:   -  bs.z bs.y bs.x"
			movups xmm4, [edx]      ;"xmm4:   -  bh.z bh.y bh.x"
			mov eax, b_f
			mov edx, b_o
			movups xmm6, [eax]      ;"xmm6:   -  bf.z bf.y bf.x"
			movups xmm3, [edx]      ;"xmm3:   -  bo.z bo.y bo.x"

			mov eax, a_s
			mov edx, a_h

			movaps xmm2, xmm0       ;"xmm2:   -  bs.z bs.y bs.x"
			movaps xmm5, xmm6       ;"xmm5:   -  bf.z bf.y bf.x"
			unpcklps xmm0, xmm4     ;"xmm0: bh.y bs.y bh.x bs.x"
			unpcklps xmm6, xmm3     ;"xmm6: bo.y bf.y bo.x bf.x"
			movhlps xmm1, xmm0      ;"xmm1:   -    -  bh.y bs.y"
			movhlps xmm7, xmm6      ;"xmm7:   -    -  bo.y bf.y"
			unpckhps xmm2, xmm4     ;"xmm2:   -    -  bh.z bs.z"
			unpckhps xmm5, xmm3     ;"xmm5:   -    -  bo.z bf.z"
			movlhps xmm0, xmm6      ;"xmm0: bo.x bf.x bh.x bs.x"
			movlhps xmm1, xmm7      ;"xmm1: bo.y bf.y bh.y bs.y"
			movlhps xmm2, xmm5      ;"xmm2: bo.z bf.z bh.z bs.z"

			movss xmm3, [eax]
			shufps xmm3, xmm3, 0
			movss xmm4, [eax+4]
			shufps xmm4, xmm4, 0
			movss xmm5, [eax+8]
			shufps xmm5, xmm5, 0
			mulps xmm3, xmm0
			mulps xmm4, xmm0
			mulps xmm5, xmm0

			mov eax, a_f

			movss xmm6, [edx]
			shufps xmm6, xmm6, 0
			movss xmm7, [edx+4]
			shufps xmm7, xmm7, 0
			movss xmm0, [edx+8]
			shufps xmm0, xmm0, 0
			mulps xmm6, xmm1
			mulps xmm7, xmm1
			mulps xmm0, xmm1
			addps xmm3, xmm6
			addps xmm4, xmm7
			addps xmm5, xmm0

			mov edx, c_s

			movss xmm6, [eax]
			shufps xmm6, xmm6, 0
			movss xmm7, [eax+4]
			shufps xmm7, xmm7, 0
			movss xmm0, [eax+8]
			shufps xmm0, xmm0, 0
			mulps xmm6, xmm2
			mulps xmm7, xmm2
			mulps xmm0, xmm2
			addps xmm3, xmm6        ;xmm3: to.x tf.x th.x ts.x
			addps xmm4, xmm7        ;xmm4: to.y tf.y th.y ts.y
			addps xmm5, xmm0        ;xmm5: to.z tf.z th.z ts.z

			mov eax, c_f

			movss [edx], xmm3
			movhlps xmm0, xmm3
			movss [edx+4], xmm4
			movhlps xmm1, xmm4
			movss [edx+8], xmm5
			movhlps xmm2, xmm5
			mov edx, c_h
			movss [eax], xmm0
			movss [eax+4], xmm1
			movss [eax+8], xmm2
			shufps xmm3, xmm3, 0xb1 ;xmm3:   -  to.x   -  th.x
			shufps xmm4, xmm4, 0xb1 ;xmm4:   -  to.y   -  th.y
			shufps xmm5, xmm5, 0xb1 ;xmm5:   -  to.z   -  th.z
			mov eax, a_o
			movss [edx], xmm3
			movss [edx+4], xmm4
			movss [edx+8], xmm5
			mov edx, c_o
			movhlps xmm0, xmm3
			addss xmm0, [eax]
			movhlps xmm1, xmm4
			addss xmm1, [eax+4]
			movhlps xmm2, xmm5
			addss xmm2, [eax+8]
			movss [edx], xmm0
			movss [edx+4], xmm1
			movss [edx+8], xmm2
		}
	}
	else
	{
	#else
		point3d ts, th, tf, to;
		ts.x = a_s->x*b_s->x + a_h->x*b_s->y + a_f->x*b_s->z;
		ts.y = a_s->y*b_s->x + a_h->y*b_s->y + a_f->y*b_s->z;
		ts.z = a_s->z*b_s->x + a_h->z*b_s->y + a_f->z*b_s->z;
		th.x = a_s->x*b_h->x + a_h->x*b_h->y + a_f->x*b_h->z;
		th.y = a_s->y*b_h->x + a_h->y*b_h->y + a_f->y*b_h->z;
		th.z = a_s->z*b_h->x + a_h->z*b_h->y + a_f->z*b_h->z;
		tf.x = a_s->x*b_f->x + a_h->x*b_f->y + a_f->x*b_f->z;
		tf.y = a_s->y*b_f->x + a_h->y*b_f->y + a_f->y*b_f->z;
		tf.z = a_s->z*b_f->x + a_h->z*b_f->y + a_f->z*b_f->z;
		to.x = a_s->x*b_o->x + a_h->x*b_o->y + a_f->x*b_o->z + a_o->x;
		to.y = a_s->y*b_o->x + a_h->y*b_o->y + a_f->y*b_o->z + a_o->y;
		to.z = a_s->z*b_o->x + a_h->z*b_o->y + a_f->z*b_o->z + a_o->z;
		(*c_s) = ts; (*c_h) = th; (*c_f) = tf; (*c_o) = to;
	#endif
	#if defined(_WIN32) && defined(__i386__) && !defined(NOASM)
	}
	#endif
}

static void setlimb (kfatype *kfa, long i, long p, long trans_type, short val)
{
	point3d ps, ph, pf, pp;
	point3d qs, qh, qf, qp;
	float r[2];

		//Generate orthonormal matrix in world space for child limb
	qp = kfa->hinge[i].p[0]; qs = kfa->hinge[i].v[0]; genperp(&qs,&qh,&qf);

	switch (trans_type)
	{
		case 0: //Hinge rotate!
			//fcossin(((float)val)*(PI/32768.0),&c,&s);
			ucossin(((long)val)<<16,r);
			ph = qh; pf = qf;
			qh.x = ph.x*r[0] - pf.x*r[1]; qf.x = ph.x*r[1] + pf.x*r[0];
			qh.y = ph.y*r[0] - pf.y*r[1]; qf.y = ph.y*r[1] + pf.y*r[0];
			qh.z = ph.z*r[0] - pf.z*r[1]; qf.z = ph.z*r[1] + pf.z*r[0];
			break;
		default: _gtfo(); //tells MSVC default can't be reached
	}

		//Generate orthonormal matrix in world space for parent limb
	pp = kfa->hinge[i].p[1]; ps = kfa->hinge[i].v[1]; genperp(&ps,&ph,&pf);

		//mat0(rotrans, loc_velcro, par_velcro)
	mat0(&qs,&qh,&qf,&qp, &qs,&qh,&qf,&qp, &ps,&ph,&pf,&pp);
		//mat2(par, rotrans, parent * par_velcro * (loc_velcro x rotrans)^-1)
	mat2(&kfa->spr[p].s,&kfa->spr[p].h,&kfa->spr[p].f,&kfa->spr[p].p,
		  &qs,&qh,&qf,&qp,
		  &kfa->spr[i].s,&kfa->spr[i].h,&kfa->spr[i].f,&kfa->spr[i].p);
}

	//Uses binary search to find sequence index at time "tim"
static long kfatime2seq (kfatype *kfa, long tim)
{
	long i, a, b;

	for(a=0,b=(kfa->seqnum)-1;b-a>=2;)
		{ i = ((a+b)>>1); if (tim >= kfa->seq[i].tim) a = i; else b = i; }
	return(a);
}

void animsprite (vx5sprite *s, long ti)
{
	kfatype *kfa;
	long i, j, k, x, y, z, zz, trat;
	long trat2, z0, zz0, frm0;

	if (!(s->flags&2)) return;
	kfa = s->kfaptr; if (!kfa) return;

	z = kfatime2seq(kfa,s->kfatim);
	while (ti > 0)
	{
		z++; if (z >= kfa->seqnum) break;
		i = kfa->seq[z].tim-s->kfatim; if (i <= 0) break;
		if (i > ti) { s->kfatim += ti; break; }
		ti -= i;
		zz = ~kfa->seq[z].frm; if (zz >= 0) { if (z == zz) break; z = zz; }
		s->kfatim = kfa->seq[z].tim;
	}

// --------------------------------------------------------------------------

	z = kfatime2seq(kfa,s->kfatim); zz = z+1;
	if ((zz < kfa->seqnum) && (kfa->seq[zz].frm != ~zz))
	{
		trat = kfa->seq[zz].tim-kfa->seq[z].tim;
		if (trat) trat = shldiv16(s->kfatim-kfa->seq[z].tim,trat);
		i = kfa->seq[zz].frm; if (i < 0) zz = kfa->seq[~i].frm; else zz = i;
	} else trat = 0;
	z = kfa->seq[z].frm;
	if (z < 0)
	{
		z0 = kfatime2seq(kfa,s->okfatim); zz0 = z0+1;
		if ((zz0 < kfa->seqnum) && (kfa->seq[zz0].frm != ~zz0))
		{
			trat2 = kfa->seq[zz0].tim-kfa->seq[z0].tim;
			if (trat2) trat2 = shldiv16(s->okfatim-kfa->seq[z0].tim,trat2);
			i = kfa->seq[zz0].frm; if (i < 0) zz0 = kfa->seq[~i].frm; else zz0 = i;
		} else trat2 = 0;
		z0 = kfa->seq[z0].frm; if (z0 < 0) { z0 = zz0; trat2 = 0; }
	} else trat2 = -1;

	for(i=(kfa->numhin)-1;i>=0;i--)
	{
		if (kfa->hinge[i].parent < 0) continue;

		if (trat2 < 0) frm0 = (long)kfa->frmval[z*(kfa->numhin)+i];
		else
		{
			frm0 = (long)kfa->frmval[z0*(kfa->numhin)+i];
			if (trat2 > 0)
			{
				x = (((long)(kfa->frmval[zz0*(kfa->numhin)+i]-frm0))&65535);
				if (kfa->hinge[i].vmin == kfa->hinge[i].vmax) x = ((x<<16)>>16);
				else if ((((long)(kfa->frmval[zz0*(kfa->numhin)+i]-kfa->hinge[i].vmin))&65535) <
							(((long)(frm0-kfa->hinge[i].vmin))&65535))
					x -= 65536;
				frm0 += mulshr16(x,trat2);
			}
		}
		if (trat > 0)
		{
			x = (((long)(kfa->frmval[zz*(kfa->numhin)+i]-frm0))&65535);
			if (kfa->hinge[i].vmin == kfa->hinge[i].vmax) x = ((x<<16)>>16);
			else if ((((long)(kfa->frmval[zz*(kfa->numhin)+i]-kfa->hinge[i].vmin))&65535) <
						(((long)(frm0-kfa->hinge[i].vmin))&65535))
				x -= 65536;
			frm0 += mulshr16(x,trat);
		}
		vx5.kfaval[i] = frm0;
	}
}

static void kfadraw (vx5sprite *s)
{
	point3d tp;
	kfatype *kfa;
	long i, j, k;

	kfa = s->kfaptr; if (!kfa) return;

	for(i=(kfa->numhin)-1;i>=0;i--)
	{
		j = kfa->hingesort[i]; k = kfa->hinge[j].parent;
		if (k >= 0) setlimb(kfa,j,k,kfa->hinge[j].htype,vx5.kfaval[j]);
		else
		{
			kfa->spr[j].s = s->s;
			kfa->spr[j].h = s->h;
			kfa->spr[j].f = s->f;
			//kfa->spr[j].p = s->p;
			tp.x = kfa->hinge[j].p[0].x;
			tp.y = kfa->hinge[j].p[0].y;
			tp.z = kfa->hinge[j].p[0].z;
			kfa->spr[j].p.x = s->p.x - tp.x*s->s.x - tp.y*s->h.x - tp.z*s->f.x;
			kfa->spr[j].p.y = s->p.y - tp.x*s->s.y - tp.y*s->h.y - tp.z*s->f.y;
			kfa->spr[j].p.z = s->p.z - tp.x*s->s.z - tp.y*s->h.z - tp.z*s->f.z;
		}
		if (j < kfa->numspr) kv6draw(&kfa->spr[j]);
	}
}

//--------------------------- KFA sprite code ends ---------------------------
/** Draw a .KV6/.KFA voxel sprite to the screen.
 *  Position & orientation are specified in the vx5sprite structure. 
 *  See VOXLAP5.H for details on the structure.
 *
 *  @param spr pointer to vx5sprite
 */
void drawsprite (vx5sprite *spr)
{
	if (spr->flags&4) return;
	if (!(spr->flags&2)) kv6draw(spr); else kfadraw(spr);
}

#if 0

void setkv6 (vx5sprite *spr)
{
	point3d r0, r1;
	long x, y, vx, vy, vz;
	kv6data *kv;
	kv6voxtype *v, *ve;

	if (spr->flags&2) return;
	kv = spr->voxnum; if (!kv) return;

	vx5.minx = ?; vx5.maxx = ?+1;
	vx5.miny = ?; vx5.maxy = ?+1;
	vx5.minz = ?; vx5.maxz = ?+1;

	v = kv->vox; //.01 is to fool rounding so they aren't all even numbers
	r0.x = spr->p.x - kv->xpiv*spr->s.x - kv->ypiv*spr->h.x - kv->zpiv*spr->f.x - .01;
	r0.y = spr->p.y - kv->xpiv*spr->s.y - kv->ypiv*spr->h.y - kv->zpiv*spr->f.y - .01;
	r0.z = spr->p.z - kv->xpiv*spr->s.z - kv->ypiv*spr->h.z - kv->zpiv*spr->f.z - .01;
	vx5.colfunc = curcolfunc;
	for(x=0;x<kv->xsiz;x++)
	{
		r1 = r0;
		for(y=0;y<kv->ysiz;y++)
		{
			for(ve=&v[kv->ylen[x*kv->ysiz+y]];v<ve;v++)
			{
				ftol(spr->f.x*v->z + r1.x,&vx);
				ftol(spr->f.y*v->z + r1.y,&vy);
				ftol(spr->f.z*v->z + r1.z,&vz);
				vx5.curcol = ((v->col&0xffffff)|0x80000000);
				setcube(vx,vy,vz,-2);
			}
			r1.x += spr->h.x; r1.y += spr->h.y; r1.z += spr->h.z;
		}
		r0.x += spr->s.x; r0.y += spr->s.y; r0.z += spr->s.z;
	}
}

#else

char umost[VSID*VSID], dmost[VSID*VSID];
static kv6data *gfrezkv;
static lpoint3d gfrezx, gfrezy, gfrezz, gfrezp;
static const signed char gkv6colx[27] = {0,  0, 0, 0, 0, 1,-1, -1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,  1, 1, 1, 1,-1,-1,-1,-1};
static const signed char gkv6coly[27] = {0,  0, 0, 1,-1, 0, 0,  0, 0,-1, 1, 1, 1,-1,-1, 0, 0,-1, 1,  1, 1,-1,-1, 1, 1,-1,-1};
static const signed char gkv6colz[27] = {0,  1,-1, 0, 0, 0, 0, -1, 1, 0, 0, 1,-1, 1,-1, 1,-1, 0, 0,  1,-1, 1,-1, 1,-1, 1,-1};
long kv6colfunc (lpoint3d *p)
{
	kv6voxtype *v0, *v1, *v, *ve;
	long i, j, k, x, y, z, ox, oy, nx, ny, nz, mind, d;

	x = ((p->x*gfrezx.x + p->y*gfrezy.x + p->z*gfrezz.x + gfrezp.x)>>16);
	y = ((p->x*gfrezx.y + p->y*gfrezy.y + p->z*gfrezz.y + gfrezp.y)>>16);
	z = ((p->x*gfrezx.z + p->y*gfrezy.z + p->z*gfrezz.z + gfrezp.z)>>16);
	x = lbound0(x,gfrezkv->xsiz-1);
	y = lbound0(y,gfrezkv->ysiz-1);
	z = lbound0(z,gfrezkv->zsiz-1);

		//Process x
	v0 = gfrezkv->vox;
	if ((x<<1) < gfrezkv->xsiz) { ox = oy = j = 0; }
	else { v0 += gfrezkv->numvoxs; ox = gfrezkv->xsiz; oy = gfrezkv->ysiz; j = ox*oy; }
	v1 = v0;

	for(k=0;k<27;k++)
	{
		nx = ((long)gkv6colx[k])+x; if ((unsigned long)nx >= gfrezkv->xsiz) continue;
		ny = ((long)gkv6coly[k])+y; if ((unsigned long)ny >= gfrezkv->ysiz) continue;
		nz = ((long)gkv6colz[k])+z; if ((unsigned long)nz >= gfrezkv->zsiz) continue;

		if (nx != ox)
		{
			while (nx > ox) { v0 += gfrezkv->xlen[ox]; ox++; j += gfrezkv->ysiz; }
			while (nx < ox) { ox--; v0 -= gfrezkv->xlen[ox]; j -= gfrezkv->ysiz; }
			if ((ny<<1) < gfrezkv->ysiz) { oy = 0; v1 = v0; }
			else { oy = gfrezkv->ysiz; v1 = v0+gfrezkv->xlen[nx]; }
		}
		if (ny != oy)
		{
			while (ny > oy) { v1 += gfrezkv->ylen[j+oy]; oy++; }
			while (ny < oy) { oy--; v1 -= gfrezkv->ylen[j+oy]; }
		}

			//Process z
		for(v=v1,ve=&v1[gfrezkv->ylen[j+ny]];v<ve;v++)
			if (v->z == nz) return(v->col);
	}

		//Use brute force when all else fails.. :/
	v = gfrezkv->vox; mind = 0x7fffffff;
	for(nx=0;nx<gfrezkv->xsiz;nx++)
		for(ny=0;ny<gfrezkv->ysiz;ny++)
			for(ve=&v[gfrezkv->ylen[nx*gfrezkv->ysiz+ny]];v<ve;v++)
			{
				d = labs(x-nx)+labs(y-ny)+labs(z-v->z);
				if (d < mind) { mind = d; k = v->col; }
			}
	return(k);
}

static void kv6colfuncinit (vx5sprite *spr, float det)
{
	point3d tp, tp2;
	float f;

	gfrezkv = spr->voxnum; if (!gfrezkv) { vx5.colfunc = curcolfunc; return; }

	tp2.x = gfrezkv->xpiv + .5;
	tp2.y = gfrezkv->ypiv + .5;
	tp2.z = gfrezkv->zpiv + .5;
	tp.x = spr->p.x - spr->s.x*tp2.x - spr->h.x*tp2.y - spr->f.x*tp2.z;
	tp.y = spr->p.y - spr->s.y*tp2.x - spr->h.y*tp2.y - spr->f.y*tp2.z;
	tp.z = spr->p.z - spr->s.z*tp2.x - spr->h.z*tp2.y - spr->f.z*tp2.z;

		//spr->s.x*x + spr->h.x*y + spr->f.x*z = np.x; //Solve for x,y,z
		//spr->s.y*x + spr->h.y*y + spr->f.y*z = np.y;
		//spr->s.z*x + spr->h.z*y + spr->f.z*z = np.z;
	f = 65536.0 / det;

	tp2.x = (spr->h.y*spr->f.z - spr->h.z*spr->f.y)*f; ftol(tp2.x,&gfrezx.x);
	tp2.y = (spr->h.z*spr->f.x - spr->h.x*spr->f.z)*f; ftol(tp2.y,&gfrezy.x);
	tp2.z = (spr->h.x*spr->f.y - spr->h.y*spr->f.x)*f; ftol(tp2.z,&gfrezz.x);
	ftol(-tp.x*tp2.x - tp.y*tp2.y - tp.z*tp2.z,&gfrezp.x); gfrezp.x += 32767;

	tp2.x = (spr->f.y*spr->s.z - spr->f.z*spr->s.y)*f; ftol(tp2.x,&gfrezx.y);
	tp2.y = (spr->f.z*spr->s.x - spr->f.x*spr->s.z)*f; ftol(tp2.y,&gfrezy.y);
	tp2.z = (spr->f.x*spr->s.y - spr->f.y*spr->s.x)*f; ftol(tp2.z,&gfrezz.y);
	ftol(-tp.x*tp2.x - tp.y*tp2.y - tp.z*tp2.z,&gfrezp.y); gfrezp.y += 32767;

	tp2.x = (spr->s.y*spr->h.z - spr->s.z*spr->h.y)*f; ftol(tp2.x,&gfrezx.z);
	tp2.y = (spr->s.z*spr->h.x - spr->s.x*spr->h.z)*f; ftol(tp2.y,&gfrezy.z);
	tp2.z = (spr->s.x*spr->h.y - spr->s.y*spr->h.x)*f; ftol(tp2.z,&gfrezz.z);
	ftol(-tp.x*tp2.x - tp.y*tp2.y - tp.z*tp2.z,&gfrezp.z); gfrezp.z += 32768;
}

#define LSC3 8 //2 for testing, 8 is normal
typedef struct
{
	long xo, yo, zo, xu, yu, zu, xv, yv, zv, d, msk, pzi;
	long xmino, ymino, xmaxo, ymaxo, xusc, yusc, xvsc, yvsc;
} gfrezt;
static gfrezt gfrez[6];
typedef struct { char z[2]; long n; } slstype;
void setkv6 (vx5sprite *spr, long dacol)
{
	point3d tp, tp2; float f, det;
	long i, j, k, x, y, z, c, d, x0, y0, z0, x1, y1, z1, xi, yi, zi;
	long xo, yo, zo, xu, yu, zu, xv, yv, zv, stu, stv, tu, tv;
	long xx, yy, xmin, xmax, ymin, ymax, isrhs, ihxi, ihyi, ihzi, syshpit;
	long isx, isy, isz, ihx, ihy, ihz, ifx, ify, ifz, iox, ioy, ioz;
	long sx, sy, sx0, sy0, sz0, rx, ry, rz, pz, dcnt, dcnt2, vismask, xysiz;
	long bx0, by0, bz0, bx1, by1, bz1, *lptr, *shead, shpit, scnt, sstop;
	gfrezt *gf;
	slstype *slst;
	kv6data *kv;
	kv6voxtype *v0, *v1, *v2, *v3;
	void (*modslab)(long *, long, long);

	if (spr->flags&2) return;
	kv = spr->voxnum; if (!kv) return;

		//Calculate top-left-up corner in VXL world coordinates
	tp.x = kv->xpiv + .5;
	tp.y = kv->ypiv + .5;
	tp.z = kv->zpiv + .5;
	tp2.x = spr->p.x - spr->s.x*tp.x - spr->h.x*tp.y - spr->f.x*tp.z;
	tp2.y = spr->p.y - spr->s.y*tp.x - spr->h.y*tp.y - spr->f.y*tp.z;
	tp2.z = spr->p.z - spr->s.z*tp.x - spr->h.z*tp.y - spr->f.z*tp.z;

		//Get bounding x-y box of entire freeze area:
	bx0 = VSID; by0 = VSID; bz0 = MAXZDIM; bx1 = 0; by1 = 0; bz1 = 0;
	for(z=kv->zsiz;z>=0;z-=kv->zsiz)
		for(y=kv->ysiz;y>=0;y-=kv->ysiz)
			for(x=kv->xsiz;x>=0;x-=kv->xsiz)
			{
				ftol(spr->s.x*(float)x + spr->h.x*(float)y + spr->f.x*(float)z + tp2.x,&i);
				if (i < bx0) bx0 = i;
				if (i > bx1) bx1 = i;
				ftol(spr->s.y*(float)x + spr->h.y*(float)y + spr->f.y*(float)z + tp2.y,&i);
				if (i < by0) by0 = i;
				if (i > by1) by1 = i;
				ftol(spr->s.z*(float)x + spr->h.z*(float)y + spr->f.z*(float)z + tp2.z,&i);
				if (i < bz0) bz0 = i;
				if (i > bz1) bz1 = i;
			}
	bx0 -= 2; if (bx0 < 0) bx0 = 0;
	by0 -= 2; if (by0 < 0) by0 = 0;
	bz0 -= 2; if (bz0 < 0) bz0 = 0;
	bx1 += 2; if (bx1 > VSID) bx1 = VSID;
	by1 += 2; if (by1 > VSID) by1 = VSID;
	bz1 += 2; if (bz1 > MAXZDIM) bz1 = MAXZDIM;
	vx5.minx = bx0; vx5.maxx = bx1;
	vx5.miny = by0; vx5.maxy = by1;
	vx5.minz = bz0; vx5.maxz = bz1;

	shpit = bx1-bx0; i = (by1-by0)*shpit*sizeof(shead[0]);
		//Make sure to use array that's big enough: umost is 1MB
	shead = (long *)(((long)umost) - (by0*shpit+bx0)*sizeof(shead[0]));
	slst = (slstype *)(((long)umost)+i);
	scnt = 1; sstop = (sizeof(umost)-i)/sizeof(slstype);
	memset(umost,0,i);

	f = (float)(1<<LSC3);
	ftol(spr->s.x*f,&isx); ftol(spr->s.y*f,&isy); ftol(spr->s.z*f,&isz);
	ftol(spr->h.x*f,&ihx); ftol(spr->h.y*f,&ihy); ftol(spr->h.z*f,&ihz);
	ftol(spr->f.x*f,&ifx); ftol(spr->f.y*f,&ify); ftol(spr->f.z*f,&ifz);
	ftol(tp2.x*f,&iox);
	ftol(tp2.y*f,&ioy);
	ftol(tp2.z*f,&ioz);

		//Determine whether sprite is RHS(1) or LHS(0)
	det = (spr->h.y*spr->f.z - spr->h.z*spr->f.y)*spr->s.x +
			(spr->h.z*spr->f.x - spr->h.x*spr->f.z)*spr->s.y +
			(spr->h.x*spr->f.y - spr->h.y*spr->f.x)*spr->s.z;
	if ((*(long *)&det) > 0) isrhs = 1;
	else if ((*(long *)&det) < 0) isrhs = 0;
	else return;

	xi = (((ifx*ihy-ihx*ify)>>31)|1);
	yi = (((isx*ify-ifx*isy)>>31)|1);
	zi = (((ihx*isy-isx*ihy)>>31)|1);
	if (xi > 0) { x0 = 0; x1 = kv->xsiz; } else { x0 = kv->xsiz-1; x1 = -1; }
	if (yi > 0) { y0 = 0; y1 = kv->ysiz; } else { y0 = kv->ysiz-1; y1 = -1; }
	if (zi > 0) { z0 = 0; z1 = kv->zsiz; } else { z0 = kv->zsiz-1; z1 = -1; }

	vismask = (zi<<3)+24 + (yi<<1)+6 + (xi>>1)+2;

	dcnt = 0;
	for(j=2;j;j--)
	{
		dcnt2 = dcnt;
		vismask = ~vismask;
		for(i=1;i<64;i+=i)
		{
			if (!(vismask&i)) continue;

			if (i&0x15) { xo = yo = zo = 0; }
			else if (i == 2) { xo = isx; yo = isy; zo = isz; }
			else if (i == 8) { xo = ihx; yo = ihy; zo = ihz; }
			else             { xo = ifx; yo = ify; zo = ifz; }

				  if (i&3)  { xu = ihx; yu = ihy; zu = ihz; xv = ifx; yv = ify; zv = ifz; }
			else if (i&12) { xu = isx; yu = isy; zu = isz; xv = ifx; yv = ify; zv = ifz; }
			else           { xu = isx; yu = isy; zu = isz; xv = ihx; yv = ihy; zv = ihz; }

			if ((yu < 0) || ((!yu) && (xu < 0)))
				{ xo += xu; yo += yu; zo += zu; xu = -xu; yu = -yu; zu = -zu; }
			if ((yv < 0) || ((!yv) && (xv < 0)))
				{ xo += xv; yo += yv; zo += zv; xv = -xv; yv = -yv; zv = -zv; }
			d = xv*yu - xu*yv; if (!d) continue;
			if (d < 0)
			{
				k = xu; xu = xv; xv = k;
				k = yu; yu = yv; yv = k;
				k = zu; zu = zv; zv = k; d = -d;
			}
			xmin = ymin = xmax = ymax = 0;
			if (xu < 0) xmin += xu; else xmax += xu;
			if (yu < 0) ymin += yu; else ymax += yu;
			if (xv < 0) xmin += xv; else xmax += xv;
			if (yv < 0) ymin += yv; else ymax += yv;

			gf = &gfrez[dcnt];
			gf->xo = xo; gf->yo = yo; gf->zo = zo;
			gf->xu = xu; gf->yu = yu;
			gf->xv = xv; gf->yv = yv;
			gf->xmino = xmin; gf->ymino = ymin;
			gf->xmaxo = xmax; gf->ymaxo = ymax;
			gf->xusc = (xu<<LSC3); gf->yusc = (yu<<LSC3);
			gf->xvsc = (xv<<LSC3); gf->yvsc = (yv<<LSC3);
			gf->d = d; gf->msk = i;

			f = 1.0 / (float)d;
			ftol(((float)gf->yusc * (float)zv - (float)gf->yvsc * (float)zu) * f,&gf->pzi);
			f *= 4194304.0;
			ftol((float)zu*f,&gf->zu);
			ftol((float)zv*f,&gf->zv);

			dcnt++;
		}
	}

	ihxi = ihx*yi;
	ihyi = ihy*yi;
	ihzi = ihz*yi;

	if (xi < 0) v0 = kv->vox+kv->numvoxs; else v0 = kv->vox;
	for(x=x0;x!=x1;x+=xi)
	{
		i = (long)kv->xlen[x];
		if (xi < 0) v0 -= i;
		if (yi < 0) v1 = v0+i; else v1 = v0;
		if (xi >= 0) v0 += i;
		xysiz = x*kv->ysiz;
		sx0 = isx*x + ihx*y0 + iox;
		sy0 = isy*x + ihy*y0 + ioy;
		sz0 = isz*x + ihz*y0 + ioz;
		for(y=y0;y!=y1;y+=yi)
		{
			i = (long)kv->ylen[xysiz+y];
			if (yi < 0) v1 -= i;
			if (zi < 0) { v2 = v1+i-1; v3 = v1-1; }
					 else { v2 = v1; v3 = v1+i; }
			if (yi >= 0) v1 += i;
			while (v2 != v3)
			{
				z = v2->z; //c = v2->col;
				rx = ifx*z + sx0;
				ry = ify*z + sy0;
				rz = ifz*z + sz0;
				for(i=0;i<dcnt;i++)
				{
					gf = &gfrez[i]; if (!(v2->vis&gf->msk)) continue;
					xo = gf->xo + rx;
					yo = gf->yo + ry;
					zo = gf->zo + rz;
					xmin = ((gf->xmino + xo)>>LSC3); if (xmin < 0) xmin = 0;
					ymin = ((gf->ymino + yo)>>LSC3); if (ymin < 0) ymin = 0;
					xmax = ((gf->xmaxo + xo)>>LSC3)+1; if (xmax > VSID) xmax = VSID;
					ymax = ((gf->ymaxo + yo)>>LSC3)+1; if (ymax > VSID) ymax = VSID;
					xx = (xmin<<LSC3) - xo;
					yy = (ymin<<LSC3) - yo;
					stu = yy*gf->xu - xx*gf->yu;
					stv = yy*gf->xv - xx*gf->yv - gf->d;
					syshpit = ymin*shpit;
					for(sy=ymin;sy<ymax;sy++,syshpit+=shpit)
					{
						tu = stu; stu += gf->xusc;
						tv = stv; stv += gf->xvsc;
						sx = xmin;
						while ((tu&tv) >= 0)
						{
							sx++; if (sx >= xmax) goto freezesprcont;
							tu -= gf->yusc; tv -= gf->yvsc;
						}
						tu = ~tu; tv += gf->d;
						pz = dmulshr22(tu,gf->zv,tv,gf->zu) + zo; j = syshpit+sx;
						tu -= gf->d; tv = ~tv;
						while ((tu&tv) < 0)
						{
							if (i < dcnt2)
							{
								if (scnt >= sstop) return; //OUT OF BUFFER SPACE!
								slst[scnt].z[0] = (char)lbound0(pz>>LSC3,255);
								slst[scnt].n = shead[j]; shead[j] = scnt; scnt++;
							}
							else slst[shead[j]].z[1] = (char)lbound0(pz>>LSC3,255);
							tu += gf->yusc; tv += gf->yvsc; pz += gf->pzi; j++;
						}
freezesprcont:;
					}
				}
				v2 += zi;
			}
			sx0 += ihxi; sy0 += ihyi; sz0 += ihzi;
		}
	}

	if (dacol == -1) modslab = delslab; else modslab = insslab;

	if (vx5.colfunc == kv6colfunc) kv6colfuncinit(spr,det);

	j = by0*shpit+bx0;
	for(sy=by0;sy<by1;sy++)
		for(sx=bx0;sx<bx1;sx++,j++)
		{
			i = shead[j]; if (!i) continue;
			lptr = scum2(sx,sy);
			do
			{
				modslab(lptr,(long)slst[i].z[isrhs],(long)slst[i].z[isrhs^1]);
				i = slst[i].n;
			} while (i);
		}
	scum2finish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
}

#endif

	//Sprite structure is already allocated
	//kv6, vox, xlen, ylen are all malloced in here!

static void setlighting (long x0, long y0, long z0, long x1, long y1, long z1, long lval)
{
	long i, x, y;
	char *v;

	x0 = max(x0,0); x1 = min(x1,VSID);
	y0 = max(y0,0); y1 = min(y1,VSID);
	z0 = max(z0,0); z1 = min(z1,MAXZDIM);

	lval <<= 24;

		//Set 4th byte of colors to full intensity
	for(y=y0;y<y1;y++)
		for(x=x0;x<x1;x++)
		{
			for(v=sptr[y*VSID+x];v[0];v+=v[0]*4)
				for(i=1;i<v[0];i++)
					(*(long *)&v[i<<2]) = (((*(long *)&v[i<<2])&0xffffff)|lval);
			for(i=1;i<=v[2]-v[1]+1;i++)
				(*(long *)&v[i<<2]) = (((*(long *)&v[i<<2])&0xffffff)|lval);
		}
}

	//Updates Lighting, Mip-mapping, and Floating objects list


#define BBOXSIZ 256
static bboxtyp bbox[BBOXSIZ];
static long bboxnum = 0;
void updatevxl ()
{
	long i;

	for(i=bboxnum-1;i>=0;i--)
	{
		if (vx5.lightmode)
			updatelighting( bbox[i].x0,bbox[i].y0,bbox[i].z0,bbox[i].x1,bbox[i].y1,bbox[i].z1 ) ;
		if (vx5.vxlmipuse > 1)
			genmipvxl(bbox[i].x0,bbox[i].y0,bbox[i].x1,bbox[i].y1);
		if ((vx5.fallcheck) && (bbox[i].csgdel))
			checkfloatinbox(bbox[i].x0,bbox[i].y0,bbox[i].z0,bbox[i].x1,bbox[i].y1,bbox[i].z1);
	}
	bboxnum = 0;
}

void updatebbox (long x0, long y0, long z0, long x1, long y1, long z1, long csgdel)
{
	long i;

	if ((x0 >= x1) || (y0 >= y1) || (z0 >= z1)) return;
	for(i=bboxnum-1;i>=0;i--)
	{
		if ((x0 >= bbox[i].x1) || (bbox[i].x0 >= x1)) continue;
		if ((y0 >= bbox[i].y1) || (bbox[i].y0 >= y1)) continue;
		if ((z0 >= bbox[i].z1) || (bbox[i].z0 >= z1)) continue;
		if (bbox[i].x0 < x0) x0 = bbox[i].x0;
		if (bbox[i].y0 < y0) y0 = bbox[i].y0;
		if (bbox[i].z0 < z0) z0 = bbox[i].z0;
		if (bbox[i].x1 > x1) x1 = bbox[i].x1;
		if (bbox[i].y1 > y1) y1 = bbox[i].y1;
		if (bbox[i].z1 > z1) z1 = bbox[i].z1;
		csgdel |= bbox[i].csgdel;
		bboxnum--; bbox[i] = bbox[bboxnum];
	}
	bbox[bboxnum].x0 = x0; bbox[bboxnum].x1 = x1;
	bbox[bboxnum].y0 = y0; bbox[bboxnum].y1 = y1;
	bbox[bboxnum].z0 = z0; bbox[bboxnum].z1 = z1;
	bbox[bboxnum].csgdel = csgdel; bboxnum++;
	if (bboxnum >= BBOXSIZ) updatevxl();
}

static long lightlst[MAXLIGHTS];
static float lightsub[MAXLIGHTS];
/** Re-calculates lighting byte #4 of all voxels inside bounding box
 *  Takes 6 longs ( which should be 2 vectors ) representing min and max extents
 *  of an AABB.
 *  @param x0
 *  @param y0
 *  @param z0
 *  @param x1
 *  @param y1
 *  @param z1
*/
void updatelighting (long x0, long y0, long z0, long x1, long y1, long z1)
{
	point3d tp;
	float f, g, h, fx, fy, fz;
	long i, j, x, y, z, sz0, sz1, offs, cstat, lightcnt;
	long x2, y2, x3, y3;
	char *v;

	if (!vx5.lightmode) return;
	xbsox = -17;

	x0 = MAX(x0-ESTNORMRAD,0); x1 = MIN(x1+ESTNORMRAD,VSID);
	y0 = MAX(y0-ESTNORMRAD,0); y1 = MIN(y1+ESTNORMRAD,VSID);
	z0 = MAX(z0-ESTNORMRAD,0); z1 = MIN(z1+ESTNORMRAD,MAXZDIM);

	x2 = x0; y2 = y0;
	x3 = x1; y3 = y1;
	for(y0=y2;y0<y3;y0=y1)
	{
		y1 = MIN(y0+64,y3);  //"justfly -" (256 lights): +1024:41sec 512:29 256:24 128:22 64:21 32:21 16:21
		for(x0=x2;x0<x3;x0=x1)
		{
			x1 = MIN(x0+64,x3);


			if (vx5.lightmode == 2)
			{
				lightcnt = 0; //Find which lights are close enough to affect rectangle
				for(i=vx5.numlights-1;i>=0;i--)
				{
					ftol(vx5.lightsrc[i].p.x,&x);
					ftol(vx5.lightsrc[i].p.y,&y);
					ftol(vx5.lightsrc[i].p.z,&z);
					if (x < x0) x -= x0; else if (x > x1) x -= x1; else x = 0;
					if (y < y0) y -= y0; else if (y > y1) y -= y1; else y = 0;
					if (z < z0) z -= z0; else if (z > z1) z -= z1; else z = 0;
					f = vx5.lightsrc[i].r2;
					if ((float)(x*x+y*y+z*z) < f)
					{
						lightlst[lightcnt] = i;
						lightsub[lightcnt] = (1.0f*(f_rsqrt(f)))/f;
						lightcnt++;
					}
				}
			}

			for(y=y0;y<y1;y++)
				for(x=x0;x<x1;x++)
				{
					v = sptr[y*VSID+x]; cstat = 0;
					while (1)
					{
						if (!cstat)
						{
							sz0 = ((long)v[1]); sz1 = ((long)v[2])+1; offs = 7-(sz0<<2);
							cstat = 1;
						}
						else
						{
							sz0 = ((long)v[2])-((long)v[1])-((long)v[0])+2;
							if (!v[0]) break; v += v[0]*4;
							sz1 = ((long)v[3]); sz0 += sz1; offs = 3-(sz1<<2);
							cstat = 0;
						}
						if (z0 > sz0) sz0 = z0;
						if (z1 < sz1) sz1 = z1;
						if (vx5.lightmode < 2)
						{
							for(z=sz0;z<sz1;z++)
							{
								estnorm(x,y,z,&tp);
								ftol((tp.y*.5+tp.z)*64.f+103.5f,&i);
								v[(z<<2)+offs] = *(char *)&i;
							}
						}
						else
						{
							for(z=sz0;z<sz1;z++)
							{
								estnorm(x,y,z,&tp);
								f = (tp.y*.5+tp.z)*16+47.5;
								for(i=lightcnt-1;i>=0;i--)
								{
									j = lightlst[i];
									fx = vx5.lightsrc[j].p.x-(float)x;
									fy = vx5.lightsrc[j].p.y-(float)y;
									fz = vx5.lightsrc[j].p.z-(float)z;
									h = tp.x*fx+tp.y*fy+tp.z*fz; if (*(long *)&h >= 0) continue;
									g = fx*fx+fy*fy+fz*fz; if (g >= vx5.lightsrc[j].r2) continue;


									#if !defined(NOASM) // equiv g = 1.0/(g*sqrt(g))-lightsub[i]; //1.0/g;
									if (cputype&(1<<25))
									{
										#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
										__asm__ __volatile__
										(
											"rcpss	%[x0], %[x1]\n"     //xmm1=1/g
											"rsqrtss	%[x0], %[x0]\n" //xmm0=1/sqrt(g)
											"mulss	%[x0], %[x1]\n"     //xmm1=1/(g*sqrt(g))
											"subss	%c[lsub](,%[i],4), %[x1]\n"
											: [x1] "=x" (g)
											: [x0] "x" (g), [i] "r" (i), [lsub] "p" (lightsub)
											:
										);
										#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
										_asm
										{
											movss	xmm0, g        //xmm0=g
											rcpss	xmm1, xmm0     //xmm1=1/g
											rsqrtss	xmm0, xmm0     //xmm0=1/sqrt(g)
											mulss	xmm1, xmm0     //xmm1=1/(g*sqrt(g))
											mov	eax, i
											subss	xmm1, lightsub[eax*4]
											movss	g, xmm1
										}
										#endif
									}
									else
									{
										#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
										__asm__ __volatile__
										(
											"pfrcp	%[y0], %%mm1\n"
											"pfrsqrt	%[y0], %[y0]\n"
											"pfmul	%%mm1, %[y0]\n"
											"pfsub	%c[lsub](,%[i],4), %[y0]\n"
											"femms\n"
											: [y0] "+y" (g)
											: [i]   "r" (i), [lsub] "p" (lightsub)
											: "mm1"
										);
										#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
										_asm
										{
											movd mm0, g
											pfrcp	mm1, mm0
											pfrsqrt	mm0, mm0
											pfmul	mm0, mm1
											mov	eax, i
											pfsub	mm0, lightsub[eax*4]
											movd	g, xmm0
											femms
										}
										#endif
									}
									#else //C Default
									g = 1.0/(g*sqrt(g))-lightsub[i]; //1.0/g;
									#endif
									f -= g*h*vx5.lightsrc[j].sc;
								}
								if (*(long *)&f > 0x437f0000) f = 255; //0x437f0000 is 255.0
								ftol(f,&i);
								v[(z<<2)+offs] = *(char *)&i;
							}
						}
					}
				}
		}
	}
}

//----------------------------------------------------------------------------
/** @warning This function can evilquit on bad malloc
 * Since voxlap is currently a software renderer and I don't have any system
 * dependent code in it, you must provide it with the frame buffer. You
 * MUST call this once per frame, AFTER startdirectdraw(), but BEFORE any
 * functions that access the frame buffer.
 * @param p pointer to the top-left corner of the frame
 * @param b pitch (bytes per line)
 * @param x frame width
 * @param y frame height
 */
void voxsetframebuffer (long p, long b, long x, long y)
{
	long i;

	frameplace = p;
	if (x > MAXXDIM) x = MAXXDIM; //This sucks, but it crashes without it
	if (y > MAXYDIM) y = MAXYDIM;

		//Set global variables used by kv6draw's PIII asm (drawboundcube)
	qsum1[3] = qsum1[1] = 0x7fff-y; qsum1[2] = qsum1[0] = 0x7fff-x;
	kv6bytesperline = qbplbpp[1] = b; qbplbpp[0] = 4;
	kv6frameplace = p - (qsum1[0]*qbplbpp[0] + qsum1[1]*qbplbpp[1]);

	if ((b != ylookup[1]) || (x != xres) || (y != yres))
	{
		bytesperline = b; xres = x; yres = y; xres4 = (xres<<2);
		ylookup[0] = 0; for(i=0;i<yres;i++) ylookup[i+1] = ylookup[i]+bytesperline;
		//gihx = gihz = (float)xres*.5f; gihy = (float)yres*.5f; //BAD!!!
#if (USEZBUFFER == 1)
		if ((ylookup[yres]+256 > zbuffersiz) || (!zbuffermem))  //Increase Z buffer size if too small
		{
			if (zbuffermem) { free(zbuffermem); zbuffermem = 0; }
			zbuffersiz = ylookup[yres]+256;
			if (!(zbuffermem = (long *)malloc(zbuffersiz))) evilquit("voxsetframebuffer: allocation too big");
		}
#endif
	}
#if (USEZBUFFER == 1)
		//zbuffer aligns its memory to the same pixel boundaries as the screen!
		//WARNING: Pentium 4's L2 cache has severe slowdowns when 65536-64 <= (zbufoff&65535) < 64
	zbufoff = (((((long)zbuffermem)-frameplace-128)+255)&~255)+128;
#endif
	uurend = &uurendmem[((frameplace&4)^(((long)uurendmem)&4))>>2];

	if (vx5.fogcol >= 0)
	{
		fogcol = (((int64_t)(vx5.fogcol&0xff0000))<<16) +
					(((int64_t)(vx5.fogcol&0x00ff00))<< 8) +
					(((int64_t)(vx5.fogcol&0x0000ff))    );

		if (vx5.maxscandist > 2047) vx5.maxscandist = 1024;
		if ((vx5.maxscandist != ofogdist) && (vx5.maxscandist > 0))
		{
			ofogdist = vx5.maxscandist;

			//foglut[?>>20] = min(?*32767/vx5.maxscandist,32767)
#if 0
			long j, k, l;
			j = 0; l = 0x7fffffff/vx5.maxscandist;
			for(i=0;i<2048;i++)
			{
				k = (j>>16); j += l;
				if (k < 0) break;
				foglut[i] = (((__int64)k)<<32)+(((__int64)k)<<16)+((__int64)k);
			}
			while (i < 2048) foglut[i++] = all32767;
#else
			i = 0x7fffffff/vx5.maxscandist;
			#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM) //gcc inline asm

			__asm__ __volatile__
			(
				".intel_syntax noprefix\n"
				"xor	eax, eax\n"
				"mov	ecx, -2048*8\n"
			".Lfogbeg:\n"
				"movd	mm0, eax\n"
				"add	eax, edx\n"
				"jo	short .Lfogend\n"
				"pshufw	mm0, mm0, 0x55\n"
				"movq	foglut[ecx+2048*8], mm0\n"
				"add	ecx, 8\n"
				"js	short .Lfogbeg\n"
				"jmp	short .Lfogend2\n"
			".Lfogend:\n"
				"movq	mm0, all32767\n"
			".Lfogbeg2:\n"
				"movq	foglut[ecx+2048*8], mm0\n"
				"add	ecx, 8\n"
				"js	short .Lfogbeg2\n"
			".Lfogend2:\n"
				"emms\n"
				".att_syntax prefix\n"
				:
				: [i] "d" (i)
				: "eax", "ecx", "mm0"
			);
			#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
			_asm
			{
				xor	eax, eax
				mov	ecx, -2048*8
				mov	edx, i
			fogbeg:
				movd	mm0, eax
				add	eax, edx
				jo	short fogend
				pshufw	mm0, mm0, 0x55
				movq	foglut[ecx+2048*8], mm0
				add	ecx, 8
				js	short fogbeg
				jmp	short fogend2
			fogend:
				movq	mm0, all32767
			fogbeg2:
				movq	foglut[ecx+2048*8], mm0
				add	ecx, 8
				js	short fogbeg2
			fogend2:
				emms
			}
			#else // C Default
			#error No C Default yet defined
			#endif
#endif // 0
		}
	} else ofogdist = -1;

	if (cputype&(1<<25)) drawboundcubesseinit(); else drawboundcube3dninit();
}


/** same as screencapture32bit, but with zbuffer */
long capture_zbuffer_32bit (const char *fname)
{
#if USEZBUFFER ==1
	long p, x, y;

	pngoutopenfile(fname,xres,yres);
	p = *(long *)&zbuffermem;
	for(y=0;y<yres;y++,p+=bytesperline)
		for(x=0;x<xres;x++)
			pngoutputpixel(*(long *)(p+(x<<2)));
#endif
	return(0);
}


#if 0
/** @note This doesn't speed it up and it only makes it crash on some computers :/ */
static inline void fixsse ()
{
	static long asm32;
	#if (defined(__GNUC__) && defined(__i386__) && (USEV5ASM != 0))
	__asm__ __volatile__
	(
		".intel_syntax noprefix\n"
		"stmxcsr	[asm32]\n" //Default is:0x1f80
		"or	asm32, 0x8040\n"   //enable ftz&daz to prevent slow denormals!
		"ldmxcsr	[asm32]\n"
		".att_syntax prefix\n"
	);
	#elif (defined(_MSC_VER) && defined(__i386__) && (USEV5ASM != 0))
	_asm
	{
		stmxcsr	[asm32]   //Default is:0x1f80
		or	asm32, 0x8040 //enable ftz&daz to prevent slow denormals!
		ldmxcsr	[asm32]
	}
	#endif
}
#endif

void freekv6 (kv6data *kv6)
{
	if (kv6->lowermip) freekv6(kv6->lowermip); //NOTE: dangerous - recursive!
	free((void *)kv6);
}

void uninitvoxlap ()
{
	if (sxlbuf) { free(sxlbuf); sxlbuf = 0; }

	if (vbuf) { free(vbuf); vbuf = 0; }
	if (vbit) { free(vbit); vbit = 0; }

	if (khashbuf)
	{     //Free all KV6&KFA on hash list
		long i, j;
		kfatype *kfp;
		for(i=0;i<khashpos;i+=strlen(&khashbuf[i+9])+10)
		{
			switch (khashbuf[i+8])
			{
				case 0: //KV6
					freekv6(*(kv6data **)&khashbuf[i+4]);
					break;
				case 1: //KFA
					kfp = *(kfatype **)&khashbuf[i+4];
					if (!kfp) continue;
					if (kfp->seq) free((void *)kfp->seq);
					if (kfp->frmval) free((void *)kfp->frmval);
					if (kfp->hingesort) free((void *)kfp->hingesort);
					if (kfp->hinge) free((void *)kfp->hinge);
					if (kfp->spr)
					{
						for(j=kfp->numspr-1;j>=0;j--)
							if (kfp->spr[j].voxnum)
								freekv6((kv6data *)kfp->spr[j].voxnum);
						free((void *)kfp->spr);
					}
					free((void *)kfp);
					break;
				default: _gtfo(); //tells MSVC default can't be reached
			}
		}
		free(khashbuf); khashbuf = 0; khashpos = khashsiz = 0;
	}

	if (skylng) { free((void *)skylng); skylng = 0; }
	if (skylat) { free((void *)skylat); skylat = 0; }
	if (skypic) { free((void *)skypic); skypic = skyoff = 0; }

	if (vx5.pic) { free(vx5.pic); vx5.pic = 0; }
#if (USEZBUFFER == 1)
	if (zbuffermem) { free(zbuffermem); zbuffermem = 0; }
#endif
	if (radarmem) { free(radarmem); radarmem = 0; radar = 0; }
}

long initvoxlap ()
{
	int64_t q;
	long i, j, k, z, zz;
	float f, ff;

	v5_asm_dep_unlock();

		//CPU Must have: FPU,RDTSC,CMOV,MMX,MMX+
	if ((cputype&((1<<0)|(1<<4)|(1<<15)|(1<<22)|(1<<23))) !=
					 ((1<<0)|(1<<4)|(1<<15)|(1<<22)|(1<<23))) return(-1);
		//CPU UNSUPPORTED!
	if ((!(cputype&(1<<25))) && //SSE
		(!((cputype&((1<<30)|(1<<31))) == ((1<<30)|(1<<31))))) //3DNow!+
		return(-1);
	//if (cputype&(1<<25)) fixsse(); //SSE

	  //WARNING: xres&yres are local to VOXLAP5.C so don't rely on them here!
	if (!(radarmem = (long *)malloc(max((((MAXXDIM*MAXYDIM*27)>>1)+7)&~7,(VSID+4)*3*SCPITCH*4+8))))
		return(-1);
	radar = (long *)((((long)radarmem)+7)&~7);

	for(i=0;i<32;i++) { xbsflor[i] = (-1<<i); xbsceil[i] = ~xbsflor[i]; }

		//Setsphere precalculations (factr[] tables) (Derivation in POWCALC.BAS)
		//   if (!factr[z][0]) z's prime else factr[z][0]*factr[z][1] == z
	factr[2][0] = 0; i = 1; j = 9; k = 0;
	for(z=3;z<SETSPHMAXRAD;z+=2)
	{
		if (z == j) { j += (i<<2)+12; i += 2; }
		factr[z][0] = 0; factr[k][1] = z;
		for(zz=3;zz<=i;zz=factr[zz][1])
			if (!(z%zz)) { factr[z][0] = zz; factr[z][1] = z/zz; break; }
		if (!factr[z][0]) k = z;
		factr[z+1][0] = ((z+1)>>1); factr[z+1][1] = 2;
	}
	for(z=1;z<SETSPHMAXRAD;z++) logint[z] = log((double)z);

#if (ESTNORMRAD == 2)
		//LUT for ESTNORM
	fsqrecip[0] = 0.f; fsqrecip[1] = 1.f;
	fsqrecip[2] = (float)(1.f/sqrt(2.f)); fsqrecip[3] = (float)1.f/sqrt(3.f);
	for(z=4,i=3;z<sizeof(fsqrecip)/sizeof(fsqrecip[0]);z+=6) //fsqrecip[z] = 1/sqrt(z);
	{
		fsqrecip[z+0] = fsqrecip[(z+0)>>1]*fsqrecip[2];
		fsqrecip[z+2] = fsqrecip[(z+2)>>1]*fsqrecip[2];
		fsqrecip[z+4] = fsqrecip[(z+4)>>1]*fsqrecip[2];
		fsqrecip[z+5] = fsqrecip[i]*fsqrecip[3]; i += 2;

		f = (fsqrecip[z+0]+fsqrecip[z+2])*.5f;
		if (z <= 22) f = (1.5f-(.5f*((float)(z+1))) * f*f)*f;
		fsqrecip[z+1] = (1.5f-(.5f*((float)(z+1))) * f*f)*f;

		f = (fsqrecip[z+2]+fsqrecip[z+4])*.5f;
		if (z <= 22) f = (1.5f-(.5f*((float)(z+3))) * f*f)*f;
		fsqrecip[z+3] = (1.5f-(.5f*((float)(z+3))) * f*f)*f;
	}
#endif

		//Lookup table to save 1 divide for gline()
	for(i=1;i<CMPRECIPSIZ;i++) cmprecip[i] = CMPPREC/(float)i;

		//Flashscan equal-angle compare table
	for(i=0;i<(1<<LOGFLASHVANG)*8;i++)
	{
		if (!(i&((1<<LOGFLASHVANG)-1)))
			j = (gfclookup[i>>LOGFLASHVANG]<<4)+8 - (1<<LOGFLASHVANG)*64;
		gfc[i].y = j; j += 64*2;
		ftol(sqrt((1<<(LOGFLASHVANG<<1))*64.f*64.f-gfc[i].y*gfc[i].y),&gfc[i].x);
	}

		//Init norm flash variables:
	ff = (float)GSIZ*.5f; // /(1);
	for(z=1;z<(GSIZ>>1);z++)
	{
		ffxptr = &ffx[(z+1)*z-1];
		f = ff; ff = (float)GSIZ*.5f/((float)z+1);
		for(zz=-z;zz<=z;zz++)
		{
			if (zz <= 0) i = (long)(((float)zz-.5f)*f); else i = (long)(((float)zz-.5f)*ff);
			if (zz >= 0) j = (long)(((float)zz+.5f)*f); else j = (long)(((float)zz+.5f)*ff);
			ffxptr[zz].x = (unsigned short)max(i+(GSIZ>>1),0);
			ffxptr[zz].y = (unsigned short)min(j+(GSIZ>>1),GSIZ);
		}
	}
	for(i=0;i<=25*5;i+=5) xbsbuf[i] = 0x00000000ffffffff;
	for(z=0;z<32;z++) { p2c[z] = (1<<z); p2m[z] = p2c[z]-1; }

		//Drawtile lookup table:
	//q = 0;
	//for(i=0;i<256;i++) { alphalookup[i] = q; q += 0x1000100010; }

		//Initialize univec normals (for KV6 lighting)
	equivecinit(255);
	//for(i=0;i<255;i++)
	//{
	//   univec[i].z = ((float)((i<<1)-254))/255.0;
	//   f = sqrt(1.0 - univec[i].z*univec[i].z);
	//   fcossin((float)i*(GOLDRAT*PI*2),&univec[i].x,&univec[i].y);
	//   univec[i].x *= f; univec[i].y *= f;
	//}
	//univec[255].x = univec[255].y = univec[255].z = 0;
	for(i=0;i<256;i++)
	{
		iunivec[i][0] = (short)(univec[i].x*4096.0);
		iunivec[i][1] = (short)(univec[i].y*4096.0);
		iunivec[i][2] = (short)(univec[i].z*4096.0);
		iunivec[i][3] = 4096;
	}
	ucossininit();

	memset(mixn,0,sizeof(mixn));

		//Initialize hash table for getkv6()
	memset(khashead,-1,sizeof(khashead));
	if (!(khashbuf = (char *)malloc(KHASHINITSIZE))) return(-1);
	khashsiz = KHASHINITSIZE;

	vx5.anginc = 1; //Higher=faster (1:full,2:half)
	vx5.sideshademode = 0; setsideshades(0,0,0,0,0,0);
	vx5.mipscandist = 128;
	vx5.maxscandist = 256; //must be <= 2047
	vx5.colfunc = curcolfunc; //This prevents omission bugs from crashing voxlap5
	vx5.curcol = 0x80804c33;
	vx5.currad = 8;
	vx5.curhei = 0;
	vx5.curpow = 2.0;
	vx5.amount = 0x70707;
	vx5.pic = 0;
	vx5.cliphitnum = 0;
	vx5.xplanemin = 0;
	vx5.xplanemax = 0x7fffffff;
	vx5.flstnum = 0;
	vx5.lightmode = 0;
	vx5.numlights = 0;
	vx5.kv6mipfactor = 96;
	vx5.kv6col = 0x808080;
	vx5.vxlmipuse = 1;
	vx5.fogcol = -1;
	vx5.fallcheck = 0;

	gmipnum = 0;

	return(0);
}

#if 0 
	long i, j, k, l;
	char *v;

	j = k = l = 0;
	for(i=0;i<VSID*VSID;i++)
	{
		for(v=sptr[i];v[0];v+=v[0]*4) { j++; k += v[2]-v[1]+1; l += v[0]-1; }
		k += v[2]-v[1]+1; l += v[2]-v[1]+1;
	}

	printf("VOXLAP5 programmed by Ken Silverman (www.advsys.net/ken)\n");
	printf("Please DO NOT DISTRIBUTE! If this leaks, I will not be happy.\n\n");
	//printf("This copy licensed to:  \n\n");
	printf("Memory statistics upon exit: (all numbers in bytes)");
	printf("\n");
	if (screen) printf("   screen: %8ld\n",imageSize);
	printf("    radar: %8ld\n",max((((MAXXDIM*MAXYDIM*27)>>1)+7)&~7,(VSID+4)*3*SCPITCH*4+8));
	printf("  bacsptr: %8ld\n",sizeof(bacsptr));
	printf("     sptr: %8ld\n",(VSID*VSID)<<2);
	printf("     vbuf: %8ld(%8ld)\n",(j+VSID*VSID+l)<<2,VOXSIZ);
	printf("     vbit: %8ld(%8ld)\n",VOXSIZ>>5);
	printf("\n");
	printf("vbuf head: %8ld\n",(j+VSID*VSID)<<2);
	printf("vbuf cols: %8ld\n",l<<2);
	printf("     fcol: %8ld\n",k<<2);
	printf("     ccol: %8ld\n",(l-k)<<2);
	printf("\n");
	printf("%.2f bytes/column\n",(float)((j+VSID*VSID+l)<<2)/(float)(VSID*VSID));
#endif
