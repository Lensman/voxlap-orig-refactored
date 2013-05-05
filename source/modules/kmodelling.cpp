#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <xmmintrin.h>

#include "../../include/sysmain.h"
#include "../../include/porthacks.h"
#include "../../include/voxlap5.h"
#include "../../include/kglobals.h"
#include "../../include/cpu_detect.h"
#include "../../include/ksnippits.h"

extern long hpngcolfunc (point3d *);
extern char umost[], dmost[];
extern long templongbuf[];
extern long *radar, *radarmem;
extern char tbuf[];
extern long tbuf2[];
extern long *vbuf, *vbit;
	//Look tables for expandbitstack256:
long xbsceil[32], xbsflor[32]; // Prefilled by initvoxlap

#define SETSPHMAXRAD 256
double logint[SETSPHMAXRAD];
float tempfloatbuf[SETSPHMAXRAD];
long factr[SETSPHMAXRAD][2];

void expandrle (long x, long y, long *uind)
{
	unsigned long i;
	char *v;

	if ((x|y)&(~(VSID-1))) { uind[0] = 0; uind[1] = MAXZDIM; return; }

	v = sptr[y*VSID+x]; uind[0] = v[1]; i = 2;
	while (v[0])
	{
		v += v[0]*4; if (v[3] >= v[1]) continue;
		uind[i-1] = v[3]; uind[i] = v[1]; i += 2;
	}
	uind[i-1] = MAXZDIM;
}

	//Inputs:  n0[<=MAXZDIM]: rle buffer of column to compress
	//         n1-4[<=MAXZDIM]: neighboring rle buffers
	//         top,bot,top,bot,... (ends when bot == MAXZDIM)
	//         px,py: takes color from original column (before modification)
	//            If originally unexposed, calls vx5.colfunc(.)
	//Outputs: cbuf[MAXCSIZ]: compressed output buffer
	//Returns: n: length of compressed buffer (in bytes)
long compilerle (long *n0, long *n1, long *n2, long *n3, long *n4, char *cbuf, long px, long py)
{
	long i, ia, ze, zend, onext, dacnt, n, *ic;
	lpoint3d p;
	char *v;

	p.x = px; p.y = py;

		//Generate pointers to color slabs in this format:
		//   0:z0,  1:z1,  2:(pointer to z0's color)-z0
	v = sptr[py*VSID+px]; ic = tbuf2;
	while (1)
	{
		ia = v[1]; p.z = v[2];
		ic[0] = ia; ic[1] = p.z+1; ic[2] = ((long)v)-(ia<<2)+4; ic += 3;
		i = v[0]; if (!i) break;
		v += i*4; ze = v[3];
		ic[0] = ze+p.z-ia-i+2; ic[1] = ze; ic[2] = ((long)v)-(ze<<2); ic += 3;
	}
	ic[0] = MAXZDIM; ic[1] = MAXZDIM;

	p.z = n0[0]; cbuf[1] = n0[0];
	ze = n0[1]; cbuf[2] = ze-1;
	cbuf[3] = 0;
	i = onext = 0; ic = tbuf2; ia = 15; n = 4;
	if (ze != MAXZDIM) zend = ze-1; else zend = -1;
	while (1)
	{
		dacnt = 0;
		while (1)
		{
			do
			{
				while (p.z >= ic[1]) ic += 3;
				if (p.z >= ic[0]) *(long *)&cbuf[n] = *(long *)(ic[2]+(p.z<<2));
								 else *(long *)&cbuf[n] = vx5.colfunc(&p);
				n += 4; p.z++; if (p.z >= ze) goto rlendit2;
				while (p.z >= n1[0]) { n1++; ia ^= 1; }
				while (p.z >= n2[0]) { n2++; ia ^= 2; }
				while (p.z >= n3[0]) { n3++; ia ^= 4; }
				while (p.z >= n4[0]) { n4++; ia ^= 8; }
			} while ((ia) || (p.z == zend));

			if (!dacnt) { cbuf[onext+2] = p.z-1; dacnt = 1; }
			else
			{
				cbuf[onext] = ((n-onext)>>2); onext = n;
				cbuf[n+1] = p.z; cbuf[n+2] = p.z-1; cbuf[n+3] = p.z; n += 4;
			}

			if ((n1[0] < n2[0]) && (n1[0] < n3[0]) && (n1[0] < n4[0]))
				{ if (n1[0] >= ze) { p.z = ze-1; } else { p.z = *n1++; ia ^= 1; } }
			else if ((n2[0] < n3[0]) && (n2[0] < n4[0]))
				{ if (n2[0] >= ze) { p.z = ze-1; } else { p.z = *n2++; ia ^= 2; } }
			else if (n3[0] < n4[0])
				{ if (n3[0] >= ze) { p.z = ze-1; } else { p.z = *n3++; ia ^= 4; } }
			else
				{ if (n4[0] >= ze) { p.z = ze-1; } else { p.z = *n4++; ia ^= 8; } }

			if (p.z == MAXZDIM-1) goto rlenditall;
		}
rlendit2:;
		if (ze >= MAXZDIM) break;

		i += 2;
		cbuf[onext] = ((n-onext)>>2); onext = n;
		p.z = n0[i]; cbuf[n+1] = n0[i]; cbuf[n+3] = ze;
		ze = n0[i+1]; cbuf[n+2] = ze-1;
		n += 4;
	}
rlenditall:;
	cbuf[onext] = 0;
	return(n);
}

	//Delete everything on b2() in y0<=y<y1
void delslab (long *b2, const long y0, long y1)
{
	long i, j, z;

	if (y1 >= MAXZDIM) y1 = MAXZDIM-1;
	if ((y0 >= y1) || (!b2)) return;
	for(z=0;y0>=b2[z+1];z+=2);
	if (y0 > b2[z])
	{
		if (y1 < b2[z+1])
		{
			for(i=z;b2[i+1]<MAXZDIM;i+=2);
			while (i > z) { b2[i+3] = b2[i+1]; b2[i+2] = b2[i]; i -= 2; }
			b2[z+3] = b2[z+1]; b2[z+1] = y0; b2[z+2] = y1; return;
		}
		b2[z+1] = y0; z += 2;
	}
	if (y1 >= b2[z+1])
	{
		for(i=z+2;y1>=b2[i+1];i+=2);
		j = z-i; b2[z] = b2[i]; b2[z+1] = b2[i+1];
		while (b2[i+1] < MAXZDIM)
			{ i += 2; b2[i+j] = b2[i]; b2[i+j+1] = b2[i+1]; }
	}
	if (y1 > b2[z]) b2[z] = y1;
}
	//Insert everything on b2() in y0<=y<y1
void insslab (long *b2, const long y0, long y1)
{
	long i, j, z;

	if ((y0 >= y1) || (!b2)) return;
	for(z=0;y0>b2[z+1];z+=2);
	if (y1 < b2[z])
	{
		for(i=z;b2[i+1]<MAXZDIM;i+=2);
		do { b2[i+3] = b2[i+1]; b2[i+2] = b2[i]; i -= 2; } while (i >= z);
		b2[z+1] = y1; b2[z] = y0; return;
	}
	if (y0 < b2[z]) b2[z] = y0;
	if ((y1 >= b2[z+2]) && (b2[z+1] < MAXZDIM))
	{
		for(i=z+2;(y1 >= b2[i+2]) && (b2[i+1] < MAXZDIM);i+=2);
		j = z-i; b2[z+1] = b2[i+1];
		while (b2[i+1] < MAXZDIM)
			{ i += 2; b2[i+j] = b2[i]; b2[i+j+1] = b2[i+1]; }
	}
	if (y1 > b2[z+1]) b2[z+1] = y1;
}

#ifdef _MSC_VER

void expandbit256 (const void * __restrict s, void *d)
{
 #if defined(__GNUC__) && defined(__i386__) && ! defined(NOASM)
	__asm__ __volatile__
	(
		".intel_syntax noprefix\n"
		"mov	ecx, 32\n"   //current bit index
		"xor	edx, edx\n"  //value of current 32-bit bits
		"jmp	short 0f\n"
	"1:\n" //begit
		"lea	esi, [esi+eax*4]\n"
		"movzx	eax, byte ptr [esi+3]\n"
		"sub	eax, ecx\n"                //xor mask [eax] for ceiling begins
		"jl	short 3f\n"
	"2:\n" //xdoc
		"mov	[edi], edx\n"
		"add	edi, 4\n"
		"mov	edx, -1\n"
		"add	ecx, 32\n"
		"sub	eax, 32\n"
		"jge	short 2b\n"
	"3:\n" //xskpc
		".att_syntax prefix\n"
		"andl	%c[ceil]+0x80(,%%eax,4), %%edx\n" //~(-1<<eax)
		".intel_syntax noprefix\n"

	"0:\n" //in2it
		"movzx	eax, byte ptr [esi+1]\n"
		"sub	eax, ecx\n"                //xor mask [eax] for floor begins
		"jl	short 5f\n"
	"4:\n" //xdof
		"mov	[edi], edx\n"
		"add	edi, 4\n"
		"xor	edx, edx\n"
		"add	ecx, 32\n"
		"sub	eax, 32\n"
		"jge	short 4b\n"
	"5:\n" //xskpf
		".att_syntax prefix\n"
		"orl	%c[floor]+0x80(,%%eax,4), %%edx\n" //(-1<<eax)
		".intel_syntax noprefix\n"

		"movzx	eax, byte ptr [esi]\n"
		"test	eax, eax\n"
		"jnz	short 1b\n"
		"sub	ecx, 256\n"              //finish writing buffer to [edi]
		"jg	short 7f\n"
	"6:\n" //xdoe
		"mov	[edi], edx\n"
		"add	edi, 4\n"
		"mov	edx, -1\n"
		"add	ecx, 32\n"
		"jle	short 6b\n"
	"7:\n" //xskpe
		".att_syntax prefix\n"
		:
		: "S" (s), "D" (d), [ceil] "p" (xbsceil), [floor] "p" (xbsflor)
		: "cc", "eax", "ecx", "edx"
	);
	#elif defined(_MSC_VER) && defined(__i386__) && ! defined(NOASM) //msvc inline asm
	_asm
	{
		push	esi
		push	edi
		mov	esi, s
		mov	edi, d
		mov	ecx, 32   //current bit index
		xor	edx, edx  //value of current 32-bit bits
		jmp	short in2it
	begit:
		lea	esi, [esi+eax*4]
		movzx	eax, byte ptr [esi+3]
		sub	eax, ecx                //xor mask [eax] for ceiling begins
		jl	short xskpc
	xdoc:
		mov	[edi], edx
		add	edi, 4
		mov	edx, -1
		add	ecx, 32
		sub	eax, 32
		jge	short xdoc
	xskpc:
		and	edx, xbsceil[eax*4+128] //~(-1<<eax); xor mask [eax] for ceiling ends
	in2it:
		movzx	eax, byte ptr [esi+1]
		sub	eax, ecx                //xor mask [eax] for floor begins
		jl	short xskpf
	xdof:
		mov	[edi], edx
		add	edi, 4
		xor	edx, edx
		add	ecx, 32
		sub	eax, 32
		jge	short xdof
	xskpf:
		or	edx, xbsflor[eax*4+128] //(-1<<eax); xor mask [eax] for floor ends
		movzx	eax, byte ptr [esi]
		test	eax, eax
		jnz	short begit
		sub	ecx, 256                //finish writing buffer to [edi]
		jg	short xskpe
	xdoe:
		mov	[edi], edx
		add	edi, 4
		mov	edx, -1
		add	ecx, 32
		jle	short xdoe
	xskpe:
		pop	edi
		pop	esi
	}
	#else // C Default
	int32_t eax;
	int32_t ecx = 32; //current bit index
	uint32_t edx = 0; //value of current 32-bit bits

	goto in2it;

	while (eax != 0)
	{
		s += eax * 4;
		eax = ((uint8_t*)s)[3];

		if ((eax -= ecx) >= 0) //xor mask [eax] for ceiling begins
		{
			do
			{
				*(uint32_t*)d = edx;
				d += 4;
				edx = -1;
				ecx += 32;
			}
			while ((eax -= 32) >= 0);
		}

		edx &= xbsceil[32+eax];
		//no jump

	in2it:
		eax = ((uint8_t*)s)[1];

		if ((eax -= ecx) > 0) //xor mask [eax] for floor begins
		{
			do
			{
				*(uint32_t*)d = edx;
				d += 4;
				edx = 0;
				ecx += 32;
			}
			while ((eax -= 32) >= 0);
		}

		edx |= xbsflor[32+eax];
		eax = *(uint8_t*)s;
	}

	if ((ecx -= 256) <= 0)
	{
		do
		{
			*(uint32_t*)d = edx;
			d += 4;
			edx = -1;
		}
		while ((ecx += 32) <= 0);
	}
	#endif
}

#endif
void expandbitstack ( const long &x, const long &y, int64_t *bind)
{
	if ((x|y)&(~(VSID-1))) { clearbuf((void *)bind,8,0L); return; }
	expandbit256(sptr[y*VSID+x],(void *)bind);
}

void expandstack (const long &x, const long &y, long *uind)
{
	long z, topz;
	char *v, *v2;

	if ((x|y)&(~(VSID-1))) { clearbuf((void *)uind,MAXZDIM,0); return; }

		//Expands compiled voxel info to 32-bit uind[?]
	v = sptr[y*VSID+x]; z = 0;
	while (1)
	{
		while (z < v[1]) { uind[z] = -1; z++; }
		while (z <= v[2]) { uind[z] = (*(long *)&v[(z-v[1])*4+4]); z++; }
		v2 = &v[(v[2]-v[1]+1)*4+4];

		if (!v[0]) break;
		v += v[0]*4;

		topz = v[3]+(((long)v2-(long)v)>>2);
		while (z < topz) { uind[z] = -2; z++; }
		while (z < v[3]) { uind[z] = *(long *)v2; z++; v2 += 4; }
	}
	while (z < MAXZDIM) { uind[z] = -2; z++; }
}
	// Inputs: uind[MAXZDIM]: uncompressed 32-bit color buffer (-1: air)
	//         nind?[MAXZDIM]: neighbor buf:
	//            -2: unexposed solid
	//            -1: air
	//    0-16777215: exposed solid (color)
	//         px,py: parameters for setting unexposed voxel colors
	//Outputs: cbuf[MAXCSIZ]: compressed output buffer
	//Returns: n: length of compressed buffer (in bytes)
long compilestack (long *uind, const long *n0, long *n1, long *n2, long *n3, char *cbuf, const long px, const long py)
{
	long oz, onext, n, cp2, cp1, cp0, rp1, rp0;
	lpoint3d p;

	p.x = px; p.y = py;

		//Do top slab (sky)
	oz = -1;
	p.z = -1; while (uind[p.z+1] == -1) p.z++;
	onext = 0;
	cbuf[1] = p.z+1;
	cbuf[2] = p.z+1;
	cbuf[3] = 0;  //Top z0 (filler, not used yet)
	n = 4;
	cp1 = 1; cp0 = 0;
	rp1 = -1; rp0 = -1;

	do
	{
			//cp2 = state at p.z-1 (0 = air, 1 = next2air, 2 = solid)
			//cp1 = state at p.z   (0 = air, 1 = next2air, 2 = solid)
			//cp0 = state at p.z+1 (0 = air, 1 = next2air, 2 = solid)
		cp2 = cp1; cp1 = cp0; cp0 = 2;
		if (p.z < MAXZDIM-2)  //Bottom must be solid!
		{
			if (uind[p.z+1] == -1)
				cp0 = 0;
			else if ((n0[p.z+1] == -1) || (n1[p.z+1] == -1) ||
						(n2[p.z+1] == -1) || (n3[p.z+1] == -1))
				cp0 = 1;
		}

			//Add slab
		if (cp1 != rp0)
		{
			if ((!cp1) && (rp0 > 0)) { oz = p.z; }
			else if ((rp0 < cp1) && (rp0 < rp1))
			{
				if (oz < 0) oz = p.z;
				cbuf[onext] = ((n-onext)>>2); onext = n;
				cbuf[n+1] = p.z;
				cbuf[n+2] = p.z-1;
				cbuf[n+3] = oz;
				n += 4; oz = -1;
			}
			rp1 = rp0; rp0 = cp1;
		}

			//Add color
		if ((cp1 == 1) || ((cp1 == 2) && ((!cp0) || (!cp2))))
		{
			if (cbuf[onext+2] == p.z-1) cbuf[onext+2] = p.z;
			if (uind[p.z] == -2) *(long *)&cbuf[n] = vx5.colfunc(&p);
								 else *(long *)&cbuf[n] = uind[p.z];
			n += 4;
		}

		p.z++;
	} while (p.z < MAXZDIM);
	cbuf[onext] = 0;
	return(n);
}

static long scx0, scx1, scox0, scox1, scoox0, scoox1;
static long scex0, scex1, sceox0, sceox1, scoy = 0x80000000, *scoym3;


void scumline ()
{
	long i, j, k, x, y, x0, x1, *mptr, *uptr;
	char *v;

	x0 = min(scox0-1,min(scx0,scoox0)); scoox0 = scox0; scox0 = scx0;
	x1 = max(scox1+1,max(scx1,scoox1)); scoox1 = scox1; scox1 = scx1;

	uptr = &scoym3[SCPITCH]; if (uptr == &radar[SCPITCH*9]) uptr = &radar[SCPITCH*6];
	mptr = &uptr[SCPITCH];   if (mptr == &radar[SCPITCH*9]) mptr = &radar[SCPITCH*6];
	long SCPITCH_MUL3 = SCPITCH * 3;
	if ((x1 < sceox0) || (x0 > sceox1))
	{
		for(x=x0;x<=x1;x++) expandstack(x,scoy-2,&uptr[x*SCPITCH_MUL3]);
	}
	else
	{
		for(x=x0;x<sceox0;x++) expandstack(x,scoy-2,&uptr[x*SCPITCH_MUL3]);
		for(x=x1;x>sceox1;x--) expandstack(x,scoy-2,&uptr[x*SCPITCH_MUL3]);
	}

	if ((scex1|x1) >= 0)
	{
		for(x=x1+2;x<scex0;x++) expandstack(x,scoy-1,&mptr[x*SCPITCH_MUL3]);
		for(x=x0-2;x>scex1;x--) expandstack(x,scoy-1,&mptr[x*SCPITCH_MUL3]);
	}
	if ((x1+1 < scex0) || (x0-1 > scex1))
	{
		for(x=x0-1;x<=x1+1;x++) expandstack(x,scoy-1,&mptr[x*SCPITCH_MUL3]);
	}
	else
	{
		for(x=x0-1;x<scex0;x++) expandstack(x,scoy-1,&mptr[x*SCPITCH_MUL3]);
		for(x=x1+1;x>scex1;x--) expandstack(x,scoy-1,&mptr[x*SCPITCH_MUL3]);
	}
	sceox0 = min(x0-1,scex0);
	sceox1 = max(x1+1,scex1);

	if ((x1 < scx0) || (x0 > scx1))
	{
		for(x=x0;x<=x1;x++) expandstack(x,scoy,&scoym3[x*SCPITCH_MUL3]);
	}
	else
	{
		for(x=x0;x<scx0;x++) expandstack(x,scoy,&scoym3[x*SCPITCH_MUL3]);
		for(x=x1;x>scx1;x--) expandstack(x,scoy,&scoym3[x*SCPITCH_MUL3]);
	}
	scex0 = x0;
	scex1 = x1;

	y = scoy-1; if (y&(~(VSID-1))) return;
	if (x0 < 0) x0 = 0;
	if (x1 >= VSID) x1 = VSID-1;
	i = y*VSID+x0; k = x0*SCPITCH_MUL3;
	for(x=x0;x<=x1;x++,i++,k+=SCPITCH_MUL3)
	{
		v = sptr[i]; vx5.globalmass += v[1];
		while (v[0]) { v += v[0]*4; vx5.globalmass += v[1]-v[3]; }

			//De-allocate column (x,y)
		voxdealloc(sptr[i]);

		j = compilestack(&mptr[k],&mptr[k-SCPITCH_MUL3],&mptr[k+SCPITCH_MUL3],&uptr[k],&scoym3[k],
							  tbuf,x,y);

			//Allocate & copy to new column (x,y)
		sptr[i] = v = voxalloc(j); copybuf((void *)tbuf,(void *)v,j>>2);

		vx5.globalmass -= v[1];
		while (v[0]) { v += v[0]*4; vx5.globalmass += v[3]-v[1]; }
	}
}

	//x: x on voxel map
	//y: y on voxel map
	//z0: highest z on column
	//z1: lowest z(+1) on column
	//nbuf: buffer of color data from nbuf[z0] to nbuf[z1-1];
	//           -3: don't modify voxel
	//           -2: solid voxel (unexposed): to be calculated in compilestack
	//           -1: write air voxel
	//   0-16777215: write solid voxel (exposed)
void scum (const long &x, const long &y, const long &z0, const long &z1, long *nbuf)
{
	long z, *mptr;

	if ((x|y)&(~(VSID-1))) return;

	if (y != scoy)
	{
		if (scoy >= 0)
		{
			scumline();
			while (scoy < y-1)
			{
				scx0 = 0x7fffffff; scx1 = 0x80000000;
				scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
				scumline();
			}
			scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
		}
		else
		{
			scoox0 = scox0 = 0x7fffffff;
			sceox0 = scex0 = x+1;
			sceox1 = scex1 = x;
			scoy = y; scoym3 = &radar[SCPITCH*6];
		}
		scx0 = x;
	}
	else
	{
		while (scx1 < x-1) { scx1++; expandstack(scx1,y,&scoym3[scx1*SCPITCH*3]); }
	}

	mptr = &scoym3[x*SCPITCH*3]; scx1 = x; expandstack(x,y,mptr);

		//Modify column (x,y):
	if (nbuf[MAXZDIM-1] == -1) nbuf[MAXZDIM-1] = -2; //Bottom voxel must be solid
	for(z=z0;z<z1;z++)
		if (nbuf[z] != -3) mptr[z] = nbuf[z];
}

void scumfinish ()
{
	long i;

	if (scoy == 0x80000000) return;
	for(i=2;i;i--)
	{
		scumline(); scx0 = 0x7fffffff; scx1 = 0x80000000;
		scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
	}
	scumline(); scoy = 0x80000000;
}

/**
	//Example of how to use this code:
	//vx5.colfunc = curcolfunc; //0<x0<x1<VSID, 0<y0<y1<VSID, 0<z0<z1<256,
	//clearbuf((void *)&templongbuf[z0],z1-z0,-1); //Ex: set all voxels to air
	//for(y=y0;y<y1;y++) //MUST iterate x&y in this order, but can skip around
	//   for(x=x0;x<x1;x++)
	//      if (rand()&8) scum(x,y,z0,z1,templongbuf));
	//scumfinish(); //MUST call this when done!
*/
void scum2line ()
{
	long i, j, k, x, y, x0, x1, *mptr, *uptr;
	char *v;

	x0 = min(scox0-1,min(scx0,scoox0)); scoox0 = scox0; scox0 = scx0;
	x1 = max(scox1+1,max(scx1,scoox1)); scoox1 = scox1; scox1 = scx1;

	uptr = &scoym3[SCPITCH]; if (uptr == &radar[SCPITCH*9]) uptr = &radar[SCPITCH*6];
	mptr = &uptr[SCPITCH];   if (mptr == &radar[SCPITCH*9]) mptr = &radar[SCPITCH*6];

	if ((x1 < sceox0) || (x0 > sceox1))
	{
		for(x=x0;x<=x1;x++) expandrle(x,scoy-2,&uptr[x*SCPITCH*3]);
	}
	else
	{
		for(x=x0;x<sceox0;x++) expandrle(x,scoy-2,&uptr[x*SCPITCH*3]);
		for(x=x1;x>sceox1;x--) expandrle(x,scoy-2,&uptr[x*SCPITCH*3]);
	}

	if ((scex1|x1) >= 0)
	{
		for(x=x1+2;x<scex0;x++) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
		for(x=x0-2;x>scex1;x--) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
	}
	if ((x1+1 < scex0) || (x0-1 > scex1))
	{
		for(x=x0-1;x<=x1+1;x++) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
	}
	else
	{
		for(x=x0-1;x<scex0;x++) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
		for(x=x1+1;x>scex1;x--) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
	}
	sceox0 = min(x0-1,scex0);
	sceox1 = max(x1+1,scex1);

	if ((x1 < scx0) || (x0 > scx1))
	{
		for(x=x0;x<=x1;x++) expandrle(x,scoy,&scoym3[x*SCPITCH*3]);
	}
	else
	{
		for(x=x0;x<scx0;x++) expandrle(x,scoy,&scoym3[x*SCPITCH*3]);
		for(x=x1;x>scx1;x--) expandrle(x,scoy,&scoym3[x*SCPITCH*3]);
	}
	scex0 = x0;
	scex1 = x1;

	y = scoy-1; if (y&(~(VSID-1))) return;
	if (x0 < 0) x0 = 0;
	if (x1 >= VSID) x1 = VSID-1;
	i = y*VSID+x0; k = x0*SCPITCH*3;
	for(x=x0;x<=x1;x++,i++,k+=SCPITCH*3)
	{
		j = compilerle(&mptr[k],&mptr[k-SCPITCH*3],&mptr[k+SCPITCH*3],&uptr[k],&scoym3[k],
							tbuf,x,y);

		v = sptr[i]; vx5.globalmass += v[1];
		while (v[0]) { v += v[0]*4; vx5.globalmass += v[1]-v[3]; }

			//De-allocate column (x,y)  Note: Must be AFTER compilerle!
		voxdealloc(sptr[i]);

			//Allocate & copy to new column (x,y)
		sptr[i] = v = voxalloc(j); copybuf((void *)tbuf,(void *)v,j>>2);

		vx5.globalmass -= v[1];
		while (v[0]) { v += v[0]*4; vx5.globalmass += v[3]-v[1]; }
	}
}

/**
 *  @param x: x on voxel map
 *  @param y: y on voxel map
 *  @return pointer to rle column (x,y)
 */
long *scum2 (const long &x, const long &y)
{
	long *mptr;

	if ((x|y)&(~(VSID-1))) return(0);

	if (y != scoy)
	{
		if (scoy >= 0)
		{
			scum2line();
			while (scoy < y-1)
			{
				scx0 = 0x7fffffff; scx1 = 0x80000000;
				scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
				scum2line();
			}
			scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
		}
		else
		{
			scoox0 = scox0 = 0x7fffffff;
			sceox0 = scex0 = x+1;
			sceox1 = scex1 = x;
			scoy = y; scoym3 = &radar[SCPITCH*6];
		}
		scx0 = x;
	}
	else
	{
		while (scx1 < x-1) { scx1++; expandrle(scx1,y,&scoym3[scx1*SCPITCH*3]); }
	}

	mptr = &scoym3[x*SCPITCH*3]; scx1 = x; expandrle(x,y,mptr);
	return(mptr);
}

void scum2finish ()
{
	long i;

	if (scoy == 0x80000000) return;
	for(i=2;i;i--)
	{
		scum2line(); scx0 = 0x7fffffff; scx1 = 0x80000000;
		scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
	}
	scum2line(); scoy = 0x80000000;
}

	//WARNING! Make sure to set vx5.colfunc before calling this function!
	//This function is here for simplicity only - it is NOT optimal.
	//
	//   -1: set air
	//   -2: use vx5.colfunc
void setcube (const long &px, const long &py, const long &pz, const long &col)
{
	long bakcol, (*bakcolfunc)(lpoint3d *), *lptr;
    bbox_t t = *(bbox_t*)(long*)&vx5.minx;
    t.bmin.x = px;
    t.bmin.y = py;
    t.bmin.z = pz;
    t.bmax.x = px+1;
    t.bmax.y = py+1;
    t.bmax.z = pz+1;
	if ((unsigned long)pz >= MAXZDIM) return;
	if ((unsigned long)col >= (unsigned long)0xfffffffe) //-1 or -2
	{
		lptr = scum2(px,py);
		if (col == -1) delslab(lptr,pz,pz+1); else insslab(lptr,pz,pz+1);
		scum2finish();
		updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,col);
		return;
	}

	bakcol = getcube(px,py,pz);
	if (bakcol == 1) return; //Unexposed solid
	if (bakcol != 0) //Not 0 (air)
		*(long *)bakcol = col;
	else
	{
		bakcolfunc = vx5.colfunc; bakcol = vx5.curcol;
		vx5.colfunc = curcolfunc; vx5.curcol = col;
		insslab(scum2(px,py),pz,pz+1); scum2finish();
		vx5.colfunc = bakcolfunc; vx5.curcol = bakcol;
	}
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,0);
}

void setsphere (const lpoint3d *hit, long hitrad, const long &dacol)
{
	void (*modslab)(long *, long, long);
	long i, x, y, xs, ys, zs, xe, ye, ze, sq;
	float f, ff;

	xs = max(hit->x-hitrad,0); xe = min(hit->x+hitrad,VSID-1);
	ys = max(hit->y-hitrad,0); ye = min(hit->y+hitrad,VSID-1);
	zs = max(hit->z-hitrad,0); ze = min(hit->z+hitrad,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;

	if (vx5.colfunc == sphcolfunc)
	{
		vx5.cen = hit->x+hit->y+hit->z;
		vx5.daf = 1.f/(hitrad*sqrt(3.f));
	}

	if (hitrad >= SETSPHMAXRAD-1) hitrad = SETSPHMAXRAD-2;
	if (dacol == -1) modslab = delslab; else modslab = insslab;

	tempfloatbuf[0] = 0.0f;
#if 0
		//Totally unoptimized
	for(i=1;i<=hitrad;i++) tempfloatbuf[i] = pow((float)i,vx5.curpow);
#else
	tempfloatbuf[1] = 1.0f;
	for(i=2;i<=hitrad;i++)
	{
		if (!factr[i][0]) tempfloatbuf[i] = exp(logint[i]*vx5.curpow);
		else tempfloatbuf[i] = tempfloatbuf[factr[i][0]]*tempfloatbuf[factr[i][1]];
	}
#endif
	*(long *)&tempfloatbuf[hitrad+1] = 0x7f7fffff; //3.4028235e38f; //Highest float

	sq = 0; //pow(fabs(x-hit->x),vx5.curpow) + "y + "z < pow(vx5.currad,vx5.curpow)
	for(y=ys;y<=ye;y++)
	{
		ff = tempfloatbuf[hitrad]-tempfloatbuf[labs(y-hit->y)];
		if (*(long *)&ff <= 0) continue;
		for(x=xs;x<=xe;x++)
		{
			f = ff-tempfloatbuf[labs(x-hit->x)]; if (*(long *)&f <= 0) continue;
			while (*(long *)&tempfloatbuf[sq] <  *(long *)&f) sq++;
			while (*(long *)&tempfloatbuf[sq] >= *(long *)&f) sq--;
			modslab(scum2(x,y),max(hit->z-sq,zs),min(hit->z+sq+1,ze));
		}
	}
	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
}

void setellipsoid (lpoint3d *hit, lpoint3d *hit2, long hitrad, long dacol, long bakit)
{
	void (*modslab)(long *, long, long);
	long x, y, xs, ys, zs, xe, ye, ze;
	float a, b, c, d, e, f, g, h, r, t, u, Za, Zb, fx0, fy0, fz0, fx1, fy1, fz1;

	xs = min(hit->x,hit2->x)-hitrad; xs = max(xs,0);
	ys = min(hit->y,hit2->y)-hitrad; ys = max(ys,0);
	zs = min(hit->z,hit2->z)-hitrad; zs = max(zs,0);
	xe = max(hit->x,hit2->x)+hitrad; xe = min(xe,VSID-1);
	ye = max(hit->y,hit2->y)+hitrad; ye = min(ye,VSID-1);
	ze = max(hit->z,hit2->z)+hitrad; ze = min(ze,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze))
		{ if (bakit) voxbackup(xs,ys,xs,ys,bakit); return; }

	fx0 = (float)hit->x; fy0 = (float)hit->y; fz0 = (float)hit->z;
	fx1 = (float)hit2->x; fy1 = (float)hit2->y; fz1 = (float)hit2->z;

	r = (fx1-fx0)*(fx1-fx0) + (fy1-fy0)*(fy1-fy0) + (fz1-fz0)*(fz1-fz0);
	r = sqrt((float)hitrad*(float)hitrad + r*.25);
	c = fz0*fz0 - fz1*fz1; d = r*r*-4; e = d*4;
	f = c*c + fz1*fz1 * e; g = c + c; h = (fz1-fz0)*2; c = c*h - fz1*e;
	Za = -h*h - e; if (Za <= 0) { if (bakit) voxbackup(xs,ys,xs,ys,bakit); return; }
	u = 1 / Za;

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;

	if (dacol == -1) modslab = delslab; else modslab = insslab;

	if (bakit) voxbackup(xs,ys,xe+1,ye+1,bakit);

	for(y=ys;y<=ye;y++)
		for(x=xs;x<=xe;x++)
		{
			a = (x-fx0)*(x-fx0) + (y-fy0)*(y-fy0);
			b = (x-fx1)*(x-fx1) + (y-fy1)*(y-fy1);
			t = a-b+d; Zb = t*h + c;
			t = ((t+g)*t + b*e + f)*Za + Zb*Zb; if (t <= 0) continue;
			t = sqrt(t);
			ftol((Zb - t)*u,&zs); if (zs < 0) zs = 0;
			ftol((Zb + t)*u,&ze); if (ze > MAXZDIM) ze = MAXZDIM;
			modslab(scum2(x,y),zs,ze);
		}
	scum2finish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
}

	//Draws a cylinder, given: 2 points, a radius, and a color
	//Code mostly optimized - original code from CYLINDER.BAS:drawcylinder
void setcylinder (lpoint3d *p0, lpoint3d *p1, long cr, long dacol, long bakit)
{
	void (*modslab)(long *, long, long);

	float t, ax, ay, az, bx, by, bz, cx, cy, cz, ux, uy, uz, vx, vy, vz;
	float Za, Zb, Zc, tcr, xxyy, rcz, rZa;
	float fx, fxi, xof, vx0, vy0, vz0, vz0i, vxo, vyo, vzo;
	long i, j, ix, iy, ix0, ix1, iz0, iz1, minx, maxx, miny, maxy;
	long x0, y0, z0, x1, y1, z1;

		//Map generic cylinder into unit space:  (0,0,0), (0,0,1), cr = 1
		//   x*x + y*y < 1, z >= 0, z < 1
	if (p0->z > p1->z)
	{
		x0 = p1->x; y0 = p1->y; z0 = p1->z;
		x1 = p0->x; y1 = p0->y; z1 = p0->z;
	}
	else
	{
		x0 = p0->x; y0 = p0->y; z0 = p0->z;
		x1 = p1->x; y1 = p1->y; z1 = p1->z;
	}

	xxyy = (float)((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
	t = xxyy + (float)(z1-z0)*(z1-z0);
	if ((t == 0) || (cr == 0))
	{
		vx5.minx = x0; vx5.maxx = x0+1;
		vx5.miny = y0; vx5.maxy = y0+1;
		vx5.minz = z0; vx5.maxz = z0+1;
		if (bakit) voxbackup(x0,y0,x0,y0,bakit);
		return;
	}
	t = 1 / t; cx = ((float)(x1-x0))*t; cy = ((float)(y1-y0))*t; cz = ((float)(z1-z0))*t;
	t = sqrt(t); ux = ((float)(x1-x0))*t; uy = ((float)(y1-y0))*t; uz = ((float)(z1-z0))*t;

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;

	if (dacol == -1) modslab = delslab; else modslab = insslab;

	if (xxyy == 0)
	{
		iz0 = max(z0,0); iz1 = min(z1,MAXZDIM);
		minx = max(x0-cr,0); maxx = min(x0+cr,VSID-1);
		miny = max(y0-cr,0); maxy = min(y0+cr,VSID-1);

		vx5.minx = minx; vx5.maxx = maxx+1;
		vx5.miny = miny; vx5.maxy = maxy+1;
		vx5.minz = iz0; vx5.maxz = iz1;
		if (bakit) voxbackup(minx,miny,maxx+1,maxy+1,bakit);

		j = cr*cr;
		for(iy=miny;iy<=maxy;iy++)
		{
			i = j-(iy-y0)*(iy-y0);
			for(ix=minx;ix<=maxx;ix++)
				if ((ix-x0)*(ix-x0) < i) modslab(scum2(ix,iy),iz0,iz1);
		}
		scum2finish();
		updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
		return;
	}

	if (x0 < x1) { minx = x0; maxx = x1; } else { minx = x1; maxx = x0; }
	if (y0 < y1) { miny = y0; maxy = y1; } else { miny = y1; maxy = y0; }
	tcr = cr / sqrt(xxyy); vx = fabs((float)(x1-x0))*tcr; vy = fabs((float)(y1-y0))*tcr;
	t = vx*uz + vy;
	ftol((float)minx-t,&minx); if (minx < 0) minx = 0;
	ftol((float)maxx+t,&maxx); if (maxx >= VSID) maxx = VSID-1;
	t = vy*uz + vx;
	ftol((float)miny-t,&miny); if (miny < 0) miny = 0;
	ftol((float)maxy+t,&maxy); if (maxy >= VSID) maxy = VSID-1;

	vx5.minx = minx; vx5.maxx = maxx+1;
	vx5.miny = miny; vx5.maxy = maxy+1;
	vx5.minz = z0-cr; vx5.maxz = z1+cr+1;
	if (bakit) voxbackup(minx,miny,maxx+1,maxy+1,bakit);

	vx = (fabs(ux) < fabs(uy)); vy = 1.0f-vx; vz = 0;
	ax = uy*vz - uz*vy; ay = uz*vx - ux*vz; az = ux*vy - uy*vx;
	t = 1.0 / (sqrt(ax*ax + ay*ay + az*az)*cr);
	ax *= t; ay *= t; az *= t;
	bx = ay*uz - az*uy; by = az*ux - ax*uz; bz = ax*uy - ay*ux;

	Za = az*az + bz*bz; rZa = 1.0f / Za;
	if (cz != 0) { rcz = 1.0f / cz; vz0i = -rcz*cx; }
	if (y0 != y1)
	{
		t = 1.0f / ((float)(y1-y0)); fxi = ((float)(x1-x0))*t;
		fx = ((float)miny-y0)*fxi + x0; xof = fabs(tcr*xxyy*t);
	}
	else { fx = (float)minx; fxi = 0.0; xof = (float)(maxx-minx); }

	vy = (float)(miny-y0);
	vxo = vy*ay - z0*az;
	vyo = vy*by - z0*bz;
	vzo = vy*cy - z0*cz;
	for(iy=miny;iy<=maxy;iy++)
	{
		ftol(fx-xof,&ix0); if (ix0 < minx) ix0 = minx;
		ftol(fx+xof,&ix1); if (ix1 > maxx) ix1 = maxx;
		fx += fxi;

		vx = (float)(ix0-x0);
		vx0 = vx*ax + vxo; vxo += ay;
		vy0 = vx*bx + vyo; vyo += by;
		vz0 = vx*cx + vzo; vzo += cy;

		if (cz != 0)   //(vx0 + vx1*t)ý + (vy0 + vy1*t)ý = 1
		{
			vz0 *= -rcz;
			for(ix=ix0;ix<=ix1;ix++,vx0+=ax,vy0+=bx,vz0+=vz0i)
			{
				Zb = vx0*az + vy0*bz; Zc = vx0*vx0 + vy0*vy0 - 1;
				t = Zb*Zb - Za*Zc; if (*(long *)&t <= 0) continue; t = sqrt(t);
				ftol(max((-Zb-t)*rZa,vz0    ),&iz0); if (iz0 < 0) iz0 = 0;
				ftol(min((-Zb+t)*rZa,vz0+rcz),&iz1); if (iz1 > MAXZDIM) iz1 = MAXZDIM;
				modslab(scum2(ix,iy),iz0,iz1);
			}
		}
		else
		{
			for(ix=ix0;ix<=ix1;ix++,vx0+=ax,vy0+=bx,vz0+=cx)
			{
				if (*(unsigned long *)&vz0 >= 0x3f800000) continue; //vz0<0||vz0>=1
				Zb = vx0*az + vy0*bz; Zc = vx0*vx0 + vy0*vy0 - 1;
				t = Zb*Zb - Za*Zc; if (*(long *)&t <= 0) continue; t = sqrt(t);
				ftol((-Zb-t)*rZa,&iz0); if (iz0 < 0) iz0 = 0;
				ftol((-Zb+t)*rZa,&iz1); if (iz1 > MAXZDIM) iz1 = MAXZDIM;
				modslab(scum2(ix,iy),iz0,iz1);
			}
		}
	}
	scum2finish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz, dacol );
}

	//Draws a rectangle, given: 2 points as opposite corners, and a color
void setrect (lpoint3d *hit, lpoint3d *hit2, long dacol)
{
	long x, y, xs, ys, zs, xe, ye, ze;

		//WARNING: do NOT use lbound because 'c' not guaranteed to be >= 'b'
	xs = max(min(hit->x,hit2->x),0); xe = min(max(hit->x,hit2->x),VSID-1);
	ys = max(min(hit->y,hit2->y),0); ye = min(max(hit->y,hit2->y),VSID-1);
	zs = max(min(hit->z,hit2->z),0); ze = min(max(hit->z,hit2->z),MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;

	ze++;
	if (dacol == -1)
	{
		for(y=ys;y<=ye;y++)
			for(x=xs;x<=xe;x++)
				delslab(scum2(x,y),zs,ze);
	}
	else
	{
		for(y=ys;y<=ye;y++)
			for(x=xs;x<=xe;x++)
				insslab(scum2(x,y),zs,ze);
	}
	scum2finish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
}

	//Does CSG using pre-sorted spanlist
void setspans (vspans *lst, long lstnum, lpoint3d *offs, long dacol)
{
	void (*modslab)(long *, long, long);
	long i, j, x, y, z0, z1, *lptr;
	char ox, oy;

	if (lstnum <= 0) return;
	if (dacol == -1) modslab = delslab; else modslab = insslab;
	vx5.minx = vx5.maxx = ((long)lst[0].x)+offs->x;
	vx5.miny = ((long)lst[       0].y)+offs->y;
	vx5.maxy = ((long)lst[lstnum-1].y)+offs->y+1;
	vx5.minz = vx5.maxz = ((long)lst[0].z0)+offs->z;

	i = 0; goto in2setlist;
	do
	{
		if ((ox != lst[i].x) || (oy != lst[i].y))
		{
in2setlist:;
			ox = lst[i].x; oy = lst[i].y;
			x = ((long)lst[i].x)+offs->x;
			y = ((long)lst[i].y)+offs->y;
				  if (x < vx5.minx) vx5.minx = x;
			else if (x > vx5.maxx) vx5.maxx = x;
			lptr = scum2(x,y);
		}
		if ((x|y)&(~(VSID-1))) { i++; continue; }
		z0 = ((long)lst[i].z0)+offs->z;   if (z0 < 0) z0 = 0;
		z1 = ((long)lst[i].z1)+offs->z+1; if (z1 > MAXZDIM) z1 = MAXZDIM;
		if (z0 < vx5.minz) vx5.minz = z0;
		if (z1 > vx5.maxz) vx5.maxz = z1;
		modslab(lptr,z0,z1);
		i++;
	} while (i < lstnum);
	vx5.maxx++; vx5.maxz++;
	if (vx5.minx < 0) vx5.minx = 0;
	if (vx5.miny < 0) vx5.miny = 0;
	if (vx5.maxx > VSID) vx5.maxx = VSID;
	if (vx5.maxy > VSID) vx5.maxy = VSID;

	scum2finish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
}

void setheightmap (const unsigned char *hptr, long hpitch, long hxdim, long hydim,
						 long x0, long y0, long x1, long y1)
{
	long x, y, su, sv, u, v;

	if (x0 < 0) x0 = 0;
	if (y0 < 0) y0 = 0;
	if (x1 > VSID) x1 = VSID;
	if (y1 > VSID) y1 = VSID;
	vx5.minx = x0; vx5.maxx = x1;
	vx5.miny = y0; vx5.maxy = y1;
	vx5.minz = 0; vx5.maxz = MAXZDIM;
	if ((x0 >= x1) || (y0 >= y1)) return;

	su = x0%hxdim; sv = y0%hydim;
	for(y=y0,v=sv;y<y1;y++)
	{
		for(x=x0,u=su;x<x1;x++)
		{
			insslab(scum2(x,y),hptr[v*hpitch+u],MAXZDIM);
			u++; if (u >= hxdim) u = 0;
		}
		v++; if (v >= hydim) v = 0;
	}
	scum2finish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,0 );
}

#define MAXPOINTS (256 *2) //Leave the *2 here for safety!
point3d nm[MAXPOINTS*2+2];
float nmc[MAXPOINTS*2+2];
long tri[MAXPOINTS*8+8], lnk[MAXPOINTS*8+8], tricnt;


void initetrasid (point3d *pt, long z)
{
	long i, j, k;
	float x0, y0, z0, x1, y1, z1;

	i = tri[z*4]; j = tri[z*4+1]; k = tri[z*4+2];
	x0 = pt[i].x-pt[k].x; y0 = pt[i].y-pt[k].y; z0 = pt[i].z-pt[k].z;
	x1 = pt[j].x-pt[k].x; y1 = pt[j].y-pt[k].y; z1 = pt[j].z-pt[k].z;
	nm[z].x = y0*z1 - z0*y1;
	nm[z].y = z0*x1 - x0*z1;
	nm[z].z = x0*y1 - y0*x1;
	nmc[z] = nm[z].x*pt[k].x + nm[z].y*pt[k].y + nm[z].z*pt[k].z;
}

void inithull3d (point3d *pt, long nump)
{
	float px, py, pz;
	long i, k, s, z, szz, zz, zx, snzz, nzz, zzz, otricnt;

	tri[0] = 0; tri[4] = 0; tri[8] = 0; tri[12] = 1;
	tri[1] = 1; tri[2] = 2; initetrasid(pt,0);
	if (nm[0].x*pt[3].x + nm[0].y*pt[3].y + nm[0].z*pt[3].z >= nmc[0])
	{
		tri[1] = 1; tri[2] = 2; lnk[0] = 10; lnk[1] = 14; lnk[2] = 4;
		tri[5] = 2; tri[6] = 3; lnk[4] = 2; lnk[5] = 13; lnk[6] = 8;
		tri[9] = 3; tri[10] = 1; lnk[8] = 6; lnk[9] = 12; lnk[10] = 0;
		tri[13] = 3; tri[14] = 2; lnk[12] = 9; lnk[13] = 5; lnk[14] = 1;
	}
	else
	{
		tri[1] = 2; tri[2] = 1; lnk[0] = 6; lnk[1] = 12; lnk[2] = 8;
		tri[5] = 3; tri[6] = 2; lnk[4] = 10; lnk[5] = 13; lnk[6] = 0;
		tri[9] = 1; tri[10] = 3; lnk[8] = 2; lnk[9] = 14; lnk[10] = 4;
		tri[13] = 2; tri[14] = 3; lnk[12] = 1; lnk[13] = 5; lnk[14] = 9;
	}
	tricnt = 4*4;

	for(z=0;z<4;z++) initetrasid(pt,z);

	for(z=4;z<nump;z++)
	{
		px = pt[z].x; py = pt[z].y; pz = pt[z].z;
		for(zz=tricnt-4;zz>=0;zz-=4)
		{
			i = (zz>>2);
			if (nm[i].x*px + nm[i].y*py + nm[i].z*pz >= nmc[i]) continue;

			s = 0;
			for(zx=2;zx>=0;zx--)
			{
				i = (lnk[zz+zx]>>2);
				s += (nm[i].x*px + nm[i].y*py + nm[i].z*pz < nmc[i]) + s;
			}
			if (s == 7) continue;

			nzz = ((0x4a4>>(s+s))&3); szz = zz; otricnt = tricnt;
			do
			{
				snzz = nzz;
				do
				{
					zzz = nzz+1; if (zzz >= 3) zzz = 0;

						//Insert triangle tricnt: (p0,p1,z)
					tri[tricnt+0] = tri[zz+nzz];
					tri[tricnt+1] = tri[zz+zzz];
					tri[tricnt+2] = z;
					initetrasid(pt,tricnt>>2);
					k = lnk[zz+nzz]; lnk[tricnt] = k; lnk[k] = tricnt;
					lnk[tricnt+1] = tricnt+6;
					lnk[tricnt+2] = tricnt-3;
					tricnt += 4;

						//watch out for loop inside single triangle
					if (zzz == snzz) goto endit;
					nzz = zzz;
				} while (!(s&(1<<zzz)));
				do
				{
					i = zz+nzz;
					zz = (lnk[i]&~3);
					nzz = (lnk[i]&3)+1; if (nzz == 3) nzz = 0;
					s = 0;
					for(zx=2;zx>=0;zx--)
					{
						i = (lnk[zz+zx]>>2);
						s += (nm[i].x*px + nm[i].y*py + nm[i].z*pz < nmc[i]) + s;
					}
				} while (s&(1<<nzz));
			} while (zz != szz);
endit:;  lnk[tricnt-3] = otricnt+2; lnk[otricnt+2] = tricnt-3;

			for(zz=otricnt-4;zz>=0;zz-=4)
			{
				i = (zz>>2);
				if (nm[i].x*px + nm[i].y*py + nm[i].z*pz < nmc[i])
				{
					tricnt -= 4; //Delete triangle zz%
					nm[i] = nm[tricnt>>2]; nmc[i] = nmc[tricnt>>2];
					for(i=0;i<3;i++)
					{
						tri[zz+i] = tri[tricnt+i];
						lnk[zz+i] = lnk[tricnt+i];
						lnk[lnk[zz+i]] = zz+i;
					}
				}
			}
			break;
		}
	}
	tricnt >>= 2;
}

static long incmod3[3];
void tmaphulltrisortho (point3d *pt)
{
	point3d *i0, *i1;
	float r, knmx, knmy, knmc, xinc;
	long i, k, op, p, pe, y, yi, z, zi, sy, sy1, itop, ibot, damost;

	for(k=0;k<tricnt;k++)
	{
		if (nm[k].z >= 0)
			{ damost = (long)umost; incmod3[0] = 1; incmod3[1] = 2; incmod3[2] = 0; }
		else
			{ damost = (long)dmost; incmod3[0] = 2; incmod3[1] = 0; incmod3[2] = 1; }

		itop = (pt[tri[(k<<2)+1]].y < pt[tri[k<<2]].y); ibot = 1-itop;
			  if (pt[tri[(k<<2)+2]].y < pt[tri[(k<<2)+itop]].y) itop = 2;
		else if (pt[tri[(k<<2)+2]].y > pt[tri[(k<<2)+ibot]].y) ibot = 2;

			//Pre-calculations
		if (fabs(nm[k].z) < .000001) r = 0; else r = -65536.0 / nm[k].z;
		knmx = nm[k].x*r; knmy = nm[k].y*r;
		//knmc = 65536.0-nmc[k]*r-knmx-knmy;
		//knmc = -nmc[k]*r-(knmx+knmy)*.5f;
		//knmc = /*65536.0-nmc[k]*r+knmx;
		knmc = -nmc[k]*r+knmx;
		ftol(knmx,&zi);

		i = ibot;
		do
		{
			i1 = &pt[tri[(k<<2)+i]]; ftol(i1->y,&sy1); i = incmod3[i];
			i0 = &pt[tri[(k<<2)+i]]; ftol(i0->y,&sy); if (sy == sy1) continue;
			xinc = (i1->x-i0->x)/(i1->y-i0->y);
			ftol((((float)sy-i0->y)*xinc+i0->x)*65536,&y); ftol(xinc*65536,&yi);
			for(;sy<sy1;sy++,y+=yi) lastx[sy] = (y>>16);
		} while (i != itop);
		do
		{
			i0 = &pt[tri[(k<<2)+i]]; ftol(i0->y,&sy); i = incmod3[i];
			i1 = &pt[tri[(k<<2)+i]]; ftol(i1->y,&sy1); if (sy == sy1) continue;
			xinc = (i1->x-i0->x)/(i1->y-i0->y);
			ftol((((float)sy-i0->y)*xinc+i0->x)*65536,&y); ftol(xinc*65536,&yi);
			op = sy*VSID+damost;
			for(;sy<sy1;sy++,y+=yi,op+=VSID)
			{
				ftol(knmx*(float)lastx[sy] + knmy*(float)sy + knmc,&z);
				pe = (y>>16)+op; p = lastx[sy]+op;
				for(;p<pe;p++,z+=zi) *(char *)p = (z>>16);
			}
		} while (i != ibot);
	}
}

void sethull3d (point3d *pt, long nump, long dacol, long bakit)
{
	void (*modslab)(long *, long, long);
	float fminx, fminy, fminz, fmaxx, fmaxy, fmaxz;
	long i, x, y, xs, ys, xe, ye, z0, z1;

	if (nump > (MAXPOINTS>>1)) nump = (MAXPOINTS>>1); //DANGER!!!

	fminx = fminy = VSID; fminz = MAXZDIM; fmaxx = fmaxy = fmaxz = 0;
	for(i=0;i<nump;i++)
	{
		pt[i].x = min(max(pt[i].x,0),VSID-1);
		pt[i].y = min(max(pt[i].y,0),VSID-1);
		pt[i].z = min(max(pt[i].z,0),MAXZDIM-1);

		if (pt[i].x < fminx) fminx = pt[i].x;
		if (pt[i].y < fminy) fminy = pt[i].y;
		if (pt[i].z < fminz) fminz = pt[i].z;
		if (pt[i].x > fmaxx) fmaxx = pt[i].x;
		if (pt[i].y > fmaxy) fmaxy = pt[i].y;
		if (pt[i].z > fmaxz) fmaxz = pt[i].z;
	}

	ftol(fminx,&xs); if (xs < 0) xs = 0;
	ftol(fminy,&ys); if (ys < 0) ys = 0;
	ftol(fmaxx,&xe); if (xe >= VSID) xe = VSID-1;
	ftol(fmaxy,&ye); if (ye >= VSID) ye = VSID-1;
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	ftol(fminz-.5,&vx5.minz); ftol(fmaxz+.5,&vx5.maxz);
	if ((xs > xe) || (ys > ye))
		{ if (bakit) voxbackup(xs,ys,xs,ys,bakit); return; }
	if (bakit) voxbackup(xs,ys,xe,ye,bakit);

	i = ys*VSID+(xs&~3); x = ((((xe+3)&~3)-(xs&~3))>>2)+1;
	for(y=ys;y<=ye;y++,i+=VSID)
		{ clearbuf((void *)&umost[i],x,-1); clearbuf((void *)&dmost[i],x,0); }

	inithull3d(pt,nump);
	tmaphulltrisortho(pt);

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;

	if (dacol == -1) modslab = delslab; else modslab = insslab;

	for(y=ys;y<=ye;y++)
		for(x=xs;x<=xe;x++)
			modslab(scum2(x,y),(long)umost[y*VSID+x],(long)dmost[y*VSID+x]);
	scum2finish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
}

// ------------------------- CONVEX 3D HULL CODE ENDS -------------------------

	//Old&Slow sector code, but only this one supports the 3D bumpmapping :(
static void setsectorb (point3d *p, long *point2, long n, float thick, long dacol, long bakit, long bumpmap)
{
	point3d norm, p2;
	float d, f, x0, y0, x1, y1;
	long i, j, k, got, x, y, z, xs, ys, zs, xe, ye, ze, maxis, ndacol;

	norm.x = 0; norm.y = 0; norm.z = 0;
	for(i=0;i<n;i++)
	{
		j = point2[i]; k = point2[j];
		norm.x += (p[i].y-p[j].y)*(p[k].z-p[j].z) - (p[i].z-p[j].z)*(p[k].y-p[j].y);
		norm.y += (p[i].z-p[j].z)*(p[k].x-p[j].x) - (p[i].x-p[j].x)*(p[k].z-p[j].z);
		norm.z += (p[i].x-p[j].x)*(p[k].y-p[j].y) - (p[i].y-p[j].y)*(p[k].x-p[j].x);
	}
	f = 1.0 / sqrt(norm.x*norm.x + norm.y*norm.y + norm.z*norm.z);
	norm.x *= f; norm.y *= f; norm.z *= f;

	if ((fabs(norm.z) >= fabs(norm.x)) && (fabs(norm.z) >= fabs(norm.y)))
		maxis = 2;
	else if (fabs(norm.y) > fabs(norm.x))
		maxis = 1;
	else
		maxis = 0;

	xs = xe = p[0].x;
	ys = ye = p[0].y;
	zs = ze = p[0].z;
	for(i=n-1;i;i--)
	{
		if (p[i].x < xs) xs = p[i].x;
		if (p[i].y < ys) ys = p[i].y;
		if (p[i].z < zs) zs = p[i].z;
		if (p[i].x > xe) xe = p[i].x;
		if (p[i].y > ye) ye = p[i].y;
		if (p[i].z > ze) ze = p[i].z;
	}
	xs = max(xs-thick-bumpmap,0); xe = min(xe+thick+bumpmap,VSID-1);
	ys = max(ys-thick-bumpmap,0); ye = min(ye+thick+bumpmap,VSID-1);
	zs = max(zs-thick-bumpmap,0); ze = min(ze+thick+bumpmap,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;
	if (bakit) voxbackup(xs,ys,xe+1,ye+1,bakit);

	clearbuf((void *)&templongbuf[zs],ze-zs+1,-3);

	ndacol = (dacol==-1)-2;

	for(y=ys;y<=ye;y++)
		for(x=xs;x<=xe;x++)
		{
			got = 0;
			d = ((float)x-p[0].x)*norm.x + ((float)y-p[0].y)*norm.y + ((float)zs-p[0].z)*norm.z;
			for(z=zs;z<=ze;z++,d+=norm.z)
			{
				if (bumpmap)
				{
					if (d < -thick) continue;
					p2.x = (float)x - d*norm.x;
					p2.y = (float)y - d*norm.y;
					p2.z = (float)z - d*norm.z;
					if (d > (float)hpngcolfunc(&p2)+thick) continue;
				}
				else
				{
					if (fabs(d) > thick) continue;
					p2.x = (float)x - d*norm.x;
					p2.y = (float)y - d*norm.y;
					p2.z = (float)z - d*norm.z;
				}

				k = 0;
				for(i=n-1;i>=0;i--)
				{
					j = point2[i];
					switch(maxis)
					{
						case 0: x0 = p[i].z-p2.z; x1 = p[j].z-p2.z;
								  y0 = p[i].y-p2.y; y1 = p[j].y-p2.y; break;
						case 1: x0 = p[i].x-p2.x; x1 = p[j].x-p2.x;
								  y0 = p[i].z-p2.z; y1 = p[j].z-p2.z; break;
						case 2: x0 = p[i].x-p2.x; x1 = p[j].x-p2.x;
								  y0 = p[i].y-p2.y; y1 = p[j].y-p2.y; break;
						default: _gtfo(); //tells MSVC default can't be reached
					}
					if (((*(long *)&y0)^(*(long *)&y1)) < 0)
					{
						if (((*(long *)&x0)^(*(long *)&x1)) >= 0) k ^= (*(long *)&x0);
						else { f = (x0*y1-x1*y0); k ^= (*(long *)&f)^(*(long *)&y1); }
					}
				}
				if (k >= 0) continue;

				templongbuf[z] = ndacol; got = 1;
			}
			if (got)
			{
				scum(x,y,zs,ze+1,templongbuf);
				clearbuf((void *)&templongbuf[zs],ze-zs+1,-3);
			}
		}
	scumfinish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
}

	//This is for ordfillpolygon&splitpoly
typedef struct { long p, i, t; } raster;
#define MAXCURS 100 //THIS IS VERY EVIL... FIX IT!!!
static raster rst[MAXCURS];
static long slist[MAXCURS];

	//Code taken from POLYOLD\POLYSPLI.BAS:splitpoly (06/09/2001)
void splitpoly (float *px, float *py, long *point2, long *bakn,
					 float x0, float y0, float dx, float dy)
{
	long i, j, s2, n, sn, splcnt, z0, z1, z2, z3;
	float t, t1;

	n = (*bakn); if (n < 3) return;
	i = 0; s2 = sn = n; splcnt = 0;
	do
	{
		t1 = (px[i]-x0)*dy - (py[i]-y0)*dx;
		do
		{
			j = point2[i]; point2[i] |= 0x80000000;
			t = t1; t1 = (px[j]-x0)*dy - (py[j]-y0)*dx;
			if ((*(long *)&t) < 0)
				{ px[n] = px[i]; py[n] = py[i]; point2[n] = n+1; n++; }
			if (((*(long *)&t) ^ (*(long *)&t1)) < 0)
			{
				if ((*(long *)&t) < 0) slist[splcnt++] = n;
				t /= (t-t1);
				px[n] = (px[j]-px[i])*t + px[i];
				py[n] = (py[j]-py[i])*t + py[i];
				point2[n] = n+1; n++;
			}
			i = j;
		} while (point2[i] >= 0);
		if (n > s2) { point2[n-1] = s2; s2 = n; }
		for(i=sn-1;(i) && (point2[i] < 0);i--);
	} while (i > 0);

	if (fabs(dx) > fabs(dy))
	{
		for(i=1;i<splcnt;i++)
		{
			z0 = slist[i];
			for(j=0;j<i;j++)
			{
				z1 = point2[z0]; z2 = slist[j]; z3 = point2[z2];
				if (fabs(px[z0]-px[z3])+fabs(px[z2]-px[z1]) < fabs(px[z0]-px[z1])+fabs(px[z2]-px[z3]))
					{ point2[z0] = z3; point2[z2] = z1; }
			}
		}
	}
	else
	{
		for(i=1;i<splcnt;i++)
		{
			z0 = slist[i];
			for(j=0;j<i;j++)
			{
				z1 = point2[z0]; z2 = slist[j]; z3 = point2[z2];
				if (fabs(py[z0]-py[z3])+fabs(py[z2]-py[z1]) < fabs(py[z0]-py[z1])+fabs(py[z2]-py[z3]))
					{ point2[z0] = z3; point2[z2] = z1; }
			}
		}
	}

	for(i=sn;i<n;i++)
		{ px[i-sn] = px[i]; py[i-sn] = py[i]; point2[i-sn] = point2[i]-sn; }
	(*bakn) = n-sn;
}

void ordfillpolygon (float *px, float *py, long *point2, long n, long day, long xs, long xe, void (*modslab)(long *, long, long))
{
	float f;
	long k, i, z, zz, z0, z1, zx, sx0, sy0, sx1, sy1, sy, nsy, gap, numrst;
	long np, ni;

	if (n < 3) return;

	for(z=0;z<n;z++) slist[z] = z;

		//Sort points by y's
	for(gap=(n>>1);gap;gap>>=1)
		for(z=0;z<n-gap;z++)
			for(zz=z;zz>=0;zz-=gap)
			{
				if (py[point2[slist[zz]]] <= py[point2[slist[zz+gap]]]) break;
				z0 = slist[zz]; slist[zz] = slist[zz+gap]; slist[zz+gap] = z0;
			}

	ftol(py[point2[slist[0]]]+.5,&sy); if (sy < xs) sy = xs;

	numrst = 0; z = 0; n--; //Note: n is local variable!
	while (z < n)
	{
		z1 = slist[z]; z0 = point2[z1];
		for(zx=0;zx<2;zx++)
		{
			ftol(py[z0]+.5,&sy0); ftol(py[z1]+.5,&sy1);
			if (sy1 > sy0) //Insert raster (z0,z1)
			{
				f = (px[z1]-px[z0]) / (py[z1]-py[z0]);
				ftol(((sy-py[z0])*f + px[z0])*65536.0 + 65535.0,&np);
				if (sy1-sy0 >= 2) ftol(f*65536.0,&ni); else ni = 0;
				k = (np<<1)+ni;
				for(i=numrst;i>0;i--)
				{
					if ((rst[i-1].p<<1)+rst[i-1].i < k) break;
					rst[i] = rst[i-1];
				}
				rst[i].i = ni; rst[i].p = np; rst[i].t = (z0<<16)+z1;
				numrst++;
			}
			else if (sy1 < sy0) //Delete raster (z1,z0)
			{
				numrst--;
				k = (z1<<16)+z0; i = 0;
				while ((i < numrst) && (rst[i].t != k)) i++;
				while (i < numrst) { rst[i] = rst[i+1]; i++; }
			}
			z1 = point2[z0];
		}

		z++;
		ftol(py[point2[slist[z]]]+.5,&nsy); if (nsy > xe) nsy = xe;
		for(;sy<nsy;sy++)
			for(i=0;i<numrst;i+=2)
			{
				modslab(scum2(sy,day),max(rst[i].p>>16,0),min(rst[i+1].p>>16,MAXZDIM));
				rst[i].p += rst[i].i; rst[i+1].p += rst[i+1].i;
			}
	}
}

	//Draws a flat polygon
	//given: p&point2: 3D points, n: # points, thick: thickness, dacol: color
static float ppx[MAXCURS*4], ppy[MAXCURS*4];
static long npoint2[MAXCURS*4];
void setsector (point3d *p, long *point2, long n, float thick, long dacol, long bakit)
{
	void (*modslab)(long *, long, long);
	point3d norm;
	float f, rnormy, xth, zth, dax, daz, t, t1;
	long i, j, k, x, y, z, sn, s2, nn, xs, ys, zs, xe, ye, ze;

	norm.x = 0; norm.y = 0; norm.z = 0;
	for(i=0;i<n;i++)
	{
		j = point2[i]; k = point2[j];
		norm.x += (p[i].y-p[j].y)*(p[k].z-p[j].z) - (p[i].z-p[j].z)*(p[k].y-p[j].y);
		norm.y += (p[i].z-p[j].z)*(p[k].x-p[j].x) - (p[i].x-p[j].x)*(p[k].z-p[j].z);
		norm.z += (p[i].x-p[j].x)*(p[k].y-p[j].y) - (p[i].y-p[j].y)*(p[k].x-p[j].x);
	}
	f = 1.0 / sqrt(norm.x*norm.x + norm.y*norm.y + norm.z*norm.z);
	norm.x *= f; norm.y *= f; norm.z *= f;

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;
	else if ((vx5.colfunc == pngcolfunc) && (vx5.pic) && (vx5.xsiz > 0) && (vx5.ysiz > 0) && (vx5.picmode == 3))
	{
			//Find biggest height offset to minimize bounding box size
		j = k = vx5.pic[0];
		for(y=vx5.ysiz-1;y>=0;y--)
		{
			i = y*(vx5.bpl>>2);
			for(x=vx5.xsiz-1;x>=0;x--)
			{
				if (vx5.pic[i+x] < j) j = vx5.pic[i+x];
				if (vx5.pic[i+x] > k) k = vx5.pic[i+x];
			}
		}
		if ((j^k)&0xff000000) //If high bytes are !=, then use bumpmapping
		{
			setsectorb(p,point2,n,thick,dacol,bakit,max(labs(j>>24),labs(k>>24)));
			return;
		}
	}

	xs = xe = p[0].x;
	ys = ye = p[0].y;
	zs = ze = p[0].z;
	for(i=n-1;i;i--)
	{
		if (p[i].x < xs) xs = p[i].x;
		if (p[i].y < ys) ys = p[i].y;
		if (p[i].z < zs) zs = p[i].z;
		if (p[i].x > xe) xe = p[i].x;
		if (p[i].y > ye) ye = p[i].y;
		if (p[i].z > ze) ze = p[i].z;
	}
	xs = max(xs-thick,0); xe = min(xe+thick,VSID-1);
	ys = max(ys-thick,0); ye = min(ye+thick,VSID-1);
	zs = max(zs-thick,0); ze = min(ze+thick,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;
	if (bakit) voxbackup(xs,ys,xe+1,ye+1,bakit);

	if (dacol == -1) modslab = delslab; else modslab = insslab;

	if (fabs(norm.y) >= .001)
	{
		rnormy = 1.0 / norm.y;
		for(y=ys;y<=ye;y++)
		{
			nn = n;
			for(i=0;i<n;i++)
			{
				f = ((float)y-p[i].y) * rnormy;
				ppx[i] = norm.z*f + p[i].z;
				ppy[i] = norm.x*f + p[i].x;
				npoint2[i] = point2[i];
			}
			if (fabs(norm.x) > fabs(norm.z))
			{
				splitpoly(ppx,ppy,npoint2,&nn,p[0].z,((p[0].y-(float)y)*norm.y-thick)/norm.x+p[0].x,norm.x,-norm.z);
				splitpoly(ppx,ppy,npoint2,&nn,p[0].z,((p[0].y-(float)y)*norm.y+thick)/norm.x+p[0].x,-norm.x,norm.z);
			}
			else
			{
				splitpoly(ppx,ppy,npoint2,&nn,((p[0].y-(float)y)*norm.y-thick)/norm.z+p[0].z,p[0].x,norm.x,-norm.z);
				splitpoly(ppx,ppy,npoint2,&nn,((p[0].y-(float)y)*norm.y+thick)/norm.z+p[0].z,p[0].x,-norm.x,norm.z);
			}
			ordfillpolygon(ppx,ppy,npoint2,nn,y,xs,xe,modslab);
		}
	}
	else
	{
		xth = norm.x*thick; zth = norm.z*thick;
		for(y=ys;y<=ye;y++)
		{
			for(z=0;z<n;z++) slist[z] = 0;
			nn = 0; i = 0; sn = n;
			do
			{
				s2 = nn; t1 = p[i].y-(float)y;
				do
				{
					j = point2[i]; slist[i] = 1; t = t1; t1 = p[j].y-(float)y;
					if (((*(long *)&t) ^ (*(long *)&t1)) < 0)
					{
						k = ((*(unsigned long *)&t)>>31); t /= (t-t1);
						daz = (p[j].z-p[i].z)*t + p[i].z;
						dax = (p[j].x-p[i].x)*t + p[i].x;
						ppx[nn+k] = daz+zth; ppx[nn+1-k] = daz-zth;
						ppy[nn+k] = dax+xth; ppy[nn+1-k] = dax-xth;
						npoint2[nn] = nn+1; npoint2[nn+1] = nn+2; nn += 2;
					}
					i = j;
				} while (!slist[i]);
				if (nn > s2) { npoint2[nn-1] = s2; s2 = nn; }
				for(i=sn-1;(i) && (slist[i]);i--);
			} while (i);
			ordfillpolygon(ppx,ppy,npoint2,nn,y,xs,xe,modslab);
		}
	}
	scum2finish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
}

	//Given: p[>=3]: points 0,1 are the axis of rotation, others make up shape
	//      numcurs: number of points
	//        dacol: color
void setlathe (point3d *p, long numcurs, long dacol, long bakit)
{
	point3d norm, ax0, ax1, tp0, tp1;
	float d, f, x0, y0, x1, y1, px, py, pz;
	long i, j, cnt, got, x, y, z, xs, ys, zs, xe, ye, ze, maxis, ndacol;

	norm.x = (p[0].y-p[1].y)*(p[2].z-p[1].z) - (p[0].z-p[1].z)*(p[2].y-p[1].y);
	norm.y = (p[0].z-p[1].z)*(p[2].x-p[1].x) - (p[0].x-p[1].x)*(p[2].z-p[1].z);
	norm.z = (p[0].x-p[1].x)*(p[2].y-p[1].y) - (p[0].y-p[1].y)*(p[2].x-p[1].x);
	f = 1.0 / sqrt(norm.x*norm.x + norm.y*norm.y + norm.z*norm.z);
	norm.x *= f; norm.y *= f; norm.z *= f;

	ax0.x = p[1].x-p[0].x; ax0.y = p[1].y-p[0].y; ax0.z = p[1].z-p[0].z;
	f = 1.0 / sqrt(ax0.x*ax0.x + ax0.y*ax0.y + ax0.z*ax0.z);
	ax0.x *= f; ax0.y *= f; ax0.z *= f;

	ax1.x = ax0.y*norm.z - ax0.z*norm.y;
	ax1.y = ax0.z*norm.x - ax0.x*norm.z;
	ax1.z = ax0.x*norm.y - ax0.y*norm.x;

	x0 = 0; //Cylindrical thickness: Perp-dist from line (p[0],p[1])
	y0 = 0; //Cylindrical min dot product from line (p[0],p[1])
	y1 = 0; //Cylindrical max dot product from line (p[0],p[1])
	for(i=numcurs-1;i;i--)
	{
		d = (p[i].x-p[0].x)*ax0.x + (p[i].y-p[0].y)*ax0.y + (p[i].z-p[0].z)*ax0.z;
		if (d < y0) y0 = d;
		if (d > y1) y1 = d;
		px = (p[i].x-p[0].x) - d*ax0.x;
		py = (p[i].y-p[0].y) - d*ax0.y;
		pz = (p[i].z-p[0].z) - d*ax0.z;
		f = px*px + py*py + pz*pz;     //Note: f is thickness SQUARED
		if (f > x0) x0 = f;
	}
	x0 = sqrt(x0)+1.0;
	tp0.x = ax0.x*y0 + p[0].x; tp1.x = ax0.x*y1 + p[0].x;
	tp0.y = ax0.y*y0 + p[0].y; tp1.y = ax0.y*y1 + p[0].y;
	tp0.z = ax0.z*y0 + p[0].z; tp1.z = ax0.z*y1 + p[0].z;
	xs = max(min(tp0.x,tp1.x)-x0,0); xe = min(max(tp0.x,tp1.x)+x0,VSID-1);
	ys = max(min(tp0.y,tp1.y)-x0,0); ye = min(max(tp0.y,tp1.y)+x0,VSID-1);
	zs = max(min(tp0.z,tp1.z)-x0,0); ze = min(max(tp0.z,tp1.z)+x0,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;
	if (bakit) voxbackup(xs,ys,xe,ye,bakit);

	if ((fabs(norm.z) >= fabs(norm.x)) && (fabs(norm.z) >= fabs(norm.y)))
		maxis = 2;
	else if (fabs(norm.y) > fabs(norm.x))
		maxis = 1;
	else
		maxis = 0;

	clearbuf((void *)&templongbuf[zs],ze-zs+1,-3);

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;

	ndacol = (dacol==-1)-2;

	for(y=ys;y<=ye;y++)
		for(x=xs;x<=xe;x++)
		{
			got = 0;
			d = ((float)x-p[0].x)*ax0.x + ((float)y-p[0].y)*ax0.y + ((float)zs-p[0].z)*ax0.z;
			for(z=zs;z<=ze;z++,d+=ax0.z)
			{
					//Another way: p = sqrt((xyz dot ax1)^2 + (xyz dot norm)^2)
				px = ((float)x-p[0].x) - d*ax0.x;
				py = ((float)y-p[0].y) - d*ax0.y;
				pz = ((float)z-p[0].z) - d*ax0.z;
				f = sqrt(px*px + py*py + pz*pz);

				px = ax0.x*d + ax1.x*f + p[0].x;
				py = ax0.y*d + ax1.y*f + p[0].y;
				pz = ax0.z*d + ax1.z*f + p[0].z;

				cnt = j = 0;
				for(i=numcurs-1;i>=0;i--)
				{
					switch(maxis)
					{
						case 0: x0 = p[i].z-pz; x1 = p[j].z-pz;
								  y0 = p[i].y-py; y1 = p[j].y-py; break;
						case 1: x0 = p[i].x-px; x1 = p[j].x-px;
								  y0 = p[i].z-pz; y1 = p[j].z-pz; break;
						case 2: x0 = p[i].x-px; x1 = p[j].x-px;
								  y0 = p[i].y-py; y1 = p[j].y-py; break;
						default: _gtfo(); //tells MSVC default can't be reached
					}
					if (((*(long *)&y0)^(*(long *)&y1)) < 0)
					{
						if (((*(long *)&x0)^(*(long *)&x1)) >= 0) cnt ^= (*(long *)&x0);
						else { f = (x0*y1-x1*y0); cnt ^= (*(long *)&f)^(*(long *)&y1); }
					}
					j = i;
				}
				if (cnt >= 0) continue;

				templongbuf[z] = ndacol; got = 1;
			}
			if (got)
			{
				scum(x,y,zs,ze+1,templongbuf);
				clearbuf((void *)&templongbuf[zs],ze-zs+1,-3);
			}
		}
	scumfinish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
}

	//Given: p[>=1]: centers
	//   vx5.currad: cutoff value
	//      numcurs: number of points
	//        dacol: color
void setblobs (point3d *p, long numcurs, long dacol, long bakit)
{
	float dx, dy, dz, v, nrad;
	long i, got, x, y, z, xs, ys, zs, xe, ye, ze, ndacol;

	if (numcurs <= 0) return;

		//Boundaries are quick hacks - rewrite this code!!!
	xs = max(p[0].x-64,0); xe = min(p[0].x+64,VSID-1);
	ys = max(p[0].y-64,0); ye = min(p[0].y+64,VSID-1);
	zs = max(p[0].z-64,0); ze = min(p[0].z+64,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;
	if (bakit) voxbackup(xs,ys,xe,ye,bakit);

	clearbuf((void *)&templongbuf[zs],ze-zs+1,-3);

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;

	ndacol = (dacol==-1)-2;

	if (cputype&(1<<25))
	{
		_asm
		{
			mov eax, 256
			xorps xmm7, xmm7    ;xmm7: 0,0,0,0
			cvtsi2ss xmm6, eax  ;xmm6: ?,?,?,256
			movlhps xmm6, xmm6  ;xmm6: ?,256,?,256
		}
	}

	nrad = (float)numcurs / ((float)vx5.currad*(float)vx5.currad + 256.0);
	for(y=ys;y<=ye;y++)
		for(x=xs;x<=xe;x++)
		{
			if (cputype&(1<<25))
			{
				_asm
				{
					cvtsi2ss xmm0, x        ;xmm0:?,?,?,x
					cvtsi2ss xmm7, y        ;xmm7:0,0,0,y
					movlhps xmm0, xmm7      ;xmm0:0,y,?,x
					shufps xmm0, xmm0, 0x08 ;xmm0:x,x,y,x
				}
			}

			got = 0;
			for(z=zs;z<=ze;z++)
			{
				if (cputype&(1<<25))
				{
					_asm
					{
						movhlps xmm3, xmm7       ;"xmm3:?,?,0,0"
						cvtsi2ss xmm7, z         ;"xmm7:0,0,0,z"
						movlhps xmm0, xmm7       ;"xmm0:0,z,y,x"
						mov eax, numcurs
						mov edx, p
						lea eax, [eax+eax*2-3]
				 beg: movups xmm1, [edx+eax*4]   ;"xmm1: ?,pz,py,pz"
						subps xmm1, xmm0         ;"xmm1: ?,dz,dy,dx"
						mulps xmm1, xmm1         ;"xmm1: ?,dzý,dyý,dxý"
						movhlps xmm6, xmm1       ;"xmm6: ?,256,?,dzý"
						shufps xmm1, xmm6, 0x84  ;"xmm1: 256,dzý,dyý,dxý"
						movhlps xmm2, xmm1       ;"xmm2: ?,?,256,dzý"
						addps xmm1, xmm2         ;"xmm1: ?,?,dyý+256,dxý+dzý"
						movss xmm2, xmm1         ;"xmm2: ?,?,256,dxý+dzý"
						shufps xmm1, xmm1, 0x1   ;"xmm1: dxý+dzý,dxý+dzý,dxý+dzý,dyý+256"
						addss xmm1, xmm2         ;"xmm1: ?,?,?,dxý+dyý+dzý+256"
						rcpss xmm1, xmm1         ;"xmm1: ?,?,?,1/(dxý+dyý+dzý+256)"
						addss xmm3, xmm1
						sub eax, 3
						jnc short beg
						movss v, xmm3
					}
				}
				else
				{
					v = 0;
					for(i=numcurs-1;i>=0;i--)
					{
						dx = p[i].x-(float)x;
						dy = p[i].y-(float)y;
						dz = p[i].z-(float)z;
						v += 1.0f / (dx*dx + dy*dy + dz*dz + 256.0f);
					}
				}
				if (*(long *)&v > *(long *)&nrad) { templongbuf[z] = ndacol; got = 1; }
			}
			if (got)
			{
				scum(x,y,zs,ze+1,templongbuf);
				clearbuf((void *)&templongbuf[zs],ze-zs+1,-3);
			}
		}
	scumfinish();
	updatebbox( vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,dacol );
}

//FLOODFILL3D begins --------------------------------------------------------

#define FILLBUFSIZ 16384 //Use realloc instead!
typedef struct { unsigned short x, y, z0, z1; } spoint4d; //128K
static spoint4d fbuf[FILLBUFSIZ];

long dntil0 (long x, long y, long z)
{
	char *v = sptr[y*VSID+x];
	while (1)
	{
		if (z < v[1]) break;
		if (!v[0]) return(MAXZDIM);
		v += v[0]*4;
		if (z < v[3]) return(v[3]);
	}
	return(z);
}

long dntil1 (long x, long y, long z)
{
	char *v = sptr[y*VSID+x];
	while (1)
	{
		if (z <= v[1]) return(v[1]);
		if (!v[0]) break;
		v += v[0]*4;
		if (z < v[3]) break;
	}
	return(z);
}

long uptil1 (long x, long y, long z)
{
	char *v = sptr[y*VSID+x];
	if (z < v[1]) return(0);
	while (v[0])
	{
		v += v[0]*4;
		if (z < v[3]) break;
		if (z < v[1]) return(v[3]);
	}
	return(z);
}

	//Conducts on air and writes solid
void setfloodfill3d (long x, long y, long z, long minx, long miny, long minz,
															long maxx, long maxy, long maxz)
{
	long wholemap, j, z0, z1, nz1, i0, i1, (*bakcolfunc)(lpoint3d *);
	spoint4d a;

	if (minx < 0) minx = 0;
	if (miny < 0) miny = 0;
	if (minz < 0) minz = 0;
	maxx++; maxy++; maxz++;
	if (maxx > VSID) maxx = VSID;
	if (maxy > VSID) maxy = VSID;
	if (maxz > MAXZDIM) maxz = MAXZDIM;
	vx5.minx = minx; vx5.maxx = maxx;
	vx5.miny = miny; vx5.maxy = maxy;
	vx5.minz = minz; vx5.maxz = maxz;
	if ((minx >= maxx) || (miny >= maxy) || (minz >= maxz)) return;

	if ((x < minx) || (x >= maxx) ||
		 (y < miny) || (y >= maxy) ||
		 (z < minz) || (z >= maxz)) return;

	if ((minx != 0) || (miny != 0) || (minz != 0) || (maxx != VSID) || (maxy != VSID) || (maxz != VSID))
		wholemap = 0;
	else wholemap = 1;

	if (isvoxelsolid(x,y,z)) return;

	bakcolfunc = vx5.colfunc; vx5.colfunc = curcolfunc;

	a.x = x; a.z0 = uptil1(x,y,z); if (a.z0 < minz) a.z0 = minz;
	a.y = y; a.z1 = dntil1(x,y,z+1); if (a.z1 > maxz) a.z1 = maxz;
	if (((!a.z0) && (wholemap)) || (a.z0 >= a.z1)) { vx5.colfunc = bakcolfunc; return; } //oops! broke free :/
	insslab(scum2(x,y),a.z0,a.z1); scum2finish();
	i0 = i1 = 0; goto floodfill3dskip;
	do
	{
		a = fbuf[i0]; i0 = ((i0+1)&(FILLBUFSIZ-1));
floodfill3dskip:;
		for(j=3;j>=0;j--)
		{
			if (j&1) { x = a.x+(j&2)-1; if ((x < minx) || (x >= maxx)) continue; y = a.y; }
				 else { y = a.y+(j&2)-1; if ((y < miny) || (y >= maxy)) continue; x = a.x; }

			if (isvoxelsolid(x,y,a.z0)) { z0 = dntil0(x,y,a.z0); z1 = z0; }
										  else { z0 = uptil1(x,y,a.z0); z1 = a.z0; }
			if ((!z0) && (wholemap)) { vx5.colfunc = bakcolfunc; return; } //oops! broke free :/
			while (z1 < a.z1)
			{
				z1 = dntil1(x,y,z1);

				if (z0 < minz) z0 = minz;
				nz1 = z1; if (nz1 > maxz) nz1 = maxz;
				if (z0 < nz1)
				{
					fbuf[i1].x = x; fbuf[i1].y = y;
					fbuf[i1].z0 = z0; fbuf[i1].z1 = nz1;
					i1 = ((i1+1)&(FILLBUFSIZ-1));
					//if (i0 == i1) floodfill stack overflow!
					insslab(scum2(x,y),z0,nz1); scum2finish();
				}
				z0 = dntil0(x,y,z1); z1 = z0;
			}
		}
	} while (i0 != i1);

	vx5.colfunc = bakcolfunc;

	updatebbox(vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,0);
}

void hollowfillstart (long x, long y, long z)
{
	spoint4d a;
	char *v;
	long i, j, z0, z1, i0, i1;

	a.x = x; a.y = y;

	v = sptr[y*VSID+x]; j = ((((long)v)-(long)vbuf)>>2); a.z0 = 0;
	while (1)
	{
		a.z1 = (long)(v[1]);
		if ((a.z0 <= z) && (z < a.z1) && (!(vbit[j>>5]&(1<<j)))) break;
		if (!v[0]) return;
		v += v[0]*4; j += 2;
		a.z0 = (long)(v[3]);
	}
	vbit[j>>5] |= (1<<j); //fill a.x,a.y,a.z0<=?<a.z1

	i0 = i1 = 0; goto floodfill3dskip2;
	do
	{
		a = fbuf[i0]; i0 = ((i0+1)&(FILLBUFSIZ-1));
floodfill3dskip2:;
		for(i=3;i>=0;i--)
		{
			if (i&1) { x = a.x+(i&2)-1; if ((unsigned long)x >= VSID) continue; y = a.y; }
				 else { y = a.y+(i&2)-1; if ((unsigned long)y >= VSID) continue; x = a.x; }

			v = sptr[y*VSID+x]; j = ((((long)v)-(long)vbuf)>>2); z0 = 0;
			while (1)
			{
				z1 = (long)(v[1]);
				if ((z0 < a.z1) && (a.z0 < z1) && (!(vbit[j>>5]&(1<<j))))
				{
					fbuf[i1].x = x; fbuf[i1].y = y;
					fbuf[i1].z0 = z0; fbuf[i1].z1 = z1;
					i1 = ((i1+1)&(FILLBUFSIZ-1));
					if (i0 == i1) return; //floodfill stack overflow!
					vbit[j>>5] |= (1<<j); //fill x,y,z0<=?<z1
				}
				if (!v[0]) break;
				v += v[0]*4; j += 2;
				z0 = (long)(v[3]);
			}
		}
	} while (i0 != i1);
}

	//hollowfill
void sethollowfill ()
{
	long i, j, l, x, y, z0, z1, *lptr, (*bakcolfunc)(lpoint3d *);
	char *v;

	vx5.minx = 0; vx5.maxx = VSID;
	vx5.miny = 0; vx5.maxy = VSID;
	vx5.minz = 0; vx5.maxz = MAXZDIM;

	for(i=0;i<VSID*VSID;i++)
	{
		j = ((((long)sptr[i])-(long)vbuf)>>2);
		for(v=sptr[i];v[0];v+=v[0]*4) { vbit[j>>5] &= ~(1<<j); j += 2; }
		vbit[j>>5] &= ~(1<<j);
	}

	for(y=0;y<VSID;y++)
		for(x=0;x<VSID;x++)
			hollowfillstart(x,y,0);

	bakcolfunc = vx5.colfunc; vx5.colfunc = curcolfunc;
	i = 0;
	for(y=0;y<VSID;y++)
		for(x=0;x<VSID;x++,i++)
		{
			j = ((((long)sptr[i])-(long)vbuf)>>2);
			v = sptr[i]; z0 = MAXZDIM;
			while (1)
			{
				z1 = (long)(v[1]);
				if ((z0 < z1) && (!(vbit[j>>5]&(1<<j))))
				{
					vbit[j>>5] |= (1<<j);
					insslab(scum2(x,y),z0,z1);
				}
				if (!v[0]) break;
				v += v[0]*4; j += 2;
				z0 = (long)(v[3]);
			}
		}
	scum2finish();
	vx5.colfunc = bakcolfunc;
	updatebbox(vx5.minx, vx5.miny, vx5.minz, vx5.maxx, vx5.maxy, vx5.maxz,0);
}

//FLOODFILL3D ends ----------------------------------------------------------