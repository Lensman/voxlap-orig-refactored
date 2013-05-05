/**************************************************************************************************
 * ksnippits.h: Bit's of inline assembly Ken commonly used because the C compiler sucked          *
 **************************************************************************************************/

#pragma once
#include <xmmintrin.h>
#include <stdint.h>

	//Ericson2314's dirty porting tricks
#include "porthacks.h"
#include "voxlap5.h"
#include "kglobals.h"
#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifdef __WATCOMC__

void ftol (float, long *);
#pragma aux ftol =\
	"fistp dword ptr [eax]"\
	parm [8087][eax]\

void dcossin (double, double *, double *);
#pragma aux dcossin =\
	"fsincos"\
	"fstp qword ptr [eax]"\
	"fstp qword ptr [ebx]"\
	parm [8087][eax][ebx]\

void clearbuf (void *, long, long);
#pragma aux clearbuf =\
	"rep stosd"\
	parm [edi][ecx][eax]\
	modify exact [edi ecx]\
	value

long mulshr24 (long, long);
#pragma aux mulshr24 =\
	"imul edx",\
	"shrd eax, edx, 24",\
	parm nomemory [eax][edx]\
	modify exact [eax edx]\
	value [eax]

long umulshr32 (long, long);
#pragma aux umulshr32 =\
	"mul edx"\
	parm nomemory [eax][edx]\
	modify exact [eax edx]\
	value [edx]

long bitrev (long, long);
#pragma aux bitrev =\
	"xor eax, eax"\
	"beg: shr ebx, 1"\
	"adc eax, eax"\
	"dec ecx"\
	"jnz short beg"\
	parm [ebx][ecx]\
	modify nomemory exact [eax ebx ecx]\
	value [eax]

void clearMMX ();
#pragma aux emms =\
	".686"\
	"emms"\
	parm nomemory []\
	modify exact []\
	value

#pragma aux bsf = "bsf eax, eax" parm [eax] modify nomemory exact [eax] value [eax]
#pragma aux bsr = "bsr eax, eax" parm [eax] modify nomemory exact [eax] value [eax]

#else // ! __WATCOMC__

static inline void fcossin (float a, float *c, float *s)
{

	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"fsincos\n"
		: "=t" (*c), "=u" (*s)
		:  "0" (a)
		:
	);
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		fld a
		fsincos
		mov	eax, c
		fstp	dword ptr [eax]
		mov	eax, s
		fstp	dword ptr [eax]
	}
	#else // C Default
	*c = cosf(a);
	*s = cosf(a);
	#endif
}

static inline void dcossin (double a, double *c, double *s)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"fsincos\n"
		: "=t" (*c), "=u" (*s)
		:  "0" (a)
		:
	);
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		fld	a
		fsincos
		mov	eax, c
		fstp	qword ptr [eax]
		mov	eax, s
		fstp	qword ptr [eax]
	}
	#else // C Default
	*c = cos(a);
	*s = sin(a);
	#endif
}

static inline void ftol (float f, long *a)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"fistpl	(%[a])"
		:
		: "t" (f), [a] "r" (a)
		:
	);
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		fld	f
		fistp	dword ptr [eax]
	}
	#else // C Default
	*a = (long) f;
	#endif
}

static inline void dtol (double d, long *a)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"fistpl	(%[a])"
		:
		: "t" (d), [a] "r" (a)
		:
	);
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		fld	qword ptr d
		fistp	dword ptr [eax]
	}
	#else // C Default
	*a = (long) d;
	#endif
}


static inline double dbound (double d, double dmin, double dmax)
{
    #if !defined(NOASM)
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	//WARNING: This ASM code requires >= PPRO
	__asm__ __volatile__
	(
		"fucomi	%[dmin], %[d]\n"      //if (d < dmin)
		"fcmovb	%[dmin], %[d]\n"      //    d = dmin;
		"fucomi	%[dmax], %[d]\n"      //if (d > dmax)
		"fcmovnb	%[dmin], %[d]\n"  //    d = dmax;
		"fucom	%[dmax]\n"
		: [d] "=t" (d)
		:     "0"  (d), [dmin] "f" (dmin), [dmax] "f" (dmax)
		:
	);
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		fld	dmin
		fld	d
		fucomi	st, st(1)   //if (d < dmin)
		fcmovb	st, st(1)   //    d = dmin;
		fld	dmax
		fxch	st(1)
		fucomi	st, st(1)   //if (d > dmax)
		fcmovnb	st, st(1)   //    d = dmax;
		fstp	d
		fucompp
	}
	#endif
	return(d);
    #else // C Default
    return BOUND(d, dmin, dmax);
	#endif
}

static inline long mulshr16 (long a, long d)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"imul	%[d]\n"
		"shrd	$16, %%edx, %[a]\n"
		: [a] "+a" (a)
		: [d]  "r" (d)
		: "edx"
	);
	return a;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		mov	edx, d
		imul	edx
		shrd	eax, edx, 16
	}
	#else // C Default
	return (long)(((int64_t)a * (int64_t)d) >> 16);
	#endif
}

static inline long mulshr24 (long a, long d)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"imul	%[d]\n"
		"shrd	$24, %%edx, %[a]\n"
		: [a] "+a" (a)
		: [d]  "r" (d)
		: "edx"
	);
	return a;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov eax, a
		mov edx, d
		imul edx
		shrd eax, edx, 24
	}
	#else
	return (long)(((int64_t)a * (int64_t)d) >> 24);
	#endif
}

static inline long mulshr32 (long a, long d)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"imul %%edx"
		: "+d" (d)
		:  "a" (a)
	);
	return d;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov eax, a
		imul d
		mov eax, edx
	}
	#else
	return (long)(((int64_t)a * (int64_t)d) >> 32);
	#endif
}

static inline int64_t mul64 (long a, long d)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	int64_t out64;
	__asm__ __volatile__
	(
		"imul	%[d] \n"
		: "=A" (out64)
		:  "a" (a),    [d] "r" (d)
		:
	);
	return out64;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		imul	d
	}
	#else // C Default
	return (int64_t)a * (int64_t)d;
	#endif
}

static inline long shldiv16 (long a, long b)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"mov	%[a], %%edx\n"
		"shl	$16, %[a]\n"
		"sar	$16, %%edx\n"
		"idiv	%[b]\n"
		: [a] "+a" (a)
		: [b] "r" (b)
		: "edx"
	);
	return a;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		mov	edx, eax
		shl	eax, 16
		sar	edx, 16
		idiv	b
	}
	#else // C Default
	return (long)(((int64_t)a << 16) / (int64_t)b);
	#endif
}

static inline long isshldiv16safe (long a, long b)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		".intel_syntax prefix\n"
		"test	%[a], %[a]\n"
		"js	short .Lskipneg0\n"
		"neg	%[a]\n"
	".Lskipneg0:\n"
		"sar	%[a], 14\n"

		"test	%[b], %[b]\n"
		"js	short .Lskipneg1\n"
		"neg	%[b]\n"
	".Lskipneg1:\n"
			//abs((a<<16)/b) < (1<<30) //1 extra for good luck!
			//-abs(a)>>14 > -abs(b)    //use -abs because safe for 0x80000000
			//eax-edx < 0
		"sub	%[b], %[a]\n"
		"shr	%[b], 31\n"
		".att_syntax prefix\n"
		:               [b] "=r" (b)
		: [a]  "r" (a),      "0" (b)
		: "cc"
	);
	return b;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	edx, a
		test	edx, edx
		js	short skipneg0
		neg	edx
	skipneg0:
		sar	edx, 14

		mov	eax, b
		test	eax, eax
		js	short skipneg1
		neg	eax
	skipneg1:
			//abs((a<<16)/b) < (1<<30) //1 extra for good luck!
			//-abs(a)>>14 > -abs(b)    //use -abs because safe for 0x80000000
			//eax-edx	< 0
		sub	eax, edx
		shr	eax, 31
	}
    #else // C Default
	return ((uint32_t)((-abs(b) - ((-abs(a)) >> 14)))) >> 31;
	#endif
}

static inline long umulshr32 (long a, long d)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"mul	%%edx\n" //dword ptr
		: "+d" (d)
		:  "a" (a)
		:
	);
	return d;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		mul	d
		mov	eax, edx
	}
	#else // C Default
	return (long)(((uint64_t)a * (uint64_t)d) >> 32);
	#endif
}

static inline long scale (long a, long d, long c)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"imul	%[d]\n"
		"idiv	%[c]\n"
		: "+a" (a)
		: [c] "r" (c), [d] "r" (d)
		:
	);
	return a;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		imul	d
		idiv	c
	}
	#else // C Default
	return (long)((int64_t)a * (int64_t)d / (int64_t)c);
	#endif
}

static inline long dmulshr0 (long a, long d, long s, long t)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"imul	%[d]\n"
		:    "+a" (a)
		: [d] "r" (d)
		:
	);
	__asm__ __volatile__
	(
		"imul	%[d]\n"
		:    "+a" (s)
		: [d] "r" (t)
		:
	);
	return a + s;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		imul	d
		mov	ecx, eax
		mov	eax, s
		imul	t
		add	eax, ecx
	}
	#else // C default
	return (long)((int64_t)a*(int64_t)d + (int64_t)s*(int64_t)t);
	#endif
}

static inline long dmulshr22 (long a, long b, long c, long d)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"imul	%[b]\n"
		:    "+a,a" (a), "=d,d" (b)
		: [b] "r,1" (b)
		:
	);
	__asm__ __volatile__
	(
		"imul	%[d]\n"
		:    "+a,a" (c), "=d,d" (d)
		: [d] "r,1" (d)
		:
	);
	__asm__ __volatile__
	(
		"add	%[a], %[c]\n"
		"adc	%[b], %[d]\n"
		"shrd	$22,  %[d], %[c]\n"
		: [c] "+r" (c)
		: [a]  "r" (a), [b] "r" (b), [d] "r" (d)
		:
	);
	return c;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		imul	b
		mov	ecx, eax
		push	edx
		mov	eax, c
		imul	d
		add	eax, ecx
		pop	ecx
		adc	edx, ecx
		shrd	eax, edx, 22
	}
    #else // C Default
	return (long)(((((int64_t)a)*((int64_t)b)) + (((int64_t)c)*((int64_t)d))) >> 22);
	#endif

}

static inline long dmulrethigh (long b, long c, long a, long d)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"imul	%[d]\n"
		:    "+a,a" (a), "=d,d" (d)
		: [d] "r,1" (d)
		:
	);
	__asm__ __volatile__
	(
		"imul	%[c]\n"
		:    "+a,a" (b), "=d,d" (c)
		: [c] "r,1" (c)
		:
	);
	__asm__ __volatile__
	(
		"sub	%[a], %[b]\n"
		"sbb	%[d], %[c]\n"
		: [c] "+r" (c)
		: [a]  "r" (a), [b] "r" (b), [d] "r" (d)
		:
	);
	return c;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		imul	d
		mov	ecx, eax
		push	edx
		mov	eax, b
		imul	c
		sub	eax, ecx
		pop	ecx
		sbb	edx, ecx
		mov	eax, edx
	}
	#else // C Default
	return (long)(((int64_t)b*(int64_t)c - (int64_t)a*(int64_t)d) >> 32);
	#endif
}

static inline void copybuf (void *s, void *d, long c)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"rep	movsl\n"
		:
		: "S" (s), "D" (d), "c" (c)
		:
	);
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		push	esi
		push	edi
		mov	esi, s
		mov	edi, d
		mov	ecx, c
		rep	movsd
		pop	edi
		pop	esi
	}
	#else // C Default
	int i;
	for (i = 0;	i < c; i++)	((long *)d)[i] = ((long *)s)[i];
	#endif
}

static inline void clearbuf (void *d, long c, long a)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"rep	stosl\n"
		:
		: "D" (d), "c" (c), "a" (a)
		:
	);
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		push	edi
		mov	edi, d
		mov	ecx, c
		mov	eax, a
		rep	stosd
		pop	edi
	}
	#else // C Default
	int i;
	for (i = 0;	i < c; i++)
		((long *)d)[i] = a;
	#endif
}

static inline unsigned long bswap (unsigned long a)
{
	#if defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"bswap	%[a]\n"
		: [a] "+r" (a)
		:
		:
	);
	return a;
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		bswap	eax
	}
	#else // C Default

	#if defined(__GNUC__)
	return __builtin_bswap32(a);
	#elif defined(_MSC_VER)
	return _byteswap_ulong(a);
    #endif

	#endif
}

static inline long bitrev (long b, long c)
{
    #if defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm
	{
		mov edx, b
		mov ecx, c
		xor eax, eax
    beg: shr edx, 1
		adc eax, eax
		sub ecx, 1
		jnz short beg
	}
    #elif defined(__GNUC__) && defined(__i386__) && !defined(NOASM)
	long a;
	__asm__ __volatile__ (
		"xorl %%eax, %%eax\n\t0:\n\tshrl $1, %%ebx\n\tadcl %%eax, %%eax\n\tsubl $1, %%ecx\n\tjnz 0b"
		: "+a" (a), "+b" (b), "+c" (c) : : "cc");
	return a;
    #else
	long i, j;
	for(i=1,j=0,c=(1<<c);i<c;i+=i) { j += j; if (b&i) j++; }
	return(j);
    #endif
}


#pragma warning( disable : 4799 )

/** Add two colors, both represented by 4 bytes, with saturation
 *  Modifies *color
 *  @param color Pointer to color to modify
*/
static inline void mmxcoloradd (long *color)
{
    
    #if defined( __GNUC__ ) && defined(__i386__) && !defined(NOASM)
    #elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	__asm
	{
		mov eax, color
		movd mm0, [eax]
		paddusb mm0, flashbrival
		movd [eax], mm0
	}
    #elif defined(_MSC_VER) && defined(__i386__) && defined(NOASM)
	*color = _m_to_int(_mm_adds_pu8(_m_from_int(*color), _m_from_int(flashbrival)));
	#endif

}
/** Subtract two colors, both represented by 4 bytes, with saturation
 *  Modifies *color
 *  @param color Pointer to color to modify
*/
static inline void mmxcolorsub (long *color)
{
    #if defined( __GNUC__ ) && defined(__i386__) && !defined(NOASM)
    #elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	__asm
	{
		mov eax, color
		movd mm0, [eax]
		psubusb mm0, flashbrival
		movd [eax], mm0
	}
    #elif defined(_MSC_VER) && defined(__i386__) && defined(NOASM)
	*color = _m_to_int(_mm_subs_pu8(_m_from_int(*color), _m_from_int(flashbrival)));
	#endif
}

static inline void clearMMX () 
{ 
	#if defined( __GNUC__ ) && defined(__i386__) && !defined(NOASM)
    #elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_mm_empty(); 
	#endif
}

/** Bitscan forward */
static inline long bsf (long a)
{
	#if defined( __GNUC__ ) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"bsf %0, %%eax\n"
		:
		: "m" (a)
		:
	);
	#elif defined(__GNUC__) && defined(__LP64__) && !defined(NOASM)
	assert (a != 0);
	__asm__ __volatile__ 
	(
		"bsfq %0, %0" 
		: "=r" (a) 
		: "0" (a)
	);
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm bsf eax, a
	#elif defined(_MSC_VER) && defined(__i386__) && defined(NOASM)
	#pragma intrinsic(_BitScanForward)
	unsigned long index;
	_BitScanForward(&index, a);
	return index;
	#elif defined(_MSC_VER) && defined(__LP64__) && defined(NOASM)
	#pragma intrinsic(_BitScanForward64)
	unsigned long index;
	_BitScanForward64(&index, (__int64)a);
	return index;
    #endif
}

/** Bitscan reverse */
static inline long bsr (long a)
{
	#if defined( __GNUC__ ) && defined(__i386__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"bsr %0, %%eax\n"
		:
		: "m" (a)
		:
	);
	#elif defined(__GNUC__) && defined(__LP64__) && !defined(NOASM)
	assert (a != 0);
	__asm__ __volatile__ 
	(
		"bsrq %0, %0" 
		: "=r" (a) 
		: "0" (a)
	);
	#elif defined(_MSC_VER) && defined(__i386__) && !defined(NOASM)
	_asm bsr eax, a
	#elif defined(_MSC_VER) && defined(__i386__) && defined(NOASM)
	#pragma intrinsic(_BitScanReverse)
	unsigned long index;
	_BitScanReverse(&index, a);
	return index;
	#elif defined(_MSC_VER) && defined(__LP64__) && defined(NOASM)
	#pragma intrinsic(_BitScanReverse64)
	unsigned long index;
	_BitScanReverse64(&index, (__int64)a);
	return index;
    #endif
}
/** @shared */
static inline void movps (point4d *dest, point4d *src)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		""
		: "=x" (dest->vec)
		:  "0" (src->vec)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, src
		movaps	xmm7, [eax]
		mov	eax, dest
		movaps	[eax], xmm7
	}
	#else // C Default
	*dest = *src;
	#endif
}
/** @shared */
static inline void intss (point4d *dest, long src)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"cvtsi2ss	%[src], %[dest]\n"
		"shufps	$0, %[dest], %[dest]\n"
		: [dest] "=x" (dest->vec)
		: [src]   "g" (src)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, dest
		cvtsi2ss	xmm7, src
		shufps	xmm7, xmm7, 0
		movaps	[eax], xmm7
	}
	#else // C Default
	dest->x = dest->y = dest->z = dest->z2 = (float)src;
	#endif
}
/** @shared */
static inline void addps (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"addps	%[b], %[a]\n"
		: [a] "=x" (sum->vec)
		:      "0" (a->vec), [b] "x" (b->vec)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movaps	xmm7, [eax]
		mov	eax, b
		addps	xmm7, [eax]
		mov	eax, sum
		movaps	[eax], xmm7
	}
	#else // C Default
	sum->x  =  a->x  +  b->x;
	sum->y  =  a->y  +  b->y;
	sum->z  =  a->z  +  b->z;
	sum->z2 =  a->z2 +  b->z2;
	#endif
}
/** @shared */
static inline void mulps (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"mulps	%[b], %[a]\n"
		: [a] "=x" (sum->vec)
		:      "0" (a->vec), [b] "x" (b->vec)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movaps	xmm7, [eax]
		mov	eax, b
		mulps	xmm7, [eax]
		mov	eax, sum
		movaps	[eax], xmm7
	}
	#else // C Default
	sum->x  =  a->x  *  b->x;
	sum->y  =  a->y  *  b->y;
	sum->z  =  a->z  *  b->z;
	sum->z2 =  a->z2 *  b->z2;
	#endif
}
/** @shared */
static inline void subps (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"subps	%[b], %[a]\n"
		: [a] "=x" (sum->vec)
		:      "0" (a->vec), [b] "x" (b->vec)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movaps	xmm7, [eax]
		mov	eax, b
		subps	xmm7, [eax]
		mov	eax, sum
		movaps	[eax], xmm7
	}
	#else // C Default
	sum->x  =  a->x  -  b->x;
	sum->y  =  a->y  -  b->y;
	sum->z  =  a->z  -  b->z;
	sum->z2 =  a->z2 -  b->z2;
	#endif

}
/** @shared */
static inline void minps (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"minps	%[b], %[a]\n"
		: [a] "=x" (sum->vec)
		:      "0" (a->vec), [b] "x" (b->vec)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movaps	xmm7, [eax]
		mov	eax, b
		minps	xmm7, [eax]
		mov	eax, sum
		movaps	[eax], xmm7
	}
	#else // C Default
	sum->x  =  MIN(a->x,  b->x);
	sum->y  =  MIN(a->y,  b->y);
	sum->z  =  MIN(a->z,  b->z);
	sum->z2 =  MIN(a->z2, b->z2);
	#endif
}
/** @shared */
static inline void maxps (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
		__asm__ __volatile__
	(
		"maxps	%[b], %[a]\n"
		: [a] "=x" (sum->vec)
		:      "0" (a->vec), [b] "x" (b->vec)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movaps	xmm7, [eax]
		mov	eax, b
		maxps	xmm7, [eax]
		mov	eax, sum
		movaps	[eax], xmm7
	}
    #else // C Default
    sum->x  =  MAX(a->x,  b->x);
	sum->y  =  MAX(a->y,  b->y);
	sum->z  =  MAX(a->z,  b->z);
	sum->z2 =  MAX(a->z2, b->z2);
	#endif
}
/** @shared */
static inline void movps_3dn (point4d *dest, point4d *src)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		""
		: "=y" (dest->svec[0]), "=y" (dest->svec[1])
		:  "0" (src ->svec[0]),  "1" (src ->svec[1])
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, src
		movq	mm0, [eax]
		movq	mm1, [eax+8]
		mov	eax, dest
		movq	[eax], mm0
		movq	[eax+8], mm1
	}
	#else // C Default
	*dest = *src;
	#endif
}
/** @shared */
static inline void intss_3dn (point4d *dest, long src)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"pi2fd	%[y1], %[y1]\n"
		"punpckldq	%[y1], %[y1]\n"
		"movq	%[y1], (%[adrs])\n"
		: [y1] "=y" (dest->svec[0])
		:       "0" (src), [adrs] "r" (&(dest->svec[1]))
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, dest
		movd	mm0, src
		pi2fd	mm0, mm0
		punpckldq	mm0, mm0
		movq	[eax], mm0
		movq	[eax+8], mm0
	}
	#else // C default
    dest->x = dest->y = dest->z = dest->z2 = (float)src;
	#endif
}
/** @shared */
static inline void addps_3dn (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"pfadd	(%[b]),  %[a1]\n"
		"pfadd	8(%[b]), %[a2]\n"
		: [a1] "=y" (sum->svec[0]), [a2] "=y" (sum->svec[1])
		:       "0" (a  ->svec[0]),       "1" (a  ->svec[1]), [b] "r" (&b)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movq	mm0, [eax]
		movq	mm1, [eax+8]
		mov	eax, b
		pfadd	mm0, [eax]
		pfadd	mm1, [eax+8]
		mov	eax, sum
		movq	[eax], mm0
		movq	[eax+8], mm1
	}
	#else // C Default
    sum->x  =  a->x  +  b->x;
	sum->y  =  a->y  +  b->y;
	sum->z  =  a->z  +  b->z;
	sum->z2 =  a->z2 +  b->z2;
	#endif
}
/** @shared */
static inline void mulps_3dn (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"pfmul	(%[b]),  %[a1]\n"
		"pfmul	8(%[b]), %[a2]\n"
		: [a1] "=y" (sum->svec[0]), [a2] "=y" (sum->svec[1])
		:       "0" (a  ->svec[0]),       "1" (a  ->svec[1]), [b] "r" (&b)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movq	mm0, [eax]
		movq	mm1, [eax+8]
		mov	eax, b
		pfmul	mm0, [eax]
		pfmul	mm1, [eax+8]
		mov	eax, sum
		movq	[eax], mm0
		movq	[eax+8], mm1
	}
	#else // C Default
    sum->x  =  a->x  *  b->x;
	sum->y  =  a->y  *  b->y;
	sum->z  =  a->z  *  b->z;
	sum->z2 =  a->z2 *  b->z2;
	#endif
}
/** @shared */
static inline void subps_3dn (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"pfsub	(%[b]),  %[a1]\n"
		"pfsub	8(%[b]), %[a2]\n"
		: [a1] "=y" (sum->svec[0]), [a2] "=y" (sum->svec[1])
		:       "0" (a  ->svec[0]),       "1" (a  ->svec[1]), [b] "r" (&b)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movq	mm0, [eax]
		movq	mm1, [eax+8]
		mov	eax, b
		pfsub	mm0, [eax]
		pfsub	mm1, [eax+8]
		mov	eax, sum
		movq	[eax], mm0
		movq	[eax+8], mm1
	}
	#else // C Default
    sum->x  =  a->x  -  b->x;
	sum->y  =  a->y  -  b->y;
	sum->z  =  a->z  -  b->z;
	sum->z2 =  a->z2 -  b->z2;
	#endif
}
/** @shared */
static inline void minps_3dn (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"pfmin	(%[b]),  %[a1]\n"
		"pfmin	8(%[b]), %[a2]\n"
		: [a1] "=y" (sum->svec[0]), [a2] "=y" (sum->svec[1])
		:       "0" (a  ->svec[0]),       "1" (a  ->svec[1]), [b] "r" (&b)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movq	mm0, [eax]
		movq	mm1, [eax+8]
		mov	eax, b
		pfmin	mm0, [eax]
		pfmin	mm1, [eax+8]
		mov	eax, sum
		movq	[eax], mm0
		movq	[eax+8], mm1
	}
    #else // C Default
    sum->x  =  MIN(a->x,  b->x);
	sum->y  =  MIN(a->y,  b->y);
	sum->z  =  MIN(a->z,  b->z);
	sum->z2 =  MIN(a->z2, b->z2);
	#endif
}
/** @shared */
static inline void maxps_3dn (point4d *sum, point4d *a, point4d *b)
{
    #if defined(__GNUC__) && !defined(NOASM)
	__asm__ __volatile__
	(
		"pfmax	(%[b]),  %[a1]\n"
		"pfmax	8(%[b]), %[a2]\n"
		: [a1] "=y" (sum->svec[0]), [a2] "=y" (sum->svec[1])
		:       "0" (a  ->svec[0]),       "1" (a  ->svec[1]), [b] "r" (&b)
		:
	);
    #elif defined(_MSC_VER) && !defined(NOASM)
	_asm
	{
		mov	eax, a
		movq	mm0, [eax]
		movq	mm1, [eax+8]
		mov	eax, b
		pfmax	mm0, [eax]
		pfmax	mm1, [eax+8]
		mov	eax, sum
		movq	[eax], mm0
		movq	[eax+8], mm1
	}
	#else // C Default
    sum->x  =  MAX(a->x,  b->x);
	sum->y  =  MAX(a->y,  b->y);
	sum->z  =  MAX(a->z,  b->z);
	sum->z2 =  MAX(a->z2, b->z2);
	#endif
}
static inline long fstcw ()
{
	long fpumode;
	#ifdef __GNUC__ //gcc inline asm
	__asm__ __volatile__
	(
		"fstcw %0\n"
		:
		: "m" (fpumode)
		:
	);
	#endif
	#ifdef _MSC_VER //msvc inline asm
	_asm fstcw fpumode
	#endif
	return(fpumode);
}

static inline void fldcw (long fpumode)
{
	#ifdef __GNUC__ //gcc inline asm
	__asm__ __volatile__
	(
		"fldcw %0\n"
		:
		: "m" (fpumode)
		:
	);
	#endif
	#ifdef _MSC_VER //msvc inline asm
	_asm fldcw fpumode
	#endif
}
#endif // ! __WATCOMC__
