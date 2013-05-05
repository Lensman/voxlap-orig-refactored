#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <emmintrin.h>

typedef float v4sf __attribute__ ((aligned(16))) __attribute__ ((vector_size (16)));// vector of four single floats
typedef int v4si __attribute__ ((aligned(16))) __attribute__ ((vector_size (16)));// vector of four signed int
typedef int v8ss __attribute__ ((aligned(16))) __attribute__ ((vector_size (16)));// vector of 8 signed shorts
typedef int v4ss __attribute__ ((aligned(16))) __attribute__ ((vector_size (8)));// vector of 4 signed shorts


union v4si_t {
	v4si v_vec; 
	__m128 v_xmm;
	struct { signed int x,y,z,w; };
	void dump(){
		printf( "vec4ss x=%i y=%i z=%i w=%i\n",x,y,z,w );
	}
};

union v4ss_t {
	v4ss v_vec; 
	__m64 v_xmm;
	struct { short x,y,z,w;	};
	void dump(){
		printf( "vec4ss x=%i y=%i z=%i w=%i\n",x,y,z,w );
	}
};

union v8ss_t 
{
	struct { v4ss_t v1, v2; };
	v8ss v_vec;
	void dump(){
		printf( "vec8ss v1.x=%i v1.y=%i v1.z=%i v1.w=%i v2.x=%i v2.y=%i v2.z=%i v2.w=%i\n" ,v1.x,v1.y,v1.z,v1.w,v2.x,v2.y,v2.z,v2.w );
	}
};

union ss_bbox_t {
	struct { short x1,y1,z1,pad,x2,y2,z2,csgdel; };
	struct { v4ss_t bmin, bmax; };
	v8ss_t v8_vec;
	__m128i v_xmm;
};

union v4f_t 
{
  	struct{float xf,yf,zf,wf;};
  	v4sf vf;
  	__m128 v_xmm;
  	__m128i iv_xmm;
  	struct{int xi,yi,zi,wi;};
};

// "Insert" a 0 bit after each of the 16 low bits of x
uint32_t Part1By1(uint32_t x)
{
  x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x <<  8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x <<  4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x <<  2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x <<  1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  return x;
}

// "Insert" two 0 bits after each of the 10 low bits of x
uint32_t Part1By2(uint32_t x)
{
  x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
  x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
  x = (x ^ (x <<  8)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
  x = (x ^ (x <<  4)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
  x = (x ^ (x <<  2)) & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  return x;
}
// Inverse of Part1By1 - "delete" all odd-indexed bits
uint32_t Compact1By1(uint32_t x)
{
  x &= 0x55555555;                  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  x = (x ^ (x >>  1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x >>  2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x >>  4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x >>  8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
  return x;
}

// Inverse of Part1By2 - "delete" all bits not at positions divisible by 3
uint32_t Compact1By2(uint32_t x)
{
  x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  x = (x ^ (x >>  2)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
  x = (x ^ (x >>  4)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
  x = (x ^ (x >>  8)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
  x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
  return x;
}
uint32_t EncodeMorton2(uint32_t x, uint32_t y)
{
  return (Part1By1(y) << 1) + Part1By1(x);
}

uint32_t EncodeMorton3(uint32_t x, uint32_t y, uint32_t z)
{
  return (Part1By2(z) << 2) + (Part1By2(y) << 1) + Part1By2(x);
}


uint32_t DecodeMorton3X(uint32_t code)
{
  return Compact1By2(code >> 0);
}

uint32_t DecodeMorton3Y(uint32_t code)
{
  return Compact1By2(code >> 1);
}

uint32_t DecodeMorton3Z(uint32_t code)
{
  return Compact1By2(code >> 2);
}
float f_rsqrt( float number )
{
        long i;
        float x2, y;
        const float threehalfs = 1.5F;
 
        x2 = number * 0.5F;
        y  = number;
        i  = * ( long * ) &y;                       // evil floating point bit level hacking
        i  = 0x5f3759df - ( i >> 1 );               // what the fuck?
        y  = * ( float * ) &i;
        y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
      	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
        return y;
}

inline __m128 _mm_hadd4_ps(__m128 i)
{   
	__m128   t;
	t = _mm_movehl_ps(t, i);
	i = _mm_add_ps(i, t);
	t = _mm_shuffle_ps(i, i, 0x55);
	i = _mm_add_ps(i, t);
	return i;
} 

// SSE scalar reciprocal sqrt using rsqrt op, plus one Newton-Rhaphson iteration
inline __m128 _mm_rsqrt( const __m128 x )
{
	__m128 recip = _mm_rsqrt_ss( x );  // "estimate" opcode
	const static __m128 three = { 3, 3, 3, 3 }; // aligned consts for fast load
	const static __m128 half = { 0.5,0.5,0.5,0.5 };
	__m128 halfrecip = _mm_mul_ss( half, recip );
	__m128 threeminus_xrr = _mm_sub_ss( three, _mm_mul_ss( x, _mm_mul_ss ( recip, recip ) ) );
	return _mm_mul_ss( halfrecip, threeminus_xrr );
}
float optistrx[4] = {0.0f, 1.0f, 0.0f, 1.0f};



int main()
{

	float scalar = 100.0f;
	float x = 1670.0f;
	float y = 19.0f;
	float d = 300.0f;
	__attribute__ ((__aligned__(16)) int rgba1[4] = { 0x0001, 0x0002, 0, 0 };
	__attribute__ ((__aligned__(16)) int rgba2[4] = { 0x0003, 0x0004, 0, 0 };
	__attribute__ ((__aligned__(16)) int rgba3[4] = { 0x0005, 0x0006, 0, 0 };
	__attribute__ ((__aligned__(16)) int rgba4[4] = { 0x0007, 0x0008, 0, 0 };

	v4f_t col1, col2, col3, col4;

    col3.iv_xmm = _mm_unpacklo_epi32 (*(__m128i*)rgba1, *(__m128i*)rgba2 );
	col4.iv_xmm = _mm_unpacklo_epi32 (*(__m128i*)rgba3, *(__m128i*)rgba4 );
	col4.iv_xmm = _mm_unpacklo_epi64 ( col4.iv_xmm, col3.iv_xmm );


    //col3.iv_xmm = _mm_shuffle_epi32( col3.iv_xmm , _MM_SHUFFLE (1,0,0,0) );

    
	printf("col1 %d %d %d %d %d\n", col4.xi , col4.yi, col4.zi, col4.wi, col1 );

	union v4f_t vect  = (v4f_t){{ x, y, 0.0f, 0.0f }};
	union v4f_t vect2 = (v4f_t){{ x, y, 0.0f, 0.0f }};
	long sx =100;
	long sy =200;
	
	v4f_t xmm0, xmm1;
	xmm0.v_xmm = _mm_setr_ps(sx, sy, sx, sy);//_mm_cvtsi32_ss( xmm0.v_xmm, 123 );
	xmm1.v_xmm = _mm_setr_ps(12, 3, 14, 5);

	printf("= vXMM %f %f %f %f\n",xmm0.xf, xmm0.yf, xmm0.zf, xmm0.wf);
	xmm0.v_xmm = _mm_shuffle_ps(xmm0.v_xmm, xmm0.v_xmm, 0 );
    printf("s vXMM %f %f %f %f\n",xmm0.xf, xmm0.yf, xmm0.zf, xmm0.wf);

	xmm0.v_xmm = _mm_mul_ps( xmm0.v_xmm, *(__m128*)&optistrx );
	printf("* vXMM %f %f %f %f\n",xmm0.xf, xmm0.yf, xmm0.zf, xmm0.wf);
	xmm0.v_xmm = _mm_add_ps( xmm0.v_xmm, xmm1.v_xmm );
    printf("+ vXMM %f %f %f %f\n",xmm0.xf, xmm0.yf, xmm0.zf, xmm0.wf);


	__m128 dist = _mm_load_ss( &d );
	float hadd4result =0.0f;
	_mm_store_ps ( &hadd4result, _mm_mul_ps( _mm_rsqrt( _mm_hadd4_ps( _mm_mul_ps( vect.v_xmm, vect2.v_xmm) ) ), dist ) );
	uint32_t morton = EncodeMorton3(916,(1<<10)-2,617);
	printf ( "sqrt (100)= %f ; fsqrt (100)= %f hadd4result=%f\n" , 300.0f/sqrt(x*x+y*y), f_rsqrt(100.0f)*300.0f , hadd4result );
	printf("%d %d %d ", DecodeMorton3X(morton), DecodeMorton3Y(morton), DecodeMorton3Z(morton));
	//b = a * scalar;
	// assert( b != a );
	//assert( b.x == a.x * scalar && b.y == a.y * scalar && b.z == a.z * scalar && b.w == a.w * scalar );

	ss_bbox_t bbmax, bbmin;
	ss_bbox_t bba = (ss_bbox_t){{ 10, 20, 30, 40, 50, 60, 70, 80 }};
	ss_bbox_t bbb = (ss_bbox_t){{ 10, 1442, 152, 56, 50, (1<<5)-1, -152, 1442 }};
	
  	// Minima and maxima of 8 signed short values
  	bbmax.v_xmm = _mm_max_epi16 ( bba.v_xmm , bbb.v_xmm );
  	bbmin.v_xmm = _mm_min_epi16 ( bba.v_xmm , bbb.v_xmm );

   	// Minima and maxima of 4 signed short values
   	bbmin.bmin.v_vec = _mm_min_pi16( bba.bmin.v_vec , bbb.bmin.v_vec );
   	bbmin.bmax.v_vec = _mm_min_pi16( bba.bmax.v_vec , bbb.bmax.v_vec );

	bbmin.bmax.dump();
	bba.v8_vec.dump();
	bbb.v8_vec.dump();
	bbmax.v8_vec.dump();
	bbmin.v8_vec.dump();

	/*

	printf("MAX(%hi,%hi) == %hi, MAX(%hi,%hi) == %hi, MAX(%hi,%hi) == %hi, MAX(%hi,%hi) == %hi MAX(%hi,%hi) == %hi, MAX(%hi,%hi) == %hi, MAX(%hi,%hi) == %hi, MAX(%hi,%hi) == %hi\n"
		, bba.x1, bbb.x1, bbmax.x1, bba.y1, bbb.y1, bbmax.y1, bba.z1, bbb.z1, bbmax.z1, bba.pad, bbb.pad, bbmax.pad
		, bba.x2, bbb.x2, bbmax.x2, bba.y2, bbb.y2, bbmax.y2, bba.z2, bbb.z2, bbmax.z2, bba.csgdel, bbb.csgdel, bbmax.csgdel
	);

	printf("MIN(%hi,%hi) == %hi, MIN(%hi,%hi) == %hi, MIN(%hi,%hi) == %hi, MIN(%hi,%hi) == %hi MIN(%hi,%hi) == %hi, MIN(%hi,%hi) == %hi, MIN(%hi,%hi) == %hi, MIN(%hi,%hi) == %hi\n"
		, bba.x1, bbb.x1, bbmin.x1, bba.y1, bbb.y1, bbmin.y1, bba.z1, bbb.z1, bbmin.z1, bba.pad, bbb.pad, bbmin.pad
		, bba.x2, bbb.x2, bbmin.x2, bba.y2, bbb.y2, bbmin.y2, bba.z2, bbb.z2, bbmin.z2, bba.csgdel, bbb.csgdel, bbmin.csgdel
	);
	printf("bbmin %hi\n"
		, bbmin.z1
	);
	*/
 // }
  
}