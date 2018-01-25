#ifndef __WISARD_SIMD_H__
#define __WISARD_SIMD_H__
#pragma once

#ifndef _WIN32
#	include <xmmintrin.h>
#endif
#include "global/common.h"
#define USE_SSE2

namespace ONETOOL {

/* Disable of Windows-dependent warnings */
#ifdef WIN32
#	ifdef USE_SSE
#		pragma warning(disable:4003)
#	endif
#endif

#ifdef USE_AVX
#	ifdef _M_ARM
#		error			"ARM architecture does not support AVX"
#	endif
#	include <immintrin.h>
#	define szSIMD(v)			((v)<<1)
#	define funcSIMD(t)			_mm256_##t
#	define func128SIMD(t)		_mm_##t
#	define typeSIMD(t)			__m256##t
#else
#	ifdef _M_ARM
#		define szSIMD(v)		(v)
#		define funcSIMD(t)		v##t
#		define typeSIMD()		float32x4_t
#	else
#		define szSIMD(v)		(v)
#		define funcSIMD(t)		_mm_##t
#		define func128SIMD(t)	_mm_##t
#		define typeSIMD(t)		__m128##t
#	endif
#endif

#ifdef USE_SSE
#	define sse128store		func128SIMD(store_si128)
#	define sse128zero		func128SIMD(set1_epi8)(0)
#	define sse128and		func128SIMD(and_si128)
#	define sse128xor		func128SIMD(xor_si128)
#	define sse128or			func128SIMD(or_si128)
#	ifdef _M_ARM
#		include <arm_neon.h>
#		define getMed(a)		(((a)>>szSIMD(2))<<szSIMD(2))
#		define sseJmp			szSIMD(4)
#		define sse_t			typeSIMD()
#		define sseAdd			funcSIMD(addq_f32)
#		define sseSub			funcSIMD(subq_f32)
#		define sseMul			funcSIMD(mulq_f32)
#		define sseDiv(a,b)		neon_div(a,b)
#		define sseStore			funcSIMD(st1q_f32)
#		define sseSet(a)		funcSIMD(dupq_n_f32)((float)a)
#		define sseNot			funcSIMD(mvnq_u32)
#		define sseAndNot(a,b)	funcSIMD(andq_u32)(sseNot(a), b)
#		define sseNeq(a,b)		sseNot(funcSIMD(ceq_f32)(a, b))
#		define sseAnd(a,b)		neon_and(a,b)
#		define sseSqrt
#		define sseXor(a,b)		neon_xor(a,b)
#		ifdef USE_DBL
#			undef USE_DBL
#		endif
sse_t	neon_div(sse_t a, sse_t b);
sse_t	neon_and(sse_t a, sse_t b);
sse_t	neon_xor(sse_t a, sse_t b);
#		define MASK_ON			(wsReal)(0xffffff)
#		define MASK_OFF			(wsReal)(0x000000)
#		if defined(_M_X64)
#			error		"ARM architecture does not support x64"
#		endif
#	else
#		include <xmmintrin.h>
#		include <emmintrin.h>
#		include <tmmintrin.h>
#		ifdef _WIN32
#			include <nmmintrin.h>
#			ifdef _M_X64
#				define WISARD_ARCH u64
#			else
#				define WISARD_ARCH u32
#			endif
#		elif defined(USE_SSE4)
#			include <smmintrin.h>
#		endif
#		ifdef USE_DBL
#			define getMed(a)		(((a)>>szSIMD(1))<<szSIMD(1))
#			define sseJmp			szSIMD(2)
#			define sse_t			typeSIMD(d)
#			define sseAdd			funcSIMD(add_pd)
#			define sseSub			funcSIMD(sub_pd)
#			define sseMul			funcSIMD(mul_pd)
#			define sseDiv			funcSIMD(div_pd)
#			define sseStore			funcSIMD(store_pd)
#			define sseLoad			funcSIMD(load_pd)
#			define sseSet(a)		funcSIMD(set1_pd)((wsReal)a)
#			define sseAndNot		funcSIMD(andnot_pd)
#			ifdef _CMP_NEQ_UQ
#				define sseNeq(a,b)	funcSIMD(cmp_pd)(a, b, _CMP_NEQ_UQ)
#			else
#				define sseNeq		funcSIMD(cmpneq_pd)
#			endif
#			define sseLes			funcSIMD(cmplt_pd)
#			define sseAnd			funcSIMD(and_pd)
#			define sseSqrt			funcSIMD(sqrt_pd)
#			define sseXor			funcSIMD(xor_pd)
#			define MASK_ON			(wsReal)(0xffffffffffffffff)
#			define MASK_OFF			(wsReal)(0x0000000000000000)
#		else
#			define getMed(a)		(((a)>>szSIMD(2))<<szSIMD(2))
#			define sseJmp			szSIMD(4)
#			define sse_t			typeSIMD()
#			define sseAdd			funcSIMD(add_ps)
#			define sseSub			funcSIMD(sub_ps)
#			define sseMul			funcSIMD(mul_ps)
#			define sseDiv			funcSIMD(div_ps)
#			define sseStore			funcSIMD(store_ps)
#			define sseLoad			funcSIMD(load_ps)
#			define sseSet(a)		funcSIMD(set1_ps)((wsReal)a)
#			define sseAndNot		funcSIMD(andnot_ps)
#			ifdef _CMP_NEQ_UQ
#				define sseNeq(a,b)	funcSIMD(cmp_ps)(a, b, _CMP_NEQ_UQ)
#			else
#				define sseNeq		funcSIMD(cmpneq_ps)
#			endif
#			define sseLes			funcSIMD(cmplt_ps)
#			define sseAnd			funcSIMD(and_ps)
#			define sseSqrt			funcSIMD(sqrt_ps)
#			define sseXor			funcSIMD(xor_ps)
#			define MASK_ON			(wsReal)(0xffffffff)
#			define MASK_OFF			(wsReal)(0x00000000)
#		endif
#	endif
#	define get16(a)		(((a)>>4)<<4)
#endif

#ifdef _WIN32
#	if defined(_M_X64)
#		define WISARD_POPCNT _mm_popcnt_u64
#		define WISARD_pc	unsigned __int64
#		define WISARD_szpc	64
#		define WISARD_sftpc	6
#		define WISARD32		0
#	elif defined(_M_ARM)
#		define WISARD_POPCNT _mm_popcnt_u32
int _mm_popcnt_u32(unsigned int);
#		define WISARD_pc	unsigned int
#		define WISARD_szpc	32
#		define WISARD_sftpc	5
#		define WISARD32		1
#	else
#		define WISARD_POPCNT _mm_popcnt_u32
#		define WISARD_pc	unsigned int
#		define WISARD_szpc	32
#		define WISARD_sftpc	5
#		define WISARD32		1
#	endif
#else
#	if defined(__alpha__) || defined(__ia64__) || defined(__x86_64__)
#		define WISARD_POPCNT __builtin_popcountll
#		define WISARD_pc	unsigned long long int
#		define WISARD_szpc	64
#		define WISARD_sftpc	6
#		define WISARD32		0
#	else
#		define WISARD_POPCNT __builtin_popcount
#		define WISARD_pc	unsigned int
#		define WISARD_szpc	32
#		define WISARD_sftpc	5
#		define WISARD32		1
#	endif
#endif

#ifdef USE_AVX
#	define alignSSE(a) ((char *)(a) + (0x20-(0x0000001f&(size_t)a)))
#else
#	define alignSSE(a) ((char *)(a) + (0x10-(0x0000000f&(size_t)a)))
#endif

#define sseSum(s, d)		sseSumType(s, d, wsReal, sseJmp, sseStore, wsReal)
#define sseProd(s, d)		sseProdType(s, d, wsReal, sseJmp, sseStore, wsReal)
#define sseUsum(s, d)		sseSumType(s, d, wsUint, sseJmp, sseStore, wsReal)
#define sseAppendSum(s, d)	sseAppendSumType(s, d, wsReal, sseJmp, sseStore, wsReal)
#define sseAppendUsum(s, d)	sseAppendSumType(s, d, wsUint, sseJmp, sseStore, wsReal)
#define sse8SHORTasum(s, d)	{ unsigned short Na_sseSumBuf[32], *Rp_sseStorePtr; \
	Rp_sseStorePtr = (unsigned short *)alignSSE(Na_sseSumBuf); \
	_mm_store_si128((__m128i *)Rp_sseStorePtr, s); \
	d += Rp_sseStorePtr[0]; \
	d += Rp_sseStorePtr[1]; \
	d += Rp_sseStorePtr[2]; \
	d += Rp_sseStorePtr[3]; \
	d += Rp_sseStorePtr[4]; \
	d += Rp_sseStorePtr[5]; \
	d += Rp_sseStorePtr[6]; \
	d += Rp_sseStorePtr[7]; }
#define sse16CHARasum(s, d)	{ unsigned char Na_sseSumBuf[64], *Rp_sseStorePtr; \
	Rp_sseStorePtr = (unsigned char *)alignSSE(Na_sseSumBuf); \
	_mm_store_si128((__m128i *)Rp_sseStorePtr, s); \
	d += Rp_sseStorePtr[0]; \
	d += Rp_sseStorePtr[1]; \
	d += Rp_sseStorePtr[2]; \
	d += Rp_sseStorePtr[3]; \
	d += Rp_sseStorePtr[4]; \
	d += Rp_sseStorePtr[5]; \
	d += Rp_sseStorePtr[6]; \
	d += Rp_sseStorePtr[7]; \
	d += Rp_sseStorePtr[8]; \
	d += Rp_sseStorePtr[9]; \
	d += Rp_sseStorePtr[10]; \
	d += Rp_sseStorePtr[11]; \
	d += Rp_sseStorePtr[12]; \
	d += Rp_sseStorePtr[13]; \
	d += Rp_sseStorePtr[14]; \
	d += Rp_sseStorePtr[15]; }

#define sseSumType(s, d, t, j, x, p) { p Ra_sseSumBuf[j<<1], *Rp_sseStorePtr; \
	Rp_sseStorePtr = (p *)alignSSE(Ra_sseSumBuf); \
	x(Rp_sseStorePtr, s); \
	d = (t)Rp_sseStorePtr[0]; \
	for (wsUint X=1 ; X<j ; X++) \
	d += (t)Rp_sseStorePtr[X]; }
#define sseProdType(s, d, t, j, x, p) { p Ra_sseSumBuf[j<<1], *Rp_sseStorePtr; \
	Rp_sseStorePtr = (p *)alignSSE(Ra_sseSumBuf); \
	x(Rp_sseStorePtr, s); \
	d = (t)Rp_sseStorePtr[0]; \
	for (wsUint X=1 ; X<j ; X++) \
	d *= (t)Rp_sseStorePtr[X]; }
#define sseAppendSumType(s, d, t, j, x, p) { p Ra_sseSumBuf[j<<1], *Rp_sseStorePtr; \
	Rp_sseStorePtr = (p *)alignSSE(Ra_sseSumBuf); \
	x(Rp_sseStorePtr, s); \
	d += (t)Rp_sseStorePtr[0]; \
	for (wsUint X=1 ; X<j ; X++) \
	d += (t)Rp_sseStorePtr[X]; }

typedef float wsFloat;
typedef	const wsFloat wsFloatCst;
typedef wsFloat* wsFvec;
typedef wsFvec* wsFmat;
#define	CORR_CONST(X_l)	X_l##f
#define	corrsse_t				typeSIMD()
#ifdef _M_ARM
#	define corrsseSum(s, d)		sseSumType(s, d, wsFloat, corrsseJmp, funcSIMD(st1q_f32), wsFloat)
#	define corrsseUsum(s, d)	sseSumType(s, d, wsUint, corrsseJmp, funcSIMD(st1q_f32), wsFloat)
#else
#	define corrsseSum(s, d)		sseSumType(s, d, wsFloat, corrsseJmp, funcSIMD(store_ps), wsFloat)
#	define corrsseUsum(s, d)	sseSumType(s, d, wsUint, corrsseJmp, funcSIMD(store_ps), wsFloat)
#endif
#define	corrGetMed(a)			(((a)>>szSIMD(2))<<szSIMD(2))
#define	corrsseJmp				szSIMD(4)
#ifdef _M_ARM
#	define	corrsseSet(a)		funcSIMD(dupq_n_f32)((wsFloat)a)
#	define	corrsseAdd			funcSIMD(addq_f32)
#	define	corrsseMul			funcSIMD(mulq_f32)
#	define	corrsseNeq(a,b)		sseNot(funcSIMD(ceq_f32)(a, b))
#	define	corrsseAnd(a,b)		neon_and(a,b)
#else
#	define	corrsseSet(a)		funcSIMD(set1_ps)((wsFloat)a)
#	define	corrsseAdd			funcSIMD(add_ps)
#	define	corrsseMul			funcSIMD(mul_ps)
#	ifdef _CMP_NEQ_UQ
#		define corrsseNeq(a,b)	funcSIMD(cmp_ps)(a, b, _CMP_NEQ_UQ)
#	else
#		define corrsseNeq		funcSIMD(cmpneq_ps)
#	endif
#	define corrsseAnd			funcSIMD(and_ps)
inline void sse2corr(wsReal *Ra_sse, wsFloat *Ra_cor, wsUint N_sz) {
	wsUint N_med = (N_sz >> 2) << 2;
	/* Should jump per 32 bytes (8 bytes * 4 elements) */
	for (wsUint i=0 ; i<N_med ; i+=4) {
		sse_t *v0 = (sse_t *)(Ra_sse + i);
		sse_t *v1 = (sse_t *)(Ra_sse + i + 2);
		__m128 *vt = (__m128 *)(Ra_cor + i);

		*vt = _mm_movelh_ps (_mm_cvtpd_ps(*v0), _mm_cvtpd_ps(*v1));
	}
	for (wsUint i=N_med ; i<N_sz ; i++)
		Ra_cor[i] = (wsFloat)Ra_sse[i];
}
#endif

/* SIMD (SSE1+MMX or SSE2) implementation of sin, cos, exp and log

   Inspired by Intel Approximate Math library, and based on the
   corresponding algorithms of the cephes math library

   The default is to use the SSE1 version. If you define USE_SSE2 the
   the SSE2 intrinsics will be used in place of the MMX intrinsics. Do
   not expect any significant performance improvement with SSE2.
*/

/* Copyright (C) 2007  Julien Pommier

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  (this is the zlib license)
*/

/* yes I know, the top of this file is quite ugly */

#ifdef _MSC_VER /* visual c++ */
#	define ALIGN16_BEG __declspec(align(16))
#	define ALIGN16_END 
#else /* gcc or icc */
#	define ALIGN16_BEG
#	define ALIGN16_END __attribute__((aligned(16)))
#endif

/* __m128 is ugly to write */
typedef __m128 v4sf;  // vector of 4 float (sse1)

#ifdef USE_SSE2
#	include <emmintrin.h>
typedef __m128i v4si; // vector of 4 int (sse2)
#else
typedef __m64 v2si;   // vector of 2 int (mmx)
#endif

/* declare some SSE constants -- why can't I figure a better way to do that? */
#define _PS_CONST(Name, Val)                                            \
  static const ALIGN16_BEG float _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PI32_CONST(Name, Val)                                            \
  static const ALIGN16_BEG int _pi32_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PS_CONST_TYPE(Name, Type, Val)                                 \
  static const ALIGN16_BEG Type _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }

_PS_CONST(1  , 1.0f);
_PS_CONST(0p5, 0.5f);
/* the smallest non denormalized float number */
_PS_CONST_TYPE(min_norm_pos, int, 0x00800000);
_PS_CONST_TYPE(mant_mask, int, 0x7f800000);
_PS_CONST_TYPE(inv_mant_mask, int, ~0x7f800000);

_PS_CONST_TYPE(sign_mask, int, (int)0x80000000);
_PS_CONST_TYPE(inv_sign_mask, int, ~0x80000000);

_PI32_CONST(1, 1);
_PI32_CONST(inv1, ~1);
_PI32_CONST(2, 2);
_PI32_CONST(4, 4);
_PI32_CONST(0x7f, 0x7f);

_PS_CONST(cephes_SQRTHF, 0.707106781186547524f);
_PS_CONST(cephes_log_p0, 7.0376836292E-2f);
_PS_CONST(cephes_log_p1, - 1.1514610310E-1f);
_PS_CONST(cephes_log_p2, 1.1676998740E-1f);
_PS_CONST(cephes_log_p3, - 1.2420140846E-1f);
_PS_CONST(cephes_log_p4, + 1.4249322787E-1f);
_PS_CONST(cephes_log_p5, - 1.6668057665E-1f);
_PS_CONST(cephes_log_p6, + 2.0000714765E-1f);
_PS_CONST(cephes_log_p7, - 2.4999993993E-1f);
_PS_CONST(cephes_log_p8, + 3.3333331174E-1f);
_PS_CONST(cephes_log_q1, -2.12194440e-4f);
_PS_CONST(cephes_log_q2, 0.693359375);

#ifndef USE_SSE2
typedef union xmm_mm_union {
  __m128 xmm;
  __m64 mm[2];
} xmm_mm_union;

#define COPY_XMM_TO_MM(xmm_, mm0_, mm1_) {          \
    xmm_mm_union u; u.xmm = xmm_;                   \
    mm0_ = u.mm[0];                                 \
    mm1_ = u.mm[1];                                 \
}

#define COPY_MM_TO_XMM(mm0_, mm1_, xmm_) {                         \
    xmm_mm_union u; u.mm[0]=mm0_; u.mm[1]=mm1_; xmm_ = u.xmm;      \
  }

#endif // USE_SSE2

/* natural logarithm computed for 4 simultaneous float 
   return NaN for x <= 0
*/
inline v4sf log_ps(v4sf x) {
#ifdef USE_SSE2
	v4si emm0;
#else
	v2si mm0, mm1;
#endif
	v4sf one = *(v4sf*)_ps_1;
	v4sf invalid_mask = _mm_cmple_ps(x, _mm_setzero_ps());

	x = _mm_max_ps(x, *(v4sf*)_ps_min_norm_pos);  /* cut off denormalized stuff */

#ifndef USE_SSE2
	/* part 1: x = frexpf(x, &e); */
	COPY_XMM_TO_MM(x, mm0, mm1);
	mm0 = _mm_srli_pi32(mm0, 23);
	mm1 = _mm_srli_pi32(mm1, 23);
#else
	emm0 = _mm_srli_epi32(_mm_castps_si128(x), 23);
#endif
	/* keep only the fractional part */
	x = _mm_and_ps(x, *(v4sf*)_ps_inv_mant_mask);
	x = _mm_or_ps(x, *(v4sf*)_ps_0p5);

#ifndef USE_SSE2
	/* now e=mm0:mm1 contain the really base-2 exponent */
	mm0 = _mm_sub_pi32(mm0, *(v2si*)_pi32_0x7f);
	mm1 = _mm_sub_pi32(mm1, *(v2si*)_pi32_0x7f);
	v4sf e = _mm_cvtpi32x2_ps(mm0, mm1);
	_mm_empty(); /* bye bye mmx */
#else
	emm0 = _mm_sub_epi32(emm0, *(v4si*)_pi32_0x7f);
	v4sf e = _mm_cvtepi32_ps(emm0);
#endif

	e = _mm_add_ps(e, one);

	/* part2: 
	if( x < SQRTHF ) {
		e -= 1;
		x = x + x - 1.0;
	} else { x = x - 1.0; }
	*/
	v4sf mask	= _mm_cmplt_ps(x, *(v4sf*)_ps_cephes_SQRTHF);
	v4sf tmp	= _mm_and_ps(x, mask);
	x			= _mm_sub_ps(x, one);
	e			= _mm_sub_ps(e, _mm_and_ps(one, mask));
	x			= _mm_add_ps(x, tmp);


	v4sf z		= _mm_mul_ps(x,x);
	v4sf y		= *(v4sf*)_ps_cephes_log_p0;
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_log_p1);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_log_p2);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_log_p3);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_log_p4);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_log_p5);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_log_p6);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_log_p7);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_log_p8);
	y = _mm_mul_ps(y, x);

	y = _mm_mul_ps(y, z);
  

	tmp = _mm_mul_ps(e, *(v4sf*)_ps_cephes_log_q1);
	y = _mm_add_ps(y, tmp);


	tmp = _mm_mul_ps(z, *(v4sf*)_ps_0p5);
	y = _mm_sub_ps(y, tmp);

	tmp = _mm_mul_ps(e, *(v4sf*)_ps_cephes_log_q2);
	x = _mm_add_ps(x, y);
	x = _mm_add_ps(x, tmp);
	x = _mm_or_ps(x, invalid_mask); // negative arg will be NAN
	return x;
}

_PS_CONST(exp_hi,	88.3762626647949f);
_PS_CONST(exp_lo,	-88.3762626647949f);

_PS_CONST(cephes_LOG2EF, 1.44269504088896341f);
_PS_CONST(cephes_exp_C1, 0.693359375);
_PS_CONST(cephes_exp_C2, -2.12194440e-4f);

_PS_CONST(cephes_exp_p0, 1.9875691500E-4f);
_PS_CONST(cephes_exp_p1, 1.3981999507E-3f);
_PS_CONST(cephes_exp_p2, 8.3334519073E-3f);
_PS_CONST(cephes_exp_p3, 4.1665795894E-2f);
_PS_CONST(cephes_exp_p4, 1.6666665459E-1f);
_PS_CONST(cephes_exp_p5, 5.0000001201E-1f);

inline v4sf exp_ps(v4sf x) {
	v4sf tmp = _mm_setzero_ps(), fx;
#ifdef USE_SSE2
	v4si emm0;
#else
	v2si mm0, mm1;
#endif
	v4sf one = *(v4sf*)_ps_1;

	x = _mm_min_ps(x, *(v4sf*)_ps_exp_hi);
	x = _mm_max_ps(x, *(v4sf*)_ps_exp_lo);

	/* express exp(x) as exp(g + n*log(2)) */
	fx = _mm_mul_ps(x, *(v4sf*)_ps_cephes_LOG2EF);
	fx = _mm_add_ps(fx, *(v4sf*)_ps_0p5);

	/* how to perform a floorf with SSE: just below */
#ifndef USE_SSE2
	/* step 1 : cast to int */
	tmp = _mm_movehl_ps(tmp, fx);
	mm0 = _mm_cvttps_pi32(fx);
	mm1 = _mm_cvttps_pi32(tmp);
	/* step 2 : cast back to float */
	tmp = _mm_cvtpi32x2_ps(mm0, mm1);
#else
	emm0 = _mm_cvttps_epi32(fx);
	tmp  = _mm_cvtepi32_ps(emm0);
#endif
	/* if greater, substract 1 */
	v4sf mask = _mm_cmpgt_ps(tmp, fx);    
	mask = _mm_and_ps(mask, one);
	fx = _mm_sub_ps(tmp, mask);

	tmp = _mm_mul_ps(fx, *(v4sf*)_ps_cephes_exp_C1);
	v4sf z = _mm_mul_ps(fx, *(v4sf*)_ps_cephes_exp_C2);
	x = _mm_sub_ps(x, tmp);
	x = _mm_sub_ps(x, z);

	z = _mm_mul_ps(x,x);
  
	v4sf y = *(v4sf*)_ps_cephes_exp_p0;
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p1);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p2);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p3);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p4);
	y = _mm_mul_ps(y, x);
	y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p5);
	y = _mm_mul_ps(y, z);
	y = _mm_add_ps(y, x);
	y = _mm_add_ps(y, one);

  /* build 2^n */
#ifndef USE_SSE2
	z = _mm_movehl_ps(z, fx);
	mm0 = _mm_cvttps_pi32(fx);
	mm1 = _mm_cvttps_pi32(z);
	mm0 = _mm_add_pi32(mm0, *(v2si*)_pi32_0x7f);
	mm1 = _mm_add_pi32(mm1, *(v2si*)_pi32_0x7f);
	mm0 = _mm_slli_pi32(mm0, 23); 
	mm1 = _mm_slli_pi32(mm1, 23);
  
	v4sf pow2n; 
	COPY_MM_TO_XMM(mm0, mm1, pow2n);
	_mm_empty();
#else
	emm0 = _mm_cvttps_epi32(fx);
	emm0 = _mm_add_epi32(emm0, *(v4si*)_pi32_0x7f);
	emm0 = _mm_slli_epi32(emm0, 23);
	v4sf pow2n = _mm_castsi128_ps(emm0);
#endif
	y = _mm_mul_ps(y, pow2n);
	return y;
}

_PS_CONST(minus_cephes_DP1, -0.78515625);
_PS_CONST(minus_cephes_DP2, -2.4187564849853515625e-4);
_PS_CONST(minus_cephes_DP3, -3.77489497744594108e-8f);
_PS_CONST(sincof_p0, -1.9515295891E-4f);
_PS_CONST(sincof_p1,  8.3321608736E-3f);
_PS_CONST(sincof_p2, -1.6666654611E-1f);
_PS_CONST(coscof_p0,  2.443315711809948E-005f);
_PS_CONST(coscof_p1, -1.388731625493765E-003f);
_PS_CONST(coscof_p2,  4.166664568298827E-002f);
_PS_CONST(cephes_FOPI, 1.27323954473516f); // 4 / M_PI


/* evaluation of 4 sines at onces, using only SSE1+MMX intrinsics so
   it runs also on old athlons XPs and the pentium III of your grand
   mother.

   The code is the exact rewriting of the cephes sinf function.
   Precision is excellent as long as x < 8192 (I did not bother to
   take into account the special handling they have for greater values
   -- it does not return garbage for arguments over 8192, though, but
   the extra precision is missing).

   Note that it is such that sinf((float)M_PI) = 8.74e-8, which is the
   surprising but correct result.

   Performance is also surprisingly good, 1.33 times faster than the
   macos vsinf SSE2 function, and 1.5 times faster than the
   __vrs4_sinf of amd's ACML (which is only available in 64 bits). Not
   too bad for an SSE1 function (with no special tuning) !
   However the latter libraries probably have a much better handling of NaN,
   Inf, denormalized and other special arguments..

   On my core 1 duo, the execution of this function takes approximately 95 cycles.

   From what I have observed on the experiments with Intel AMath lib, switching to an
   SSE2 version would improve the perf by only 10%.

   Since it is based on SSE intrinsics, it has to be compiled at -O2 to
   deliver full speed.
*/
inline v4sf sin_ps(v4sf x) { // any x
	v4sf xmm1, xmm2 = _mm_setzero_ps(), xmm3, sign_bit, y;

#ifdef USE_SSE2
	v4si emm0, emm2;
#else
	v2si mm0, mm1, mm2, mm3;
#endif
	sign_bit = x;
	/* take the absolute value */
	x = _mm_and_ps(x, *(v4sf*)_ps_inv_sign_mask);
	/* extract the sign bit (upper one) */
	sign_bit = _mm_and_ps(sign_bit, *(v4sf*)_ps_sign_mask);
  
	/* scale by 4/Pi */
	y = _mm_mul_ps(x, *(v4sf*)_ps_cephes_FOPI);

#ifdef USE_SSE2
	/* store the integer part of y in mm0 */
	emm2 = _mm_cvttps_epi32(y);
	/* j=(j+1) & (~1) (see the cephes sources) */
	emm2 = _mm_add_epi32(emm2, *(v4si*)_pi32_1);
	emm2 = _mm_and_si128(emm2, *(v4si*)_pi32_inv1);
	y = _mm_cvtepi32_ps(emm2);

	/* get the swap sign flag */
	emm0 = _mm_and_si128(emm2, *(v4si*)_pi32_4);
	emm0 = _mm_slli_epi32(emm0, 29);
	/* get the polynom selection mask 
		there is one polynom for 0 <= x <= Pi/4
		and another one for Pi/4<x<=Pi/2

		Both branches will be computed.
	*/
	emm2 = _mm_and_si128(emm2, *(v4si*)_pi32_2);
	emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
  
	v4sf swap_sign_bit = _mm_castsi128_ps(emm0);
	v4sf poly_mask = _mm_castsi128_ps(emm2);
	sign_bit = _mm_xor_ps(sign_bit, swap_sign_bit);
  
#else
	/* store the integer part of y in mm0:mm1 */
	xmm2 = _mm_movehl_ps(xmm2, y);
	mm2 = _mm_cvttps_pi32(y);
	mm3 = _mm_cvttps_pi32(xmm2);
	/* j=(j+1) & (~1) (see the cephes sources) */
	mm2 = _mm_add_pi32(mm2, *(v2si*)_pi32_1);
	mm3 = _mm_add_pi32(mm3, *(v2si*)_pi32_1);
	mm2 = _mm_and_si64(mm2, *(v2si*)_pi32_inv1);
	mm3 = _mm_and_si64(mm3, *(v2si*)_pi32_inv1);
	y = _mm_cvtpi32x2_ps(mm2, mm3);
	/* get the swap sign flag */
	mm0 = _mm_and_si64(mm2, *(v2si*)_pi32_4);
	mm1 = _mm_and_si64(mm3, *(v2si*)_pi32_4);
	mm0 = _mm_slli_pi32(mm0, 29);
	mm1 = _mm_slli_pi32(mm1, 29);
	/* get the polynom selection mask */
	mm2 = _mm_and_si64(mm2, *(v2si*)_pi32_2);
	mm3 = _mm_and_si64(mm3, *(v2si*)_pi32_2);
	mm2 = _mm_cmpeq_pi32(mm2, _mm_setzero_si64());
	mm3 = _mm_cmpeq_pi32(mm3, _mm_setzero_si64());
	v4sf swap_sign_bit, poly_mask;
	COPY_MM_TO_XMM(mm0, mm1, swap_sign_bit);
	COPY_MM_TO_XMM(mm2, mm3, poly_mask);
	sign_bit = _mm_xor_ps(sign_bit, swap_sign_bit);
	_mm_empty(); /* good-bye mmx */
#endif
  
	/* The magic pass: "Extended precision modular arithmetic" 
		x = ((x - y * DP1) - y * DP2) - y * DP3; */
	xmm1 = *(v4sf*)_ps_minus_cephes_DP1;
	xmm2 = *(v4sf*)_ps_minus_cephes_DP2;
	xmm3 = *(v4sf*)_ps_minus_cephes_DP3;
	xmm1 = _mm_mul_ps(y, xmm1);
	xmm2 = _mm_mul_ps(y, xmm2);
	xmm3 = _mm_mul_ps(y, xmm3);
	x = _mm_add_ps(x, xmm1);
	x = _mm_add_ps(x, xmm2);
	x = _mm_add_ps(x, xmm3);

	/* Evaluate the first polynom  (0 <= x <= Pi/4) */
	y = *(v4sf*)_ps_coscof_p0;
	v4sf z = _mm_mul_ps(x,x);

	y = _mm_mul_ps(y, z);
	y = _mm_add_ps(y, *(v4sf*)_ps_coscof_p1);
	y = _mm_mul_ps(y, z);
	y = _mm_add_ps(y, *(v4sf*)_ps_coscof_p2);
	y = _mm_mul_ps(y, z);
	y = _mm_mul_ps(y, z);
	v4sf tmp = _mm_mul_ps(z, *(v4sf*)_ps_0p5);
	y = _mm_sub_ps(y, tmp);
	y = _mm_add_ps(y, *(v4sf*)_ps_1);
  
	/* Evaluate the second polynom  (Pi/4 <= x <= 0) */

	v4sf y2 = *(v4sf*)_ps_sincof_p0;
	y2 = _mm_mul_ps(y2, z);
	y2 = _mm_add_ps(y2, *(v4sf*)_ps_sincof_p1);
	y2 = _mm_mul_ps(y2, z);
	y2 = _mm_add_ps(y2, *(v4sf*)_ps_sincof_p2);
	y2 = _mm_mul_ps(y2, z);
	y2 = _mm_mul_ps(y2, x);
	y2 = _mm_add_ps(y2, x);

	/* select the correct result from the two polynoms */  
	xmm3 = poly_mask;
	y2 = _mm_and_ps(xmm3, y2); //, xmm3);
	y = _mm_andnot_ps(xmm3, y);
	y = _mm_add_ps(y,y2);
	/* update the sign */
	y = _mm_xor_ps(y, sign_bit);
	return y;
}

/* almost the same as sin_ps */
inline v4sf cos_ps(v4sf x) { // any x
  v4sf xmm1, xmm2 = _mm_setzero_ps(), xmm3, y;
#ifdef USE_SSE2
  v4si emm0, emm2;
#else
  v2si mm0, mm1, mm2, mm3;
#endif
  /* take the absolute value */
  x = _mm_and_ps(x, *(v4sf*)_ps_inv_sign_mask);
  
  /* scale by 4/Pi */
  y = _mm_mul_ps(x, *(v4sf*)_ps_cephes_FOPI);
  
#ifdef USE_SSE2
  /* store the integer part of y in mm0 */
  emm2 = _mm_cvttps_epi32(y);
  /* j=(j+1) & (~1) (see the cephes sources) */
  emm2 = _mm_add_epi32(emm2, *(v4si*)_pi32_1);
  emm2 = _mm_and_si128(emm2, *(v4si*)_pi32_inv1);
  y = _mm_cvtepi32_ps(emm2);

  emm2 = _mm_sub_epi32(emm2, *(v4si*)_pi32_2);
  
  /* get the swap sign flag */
  emm0 = _mm_andnot_si128(emm2, *(v4si*)_pi32_4);
  emm0 = _mm_slli_epi32(emm0, 29);
  /* get the polynom selection mask */
  emm2 = _mm_and_si128(emm2, *(v4si*)_pi32_2);
  emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
  
  v4sf sign_bit = _mm_castsi128_ps(emm0);
  v4sf poly_mask = _mm_castsi128_ps(emm2);
#else
  /* store the integer part of y in mm0:mm1 */
  xmm2 = _mm_movehl_ps(xmm2, y);
  mm2 = _mm_cvttps_pi32(y);
  mm3 = _mm_cvttps_pi32(xmm2);

  /* j=(j+1) & (~1) (see the cephes sources) */
  mm2 = _mm_add_pi32(mm2, *(v2si*)_pi32_1);
  mm3 = _mm_add_pi32(mm3, *(v2si*)_pi32_1);
  mm2 = _mm_and_si64(mm2, *(v2si*)_pi32_inv1);
  mm3 = _mm_and_si64(mm3, *(v2si*)_pi32_inv1);

  y = _mm_cvtpi32x2_ps(mm2, mm3);


  mm2 = _mm_sub_pi32(mm2, *(v2si*)_pi32_2);
  mm3 = _mm_sub_pi32(mm3, *(v2si*)_pi32_2);

  /* get the swap sign flag in mm0:mm1 and the 
     polynom selection mask in mm2:mm3 */

  mm0 = _mm_andnot_si64(mm2, *(v2si*)_pi32_4);
  mm1 = _mm_andnot_si64(mm3, *(v2si*)_pi32_4);
  mm0 = _mm_slli_pi32(mm0, 29);
  mm1 = _mm_slli_pi32(mm1, 29);

  mm2 = _mm_and_si64(mm2, *(v2si*)_pi32_2);
  mm3 = _mm_and_si64(mm3, *(v2si*)_pi32_2);

  mm2 = _mm_cmpeq_pi32(mm2, _mm_setzero_si64());
  mm3 = _mm_cmpeq_pi32(mm3, _mm_setzero_si64());

  v4sf sign_bit, poly_mask;
  COPY_MM_TO_XMM(mm0, mm1, sign_bit);
  COPY_MM_TO_XMM(mm2, mm3, poly_mask);
  _mm_empty(); /* good-bye mmx */
#endif
  /* The magic pass: "Extended precision modular arithmetic" 
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
  xmm1 = *(v4sf*)_ps_minus_cephes_DP1;
  xmm2 = *(v4sf*)_ps_minus_cephes_DP2;
  xmm3 = *(v4sf*)_ps_minus_cephes_DP3;
  xmm1 = _mm_mul_ps(y, xmm1);
  xmm2 = _mm_mul_ps(y, xmm2);
  xmm3 = _mm_mul_ps(y, xmm3);
  x = _mm_add_ps(x, xmm1);
  x = _mm_add_ps(x, xmm2);
  x = _mm_add_ps(x, xmm3);
  
  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  y = *(v4sf*)_ps_coscof_p0;
  v4sf z = _mm_mul_ps(x,x);

  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(v4sf*)_ps_coscof_p1);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(v4sf*)_ps_coscof_p2);
  y = _mm_mul_ps(y, z);
  y = _mm_mul_ps(y, z);
  v4sf tmp = _mm_mul_ps(z, *(v4sf*)_ps_0p5);
  y = _mm_sub_ps(y, tmp);
  y = _mm_add_ps(y, *(v4sf*)_ps_1);
  
  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

  v4sf y2 = *(v4sf*)_ps_sincof_p0;
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(v4sf*)_ps_sincof_p1);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(v4sf*)_ps_sincof_p2);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_mul_ps(y2, x);
  y2 = _mm_add_ps(y2, x);

  /* select the correct result from the two polynoms */  
  xmm3 = poly_mask;
  y2 = _mm_and_ps(xmm3, y2); //, xmm3);
  y = _mm_andnot_ps(xmm3, y);
  y = _mm_add_ps(y,y2);
  /* update the sign */
  y = _mm_xor_ps(y, sign_bit);

  return y;
}

/* since sin_ps and cos_ps are almost identical, sincos_ps could replace both of them..
   it is almost as fast, and gives you a free cosine with your sine */
inline void sincos_ps(v4sf x, v4sf *s, v4sf *c) {
  v4sf xmm1, xmm2, xmm3 = _mm_setzero_ps(), sign_bit_sin, y;
#ifdef USE_SSE2
  v4si emm0, emm2, emm4;
#else
  v2si mm0, mm1, mm2, mm3, mm4, mm5;
#endif
  sign_bit_sin = x;
  /* take the absolute value */
  x = _mm_and_ps(x, *(v4sf*)_ps_inv_sign_mask);
  /* extract the sign bit (upper one) */
  sign_bit_sin = _mm_and_ps(sign_bit_sin, *(v4sf*)_ps_sign_mask);
  
  /* scale by 4/Pi */
  y = _mm_mul_ps(x, *(v4sf*)_ps_cephes_FOPI);
    
#ifdef USE_SSE2
  /* store the integer part of y in emm2 */
  emm2 = _mm_cvttps_epi32(y);

  /* j=(j+1) & (~1) (see the cephes sources) */
  emm2 = _mm_add_epi32(emm2, *(v4si*)_pi32_1);
  emm2 = _mm_and_si128(emm2, *(v4si*)_pi32_inv1);
  y = _mm_cvtepi32_ps(emm2);

  emm4 = emm2;

  /* get the swap sign flag for the sine */
  emm0 = _mm_and_si128(emm2, *(v4si*)_pi32_4);
  emm0 = _mm_slli_epi32(emm0, 29);
  v4sf swap_sign_bit_sin = _mm_castsi128_ps(emm0);

  /* get the polynom selection mask for the sine*/
  emm2 = _mm_and_si128(emm2, *(v4si*)_pi32_2);
  emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
  v4sf poly_mask = _mm_castsi128_ps(emm2);
#else
  /* store the integer part of y in mm2:mm3 */
  xmm3 = _mm_movehl_ps(xmm3, y);
  mm2 = _mm_cvttps_pi32(y);
  mm3 = _mm_cvttps_pi32(xmm3);

  /* j=(j+1) & (~1) (see the cephes sources) */
  mm2 = _mm_add_pi32(mm2, *(v2si*)_pi32_1);
  mm3 = _mm_add_pi32(mm3, *(v2si*)_pi32_1);
  mm2 = _mm_and_si64(mm2, *(v2si*)_pi32_inv1);
  mm3 = _mm_and_si64(mm3, *(v2si*)_pi32_inv1);

  y = _mm_cvtpi32x2_ps(mm2, mm3);

  mm4 = mm2;
  mm5 = mm3;

  /* get the swap sign flag for the sine */
  mm0 = _mm_and_si64(mm2, *(v2si*)_pi32_4);
  mm1 = _mm_and_si64(mm3, *(v2si*)_pi32_4);
  mm0 = _mm_slli_pi32(mm0, 29);
  mm1 = _mm_slli_pi32(mm1, 29);
  v4sf swap_sign_bit_sin;
  COPY_MM_TO_XMM(mm0, mm1, swap_sign_bit_sin);

  /* get the polynom selection mask for the sine */

  mm2 = _mm_and_si64(mm2, *(v2si*)_pi32_2);
  mm3 = _mm_and_si64(mm3, *(v2si*)_pi32_2);
  mm2 = _mm_cmpeq_pi32(mm2, _mm_setzero_si64());
  mm3 = _mm_cmpeq_pi32(mm3, _mm_setzero_si64());
  v4sf poly_mask;
  COPY_MM_TO_XMM(mm2, mm3, poly_mask);
#endif

  /* The magic pass: "Extended precision modular arithmetic" 
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
  xmm1 = *(v4sf*)_ps_minus_cephes_DP1;
  xmm2 = *(v4sf*)_ps_minus_cephes_DP2;
  xmm3 = *(v4sf*)_ps_minus_cephes_DP3;
  xmm1 = _mm_mul_ps(y, xmm1);
  xmm2 = _mm_mul_ps(y, xmm2);
  xmm3 = _mm_mul_ps(y, xmm3);
  x = _mm_add_ps(x, xmm1);
  x = _mm_add_ps(x, xmm2);
  x = _mm_add_ps(x, xmm3);

#ifdef USE_SSE2
  emm4 = _mm_sub_epi32(emm4, *(v4si*)_pi32_2);
  emm4 = _mm_andnot_si128(emm4, *(v4si*)_pi32_4);
  emm4 = _mm_slli_epi32(emm4, 29);
  v4sf sign_bit_cos = _mm_castsi128_ps(emm4);
#else
  /* get the sign flag for the cosine */
  mm4 = _mm_sub_pi32(mm4, *(v2si*)_pi32_2);
  mm5 = _mm_sub_pi32(mm5, *(v2si*)_pi32_2);
  mm4 = _mm_andnot_si64(mm4, *(v2si*)_pi32_4);
  mm5 = _mm_andnot_si64(mm5, *(v2si*)_pi32_4);
  mm4 = _mm_slli_pi32(mm4, 29);
  mm5 = _mm_slli_pi32(mm5, 29);
  v4sf sign_bit_cos;
  COPY_MM_TO_XMM(mm4, mm5, sign_bit_cos);
  _mm_empty(); /* good-bye mmx */
#endif

  sign_bit_sin = _mm_xor_ps(sign_bit_sin, swap_sign_bit_sin);

  
  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  v4sf z = _mm_mul_ps(x,x);
  y = *(v4sf*)_ps_coscof_p0;

  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(v4sf*)_ps_coscof_p1);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, *(v4sf*)_ps_coscof_p2);
  y = _mm_mul_ps(y, z);
  y = _mm_mul_ps(y, z);
  v4sf tmp = _mm_mul_ps(z, *(v4sf*)_ps_0p5);
  y = _mm_sub_ps(y, tmp);
  y = _mm_add_ps(y, *(v4sf*)_ps_1);
  
  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

  v4sf y2 = *(v4sf*)_ps_sincof_p0;
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(v4sf*)_ps_sincof_p1);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, *(v4sf*)_ps_sincof_p2);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_mul_ps(y2, x);
  y2 = _mm_add_ps(y2, x);

  /* select the correct result from the two polynoms */  
  xmm3 = poly_mask;
  v4sf ysin2 = _mm_and_ps(xmm3, y2);
  v4sf ysin1 = _mm_andnot_ps(xmm3, y);
  y2 = _mm_sub_ps(y2,ysin2);
  y = _mm_sub_ps(y, ysin1);

  xmm1 = _mm_add_ps(ysin1,ysin2);
  xmm2 = _mm_add_ps(y,y2);
 
  /* update the sign */
  *s = _mm_xor_ps(xmm1, sign_bit_sin);
  *c = _mm_xor_ps(xmm2, sign_bit_cos);
}

inline sse_t sseAbs(const sse_t &x) {
	static const sse_t sign_mask = sseSet(REAL_CONST(-0.)); // -0.f = 1 << 31
	return sseAndNot(sign_mask, x);
}

#endif

} // End namespace ONETOOL
