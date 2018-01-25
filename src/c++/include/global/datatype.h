#pragma once
#ifndef __WISARD_DATATYPE_H__
#define __WISARD_DATATYPE_H__

namespace ONETOOL {

/* Basic typedefs */
typedef unsigned short	USHORT_t;
typedef unsigned int	wsUint;
typedef wsUint*			wsUintPtr;
typedef const int		wsIntCst;
typedef const char*		wsStrCst;
typedef const wsUint	wsUintCst;
#define UINT_s			static wsUint

/* System dependent typedefs */
#ifdef __GNUC__
#	define FMT_INT64	"%lld"
#else
#	define FMT_INT64	"%I64d"
#endif

#ifdef _WIN32
typedef __int64 INT64_t;
#else
typedef long long INT64_t;
#endif

/* Floating-point data type will be "double" when it is 1,
 * or be "wsReal" otherwise.
 */
#ifdef USE_DBL
typedef double wsReal; /* Type of floating point data */
#	define REAL_CONST(X_l)	X_l
#	define REAL_FMT			"%lf"
#	ifdef WIN32
#		define REAL_Ut		unsigned INT64_t
#	else
#		define REAL_Ut		unsigned INT64_t
#	endif
#else
typedef float wsReal; /* Type of floating point data */
#	define REAL_CONST(X_l)	X_l##f
#	define REAL_FMT			"%f"
#	define REAL_UINT		unsigned int
#endif

typedef const wsReal	wsRealCst;
typedef wsReal*			wsRealPtr;
typedef	wsReal*			wsVec;
typedef	const wsReal*	wsVecCst;
typedef	const wsReal**	wsMatCst;
typedef	const wsReal**	wsSymCst;
typedef	wsReal**		wsMat;
typedef wsReal**		wsSym;
typedef wsUint**		wsUmat;
typedef wsReal*			wsDiag;
#define SPARSE_t		wsReal*
#define DMAT_t			double**	/* Explicit double matrix */

} // End namespace ONETOOL

#endif
