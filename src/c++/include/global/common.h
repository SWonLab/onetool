#pragma once
#ifndef __WISARD_COMMON_H__
#define __WISARD_COMMON_H__

#define USE_R		// Force usage of R in analysis
#define USE_GZ		// Force readability for gzipped file
//#define USE_NET		// Force use network interface
#define USE_OPTCHECK	// Force check option elements' integrity
//#define USE_LONGFUNC	// Force to use long function name for logging
#define USE_CORRv2	// Force to use GCTA style correlation
#define USE_SSE		// Force to use SSE-enabled instruction
#define USE_DBL		// Force to use double-precision real
#define USE_SYM
#define USE_ED2
//#define USE_FTP
//#define USE_SSE4
//#define USE_AVX
//#define USE_MEMDBG
//#define USE_CUDA	// Force to use CUDA when CUDA-logic is applicable

#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <list>

#include "program.h"
/*
 * datatype.h should be included at this point since
 * the definition is affected by USE_* directives
 */
#include "datatype.h"

/* Directory letter */
#ifdef _WIN32
#	include <direct.h>
#	define GetCurrentDir _getcwd
#	define WISARD_DIR_LETTER "\\"
#	define WISARD_DIR_CHAR '\\'
#else
#	include <unistd.h>
#	define GetCurrentDir getcwd
#	define WISARD_DIR_LETTER "/"
#	define WISARD_DIR_CHAR '/'
#endif

/*
 * USE_CUDA
 */
#ifdef USE_CUDA
#	include "cuda/cucommon.h"
#endif

/*
 * USE_GZ
 */
#ifdef USE_GZ
 #	include <zlib.h>
#	ifdef _WIN32
#		ifdef _M_X64
#			pragma comment(lib, "zlibstatic.lib")
#		else
#			pragma comment(lib, "zlibstatic.lib")
#		endif
#	endif
#	define myGzClose			gzclose
#	define myGzOpen				gzopen
#	define myGzHandler			gzFile
#	define myGzRewind			gzrewind
#	define myGzGets				gzgets
#	define myGzPrint			gzprintf
#	define myGzPuts(s,p)		gzputs(p,s)
#	define myGzWrite(b,l,c,p)	gzwrite(p,b,(l)*(c))
#else
#	define myGzClose			fclose
#	define myGzOpen				fopen
#	define myGzHandler			FILE *
#	define myGzRewind			::rewind
#	define myGzGets				fgets
#	define myGzPrint			vfprintf
#	define myGzPuts				fputs
#	define myGzWrite			fwrite
#endif

/*
 * USE_R
 */
/* REMOVED 170628
#ifdef USE_R
#	ifdef _WIN32
#		ifdef _M_X64
#			pragma comment(lib, "Rdll_x64.lib")
#		else
#			pragma comment(lib, "Rdll_x86.lib")
#		endif
#	endif
#endif
*/

namespace ONETOOL {

/*
 * Extended system feature
 */
typedef struct _xSysSpec {
	bool	SSSE3;
} xSysSpec;
extern xSysSpec G;
#define OK(a) (G.a)

/*
 * Internal constants
 */

extern const wsReal WISARD_NAN;
extern const wsReal WISARD_INF;
extern const wsReal WISARD_PI;
extern const wsReal WISARD_NA_REAL;
extern const wsReal W0;
extern const wsReal W1;
extern const wsReal W2;

#define WIS_I64NA			0xffffffffffffffff
#define WIS_I32NA			0xffffffff
#define WISARD_NA			-9
#define WISARD_AFFECTED		1
#define WISARD_UNAFFECTED	0
#define NA(v)				((v) != (v))

/*
 * typedefs of globally used vectors and maps
 */
typedef std::vector<char>					vBool;
typedef vBool::iterator						vBool_it;
typedef std::map<std::string, int>			strMap;
typedef strMap::iterator					strMap_it;
typedef std::map<float, int>				realMap;
typedef realMap::iterator					realMap_it;
typedef std::vector<std::string>			vStr;
typedef vStr::iterator						vStr_it;
typedef std::vector<int>					vInt;		/* Explicit 32-bit integer vector */
typedef vInt::iterator						vInt_it;
typedef	std::vector<vInt>					vvInt;
typedef	vvInt::iterator						vvInt_it;
typedef	std::vector<double>					vDbl;		/* Explicit double vector */
typedef	vDbl::iterator						vDbl_it;
typedef std::vector<wsReal>					vReal;		/* Implicit real vector */
typedef vReal::iterator						vReal_it;
typedef std::vector<short *>				vArrWord;
typedef vArrWord::iterator					vArrWord_it;
typedef	std::map<std::string,vStr>			mvStr;
typedef	mvStr::iterator						mvStr_it;
typedef std::map<std::string,wsReal>		mStrReal;
typedef mStrReal::iterator					mStrReal_it;
typedef std::map<std::string, std::string>	mStr;
typedef mStr::iterator						mStr_it;

#ifdef _WIN32
#	define __LONGFUNC__	__FUNCDNAME__
/* Non-GCC doesn't support this, replace to __FUNCTION__ */
#	ifdef USE_LONGFUNC
#		define __FUNC__ __FUNCDNAME__
#	else
#		define __FUNC__ __FUNCTION__
#	endif
#	pragma warning(disable:4996)
#	pragma warning(disable:4503)
#	pragma warning(disable:4706)
#	pragma warning(disable:4819)
#	pragma warning(disable:4800)
/* 4996(deprecated) disable */
#	pragma warning(disable:4996)
/* 4819(Unpresentablechar) disable */
#	pragma warning(disable:4819)
#else
#	define __LONGFUNC__	__PRETTY_FUNCTION__
#	ifdef USE_LONGFUNC
#		define __FUNC__ __PRETTY_FUNCTION__
#	else
#		define __FUNC__ __FUNCTION__
#	endif

#	define __int64 long int
#	define stricmp strcasecmp
#	ifdef __APPLE__
#		define _aligned_malloc(s, a) malloc(s)
#	else
#		define _aligned_malloc(s, a) memalign(a, s)
#	endif
#endif

/* Arithmetic preprocessors */
#define SQR(a)	((a)*(a))
#define CUBE(a)	((a)*(a)*(a))

/* Get the # of entries */
#define len(v,t)	(sizeof(v)/sizeof(t))

} // End namespace ONETOOL

#endif
