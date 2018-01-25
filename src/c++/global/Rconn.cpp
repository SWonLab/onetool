#ifdef _WIN32
#	pragma warning(disable:4068)
#	pragma warning(disable:4141)
#	include <windows.h>
#endif
#include <iostream>
#include <stdarg.h>
#include "global/common.h"
#include "utils/util.h"
#include "global/Rconn.h"
#include "global/option.h"
#include "utils/data.h"

#ifdef USE_R
/* 140509 To avoid the confliction between wingdi and R */
#       ifdef ERROR
#       undef ERROR 
#       endif
#	ifndef _WIN32
#		include <dlfcn.h>
#	endif
        #include <R.h>
        #include <Rdefines.h>
        #include <Rembedded.h>
        #include <R_ext/Parse.h>
#endif

extern char S_exePath[MAX_PATH];

namespace ONETOOL {

#ifdef USE_R

	typedef SEXP(*fn_allocVector)	(SEXPTYPE, R_xlen_t);
	typedef SEXP(*fn_ParseVector)	(SEXP, int, ParseStatus*, SEXP file);
	typedef SEXP(*fn_protect)		(SEXP);
	typedef int(*fn_initEmbeddedR)	(int, char**);
	typedef void(*fn_endEmbeddedR)	(size_t);
	typedef SEXP(*fn_eval)			(SEXP, SEXP);
	typedef SEXP(*fn_tryEval)		(SEXP, SEXP, int *);
	typedef SEXP(*fn_ScalarInteger)	(int);
	typedef SEXP(*fn_install)		(const char *);
	typedef SEXP(*fn_lang2)			(SEXP, SEXP);
	typedef SEXP(*fn_allocMatrix)	(SEXPTYPE, int, int);
	typedef void(*fn_defineVar)		(SEXP, SEXP, SEXP);
	typedef void(*fn_unprotect)		(int);
	typedef SEXP(*fn_mkChar)		(const char *);
	typedef SEXP(*fn_findVar)		(SEXP, SEXP);
	typedef R_len_t(*fn_length)		(SEXP);
	typedef SEXP(*fn_STRING_ELT)	(SEXP x, R_xlen_t i);
	typedef SEXP(*fn_VECTOR_ELT)	(SEXP x, R_xlen_t i);
	typedef	void(*fn_SET_STRING_ELT)(SEXP x, R_xlen_t i, SEXP v);
	typedef	SEXP(*fn_SET_VECTOR_ELT)(SEXP x, R_xlen_t i, SEXP v);
	typedef	SEXP*(*fn_STRING_PTR)	(SEXP x);
	typedef	int*(*fn_INTEGER)		(SEXP x);
	typedef	double*(*fn_REAL)		(SEXP x);


SEXP*				_R_GlobalEnv		= NULL;
SEXP*				_R_NilValue			= NULL;
SEXP*				_R_UnboundValue		= NULL;
fn_allocVector		_R_allocVector		= NULL;
fn_ParseVector		_R_ParseVector		= NULL;
fn_protect			_R_protect			= NULL;
fn_initEmbeddedR	_R_initEmbeddedR	= NULL;
fn_endEmbeddedR		_R_endEmbeddedR		= NULL;
fn_eval				_R_eval				= NULL;
fn_tryEval			_R_tryEval			= NULL;
fn_ScalarInteger	_R_ScalarInteger	= NULL;
fn_install			_R_install			= NULL;
fn_lang2			_R_lang2			= NULL;
fn_allocMatrix		_R_allocMatrix		= NULL;
fn_defineVar		_R_defineVar		= NULL;
fn_unprotect		_R_unprotect		= NULL;
fn_mkChar			_R_mkChar			= NULL;
fn_findVar			_R_findVar			= NULL;
fn_length			_R_length			= NULL;
fn_STRING_ELT		_R_STRING_ELT		= NULL;
fn_VECTOR_ELT		_R_VECTOR_ELT		= NULL;
fn_SET_STRING_ELT	_R_SET_STRING_ELT	= NULL;
fn_SET_VECTOR_ELT	_R_SET_VECTOR_ELT	= NULL;
fn_STRING_PTR		_R_STRING_PTR		= NULL;
fn_INTEGER			_R_INTEGER			= NULL;
fn_REAL				_R_REAL				= NULL;


cRconn X_R;
cRconn& R() { return X_R; }

cRconn::cRconn() {
	B_init = 0;
}

#ifdef _WIN32
#	ifdef _M_X64
#		define R_BIN_PATH	WISARD_DIR_LETTER "bin" WISARD_DIR_LETTER "x64"
#	else
#		define R_BIN_PATH	WISARD_DIR_LETTER "bin" WISARD_DIR_LETTER "i386"
#	endif
#else
#		define R_BIN_PATH	"/bin/exec"
#endif

	char* _unixPath(char *Sp_path)
	{
#ifdef _WIN32
		char *Sp_ret = NULL, *q = NULL;
		wsAlloc(Sp_ret, char, strlen(Sp_path)+1);
		q = Sp_ret;
		for (char *p=Sp_path ; *p ; p++) {
			if (*p == '\\') {
				*(q++) = '/';
			} else *(q++) = *p;
		}
		*q = '\0';
		return Sp_ret;
#else
		return strdup(Sp_path);
#endif
	}

void _findRinstallation(char *Sp_Rpath, char &B_Rfound)
{
#ifdef _WIN32
	/* Find 'program files' */
	char	S_pf[MAX_PATH]	= { 0, };
	wsStrCst	Sp_pf			= getenv("ProgramW6432");
	if (!Sp_pf)
		Sp_pf				= getenv("ProgramFiles");
	if (Sp_pf && Sp_pf[0])
		sprintf(S_pf, "%s" WISARD_DIR_LETTER "R" WISARD_DIR_LETTER "*", Sp_pf);
	else
		strcpy(S_pf, "C:" WISARD_DIR_LETTER "Program Files" WISARD_DIR_LETTER "R" WISARD_DIR_LETTER "*");

	/* Find R installation */
	WIN32_FIND_DATA	X_fEnt		= { 0, };
	HANDLE			H_findFile = FindFirstFile(S_pf, &X_fEnt);
	if (H_findFile) {
		do {
			if (X_fEnt.cFileName[0] == '.' ||
				memcmp(X_fEnt.cFileName, "R-", 2)) continue;
			LOG("R installation [%s] found\n", X_fEnt.cFileName);
			/* Remove * */
			char* S_pff = new char[512];
			strcpy(S_pff, S_pf);
			S_pff[strlen(S_pf) - 2] = '\0';
			/* Copy with URL */
			sprintf(Sp_Rpath, "%s" WISARD_DIR_LETTER "%s", S_pff,
				X_fEnt.cFileName);
			delete[] S_pff;
			B_Rfound = 1;
		} while (FindNextFile(H_findFile, &X_fEnt));
		FindClose(H_findFile);
	}
#else
	/* Try whereis */
	FILE* H_whereis = popen("whereis R", "r");
	if (!H_whereis) return;
	char* S_buf = new char[65536];
	while (!feof(H_whereis)) {
		S_buf[0] = '\0';
		fgets(S_buf, 65536, H_whereis);
		lverbose("R whereis [%s]\n", S_buf);
		char *t = strtok(S_buf, " ");
		if (!t || strcmp("R:", t)) break;
		char *t2 = strtok(NULL, " ");
		if (!t2 || !t2[0]) break;
		while (t2 && !B_Rfound) {
			lverbose("Check R directory [%s]\n", t2);
			FILE *fp = fopen(t2, "r");
			while (fgets(S_buf, 65536, fp)) {
				if (strlen(S_buf) > 11 && !memcmp(S_buf, "R_HOME_DIR=", 11)) {
					for (char* t=S_buf+strlen(S_buf)-1 ; (*t == '\r' || *t == ' ' || *t == '\t' || *t == '\n') && t>=S_buf ; t--)
						*t = '\0';
					strcpy(Sp_Rpath, S_buf+11);
					LOG("R installation [%s] found\n", S_buf+11);
					B_Rfound = 1;
					break;
				}
			}
			t2 = strtok(NULL, " ");
			fclose(fp);
		}
	}
	fclose(H_whereis);
	delete [] S_buf;
#endif
}

char S_Rpath[512];
void cRconn::_setRpath(char &B_init)
{
	char *S_path	= getenv("PATH");
	/* Set R path variable */
	if (IS_ASSIGNED(rpath))
		strcpy(S_Rpath, OPT_STRING(rpath));
	else S_Rpath[0] = '\0';
	char *S_npath	= NULL;

	if (S_Rpath[0]) {
		wsAlloc(S_npath, char, strlen(S_path) + strlen(S_Rpath) + 16);
#ifndef _WIN32
		sprintf(S_npath, "R_HOME=%s", S_Rpath);
		setenv("R_HOME", S_Rpath, 1);
		LOG("R_HOME set to [%s]\n", getenv("R_HOME"));
		sprintf(S_npath, "%s:%s" R_BIN_PATH, S_path, S_Rpath);
#else
		sprintf(S_npath, "%s;%s" R_BIN_PATH, S_path, S_Rpath);
#endif
	} else {
		/* Find R path */
		char B_Rfound = 0;
		wsAlloc(S_npath, char, strlen(S_path) + 256 + 16);
		_findRinstallation(S_Rpath, B_Rfound);

		/* Finally, set PATH */
		if (B_Rfound) {
			LOG("Set R_HOME [%s]\n", S_Rpath);
#ifdef _WIN32
			char tmp[512];
			sprintf(tmp, "R_HOME=%s", S_Rpath);
			putenv(tmp);
			sprintf(S_npath, "PATH=%s;%s" R_BIN_PATH, S_path, S_Rpath);
#else
			setenv("R_HOME", S_Rpath, 1);
			sprintf(S_npath, "%s:%s" R_BIN_PATH, S_path, S_Rpath);
#endif
		} else {
			LOGwarn("Failed to find R, R compatibility disabled\n");
			B_init = 0;
		}
	}

	/* Set path only if it is still ok to initiate R compatibility */
	if (B_init)
#ifdef _WIN32
		putenv(S_npath);
#else
		setenv("PATH", S_npath, 1);
#endif

	DEALLOC(S_npath);
}

cRconn::~cRconn() {
	/* Cleanup embedded-R environment */
	if (B_init) {
		Rparse("options(warn=-1);q(save='no')");
		_R_endEmbeddedR(0);
	}
}

// #define R_VAL(n,t)		t _ ## n
// #define R_DEF(r,n,t)	r (* _ ## n)t
// #define R_IMP(r,n,t)	_ ## n = (r (*)t)GetProcAddress(H_dll, #n)
// #define R_VIM(n,t)		_ ## n = (t)GetProcAddress(H_dll, #n)

SEXP callrnorm(){ 
	int err = 0;
	SEXP call = _R_protect( _R_lang2( _R_install( "rnorm"), _R_ScalarInteger(10) ) ) ; 
	SEXP res = _R_protect( _R_tryEval( call, *_R_GlobalEnv, &err ) ) ;
	if (err) LOGwarn("R command was failed");
	_R_unprotect(2);
	return res;
} 

SEXP cRconn::Rparse(const char *S_fmt, ...)
{
	/* Set va */
	int		err = 0;
	char	S_buf[4096];
	va_list	H_varList;
	va_start(H_varList, S_fmt);
	vsprintf(S_buf, S_fmt, H_varList);

	SEXP cmdSexp, cmdexpr, ans = *_R_NilValue;
	ParseStatus status;
	_R_protect(cmdSexp = _R_allocVector(STRSXP, 1));
	_R_SET_STRING_ELT(cmdSexp, 0, _R_mkChar(S_buf));
	cmdexpr = _R_protect(_R_ParseVector(cmdSexp, -1, &status, *_R_NilValue));
	if (status != PARSE_OK) {
		_R_unprotect(2);
		halt("Invalid call %s", S_buf);
	}
	/* Loop is needed here as EXPSEXP will be of length > 1 */
	for(int i=0 ; i<_R_length(cmdexpr) ; i++) {
		ans = _R_tryEval(_R_VECTOR_ELT(cmdexpr, i), *_R_GlobalEnv, &err);
		if (err) {
			LOGwarn("R command was failed, check the above error message!\n");
			break;
		}
	}
	_R_unprotect(2);

	va_end(H_varList);
	return ans;
}

SEXP cRconn::Reval(const char *S_buf)
{
	/* Set va */
	int err = 0;
	SEXP cmdSexp, cmdexpr, ans = *_R_NilValue;
	ParseStatus status;
	_R_protect(cmdSexp = _R_allocVector(STRSXP, 1));
	_R_SET_STRING_ELT(cmdSexp, 0, _R_mkChar(S_buf));
	cmdexpr = _R_protect(_R_ParseVector(cmdSexp, -1, &status, *_R_NilValue));
	if (status != PARSE_OK) {
		_R_unprotect(2);
		halt("Invalid call %s", S_buf);
	}
	/* Loop is needed here as EXPSEXP will be of length > 1 */
	for(int i=0 ; i<_R_length(cmdexpr) ; i++) {
		ans = _R_tryEval(_R_VECTOR_ELT(cmdexpr, i), *_R_GlobalEnv, &err);
		if (err) {
			LOGwarn("R command was failed, check the above error message!\n");
			break;
		}
	}
	_R_unprotect(2);

	return ans;
}

void cRconn::dim(wsStrCst S_varName, wsUint *Np_R, wsUint *Np_C)
{
	char S[1024];
	sprintf(S, "dim(%s)", S_varName);
	SEXP res = Rparse(S);
	int	*Rp_Rvar = _R_INTEGER(res);
	*Np_R = (wsUint)Rp_Rvar[0];
	*Np_C = (wsUint)Rp_Rvar[1];
}

void cRconn::sendOption(cOption &C_opt) {
	char *S_buf = new char[65536], *S_ptr = S_buf;
	S_ptr += sprintf(S_ptr, "OPTS <- list(\n");

	/* Get assigned options */
	vOptPtr& X_ops = C_opt.getAsgnOpts();
	FOREACH (vOptPtr_it, X_ops, i) {
		switch ((*i)->E_type) {
		case OT_ONOFF:	S_ptr += sprintf(S_ptr, "'%s'=%s,\n", (*i)->S_longName+2,
							(*i)->N_intVal?"TRUE":"FALSE"); break;
		case OT_NUMBER:	S_ptr += sprintf(S_ptr, "'%s'=%d,\n", (*i)->S_longName+2,
							(*i)->N_intVal); break;
		case OT_RANGE:	S_ptr += sprintf(S_ptr, "'%s'=c(%g,%g),\n", (*i)->S_longName+2,
							(*i)->X_rngVal.R_s, (*i)->X_rngVal.R_e); break;
		case OT_STRING:	S_ptr += sprintf(S_ptr, "'%s'='%s',\n", (*i)->S_longName+2,
							(*i)->S_strVal); break;
		case OT_REAL:	S_ptr += sprintf(S_ptr, "'%s'=%g,\n", (*i)->S_longName+2,
							(*i)->R_realVal); break;
		case OT_EXPR:
		default: break;
		}
	}
	*(S_ptr-2) = '\n';
	sprintf(S_ptr, ")");
	Rparse(S_buf);
}

void cRconn::sendVector(const char *S_varName, wsReal *Rp_data, wsUint N_sz) {
	SEXP	X_Rvar;

	_R_protect(X_Rvar = _R_allocVector(REALSXP, N_sz));
	double	*Rp_Rvar = _R_REAL(X_Rvar);
	for (wsUint i=0 ; i<N_sz ; i++)
		Rp_Rvar[i] = Rp_data[i];
	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
	_R_defineVar(_R_install(S_varName), X_Rvar, *_R_GlobalEnv);
}

void cRconn::sendVector(const char *S_varName, int* Np_data, wsUint N_sz) {
	SEXP	X_Rvar;

	_R_protect(X_Rvar = _R_allocVector(INTSXP, N_sz));
	int *Np_Rvar = _R_INTEGER(X_Rvar);
	memcpy(Np_Rvar, Np_data, sizeof(int)*N_sz);
	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
	_R_defineVar(_R_install(S_varName), X_Rvar, *_R_GlobalEnv);
}

void cRconn::sendVector(const char *S_varName, char *Np_data, wsUint N_sz, char N_varNA/*=-9*/) {
	SEXP	X_Rvar;

	_R_protect(X_Rvar = _R_allocVector(REALSXP, N_sz));
	double	*Rp_Rvar = _R_REAL(X_Rvar);
	for (wsUint i=0 ; i<N_sz ; i++)
		Rp_Rvar[i] = Np_data[i] == N_varNA ? WISARD_NAN : (wsReal)Np_data[i];
	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
	_R_defineVar(_R_install(S_varName), X_Rvar, *_R_GlobalEnv);
	_R_unprotect(1);
}

void cRconn::sendVector(const char *S_varName, char **Sp_strs, wsUint N_sz) {
	SEXP	X_Rvar;

	_R_protect(X_Rvar = _R_allocVector(STRSXP, N_sz));
	SEXP	*Xp_Rvar = _R_STRING_PTR(X_Rvar);
	for (wsUint i=0 ; i<N_sz ; i++)
		_R_SET_STRING_ELT(Xp_Rvar[i], 0, _R_mkChar(Sp_strs[i]));
	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
	_R_defineVar(_R_install(S_varName), X_Rvar, *_R_GlobalEnv);
	_R_unprotect(1);
}

void cRconn::sendVector(const char *S_varName, vStr &Xv_strs) {
	SEXP	X_Rvar;
	wsUint	N_sz = (wsUint)Xv_strs.size();
	if (N_sz == 0) return;

	_R_protect(X_Rvar = _R_allocVector(STRSXP, N_sz));
	SEXP	*Xp_Rvar = _R_STRING_PTR(X_Rvar);
	for (wsUint i=0 ; i<N_sz ; i++)
		_R_SET_STRING_ELT(Xp_Rvar[i], 0, _R_mkChar(Xv_strs[i].c_str()));
	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
	_R_defineVar(_R_install(S_varName), X_Rvar, *_R_GlobalEnv);
	_R_unprotect(1);
}

void cRconn::sendMatrix(const char *S_varName, wsReal **Rp_data, wsUint N_row, wsUint N_col) {
	SEXP	X_Rvar;

	_R_protect(X_Rvar = _R_allocMatrix(REALSXP, N_row, N_col));
	double	*Rp_Rvar = _R_REAL(X_Rvar);
	for (wsUint i=0,k=0 ; i<N_row ; k=++i) {
		for (wsUint j=0 ; j<N_col ; j++,k+=N_row) {
			Rp_Rvar[k] = Rp_data[i][j];
//			printf("Allocate %d,%d to %d\n", i, j, k);
		}
	}
	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
	SEXP X_var = _R_protect(_R_install(S_varName));
	_R_defineVar(X_var, X_Rvar, *_R_GlobalEnv);
	_R_unprotect(2);
}

void cRconn::recvVector(const char *S_varName, wsReal *Rp_data, wsUint N_sz) {
	SEXP	X_Rvar;

	_R_protect(X_Rvar = _R_findVar(_R_install(S_varName), *_R_GlobalEnv));
	double	*Rp_Rvar = _R_REAL(X_Rvar);
	for (wsUint i=0 ; i<N_sz ; i++)
		Rp_data[i] = (wsReal)Rp_Rvar[i];
	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
	_R_unprotect(1);
}

void cRconn::recvMatrix(const char *S_varName, wsReal **Rp_data, wsUint N_row, wsUint N_col) {
	SEXP	X_Rvar;

	_R_protect(X_Rvar = _R_findVar(_R_install(S_varName), *_R_GlobalEnv));
	double	*Rp_Rvar = _R_REAL(X_Rvar);
	for (wsUint j=0,k=0 ; j<N_col ; j++)
		for (wsUint i=0 ; i<N_row ; i++,k++)
			Rp_data[i][j] = (wsReal)Rp_Rvar[k];
	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
	_R_unprotect(1);
}

#ifdef _WIN32
typedef HMODULE wsDllHandle;
#	define wsDllLoad(a) LoadLibrary(a)
#	define wsDllGet(h,a) GetProcAddress(h,a)
#	define R_DLL "R.dll"
#else
typedef void* wsDllHandle;
#	define wsDllLoad(a) dlopen(a, RTLD_NOW)
#	define wsDllGet(h,a) dlsym(h,a)
#	define R_DLL "libR.so"
#endif

class cDynLibLoader
{
	wsDllHandle	H_dll;
	string		S_path;
public:
	cDynLibLoader(const char* Sp_path) {
		S_path = Sp_path;
		LOG("Load dynamic library [%s]\n", Sp_path);
		H_dll = wsDllLoad(Sp_path);
	}
	char sane() {
		return H_dll != NULL;
	}
	char tryFrom(const char* Sp_path) {
		S_path = Sp_path;
		H_dll = wsDllLoad(Sp_path);
		return sane();
	}
	void* get(const char* Sp_symbol) {
		return wsDllGet(H_dll, Sp_symbol);
	}
};

void cRconn::init() {
	B_init = 1;
	const char *tmp;
	// FIXME:  if per-session temp directory is used (as R does) then return
	tmp = getenv("TMPDIR");
	if (tmp == NULL) {
		tmp = getenv("TMP");
		if (tmp == NULL) {
			tmp = getenv("TEMP");
			if (tmp == NULL)
				tmp = "/tmp";
			}
	}
#ifndef _WIN32
	if (setenv("R_SESSION_TMPDIR",tmp,1) != 0)
		halt("Could not set / replace R_SESSION_TMPDIR to %s", tmp);
#endif
	/* Set PATH automatically */
	_setRpath(B_init);

	/* Initialize R compatibility only if it is allowed */
	if (B_init) {
		const char *R_argv[] = { "rconn", "--no-readline", "--no-save", "--silent" };
		int R_argc = (sizeof(R_argv) ) / sizeof(R_argv[0]);

		cDynLibLoader C_libR(R_DLL);
		if (!C_libR.sane()) {
#if 0
			char S_dirSO[512];
			sprintf(S_dirSO, "%s/lib/" R_DLL, S_Rpath);
			if (!C_libR.tryFrom(S_dirSO)) {
				LOGwarn("Failed to load R dynamic library, is it exists in the directories of LD_LIBRARY_PATH?\n");
#else
			char S_dirSO[512];
			sprintf(S_dirSO, "%s/lib/", S_Rpath);
			LOGwarn("Failed to load R dynamic library, is it exists in the directories of LD_LIBRARY_PATH?\n       Possible solution: Execute 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s' and re-run.\n", S_dirSO);
#endif
				B_init = 0;
				return;
#if 0
			}
			sprintf(S_dirSO, "%s/lib/", S_Rpath);
			char *Sp_curLDLP = getenv("LD_LIBRARY_PATH");
			char S_newLDLP[2048];
			sprintf(S_newLDLP, "%s:%s", Sp_curLDLP, S_dirSO);
			setenv("LD_LIBRARY_PATH", S_dirSO, 1);
			setenv("LIBPATH", S_dirSO, 1);
#endif
		}
		_R_GlobalEnv		= (SEXP *)C_libR.get("R_GlobalEnv");
		_R_NilValue			= (SEXP *)C_libR.get("R_NilValue");
		_R_UnboundValue		= (SEXP *)C_libR.get("R_UnboundValue");
		_R_allocVector		= (fn_allocVector)C_libR.get("Rf_allocVector");
		_R_ParseVector		= (fn_ParseVector)C_libR.get("R_ParseVector");
		_R_protect			= (fn_protect)C_libR.get("Rf_protect");
		_R_initEmbeddedR	= (fn_initEmbeddedR)C_libR.get("Rf_initEmbeddedR");
		_R_endEmbeddedR		= (fn_endEmbeddedR)C_libR.get("Rf_endEmbeddedR");
		_R_eval				= (fn_eval)C_libR.get("Rf_eval");
		_R_tryEval			= (fn_tryEval)C_libR.get("R_tryEval");
		_R_ScalarInteger	= (fn_ScalarInteger)C_libR.get("Rf_ScalarInteger");
		_R_install			= (fn_install)C_libR.get("Rf_install");
		_R_lang2			= (fn_lang2)C_libR.get("Rf_lang2");
		_R_allocMatrix		= (fn_allocMatrix)C_libR.get("Rf_allocMatrix");
		_R_defineVar		= (fn_defineVar)C_libR.get("Rf_defineVar");
		_R_unprotect		= (fn_unprotect)C_libR.get("Rf_unprotect");
		_R_mkChar			= (fn_mkChar)C_libR.get("Rf_mkChar");
		_R_findVar			= (fn_findVar)C_libR.get("Rf_findVar");
		_R_length			= (fn_length)C_libR.get("Rf_length");
		_R_STRING_ELT		= (fn_STRING_ELT)C_libR.get("STRING_ELT");
		_R_VECTOR_ELT		= (fn_VECTOR_ELT)C_libR.get("VECTOR_ELT");
		_R_SET_STRING_ELT	= (fn_SET_STRING_ELT)C_libR.get("SET_STRING_ELT");
		_R_SET_VECTOR_ELT	= (fn_SET_VECTOR_ELT)C_libR.get("SET_VECTOR_ELT");
		_R_STRING_PTR		= (fn_STRING_PTR)C_libR.get("STRING_PTR");


		/* Always returns 1 */
		LOG("Trying to initialize R (with %d args)...\n", R_argc);

		_R_initEmbeddedR(R_argc, (char **)R_argv);

		/* Loading base */
		if (S_exePath[0] == '\0' && getItselfPath(S_exePath) == 0)
			halt("SYSERR: Failed to get the path of " TOOLSET_NAME);
		lverbose("WISARD executable path [%s]\n", S_exePath);
		char *S_RstyleBase = _unixPath(S_exePath);
		lverbose("Rbase path [%s]\n", S_RstyleBase);

		/* Check file existance */
		char S_baseRpath[MAX_PATH] = { 0, };
		if (S_exePath[strlen(S_exePath)-1] == WISARD_DIR_CHAR)
			sprintf(S_baseRpath, "%s" "onetool.R", S_exePath);
		else
			sprintf(S_baseRpath, "%s" WISARD_DIR_LETTER "onetool.R", S_exePath);
		FILE *f = fopen(S_baseRpath, "r");
		if (f != NULL) {
			fclose(f);
			Rparse("source('%s/onetool.R')", S_RstyleBase);
			DEALLOC(S_RstyleBase);
			LOGnote("Base R script loaded\n");
		} else {
			LOGwarn("Base R script [%s] not found\n", S_baseRpath);
			SEXP load = _R_protect( Rparse("attr(try(library(kinship2), T), 'class')") );
			if (load == *_R_NilValue)
				LOGnote("[kinship2] R package successfully loaded\n");
			else
				LOGwarn("Failed to load [kinship2] R package!\n");
			_R_unprotect(1);
		}
			
		LOG("R initialized\n");
	}
}

cRvar& cRconn::operator[]( const char *S_varName )
{
	cRvar *C_ret;

	SEXP X_ret = _R_findVar(_R_install(S_varName), *_R_GlobalEnv);
	if (X_ret == *_R_UnboundValue)
		C_ret = new cRvar(S_varName, *_R_NilValue);
	else
		C_ret = new cRvar(S_varName, X_ret);

	return *C_ret;
}

// void cRconn::operator[](const char *S_varName)
// {
// 	SEXP			X_Rvar;
// 	MYuint	N_row = C_mat.getRow();
// 	MYuint	N_col = C_mat.getCol();
// 	wsReal			**Rp_data = C_mat.getData();
// 
// 	_R_protect(X_Rvar = _R_allocMatrix(REALSXP, N_row, N_col));
// 	double	*Rp_Rvar = _R_REAL(X_Rvar);
// 	for (MYuint i=0,k=0 ; i<N_row ; k=++i) {
// 		for (MYuint j=0 ; j<N_col ; j++,k+=N_row) {
// 			Rp_Rvar[k] = Rp_data[i][j];
// 			//			printf("Allocate %d,%d to %d\n", i, j, k);
// 		}
// 	}
// 	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
// 	SEXP X_var = _R_protect(_R_install(S_varName));
// 	_R_defineVar(X_var, X_Rvar, *_R_GlobalEnv);
// }

void cRvar::operator=(cMatrix &C_mat)
{
	SEXP			X_Rvar;
	wsUint	N_row = C_mat.row();
	wsUint	N_col = C_mat.col();
	wsReal			**Rp_data = C_mat.get();

	_R_protect(X_Rvar = _R_allocMatrix(REALSXP, N_row, N_col));
	double	*Rp_Rvar = _R_REAL(X_Rvar);
	for (wsUint i=0,k=0 ; i<N_row ; k=++i) {
		for (wsUint j=0 ; j<N_col ; j++,k+=N_row) {
			Rp_Rvar[k] = Rp_data[i][j];
			//			printf("Allocate %d,%d to %d\n", i, j, k);
		}
	}
	//	memcpy(Rp_Rvar, Rp_data, sizeof(double)*N_sz);
	X_var = _R_protect(_R_install(S_varName));
	_R_defineVar(X_var, X_Rvar, *_R_GlobalEnv);
}

cRvar::~cRvar()
{
	DEALLOC(S_varName);
}

void cRvar::set(const char *S_newVarName, SEXP X_exp)
{
	S_varName = strdup(S_newVarName);
	X_var = X_exp;
}

#else

void cRconn::init() {
	/* Do nothing */
}

cRconn::~cRconn() {
}

#endif

} // End namespace ONETOOL
