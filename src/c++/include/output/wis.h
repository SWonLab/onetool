#pragma once
#include "utils/util.h"

namespace ONETOOL {

typedef enum _xCorType {
	ICOR_UNDEF,
	ICOR_BNFNDERS,
	ICOR_BNALL,
	ICOR_BNXFNDERS,
	ICOR_BNXALL,
	ICOR_IBS
} xCorType;

/*

FOURCC

WISD	Main chunk


 */

typedef struct _xAttrDeft {
	wsUint	N_ver;
	wsUint	N_subChunk;
	wsUint	N_szData;
} xAttrDeft;

typedef struct _xAttrWISD : xAttrDeft {
	long	N_time;		/* File created time */
} xAttrWISD;

typedef struct _xAttrICOR : xAttrDeft {
	wsUint		B_header;
	wsUint		N_samp, N_marker;
	xCorType	X_type;
} xAttrICOR;

typedef struct _xAttrVECT : xAttrDeft {
	wsUint		N_sz;
} xAttrVECT;

typedef struct _xAttrMTRX : xAttrDeft {
	wsUint		N_row, N_col, N_szElem;
} xAttrMTRX;

typedef struct _xAttrGENO : xAttrDeft {
	wsUint		N_samp, N_marker;
	char		B_sampWise;
} xAttrGENO;

typedef enum _xChunkType {
	CT_WISD,
	CT_ICOR,
	CT_STRV,
	CT_MTRX,
	CT_GENO,
	CT_NONE=-1,
} xChunkType;

typedef struct _xWisData
{
	xChunkType	X_type;
	char		S_chunk[5];
	void*		Sp_attr;
	wsUint		N_szdata;
	void*		Sp_data;
	wsUint		N_sub;
	struct _xWisData **Xp_subs;

	struct _xWisData*	add(xChunkType X_type, ...);
} xWisData;

#define	CHUNK_LOAD_PROC(a)	void* ckload_##a(cStream &C_stream, xWisData *Xp_par, xWisData *Xp_d)
#define	CHUNK_ADD_PROC(a)	void ckadd_##a(xWisData *Xp_par, xWisData *Xp_d, va_list X_args)
#define	CHUNK_STORE_PROC(a)	void ckstore_##a(cExporter *Cp_e, xWisData *Xp_d)

class cStream;
class cExporter;
typedef void* (*hChunkLoadProc)(cStream &C_stream, xWisData *Xp_par, xWisData *Xp_ret);
typedef void (*hChunkAddProc)(xWisData *Xp_par, xWisData *Xp_d, va_list X_args);
typedef void (*hChunkStoreProc)(cExporter *Cp_e, xWisData *Xp_d);

typedef struct _xChunkDef {
	xChunkType	X_type;
	wsStrCst		S_cname;
	bool		B_recursivable;
	wsUint		N_ver;
	wsUint		N_szAttr;
	wsStrCst		S_parent;
	hChunkLoadProc	H_load;
	hChunkAddProc	H_add;
	hChunkStoreProc	H_store;
} xChunkDef;

class cWisReader
{
	xWisData *Xp_d;
	xWisData* _read(cStream &S, xWisData *Xp_par, xChunkType X_expect=CT_NONE);
public:
	cWisReader(wsStrCst S_path);
	xWisData* fetch(wsStrCst S_path) {
		char *p = strdup(S_path);

		/* Divide by / */
		char		*S = strtok(p, "/");
		xWisData	*X = Xp_d;
		while (S) {
			xWisData *oldX = X;
			for (wsUint i=0 ; i<X->N_sub ; i++)
				if (!stricmp(X->Xp_subs[i]->S_chunk, S)) {
					X = X->Xp_subs[i];
					break;
				}
			if (oldX == X)
				return NULL;
			S = strtok(NULL, "/");
		}

		return X;
	}
};

class cWisWriter
{
	xWisData *Xp_d;
	void _write(cExporter *Cp_exp, xWisData *Xp_d);
public:
	cWisWriter();
	xWisData*	root();
	void		write(wsStrCst S_ext);
};

} // End namespace ONETOOL
