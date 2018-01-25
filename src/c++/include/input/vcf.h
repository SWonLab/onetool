#ifndef __WISARD_VCF_H__
#define __WISARD_VCF_H__
#pragma once
#include "global/io.h"
#include "input/stream.h"

namespace ONETOOL {

class cFa;
struct xSeq
{
	cFa		*Cp_fa;
	string	S_desc;
	string	S_buf;
	size_t	N_pos;
	void load();
	void unload();
};

typedef vector<xSeq>	vSeq;
typedef	vSeq::iterator	vSeq_it;

class cFa
{
	bool			B_init;
	vSeq			Xv_seq;
	cStrFile		X_f;
public:
	cFa() {
		B_init = false;
	}
	cFa(wsStrCst S_fn, wsStrCst S_desc=NULL) {
		B_init = init(S_fn, S_desc);
	}
	~cFa() {
	}
	cStrFile&	getHandle();
	vSeq&		getSeqs();
	bool		init(wsStrCst S_fn, wsStrCst S_desc=NULL);
};

typedef enum _xVcfLoadState
{
	VCF_LOAD_META,
	VCF_LOAD_HEADER,
	VCF_LOAD_SNP,
	VCF_LOAD_OK
} xVcfLoadState;

typedef struct _xVcfRecord
{
	wsUint	L_shared;
	wsUint	L_indiv;
	int		N_chrom;
	int		N_rlen;
	float	R_qual;
	wsUint	N_info:16;
	wsUint	N_allele:16;
	wsUint	N_sample:24;
	wsUint	N_fmt:8;
} xVcfRecord;

#define BCF_MISSING 0x7F800001

typedef vector<char> vChar;
typedef vector<short> vShort;

class cVcfIO : public cIO
{
	xFileType		X_fType;
	vStr			Sa_meta;
	xVcfLoadState	N_loadState;
	vChar*			Xa_geno;
	vShort*			Xa_dosage;
	wsUint			N_idxGeno;
	size_t			N_dataPos;
	wsUint			N_idxSNP;
	char			B_bin;
	char			N_version;

	// For BCF - requires INFO / CONTIG sequence in VCF/BCF
public:
	mDataIdx		Xm_seqINFO, Xm_seqCONTIG;

	cVcfIO(char *S_fn, xFileType X_type, char B_inDry);
	void		init(char *S_fn);
	wsStrCst	getFormat() { return B_bin ? "io.bcf" : "io.vcf"; }
	int			procVcfLine(char *Sp_buf, wsUint N_line);
	void		_procVcfLine_header(char *Sp_buf, wsUint N_line);
	void		_procVcfLine_meta(char *Sp_buf, wsUint N_line);
	vStr&		getMeta() { return Sa_meta; }
};

} // End namespace ONETOOL

#endif
