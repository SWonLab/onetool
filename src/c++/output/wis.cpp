#include <stdarg.h>
#include "global/common.h"
#include "output/wis.h"
#include "input/stream.h"
#include "utils/util.h"

namespace ONETOOL {

unsigned char cs(unsigned char *Na_buf, wsUint N_sz)
{
	unsigned char N = 0x00;
	wsUint N_med = 0;
#ifdef USE_SSE4
	N_med = getMed(N_sz); 
	sse_t sse_N = sseSet(0.0);
	wsReal *sse_P = (wsReal *)Na_buf;
	for (wsUint i=0 ; i<N_med ; i+=sseJmp)
		sse_N = sseXor(sse_N, *((sse_t *)(sse_P+i)));
#else
	for (wsUint i=N_med ; i<N_sz ; i++)
		N ^= Na_buf[i];
#endif
	return N;
}

CHUNK_LOAD_PROC(MTRX)
{
	xAttrMTRX	*Xp_ai	= (xAttrMTRX *)(Xp_d->Sp_attr);
	wsUint		N_s		= Xp_ai->N_row;
	wsMat		Ra_ret	= sseMatrix(N_s, N_s);

	if (sizeof(wsReal) != sizeof(float)) {
		float		*Ra_buf	= NULL;
		wsAlloc(Ra_buf, float, N_s);
		for (wsUint i=0 ; i<N_s ; i++) {
			/* Read N (float) */
			if (C_stream.read(Ra_buf, N_s) != N_s)
				return NULL;
			for (wsUint j=0 ; j<N_s ; j++)
				Ra_ret[i][j] = (wsReal)Ra_buf[j];
		}
		DEALLOC(Ra_buf);
	} else {
		for (wsUint i=0 ; i<N_s ; i++) {
			/* Read N (float) */
			if (C_stream.read(Ra_ret[i], N_s) != N_s)
				return NULL;
		}
	}

	return Ra_ret;
}

CHUNK_ADD_PROC(WISD)
{
	timeval t;
	gettimeofday(&t, NULL);
	xAttrWISD *Xp_am = (xAttrWISD *)(Xp_d->Sp_attr);
	Xp_am->N_ver	= 1;
	Xp_am->N_time	= t.tv_sec;
}

CHUNK_ADD_PROC(MTRX)
{
	xAttrMTRX *Xp_am = (xAttrMTRX *)(Xp_d->Sp_attr);
	Xp_d->Sp_data	= va_arg(X_args, wsMat);
	Xp_am->N_ver	= 1;
	Xp_am->N_row	= va_arg(X_args, wsUint);
	Xp_am->N_col	= va_arg(X_args, wsUint);
	Xp_am->N_szElem	= va_arg(X_args, wsUint);
	Xp_d->N_szdata	= Xp_am->N_row * Xp_am->N_col * Xp_am->N_szElem;
}

CHUNK_STORE_PROC(MTRX)
{
	xAttrMTRX	*Xp_am	= (xAttrMTRX *)(Xp_d->Sp_attr);
	wsMat		Ra_m	= (wsMat)Xp_d->Sp_data;

	if (sizeof(wsReal) == Xp_am->N_szElem)
		for (wsUint i=0 ; i<Xp_am->N_row ; i++)
			Cp_e->write(Ra_m[i], sizeof(wsReal)*Xp_am->N_col);
	else if (Xp_am->N_szElem == 4) {
		float *Ra_v = NULL;
		wsAlloc(Ra_v, float, Xp_am->N_col);

		for (wsUint i=0 ; i<Xp_am->N_row ; i++) {
			for (wsUint j=0 ; j<Xp_am->N_col ; j++)
				Ra_v[j] = (float)Ra_m[i][j];
			Cp_e->write(Ra_v, sizeof(float)*Xp_am->N_col);
		}
		DEALLOC(Ra_v);
	} else if (Xp_am->N_szElem == 8) {
		double *Ra_v = NULL;
		wsAlloc(Ra_v, double, Xp_am->N_col);

		for (wsUint i=0 ; i<Xp_am->N_row ; i++) {
			for (wsUint j=0 ; j<Xp_am->N_col ; j++)
				Ra_v[j] = (double)Ra_m[i][j];
			Cp_e->write(Ra_v, sizeof(double)*Xp_am->N_col);
		}
		DEALLOC(Ra_v);
	}
}

CHUNK_ADD_PROC(STRV)
{
	xAttrVECT *Xp_av = (xAttrVECT *)(Xp_d->Sp_attr);
	Xp_d->Sp_data	= va_arg(X_args, char *);
	Xp_av->N_ver	= 1;
	Xp_av->N_sz		= va_arg(X_args, wsUint);
}

CHUNK_ADD_PROC(ICOR)
{
	xAttrICOR *Xp_ai = (xAttrICOR *)(Xp_d->Sp_attr);
	Xp_ai->N_ver	= 1;
	Xp_ai->B_header	= va_arg(X_args, wsUint);
	Xp_ai->N_samp	= va_arg(X_args, wsUint);
	Xp_ai->N_marker	= va_arg(X_args, wsUint);
	Xp_ai->X_type	= (xCorType)va_arg(X_args, int);
}

CHUNK_LOAD_PROC(STRV)
{
	xAttrVECT	*Xp_av	= (xAttrVECT *)(Xp_d->Sp_attr);
	wsUint		N_s		= Xp_av->N_sz;
	char**		Sa_ret	= NULL;
	wsAlloc(Sa_ret, char*, N_s);
	wsUint		i=0;

	/* Retrieve */
	char*		S_tmp = NULL;
	wsAlloc(S_tmp, char, Xp_av->N_szData);
	if (C_stream.read(S_tmp, Xp_av->N_szData) != Xp_av->N_szData)
		goto _term;

	/* Get strings */
	for (char *a=S_tmp,*b=NULL ; a ; a=b,i++) {
		getString(&a, &b);
		if (i >= N_s)
			goto _term;
		Sa_ret[i] = a;
	}
	if (i != N_s)
		goto _term;

	return Sa_ret;
_term:
	DEALLOC(Sa_ret);
	DEALLOC(S_tmp);
	return NULL;
}

CHUNK_LOAD_PROC(GENO)
{
	xAttrGENO *Xp_ag = (xAttrGENO *)(Xp_d->Sp_attr);

	/* Allocate memory */
// 	char**	Na_data = NULL;
// 	MULTI_MALLOC(Na_data, char*, Xp_ag->N_samp);

	/* Determine readsize */
	wsUint N_rsz = Xp_ag->B_sampWise ? ((Xp_ag->N_marker+3) >> 2) << 2 :
		((Xp_ag->N_samp+3) >> 2) << 2;
	wsUint N_loop = Xp_ag->B_sampWise ? Xp_ag->N_samp : Xp_ag->N_marker;
	unsigned char *Na_buf = NULL;
	wsAlloc(Na_buf, unsigned char, N_rsz);

	/* Read each block, and validate */
	for (wsUint i=0 ; i<N_loop ; i++) {
		/* Read */
		C_stream.read(Na_buf, N_rsz);
		/* Checksum */
		//unsigned char N_cs = cs(Na_buf, N_rsz);
	}

	return NULL;
}

xChunkDef Xa_chunks[] = {
	{ CT_WISD, "WISD", true , 1, sizeof(xAttrWISD), NULL, NULL, ckadd_WISD, NULL },
	{ CT_ICOR, "ICOR", false, 1, sizeof(xAttrICOR), NULL, NULL, ckadd_ICOR, NULL },
	{ CT_STRV, "STRV", false, 1, sizeof(xAttrVECT), NULL, ckload_STRV, ckadd_STRV, NULL },
	{ CT_MTRX, "MTRX", false, 1, sizeof(xAttrMTRX), NULL, ckload_MTRX, ckadd_MTRX, ckstore_MTRX },
	{ CT_GENO, "GENO", false, 1, sizeof(xAttrGENO), NULL, ckload_GENO, NULL, NULL },
};	

void* load_BINV(cStream *S, size_t N_sz)
{
	char *Sp_ret = new char[N_sz+1];

	/* Read specified bytes from stream */
	if (S->read(Sp_ret, (wsUint)N_sz+1) != (size_t)(N_sz+1))
		return NULL;

	/* Get hash */
	char H = 0;
	for (size_t i=0 ; i<N_sz ; i++)
		H ^= Sp_ret[i];
	/* Hash check */
	if (H != Sp_ret[N_sz])
		return NULL;

	return Sp_ret;
}

xWisData* cWisReader::_read(cStream &S, xWisData *Xp_par,
							xChunkType X_expect/*=CT_NONE*/)
{
	xWisData	*Xp_ret	= new xWisData;
	memset(Xp_ret, 0x00, sizeof(xWisData));

	/* Read 4 bytes to check chunk */
	if (S.read(Xp_ret->S_chunk, 4) != 4)
		return NULL;

	/* If something is expect, get that */
	xChunkDef *Xp_cd = NULL;
	if (X_expect != CT_NONE) {
		Xp_cd = Xa_chunks + X_expect;

		/* Check is this correct chunk */
		if (strcmp(Xp_cd->S_cname, Xp_ret->S_chunk))
			return NULL;
	} else for (wsUint i=0 ; i<len(Xa_chunks, xChunkDef) ; i++) {
		if (!strcmp(Xa_chunks[i].S_cname, Xp_ret->S_chunk)) {
			Xp_cd = Xa_chunks + i;

			break;
		}

	}

	/* If chunk cannot found */
	if (Xp_cd == NULL)
		return NULL;

	/* Now we have chunk, read attribute field */
	Xp_ret->X_type	= Xp_cd->X_type;
	Xp_ret->Sp_attr	= new char[Xp_cd->N_szAttr];
	if (S.read(Xp_ret->Sp_attr, Xp_cd->N_szAttr) != Xp_cd->N_szAttr)
		return NULL;

	xAttrDeft *Xp_attr = (xAttrDeft *)Xp_ret->Sp_attr;
	/* According to the number of chunks have, read */
	if (Xp_attr->N_subChunk) {
		wsAlloc(Xp_ret->Xp_subs, xWisData*, Xp_attr->N_subChunk);
		for (wsUint i=0 ; i<Xp_attr->N_subChunk ; i++) {
			Xp_ret->Xp_subs[i] = _read(S, Xp_ret);
			if (Xp_ret->Xp_subs[i] == NULL)
				return NULL;
		}
		Xp_ret->N_sub = Xp_attr->N_subChunk;
	}

	/* Now read the data */
	if (Xp_cd->H_load)
		Xp_ret->Sp_data = Xp_cd->H_load(S, Xp_par, Xp_ret);
	else {
		Xp_ret->Sp_data = new char[Xp_attr->N_szData];
		if (S.read(Xp_ret->Sp_data, Xp_attr->N_szData) != Xp_attr->N_szData)
			return NULL;
	}

	/* Now OK */
	return Xp_ret;
}

cWisReader::cWisReader(wsStrCst S_path) {
	/* Init stream */
	cStrFile S(S_path);
	Xp_d = _read(S, NULL, CT_WISD);
	if (Xp_d == NULL)
		halt("Failed to load WIS file [%s]", S_path);
}

cWisWriter::cWisWriter()
{
	Xp_d = new xWisData;
	memset(Xp_d, 0x00, sizeof(xWisData));
	strcpy(Xp_d->S_chunk, "WISD");
	Xp_d->X_type = CT_WISD;
	Xp_d->Sp_attr = new xAttrWISD;
	((xAttrWISD *)Xp_d->Sp_attr)->N_ver = 1;
	((xAttrWISD *)Xp_d->Sp_attr)->N_subChunk = 0;
	((xAttrWISD *)Xp_d->Sp_attr)->N_szData = 0;
	timeval v;
	gettimeofday(&v, NULL);
	((xAttrWISD *)Xp_d->Sp_attr)->N_time = v.tv_sec;
}

xWisData* cWisWriter::root()
{
	return Xp_d;
}

void cWisWriter::write(wsStrCst S_ext)
{
	cExporter *C_exp = cExporter::summon(S_ext, 0, ET_BIN);
	_write(C_exp, Xp_d);
	delete C_exp;
}

void cWisWriter::_write(cExporter *Cp_exp, xWisData *Xp_d)
{
	/* Write out FOURCC */
	Cp_exp->write(Xp_d->S_chunk, 4);
	/* Write out attributes */
	Cp_exp->write(Xp_d->Sp_attr, Xa_chunks[Xp_d->X_type].N_szAttr);
	/* Write out subchunks */
	for (wsUint i=0 ; i<Xp_d->N_sub ; i++)
		_write(Cp_exp, Xp_d->Xp_subs[i]);
	/* Write out data */
	if (Xa_chunks[Xp_d->X_type].H_store)
		Xa_chunks[Xp_d->X_type].H_store(Cp_exp, Xp_d);
	else
		Cp_exp->write(Xp_d->Sp_data, Xp_d->N_szdata);
}

xWisData* xWisData::add(xChunkType X_type, ...)
{
	va_list args;
	va_start(args, X_type);

	xWisData *Xp_wd = new xWisData;
	memset(Xp_wd, 0x00, sizeof(xWisData));
	strcpy(Xp_wd->S_chunk, Xa_chunks[X_type].S_cname);
	Xp_wd->Sp_attr	= new unsigned char[Xa_chunks[X_type].N_szAttr];
	memset(Xp_wd->Sp_attr, 0x00, Xa_chunks[X_type].N_szAttr);
	Xp_wd->X_type	= X_type;

	xAttrDeft *Xp_attr = (xAttrDeft *)Sp_attr;
	/* Insert newly to the original vector */
	xWisData **Xp_new = NULL;
	wsAlloc(Xp_new, xWisData*, this->N_sub+1);
	if (this->Xp_subs) {
		memcpy(Xp_new, this->Xp_subs, sizeof(xWisData*)*this->N_sub);
		Xp_new[this->N_sub] = Xp_wd;
		this->N_sub++;
		DEALLOC(this->Xp_subs);
		Xp_subs = Xp_new;
	} else {
		Xp_new[this->N_sub] = Xp_wd;
		this->N_sub++;
		Xp_subs = Xp_new;
	}
	Xp_attr->N_subChunk++;

	/* Perform add command */
	if (Xa_chunks[X_type].H_add == NULL)
		halt("SYSERR : [%s] does not have add function", Xa_chunks[X_type].S_cname);
	Xa_chunks[X_type].H_add(this, Xp_wd, args);

	va_end(args);

	return Xp_wd;
}

#if 0
void exportCorr()
{
	/* Load samples */
	vSampPtr &Xa_samp = Cp_IO->getSamplesV2();

	/* Make FID/IID vector string */
	char *Sp_strv = NULL, *p = NULL;
	MULTI_MALLOC(Sp_strv, 256 * N_samp);
	FOREACH (vSampPtr_it, Xa_samp, i) {
		size_t q;

		sprintf(p, "%s%n", (*i)->S_FID.c_str(), &q);
		p += q+1;
		sprintf(p, "%s%n", (*i)->S_IID.c_str(), &q);
		p += q+1;
	}

	/* Prepare writer */
	cWisWriter X_cor;
	xWisData	*Xp_r	= X_cor.getRoot();
	/* Add ICOR */
	xWisData	*Xp_i	= Xp_r->add(CT_ICOR, N_samp, N_snp);
	/* Add STRV and MTRX */
	Xp_i->add(CT_STRV, Sp_strv, N_samp*2);
	Xp_i->add(CT_MTRX, getFullCorMat(), N_samp, N_samp);

	/* Write out */
	X_cor.write("");

}
#endif

} // End namespace ONETOOL
