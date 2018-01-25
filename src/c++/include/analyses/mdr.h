#pragma once
#ifndef __WISARD_MDR_H__
#define __WISARD_MDR_H__
#include "global/analysis.h"

namespace ONETOOL {

class cResult;

typedef struct _xHmdrEntry {
	vInt	Nv_comb;
	float	R_mea;
	vInt	Nv_cellIdx;
} xHmdrEntry;
typedef vector<xHmdrEntry>		vHmdrEntry;
typedef vHmdrEntry::iterator	vHmdrEntry_it;

typedef enum _xMdrType {
	MT_ERR	= 0x00,
	MT_MDR	= 0x01,
	MT_GMDR	= 0x02,
	MT_HMDR	= 0x04,
	MT_FUZZYMDR = 0x08,
} xMdrType;

/* !!! HIMINI-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE & FUNCTION_HIMINI)

#include <stdint.h>
typedef struct _xBoostGeno {
	WISARD_pc *ca, *ct;
} xBoostGeno;

typedef struct _xFastEpi {
	vVariant&		Xa_snp;
	uint64_t	nepi;
	wsUint*		summary_good;
	wsUint*		summary_sig;
	wsReal*		best_score;
	wsUint*		best_partner;
	cExporter*	C_fe;
} xFastEpi;

class cBoostAnalysis : public cAnalysis
{
	wsUint		N_marker;
	xBoostGeno	*Xa_cnt;
	wsUint *N_ca, *N_ct, *N_caL, *N_ctL;
	void _comp(wsUint I, wsUint N_samp, wsUint N_snp,
		/*int *DistrCollection, */wsReal &minInteraction, wsReal &maxInteraction,
		wsUint &buffersize, wsReal thresholdRecord, wsUint &InteractionCount,
		wsUint **InteractionSNPpairs, wsReal **InteractionMeasureSNPpair,
		wsUint *N_ca, wsUint *N_ct, wsUint *N_caL, wsUint *N_ctL,
		wsUint **Na_genoMD, wsUint **Na_genoMDs, wsReal *Ra_phe, char **Na_data,
		xFastEpi *Xp_e);
public:
	cBoostAnalysis(cIO *Cp_inpIO);
	~cBoostAnalysis();
	void run();
};

class cMdrAnalysis : public cAnalysis
{
	wsUint				N_anaSamp;
	wsReal*				Ra_anaY;
	wsReal*				Ra_sc;
	cResult**			Cp_res;
	xMdrType			X_mdrType;
	cSetManagerAnalysis	*Cp_anaSetMgr;
	/* --cv */
	wsUint				N_anaCase, N_anaCtrl;
	wsUint*				Na_idxCase;
	wsUint*				Na_idxCtrl;
	wsMat				Ra_permY;
public:
	cMdrAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSetMgr);
	~cMdrAnalysis();
	void run();
	cResult*	getResult(wsUint N_idx=0);
	wsRealCst*	getAnaY() { return Ra_anaY; };
	xMdrType	getMdrType() { return X_mdrType; };
	wsRealCst*	getScore() { return Ra_sc; };
	wsMat		getPermY() { return Ra_permY; };
};

#else

/* DUMMY DEFINITION */
class cMdrAnalysis : public cAnalysis
{
public:
	cMdrAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSetMgr) : cAnalysis(Cp_inpIO) {}
	~cMdrAnalysis() {}
	void		run() {}
	cResult*	getResult() { return NULL; }
	xMdrType	getMdrType() { return MT_ERR; }
	wsRealCst*	getScore() { return NULL; };
};

/* DUMMY DEFINITION */
class cBoostAnalysis : public cAnalysis
{
	void _comp(wsUint i, wsUint N_snp);
public:
	cBoostAnalysis(cIO *Cp_inpIO) : cAnalysis(Cp_inpIO) {}
	~cBoostAnalysis() {}
	void run() {}
	cResult* getResult() { return NULL; }
};

#endif

} // End namespace ONETOOL

#endif
