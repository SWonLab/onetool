#ifndef __WISARD_SETMGR_H__
#define __WISARD_SETMGR_H__
#pragma once
#include "global/analysis.h"

namespace ONETOOL {

typedef map<string,vector<int> >	mGeneDef;
typedef mGeneDef::iterator			mGeneDef_it;
typedef map<string,vStr>			mGeneSet;
typedef mGeneSet::iterator			mGeneSet_it;

class cVariantMap;
class cSetManagerAnalysis : public cAnalysis
{
	mGeneDef	Xa_gDef;
	mGeneSet	Xa_gsDef;
	mDataIdx	Xa_gsSize;
	FILE		**Ha_fpGene;
	FILE		*H_fpGset;
	wsUint		N_defGene;
	vStr		Sa_gDefPath;
protected:
	int			_initRangedGeneDef(FILE *H_fp, vVariant &Xa_snp);
	int			_initRefSeqGeneDef(FILE *H_fp, vVariant &Xa_snp);
	int			_initPairedOrListGeneDef(FILE *H_fp, vVariant &Xa_snp,
		char B_mustListForm=0);
	int			_initPLINKGeneDef(FILE *H_fp, vVariant &Xa_snp);
public:
	cSetManagerAnalysis(cIO *Cp_inpIO);
	cSetManagerAnalysis(cIO *Cp_inpIO, cVariantMap *Cp_inpMM);
	~cSetManagerAnalysis();
	
	void		run();
	void		summary();
	void		exportSet();
	mGeneDef&	getGeneDef();
	mGeneSet&	getGeneSetDef();
	mDataIdx&	getGeneSetOrigSize();
	vInt&		getSNPindicesByName(wsStrCst S_gName);
};

} // End namespace ONETOOL

#endif
