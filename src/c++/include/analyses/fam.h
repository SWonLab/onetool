#pragma once
#include "global/analysis.h"

namespace ONETOOL {

struct xNucFamily
{
	xSample		*Xp_pat, *Xp_mat;	///< Pointer of xSample for parents
	vSampPtr	Xp_childs;			///< Pointer of xSample for childs

	wsUint getMissingChildSize()
	{
		wsUint N_missing = 0;
		FOREACH (vSampPtr_it, Xp_childs, it)
			if ((*it)->B_isMissing)
				N_missing++;

		return N_missing;
	}
};

typedef vector<xNucFamily*>		vNucFamPtr;
typedef vNucFamPtr::iterator	vNucFamPtr_it;
typedef map<string,xNucFamily>	mNucFam;
typedef mNucFam::iterator		mNucFam_it;

class cFamStrAnalysis : public cAnalysis
{
	mNucFam		Xm_couple2data;
	mSampPtr	Xm_missingFounderIID2data;
	wsUint*		Na_famIndices;				///< [#samp] Index of belonged family
	char		B_sameFamStr;
public:
	cFamStrAnalysis(cIO *Cp_inpIO);
	~cFamStrAnalysis();

	void		run();
	void		_buildNucFam(xSample &X_husband, xSample &X_wife);
	wsUint*		getFamIndices();
	mNucFam&	getNucFamData();
	mSampPtr&	getMissingFounderData();
};

class cTDTanalysis : public cAnalysis
{
	cFamStrAnalysis *Cp_anaFS;
public:
	static void	_test(wsUint N_pheno, wsReal **Ra_pheno, char **Na_data,
		mNucFam &Xm_nf, xVariant &X_snp, wsUint idx, wsReal *Rp_rets);
	cTDTanalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS);
	~cTDTanalysis();
	void	run();
	cFamStrAnalysis*
			getFamStrAna() { return Cp_anaFS; }
};

class cSDTanalysis : public cAnalysis
{
	cFamStrAnalysis *Cp_anaFS;
public:
	static void	_test(wsUint N_pheno, wsReal **Ra_pheno, char **Na_data,
		mNucFam &Xm_nf, xVariant &X_snp, wsUint idx, wsReal *Rp_rets);
	cSDTanalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS);
	~cSDTanalysis();
	void	run();
	cFamStrAnalysis*
			getFamStrAna() { return Cp_anaFS; }
};

} // End namespace ONETOOL
