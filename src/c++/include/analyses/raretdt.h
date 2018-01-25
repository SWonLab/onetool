#pragma once
#ifndef __WISARD_RARETDT_H__
#define __WISARD_RARETDT_H__

#include "global/analysis.h"
#include "utils/stat.h"

namespace ONETOOL {

typedef struct _gsl_rng {
	wsReal N_s;
	wsReal N_e;
} gsl_rng;

/** Others **/
inline void trimNonSenseSites(vector<vReal>& genotype, const vBool& criteria)
{
	for(wsUint i=0; i<genotype.size(); i++){
		for(wsUint j=0; j<genotype[i].size(); j++){
			if(!criteria[j]) genotype[i][j] = WISARD_NA;
		}
	}
	return;
}

struct shuffleOutput
{
	shuffleOutput() {}
	shuffleOutput(const vector<vReal>& table, const vector<vReal>& untrans): tdtTable(table), untransmitted(untrans) {}
	vector<vReal> tdtTable, untransmitted;
};

vector<vReal> phasedToUnphased(const vector<vReal>& phased);

/** Base Trio **/
class  trio
{
public:
	trio() {}
	virtual void load() {}
	virtual ~trio() {}

	bool checkInformTrio() const {return __isInformative;}
	vector<vReal> getTdtTable() const {return __tdtTable;}
	vector<vReal> getNuTransPatChrs() const {return __unTransParentalChrs;}
	vBool getAnalyzed() const {return __tobeAnalyzed;}
	vBool getDenovo() const {return __denovoSite;}
	wsMat	getGeno() { return Ra_geno; }

	virtual void tdtTableCount() {}
	virtual shuffleOutput shuffle(std::string GenoHapo, gsl_rng* gslr);
	virtual shuffleOutput shuffleWithGenos(std::string GenoHapo, gsl_rng* gslr, vector<vReal>& genos);

protected:
	wsMat		Ra_geno;
	wsUint		N_pat, N_mat, N_chi;
	vInt*		Xv_set;
	cIO*		Cp_IO;
	vBool __tobeAnalyzed;
	bool __skipMiss;
	vector<vReal> __unTransParentalChrs;
	// four vectors for phased trio, first two vector is for father: v[0] is B count; v[1] is C count; Two vector for unphased trio
	vector<vReal> __tdtTable;
	vBool __denovoSite;
	bool __PesudoGenoReady;
	bool __isInformative; // check if this trio informative
	vector<vReal> __pesudoGeno; // genotype for shuffle

	vector<vReal> m_GenoShuffle(gsl_rng* gslr);
	vector<vReal> m_HapoShuffle(gsl_rng* gslr);
//	virtual void m_setUpPesudoGeno(gsl_rng* gslr) {} // For shuffle: convert unphased to phased, guess missing hapotype
};

/** Phased Trio **/
// class phasedTrio: public trio
// {
// public:
// 	phasedTrio() {}
// 	//~phasedTrio() {}
// 	void load(cIO *Cp_IO, vInt& Xv_set, const vectorL& tobeAnalyzed, bool skipMiss);
// 	void tdtTableCount();
// 	shuffleOutput shuffle(std::string GenoHapo, gsl_rng* gslr);
// 
// protected:
// 	vectorUI __chrOrigin;
// 	bool m_ChromOrigin(bool allowMissing, bool allowDenovo);
// 	void m_phasedTdtCount(vectorF& transParent, vectorF& unTransParent, vectorF& Child, wsUint base);
// 	void m_setUpPesudoGeno(gsl_rng* gslr);
// };

bool ChromEqual(const vReal& ChromOne, const vReal& ChromTwo, bool Missing, bool Denovo);

/** unPhased Trio **/
class unPhasedTrio: public trio
{
public:
	unPhasedTrio() {}
	void load(cIO *Cp_IO, vInt& Xv_set, const vBool& tobeAnalyzed, bool skipMiss,
		wsUint N_idxPat, wsUint N_idxMat, wsUint N_idxChi);
	void tdtTableCount();
	shuffleOutput shuffle(std::string GenoHapo, gsl_rng* gslr);
	shuffleOutput shuffleWithGenos(std::string GenoHapo, gsl_rng* gslr, vector<vReal>& genos);

protected:
	void m_unPhasedTdtCount(double fat, double mot, double kid, wsUint index);
	void m_setUpPesudoGeno(gsl_rng* gslr);
	double guessMissing(double knownParent, double kid, gsl_rng* gslr);
};

typedef trio* trioPtr;
//typedef phasedTrio* phasedTrioPtr;
typedef unPhasedTrio* unPhasedTrioPtr;
typedef std::vector<trioPtr>::iterator trioIterator;

struct tdtResult
{
	tdtResult() {}
	tdtResult(double C_MZ, double C_CMC, double P_MZ, double P_CMC){
		mz_chi2=C_MZ;
		cmc_chi2=C_CMC;
		mz_pval=P_MZ;
		cmc_pval=P_CMC;
	}
	void print() {printf(">> %s=%+-4.2f %s=%+-4.2f %s=%+-4.2f %s=%+-4.2f\n", "mz_chi2", mz_chi2,"cmc_chi2", cmc_chi2,"mz_pval", mz_pval,"cmc_pval", cmc_pval);}
	double mz_chi2, cmc_chi2, mz_pval, cmc_pval;
};

struct rvTDTshuffleOutput
{
	rvTDTshuffleOutput() {}
	rvTDTshuffleOutput(const vector<vReal>& __unTrans, const vector<vReal>& __tdtBtable, const vector<vReal>& __tdtCtable){
		unTrans=__unTrans;
		tdtBtable=__tdtBtable;
		tdtCtable=__tdtCtable;
	}
	vector<vReal> unTrans, tdtBtable, tdtCtable;
};

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cFamStrAnalysis;
class cRvtdtAnalysis : public cAnalysis
{
	//	Json::Value logRoot;
	std::vector<trioPtr> __trios;
	/* FB-SKAT */
	wsReal		Echp1p2[10][10];
	bool		__skipMiss;
	double		__missCutoff;
	vReal		__samMafs;
	wsUint		__nInformTrios;
	vReal		__missingRatio;
	vBool		__varToBeAnalyzed;
	vector<wsUint>	__denovoCount;
	vector<wsUint>	__effTriosBySite;
	vReal		__wssWeight;
	cSetManagerAnalysis
				*Cp_anaSM;
	cFamStrAnalysis
				*Cp_anaFS;

	void	m_preprocess(vInt& Xv_set);
	wsUint	__getInformVarNum();
	char	__load(std::map<std::string, double>& pvalues, vInt& Xv_set);
	void	m_tdtTableLoad(vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable);
	vReal	m_wssWeight(const vector<vReal>& unTrans);
	vReal	m_singleSiteP(const vector<vReal> &__tdtBtable, const vector<vReal> &__tdtCtable) const;

	rvTDTshuffleOutput m_shuffle(std::string shuffleMethod, gsl_rng* gslr);
	tdtResult m_rvAggreate(const vector<vReal>& tdtBtable, const vector<vReal>& tdtCtable, const vBool& tobeAnalyzed, const vReal& Weight, bool allowDenovo) const;
	tdtResult m_VT(const vector<vReal>& tdtBtable, const vector<vReal>& tdtCtable, const vBool& tobeAnalyzed, const vReal& weight, vInt& Xv_set) const;

	vReal m_noPermut(const vBool& tobeAnalyzed,
		vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable, vInt& Xv_set);
	vReal m_Permut(std::string shuffleMethod,
		const vBool& tobeAnalyzed, wsUint adaptive, wsUint N_perm,
		double R_alpha, vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable,
		gsl_rng* gslr, vInt& Xv_set);
	vReal m_MzCmcPermut(vBool& tobeAnalyzed, wsUint adaptive, wsUint PermutateTimes, double alpha, gsl_rng* gslr, double upperMaf,
		vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable, vInt& Xv_set);
	vReal m_CmpPermut(std::string GenoHapo, const vBool& tobeAnalyzed, wsUint adaptive, wsUint PermutateTimes, double alpha, gsl_rng* gslr, double lowerMaf, double upperMaf);
public:
	cRvtdtAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSM, cFamStrAnalysis *Cp_inpAnaFS, bool phased, const vReal& popMafs, bool skipMiss, double missCutoff);
	cRvtdtAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSM, cFamStrAnalysis *Cp_inpAnaFS, bool phased, bool skipMiss, double missCutoff);
	~cRvtdtAnalysis(){}
	void run();
	void tdtTest(std::string GenoHapo, wsReal R_mafL, wsReal R_mafU,
		wsUint adaptive, wsUint N_perm, wsReal R_alpha, std::map<std::string, double>& pvalues,
		std::vector<std::string>& tests, vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable,
		vInt& Xv_set);
//	Json::Value getJsonLog() const {return logRoot;}
};

#else

/* DUMMY DEFINITION */
class cFamStrAnalysis;
class cRvtdtAnalysis : public cAnalysis
{
public:
	cRvtdtAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSM, cFamStrAnalysis *Cp_inpAnaFS, bool phased, const vReal& popMafs, bool skipMiss, double missCutoff)
		: cAnalysis(Cp_inpIO) {}
	cRvtdtAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSM, cFamStrAnalysis *Cp_inpAnaFS, bool phased, bool skipMiss, double missCutoff)
		: cAnalysis(Cp_inpIO) {}
	~cRvtdtAnalysis(){}
	void run(){};
	void tdtTest(std::string GenoHapo, wsReal R_mafL, wsReal R_mafU,
		wsUint adaptive, wsUint N_perm, wsReal R_alpha, std::map<std::string, double>& pvalues,
		std::vector<std::string>& tests, vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable,
		vInt& Xv_set){};
	//	Json::Value getJsonLog() const {return logRoot;}
};

#endif

inline vReal m_colTotal(const vector<vReal>& table)
{
	vReal count(table[0].size(),0);
	for(wsUint i = 0; i < table.size(); i++){
		for(wsUint j = 0; j < table[i].size(); j++){
			count[j] += (table[i][j]>=100)?(table[i][j]-100):table[i][j];
		}
	}
	return count;
}

#define BIN(gslr) (wsUnifrand() * (gslr->N_e - gslr->N_s) + gslr->N_s < 0.5)
#define BINV(gslr, v) (wsUnifrand() * (gslr->N_e - gslr->N_s) + gslr->N_s < v)

inline void m_updateCount(vector<wsUint>& permCount, vReal& oriRes, vReal& permRes, gsl_rng* gslr)
{
	if(permCount.size() != oriRes.size() || permCount.size() != permRes.size()) {return;}
	for(wsUint i=0; i<permCount.size(); i++){
		double diff=permRes[i] - oriRes[i];
		if(diff>1e-6) permCount[i]++;
		else if(diff<=1e-6 && diff>=-1e-6) {
			if(BIN(gslr)) permCount[i]++;
		}
		else {;}
	}
	return;
}

inline void m_updatePvalue(vReal& permPval, vector<wsUint>& permCount, wsUint iPermut, double alpha)
{
	for(wsUint i=0; i<permCount.size(); i++){
		double pval = (permCount[i]+1)*1.0/((iPermut+1)*1.0);
		double sigma = sqrt(pval*(1-pval)/iPermut);
		double beta = 0.05;
		double gs = qnorm(1.0-beta/2.0, 0, sigma);
//		double gs = gsl_cdf_gaussian_Pinv(1.0-beta/2.0, sigma);

		permPval[i] = (pval-gs > alpha)? pval:9.0;
	}
	return;
}

inline void m_ifTrueThenAddOne(const vBool& logic, vector<wsUint>& count)
{
	if(count.size()==0) count.resize(logic.size(), 0);
	if(count.size() != logic.size()) return;
	for(wsUint i=0; i<logic.size(); i++){
		if(logic[i]) count[i]++;
	}
	return;
}

template<class T> std::vector<T> eraseTrue (std::vector<T>& vec, const vBool stdV)
{
	std::vector<T> output;
	for(wsUint i=0; i<vec.size(); i++){
		if(!stdV[i]) output.push_back(vec[i]); // remove element that is true in stdV
	}
	return output;
}

} // End namespace ONETOOL

#endif
