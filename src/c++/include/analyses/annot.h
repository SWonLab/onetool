#ifndef __WISARD_ANNOT_H__
#define __WISARD_ANNOT_H__

#pragma once
#include "global/analysis.h"
#include "global/io.h"
#include <string>
using namespace std;

namespace ONETOOL {

typedef struct _xAnnoGene {
	wsUint	N_start, N_end;
	char	*S_name;
} xAnnoGene;

class cAnnotAnalysis : public cAnalysis
{
	vector<xAnnoGene> *Xa_annoGene;
	void _loadAnnoGene(wsStrCst Sp_gene);
	void	anno(xVariant &X_snp);
public:
	cAnnotAnalysis(cIO *Cp_IO);
	~cAnnotAnalysis();
	void	run();
};

cAnnotAnalysis* wsAnnot(cIO *Cp_IO);
void ANNOclear();

} // End namespace ONETOOL

#endif

