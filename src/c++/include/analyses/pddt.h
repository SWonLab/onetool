#ifndef __WISARD_PDDT_H__
#define __WISARD_PDDT_H__
#pragma once

#include "global/io.h"
#include "global/analysis.h"

namespace ONETOOL {

typedef struct _xPDDTnode
{
	int			N_depth;
	vSampPtr	Xa_next;
} xPDDTnode;

class cPDDTAnalysis : public cAnalysis
{
	wsUint	N_sampPDDT;
	void	_getEmpiCorr(wsReal **Ra_pddt, wsFloat **Ra_data,
		wsUint N_1, wsUint N_2, wsUint N_SNP);
public:
	wsReal	**Ra_pddt; 	///< Array[size #sample*#sample], PDDT statistics for dataset
	cPDDTAnalysis(cIO *Cp_inpIO);
	~cPDDTAnalysis();
	void run();

	wsReal**	getPDDT() { return Ra_pddt; }
};

int buildGraph(map<string,xPDDTnode> &X_map, xSample *Xp_t, char B_fromMe,
			   int N_depth=0);

} // End namespace ONETOOL

#endif
