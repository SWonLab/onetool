#pragma once
#include "utils/matrix.h"
#include "global/common.h"

namespace ONETOOL {

class cIO;
class cAnalysis;

class cImpute
{
protected:
	wsMat	Ra_data;
	cIO*	Cp_IO;
	char	B_useNorm;
public:
	cImpute(cIO *Cp_inpIO, char B_useNorm=0);
	~cImpute();
	cIO*			getIO() { return Cp_IO; }
	wsMat			getData() { return Ra_data; }
	virtual void	impute()=0;
};

class cMAFimpute : public cImpute
{
public:
	cMAFimpute(cIO *Cp_inpIO);
	~cMAFimpute();
	void impute();
};

class cCorImpute : public cImpute
{
	cSymMatrix	M_cor;
public:
	cCorImpute(cIO *Cp_inpIO, wsSym Ra_cor);
	~cCorImpute();
	void	impute();
	wsSym	getCor() { return M_cor.get(); }
};

} // End namespace ONETOOL
