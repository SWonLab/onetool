#pragma once
#include "global/io.h"

namespace ONETOOL {

vector<vInt> _buildPropagatePath(cIO *Cp_IO, wsUint **Np_t1, wsUint **Np_t2,
	mDataIdx *Xp_insIdx=NULL, char B_noMissing=1);

class cSimTrioIO : public cIO
{
public:
	cSimTrioIO(int N_simFam, int N_simSNP);
	wsStrCst getFormat() { return "sim.trio"; }
};

class cSymMatrix;
class cSimFamIO : public cIO
{
	cSymMatrix*		Mp_V;
	vector<vInt>	Xv_prop;
	vReal			Xv_maf;
	wsUint			*Na_t1, *Na_t2;
	void			_procFAM(char *Sp_buf, wsUint i, wsUint &j);
public:
	int				_sim(vector<vInt> &Xv_prop, wsUint i, wsReal R_maf);
	int				_simsig(vector<vInt> &Xv_prop, wsUint i, wsReal R_maf,
		cSymMatrix *Mp_V, cVector &V_y, wsReal R_beta);
	cSimFamIO(wsStrCst S_fn, char B_inDry);
	~cSimFamIO();
	vector<vInt>&	getPropIndices() { return Xv_prop; }
	vReal&			getMAFs() { return Xv_maf; }
	cSymMatrix*		getV() { return Mp_V; }
	wsStrCst getFormat() { return "sim.fam"; }
};

} // End namespace ONETOOL
