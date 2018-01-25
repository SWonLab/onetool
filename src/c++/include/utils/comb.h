#pragma once
#ifndef __WISARD_COMB_H__
#define __WISARD_COMB_H__

namespace ONETOOL {

class cComb
{
	char	B_stop;
	wsUint	N_sz, N_order;
	wsUint	*Na_comb;
public:
	cComb();
	cComb(wsUint N_inpOrder, wsUint N_inpSz);
	~cComb();
	wsUint	get(wsUint N_len, wsUint *Np_comb, wsUint N_szJmp=0);
	void	init(wsUint N_inpOrder, wsUint N_inpSz);
	void	uninit();
};

} // End namespace ONETOOL

#endif
