#include "global/common.h"
#include "utils/util.h"
#include "utils/comb.h"

namespace ONETOOL {

cComb::cComb()
{
	Na_comb = NULL;
}

cComb::cComb(wsUint N_inpOrder, wsUint N_inpSz)
{
	init(N_inpOrder, N_inpSz);
}

cComb::~cComb()
{
	DEALLOC(Na_comb);
}

wsUint cComb::get(wsUint N_len, wsUint *Np_comb, wsUint N_szJmp)
{
	wsUint i;

	if (Na_comb == NULL)
		halt("cComb was not initialized!");
	if (N_szJmp == 0)
		N_szJmp = N_order;

	/* Perform */
	wsUint *com = Na_comb;
	wsUint k = N_order;
	wsUint n = N_sz;
	for (i=0 ; i<N_len && com[k - 1] < n ; i++) {
		/* Copy generated combination */
// 		pverbose("Generated comb [%d] : ", i);
// 		for (j=0 ; j<N_order ; j++)
// 			verbosenf(" %d", Na_comb[j]);
// 		verbosenf("\n");
		memcpy(Np_comb, Na_comb, sizeof(wsUint)*N_order);
		Np_comb += N_szJmp;

		/* Increase last index */
		int t = k - 1;
		while (t != 0 && com[t] == n - k + t) t--;
		com[t]++;
		for (wsUint i = t + 1; i < k; i++) com[i] = com[i - 1] + 1;
 	}
	//pverbose("Number of generated combination : %d\n", i);

	return i;
}

void cComb::init(wsUint N_inpOrder, wsUint N_inpSz)
{
	N_order	= N_inpOrder;
	N_sz	= N_inpSz;
	B_stop	= 0;

	wsAlloc(Na_comb, wsUint, N_order);
	for (wsUint i=0 ; i<N_order ; i++)
		Na_comb[i] = i;
}

void cComb::uninit()
{
	DEALLOC(Na_comb);
}

} // End namespace ONETOOL
