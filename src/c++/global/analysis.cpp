#include "global/common.h"
#include "global/option.h"
#include "global/analysis.h"

namespace ONETOOL {

/* Running only once
 * Divide SNPs ALMOST equally to each thread */
int forVariant_equal(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	int			*Np_idx	= (int *)Vp_thrData;
	xAnaThread	*Cp_ana	= (xAnaThread *)Vp_shareData;
	cIO			*Cp_IO	= Cp_ana->Cp_IO;
	wsUintCst		N_vrt	= Cp_IO->sizeVariant();
	static wsUint
				N_proc	= 0;
	wsReal		R_szDiv	= (wsReal)N_vrt/(wsReal)N_thread;
	int			N_ret	= N_proc != 0 ? 0 : 1;
	int			i=0;
	wsReal		j=W0;

	switch (N_mode) {
	case DISTFUN_INIT:
		i = 0;
		j = W0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		if (N_ret == 0)
			return 0;
		for ( ; i<N_thread ; i++,j+=R_szDiv) {
			wsUint E = (wsUint)(j+R_szDiv+REAL_CONST(0.5));
			if (E > N_vrt)
				E = N_vrt;
			Np_idx[3*i] = (wsUint)(j+REAL_CONST(0.5));
			Np_idx[3*i+1] = E;
			Np_idx[3*i+2] = 0;
		}
		N_proc = N_ret = N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d variants processed\n", N_vrt, N_vrt);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}

	return N_ret;
}

int forAllSNP_ana(int N_mode, int N_thread, void *Vp_data,
	void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx = (int *)Vp_data;
	static wsUint N_idxSNP;
	cIO *Cp_IO = ((cAnalysis *)Vp_shareData)->getIO();
	int i=0;

	switch (N_mode) {
		/* Initializes N_idxSNP <- 0 */
		case DISTFUN_INIT:
			N_idxSNP = 0;
			break;
		/* Assign thread-specific value */
		case DISTFUN_DISTRIBUTE:
			for ( ; i<N_thread ; i++) {
				if ((N_idxSNP+i) >= Cp_IO->sizeVariant())
					break;
				Np_idx[i] = N_idxSNP+i;
				if (((N_idxSNP+i)%1000) == 0)
					notice("%d/%d SNPs processed\r", N_idxSNP+i, Cp_IO->sizeVariant());
			}
			N_idxSNP += N_thread;
			break;
		case DISTFUN_UNINIT:
			N_idxSNP = 0;
			break;
		case DISTFUN_AFTERLOOP:
			break;
		default:
			halt("Unsupported distFunc command [%d]", N_mode);
			break;
	}
	return i;
}

int forAllSample_ana(int N_mode, int N_thread, void *Vp_data,
	void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx = (int *)Vp_data;
	static wsUint N_idxSample;
	cIO *Cp_IO = ((cAnalysis *)Vp_shareData)->getIO();
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxSample = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			if ((N_idxSample+i) >= Cp_IO->sizeSample())
				break;
			Np_idx[i] = N_idxSample+i;
			if (((N_idxSample+i)%10) == 0)
				notice("%d/%d samples processed\r", N_idxSample+i, Cp_IO->sizeSample());
		}
		N_idxSample += N_thread;
		break;
	case DISTFUN_UNINIT:
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return i;
}

int forAllSample_anaXthr(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx = (int *)Vp_thrData;
	static	wsUint N_idxSample;
	cIO		*Cp_IO = ((xAnaThread *)Vp_shareData)->Cp_IO;
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxSample = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			if ((N_idxSample+i) >= Cp_IO->sizeSample())
				break;
			Np_idx[i] = N_idxSample+i;
			if (((N_idxSample+i)%10) == 0)
				notice("%d/%d samples processed\r", N_idxSample+i, Cp_IO->sizeSample());
		}
		N_idxSample += N_thread;
		break;
	case DISTFUN_UNINIT:
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return i;
}

} // End namespace ONETOOL
