#ifndef __WISARD_VIS_H__
#define __WISARD_VIS_H__
#pragma once

namespace ONETOOL {

void qqVariant(wsStrCst S_title, wsStrCst S_prefix, wsMat Ra_pvals, cIO *Cp_IO);
void hist(wsStrCst S_inpTitle, wsStrCst S_prefix, wsStrCst S_xLab, wsVec Ra_obs,
	wsUint N_obs);
void mhtVariant(wsStrCst S_inpTitle, wsStrCst S_prefix, wsMat Ra_pvals,
	int** Na_chrs, cIO* Cp_IO);

} // End namespace ONETOOL

#endif