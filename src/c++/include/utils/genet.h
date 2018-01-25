#pragma once
#ifndef __WISARD_GENET_H__
#define __WISARD_GENET_H__

#include "utils/util.h"

/* Maximum number of chromosomes allowed */
#define MAX_NCHR		40
#define NCHR_SPECIES	(wsUint)(OPT_NUMBER(maxNumChr))
#define NAUTO_SPECIES	(wsUint)(OPT_NUMBER(maxNumAutoChr))
#define CHR_X			(NAUTO_SPECIES+1)

namespace ONETOOL {

/* Chr# into integer */
int		getChr(char *S_buf);
int		getChrAuto(char *S_buf);
int		getChrAuto(int N_chr);

/* Get Chr# from integer */
const char * getChrName2(int N_chr);

} // End namespace ONETOOL

#endif
