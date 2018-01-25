#pragma once
#ifndef __WISARD_SORT_H__
#define __WISARD_SORT_H__

#include <stdio.h>
#include "global/datatype.h"

namespace ONETOOL {

struct xRealSort {
	wsUint i;
	wsReal V;
};

struct xStrSort {
	wsUint i;
	char *V;
};

struct xUintSort {
	xUintSort() {
		i = 0;
		V = 0;
	};
	xUintSort(wsUint _i, wsUint _V) {
		i = _i;
		V = _V;
	};
	wsUint i;
	wsUint V;
};

int			sort_real_unorder(const void *a, const void *b);
int			sort_real(const void *a, const void *b);
int			sort_real_desc(const void *a, const void *b);
int			sort_str(const void *a, const void *b);
int			sort_str_desc(const void *a, const void *b);
int			sort_natstr(const void *a, const void *b);
int			sort_natstr_desc(const void *a, const void *b);
int			sort_uint(const void *a, const void *b);
int			sort_uint_desc(const void *a, const void *b);
int			sortabs_real(const void *a, const void *b);
int			sortabs_real_unorder(const void *a, const void *b);
xRealSort*	buildRealSort(wsVecCst Ra_vec, wsUint N_szVec, wsUint *Np_szAfter=NULL);
xRealSort*	buildArealSort(wsVecCst Ra_vec, wsUint N_szVec, wsUint *Np_szAfter=NULL);

/* Ascend-order sorting function for 3'-UTR length */
int			asc_sort_rsname(const void *X_l, const void *X_r);
/* ascend-order sorting function for double array */
int			asc_sort_real(const void *a, const void *b);

/*
 *
 * Natural sort definition
 * 
 */

typedef char nat_char;

int strnatcmp(nat_char const *a, nat_char const *b);
int strnatcasecmp(nat_char const *a, nat_char const *b);

} // End namespace ONETOOL

#endif
