#include <string.h>
#include <math.h>
#include "utils/util.h"
#include "utils/sort.h"

namespace ONETOOL {

int sort_real_unorder(const void *a, const void *b)
{
	wsReal res = (*((wsReal *)a)) - (*((wsReal *)b));
	if (res > 0)
		return -1;
	else if (res < 0)
		return 1;
	else
		return 0;
}

int sort_real(const void *a, const void *b)
{
	wsReal res = ((xRealSort *)a)->V - ((xRealSort *)b)->V;
	if (res > 0)
		return 1;
	else if (res < 0)
		return -1;
	else
		return 0;
}

int sort_real_desc(const void *a, const void *b)
{
	wsReal res = ((xRealSort *)b)->V - ((xRealSort *)a)->V;
	if (res > 0)
		return 1;
	else if (res < 0)
		return -1;
	else
		return 0;
}

int sort_str(const void *a, const void *b)
{
	return strcmp(((xStrSort *)a)->V, ((xStrSort *)b)->V);
}

int sort_str_desc(const void *a, const void *b)
{
	return strcmp(((xStrSort *)b)->V, ((xStrSort *)a)->V);
}

int sort_natstr(const void *a, const void *b)
{
	return strnatcmp(((xStrSort *)a)->V, ((xStrSort *)b)->V);
}

int sort_natstr_desc(const void *a, const void *b)
{
	return strnatcmp(((xStrSort *)b)->V, ((xStrSort *)a)->V);
}

int sort_uint(const void *a, const void *b)
{
	if (((xUintSort *)a)->V > ((xUintSort *)b)->V)
		return -1;
	else if (((xUintSort *)a)->V < ((xUintSort *)b)->V)
		return 1;
	else
		return 0;
}

int sort_uint_desc(const void *a, const void *b)
{
	if (((xUintSort *)a)->V > ((xUintSort *)b)->V)
		return 1;
	else if (((xUintSort *)a)->V < ((xUintSort *)b)->V)
		return -1;
	else
		return 0;
}

int sortabs_real(const void *a, const void *b)
{
	wsReal res = fabs(((xRealSort *)a)->V) - fabs(((xRealSort *)b)->V);
	if (res > 0)
		return -1;
	else if (res < 0)
		return 1;
	else
		return 0;
}

int sortabs_real_unorder(const void *a, const void *b)
{
	wsReal res = fabs(*((wsReal *)a)) - fabs(*((wsReal *)b));
	if (res > 0)
		return -1;
	else if (res < 0)
		return 1;
	else
		return 0;
}

xRealSort*	buildRealSort(wsVecCst Ra_vec, wsUint N_szVec,
	wsUint *Np_szAfter/*=NULL*/)
{
	xRealSort	*Xa_ret = NULL;

	wsAlloc(Xa_ret, xRealSort, N_szVec);
	if (Np_szAfter) {
		wsUint I=0;
		for (wsUint i=0 ; i<N_szVec ; i++) {
			if (isMissingReal(Ra_vec[i])) continue;
			Xa_ret[I].i = i;
			Xa_ret[I].V = Ra_vec[i];
			I++;
		}
		*Np_szAfter = I;
	} else for (wsUint i=0 ; i<N_szVec ; i++) {
		Xa_ret[i].i = i;
		Xa_ret[i].V = Ra_vec[i];
	}

	return Xa_ret;
}

xRealSort*	buildArealSort(wsVecCst Ra_vec, wsUint N_szVec,
	wsUint *Np_szAfter/*=NULL*/)
{
	xRealSort	*Xa_ret = NULL;

	wsAlloc(Xa_ret, xRealSort, N_szVec);
	if (Np_szAfter) {
		wsUint I=0;
		for (wsUint i=0 ; i<N_szVec ; i++) {
			if (isMissingReal(Ra_vec[i])) continue;
			Xa_ret[I].i = i;
			Xa_ret[I].V = fabs(Ra_vec[i]);
			I++;
		}
		*Np_szAfter = I;
	} else for (wsUint i=0 ; i<N_szVec ; i++) {
		Xa_ret[i].i = i;
		Xa_ret[i].V = fabs(Ra_vec[i]);
	}

	return Xa_ret;
}

wsReal* sseVpP(wsReal *Ra_v, wsUint N_sz, wsReal *Ra_p, wsUint N_c,
			   wsUint *Na_rEnd, wsUint *Na_cIdx, wsReal *Ra_mask/*=NULL*/)
{
	wsReal*	Ra_ret	= sseVector(N_c);
	wsUint	i		= 0;
	wsUint	N_s		= 0;
	wsUint	N_e;

	if (Ra_mask == NULL) do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e) {
			/* Sum up them */
			for (wsUint j=N_s ; j<N_e ; j++) {
				wsUint TC = Na_cIdx[j];
				Ra_ret[TC] += Ra_v[i] * Ra_p[j];
			}
		}

		/* Update start point */
		i++;
		if (i < N_sz)
			N_s = Na_rEnd[i];
		else break;
	} while (1); else do {
		/* Set end point */
		N_e = Na_rEnd[i];

		if (Ra_mask[i] && N_s != N_e) {
			/* Sum up them */
			for (wsUint j=N_s ; j<N_e ; j++) {
				wsUint TC = Na_cIdx[j];
				Ra_ret[TC] += Ra_v[i] * Ra_p[j];
			}
		}

		/* Update start point */
		i++;
		if (i < N_sz)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	return Ra_ret;
}

/* ascend-order sorting function for double array */
int asc_sort_real(const void *a, const void *b)
{
	wsReal _a = *((wsReal *)a);
	wsReal _b = *((wsReal *)b);

	return _a > _b ? 1 : (_a < _b ? -1 : 0);
}

/*
 *
 * Natural sort implementation
 * 
 */

/* These are defined as macros to make it easier to adapt this code to
 * different characters types or comparison functions. */
static inline int nat_isdigit(nat_char a)
{
	return isdigit((unsigned char) a);
}


static inline int nat_isspace(nat_char a)
{
	return isspace((unsigned char) a);
}


static inline nat_char nat_toupper(nat_char a)
{
	return (nat_char)toupper((unsigned char) a);
}



inline int compare_right(nat_char const *a, nat_char const *b)
{
	int bias = 0;
	
	/* The longest run of digits wins.  That aside, the greatest
	value wins, but we can't know that it will until we've scanned
	both numbers to know that they have the same magnitude, so we
	remember it in BIAS. */
	for ( ; ; a++,b++) {
	 if (!nat_isdigit(*a) && !nat_isdigit(*b))
		 return bias;
	 else if (!nat_isdigit(*a))
		 return -1;
	 else if (!nat_isdigit(*b))
		 return +1;
	 else if (*a < *b) {
		 if (!bias)
		   bias = -1;
	 } else if (*a > *b) {
		 if (!bias)
		   bias = +1;
	 } else if (!*a  &&  !*b)
		 return bias;
	}

	// Unreachable confirmed
	// return 0;
}


inline int compare_left(nat_char const *a, nat_char const *b)
{
	/* Compare two left-aligned numbers: the first to have a
	   different value wins. */
	for ( ; ; a++, b++) {
		if (!nat_isdigit(*a) && !nat_isdigit(*b))
			return 0;
		else if (!nat_isdigit(*a))
			return -1;
		else if (!nat_isdigit(*b))
			return +1;
		else if (*a < *b)
			return -1;
		else if (*a > *b)
			return +1;
	}

	// Unreachable confirmed
	// return 0;
}

inline int strnatcmp0(nat_char const *a, nat_char const *b, int fold_case)
{
	int ai, bi;
	nat_char ca, cb;
	int fractional, result;
	
//	assert(a && b);
	ai = bi = 0;
	while (1) {
	ca = a[ai]; cb = b[bi];

	/* skip over leading spaces or zeros */
	while (nat_isspace(ca))
		ca = a[++ai];

	while (nat_isspace(cb))
		cb = b[++bi];

	/* process run of digits */
	if (nat_isdigit(ca)  &&  nat_isdigit(cb)) {
		fractional = (ca == '0' || cb == '0');

		if (fractional) {
		  if ((result = compare_left(a+ai, b+bi)) != 0)
			return result;
		} else {
		  if ((result = compare_right(a+ai, b+bi)) != 0)
			return result;
		}
	}

	if (!ca && !cb) {
		/* The strings compare the same.  Perhaps the caller
			   will want to call strcmp to break the tie. */
		return 0;
	}

	if (fold_case) {
		ca = nat_toupper(ca);
		cb = nat_toupper(cb);
	}
	
	if (ca < cb)
		return -1;
	else if (ca > cb)
		return +1;

	++ai; ++bi;
	}
}

int strnatcmp(nat_char const *a, nat_char const *b)
{
	return strnatcmp0(a, b, 0);
}

/* Compare, recognizing numeric string and ignoring case. */
int strnatcasecmp(nat_char const *a, nat_char const *b)
{
	return strnatcmp0(a, b, 1);
}

} // End namespace ONETOOL
