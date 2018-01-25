//==========================================================================
// File:      vcf.ipp
//
// Author:    Sungyoung Lee
//
// History:   Initial implementation                              syl Dec 10
//
// Notes:     Inline functions in cVcfIO class in vcf.h
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "global/mem.h"

namespace ONETOOL {
namespace INPUT  {

inline void
cVcfIO::parseFormatField(char *Sp_fmt, xVCFfield &X_fdef)
{
    UINT_t  N_idx = 1;
    char    *a, *b;

    for (a = Sp_fmt - 1; a; N_idx++) {
        b = strchr(a + 1, ':');

        /* Check out what of a+1 have */
        if (a[1] == 'G' && a[2] == 'T')
            X_fdef.N_idxGT = N_idx;
        /* Genotype quality */
        else if (a[1] == 'G' && a[2] == 'Q')
            X_fdef.N_idxGQ = N_idx;
        /* Read depth */
        else if (a[1] == 'D' && a[2] == 'P')
            X_fdef.N_idxDP = N_idx;
        a = b;
    }

    /* Cannot found GT field */
    //  return 0xffffffff;
}

inline void
cVcfIO::loadVCFinfo(char *Sp_val, mVar &Xa_info)
{
    /* Get the structure of --prevalence */
    char        *a, *b;
    //UINT_t        N_prev      = 0;
    //char      **Sp_ret    = NULL;
    if (Sp_val == NULL) return;
    char        *Sp_prev = strdup(Sp_val);

    /* If --prevalence is not given, just return */
    if (Sp_prev[0] == '\0') return;

    /* Process prevalence string into real */
    //LOG("%d prevalences detected\n", N_prev);
    //  MULTI_MALLOC(Sp_ret, char*, N_prev);
    //  N_prev = 0;
    for (a = Sp_prev - 1; a; a = b) {
        b = strchr(a + 1, ';');
        if (b) *b = '\0';
        /* Find = */
        char *c = strchr(a + 1, '=');
        if (c) *(c++) = '\0';
        string S_name = a + 1;
        xOperand X;
        X.X_type = OTP_NA;
        if (c) str2op(c, X);
        Xa_info.insert(make_pair(S_name, X));
    }
}

inline char**
cVcfIO::loadVCFelements(char *Sp_val, UINT_t *Np_var)
{
    /* Get the structure of --prevalence */
    char        *a, *b;
    UINT_t      N_prev = 0;
    char        **Sp_ret = NULL;
    if (Sp_val == NULL) {
        *Np_var = 0;
        return NULL;
    }
    char        *Sp_prev = strdup(Sp_val);

    /* If --prevalence is not given, just return */
    if (Sp_prev[0] == '\0') {
        *Np_var = 0;
        return NULL;
    }

    /* Get the number of given prevalences */
    for (a = Sp_prev - 1; a; a = b) {
        b = strchr(a + 1, ':');
        N_prev++;
    }

    /* Process prevalence string into real */
    //LOG("%d prevalences detected\n", N_prev);
    wsAlloc(Sp_ret, char*, N_prev);
    N_prev = 0;
    for (a = Sp_prev - 1; a; a = b) {
        b = strchr(a + 1, ':');
        if (b) *b = '\0';
        Sp_ret[N_prev++] = a + 1;
    }

    *Np_var = N_prev;
    return Sp_ret;
}

inline char**
cVcfIO::loadVCFalleles(char *Sp_val, UINT_t *Np_var)
{
    /* Get the structure of --prevalence */
    char        *a, *b;
    UINT_t      N_prev = 0;
    char        **Sp_ret = NULL;
    if (Sp_val == NULL) {
        *Np_var = 0;
        return NULL;
    }
    char        *Sp_prev = strdup(Sp_val);

    /* If --prevalence is not given, just return */
    if (Sp_prev[0] == '\0') {
        *Np_var = 0;
        return NULL;
    }

    /* Get the number of given prevalences */
    for (a = Sp_prev - 1; a; a = b) {
        b = strchr(a + 1, ',');
        N_prev++;
    }

    /* Process prevalence string into real */
    //LOG("%d prevalences detected\n", N_prev);
    wsAlloc(Sp_ret, char*, N_prev);
    N_prev = 0;
    for (a = Sp_prev - 1; a; a = b) {
        b = strchr(a + 1, ',');
        if (b) *b = '\0';
        Sp_ret[N_prev++] = a + 1;
    }

    *Np_var = N_prev;
    return Sp_ret;
}

} // End namespace INPUT
} // End namespace ONETOOL
