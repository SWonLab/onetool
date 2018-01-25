#pragma once
#ifndef __WISARD_SAMPVAR_H__
#define __WISARD_SAMPVAR_H__
#include <string>

namespace ONETOOL {

class cIO;

typedef enum _xReservedCol {
	RC_PROBAND,
	RC_TWIN,
	RC_POPGROUP,
	RC_SAMPLEWEIGHT,
} xReservedCol;

typedef enum _xVarType
{
	WISARD_VAR_UNDET,
	WISARD_VAR_REAL,
	WISARD_VAR_FACTOR
} xVarType;

struct xVarRange
{
	int		N_sIdx, N_eIdx;
	//	xPheno	*Xp_sPhe, *Xp_ePhe;
	char	R_sEQ, R_eEQ;
	wsReal	R_s, R_e;
	char	S_varName[256];
	char	B_filt;
	char	*Sp_varValue;
};

typedef vector<xVarRange>	vVarRng;
typedef vVarRng::iterator	vVarRng_it;

typedef struct _xPheno
{
	std::string	S_name;
	std::string	S_fn;
	wsUint		N_idx;
	vVarRng		Xa_filters;
} xPheno;

typedef struct _xCovar
{
	xVarType	X_type;
	wsUint		N_szFactor;

	/* Variable name */
	char		*Sp_varName;
	/* Variable index IN input file */
	int			N_idxFile;
	/* Variable index IN data matrix */
	int			N_idx;
	/* Field value to be baseline if it is factor */
	char		*Sp_bl;
	/* Times */
	wsReal		R_times;
	/* A*B, pos of B in FILE */
	int			N_idxMul;
	/* File name came from */
	char		*S_fn;
	/* Value for flip */
	wsReal		R_v0, R_v1;

	vVarRng		Xa_filters;
} xCovar;

typedef vector<xPheno>		vPheno;
typedef vPheno::iterator	vPheno_it;

typedef vector<xCovar>		vCovar;
typedef vCovar::iterator	vCovar_it;

void loadSampVar(cIO* Cp_IO);

} // End namespace ONETOOL

#endif