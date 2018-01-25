#include "utils/vector.h"
#include "utils/matrix.h"
#include <math.h>
#include <stdarg.h>

namespace ONETOOL {

void exportVector(const char *S_ext, vInt& Nv_vec)
{
	wsUint N_sz = (wsUint)Nv_vec.size();
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	LOGoutput("Export vector [%s.%s] (size %d)\n", OPT_STRING(out), S_ext, N_sz);

	/* Do export */
	FOREACH(vInt_it, Nv_vec, i)
		Cp_mat->fmt("%d	", *i);
	Cp_mat->put("\n");

	delete Cp_mat;
}

void exportVector(const char *S_ext, vector<wsFloat>& Nv_vec)
{
	wsUint N_sz = (wsUint)Nv_vec.size();
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	LOGoutput("Export vector [%s.%s] (size %d)\n", OPT_STRING(out), S_ext, N_sz);

	/* Do export */
	FOREACH(vector<wsFloat>::iterator, Nv_vec, i)
		Cp_mat->fmt("%g	", *i);
	Cp_mat->put("\n");

	delete Cp_mat;
}

void exportVector(const char *S_ext, vStr& Nv_vec, wsUint N_sz)
{
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	LOGoutput("Export vector [%s.%s] (size %d)\n", OPT_STRING(out), S_ext, N_sz);

	/* Do export */
	FOREACH(vStr_it, Nv_vec, i)
		Cp_mat->fmt("%s\n", i->c_str());

	delete Cp_mat;
}

void exportVector(const char *S_ext, wsVecCst Rp_vec, wsUint N_sz, wsUint N_prec/*=0*/)
{
	if (Rp_vec == NULL) halt_fmt(WISARD_SYST_NULL_DATA, "Exported vector");
	char		S_fmt[64];
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	LOGoutput("Export vector [%s.%s] (size %d)\n", OPT_STRING(out), S_ext, N_sz);

	/* Do export */
	for (wsUint j = 0; j < N_sz; j++)
		Cp_mat->fmt(S_fmt, Rp_vec[j]);
	Cp_mat->put("\n");

	delete Cp_mat;
}

void exportVector(const char *S_ext, const wsFvec Rp_vec, wsUint N_sz, wsUint N_prec/*=0*/)
{
	if (Rp_vec == NULL) halt_fmt(WISARD_SYST_NULL_DATA, "Exported vector");
	char		S_fmt[64];
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	LOGoutput("Export vector [%s.%s] (size %d)\n", OPT_STRING(out), S_ext, N_sz);

	/* Do export */
	for (wsUint j=0 ; j<N_sz ; j++)
		Cp_mat->fmt(S_fmt, Rp_vec[j]);
	Cp_mat->put("\n");

	delete Cp_mat;
}

void exportUnitVector(const char *S_ext, wsUint N_sz)
{
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	LOGoutput("Export unit vector [%s.%s] (size %d)\n", OPT_STRING(out), S_ext, N_sz);

	/* Do export */
	for (wsUint j = 0; j < N_sz; j++)
		Cp_mat->put("1.0	");
	Cp_mat->put("\n");

	delete Cp_mat;
}

void exportValVector(const char *S_ext, wsReal R_val, wsUint N_sz, wsUint N_prec)
{
	char		S_buf[128];
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	LOGoutput("Export value[%g] vector [%s.%s] (size %d)\n", R_val,
		OPT_STRING(out), S_ext, N_sz);
	if (N_prec == 0)
		sprintf(S_buf, "%g	", R_val);
	else {
		char S_fmt[64];
		sprintf(S_fmt, "%%.%de	", N_prec);
		sprintf(S_buf, S_fmt, R_val);
	}

	/* Do export */
	for (wsUint j = 0; j < N_sz; j++)
		Cp_mat->put(S_buf);
	Cp_mat->put("\n");

	delete Cp_mat;
}

char cVector::B_report = 0;

cValVector::cValVector(wsUint N_inpSz, wsReal R_inpVal)
{
	N_sz = N_inpSz;
	R_val = R_inpVal;
	Ra_buf = NULL;
}

cValVector::cValVector(cValVector &&V_noname)
{
	if (V_noname.isDdealloc())
		setDontDealloc();
	Ra_buf = NULL;
	N_sz = V_noname.size();
	R_val = V_noname.getV();

	V_noname.setDontDealloc();
}

cValVector::~cValVector()
{
	sseFree(Ra_buf);
}

void cValVector::init(wsUint N_inpSz, wsReal *Ra_val, wsReal *Rp_val,
	char B_rand, char B_inpDDealloc)
{
	if (B_rand)
		halt("Unit vector can't be initiated with random");
	else if (Ra_val)
		halt("Unit vector can't be initiated with vector values");

	N_sz	= N_inpSz;
	R_val	= *Rp_val;
	Ra_buf	= NULL;
}

wsReal* cValVector::get()
{
	/* Realize buffer */
	if (Ra_buf == NULL) {
		Ra_buf = sseVector(N_sz);
		sseVinit(Ra_buf, N_sz, R_val);
	}
	return Ra_buf;
}

wsRealCst cValVector::sum(wsReal *Rp_sqSum)
{
	if (Rp_sqSum) *Rp_sqSum = (wsReal)size()*SQR(R_val);
	return size()*R_val;
}

wsRealCst cValVector::sum(wsVec Ra_vec, wsUint N_sz)
{
	ASSERT_SIZE(size(), N_sz, "vector-vector sum");
	return sseVsum(Ra_vec, N_sz)*R_val;
}

wsRealCst cValVector::sum(cVector &V)
{
	ASSERT_SIZE(size(), V.size(), "vector-vector sum");
	return V.sum()*R_val;
}

void cValVector::set0()
{
	halt("This command cannot be done!");
}

cVector& cValVector::inv()
{ /* Compute self-inversion */
	if (R_val == W0)
		halt("Can't take inverse to 0 in size [%d] valvector", size());
	/* Take inverse */
	R_val = W1 / R_val;
	/* Inverse of v = 1/v */
	return *this;
}

cValVector cValVector::p1p()
{ /* Compute v(1-v) */
	wsReal	R_res = R_val * (W1 - R_val);
	return cValVector(size(), R_res);
}

cStdMatrix cValVector::toMt()
{
	wsUint	N_z = size();
	wsMat	Ra_ret = sseMatrix(1, N_z);
	sseVinit(Ra_ret[0], N_z, R_val);
	return cStdMatrix(1, size(), Ra_ret);
}

/* Calculate quadratic form
	* this %*% mat %*% t(vec) */
wsReal cValVector::qf(wsMat Ra_mat, wsUint N_sz, wsReal *Ra_vec)
{
	/* column sum of mat %*% t(vec) * R_val */
	wsReal*	Ra_colSum	= sseMsumC(Ra_mat, N_sz, N_sz);
	wsReal	R_ret		= sseVV(Ra_colSum, N_sz, Ra_vec, R_val);
	sseFree(Ra_colSum);
	return R_ret;
}

/* Calculate quadratic form
	* this %*% mat %*% t(this) */
wsReal cValVector::qf(wsMat Ra_mat, wsUint N_sz)
{
	/* Just sum of Ra_mat * SQR(R_val) */
	return sseMsum(Ra_mat, N_sz) * SQR(R_val);
}

wsReal cValVector::qf(cBlkMatrix &B)
{
	if (B.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-block matrix QF", B.row(), size());
	return B.sum() * SQR(R_val);
}

wsReal cValVector::qf(cBlkMatrix &B, cVector &V)
{
	if (B.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-block matrix QF #1", B.row(), size());
	if (size() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-block matrix QF #2", size(), V.size());
	cVector VV = B.sumC();
	VV *= R_val;
	return VV.sum(V);
}

wsReal cValVector::qf(cStdMatrix &M)
{
	if (!M.sqr() || M.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix QF", size(), M.row());
	return M.sum() * SQR(R_val);
}

wsReal cValVector::qf(cStdMatrix &M, cVector &V)
{
	if (M.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix QF #1", size(), M.row());
	if (M.col() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix QF #2", M.col(), V.size());
	/* Do column sum */
	cVector VV = M.sumC();
	VV *= R_val;
	return VV.sum(V);
}

wsReal cValVector::qf(cIdtMatrix &I)
{
	if (I.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-identity matrix QF", size(), I.row());
	return I.row() * SQR(R_val);
}

wsReal cValVector::qf(cIdtMatrix &I, cVector &V)
{
	if (I.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-identity matrix QF #1", size(), I.row());
	if (V.size() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-identity matrix QF #2", size(), V.size());
	/* No need to column sum since it is same */
	/* Just summing V */
	return V.sum() * R_val;
}

wsReal cValVector::qf(cDiagMatrix &D)
{
	if (D.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-diagonal matrix QF", size(), D.row());
	return D.sum() * SQR(R_val);
}

wsReal cValVector::qf(cDiagMatrix &D, cVector &V)
{
	if (D.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-diagonal matrix QF #1", size(), D.row());
	if (size() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-diagonal matrix QF #2", size(), V.size());
	/* Just D . V */
	return sseVV(D.get()[0], size(), V.get(), R_val);
}

wsReal cValVector::qf(cSymMatrix &S)
{
	if (S.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF", size(), S.row());
	return S.sum() * SQR(R_val);
}

wsReal cValVector::qf(cSpsMatrix &P, cVector &V)
{
	if (P.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF #1", size(), P.row());
	if (size() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF #2", size(), V.size());
	/* Do column sum */
	cVector VV = P.sumC();
	VV *= R_val;
	return VV.sum(V);
}

wsReal cValVector::qf(cSpsMatrix &P)
{
	if (P.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF", size(), P.row());
	return P.sum() * SQR(R_val);
}

wsReal cValVector::qf(cSymMatrix &S, cVector &V)
{
	if (S.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF #1", size(), S.row());
	if (size() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF #2", size(), V.size());
	/* Do column sum */
	cVector VV = S.sumC();
	VV *= R_val;
	return VV.sum(V);
}

cVector& cValVector::operator=(const cVector &V)
{
	halt("This command cannot be done!");
	return *this;
}

cVector& cValVector::operator=(cVector &&V)
{
	halt("This command cannot be done!");
	return *this;
}

/* Self subtraction */
void cValVector::operator-=(wsReal R_val)
{
	halt("This command cannot be done!");
}

void cValVector::operator-=(cVector &C_vec)
{
	halt("This command cannot be done!");
}

/* Self addition */
cVector& cValVector::operator+=(wsReal R_val)
{
	halt("This command cannot be done!");
	return *this;
}

cVector& cValVector::operator+=(wsUint N_val)
{
	halt("This command cannot be done!");
	return *this;
}

cVector& cValVector::operator+=(cVector &C_vec)
{
	halt("This command cannot be done!");
	return *this;
}

/* Self production */
void cValVector::operator*=(wsRealCst R_val)
{
	halt("This command cannot be done!");
}

void cValVector::operator*=(const cVector &C_vec)
{
	halt("This command cannot be done!");
}

/* Self division */
void cValVector::operator/=(wsReal R_val)
{
	halt("This command cannot be done!");
}

void cValVector::operator/=(cVector &C_vec)
{
	halt("This command cannot be done!");
}

/* Subtract and return */
cValVector cValVector::operator-(wsRealCst R_inpVal)
{
	return cValVector(size(), R_val - R_inpVal);
}

cVector cValVector::operator-(const cVector &V)
{
	wsReal *Ra_ret = sseVector(size());
	wsReal *p = V.get();

	for (wsUint i=0 ; i<size() ; i++)
		Ra_ret[i] = R_val - p[i];

	return cVector(Ra_ret, size());
}

/* Addition and return */
cVector cValVector::operator+(wsRealCst R_inpVal)
{
	return cValVector(size(), R_val + R_inpVal);
}

cVector cValVector::operator+(const cVector &V)
{
	wsReal *Ra_ret = sseVector(size());
	wsReal *p = V.get();

	for (wsUint i=0 ; i<size() ; i++)
		Ra_ret[i] = R_val + p[i];

	return cVector(Ra_ret, size());
}

/* Production and return */
cVector cValVector::operator*(wsRealCst R_inpVal)
{
	return cValVector(size(), R_val*R_inpVal);
}

cVector cValVector::operator*(cVector &V)
{
	return V.clone() * R_val;
}

cVector cValVector::operator*(cStdMatrix &M)
{
	/* Dim check */
	if (size()!=M.row())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix multiplication", size(), M.row());
	/* Just do column sum */
	return M.sumC() * R_val;
}

cVector cValVector::operator*(cSymMatrix &S)
{
	/* Dim check */
	if (size()!=S.row())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix multiplication", size(), S.row());
	/* Just do column sum */
	return S.sumC() * R_val;
}

cVector cValVector::Mt(cStdMatrix &M)
{
	/* Dim check */
	if (size()!=M.col())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix multiplication", size(), M.col());
	/* Just do row sum */
	return M.sumR() * R_val;
}

cVector& cValVector::operator*=(cDiagMatrix &D)
{
	halt("This command cannot be done!");
	return *this;
}

cVector& cValVector::operator*=(cIdtMatrix &I)
{
	halt("This command cannot be done!");
	return *this;
}

cVector cValVector::operator*(cDiagMatrix &D)
{
	/* Dim check */
	if (size()!=D.row())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix multiplication", size(), D.row());
	/* Just copy and return as vector */
	wsReal *Ra_ret = sseVector(D.row());
	sseVpC(D.get()[0], R_val, Ra_ret, D.row());
	return cVector(Ra_ret, D.row());
}

cVector cValVector::operator*(const cIdtMatrix &I)
{
	return cValVector(size(), R_val);
}

cVector cValVector::sq()
{
	return cValVector(size(), SQR(R_val));
}

cVector cValVector::sqrt()
{
	return cValVector(size(), ::sqrt(R_val));
}

cVector cValVector::cb()
{
	return cValVector(size(), CUBE(R_val));
}

/* Division and return */
cValVector cValVector::operator/(wsRealCst R_inpVal)
{
	return cValVector(size(), R_val/R_inpVal);
}

cVector cValVector::operator/(cVector &V)
{
	return V.inv() * R_val;
}

/* Special operations */
void cValVector::file(wsStrCst S_ext)
{
	exportValVector(S_ext, R_val, size());
}

cVector& cValVector::sub(cVector &V_s, cVector &V_t)
{
	halt("This command cannot be done!");
	return *this;
}

wsRealCst cValVector::ss()
{
	/* 1'1 = size() */
	return SQR(R_val) * size();
}

cValVector cValVector::clone()
{
	return cValVector(size(), R_val);
}

cValVector cValVector::subset(wsUint N_start, wsUint N_sz)
{
	wsUint N_z = size();
	/* size check */
	if (N_start >= N_z || (N_start+N_sz) > N_z)
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector subsetting", N_z, N_start+N_sz);
	return cValVector(N_sz, R_val);
}

cSymMatrix cValVector::tV(wsReal R_prod)
{
	wsUint	N_z		= size();
	wsReal	R_v		= NA(R_prod) ? SQR(R_val) : SQR(R_val)*R_prod; 
	wsSym	Ra_ret	= sseSymMat(N_z);
	/* Init with value^2 or value^2 * R_prod */
	for (wsUint i=0 ; i<N_z ; i++)
		sseVinit(Ra_ret[i], i+1, R_v);
	return cSymMatrix(Ra_ret, N_z);
}

cStdMatrix cValVector::tV(cVector &V, wsReal R_prod)
{
	wsUint	N_z		= size();
	wsUint	N_vz	= V.size();
	/* Unit vector %*% normal vector = row is each normal vector */
	wsMat	Ra_ret	= sseMatrix(N_z, N_vz);
	cVector	V_tmp	= V * R_prod * R_val;
	wsReal*	Ra_tmp	= V_tmp.get();
	for (wsUint i=0 ; i<N_z ; i++)
		memcpy(Ra_ret[i], Ra_tmp, sizeof(wsReal)*N_vz);

	return cStdMatrix(N_z, N_vz, Ra_ret);
}

cVector& cValVector::neg(wsReal R_add)
{
	halt("This command cannot be done!");
	return *this;
}

/*
 *
 * cUnitVector definition
 * 
 */

cUnitVector::cUnitVector(wsUint N_inpSz)
{
	N_sz = N_inpSz;
	Ra_buf = NULL;
}

cUnitVector::cUnitVector(cUnitVector&& V_noname)
{
	if (V_noname.isDdealloc()) setDontDealloc();

	V_noname.setDontDealloc();
	N_sz = V_noname.size();
	Ra_buf = NULL;
}

void cUnitVector::init(wsUint N_inpSz, wsReal *Ra_val, wsReal *Rp_val,
	char B_rand, char B_inpDDealloc)
{
	if (B_rand)
		halt("Unit vector can't be initiated with random");
	else if (Rp_val)
		halt("Unit vector can't be initiated with scalar value");
	else if (Ra_val)
		halt("Unit vector can't be initiated with vector values");

	N_sz	= N_inpSz;
	Ra_buf	= NULL;
}

wsReal* cUnitVector::get() const
{
	/* Return nothing, since there is nothing */
	return NULL;
}

wsRealCst cUnitVector::sum(wsReal *Rp_sqSum)
{
	/* Squared sum of 1'1 = size */
	if (Rp_sqSum) *Rp_sqSum = (wsReal)size();
	return size();
}

wsRealCst cUnitVector::sum(wsVec Ra_vec, wsUint N_sz)
{
	ASSERT_SIZE(size(), N_sz, "vector-vector sum");
	/* 1'V = just sum of V */
	return sseVsum(Ra_vec, N_sz);
}

wsRealCst cUnitVector::sum(cVector &V)
{
	ASSERT_SIZE(size(), V.size(), "vector-vector sum");
	/* 1'V = just sum of V */
	return V.sum();
}

void cUnitVector::set0()
{
	/* Can't set with 0 since the value cannot be changed */
	halt("This command cannot be done!");
}

cVector& cUnitVector::inv()
{ /* Compute self-inversion */
	/* Inverse of unit = unit */
	return *this;
}

cVector cUnitVector::p1p()
{ /* Compute v(1-v) */
	wsUint	N_z = size();
	/* 1 * 0 == 0 */
	wsReal*	Ra_ret = sseVector(N_z);
	sseVinit(Ra_ret, N_z, W0);
	return cVector(Ra_ret, N_z);
}

cStdMatrix cUnitVector::toMt()
{
	wsUint	N_z = size();
	wsMat	Ra_ret = sseMatrix(1, N_z);
	sseVinit(Ra_ret[0], N_z, W1);
	return cStdMatrix(1, size(), Ra_ret);
}

/* Calculate quadratic form
 * this %*% mat %*% t(vec) */
wsReal cUnitVector::qf(wsMat Ra_mat, wsUint N_sz, wsReal *Ra_vec)
{
	/* column sum of mat %*% t(vec) */
	wsReal*	Ra_colSum	= sseMsumC(Ra_mat, N_sz, N_sz);
	wsReal	R_ret		= sseVV(Ra_colSum, N_sz, Ra_vec);
	sseFree(Ra_colSum);
	return R_ret;
}

/* Calculate quadratic form
 * this %*% mat %*% t(this) */
wsReal cUnitVector::qf(wsMat Ra_mat, wsUint N_sz)
{
	/* Just sum of Ra_mat */
	return sseMsum(Ra_mat, N_sz);
}

wsReal cUnitVector::qf(cBlkMatrix &B)
{
	if (B.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-block matrix QF", B.row(), size());
	return B.sum();
}

wsReal cUnitVector::qf(cBlkMatrix &B, cVector &V)
{
	if (B.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-block matrix QF #1", B.row(), size());
	if (size() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-block matrix QF #2", size(), V.size());
	cVector VV = B.sumC();
	return VV.sum(V);
}

wsReal cUnitVector::qf(cStdMatrix &M)
{
	if (!M.sqr() || M.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix QF", size(), M.row());
	return M.sum();
}

wsReal cUnitVector::qf(cStdMatrix &M, cVector &V)
{
	if (M.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix QF #1", size(), M.row());
	if (M.col() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix QF #2", M.col(), V.size());
	/* Do column sum */
	cVector VV = M.sumC();
	return VV.sum(V);
}

wsReal cUnitVector::qf(cIdtMatrix &I)
{
	if (I.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-identity matrix QF", size(), I.row());
	return I.row();
}

wsReal cUnitVector::qf(cIdtMatrix &I, cVector &V)
{
	if (I.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-identity matrix QF #1", size(), I.row());
	if (V.size() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-identity matrix QF #2", size(), V.size());
	/* No need to column sum since it is same */
	/* Just summing V */
	return V.sum();
}

wsReal cUnitVector::qf(cDiagMatrix &D)
{
	if (D.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-diagonal matrix QF", size(), D.row());
	return D.sum();
}

wsReal cUnitVector::qf(cDiagMatrix &D, cVector &V)
{
	if (D.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-diagonal matrix QF #1", size(), D.row());
	if (size() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-diagonal matrix QF #2", size(), V.size());
	/* Just D . V */
	return sseVV(D.get()[0], size(), V.get());
}

wsReal cUnitVector::qf(cSymMatrix &S)
{
	if (S.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF", size(), S.row());
	return S.sum();
}

wsReal cUnitVector::qf(cSymMatrix &S, cVector &V)
{
	if (S.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF #1", size(), S.row());
	if (size() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF #2", size(), V.size());
	/* Do column sum */
	cVector VV = S.sumC();
	return VV.sum(V);
}

wsReal cUnitVector::qf(cSpsMatrix &S)
{
	if (S.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF", size(), S.row());
	return S.sum();
}

wsReal cUnitVector::qf(cSpsMatrix &S, cVector &V)
{
	if (S.row() != size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF #1", size(), S.row());
	if (size() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-symmetric matrix QF #2", size(), V.size());
	/* Do column sum */
	cVector VV = S.sumC();
	return VV.sum(V);
}

cVector& cUnitVector::operator=(const cVector &V)
{
	halt("This command cannot be done!");
	return *this;
}

cVector& cUnitVector::operator=(cVector &&V)
{
	halt("This command cannot be done!");
	return *this;
}

/* Self subtraction */
void cUnitVector::operator-=(wsReal R_val)
{
	halt("This command cannot be done!");
}

void cUnitVector::operator-=(cVector &C_vec)
{
	halt("This command cannot be done!");
}

/* Self addition */
cVector& cUnitVector::operator+=(wsReal R_val)
{
	halt("This command cannot be done!");
	return *this;
}

cVector& cUnitVector::operator+=(wsUint N_val)
{
	halt("This command cannot be done!");
	return *this;
}

cVector& cUnitVector::operator+=(cVector &C_vec)
{
	halt("This command cannot be done!");
	return *this;
}

/* Self production */
void cUnitVector::operator*=(wsRealCst R_val)
{
	halt("This command cannot be done!");
}

void cUnitVector::operator*=(const cVector &C_vec)
{
	halt("This command cannot be done!");
}

/* Self division */
void cUnitVector::operator/=(wsReal R_val)
{
	halt("This command cannot be done!");
}

void cUnitVector::operator/=(cVector &C_vec)
{
	halt("This command cannot be done!");
}

/* Subtract and return */
cValVector cUnitVector::operator-(wsRealCst R_val)
{
	return cValVector(size(), W1 - R_val);
}

cVector cUnitVector::operator-(const cVector &C_vec)
{
	wsReal *Ra_ret = sseVector(size());
	wsReal *p = C_vec.get();

	for (wsUint i=0 ; i<size() ; i++)
		Ra_ret[i] = W1 - p[i];

	return cVector(Ra_ret, size());
}

/* Addition and return */
cVector cUnitVector::operator+(wsRealCst R_val)
{
	return cValVector(size(), W1 + R_val);
}

cVector cUnitVector::operator+(const cVector &C_vec)
{
	wsReal *Ra_ret = sseVector(size());
	wsReal *p = C_vec.get();

	for (wsUint i=0 ; i<size() ; i++)
		Ra_ret[i] = W1 + p[i];

	return cVector(Ra_ret, size());
}

/* Production and return */
cVector cUnitVector::operator*(wsRealCst R_val)
{
	return cValVector(size(), R_val);
}

cVector cUnitVector::operator*(cVector &C_vec)
{
	return C_vec.clone();
}

cVector cUnitVector::operator*(cStdMatrix &M)
{
	/* Dim check */
	if (size()!=M.row())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix multiplication", size(), M.row());
	/* Just do column sum */
	return M.sumC();
}

cVector cUnitVector::operator*(cSymMatrix &S)
{
	/* Dim check */
	if (size()!=S.row())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix multiplication", size(), S.row());
	/* Just do column sum */
	return S.sumC();
}

cVector cUnitVector::Mt(cStdMatrix &M)
{
	/* Dim check */
	if (size()!=M.col())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix multiplication", size(), M.col());
	/* Just do row sum */
	return M.sumR();
}

cVector& cUnitVector::operator*=(cDiagMatrix &D)
{
	halt("This command cannot be done!");
	return *this;
}

cVector& cUnitVector::operator*=(cIdtMatrix &I)
{
	/* No effective work */
	return *this;
}

cVector cUnitVector::operator*(cDiagMatrix &D)
{
	/* Dim check */
	if (size()!=D.row())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-matrix multiplication", size(), D.row());
	/* Just copy and return as vector */
	wsReal *Ra_ret = sseVector(D.row());
	memcpy(Ra_ret, D.get()[0], sizeof(wsReal)*D.row());
	return cVector(Ra_ret, D.row());
}

cVector cUnitVector::operator*(const cIdtMatrix &I)
{
	return cUnitVector(size());
}

cVector cUnitVector::sq()
{
	/* 1^2 = 1 */
	return cUnitVector(size());
}

cVector cUnitVector::sqrt()
{
	/* sqrt(1) = 1 */
	return cUnitVector(size());
}

cVector cUnitVector::cb()
{
	/* 1^3 = 1 */
	return cUnitVector(size());
}

/* Division and return */
cValVector cUnitVector::operator/(wsRealCst R_val)
{
	return cValVector(size(), W1/R_val);
}

cVector cUnitVector::operator/(cVector &C_vec)
{
	wsUint N = size();
	/* Dim check */
	if (N != C_vec.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector-vector division", size(), C_vec.size());
	wsReal *R = sseVector(N);
	sseVdV(get(), C_vec.get(), R, N);
	return cVector(R, N);
}

/* Special operations */
void cUnitVector::file(wsStrCst S_ext)
{
	/* Export unit vector */
	exportUnitVector(S_ext, size());
}

cVector& cUnitVector::sub(cVector &V_s, cVector &V_t)
{
	halt("This command cannot be done!");
	return *this;
}

wsRealCst cUnitVector::ss()
{
	/* 1'1 = size() */
	return size();
}

cVector cUnitVector::clone()
{
	return cUnitVector(size());
}

cVector cUnitVector::subset(wsUint N_start, wsUint N_sz)
{
	wsUint N_z = size();
	/* size check */
	if (N_start >= N_z || (N_start+N_sz) > N_z)
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "vector subsetting", N_z, N_start+N_sz);
	return cUnitVector(N_sz);
}

cSymMatrix cUnitVector::tV(wsReal R_prod)
{
	wsUint	N_z		= size();
	wsSym	Ra_ret	= sseSymMat(N_z);
	for (wsUint i=0 ; i<N_z ; i++)
		sseVinit(Ra_ret[i], i+1, W1);
	return cSymMatrix(Ra_ret, N_z);
}

cStdMatrix cUnitVector::tV(cVector &V, wsReal R_prod)
{
	wsUint	N_z		= size();
	wsUint	N_vz	= V.size();
	/* Unit vector %*% normal vector = row is each normal vector */
	wsMat	Ra_ret	= sseMatrix(N_z, N_vz);
	cVector	V_tmp	= V * R_prod;
	wsReal*	Ra_tmp	= V_tmp.get();
	for (wsUint i=0 ; i<N_z ; i++)
		memcpy(Ra_ret[i], Ra_tmp, sizeof(wsReal)*N_vz);

	return cStdMatrix(N_z, N_vz, Ra_ret);
}

cVector& cUnitVector::neg(wsReal R_add)
{
	halt("This command cannot be done!");
	return *this;
}



// Do vt = v1+v2 or vt = (v1+v2)*m
void sseVaV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint	N_med = 0;
#ifdef USE_SSE
	N_med = getMed(N_sz);
#endif

	if (NA(R_mul)) {
#ifdef USE_SSE
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseAdd(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v1[i] + Ra_v2[i];
	} else {
#ifdef USE_SSE
		sse_t sse_m = sseSet(R_mul);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseMul(sseAdd(*sse_v1, *sse_v2), sse_m);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = (Ra_v1[i] + Ra_v2[i]) * R_mul;
	}
}

/* Do tgt = v1 + S when tgt != NULL
 *    v1 += S      when tgt == NULL */
//void sseVectorAdd(wsUint N_sz, wsReal *Ra_v1, wsReal R_val,
// 	wsReal *Ra_tgt=NULL)
// {
// 	wsUint	i;
// 	wsUint	N_med = 0;
// #ifdef USE_SSE
// 	sse_t	sse_S	= sseSet(R_val);
// 	N_med	= getMed(N_sz);
// #endif
// 
// 	if (Ra_tgt != NULL) {
// #ifdef USE_SSE
// 		for (i=0 ; i<N_med ; i+=sseJmp) {
// 			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
// 			sse_t *sse_vt	= (sse_t *)(Ra_tgt + i);
// 
// 			*sse_vt = sseAdd(*sse_v1, sse_S);
// 		}
// #endif
// 		/* Do rest part */
// 		for (i=N_med ; i<N_sz ; i++) Ra_tgt[i] = Ra_v1[i] + R_val;
// 	} else {
// #ifdef USE_SSE
// 		for (i=0 ; i<N_med ; i+=sseJmp) {
// 			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
// 			*sse_v1 = sseAdd(*sse_v1, sse_S);
// 		}
// #endif
// 		/* Do rest part */
// 		for (i=N_med ; i<N_sz ; i++) Ra_v1[i] += R_val;
// 	}
// }

/* Do tgt = v1 - v2 when tgt != NULL
 *    v1 -= v2      when tgt == NULL */
void sseVectorSub(wsUint N_sz, wsReal *Ra_v1, wsReal *Ra_v2,
	wsReal *Ra_tgt=NULL)
{
	wsUint	i;
	wsUint	N_med = 0;
#ifdef USE_SSE
	N_med	= getMed(N_sz);
#endif

	if (Ra_tgt != NULL) {
#ifdef USE_SSE
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2	= (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt	= (sse_t *)(Ra_tgt + i);

			*sse_vt = sseSub(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++) Ra_tgt[i] = Ra_v1[i] - Ra_v2[i];
	} else {
#ifdef USE_SSE
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2	= (sse_t *)(Ra_v2 + i);

			*sse_v1 = sseSub(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++) Ra_v1[i] -= Ra_v2[i];
	}
}

/* Do tgt = v1 - S when tgt != NULL
 *    v1 -= S      when tgt == NULL */
void sseVectorSub(wsUint N_sz, wsReal *Ra_v1, wsReal R_val,
	wsReal *Ra_tgt=NULL)
{
	wsUint	i;
	wsUint	N_med = 0;
#ifdef USE_SSE
	sse_t	sse_S	= sseSet(R_val);
	N_med	= getMed(N_sz);
#endif

	if (Ra_tgt != NULL) {
#ifdef USE_SSE
		N_med	= getMed(N_sz);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			sse_t *sse_vt	= (sse_t *)(Ra_tgt + i);

			*sse_vt = sseSub(*sse_v1, sse_S);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_tgt[i] = Ra_v1[i] - R_val;
	} else {
		wsUint	i;
		wsUint	N_med = 0;
#ifdef USE_SSE
		N_med	= getMed(N_sz);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);

			*sse_v1 = sseSub(*sse_v1, sse_S);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_v1[i] -= R_val;
	}
}


/* Do tgt = v1 * v2 when tgt != NULL
 *    v1 *= v2      when tgt == NULL */
void sseVectorMul(wsUint N_sz, wsReal *Ra_v1, wsReal *Ra_v2,
	wsReal *Ra_tgt=NULL)
{
	wsUint	i;
	wsUint	N_med = 0;
#ifdef USE_SSE
	N_med	= getMed(N_sz);
#endif

	if (Ra_tgt != NULL) {
#ifdef USE_SSE
		N_med	= getMed(N_sz);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2	= (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt	= (sse_t *)(Ra_tgt + i);

			*sse_vt = sseMul(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_tgt[i] = Ra_v1[i] * Ra_v2[i];
	} else {
		wsUint	i;
		wsUint	N_med = 0;
#ifdef USE_SSE
		N_med	= getMed(N_sz);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2	= (sse_t *)(Ra_v2 + i);

			*sse_v1 = sseMul(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_v1[i] *= Ra_v2[i];
	}
}

/* Do tgt = v1 * S when tgt != NULL
 *    v1 *= S      when tgt == NULL */
void sseVectorMul(wsUint N_sz, wsReal *Ra_v1, wsReal R_val,
	wsReal *Ra_tgt=NULL)
{
	wsUint	i;
	wsUint	N_med = 0;
#ifdef USE_SSE
	sse_t	sse_S	= sseSet(R_val);
	N_med	= getMed(N_sz);
#endif

	if (Ra_tgt != NULL) {
#ifdef USE_SSE
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			sse_t *sse_vt	= (sse_t *)(Ra_tgt + i);

			*sse_vt = sseMul(*sse_v1, sse_S);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_tgt[i] = Ra_v1[i] * R_val;
	} else {
#ifdef USE_SSE
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			*sse_v1 = sseMul(*sse_v1, sse_S);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_v1[i] *= R_val;
	}
}

/* Do tgt = v1 / v2 when tgt != NULL
 *    v1 /= v2      when tgt == NULL */
void sseVectorDiv(wsUint N_sz, wsReal *Ra_v1, wsReal *Ra_v2,
	wsReal *Ra_tgt=NULL)
{
	wsUint	i;
	wsUint	N_med = 0;
#ifdef USE_SSE
	N_med	= getMed(N_sz);
#endif

	if (Ra_tgt != NULL) {
#ifdef USE_SSE
		N_med	= getMed(N_sz);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2	= (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt	= (sse_t *)(Ra_tgt + i);

			*sse_vt = sseDiv(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_tgt[i] = Ra_v1[i] / Ra_v2[i];
	} else {
		wsUint	i;
		wsUint	N_med = 0;
#ifdef USE_SSE
		N_med	= getMed(N_sz);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2	= (sse_t *)(Ra_v2 + i);

			*sse_v1 = sseDiv(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_v1[i] /= Ra_v2[i];
	}
}

/* Do tgt = v1 / S when tgt != NULL
 *    v1 /= S      when tgt == NULL */
void sseVectorDiv(wsUint N_sz, wsReal *Ra_v1, wsReal R_val,
	wsReal *Ra_tgt=NULL)
{
	wsUint	i;
	wsUint	N_med = 0;
#ifdef USE_SSE
//	sse_t	sse_S	= sseSet(R_val);
	sse_t	sse_S	= sseSet(W1/R_val);
	N_med	= getMed(N_sz);
#endif

	if (Ra_tgt != NULL) {
#ifdef USE_SSE
		N_med	= getMed(N_sz);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);
			sse_t *sse_vt	= (sse_t *)(Ra_tgt + i);

//			*sse_vt = sseDiv(*sse_v1, sse_S);
			*sse_vt = sseMul(*sse_v1, sse_S);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_tgt[i] = Ra_v1[i] / R_val;
	} else {
		wsUint	i;
		wsUint	N_med = 0;
#ifdef USE_SSE
		N_med	= getMed(N_sz);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1	= (sse_t *)(Ra_v1 + i);

//			*sse_v1 = sseDiv(*sse_v1, sse_S);
			*sse_v1 = sseMul(*sse_v1, sse_S);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_v1[i] /= R_val;
	}
}

/*
 *
 * cVector definition
 * 
 */

cVector::cVector()
{
	B_dontDealloc = 0;
	Ra_buf = NULL;
	N_sz = 0;
}

cVector::cVector(cVector &&V_noname)
{
	B_dontDealloc = V_noname.isDdealloc();
	Ra_buf = V_noname.get();
	N_sz = V_noname.size();

	V_noname.setDontDealloc();
}

cVector::cVector(cVector &C_vec, wsUint N_rep)
{
	Ra_buf = NULL;
	B_dontDealloc = 0;
	N_sz = C_vec.size()*N_rep;
	/* In default, set to 0 */
	sseMalloc(Ra_buf, wsReal, N_sz);

	if (N_rep == 1)
		memcpy(Ra_buf, C_vec.get(), sizeof(wsReal)*N_sz);
	else {
		wsReal *R = C_vec.get();
		for (wsUint i=0,_i=0 ; i<C_vec.size() ; i++) {
			for (wsUint j=0 ; j<N_rep ; j++)
				Ra_buf[_i++] = R[i];
		}
	}
}

// cVector::cVector(cVector &C_vec)
// {
// 	B_dontDealloc = 0;
// 	Ra_buf = C_vec.get();
// 	N_sz = C_vec.size();
//}

cVector::cVector(wsReal *Rp_inpVal, wsUint N_inpSz, char B_inpDDealloc/*=0*/)
{
	Ra_buf = NULL;
	init(N_inpSz, Rp_inpVal, NULL, 0, B_inpDDealloc);
}

cVector::cVector(wsReal R_val, wsUint N_inpSz)
{
	Ra_buf = NULL;
	init(N_inpSz, NULL, &R_val, 0, 0);
}

cVector::cVector(wsUint N_inpSz, char B_rand/*=0*/)
{
	Ra_buf = NULL;
	init(N_inpSz, NULL, NULL, B_rand);
}

cVector::cVector(xRealSort *Xa_rs, wsUint N_inpSz)
{
	Ra_buf = NULL;
	wsReal *Ra_v = sseVector(N_inpSz);
	for (wsUint i=0 ; N_inpSz ; i++)
		Ra_v[i] = Xa_rs[i].V;
	init(N_inpSz, Ra_v, NULL, 0);
}

cVector::~cVector()
{
	if (!B_dontDealloc && Ra_buf) {
		if (cVector::B_report)
			LOG("Memory deallocate [%x]\n", Ra_buf);
		sseFree(Ra_buf);
	}
	Ra_buf = NULL;
}


void cVector::init(wsUint N_inpSz, wsReal *Ra_val/*=NULL*/,
	wsReal *Rp_valInit/*=NULL*/, char B_rand/*=0*/, char B_inpDDealloc/*=0*/)
{
	if (Ra_buf)
		sseFree(Ra_buf);
	B_dontDealloc = B_inpDDealloc;
	N_sz = N_inpSz;

	if (Ra_val) {
		/* Maintain original contents */
		Ra_buf	= Ra_val;
		return;
	}

	/* In default, set to 0 */
	if (cVector::B_report)
		LOG("Vector [%d] elements assigned\n", N_sz);
	sseCalloc(Ra_buf, wsReal, N_sz);
	if (Rp_valInit)
		sseVinit(Ra_buf, N_sz, *Rp_valInit);
// 		for (wsUint i=0 ; i<N_sz ; i++)
// 			Ra_buf[i] = *Rp_val;
	else if (B_rand)
		for (wsUint i=0 ; i<N_sz ; i++)
			Ra_buf[i] = (wsReal)(rand()%1000) - REAL_CONST(500.0);
}

wsReal* cVector::get() const
{
	return Ra_buf;
}

wsUintCst cVector::size() const
{
	return N_sz;
}

wsRealCst cVector::sum(wsReal *Rp_sqSum/*=NULL*/)
{
	return sseVsum(get(), size(), Rp_sqSum);
}

wsRealCst cVector::asum(wsReal *Rp_sqSum/*=NULL*/)
{
	return sseVasum(get(), size(), Rp_sqSum);
}

wsRealCst cVector::aasum(cVector& Vsub)
{
	wsVecCst	Ra_v1	= get();
	wsVecCst	Ra_v2	= Vsub.get();
	wsUint	N_med	= 0;
	wsReal	R_sum	= W0;

#ifdef USE_SSE
	sse_t sse_sum = sseSet(0.0);
	N_med = getMed(N_sz);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
		sse_t sse_t1 = sseAbs(*sse_v1);
		sse_t sse_t2 = sseAbs(*sse_v2);

		sse_sum = sseAdd(sse_sum, sseAbs(sseSub(sse_t1, sse_t2)));
	}
	sseSum(sse_sum, R_sum);
#endif
	/* Do rest part */
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_sum += fabs(fabs(Ra_v1[i]) - fabs(Ra_v2[i]));

	return R_sum;
}

wsRealCst cVector::asum(cVector& Vsub)
{
	wsVecCst	Ra_v1	= get();
	wsVecCst	Ra_v2	= Vsub.get();
	wsUint	N_med	= 0;
	wsReal	R_sum	= W0;

#ifdef USE_SSE
	sse_t sse_sum = sseSet(0.0);
	N_med = getMed(N_sz);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);

		sse_sum = sseAdd(sse_sum, sseAbs(sseSub(*sse_v1, *sse_v2)));
	}
	sseSum(sse_sum, R_sum);
#endif
	/* Do rest part */
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_sum += fabs(Ra_v1[i] - Ra_v2[i]);

	return R_sum;
}

wsRealCst cVector::ssum(cVector& Vsub)
{
	wsVecCst	Ra_v1	= get();
	wsVecCst	Ra_v2	= Vsub.get();
	wsUint	N_med	= 0;
	wsReal	R_sum	= W0;

#ifdef USE_SSE
	sse_t sse_sum = sseSet(0.0);
	N_med = getMed(N_sz);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);

		sse_sum = sseAdd(sse_sum, sseSub(*sse_v1, *sse_v2));
	}
	sseSum(sse_sum, R_sum);
#endif
	/* Do rest part */
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_sum += Ra_v1[i] - Ra_v2[i];

	return R_sum;
}

wsRealCst cVector::sum(wsVec Ra_vec, wsUint N_sz)
{
	ASSERT_SIZE(size(), N_sz, "vector-vector sum");
	return sseVV(get(), size(), Ra_vec);
}

wsRealCst cVector::sum(cVector &C_vec)
{
	return sum(C_vec.get(), C_vec.size());
}

wsRealCst cVector::mean()
{
	return sum() / (wsReal)size();
}

cStdMatrix cVector::colP(cStdMatrix& M)
{
	wsUint	N = size();
	wsUint	C = M.col();
	ASSERT_SIZE(N, C, "Vector-matrix column-wise product");
	wsMat	Mp = M.get();
	wsMat	R = sseMatrix(M.row(), C);
	wsVec	V = get();
	LOOP (i, M.row()) sseVpC(Mp[i], V[i], R[i], C);
	return cStdMatrix(M.row(), C, R);
}

wsRealCst cVector::qf(wsMat Ra_mat, wsUint N_sz, wsVec Ra_vec)
{
	if (N_sz != this->N_sz)
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), N_sz);
	return sseVpMpV(get(), Ra_mat, N_sz, Ra_vec);
}

wsReal cVector::qf(wsReal **Ra_mat, wsUint N_sz)
{
	return qf(Ra_mat, N_sz, get());
}

wsReal cVector::qf(cStdMatrix &M)
{
	if (!M.sqr() || size()!=M.row())
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), M.row());
	return qf(M.get(), size());
}

wsReal cVector::qf(cStdMatrix &M, cVector &V)
{
	if (!M.sqr() || size()!=M.row())
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), M.row());
	return qf(M.get(), size(), V.get());
}

wsReal cVector::qf(cIdtMatrix &I)
{
	if (!I.sqr() || size()!=I.row())
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), I.row());
	return sseVV(get(), size());
	//return qf(M.get(), size());
}

wsReal cVector::qf(cIdtMatrix &I, cVector &V)
{
	if (!I.sqr() || size()!=I.row())
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), I.row());
	return sseVV(get(), size(), V.get());
	//return qf(M.get(), size(), V.get());
}

wsReal cVector::qf(cDiagMatrix &D)
{
	if (!D.sqr() || size()!=D.row()) /* Dimension check */
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), D.row());
	wsReal *Ra_D = D.get()[0];
	return sseVVV(get(), Ra_D, get(), size());
	//return qf(M.get(), size());
}

wsReal cVector::qf(cDiagMatrix &D, cVector &V)
{
	if (!D.sqr() || size()!=D.row()) /* Dimension check */
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), D.row());
	wsReal *Ra_D = D.get()[0];
	return sseVVV(get(), Ra_D, V.get(), size());
	//return qf(M.get(), size(), V.get());
}

wsReal cVector::qf(cSymMatrix &S)
{
	if (size()!=S.row()) /* Dimension check */
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), S.row());
	return sseVpSpV(get(), S.get(), size());
}

wsReal cVector::qf(cSymMatrix &S, cVector &V)
{
	if (size()!=S.row()) /* Dimension check */
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), S.row());
	return sseVpSpV(get(), S.get(), size(), V.get());
}

wsReal cVector::qf(cSpsMatrix &P)
{
	if (size()!=P.row() || size()!=P.col()) /* Dimension check */
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), P.row());
	return sseVpPpV(get(), size(), P.get()[0], P.col(), P.getRowEnds(),
		P.getColIdxs());
}

wsReal cVector::qf(cBlkMatrix &B)
{
	wsReal	R_ret	= W0;
	if (B.rep()!=size() && B.row()!=size()) /* Dimension check */
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), B.row());

	wsReal R_sumB	= B.sum();
	if (B.rep() == size()) {
		for (wsUint i=0,N=size() ; i<N ; i++)
			R_ret += SQR(Ra_buf[i])*R_sumB;
		R_ret /= (wsReal)B.rep();
	} else {
		wsUint N_b = B.blk();
		wsReal *Rp_T = NULL;
		sseCalloc(Rp_T, wsReal, size());

		for (wsUint i=0,N=B.rep() ; i<N ; i++) {
			multMpM(&Ra_buf, 1, size()-((N-i-1)*N_b), 0, i*N_b,
				B.get(), N_b, N_b, 0, 0,
				&Rp_T, 1, size(), 0, i*N_b);
		}

		/* Do Rp_T %*% this */
		R_ret = sseVV(Rp_T, size(), get());
		sseFree(Rp_T);
	}

	return R_ret;
}


wsReal cVector::qf(cBlkMatrix &B, cVector &V)
{
	wsReal R_ret	= W0;
	if (B.rep()!=size() && B.row()!=size())
		halt_fmt(WISARD_SYST_INVL_DIM, "Quadratic form", size(), B.row());

	if (B.rep() == size()) {
		wsReal R_sumB	= B.sum();
		for (wsUint i=0,N=size() ; i<N ; i++)
			R_ret += V.get()[i]*Ra_buf[i]*R_sumB;
		R_ret /= (wsReal)B.rep();
	} else {
		wsUint N_b = B.blk();
		wsReal *Rp_T = NULL;
		sseCalloc(Rp_T, wsReal, size());

		for (wsUint i=0,N=B.col() ; i<N ; i+=N_b) {
			multMpM(&Ra_buf, 1, size()-(N-i-N_b), 0, i,
				B.get(), N_b, N_b, 0, 0,
				&Rp_T, 1, size(), 0, i);
		}

		/* Do Rp_T %*% this */
		R_ret = sseVV(Rp_T, size(), V.get());
		sseFree(Rp_T);
	}

	return R_ret;
}

cVector& cVector::operator=(const cVector &V)
{
	wsReal *Ra_ret = sseVector(V.size());
	memcpy(Ra_ret, V.get(), sizeof(wsReal)*V.size());
	init(V.size(), Ra_ret);
	return *this;
}

cVector& cVector::operator=(cVector &&V)
{
	init(V.size(), V.get());
	V.setDontDealloc();
	return *this;
}

void cVector::operator-=(wsReal R_val)
{
	sseVectorSub(size(), get(), R_val);
}

void cVector::operator-=(cVector &V)
{
	/* Sanity check : Length */
	if (V.size() != size())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector subtraction", size(), V.size());

	sseVectorSub(size(), get(), V.get());
}

cVector& cVector::operator+=(wsReal R_val)
{
	sseVaC(get(), R_val, get(), size());
	return *this;
}

cVector& cVector::operator+=(wsUint N_val)
{
	sseVaC(get(), (wsReal)N_val, get(), size());
	return *this;
}

cVector& cVector::operator+=(cVector &V)
{
	/* Sanity check : Length */
	if (V.size() != size())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector addition", size(), V.size());

	sseVaV(get(), V.get(), get(), size());
	return *this;
}

void cVector::operator*=(wsRealCst R_val)
{
	sseVectorMul(size(), get(), R_val);
}

void cVector::operator*=(const cVector &V)
{
	/* Sanity check : Length */
	if (V.size() != size())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector multiplication", size(), V.size());

	sseVectorMul(size(), get(), V.get());
}

void cVector::operator/=(wsReal R_val)
{
	sseVectorDiv(size(), get(), R_val);
}

void cVector::operator/=(cVector &V)
{
	/* Sanity check : Length */
	if (V.size() != size())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector division", size(), V.size());

	sseVectorDiv(size(), get(), V.get());
}

cVector cVector::operator-(wsRealCst R_val)
{
	wsReal *Ra_ret = sseVector(size());

	sseVectorSub(size(), get(), R_val, Ra_ret);
	return cVector(Ra_ret, size());
}

cVector cVector::operator-(const cVector &V)
{
	/* Sanity check : Length */
	if (V.size() != size())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector subtraction", size(), V.size());

	wsReal *Ra_ret = sseVector(size());

	sseVectorSub(size(), get(), V.get(), Ra_ret);
	return cVector(Ra_ret, size());
}

cVector cVector::operator+(wsRealCst R_val)
{
	wsVec Ra_ret = sseVector(size());
	sseVaC(get(), R_val, Ra_ret, size());
	return cVector(Ra_ret, size());
}

cVector cVector::operator+(const cVector &V)
{
	ASSERT_SIZE(size(), V.size(), "Vector addition");
	wsVec Ra_ret = sseVector(size());
	sseVaV(get(), V.get(), Ra_ret, size());
	return cVector(Ra_ret, size());
}

cVector cVector::operator*(wsRealCst R_val)
{
	wsVec Ra_ret = sseVector(size());
	sseVectorMul(size(), get(), R_val, Ra_ret);
	return cVector(Ra_ret, size());
}

cVector cVector::operator*(const cVector &V)
{
	ASSERT_SIZE(size(), V.size(), "Vector multiplication");
	wsVec Ra_ret = sseVector(size());
	sseVectorMul(size(), get(), V.get(), Ra_ret);
	return cVector(Ra_ret, size());
}

cVector cVector::operator*(const cStdMatrix &C_mat)
{
	wsReal *V = get();
	wsReal **R = sseMpM(&V, 1, size(), C_mat.get(), C_mat.row(),
		C_mat.col());
	wsReal *RR = R[0];
	DEALLOC(R);
	return cVector(RR, C_mat.col());
}

cVector cVector::operator*(const cSymMatrix &S)
{
	if (size() != S.row())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector-matrix multiplication",
			size(), S.row());
	wsReal *V = get();
	wsReal **R = sseMpS(&V, 1, size(), S.get(), S.row());
	wsReal *Rp = R[0];
	DEALLOC(R);
	return cVector(Rp, size());
}

cVector cVector::Mt(const cStdMatrix &M)
{
	if (size() != M.col())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector-matrix multiplication",
			size(), M.col());
	wsReal *V = get();
	wsReal **R = sseMpMt(&V, 1, size(), M.get(), M.row(), M.col());
	wsReal *Rp = R[0];
	DEALLOC(R);
	return cVector(Rp, M.row());
}

cVector& cVector::operator*=(cDiagMatrix &D)
{
	if (size() != D.row())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector-matrix multiplication",
			size(), D.row());
	sseVpV(get(), D.get()[0], get(), size());

	return *this;
}

cVector& cVector::operator*=(cIdtMatrix &I)
{
	if (size() != I.row())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector-matrix multiplication",
			size(), I.row());
	return *this;
}

cVector& cVector::operator*=(cSpsMatrix &P)
{
	if (size() != P.row() || !P.sqr())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector-matrix multiplication",
			size(), P.row());
	sseVpP(get(), size(), P.get()[0], P.col(), P.getRowEnds(), P.getColIdxs());
	return *this;
}

cVector cVector::operator*(const cDiagMatrix &D)
{
	wsReal *pV = sseVector(size());
	if (size() != D.row())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector-matrix multiplication",
		size(), D.row());
	sseVpV(get(), D.get()[0], pV, size());

	return cVector(pV, size());
}

cVector cVector::operator*(const cIdtMatrix &I)
{
	wsReal *pV = sseVector(size());
	memcpy(pV, get(), sizeof(wsReal)*size());

	return cVector(pV, size());
}

cVector cVector::operator*(const cSpsMatrix &P)
{
	wsReal *pV = sseVector(size());
	if (size() != P.row())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector-matrix multiplication",
			size(), P.row());
	sseVpP(get(), size(), P.get()[0], P.col(), P.getRowEnds(), P.getColIdxs());

	return cVector(pV, size());
}

cVector cVector::sq()
{
	wsReal *Ra_ret = sseVector(size());
	sseVsquare(get(), size(), Ra_ret);
	return cVector(Ra_ret, size());
}

cVector cVector::sq(cVector& V)
{
	wsReal *Ra_ret = sseVector(size());
	sseVsquare(get(), size(), Ra_ret, V.get());
	return cVector(Ra_ret, size());
}

cVector cVector::sqrt()
{
	wsReal *Ra_ret = sseVector(size());
	sseVsqrt(get(), size(), Ra_ret);
	return cVector(Ra_ret, size());
}

void cVector::selfSqrt(char B_abs)
{
	sseVabsqrt(get(), size());
}

cVector cVector::cb()
{
	wsReal *Ra_ret = sseVector(size());
	sseVcube(get(), size(), Ra_ret);
	return cVector(Ra_ret, size());
}

void cVector::mul(wsRealCst R_multipier, cVector &V_dest)
{
	if (size() != V_dest.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "Vector-scalar multiplication",
			size(), V_dest.size());
	sseVpC(get(), R_multipier, V_dest.get(), size());
}

void cVector::mul(wsRealCst R_multipier, wsRealCst R_adder, cVector &V_dest)
{
	if (size() != V_dest.size())
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "Vector-scalar multiplication",
		size(), V_dest.size());
	sseVpC(get(), R_multipier, V_dest.get(), size(), R_adder);
}

void cVector::mul(wsRealCst R_multipier, wsRealCst R_adder)
{
	sseVpC(get(), R_multipier, get(), size(), R_adder);
}

cVector cVector::operator/(wsRealCst R_val)
{
	wsReal *Ra_ret	= sseVector(size());

	sseVectorDiv(size(), get(), R_val, Ra_ret);
	return cVector(Ra_ret, size());
}

cVector cVector::operator/(const cVector &V)
{
	/* Sanity check : Length */
	if (V.size() != size())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector division", size(), V.size());
	wsReal *Ra_ret	= sseVector(size());

	sseVectorDiv(size(), get(), V.get(), Ra_ret);
	return cVector(Ra_ret, size());
}

void cVector::file(wsStrCst S_ext)
{
	exportMatrix(S_ext, &Ra_buf, 1, N_sz);
}

void cVector::fmtfile(wsStrCst S_fmtExt, ...)
{
	va_list	H_varList;
	char	S_buf[256];
	va_start(H_varList, S_fmtExt);
	vsprintf(S_buf, S_fmtExt, H_varList);
	file(S_buf);
}

cVector& cVector::sub(cVector &V_s, cVector &V_t)
{
	if (size()!=V_s.size())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector multiplication", size(), V_s.size());
	if (size()!=V_t.size())
		halt_fmt(WISARD_SYST_INVL_DIM, "Vector multiplication", size(), V_t.size());
	sseVsV(get(), V_s.get(), V_t.get(), size());
	return V_t;
}

wsRealCst cVector::ss()
{
	return sseVV(get(), size());
}

cVector cVector::clone(wsUint N_len/*=0*/)
{
	if (N_len > size()) halt_fmt(WISARD_SYST_INVL_DIM, "Vector cloning with length", size(), N_len);
	else if (N_len == 0) N_len = size();
	wsReal *R = NULL;
	sseMalloc(R, wsReal, N_len);
	memcpy(R, get(), sizeof(wsReal)*N_len);
	return cVector(R, N_len);
}

cVector cVector::subset(wsUint N_start, wsUint N_sz)
{
	wsReal	*V	= get();
	wsReal	*R	= sseVector(N_sz);

	/* Overflow check */
	if ((N_start+N_sz) > size())
		halt("SYSERR : Subset range (%d~%d) overflow (>%d)", N_start,
		N_start+N_sz-1, size());
	/* Copy */
	memcpy(R, V+N_start, sizeof(wsReal)*N_sz);

	/* Return */
	return cVector(R, N_sz);
}

cSymMatrix cVector::tV(wsReal R_prod/*=WISARD_NAN*/)
{
	wsUint	N	= size();
	wsSym	R	= sseSymMat(N);
	wsReal	*P	= get();

	if (NA(R_prod)) for (wsUint i=0 ; i<N ; i++)
		sseVpC(P, P[i], R[i], i+1);
	else for (wsUint i=0 ; i<N ; i++)
		sseVpC(P, P[i]*R_prod, R[i], i+1);
	return cSymMatrix(R, N);
}

cStdMatrix cVector::tV(cVector &V, wsReal R_prod/*=WISARD_NAN*/)
{
	wsUint	N	= size();
	wsReal	*P	= get();

	wsMat R = sseVtV(P, N, V.get(), V.size(), R_prod);
	return cStdMatrix(N, V.size(), R);
}

cVector& cVector::neg(wsReal R_add/*=WISARD_NAN*/)
{
	sseVneg(get(), size(), NULL, R_add);
	return *this;
}

cVector cVector::reorder(wsUint *Na_seq)
{
	wsUint	N		= size();
	wsReal*	Ra_ret	= sseVector(N);
	for (wsUint i=0 ; i<N ; i++)
		Ra_ret[i] = Ra_buf[Na_seq[i]];

	return cVector(Ra_buf, N);
}

cVector cVector::krnk(cVector &V)
{
	wsUint	N1		= size();
	wsUint	N2		= V.size();
	wsUint	Na		= N1 * N2;
	wsVec	V1		= get();
	wsVec	V2		= V.get();
	wsVec	Ra_ret	= sseVector(Na);
	for (wsUint i=0,k=0 ; i<N1 ; i++) {
		wsReal E	= V1[i];
		for (wsUint j=0 ; j<N2 ; j++,k++)
			Ra_ret[k] = E * V2[j];
	}
	return cVector(Ra_ret, Na);
}

cStdMatrix cVector::krnkV(cStdMatrix &M)
{
	wsUint	R		= size() * M.row();
	wsUint	C		= M.col();
	wsVec	V		= get();
	wsMat	Mp		= M.get();
	wsMat	Ra_ret	= sseMatrix(R, C);
	for (wsUint i=0,k=0 ; i<size() ; i++) {
		wsReal E	= V[i];
		for (wsUint j=0 ; j<M.row() ; j++,k++)
			sseVpC(Mp[j], E, Ra_ret[k], C);
	}
	return cStdMatrix(R, C, Ra_ret);
}

void cVector::minSelf(wsRealCst R_val)
{
	wsUint	N_sz	= size();
	wsVec	Ra_v1	= get();
//	wsUint	N_med	= 0;

	//#ifdef USE_SSE
	//#else
	LOOP(i, N_sz)
		Ra_v1[i] = min(R_val, Ra_v1[i]);
	//#endif
}

void cVector::maxSelf(wsRealCst R_val)
{
	wsUint	N_sz	= size();
	wsVec	Ra_v1	= get();
//	wsUint	N_med	= 0;

	//#ifdef USE_SSE
	//#else
	LOOP (i, N_sz)
		Ra_v1[i] = max(R_val, Ra_v1[i]);
	//#endif
}

cVector cVector::log1exp()
{
	wsUint	N_sz	= size();
	wsVec	Ra_ret	= sseVector(N_sz);
	wsVec	Ra_v1	= get();
//	wsUint	N_med	= 0;

	//#ifdef USE_SSE
	//#else
	LOOP(i, N_sz)
		Ra_ret[i] = log(1+::exp(Ra_v1[i]));
	//#endif
	return cVector(Ra_ret, N_sz);
}

cVector cVector::logit()
{
	wsUint	N_sz	= size();
	wsVec	Ra_ret	= sseVector(N_sz);
	wsVec	Ra_v1	= get();
//	wsUint	N_med	= 0;

	//#ifdef USE_SSE
	//#else
	LOOP(i, N_sz)
		Ra_ret[i] = ::exp(Ra_v1[i])/(1+::exp(Ra_v1[i]));
	//#endif
	return cVector(Ra_ret, N_sz);
}

cVector cVector::divide(wsUint* Na_idx, wsUint N_szIdx, cVector* Vp_rest/*=NULL*/)
{
	wsVec		Ra_ret	= sseVector(N_szIdx);
	wsUintCst	N_sz	= size();
	char*		Ba_map	= NULL;
	wsCalloc(Ba_map, char, N_sz);
	LOOP (i, N_szIdx)
		Ba_map[Na_idx[i]] = 1;
	wsVec		Ra_v1	= get();
	wsUint		j		= 0;
	wsUint		k		= 0;

	if (Vp_rest) {
		Vp_rest->init(N_sz - N_szIdx);
		wsVec	Ra_rest = Vp_rest->get();
		LOOP (i, N_sz)
			if (Ba_map[i] == 1)
				Ra_ret[j++] = Ra_v1[i];
			else
				Ra_rest[k++] = Ra_v1[i];
	} else LOOP (i, N_sz)
		if (Ba_map[i] == 1) Ra_ret[j++] = Ra_v1[i];
	
	DEALLOC(Ba_map);
	return cVector(Ra_ret, N_szIdx);
}

void cVector::replace(wsUint* Na_idx, wsUint N_szIdx, wsVec Ra_rep, char B_inv/* =0 */)
{
	wsUintCst	N_sz	= size();
	wsVec		Ra_v	= get();
	char*		Ba_map	= NULL;
	wsCalloc(Ba_map, char, N_sz);
	LOOP(i, N_szIdx)
		Ba_map[Na_idx[i]] = 1;
	wsUint		j = 0;

	if (B_inv) {
		LOOP (i, N_sz) if (Ba_map[i] == 0)
			Ra_v[i] = Ra_rep[j++];
	} else LOOP (i, N_sz) if (Ba_map[i] == 1)
			Ra_v[i] = Ra_rep[j++];

	DEALLOC(Ba_map);
}

cVector cVector::exp()
{
	wsUint	N_sz	= size();
	wsVec	Ra_ret	= sseVector(N_sz);
	wsVec	Ra_v1	= get();
//	wsUint	N_med	= 0;

//#ifdef USE_SSE
//#else
	LOOP (i, N_sz)
		Ra_ret[i] = ::exp(Ra_v1[i]);
//#endif
	return cVector(Ra_ret, N_sz);
}

// V * v1 * v2
wsRealCst cVector::VV(wsReal *Ra_1, wsReal *Ra_2, wsUint N_sz)
{
	return sseVVV(get(), Ra_1, Ra_2, N_sz);
}

// V^2 * v
wsRealCst cVector::V2(wsReal *Ra_v, wsUint N_sz)
{
	wsReal	R_ret	= REAL_CONST(0.0);
	wsReal*	Ra_v1	= get();
	wsUint	N_med	= 0;

#ifdef USE_SSE
	N_med = getMed(N_sz);
	sse_t sse_ret = sseSet(0.0);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v + i);
		sse_ret = sseAdd(sse_ret, sseMul(sseMul(*sse_v1, *sse_v1), *sse_v2));
	}
	sseSum(sse_ret, R_ret);
#endif
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_ret += SQR(Ra_v1[i])*Ra_v[i];

	return R_ret;
}

char cVector::isAnyNeg()
{
	wsUint	N_sz	= size();
	wsReal*	Ra_v	= get();
	for (wsUint i=0 ; i<N_sz ;i++)
		if (Ra_v[i] < W0) return 1;
	return 0;
}

void cVector::set(wsRealCst R_v)
{
	sseVinit(get(), size(), R_v);
}

// weighted mean vec * Ra_w * N_wsz / sum(Ra_w)
wsRealCst cVector::wMean(wsReal *Ra_w, wsUint N_wsz)
{
	/* Dimension check */
	if (size() != N_wsz) halt_fmt(WISARD_SYST_INVL_DIM, size(), N_wsz);

	/* Now compute weighted mean */
	wsRealCst R_wsum = sseVsum(Ra_w, N_wsz);
	return sseVV(get(), N_wsz, Ra_w, (wsReal)N_wsz/R_wsum);
}

void cVector::pinv()
{
	wsUintCst	N = size();
	wsReal		w = W0;
	wsVec		R = get();
	dsvd(R, N, &w);
//	wsReal tol = MAX(pow(numeric_limits<wsReal>::epsilon(), 2/3), W0);

// 	wsUint n = 0;
	LOOP (i, N)
		R[i] /= w;
}

cVector& cVector::subFrom(cVector &V)
{
	/* Dimension check */
	if (size() != V.size()) halt_fmt(WISARD_SYST_INVL_DIM, size(), V.size());
	sseVsV(V.get(), get(), get(), size());
	return *this;
}

cVector cVector::permute()
{
	wsUint	N_sz		= size();
	wsVec	Ra_yCopy	= sseVectorP(N_sz, get());
	wsVec	Ra_ret		= sseVector(N_sz);
	wsUint	N_remained	= N_sz;

	/* Permute it */

	LOOP(j, N_sz) {
		wsUint N_idxPick = wsRand()%N_remained;
		Ra_ret[j] = Ra_yCopy[N_idxPick];
		memmove(Ra_yCopy+N_idxPick, Ra_yCopy+N_idxPick+1,
			sizeof(wsReal)*(N_remained-N_idxPick-1));
		N_remained--;
	}
	sseFree(Ra_yCopy);

	return cVector(Ra_ret, N_sz);
}

cVector cVector::permute(wsUint* Na_seq)
{
	wsUint	N_sz		= size();
	wsVecCst	Ra_buf		= get();
	wsVec	Ra_ret		= sseVector(N_sz);

	/* Permute it with fixed sequence*/
	LOOP(j, N_sz)
		Ra_ret[j] = Ra_buf[Na_seq[j]];

	return cVector(Ra_ret, N_sz);
}

/*
 *
 * cMask definitions
 *
 */

wsRealCst cMask::vv(cVector &V, cVector &W)
{
	if (size()!=V.size() || size()!=W.size()) {
		wsUint N_ratio = V.size()%size();
		if (V.size()!=W.size() || N_ratio)
			halt_fmt(WISARD_SYST_INVL_SAMEORMULTSZ, "Source mask size",
				"target vector", size(), V.size());
		/* If vector V and W is multiple of mask */
		N_ratio = (wsUint)(V.size()/size());
		cMask C_mulMask(*this, N_ratio);
		return sseVVsel(C_mulMask.get(), V.get(), size(), W.get());
	}
	return sseVVsel(get(), V.get(), size(), W.get());
}

wsRealCst cMask::vv(cVector &V)
{
	if (size()!=V.size()) {
		wsUint N_ratio = V.size()%size();
		if (N_ratio)
			halt_fmt(WISARD_SYST_INVL_SAMEORMULTSZ, "Source mask size",
				"target vector", size(), V.size());
		/* If vector V and W is multiple of mask */
		N_ratio = (wsUint)(V.size()/size());
		cMask C_mulMask(*this, N_ratio);
		return sseVVsel(C_mulMask.get(), V.get(), size());
	}
	return sseVVsel(get(), V.get(), size());
}

wsRealCst cMask::qf(cVector &V, cSymMatrix &M)
{
	if (!M.sqr() || size()!=V.size() || size()!=M.row()) halt("ERR");
	return sseVpSpV(V.get(), M.get(), size(), V.get(), get());
}

wsRealCst cMask::qf(cVector &V, cSymMatrix &M, cVector &W)
{
	if (!M.sqr() || size()!=V.size() || size()!=W.size() || size()!=M.row()) halt("ERR");
	return sseVpSpV(V.get(), M.get(), size(), W.get(), get());
}

wsRealCst cMask::qf(cVector &V, cIdtMatrix &M)
{
	if (!M.sqr() || size()!=V.size() || size()!=M.row()) halt("ERR");
	return sseVVsel(get(), V.get(), V.size());
}

wsRealCst cMask::qf(cVector &V, cIdtMatrix &M, cVector &W)
{
	if (!M.sqr() || size()!=V.size() || size()!=W.size() || size()!=M.row()) halt("ERR");
	return sseVVsel(get(), V.get(), V.size(), W.get());
}

void sseVinit(wsVec Ra_v, wsUintCst N_sz, wsRealCst R_val)
{
	wsUint N_med = 0;
#ifdef USE_SSE
	N_med = getMed(N_sz);
	sse_t sse_v = sseSet(R_val);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp)
		*((sse_t *)(Ra_v + i)) = sse_v;
#endif
	for (wsUint i=N_med ; i<N_sz ; i++)
		Ra_v[i] = R_val;
}

void sseVset0(wsVec Ra_v, wsUintCst N_sz)
{
	wsUint N_med = 0;
#ifdef USE_SSE
	N_med = getMed(N_sz);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_p = (sse_t *)(Ra_v + i);
		*sse_p = sseXor(*sse_p, *sse_p);
	}
#endif
	for (wsUint i=N_med ; i<N_sz ; i++)
		Ra_v[i] = REAL_CONST(0.0);
}

wsVec sseVrep(wsVecCst Ra_v, wsUintCst N_sz, wsUintCst N_times)
{
	wsVec	Ra_ret	= sseVector(N_sz*N_times);
	for (wsUint i=0,j=0 ; i<N_times ; i++,j+=N_sz)
		memcpy(Ra_ret+j, Ra_v, sizeof(wsReal)*N_sz);
	return Ra_ret;
}

void sseMinit(wsMat Ra_m, wsUintCst N_row, wsUintCst N_col, wsRealCst R_val)
{
	wsUint N_med = 0;
#ifdef USE_SSE
	N_med = getMed(N_col);
#endif
	for (wsUint i=0 ; i<N_row ; i++) {
		wsVec	Ra_v	= Ra_m[i];
#ifdef USE_SSE
		sse_t	sse_v	= sseSet(R_val);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp)
			*((sse_t *)(Ra_v + j)) = sse_v;
#endif
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_v[j] = R_val;
	}
}

wsReal** sseMrand(wsUint N_r, wsUint N_c)
{
	wsReal** Ra_ret = sseMatrix(N_r, N_c);

	for (wsUint i=0 ; i<N_r ; i++)
		for (wsUint j=0 ; j<N_c ; j++)
			Ra_ret[i][j] = (wsReal)(wsRand()%100000) / REAL_CONST(1000.0);

	return Ra_ret;
}

wsMat sseMrand01(wsUint N_r, wsUint N_c, wsReal R_prob)
{
	wsReal** Ra_ret = sseMatrix(N_r, N_c);

	for (wsUint i=0 ; i<N_r ; i++)
		for (wsUint j=0 ; j<N_c ; j++) {
			wsReal R = (wsReal)(wsRand()%1000) / REAL_CONST(999.0);
			Ra_ret[i][j] = R > R_prob ? W0 : W1;
		}

		return Ra_ret;
}

wsReal sseVvar(wsVecCst Ra_v, wsUintCst N_sz, wsRealPtr Rp_mean/*=NULL*/)
{
	wsReal R_mean = sseVsum(Ra_v, N_sz) / (wsReal)N_sz;
	if (Rp_mean) *Rp_mean = R_mean;
	return sseVsqsum(Ra_v, N_sz, -R_mean) / (wsReal)(N_sz-1);
}

wsReal sseVsum(wsVecCst Ra_v, wsUintCst N_sz, wsRealPtr Rp_sqSum)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;

	if (Rp_sqSum) {
		wsReal	R_sqRet		= W0;
#ifdef USE_SSE
		sse_t	sse_ret		= sseSet(0.0);
		sse_t	sse_sqRet	= sseSet(0.0);
		N_med				= getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v+i);
			sse_ret = sseAdd(sse_ret, *sse_v);
			sse_sqRet = sseAdd(sse_sqRet, sseMul(*sse_v, *sse_v));
		}
		sseSum(sse_ret, R_ret);
		sseSum(sse_sqRet, R_sqRet);
#endif
		for (wsUint i=N_med ; i<N_sz ; i++) {
			R_ret += Ra_v[i];
			R_sqRet += SQR(Ra_v[i]);
		}
		*Rp_sqSum = R_sqRet;
	} else {
#ifdef USE_SSE
		sse_t sse_ret = sseSet(0.0);
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp)
			sse_ret = sseAdd(sse_ret, *(sse_t *)(Ra_v+i));
		sseSum(sse_ret, R_ret);
#endif
		for (wsUint i=N_med ; i<N_sz ; i++)
			R_ret += Ra_v[i];
	}

	return R_ret;
}

wsReal sseVsqsum(wsVecCst Ra_v, wsUintCst N_sz, wsRealCst R_adder/*=WISARD_NAN*/)
{
	wsUint	N_med	= 0;
	wsReal	R_sqRet	= W0;

	if (NA(R_adder)) {
#ifdef USE_SSE
		sse_t	sse_sqRet	= sseSet(0.0);
		N_med				= getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v+i);
			sse_sqRet = sseAdd(sse_sqRet, sseMul(*sse_v, *sse_v));
		}
		sseSum(sse_sqRet, R_sqRet);
#endif
		for (wsUint i=N_med ; i<N_sz ; i++)
			R_sqRet += SQR(Ra_v[i]);
	} else {
#ifdef USE_SSE
		sse_t	sse_a		= sseSet(R_adder);
		sse_t	sse_sqRet	= sseSet(0.0);
		N_med				= getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v+i);
			sse_t sse_r = sseAdd(*sse_v, sse_a);
			sse_sqRet = sseAdd(sse_sqRet, sseMul(sse_r, sse_r));
		}
		sseSum(sse_sqRet, R_sqRet);
#endif
		for (wsUint i=N_med ; i<N_sz ; i++) {
			wsReal R_v = Ra_v[i] + R_adder;
			R_sqRet += SQR(R_v);
		}
	}

	return R_sqRet;
}

wsReal sseVasum(wsVecCst Ra_v, wsUintCst N_sz, wsRealPtr Rp_sqSum)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;

	if (Rp_sqSum) {
		wsReal	R_sqRet		= W0;
#ifdef USE_SSE
		sse_t	sse_ret		= sseSet(0.0);
		sse_t	sse_sqRet	= sseSet(0.0);
		N_med				= getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v+i);
			sse_ret = sseAdd(sse_ret, sseAbs(*sse_v));
			sse_sqRet = sseAdd(sse_sqRet, sseMul(*sse_v, *sse_v));
		}
		sseSum(sse_ret, R_ret);
		sseSum(sse_sqRet, R_sqRet);
#endif
		for (wsUint i=N_med ; i<N_sz ; i++) {
			R_ret += fabs(Ra_v[i]);
			R_sqRet += SQR(Ra_v[i]);
		}
		*Rp_sqSum = R_sqRet;
	} else {
#ifdef USE_SSE
		sse_t sse_ret = sseSet(0.0);
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp)
			sse_ret = sseAdd(sse_ret, sseAbs(*(sse_t *)(Ra_v+i)));
		sseSum(sse_ret, R_ret);
#endif
		for (wsUint i=N_med ; i<N_sz ; i++)
			R_ret += fabs(Ra_v[i]);
	}

	return R_ret;
}

wsReal sseVlogasum(wsVecCst Ra_v, wsUintCst N_sz)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;

	// #ifdef USE_SSE
	// 	sse_t sse_ret = sseSet(0.0);
	// 	N_med = getMed(N_sz);
	// 	for (wsUint i=0 ; i<N_med ; i+=sseJmp)
	// 		sse_ret = sseAdd(sse_ret, *(sse_t *)(Ra_v+i));
	// 	sseSum(sse_ret, R_ret);
	// #endif
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_ret += log(fabs(Ra_v[i]));

	return R_ret;
}

wsReal sseVlogsum(wsVecCst Ra_v, wsUintCst N_sz)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;

	for (wsUint i=N_med ; i<N_sz ; i++)
		R_ret += log(Ra_v[i]);

	return R_ret;
}

wsRealCst sseVprod(wsVecCst Ra_v, wsUintCst N_sz)
{
	wsReal	R_ret	= W1;
	wsUint	N_med	= 0;

#if 0
#ifdef USE_SSE
	sse_t sse_ret = sseSet(0.0);
	N_med = getMed(N_sz);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp)
		sse_ret = sseMul(sse_ret, *(sse_t *)(Ra_v+i));
	sseProd(sse_ret, R_ret);
#endif
#endif
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_ret *= Ra_v[i];

	return R_ret;
}

// Do sum(v1+v2)*m or sum(v1+v2)
wsReal sseVaVsum(wsVecCst Ra_v1, wsVecCst Ra_v2, wsUintCst N_sz,
	wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint	N_med	= 0;
	wsReal	R_sum	= W0;

#ifdef USE_SSE
	sse_t sse_sum = sseSet(0.0);
	N_med = getMed(N_sz);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);

		sse_sum = sseAdd(sse_sum, sseAdd(*sse_v1, *sse_v2));
	}
	sseSum(sse_sum, R_sum);
#endif
	/* Do rest part */
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_sum += Ra_v1[i] + Ra_v2[i];

	return NA(R_mul) ? R_sum : R_sum*R_mul;
}

// Do sum(v1-v2)*m or sum(v1-v2)
wsReal sseVsVsum(wsVecCst Ra_v1, wsVecCst Ra_v2, wsUintCst N_sz, char B_sq/*=0*/,
				 wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint	N_med	= 0;
	wsReal	R_sum	= W0;

	if (B_sq) {
#ifdef USE_SSE
		sse_t sse_sum = sseSet(0.0);
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t sse_tmp = sseSub(*sse_v1, *sse_v2);

			sse_sum = sseAdd(sse_sum, sseMul(sse_tmp, sse_tmp));
		}
		sseSum(sse_sum, R_sum);
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++) {
			wsReal R_tmp = Ra_v1[i] - Ra_v2[i];
			R_sum += SQR(R_tmp);
		}
	} else {
#ifdef USE_SSE
		sse_t sse_sum = sseSet(0.0);
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);

			sse_sum = sseAdd(sse_sum, sseSub(*sse_v1, *sse_v2));
		}
		sseSum(sse_sum, R_sum);
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			R_sum += Ra_v1[i] - Ra_v2[i];
	}

	return NA(R_mul) ? R_sum : R_sum*R_mul;
}

wsReal sseVsCsum(wsVecCst Ra_v1, wsRealCst R_v, wsUintCst N_sz, char B_sq/*=0*/)
{
	wsUint	N_med	= 0;
	wsReal	R_sum	= W0;

	if (B_sq) {
#ifdef USE_SSE
		sse_t sse_sum = sseSet(0.0);
		sse_t sse_V = sseSet(R_v);
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t sse_tmp = sseSub(*sse_v1, sse_V);

			sse_sum = sseAdd(sse_sum, sseMul(sse_tmp, sse_tmp));
		}
		sseSum(sse_sum, R_sum);
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++) {
			wsReal R_tmp = Ra_v1[i] - R_v;
			R_sum += SQR(R_tmp);
		}
	} else {
#ifdef USE_SSE
		sse_t sse_sum = sseSet(0.0);
		sse_t sse_V = sseSet(R_v);
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);

			sse_sum = sseAdd(sse_sum, sseSub(*sse_v1, sse_V));
		}
		sseSum(sse_sum, R_sum);
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			R_sum += Ra_v1[i] - R_v;
	}

	return R_sum;
}

wsReal sseVV(wsVecCst Ra_v1, wsUintCst N_sz, wsVec Ra_v2/*=NULL*/,
			 wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint	N_med	= 0;
	wsReal	R_ret	= W0;

	if (NA(R_mul)) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_ret = sseSet(0.0);
#endif
		if (Ra_v2 == NULL) {
#ifdef USE_SSE
			for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
				sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
				sse_ret = sseAdd(sse_ret, sseMul(*sse_v1, *sse_v1));
			}
			sseSum(sse_ret, R_ret);
#endif
			for (wsUint i=N_med ; i<N_sz ; i++)
				R_ret += SQR(Ra_v1[i]);
		} else {
#ifdef USE_SSE
			for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
				sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
				sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
				sse_ret = sseAdd(sse_ret, sseMul(*sse_v1, *sse_v2));
			}
			sseSum(sse_ret, R_ret);
#endif
			for (wsUint i=N_med ; i<N_sz ; i++)
				R_ret += Ra_v1[i]*Ra_v2[i];
		}
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_ret = sseSet(0.0);
		sse_t sse_mul = sseSet(R_mul);
#endif

		if (Ra_v2 == NULL) {
#ifdef USE_SSE
			for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
				sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
				sse_ret = sseAdd(sse_ret, sseMul(sseMul(*sse_v1, *sse_v1), sse_mul));
			}
			sseSum(sse_ret, R_ret);
#endif
			for (wsUint i=N_med ; i<N_sz ; i++)
				R_ret += SQR(Ra_v1[i]) * R_mul;
		} else {
#ifdef USE_SSE
			for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
				sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
				sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
				sse_ret = sseAdd(sse_ret, sseMul(sseMul(*sse_v1, *sse_v2), sse_mul));
			}
			sseSum(sse_ret, R_ret);
#endif
			for (wsUint i=N_med ; i<N_sz ; i++)
				R_ret += Ra_v1[i]*Ra_v2[i] * R_mul;
		}
	}

	return R_ret;
}

wsReal sseVVV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVecCst Ra_v3, wsUintCst N_sz)
{
	wsReal R_ret = W0;

#ifdef USE_SSE
	wsUint N_med = getMed(N_sz);
	sse_t sse_ret = sseSet(0.0f);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
		sse_t *sse_v3 = (sse_t *)(Ra_v3 + i);
		sse_ret = sseAdd(sse_ret, sseMul(sseMul(*sse_v1, *sse_v2), *sse_v3));
	}
	sseSum(sse_ret, R_ret);
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_ret += Ra_v1[i]*Ra_v2[i]*Ra_v3[i];
#else
	for (wsUint i=0 ; i<N_sz ; i++)
		R_ret += Ra_v1[i]*Ra_v2[i]*Ra_v3[i];
#endif

	return R_ret;
}

/* Do v <- sqrt(v) */
void sseVsqrt(wsVec Ra_V, wsUintCst N_sz, wsVec Ra_tgt/*=NULL*/)
{
	wsUint i;
	wsUint N_med = 0;

	N_med = getMed(N_sz);
	if (Ra_tgt) {
#ifdef USE_SSE
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_eVal	= (sse_t *)(Ra_V + i);
			sse_t *sse_tgt	= (sse_t *)(Ra_tgt + i);
			*sse_tgt = sseSqrt(*sse_eVal);
		}
#endif
		for (i=N_med ; i<N_sz ; i++)
			Ra_tgt[i] = sqrt(Ra_V[i]);
	} else {
#ifdef USE_SSE
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_eVal = (sse_t *)(Ra_V + i);
			*sse_eVal = sseSqrt(*sse_eVal);
		}
#endif
		for (i=N_med ; i<N_sz ; i++)
			Ra_V[i] = sqrt(Ra_V[i]);
	}
}

/* Do v <- sqrt(|v|) */
void sseVabsqrt(wsVec Ra_V, wsUintCst N_sz, wsVec Ra_tgt/*=NULL*/)
{
	wsUint i;
	wsUint N_med = 0;

	N_med = getMed(N_sz);
	if (Ra_tgt) {
#ifdef USE_SSE
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_eVal	= (sse_t *)(Ra_V + i);
			sse_t *sse_tgt	= (sse_t *)(Ra_tgt + i);
			*sse_tgt = sseSqrt(sseAbs(*sse_eVal));
		}
#endif
		for (i=N_med ; i<N_sz ; i++)
			Ra_tgt[i] = sqrt(fabs(Ra_V[i]));
	} else {
#ifdef USE_SSE
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_eVal = (sse_t *)(Ra_V + i);
			*sse_eVal = sseSqrt(sseAbs(*sse_eVal));
		}
#endif
		for (i=N_med ; i<N_sz ; i++)
			Ra_V[i] = sqrt(fabs(Ra_V[i]));
	}
}

/* Do v2 <- v*v */
void sseVsquare(wsVecCst Ra_V, wsUintCst N_sz, wsVec Ra_vt, wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint i;
	wsUint N_med = 0;

	if (NA(R_mul)) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_src	= (sse_t *)(Ra_V + i);
			sse_t *sse_tgt	= (sse_t *)(Ra_vt + i);
			*sse_tgt = sseMul(*sse_src, *sse_src);
		}
#endif
		for (i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = SQR(Ra_V[i]);
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_mul = sseSet(R_mul);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_src	= (sse_t *)(Ra_V + i);
			sse_t *sse_tgt	= (sse_t *)(Ra_vt + i);
			*sse_tgt = sseMul(sseMul(*sse_src, *sse_src), sse_mul);
		}
#endif
		for (i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = SQR(Ra_V[i]) * R_mul;
	}
}

/* Do v2 <- v*v */
void sseVsquare(wsVecCst Ra_V, wsUintCst N_sz, wsVec Ra_vt, wsVec Ra_v2)
{
	wsUint i;
	wsUint N_med = 0;

#ifdef USE_SSE
	N_med = getMed(N_sz);
	for (i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_src	= (sse_t *)(Ra_V + i);
		sse_t *sse_tgt	= (sse_t *)(Ra_vt + i);
		sse_t *sse_s2	= (sse_t *)(Ra_v2 + i);
		*sse_tgt = sseMul(sseMul(*sse_src, *sse_src), *sse_s2);
	}
#endif
	for (i=N_med ; i<N_sz ; i++)
		Ra_vt[i] = SQR(Ra_V[i]) * Ra_v2[i];
}

/* Do v3 <- v*v*v */
void sseVcube(wsVecCst Ra_V, wsUintCst N_sz, wsVec Ra_v3)
{
	wsUint i;

#ifdef USE_SSE
	wsUint N_med = getMed(N_sz);
	for (i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_src	= (sse_t *)(Ra_V + i);
		sse_t *sse_tgt	= (sse_t *)(Ra_v3 + i);
		*sse_tgt = sseMul(*sse_src, sseMul(*sse_src, *sse_src));
	}
	for (i=N_med ; i<N_sz ; i++)
		Ra_v3[i] = CUBE(Ra_V[i]);
#else
	for (i=0 ; i<N_sz ; i++)
		Ra_v3[i] = CUBE(Ra_V[i]);
#endif
}

// Do vt = v1+s
void sseVaC(wsVecCst Ra_v, wsRealCst R_const, wsVec Ra_vt, wsUintCst N_sz)
{
	wsUint	N_med = 0, i;
	sse_t	sse_s = sseSet(R_const);

#ifdef USE_SSE
	N_med = getMed(N_sz);
	for (i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v + i);
		sse_t *sse_vt = (sse_t *)(Ra_vt + i);
		*sse_vt = sseAdd(*sse_v1, sse_s);
	}
#endif
	/* Do rest part */
	for (i=N_med ; i<N_sz ; i++)
		Ra_vt[i] = Ra_v[i] + R_const;
}

// Do vt = |v1|+v2 or vt = (|v1|+v2)*m
void sseAVaV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
			 wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint	N_med = 0;
#ifdef USE_SSE
	N_med = getMed(N_sz);
#endif

	if (NA(R_mul)) {
#ifdef USE_SSE
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseAdd(sseAbs(*sse_v1), *sse_v2);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = fabs(Ra_v1[i]) + Ra_v2[i];
	} else {
		sse_t sse_m = sseSet(R_mul);

#ifdef USE_SSE
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseMul(sseAdd(sseAbs(*sse_v1), *sse_v2), sse_m);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = (fabs(Ra_v1[i]) + Ra_v2[i]) * R_mul;
	}
}

// Do vt = (v1-v2)*m or vt = v1-v2
void sseVsV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
			wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint N_med = 0;

	if (NA(R_mul)) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseSub(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v1[i] - Ra_v2[i];
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_M = sseSet(R_mul);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseMul(sseSub(*sse_v1, *sse_v2), sse_M);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = (Ra_v1[i] - Ra_v2[i]) * R_mul;
	}
}

// Do vt = (v1-v2)^2
void sseVsVsq(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz)
{
	wsUint N_med = 0;

#ifdef USE_SSE
	N_med = getMed(N_sz);
	for (wsUint i=0 ; i < N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
		sse_t *sse_vt = (sse_t *)(Ra_vt + i);
		sse_t sse_tmp = sseSub(*sse_v1, *sse_v2);
		*sse_vt = sseMul(sse_tmp, sse_tmp);
	}
#endif
	/* Do rest part */
	for (wsUint i=N_med ; i < N_sz ; i++) {
		wsRealCst R_v = Ra_v1[i] - Ra_v2[i];
		Ra_vt[i] = R_v*R_v;
	}
}

// Do vt = v1-v2*m
void sseVsVC(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
			 wsRealCst R_mul)
{
	wsUint N_med = 0;

#ifdef USE_SSE
	N_med = getMed(N_sz);
	sse_t sse_M = sseSet(R_mul);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
		sse_t *sse_vt = (sse_t *)(Ra_vt + i);
		*sse_vt = sseSub(*sse_v1, sseMul(*sse_v2, sse_M));
	}
#endif
	/* Do rest part */
	for (wsUint i=N_med ; i<N_sz ; i++)
		Ra_vt[i] = Ra_v1[i] - Ra_v2[i] * R_mul;
}

// Do vt = v1-s
void sseVsC(wsVecCst Ra_v1, wsRealCst R_s, wsVec Ra_vt, wsUintCst N_sz,
			wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint N_med = 0;

	if (NA(R_mul)) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(R_s);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseSub(*sse_v1, sse_S);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v1[i] - R_s;
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(R_s);
		sse_t sse_M = sseSet(R_mul);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseMul(sseSub(*sse_v1, sse_S), sse_M);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = (Ra_v1[i] - R_s) * R_mul;
	}
}

// Do vt = v*s
void sseVpC(wsVecCst Ra_v, wsRealCst R_s, wsVec Ra_vt, wsUintCst N_sz,
			wsRealPtr Rp_sqSum/*=NULL*/)
{
	wsUint	i;
	wsUint	N_med = 0;

	if (Rp_sqSum) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(R_s);
		sse_t sse_ss = sseSet(0.0);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);

			*sse_vt = sseMul(*sse_v, sse_S);
			sse_ss = sseAdd(sse_ss, sseMul(*sse_vt, *sse_vt));
		}
		sseSum(sse_ss, *Rp_sqSum);
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++) {
			Ra_vt[i] = Ra_v[i]*R_s;
			*Rp_sqSum += SQR(Ra_vt[i]);
		}
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(R_s);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);

			*sse_vt = sseMul(*sse_v, sse_S);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v[i]*R_s;
	}
}

// Do vt = v*s
void sseVpC(wsVecCst Ra_v, wsRealCst R_s, wsVec Ra_vt, wsUintCst N_sz,
			wsRealCst R_add)
{
	wsUint	i;
	wsUint	N_med = 0;

	if (NA(R_add)) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(R_s);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);

			*sse_vt = sseMul(*sse_v, sse_S);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v[i]*R_s;
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(R_s);
		sse_t sse_A = sseSet(R_add);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);

			*sse_vt = sseAdd(sseMul(*sse_v, sse_S), sse_A);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v[i]*R_s + R_add;
	}
}

// Do vt = v1*v2 or vt = v1*v2 * c
void sseVpV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
			wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint N_med = 0;

	if (NA(R_mul)) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseMul(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v1[i] * Ra_v2[i];
	} else {
#ifdef USE_SSE
		sse_t sse_C = sseSet(R_mul);
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseMul(sse_C, sseMul(*sse_v1, *sse_v2));
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v1[i] * Ra_v2[i] * R_mul;
	}
}

// Do vt = vt = v1*sqrt(v2)
void sseVpVsqrt(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz)
{
	wsUint N_med = getMed(N_sz);

#ifdef USE_SSE
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
		sse_t *sse_vt = (sse_t *)(Ra_vt + i);
		*sse_vt = sseMul(*sse_v1, sseSqrt(*sse_v2));
	}
	/* Do rest part */
	for (wsUint i=N_med ; i<N_sz ; i++)
		Ra_vt[i] = Ra_v1[i] * sqrt(Ra_v2[i]);
#else
	for (wsUint i=0 ; i<N_sz ; i++)
		Ra_vt[i] = Ra_v1[i] * sqrt(Ra_v2[i]);
#endif
}

// Do vt = vt = v1/v2
void sseVdV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
			wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsUint N_med = 0;

	if (NA(R_mul)) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseDiv(*sse_v1, *sse_v2);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v1[i] / Ra_v2[i];
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_M = sseSet(R_mul);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
			sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);
			*sse_vt = sseMul(sseDiv(*sse_v1, *sse_v2), sse_M);
		}
#endif
		/* Do rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = Ra_v1[i] / Ra_v2[i] * R_mul;
	}
}

// Do vt += v*s
void sseVpCadd(wsVecCst Ra_v, wsRealCst R_s, wsVec Ra_vt, wsUintCst N_sz,
			   wsRealCst R_add/*=WISARD_NAN*/)
{
	wsUint	i;
	wsUint	N_med = 0;

	if (NA(R_add)) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(R_s);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);

			*sse_vt = sseAdd(*sse_vt, sseMul(*sse_v, sse_S));
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_vt[i] += Ra_v[i]*R_s;
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(R_s);
		sse_t sse_A = sseSet(R_add);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);

			*sse_vt = sseAdd(*sse_vt, sseAdd(sseMul(*sse_v, sse_S), sse_A));
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_vt[i] += Ra_v[i]*R_s + R_add;
	}
}

// Do vt = vt = |v1|/v2
void sseAVdV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz)
{
	wsUint N_med = 0;

#ifdef USE_SSE
	N_med = getMed(N_sz);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
		sse_t *sse_vt = (sse_t *)(Ra_vt + i);
		*sse_vt = sseDiv(sseAbs(*sse_v1), *sse_v2);
	}
#endif
	/* Do rest part */
	for (wsUint i=N_med ; i<N_sz ; i++)
		Ra_vt[i] = fabs(Ra_v1[i]) / Ra_v2[i];
}

wsVec sseVp1p(wsVecCst Ra_v, wsUintCst N_sz, wsVec Ra_ret, char B_sq/*=0*/)
{
	wsUint	N_med	= 0;

	if (B_sq) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_1 = sseSet(1.0);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_s = (sse_t *)(Ra_v + i);
			sse_t *sse_r = (sse_t *)(Ra_ret + i);
			sse_t sse_tmp = sseSub(sse_1, *sse_s);
			*sse_r = sseMul(sseMul(*sse_s, *sse_s),
				sseMul(sse_tmp, sse_tmp));
		}
#endif
		/* Rest part */
		for (wsUint i=N_med ; i<N_sz ; i++) {
			wsReal R_1 = SQR(Ra_v[i]);
			wsReal R_2 = W1 - Ra_v[i];
			Ra_ret[i] = R_1 * SQR(R_2);
		}
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_1 = sseSet(1.0);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_s = (sse_t *)(Ra_v + i);
			sse_t *sse_r = (sse_t *)(Ra_ret + i);
			*sse_r = sseMul(*sse_s, sseSub(sse_1, *sse_s));
		}
#endif
		/* Rest part */
		for (wsUint i=N_med ; i<N_sz ; i++)
			Ra_ret[i] = Ra_v[i] * (W1 - Ra_v[i]);
	}

	return Ra_ret;
}

wsMat sseVtV(wsVecCst Ra_V, wsUintCst N_sz, wsRealCst R_mul/*=WISARD_NAN*/)
{
	wsMat	Ra_ret	= sseMatrix(N_sz, N_sz);
	wsUint	N_med	= 0;
#ifdef USE_SSE
	N_med = getMed(N_sz);
#endif

	for (wsUint i=0 ; i<N_sz ; i++) {
		wsVec	Ra_R	= Ra_ret[i];
		wsReal	R_q		= NA(R_mul) ? Ra_V[i] : Ra_V[i] * R_mul;
#ifdef USE_SSE
		sse_t sse_Vi = sseSet(R_q);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_R = (sse_t *)(Ra_R + j);
			sse_t *sse_V = (sse_t *)(Ra_V + j);
			*sse_R = sseMul(sse_Vi, *sse_V);
		}	
#endif
		for (wsUint j=N_med ; j<N_sz ; j++)
			Ra_R[j] = R_q * Ra_V[j];
	}

	return Ra_ret;
}

wsMat sseVtV(wsVecCst Ra_v1, wsUintCst N_sz, wsVecCst Ra_v2,
			 wsRealCst R_mul/*=WISARD_NAN*/)
{
	return sseVtV(Ra_v1, N_sz, Ra_v2, N_sz, R_mul);
}

wsMat sseVtV(wsVecCst Ra_v1, wsUintCst N_s1, wsVecCst Ra_v2, wsUintCst N_s2,
			 wsRealCst R_mul/*=WISRAD_NAN*/)
{
	wsMat	Ra_ret	= sseMatrix(N_s1, N_s2);
	wsUint	N_med	= 0;
#ifdef USE_SSE
	N_med = getMed(N_s2);
#endif

	if (NA(R_mul)) for (wsUint i=0 ; i<N_s1 ; i++) {
		wsVec	Ra_R = Ra_ret[i];
#ifdef USE_SSE
		sse_t sse_Vi = sseSet(Ra_v1[i]);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_R = (sse_t *)(Ra_R + j);
			sse_t *sse_V = (sse_t *)(Ra_v2 + j);
			*sse_R = sseMul(sse_Vi, *sse_V);
		}	
#endif
		for (wsUint j=N_med ; j<N_s2 ; j++)
			Ra_R[j] = Ra_v1[i] * Ra_v2[j];
	} else for (wsUint i=0 ; i<N_s1 ; i++) {
		wsUint	N_med	= 0;
		wsVec	Ra_R	= Ra_ret[i];
#ifdef USE_SSE
		sse_t sse_Vi	= sseSet(Ra_v1[i]);
		sse_t sse_M		= sseSet(R_mul);
		sse_t sse_ViM	= sseMul(sse_Vi, sse_M);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_R = (sse_t *)(Ra_R + j);
			sse_t *sse_V = (sse_t *)(Ra_v2 + j);
			*sse_R = sseMul(sse_ViM, *sse_V);
		}	
#endif
		for (wsUint j=N_med ; j<N_s2 ; j++)
			Ra_R[j] = Ra_v1[i] * Ra_v2[j] * R_mul;
	}

	return Ra_ret;
}

// Do vt = 1/v
void sseVinv(wsVec Ra_v, wsUintCst N_sz, wsVec Ra_vt/*=NULL*/)
{
	wsUint	i;
	wsUint	N_med = 0;

	if (Ra_vt) {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(1.0);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v + i);
			sse_t *sse_vt = (sse_t *)(Ra_vt + i);

			*sse_vt = sseDiv(sse_S, *sse_v);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = W1 / Ra_v[i];
	} else {
#ifdef USE_SSE
		N_med = getMed(N_sz);
		sse_t sse_S = sseSet(1.0);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_v = (sse_t *)(Ra_v + i);

			*sse_v = sseDiv(sse_S, *sse_v);
		}
#endif
		/* Do rest part */
		for (i=N_med ; i<N_sz ; i++)
			Ra_v[i] = W1 / Ra_v[i];
	}
}

void sseVpermSelf(wsVec Ra_v, wsUintCst N_sz, wsUintCst* Na_seq)
{
	wsVec Ra_tmp = sseVector(N_sz);

	/* Permute vector */
	LOOP (i, N_sz) Ra_tmp[i] = Ra_v[Na_seq[i]];
	/* Copy permuted result */
	memcpy(Ra_v, Ra_tmp, sizeof(wsReal)*N_sz);
}

wsVec vectorRepeat(wsVecCst Ra_vec, vInt& V_cntRep, wsUint N_sz)
{
	wsVec	Ra_ret	= sseVector(N_sz);
	wsUint	I		= 0;
	wsUint	J		= 0;
	FOREACHDO (vInt_it, V_cntRep, i, J++) {
		wsReal	R_v		= Ra_vec[J];
		wsUint	N_rep	= (wsUint)*i;
		LOOP (j, N_rep) {
			if (I>=N_sz) halt("SYSERR: This function have an error");
			Ra_ret[I++] = R_v;
		}
	}
	if (I != N_sz)
		halt("SYSERR: Repeated size[%d] != expected size[%d]", I, N_sz);
	return Ra_ret;
}

wsVec vectorRepeat(wsVecCst Ra_vec, mvStr& V_cntRep, wsUint N_sz)
{
	wsVec	Ra_ret	= sseVector(N_sz);
	wsUint	I		= 0;
	wsUint	J		= 0;
	FOREACHDO(mvStr_it, V_cntRep, i, J++) {
		wsReal	R_v		= Ra_vec[J];
		wsUint	N_rep	= (wsUint)(i->second.size());
		LOOP(j, N_rep) {
			if (I>=N_sz) halt("SYSERR: This function have an error");
			Ra_ret[I++] = R_v;
		}
	}
	if (I != N_sz)
		halt("SYSERR: Repeated size[%d] != expected size[%d]", I, N_sz);
	return Ra_ret;
}

wsVec vectorCummax(wsVecCst Ra_vec, wsUint N_sz, xRealSort* Xa_rs/*=NULL*/)
{
	wsVec Ra_ret = sseVector(N_sz);

	if (Xa_rs == NULL) {
		wsReal R_max = Ra_vec[0];
		Ra_ret[0] = R_max;
		for (wsUint i=1; i < N_sz; i++) {
			if (R_max < Ra_vec[i]) R_max = Ra_vec[i];
			Ra_ret[i] = R_max;
		}
	} else {
		wsReal R_max = Ra_vec[Xa_rs[0].i];
		Ra_ret[0] = R_max;
		for (wsUint i=1; i < N_sz; i++) {
			wsRealCst V = Ra_vec[Xa_rs[i].i];
			if (R_max < V) R_max = V;
			Ra_ret[i] = R_max;
		}
	}
	return Ra_ret;
}

wsVec vectorCummax2(wsVecCst Ra_vec, wsUint N_sz, xRealSort* Xa_rs/*=NULL*/)
{
	wsVec Ra_ret = sseVector(N_sz);

	if (Xa_rs == NULL) {
		wsReal R_max = Ra_vec[0];
		Ra_ret[0] = R_max;
		for (wsUint i=1; i < N_sz; i++) {
			if (R_max < Ra_vec[i]) R_max = Ra_vec[i];
			Ra_ret[i] = R_max;
		}
	} else {
		wsReal R_max = Ra_vec[Xa_rs[0].i];
		Ra_ret[Xa_rs[0].i] = R_max;
		for (wsUint i=1; i < N_sz; i++) {
			wsUintCst I = Xa_rs[i].i;
			wsRealCst V = Ra_vec[I];
			if (R_max < V) R_max = V;
			Ra_ret[I] = R_max;
		}
	}
	return Ra_ret;
}

wsReal sseVse(wsVec Ra_t, wsUint N_sz, wsReal *Rp_mean/*=NULL*/)
{
	wsReal R_sum = W0;
	wsReal R_mean = W0;
	R_mean = sseVsum(Ra_t, N_sz) / (wsReal)N_sz;
	if (Rp_mean)
		*Rp_mean = R_mean;
	wsUint N_med = 0;
#ifdef USE_SSE
	sse_t sse_s = sseSet(W0);
	sse_t sse_m = sseSet(R_mean);
	for (wsUint i=0; i < N_med; i+=sseJmp) {
		sse_t* sse_tt = (sse_t *)(Ra_t + i);
		sse_t sse_tmp = sseSub(*sse_tt, sse_m);
		sse_s = sseAdd(sse_s, sseMul(sse_tmp, sse_tmp));
	}
	sseSum(sse_s, R_sum);
#endif
	for (wsUint i=N_med; i < N_sz; i++) {
		wsReal R_tmp = Ra_t[i] - R_mean;
		R_sum += SQR(R_tmp);
	}
	return sqrt(R_sum / (wsReal)(N_sz - 1));
}

void sseVstdize(wsVec Ra_t, wsUint N_sz)
{
	wsReal R_mean = W0;
	wsReal R_se = sseVse(Ra_t, N_sz, &R_mean);
	sseVsC(Ra_t, R_mean, Ra_t, N_sz, W1 / R_se);
}

/* cSpsVector implementation */

cSpsVector::cSpsVector(wsVec Ra_val, wsUint N_inpSz, wsUint N_inpElem,
	wsUint* Na_inpIdx)
{
	Ra_buf	= Ra_val;
	N_sz	= N_inpSz;
	N_elem	= N_inpElem;
	Na_idx	= Na_inpIdx;
}

wsRealCst cSpsVector::sum(wsReal *Rp_sqSum/*=NULL*/)
{
	return sseVsum(get(), N_elem, Rp_sqSum);
}

wsRealCst cSpsVector::sum(wsVec Ra_vec, wsUint N_sz)
{
	ASSERT_SIZE(size(), N_sz, "vector-vector sum");
	wsReal	R_ret	= W0;
	wsVecCst	Ra_v	= get();
	LOOP (i, N_elem)
		R_ret += Ra_vec[Na_idx[i]] * Ra_v[i];
	return R_ret;
}

wsRealCst cSpsVector::sum(cVector &V)
{
	return sum(V.get(), V.size());
}

void cSpsVector::set0()
{
	/* Just set to 0 by the size of actual elements */
	sseVset0(get(), N_elem);
}

cSpsVector cSpsVector::p1p()
{
	wsVec	Ra_ret	= sseVector(N_elem);
	wsUint*	Na_ret	= NULL;
	wsAlloc(Na_ret, wsUint, N_elem);
	memcpy(Na_ret, Na_idx, sizeof(wsUint)*N_elem);
	return cSpsVector(sseVp1p(get(), N_elem, Ra_ret), size(), N_elem,
		Na_ret);
}

wsReal* cSpsVector::get()
{
	return Ra_buf;
}

} // End namespace ONETOOL
