#include "analyses/ld.h"
#include "input/stream.h"
#include "utils/data.h"
//#include "emai.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cLdAnalysis::cLdAnalysis(cIO *Cp_inpIO, wsUint *Np_sset/*=NULL*/,
	wsUint N_size/*=0*/) : cAnalysis(Cp_inpIO)
{
	/* Check */
	if (Np_sset == NULL && N_size != 0)
		halt("Subset size is given but subset indices is not given!");
	if (!OPT_ENABLED(ld))
		return;
	Na_sset		= Np_sset;
	N_set		= N_size;
	Cp_anaMM	= NULL;
}

cLdAnalysis::cLdAnalysis(cIO *Cp_inpIO, cVariantMap *Cp_inpAnaMM) :
	cAnalysis(Cp_inpIO)
{
	Cp_anaMM	= Cp_inpAnaMM;
	Na_sset		= NULL;
	N_set		= 0;

	/* --ldbin & --ldsize are mutually exclusive */
	if (IS_ASSIGNED(ldbin) && IS_ASSIGNED(ldsize))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--ldbin", "--ldsize");

	/* --ldvariant makes loop to only some variants */
	if (IS_ASSIGNED(ldvar)) {
		/* Trying to open */
		cStrFile C_lm(OPT_STRING(ldvar), "LD computing variants list", 1);
		char **Sa_lms = NULL;
		wsUint N_mk = 0;

		if (C_lm.isFailed())
			Sa_lms = loadStringValues(OPT_STRING(ldvar), &N_mk);
		else {
			/* Get the number */
			char *S_buf = NULL, *a = NULL;
			wsAlloc(S_buf, char, 512);
			for ( ; C_lm.gets(S_buf, 512) ; N_mk++) {
				getString(&S_buf, &a);
				if (!S_buf[0]) break;
			}
			C_lm.rewind();

			wsAlloc(Sa_lms, char*, N_mk);
			for (wsUint i=0 ; C_lm.gets(S_buf, 512) ; i++) {
				getString(&S_buf, &a);
				Sa_lms[i] = strdup(S_buf);
				if (!S_buf[0]) break;
			}
			C_lm.close();

			DEALLOC(S_buf);
		}

		/* Find which variants are available */
		/* Build variant map */
		vVariant&	Xa_snp	= Cp_IO->getVariant();
		mDataIdx	Xm_mmap[MAX_NCHR+1];
		wsUint		j = 0;
		FOREACHDO (vVariant_it, Xa_snp, i, j++) {
			if (i->chr > 0 && (wsUint)i->chr <= NCHR_SPECIES)
				Xm_mmap[i->chr].insert(make_pair((string)i->name, j));
			else
				Xm_mmap[0].insert(make_pair((string)i->name, j));
		}

		/* For listed variants, get indices */
		wsAlloc(Na_sset, wsUint, N_mk);
		N_set = N_mk;
		for (wsUint i=0 ; i<N_mk ; i++) {
			int N_insIdx = -1;
			/* find */
			for (wsUint N_chr=0 ; N_chr<=NCHR_SPECIES ; N_chr++) {
				mDataIdx_it X_find = Xm_mmap[N_chr].find(Sa_lms[i]);
				/* If find, record index */
				if (X_find != Xm_mmap[N_chr].end()) {
					N_insIdx = X_find->second;
					break;
				}
			}
			/* Ignore this line if no match found */
			if (N_insIdx == -1)
				halt("Variant [%s] is not exists in the dataset\n",
					Sa_lms[i]);
			Na_sset[i] = N_insIdx;
		}
	}

	/* If --ldbin */
	if (IS_ASSIGNED(ldbin) || IS_ASSIGNED(ldsize)) {
		wsUint		*Na_chr	= NULL;
		wsUint**	Xa_chr	= Cp_anaMM->getVariantPosMap(&Na_chr);
		wsUint		N_sz	= OPT_NUMBER(ldsize);
		vVariant&	Xa_snp	= Cp_IO->getVariant();
		wsUint		N_SNP	= (wsUint)Xa_snp.size();

		/* For all variants, find adjacent variants within specific range */
		wsUint j = 0;
		FOREACHDO (vVariant_it, Xa_snp, i, j++) {
			Na_SE.push_back(vInt());
			wsUint	P = i->pos;
			wsUint*	Q = Xa_chr[i->chr];
			wsUint	N = Na_chr[i->chr];

			/* Find out its position on variantmap */
			wsUint s = 0, e = N;
			while (s < e) {
				wsUint	m = (s+e)>>1;
				xVariant&	S = Xa_snp[Q[m]];
				if (S.pos < P) s = m;
				else if (S.pos > P) e = m;
				else s = e = m;
			}
			xVariant&	C = Xa_snp[Q[s]];
			if (C.pos != P) halt("SYSERR : Can't find given pos [%d]", P);

			/* Do - and + */
			if (IS_ASSIGNED(ldsize)) {
				//wsUint N_idxS = !s ? s : s-1, N_idxE = e==(N_SNP-1) ? e : e+1;
				char N_setS = !s ? 1 : 0, N_setE = e==(N_SNP-1) ? 1 : 0;
				while (!N_setS && !N_setE) {
					xVariant&	S = Xa_snp[Q[s]];
					xVariant&	E = Xa_snp[Q[e]];
					if (!N_setS) {
						if ((P - S.pos) < N_sz) {
							if (--s == 0) N_setS = 1;
						} else {
							s++;
							N_setS = 1;
						}
					}

					if (!N_setE) {
						if ((E.pos - P) < N_sz) {
							if (++e == (N-1)) N_setE = 1;
						} else {
							e--;
							N_setE = 1;
						}
					}
				}
			} else {
				int pp = (int)s - OPT_NUMBER(ldbin);

				s = pp<0 ? 0 : (wsUint)pp;
				e += OPT_NUMBER(ldbin);
				e = e>=N ? N-1 : e;
			}

			/* Insert positions */
			for (wsUint k=s ; k<=e ; k++)
				Na_SE[j].push_back(Q[k]);
		}

		if (IS_ASSIGNED(ldbin))
			LOG("Block-wise LD [bin size %d] is computed\n", OPT_NUMBER(ldbin));
		if (IS_ASSIGNED(ldsize))
			LOG("Block-wise LD [bin range +-%d] is computed\n", OPT_NUMBER(ldsize));
	} else
		LOG("All-pairwise LD is computed\n");
}

cLdAnalysis::~cLdAnalysis()
{
	Na_SE.clear();
}

wsReal calcR(char **Na_data, wsUint N_sz, wsUint i, wsUint j,
	const char *Ba_isFounder)
{
	wsUint sumXY	= 0;
	wsUint sumX		= 0;
	wsUint sumX2	= 0;
	wsUint sumY		= 0;
	wsUint sumY2	= 0;
	wsUint N		= 0;
	//wsReal R_ret	= W0;

	if (Ba_isFounder == NULL) for (wsUint k=0 ; k<N_sz ; k++) {
		if (isAvailable(Na_data[k][i]) && isAvailable(Na_data[k][j])) {
			sumX	+= Na_data[k][i];
			sumX2	+= Na_data[k][i] * Na_data[k][i];
			sumY	+= Na_data[k][j];
			sumY2	+= Na_data[k][j] * Na_data[k][j];
			sumXY	+= Na_data[k][i] * Na_data[k][j];
			N++;
		}
	}
	else for (wsUint k=0 ; k<N_sz ; k++) {
		if (!Ba_isFounder[k]) continue;
		if (isAvailable(Na_data[k][i]) && isAvailable(Na_data[k][j])) {
			sumX	+= Na_data[k][i];
			sumX2	+= Na_data[k][i] * Na_data[k][i];
			sumY	+= Na_data[k][j];
			sumY2	+= Na_data[k][j] * Na_data[k][j];
			sumXY	+= Na_data[k][i] * Na_data[k][j];
			N++;
		}
	}

	return ((wsReal)N*(wsReal)sumXY - (wsReal)sumX*(wsReal)sumY)
		/ (sqrt((wsReal)N*(wsReal)sumX2 - (wsReal)sumX*(wsReal)sumX) *
		sqrt((wsReal)N*(wsReal)sumY2 - (wsReal)sumY*(wsReal)sumY));
}

int gsl_poly_solve_cubic(double a, double b, double c, double *x0, double *x1, double *x2)
{
	return 0;
}

int phase(const int N, const char *x, const char *y, 
		  const int *diploid, double *hapfreq, double *margins, double *LLR) {
	int T[4]={0, 0, 0, 0}, G[9]={0, 0, 0, 0, 0, 0, 0, 0, 0};

	/* Defaults (monomorphic SNP) */
	for (int i=0; i<N; i++) {
		char xi = x[i], yi = y[i];

		/* Pass if missing */
		if (isMissing(xi) || isMissing(yi)) continue;
		int xyi = xi*3 + yi;
//		if (!diploid || diploid[i]) 
			G[xyi]++;
// 		else 
// 			switch (xyi) {
// 			case 0: T[0]++; break;
// 			case 2: T[1]++; break;
// 			case 6: T[2]++; break;
// 			case 8: T[3]++; break;
// 			default: return 2; /* Heterozygous haploid genotype on X */
// 		}
	}
	T[0]+= 2*G[0]+G[1]+G[3];
	T[1]+= 2*G[2]+G[1]+G[5];
	T[2]+= 2*G[6]+G[3]+G[7];
	T[3]+= 2*G[8]+G[7]+G[5];
	int Dh = G[4];
	int Nph = T[0]+T[1]+T[2]+T[3];
	if (!Nph) {
		return 2;
	}
	double Nhap = Nph +2*Dh;
	double E1[4];
	/* Allele frequencies */
	margins[0] = (double)(T[0]+T[1]+Dh)/Nhap;
	margins[1] = 1.0-margins[0];
	margins[2] = (double)(T[0]+T[2]+Dh)/Nhap;
	margins[3] = 1.0-margins[2];
	if (!(margins[0]&&margins[1]&&margins[2]&&margins[3])) {
		return 1; /* At least one SNP is monomorphic */
	}
	/* Solve cubic equation */
	double  roots[3];
	int nroot = 0;
	double w1 = (double)(T[0]+T[3]-T[1]-T[2]);
	double w2 = (double)(T[0]*T[3]);
	double w3 = (double)(T[1]*T[2]);
	if (!Dh) {
		nroot = 1;
		roots[0] = w2/(w2+w3);
	}
	else { 
		w1 /= Dh;
		double Dh2 = Dh*Dh;
		w2 /= Dh2;
		w3 /= Dh2;
		double a = (w1-3.0)/2.0;
		double b = (w2 + w3 - w1 + 1.0)/2.0;
		double c = -w2/2.0;
		nroot = gsl_poly_solve_cubic(a, b, c, roots, roots+1, roots+2);
		if (!nroot)
			return 3;
	}
	double llh = 0.0, p = -1.0;
	if (LLR || (nroot>1)) {
		for (int i=0; i<nroot; i++) {
			double pi = roots[i];
			if (pi<0.0)
				pi = 0.0;
			if (pi>1.0)
				pi = 1.0;
			double Dhpi = Dh*pi;
			double Dhqi = Dh-Dhpi;
			E1[0] = ((double)T[0] + Dhpi)/Nhap;
			E1[1] = ((double)T[1] + Dhqi)/Nhap;
			E1[2] = ((double)T[2] + Dhqi)/Nhap;
			E1[3] = ((double)T[3] + Dhpi)/Nhap;
			/* Probability of double het (divided by 2)*/
			double EDh1 = (E1[0]*E1[3] + E1[1]*E1[2])/2.0;
			double llhi = 0.0;
			if (Dh) 
				llhi += Dh*log(EDh1);
			for (int j=0; j<4; j++) 
				if (T[j]) llhi += T[j]*log(E1[j]);
			if (p<0.0 || (llhi>llh)) {
				llh = llhi;
				p = pi;
			}
		}
	}
	else 
		p = roots[0];
	if (p<0) {
		return 3;
	}

	/* Normal return */

	double Dhp = Dh*p;
	double Dhq = Dh-Dhp;
	hapfreq[0] = ((double)T[0] + Dhp)/Nhap;
	hapfreq[1] = ((double)T[1] + Dhq)/Nhap;
	hapfreq[2] = ((double)T[2] + Dhq)/Nhap;
	hapfreq[3] = ((double)T[3] + Dhp)/Nhap;
	if (LLR) {
		double llh0=0.0;
		for (int i=0; i<4; i++)
			llh0 += margins[i]*log(margins[i]);
		*LLR = llh - Nhap*llh0;
	}
	return 0;
}

void set_arrays(const double *hapfreqs, const double *margins, double LLR, 
				double **arrays, int ij)
{
	/* LLR */
	if (arrays[0]) (arrays[0])[ij] = LLR;
	/* OR */
	double ad = hapfreqs[0]*hapfreqs[3];
	double bc = hapfreqs[1]*hapfreqs[2];
	//double OR =  ad/bc;
	if (arrays[1]) (arrays[1])[ij] = ad/bc;
	/* Yules Q */
	if (arrays[2]) (arrays[2])[ij] = (ad-bc)/(ad+bc);
	/* Covariance */
	double covar = hapfreqs[0] - margins[0]*margins[2];
	if (arrays[3]) (arrays[3])[ij] = covar;
	/* D-prime */
	if (arrays[4]) {
		if (covar>0) {
			double P1Q2 = margins[0]*margins[3];
			double P2Q1 = margins[1]*margins[2];
			(arrays[4])[ij] = covar/(P1Q2<P2Q1? P1Q2: P2Q1);
		}
		else {
			double P1Q1 = margins[0]*margins[2];
			double P2Q2 = margins[1]*margins[3];
			(arrays[4])[ij] = -covar/(P1Q1<P2Q2? P1Q1: P2Q2);
		}
	}
	/* R-squared */
	double mprod = margins[0]*margins[1]*margins[2]*margins[3];
	if (arrays[5]) (arrays[5])[ij] = covar*covar/mprod;
	/* R */
	if (arrays[6]) (arrays[6])[ij] = covar/sqrt(mprod);
}

/*
def get_covD(Y,Z):
    """
    get_covD estimates pA, pB, and D w/o info on gametic phase.
    Uses the method of Rogers and Huff 2008.
    """
    pA, v0, pB, v1, cov = bivmom(Y,Z)
    pA = 0.5*pA
    pB = 0.5*pB
    qA=1-pA
    qB=1-pB
    two_1pf = sqrt((v0*v1)/(pA*qA*pB*qB)) # estimates 2(1+f)
    D = cov/two_1pf
    return pA, pB, D
*/
void get_covD(char *Y, char *Z, double *pD, double *ppA, double *ppB,
	wsUint N_cnt[3][3], wsUint N_sz)
{
	double pA = 0, pB = 0, v0 = 0, v1 = 0, cov = 0;
	int n = 0;

	for (wsUint i=0 ; i<N_sz ; i++) {
		char __X = Y[i];
		char _Y = Z[i];
		if (isMissing(__X) || isMissing(_Y)) continue;

		N_cnt[(wsUint)__X][(wsUint)_Y]++;
		double x = (double)Y[i];
		double y = (double)Z[i];

		pA += x;
		pB += y;
		v0 += x*x;
		v1 += y*y;
		cov += x*y;
		n++;
	}
	double nn = (double)n;
	pA /= nn;
	pB /= nn;
	v0 /= nn;
	v1 /= nn;
	cov /= nn;

	cov -= pA * pB;
	v0 -= pA * pA;
	v1 -= pB * pB;

	pA = 0.5*pA;
	pB = 0.5*pB;
	double qA=1-pA;
	double qB=1-pB;
	double two_1pf = sqrt((v0*v1)/(pA*qA*pB*qB)); /* estimates 2(1+f) */

	*pD = cov/two_1pf;
	*ppA = pA;
	*ppB = pB;
}

/* Rogers-Huff method */
/*
def rhesem_r(Y,Z, max_itr=1000):
    global tol
    # RH step
    pA, pB, D = get_covD(Y,Z)
    qA=1.0-pA
    qB=1.0-pB
    h = [pA*pB + D, pA*qB-D, qA*pB-D, qA*qB+D]

    # ES step
    h = rhesem(Y,Z, h, max_itr)
    pA = h[0]+h[1]
    pB = h[0]+h[2]
    qA = 1.0 - pA
    qB = 1.0 - pB
    r = h[0]*h[3] - h[1]*h[2]
    r /= sqrt(pA*qA*pB*qB)
    return r
*/
double getHfreq_RogersHuff(char *Y, char *Z, wsUint N_sz, double *Rp_Dprime, wsUint N_itr=1000)
{
	double D, pA, pB;
	double R_tol = 1e-10;
	wsUint x[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };

	/* R-H step */
	double H[4];
	get_covD(Y, Z, &D, &pA, &pB, x, N_sz);
	double R_qA = 1.0 - pA;
	double R_qB = 1.0 - pB;
	H[0] = pA*pB + D;
	H[1] = pA*R_qB - D;
	H[2] = R_qA*pB - D;
	H[3] = R_qA*R_qB + D;

	/* Do E-S step */
	double R_delta = 1.0;
	for (wsUint i=0 ; i<N_itr && R_delta>R_tol ; i++) {
		/*
		def rhesem_step(h, x):

# g is a 4X4 matrix of genotype frequencies. g[0][3]
# is the frequency of the genotype that combines gamete 0 (AB)
# with gamete 3 (ab).
			g = [[None,None,None,None],[None,None,None,None],
			[None,None,None,None],[None,None,None,None]] */

		double g[4][4] = { {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} };
		for (wsUint j=0 ; j<4 ; j++) {
			for (wsUint k=0 ; k<j ; k++)
				g[j][k] = 2.0 * H[j] * H[k];
			g[j][j] = SQR(H[j]);
		}
		/*
		for i in range(4):
		g[i][i] = h[i]*h[i]
		for j in range(i):
		g[i][j] = 2*h[i]*h[j]

# p is a 3X3 matrix of phenotype frequencies, recoded as
# described for the input matrix x.
		p = [[None,None,None],[None,None,None],[None,None,None]]
		*/
		double p[3][3];

		p[0][0] = g[0][0];
		p[0][1] = g[1][0];
		p[0][2] = g[1][1];

		p[1][0] = g[2][0];
		p[1][1] = g[3][0]+g[2][1];
		p[1][2] = g[3][1];

		p[2][0] = g[2][2];
		p[2][1] = g[3][2];
		p[2][2] = g[3][3];

		double hh[4];
		hh[0] = 2*x[0][0] + x[0][1] + x[1][0] + x[1][1]*g[3][0]/p[1][1];
		hh[1] = x[0][1] + 2*x[0][2] + x[1][1]*g[2][1]/p[1][1] + x[1][2];
		hh[2] = x[1][0] + x[1][1]*g[2][1]/p[1][1] + 2*x[2][0] + x[2][1];
		hh[3] = x[1][1]*g[3][0]/p[1][1] + x[1][2] + x[2][1] + 2*x[2][2];

		// haploid sample size
		double n = hh[0] + hh[1] + hh[2] + hh[3];

		// convert gamete counts counts to relative frequencies
		// for i in range(4):
		for (wsUint j=0 ; j<4 ; j++)
			hh[j] /= n;

		/* Update delta */
		R_delta = fabs(hh[0] - H[0]) + fabs(hh[1] - H[1]) +
			fabs(hh[2] - H[2]) + fabs(hh[3] - H[3]);
		memcpy(H, hh, sizeof(double)*4);
	}

	pA = H[0]+H[1];
	pB = H[0]+H[2];
	R_qA = 1.0 - pA;
	R_qB = 1.0 - pB;
	/* R */
	double r = (H[0]*H[3] - H[1]*H[2]) / sqrt(pA*R_qA*pB*R_qB);
	/* Dprime */
	if (Rp_Dprime) {
		double Q = H[0] - pA*pB;
		/* D-prime */
		if (Q > 0) {
			double P1Q2 = pA * R_qB;
			double P2Q1 = R_qA * pB;
			*Rp_Dprime = Q / (P1Q2<P2Q1? P1Q2: P2Q1);
		} else {
			double P1Q1 = pA * R_qA;
			double P2Q2 = R_qB * pB;
			*Rp_Dprime = -Q / (P1Q1<P2Q2? P1Q1: P2Q2);
		}
	}

	return R_delta > R_tol ? WISARD_NAN : r;
}

/*
def rhesem(Y,Z, h=[0.25,0.25,0.25,0.25], max_itr=1000):
    global tol
    x = [[0,0,0],[0,0,0],[0,0,0]]
    for y,z in zip(Y,Z):
        x[y][z] += 1
    for itr in xrange(max_itr):
        hh = rhesem_step(h, x)
        dh = 0.0
        for u,v in zip(h, hh):
            dh += abs(u-v)
        if dh <= tol:
            break
        h = hh
    if dh > tol:
        raise ConvergenceError
    return hh
*/

// int phase(const int N, const unsigned char *x, const unsigned char *y, 
// 	double *hapfreq, double *margins, double *LLR)
// {
// 	int T[4]={0, 0, 0, 0}, G[9]={0, 0, 0, 0, 0, 0, 0, 0, 0};
// 
// 	/* Defaults (monomorphic SNP) */
// 	for (int i=0; i<N; i++) {
// 		unsigned char xi = x[i], yi = y[i];
// 
// 		/* Pass if missing */
// 		if (isMissing(xi) || isMissing(yi)) continue;
// 		int xyi = xi*3 + yi;
// //		if (!diploid || diploid[i]) 
// 			G[xyi]++;
// // 		else switch (xyi) {
// // 			case 0: T[0]++; break;
// // 			case 2: T[1]++; break;
// // 			case 6: T[2]++; break;
// // 			case 8: T[3]++; break;
// // 			default: return 2; /* Heterozygous haploid genotype on X */
// // 		}
// 	}
// 	T[0]+= 2*G[0]+G[1]+G[3];
// 	T[1]+= 2*G[2]+G[1]+G[5];
// 	T[2]+= 2*G[6]+G[3]+G[7];
// 	T[3]+= 2*G[8]+G[7]+G[5];
// 	int Dh = G[4];
// 	int Nph = T[0]+T[1]+T[2]+T[3];
// 	if (!Nph) {
// 		return 2;
// 	}
// 	double Nhap = Nph +2*Dh;
// 	double E1[4];
// 
// 	/* Allele frequencies */
// 	margins[0] = (double)(T[0]+T[1]+Dh)/Nhap;
// 	margins[1] = 1.0-margins[0];
// 	margins[2] = (double)(T[0]+T[2]+Dh)/Nhap;
// 	margins[3] = 1.0-margins[2];
// 	if (!(margins[0]&&margins[1]&&margins[2]&&margins[3])) {
// 		return 1; /* At least one SNP is monomorphic */
// 	}
// 
// 	/* Solve cubic equation */
// 	double  roots[3];
// 	int nroot = 0;
// 	double w1 = (double)(T[0]+T[3]-T[1]-T[2]);
// 	double w2 = (double)(T[0]*T[3]);
// 	double w3 = (double)(T[1]*T[2]);
// 	if (!Dh) {
// 		nroot = 1;
// 		roots[0] = w2/(w2+w3);
// 	}
// 	else { 
// 		w1 /= Dh;
// 		double Dh2 = Dh*Dh;
// 		w2 /= Dh2;
// 		w3 /= Dh2;
// 		double a = (w1-3.0)/2.0;
// 		double b = (w2 + w3 - w1 + 1.0)/2.0;
// 		double c = -w2/2.0;
// 		nroot = gsl_poly_solve_cubic(a, b, c, roots, roots+1, roots+2);
// 		if (!nroot)
// 			return 3;
// 	}
// 	double llh = 0.0, p = -1.0;
// 	if (LLR || (nroot>1)) {
// 		for (int i=0; i<nroot; i++) {
// 			double pi = roots[i];
// 			if (pi<0.0)
// 				pi = 0.0;
// 			if (pi>1.0)
// 				pi = 1.0;
// 			double Dhpi = Dh*pi;
// 			double Dhqi = Dh-Dhpi;
// 			E1[0] = ((double)T[0] + Dhpi)/Nhap;
// 			E1[1] = ((double)T[1] + Dhqi)/Nhap;
// 			E1[2] = ((double)T[2] + Dhqi)/Nhap;
// 			E1[3] = ((double)T[3] + Dhpi)/Nhap;
// 			/* Probability of double het (divided by 2)*/
// 			double EDh1 = (E1[0]*E1[3] + E1[1]*E1[2])/2.0;
// 			double llhi = 0.0;
// 			if (Dh) 
// 				llhi += Dh*log(EDh1);
// 			for (int j=0; j<4; j++) 
// 				if (T[j]) llhi += T[j]*log(E1[j]);
// 			if (p<0.0 || (llhi>llh)) {
// 				llh = llhi;
// 				p = pi;
// 			}
// 		}
// 	}
// 	else 
// 		p = roots[0];
// 	if (p<0) {
// 		return 3;
// 	}
// 
// 	/* Normal return */
// 
// 	double Dhp = Dh*p;
// 	double Dhq = Dh-Dhp;
// 	hapfreq[0] = ((double)T[0] + Dhp)/Nhap;
// 	hapfreq[1] = ((double)T[1] + Dhq)/Nhap;
// 	hapfreq[2] = ((double)T[2] + Dhq)/Nhap;
// 	hapfreq[3] = ((double)T[3] + Dhp)/Nhap;
// 	if (LLR) {
// 		double llh0=0.0;
// 		for (int i=0; i<4; i++)
// 			llh0 += margins[i]*log(margins[i]);
// 		*LLR = llh - Nhap*llh0;
// 	}
// 	return 0;
// }

/* D' from MAF */
wsReal calcDp(char **Na_data, wsUint N_sz, wsUint i, wsUint j, wsReal *Rp_r,
	 const char *Ba_isFounder)
{
	wsUint N11 = 0;
	wsUint N12 = 0;
	wsUint N21 = 0;
	wsUint N22 = 0;

	wsReal R_pAB, R_pAb, R_paB, R_pab;
	wsReal R_dAB;

	if (0) {
		if (Ba_isFounder == NULL) for (wsUint k=0 ; k<N_sz ; k++) {
			if (isAvailable(Na_data[k][i]) && isAvailable(Na_data[k][j])) {
				wsUint G = (Na_data[k][i]<<1) + Na_data[k][i] + Na_data[k][j];
				switch (G) {
				/* 0 0 = 2 0 0 0 */
				case 0: N11 += 2; break;
				/* 0 1 = 1 1 0 0 */
				case 1: N11++; N12++; break;
				/* 0 2 = 0 2 0 0 */
				case 2: N12 += 2; break;
				/* 1 0 = 1 0 1 0 */
				case 3: N11++; N21++; break;
				/* 1 1 = 1 0 0 1 or 0 1 1 0 */
				case 4: rand()%2 ? (N11++,N22++) : (N12++,N21++); break;
				/* 1 2 = 0 1 0 1 */
				case 5: N12++; N22++; break;
				/* 2 0 = 0 0 2 0 */
				case 6: N21 += 2; break;
				/* 2 1 = 0 0 1 1 */
				case 7: N21++; N22++; break;
				/* 2 2 = 0 0 0 2 */
				case 8: N22 += 2; break;
				/* Impossible here */
				default: halt("Invalid genotype [%d,%d] found", Na_data[k][i],
					Na_data[k][j]);
				}
			}
		}
		else for (wsUint k=0 ; k<N_sz ; k++) {
			if (!Ba_isFounder[k]) continue;
			if (isAvailable(Na_data[k][i]) && isAvailable(Na_data[k][j])) {
				wsUint G = (Na_data[k][i]<<1) + Na_data[k][i] + Na_data[k][j];
				switch (G) {
					/* 0 0 = 2 0 0 0 */
				case 0: N11 += 2; break;
					/* 0 1 = 1 1 0 0 */
				case 1: N11++; N12++; break;
					/* 0 2 = 0 2 0 0 */
				case 2: N12 += 2; break;
					/* 1 0 = 1 0 1 0 */
				case 3: N11++; N21++; break;
					/* 1 1 = 1 0 0 1 or 0 1 1 0 */
				case 4: rand()%2 ? (N11++,N22++) : (N12++,N21++); break;
					/* 1 2 = 0 1 0 1 */
				case 5: N12++; N22++; break;
					/* 2 0 = 0 0 2 0 */
				case 6: N21 += 2; break;
					/* 2 1 = 0 0 1 1 */
				case 7: N21++; N22++; break;
					/* 2 2 = 0 0 0 2 */
				case 8: N22 += 2; break;
					/* Impossible here */
				default: halt("Invalid genotype [%d,%d] found", Na_data[k][i],
							 Na_data[k][j]);
				}
			}
		}
		wsReal R_nA = (wsReal)N11 + (wsReal)N12;
		wsReal R_na = (wsReal)N21 + (wsReal)N22;
		wsReal R_nB = (wsReal)N11 + (wsReal)N21;
		//wsReal R_nb = (wsReal)N12 + (wsReal)N22;
		wsReal R_nAll = R_nA + R_na;

		/* pAB, pAb, paB, pab */
		R_pAB = N11/R_nAll;
		R_pAb = N12/R_nAll;
		R_paB = N21/R_nAll;
		R_pab = N22/R_nAll;

		/* pA, pB, pa, pb */
		wsReal R_pA = R_nA / R_nAll;
//		wsReal R_pa = R_na / R_nAll;
		wsReal R_pB = R_nB / R_nAll;
//		wsReal R_pb = R_nb / R_nAll;

		/* dAB */
		R_dAB = R_pAB - R_pA * R_pB;
		if (Rp_r) *Rp_r = WISARD_NAN;
	} else {
		if (i == j) return 1.0;

		char *Na_X = NULL;
		char *Na_Y = NULL;
		wsAlloc(Na_X, char, N_sz);
		wsAlloc(Na_Y, char, N_sz);

		wsUint K = 0;
		if (Ba_isFounder == NULL) for (wsUint k=0; k<N_sz ; k++) {
			if (isAvailable(Na_data[k][i]) && isAvailable(Na_data[k][j])) {
				Na_X[K] = Na_data[k][i];
				Na_Y[K] = Na_data[k][j];
				K++;
			}
		} else for (wsUint k=0 ; k<N_sz ; k++) {
			if (!Ba_isFounder[k]) continue;
			if (isAvailable(Na_data[k][i]) && isAvailable(Na_data[k][j])) {
				Na_X[K] = Na_data[k][i];
				Na_Y[K] = Na_data[k][j];
				K++;
			}
		}
		/* Now do */
		double Dprime;
		double r = (wsReal)getHfreq_RogersHuff(Na_X, Na_Y, K, &Dprime);
		if (Rp_r) *Rp_r = (wsReal)r;
		return Dprime;
	}
// 	return R_dAB < W0 ? R_dAB / min(R_pA*R_pB,R_pa*R_pb) :
// 		(R_dAB > W0 ? R_dAB / min(R_pA*R_pb,R_pa*R_pB) : W0);
 	return R_dAB < W0 ? R_dAB / min(R_pAB, R_pab) :
 		(R_dAB > W0 ? R_dAB / min(R_pAb, R_paB) : W0);
}

void cLdAnalysis::run()
{
	if (!OPT_ENABLED(ldcor) && !OPT_ENABLED(ld))
		return;
	wsUint			j;
	cTableExporter*	Cp_ld;
	if (OPT_ENABLED(ldcor))
		Cp_ld = new cTableExporter("ld.res", "ssrr",
			"Result of LD estimation", 0, 4, "SNP1", "SNP2", "LD", "LD2");
	else
		Cp_ld = new cTableExporter("ld.res", "ssrrr",
			"Result of LD estimation", 0, 5, "SNP1", "SNP2", "r", "r2", "Dprime");
	char		**Na_data	= Cp_IO->getGenotype();
	wsUint		N_samp		= Cp_IO->sizeSample();
	vVariant&	Xv_vrt		= Cp_IO->getVariant();
	wsUint		N_vrt		= Na_sset ? N_set : Cp_IO->sizeVariant();

	LOOP (i, N_vrt) {
		j = i;
		wsUint N_e = N_vrt;

		/* Fetch range if --ldbin */
		if (Na_SE.size()) N_e	= (wsUint)Na_SE[i].size();

		for ( ; j<N_e ; j++) {
			wsUint	I		= Na_sset ? Na_sset[i] : i;
			wsUint	J		= Na_sset ? Na_sset[j] : (Na_SE.size() ? Na_SE[i][j] : j);
			wsReal	R_dp	= WISARD_NAN;
			wsReal	R_val	= OPT_ENABLED(ldcor) ? calcDp(Na_data, N_samp,
				I, J, &R_dp, OPT_ENABLED(founderonly)?Cp_IO->getIsFounder():NULL) :
				calcR(Na_data, N_samp, I, J,
					OPT_ENABLED(founderonly)?Cp_IO->getIsFounder():NULL);
			wsReal	R_v2 = SQR(R_val);

			if (OPT_ENABLED(ldcor)) {
				/* --remna */
				if (OPT_ENABLED(remna) && NA(R_val)) continue;

				Cp_ld->write(4, Xv_vrt[I].name, Xv_vrt[J].name, R_val, R_v2);
			} else {
				/* --remna */
				if (OPT_ENABLED(remna) && NA(R_dp) && NA(R_val)) continue;

				Cp_ld->write(5, Xv_vrt[I].name, Xv_vrt[J].name,
					R_val, R_v2, R_dp);
			}
		}
	}

	delete Cp_ld;
}

// Time-stamp: <2006-07-26 15:47:16 zaykind> (written by Dmitri Zaykin)
//
// "Composite D-prime" code
//
// After: Zaykin DV 2004 Bounds and normalization of the composite
// linkage disequilibrium coefficient. Genet Epidemiol. 27(3):252-257.
//
// R CMD SHLIB DprKK.cpp

inline int Ix(int i, int j, int r)
{
	return i + r*j;
}

// returns sample covariance estimated as Sxy/(n-1)
double Cov(int x[], int n, int ii, int jj)
{
	// this function is modified from pseudocode in
	// en.wikipedia.org/wiki/Correlation
	double ssq_xy = 0;
	double xbar = x[Ix(0,ii,n)];
	double ybar = x[Ix(0,jj,n)];
	for(int i=1; i<n; i++) {
		double dx = x[Ix(i,ii,n)] - xbar;
		double dy = x[Ix(i,jj,n)] - ybar;
		double i1 = i+1;
		double sweep = (i1 - 1.0) / i1;
		ssq_xy += dx * dy * sweep;
		xbar += dx / i1;
		ybar += dy / i1;
	}
	double Sxy = ssq_xy / (n-1.0);
	return Sxy;
}

inline int Sum(char **q, int n, int dim)
{
	int s=0;
	for(int i=0; i<n; i++) s += q[i][n]; //q[i][dim];
	return s;
}

inline int Cnt(char** q, int n, int dim, int what)
{
	int c=0;
	for(int i=0; i<n; i++) if(q[i][n] == what) ++c; //if(q[i][dim] == what) ++c;
	return c;
}

template <class T> const T& Min (const T& a, const T& b)
{
	return a <= b ? a : b;
}

template <class T> const T& Max (const T& a, const T& b)
{
	return a > b ? a : b;
}

// q is (n by k) matrix of (-1;0;1) recoded values;
// each column for each of k SNPs; n is the number of indiv
// inline double DprIJ(char** q, int n, int i, int j)
// {
// 	double naa		= Cnt(q, n, i, 0);
// 	double nbb		= Cnt(q, n, j, 0);
// 	double nAa		= Cnt(q, n, i, 1);
// 	double nBb		= Cnt(q, n, j, 1);
// 	double nAA		= Cnt(q, n, i, 2);
// 	double nBB		= Cnt(q, n, j, 2);
// 	double pa		= 0.5 - Sum(q,n,i)/(2.0*n);
// 	double pb		= 0.5 - Sum(q,n,j)/(2.0*n);
// 	double delta	= Cov(q, n, i, j);
// 	double d,s,mxd;
// 	double ld		= (n-1.0)/n * delta/2;
// 
// 	if (delta > 0) {
// 		d	= Min(nAA,nBB) + Min(naa,nbb);
// 		s	= Max(n-d-nAa-nBb, 0.0);
// 		mxd	= (d-s)/(2.0*n) - (1-2*pa)*(1-2*pb)/2;
// 	} else {
// 		d	= Min(nAA,nbb) + Min(naa,nBB);
// 		s	= Max(n-d-nAa-nBb, 0.0);
// 		mxd	= (d-s)/(2.0*n) + (1-2*pa)*(1-2*pb)/2;
// 	}
// 	double dpr = ld/mxd;
// 	return dpr;
// }

void DprKK(char** Na_data, wsSym dp, wsUint n, wsUint k)
{
	LOOP (i, k-1) LOOP (j, i) {
		wsReal r = WISARD_NAN;
		calcDp(Na_data, n, i, j, &r, NULL);
		dp[i][j] = r;
	}
}

wsSym dprc(char** Na_data, wsUint r, wsUint c, char B_usecor=0)
{
	if (B_usecor)
		return NULL;//sseCovRow(p, r, c, 1, 1);
	wsSym dpmat = sseSymMat(r);
	sseSinit(dpmat, r, W1);
	DprKK(Na_data, dpmat, c, r);
	return dpmat;
}

wsRealCst GetStat(cIO *Cp_IO, wsUint* Na_mIdx, wsUint c, wsUint n0, wsUint n1,
	char B_usecor, wsUint iter, int* y)
{
	char**	Na_data		= Cp_IO->getGenotype();
	wsUint	n			= Cp_IO->sizeSample();
	wsUint	i0			= 0;
	wsUint	i1			= 0;
	wsUint*	Na_caIdx	= NULL;
	wsUint*	Na_ctIdx	= NULL;
	wsAlloc(Na_caIdx, wsUint, n1);
	wsAlloc(Na_ctIdx, wsUint, n0);
	int*	Na_permY	= NULL;
	if (iter) {
		wsUint N_miss	= n - n0 - n1;
		wsAlloc(Na_permY, int, n);

		LOOP (i, N_miss) {
			wsUint N_idx = 0;
			do {
				N_idx = wsRand()%n;
			} while (Na_permY[N_idx] != W0);
			Na_permY[N_idx] = WISARD_NA;
		}

		LOOP (i, n1) {
			wsUint N_idx = 0;
			do {
				N_idx = wsRand()%n;
			} while (Na_permY[N_idx] != W0);
			Na_permY[N_idx] = WISARD_NA;
		}
	}
	int*	Na_y = iter ? Na_permY : y;
	
	LOOP (i, n) {
		if (Na_y[i] == WISARD_AFFECTED)
			Na_caIdx[i1++] = Na_y[i];
		else if (Na_y[i] == WISARD_UNAFFECTED)
			Na_ctIdx[i0++] = Na_y[i];
	}

	wsUint	k	= c-1;
	wsReal	t	= WISARD_NAN;

	/* Subsetting data into ca/ct */
	char**	co = NULL;
	char**	ca = NULL;
	wsAlloc(co, char*, c);
	wsAlloc(ca, char*, c);
	LOOP (i, c) {
		wsAlloc(co[i], char, n0);
		wsAlloc(ca[i], char, n1);
		LOOP (j, n1) ca[i][j] = Na_data[Na_mIdx[i]][Na_caIdx[j]];
		LOOP (j, n0) co[i][j] = Na_data[Na_mIdx[i]][Na_ctIdx[j]];
	}
	wsSym	ca2		= dprc(ca, c, n1, B_usecor);
	wsSym	co2		= dprc(co, c, n0, B_usecor);
	sseSsS(ca2, co2, ca2, c);
	wsSym	ca22	= sseSpS(ca2, c);
	wsReal	R_denom	= (wsReal)(4*(k-1)*k);
	wsRealCst	R_numer	= diagSum(ca22, c);//sum(diag(crossprod(ca2 - co2)));
	t = R_numer / R_denom;
//	ca2[row(ca2) > col(ca2)] <- co2[row(co2) > col(co2)]
//	AvD = ca2;

//	if(iter == 0) {
//		cat("\nNumber of cases and controls: ", n1, ",", n0, "\n")
//		cat("LD matrices in cases (upper part) and controls (lower)\n")
//		print(AvD)
//		cat("Statistic value: ", t, "\n")
//	}

	return t;
}

// void LDcontrast(cIO* Cp_IO)
// {
// 
// 		t0 <- GetStat(p, c, n, n0, n1, NumMiss, method, 0)
// 
// 		if(is.nan(t0) || is.na(t0) || identical(NA,t0)) {
// 			cat("\nmonomorphic loci in the input\n")
// 				return (1)
// 		}
// 
// # Obtain statistics empirical distribution
// 		cnt1 <- 0   
// 			for(i in 1:nsim) {
// 				ti <- GetStat(p, c, n, n0, n1, NumMiss, method, i)
// 					if(i %% 10 == 0) cat(".", i, ".", sep="")
// 						if(i %% 100 == 0) cat("\r")
// 							if(is.nan(ti) || is.na(ti) || identical(NA,ti)) {
// 								i <- i-1
// 							} else {
// 								if(ti >= t0) cnt1 <- cnt1+1
// 							}
// 			}
// 
// 			dyn.unload("DprKK.so")
// 				cat("\nLD intensity difference p-value = ", x1 <- cnt1/nsim, "\n")
// 				x <- list("pv1"=x1)
// 				invisible(x)
// }

#endif

} // End namespace ONETOOL
