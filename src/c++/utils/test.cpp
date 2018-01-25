#ifdef _WIN32
#	pragma warning(disable:4068)
#endif
#include <stdarg.h>
#include "global/common.h"
#include "global/option.h"
#include "utils/test.h"
#include "utils/util.h"
#include "global/io.h"
#include "global/Rconn.h"
#include "analyses/emai.h"
#include "analyses/corr.h"
#include "utils/matrix.h"
#include "utils/data.h"
#include "output/wis.h"
#include "utils/dcdflib.h"
#include "utils/logic.h"
#include "utils/ftp.h"
#include "input/stream.h"
//#include "svg.hpp"

namespace ONETOOL {

void test_matrix();
void test_matrix_function();
void test_matrix_class();
void test_function();

int test()
{
	if (0) {
		cSymMatrix X1(2u, 1);
		cSymMatrix X2(5u, 1);
		wsReal **Ra_Y = sseSkS(X1.get(), X1.row(), X2.get(), X2.row());
		X1.file("X1");
		X2.file("X2");
		exportSymMat("Y", Ra_Y, X1.row()*X2.row());
	}
	if (0) { 
		wsUint RRRR;
		wsSym MR = makeSymMatSSE("tests/res_mmt.txt", &RRRR);
		wsSym MRinv = SinvSymMat(MR, RRRR);
		wsSym MRinvT = invSymMat(MR, RRRR);
		for (wsUint i=0 ; i<RRRR ; i++) {
			for (wsUint j=0 ; j<=i ; j++)
				printf("%g	", MRinv[i][j]);
			printf("\n");
		}
		for (wsUint i=0 ; i<RRRR ; i++) {
			for (wsUint j=0 ; j<=i ; j++)
				printf("%g	", MRinvT[i][j]);
			printf("\n");
		}
	}


	// 	wsReal R_pV = 0.3;
	// 	LOG("qchisq : %g\n", chiprobP(R_pV, 2.0));
	// 	return 1;

	// 	wsUint Rxx;
	// 	wsReal **A = makeSymMatSSE("A", &Rxx);
	// 	wsReal **Rxxx = sseSS(A, Rxx, A, Rxx);
	// 	exportMatrix("RRR", Rxxx, Rxx, Rxx);
	// 	return 1;
	if (OPT_ENABLED(testmatrix))
		test_matrix();
	if (OPT_ENABLED(testmatfunc))
		test_matrix_function();
	if (OPT_ENABLED(testmatclass))
		test_matrix_class();
	if (OPT_ENABLED(testfunc))
		test_function();
	return 0;
#if 0
	// 	exportMatrix("mysse2", Ra_ret3, N_xx1r, N_xx2);
	return 1;

	OPT_RANGE(cormaf);
#define NR 120
#define NC 120
#ifdef USE_R
	REAL_t **a = new REAL_t*[NR];
	for (int i=0,k=0 ; i<NR ; i++) {
		a[i] = new REAL_t[NC];
		for (int j=0 ; j<NC ; j++,k++)
			a[i][j] = k;
	}
	R().sendMatrix("var", a, NR, NC);
	R().Rparse("svg('haha.svg');image(var);dev.off()");
	return 3;
#endif

	// 	double ncp=0.0;
	// 	printf("%g\n", chiprobP(numeric_limits<REAL_t>::infinity(), 10.701, NULL));
	// 	return 1;
	return 0;
	wsUint A, B;
	wsReal **RR = makeMatrixSSE("aa.txt", &A, &B);
	wsReal **RRvec;
	char res;
	wsReal *R_ev = eigenDecomp(RR, A, &res, &RRvec);
	exportMatrix("ev", &R_ev, 1, A);
	exportMatrix("evec", RRvec, A, A);
	return 1;

	PVchisq(37.5525, 33.3639);
	return 0;
	cComb c(2, 10);
	wsUint N_gen;
	c.get(1000, &N_gen);
	c.get(1000, &N_gen);
	return 0;
	cTimer t;

	t.start();
	wsReal **mu = makeMatrixSSE("res.ea.scr.phi", &A, &B);
	printf("Retrieval %s\n", t.getReadable());
	wsReal **mu2 = sseMatrix(A, B, mu);

	printf("Size %d*%d\n", A, B);

	t.start();	
	wsReal **Ra_x1 = sseMpM(mu, A, B, mu2, A, B);
	printf("sseMM %s\n", t.getReadable());
	t.start();	
	wsReal **Ra_x2 = sseMpMt(mu, A, B, mu2, A, B);
	printf("sseMMt %s\n", t.getReadable());
	return 1;

	compareMatrix(Ra_x1, A, B, Ra_x2, "xx");
	return 1;
	t.start();	
	//	wsReal **mu2 = sseThrRowMM(mu, A, B, mu, A, B);
	printf("Threaded %s\n", t.getReadable());
	exportMatrix("mu2", mu2, A, B);

	t.start();	
	wsReal **mu4 = sseThrRowMMt(mu, A, B, mu, A, B);
	printf("Threaded-tt %s\n", t.getReadable());
	exportMatrix("mu4", mu2, A, B);

	t.start();	
	wsReal **mu3 = sseMpM(mu, A, B, mu, A, B);
	printf("Non-hreaded %s\n", t.getReadable());
	exportMatrix("mu3", mu3, A, B);

	compareMatrix(mu2, A, B, mu3, "xx");
	compareMatrix(mu2, A, B, mu4, "xx");

	return 1;
	double R_1 = 1.0-PVchisq(5.724552, 4.790102);
	printf("%g-%g=%g\n", R_1, 0.6908296, R_1-0.6908296);
	return 1;
	wsReal **RR3 = NULL;
	SVDinverse(RR, A, &RR3);
	wsReal **RR2 = SO_invMatrix(RR, A);
	exportMatrix("phiInv", RR2, A, B);
	return 1;
	tqli1(100, RR[0], RR[1]);
	for (wsUint i=0 ; i<100 ; i++) {
		LOG("Eigen %d : %g\n", i+1, RR[0][i]);
	}
	return 0;
	wsReal **RRR = makeMatrixSSE("tests/covRowTest_NA.txt");
	t.start();
	wsReal **covR = sseMcovRow(RRR, 1000, 3000);
	printf("CovR(sse) : %s\n", t.getReadable());
	exportMatrix("covR.txt", covR, 100, 100);
	return 1;
	cMatrix my1, my2, my3;
	wsReal *Ra_ret = new wsReal[1000];
	char N_stat;
	Ra_ret = eigen(RRR, 1000, &N_stat);
	printf("Stat : %d, %s\n", N_stat, t.getReadable());
	exportMatrix("eig.tred1", &Ra_ret, 1, 1000);
	return 1;
	wsUint N_r, N_c;
	wsReal **Ra_testMat = makeMatrix((char *)"tests/mat.txt", &N_r, &N_c);

	my2.initRandom(10, 10);
	wsReal **Ra_mv2 = my2.getData();
	exportMatrix("m2", Ra_mv2, 10, 10);

	wsReal **Ra_vcv = sseMcovCol(Ra_mv2, 10, 10);
	exportMatrix("cv", Ra_vcv, 10, 10);
	return 1;

	my1.initRandom(1, 10000);
	my2.initRandom(10000, 10000);
	my3.initRandom(10000, 1);
	wsReal **Ra_mv1 = my1.getData();
	Ra_mv2 = my2.getData();
	wsReal **Ra_mv3 = my3.getData();
	t.start();
	wsReal **r1 = multMMM(Ra_mv1, 1, 10000, Ra_mv2, 10000, 10000, Ra_mv3, 10000, 1);
	printf("multMMM() %g : %s\n", r1[0][0], t.getReadable());
	t.start();
	wsReal **r2 = sseMMM(Ra_mv1, 1, 10000, Ra_mv2, 10000, 10000, Ra_mv3, 10000, 1);
	printf("sseMMM() %g : %s\n", r2[0][0], t.getReadable());
	return 1;

	t.start();
	wsReal r = multVMV(Ra_mv1[0], Ra_mv2, 20000, Ra_mv3[0]);
	printf("\nmultVMV() function test\n%g, %s\n", r, t.getReadable());
	t.start();
	r = sseVpMpV(Ra_mv1[0], Ra_mv2, 20000, Ra_mv3[0]);
	printf("\nmultVMV_SSE() function test\n%g, %s\n", r, t.getReadable());
	return 1;

	/* Test 0 */
	// 	double qfc(double* Ra_lambda, double* Ra_nonCentrals, int* Na_df,
	// 		int N_lambda, double *sigma, double *c1, int *lim1, double *acc,
	// 		double* trace, int* ifault);
	double	R_q					= { 3 };
	double	Ra_lambda[]			= { 1, 2, 3 };

	double	Qq = davies(R_q, Ra_lambda, 3);
	LOG("Qq : %g\n", Qq);
	return 0;

	/* Test 1 */
	wsReal **Ra_testMtM = multMtM(Ra_testMat, N_r, N_c);
	compareMatrix(Ra_testMtM, N_c, N_c, (char *)"tests/res_mtm.txt",
		"multMtM() function test");

	/* Test 2 */
	wsReal **Ra_testMMt = multMMt(Ra_testMat, N_r, N_c);
	compareMatrix(Ra_testMMt, N_r, N_r, (char *)"tests/res_mmt.txt",
		"multMMt() function test");
	deallocMatrix(Ra_testMMt, N_r);

	/* Test 3 */
	wsReal R_M0MtMM0 = multVMV(Ra_testMat[0], Ra_testMtM, N_c);
	wsReal** Ra_M0MtMM0 = makeMatrix((char *)"tests/res_vmv.txt");
	printf("\nmultVMV() function test\n%g\n", R_M0MtMM0-Ra_M0MtMM0[0][0]);
	deallocMatrix(Ra_M0MtMM0, 1);

	/* Test 4 */
	wsReal* Ra_MtMM0t = multMV(Ra_testMtM, N_c, N_c, Ra_testMat[0]);
	wsReal** Ra_MtMM0 = transpose(&Ra_MtMM0t, 1, N_c);
	compareMatrix(Ra_MtMM0, N_c, 1, (char *)"tests/res_mv.txt",
		"multMV() function test");
	DEALLOC(Ra_MtMM0t);
	deallocMatrix(Ra_MtMM0, N_c);

	/* Test 5 */
	wsReal **Ra_M0 = transpose(&(Ra_testMat[0]), 1, N_c);
	wsReal **Ra_MtMM0_2 = multThrRowMM(Ra_testMtM, N_c, N_c,
		Ra_M0, N_c, 1);
	wsReal **Ra_MtMM0_3 = sseThrRowMM(Ra_testMtM, N_c, N_c,
		Ra_M0, N_c, 1);
	compareMatrix(Ra_MtMM0_2, N_c, 1, (char *)"tests/res_mv.txt",
		"mult2Matrix() function test");
	compareMatrix(Ra_MtMM0_3, N_c, 1, (char *)"tests/res_mv.txt",
		"mult2MatrixSSE() function test");
	deallocMatrix(Ra_MtMM0_2, N_c);
	sseUnmat(Ra_MtMM0_3, N_c);

	/* Test 6 */
	wsReal **Ra_M0MtMM0_2 = multMMM(&(Ra_testMat[0]), 1, N_c, Ra_testMtM,
		N_c, N_c,
		Ra_M0, N_c, 1);
	compareMatrix(Ra_M0MtMM0_2, 1, 1, (char *)"tests/res_vmv.txt",
		"mult3Matrix() function test");
	deallocMatrix(Ra_M0, N_c);
	deallocMatrix(Ra_M0MtMM0_2, 1);

	/* Test 7 */
	wsUint N_r1, N_r2, N_r3, N_c1, N_c2, N_c3;
	printf("\nReading 1st large matrix...\n");
	wsReal **Ra_l1 = makeMatrixSSE((char *)"tests/lmat1.txt", &N_r1, &N_c1);
	printf("\nReading 2nd large matrix...\n");
	wsReal **Ra_l2 = makeMatrixSSE((char *)"tests/lmat2.txt", &N_r2, &N_c2);
	printf("\nReading 3rd large matrix...\n");
	wsReal **Ra_l3 = makeMatrixSSE((char *)"tests/lmat3.txt", &N_r3, &N_c3);

	t.start();
	wsReal **Ra_r1 = multThrRowMM(Ra_l1, N_r1, N_c1, Ra_l2, N_r2, N_c2);
	printf("\n%s taken without SSE\n", t.getReadable());

	t.start();
	wsReal **Ra_r2 = sseThrRowMM(Ra_l1, N_r1, N_c1, Ra_l2, N_r2, N_c2);
	printf("\n%s taken with SSE\n", t.getReadable());

	t.start();
	wsReal **Ra_r3 = multMMM(Ra_l1, N_r1, N_c1, Ra_l2, N_r2, N_c2, Ra_l3,
		N_r3, N_c3);
	printf("\n%s taken without SSE\n", t.getReadable());

	t.start();
	wsReal **Ra_r4 = mult3MatrixSSE2(Ra_l1, N_r1, N_c1, Ra_l2, N_r2, N_c2,
		Ra_l3, N_r3, N_c3);
	printf("\n%s taken with SSE\n", t.getReadable());

	compareMatrix(Ra_r1, N_r1, N_c2, Ra_r2, "diff2.txt");
	compareMatrix(Ra_r3, N_r1, N_c3, Ra_r4, "diff3.txt");
	exportMatrix("m1", Ra_r1, N_r1, N_c2);
	exportMatrix("m2", Ra_r2, N_r1, N_c2);
	exportMatrix("m3", Ra_r3, N_r1, N_c3);
	exportMatrix("m4", Ra_r4, N_r1, N_c3);
	deallocMatrix(Ra_r1, N_r1);
	sseUnmat(Ra_r2, N_r1);
	deallocMatrix(Ra_r3, N_r1);
	sseUnmat(Ra_r4, N_r1);

	sseUnmat(Ra_l1, N_r1);
	sseUnmat(Ra_l2, N_r1);
	sseUnmat(Ra_l3, N_r1);

	deallocMatrix(Ra_testMtM, N_c);
	deallocMatrix(Ra_testMat, N_r);
	return 1;
	// 	wsReal *a;
	// 	MULTI_16_MALLOC(a, wsReal, 5);
	// 	a[0] = 1.0;
	// 	a[1] = 2.0;
	// 	a[2] = 3.0;
	// 	a[3] = 5.0;
	// 	a[4] = 7.0;
	// 	a[5] = 8.0;
	// 	a[6] = 9.0;
	// 	DEALLOC_16(a);
	// 	return 1;
	cTimer t2;

	cMatrix ma1, ma2, ma4;
	ma1.initRandom(10, 5000);
	ma2.initRandom(5000, 4000);
	ma4.initRandom(4000, 10);
	ma1.file("testMat1");
	ma2.file("testMat2");
	ma4.file("testMat3");
	t2.start();
	//	cMatrix& ma3 = ma1.product(ma2, ma4);
	wsReal X = ma1.sumDiag(ma2, ma4);
	printf("Matrix product : %s, %g\n", t2.getReadable(), X);
	//	ma3.print();
	//	ma3.file("testMatRes");
	exit(1);
	// 
	// 	//	printf("%g", qnorm(0.001, 0, 1.0/sqrt((double)(516610-3)), 0, 0));
	// 	//	return 0;
	// 	//	printf("Result %d\n", corrTest(0.2, 0.4, 1000, 0.01));
	// 	//	return 0;
	// 
	// 	wsReal _mat[3][3] = {
	// 		{ 1.0f, 2.0f, 1.0f },
	// 		{ 2.0f, 5.0f, 8.0f },
	// 		{ 4.0f, 2.0f, 4.0f },
	// 	};
	// 
	// #define NR 500
	// 	FILE *fp = fopen("mat.txt", "w+");
	// 	srand((MYuint)time(NULL));
	// 	t2.start();
	// 	wsReal **mat = new wsReal*[NR];
	// 	for (int i=0 ; i<NR ; i++) {
	// 		mat[i] = new wsReal[NR];
	// 		memset(mat[i], 0, sizeof(wsReal)*NR);
	// 		mat[i][i] += 1;
	// 		for (int j=0 ; j<NR ; j++) {
	// 			mat[i][j] = (wsReal)(rand()%100);
	// 		}
	// 		fprintf(fp, "\n");
	// 	}
	// 	fclose(fp);
	// 	printf("Array make : %s\n", t2.getReadable());

	// 	printf("Input : \n");
	// 	for (int i=0 ; i<NR ; i++) {
	// 		for (int j=0 ; j<NR ; j++)
	// 			printf("%g	", mat[i][j]);
	// 		printf("\n");
	// 	}
	// 	t2.start();
	// 	wsReal **Ra_inv2 = GJ_invMatrix(mat, NR);
	// 	printf("GJ inversion : %s\n", t2.getReadable());
	// 	t2.start();
	// 	wsReal **Ra_inv = SO_invMatrix(mat, NR);
	// 	printf("SO inversion : %s\n", t2.getReadable());
	// 
	// 	t2.start();
	// 	wsReal **Ra_mult = _mult2Matrix(mat, NR, NR, Ra_inv, NR, NR);
	// 	printf("SO Matrix multiplication : %s\n", t2.getReadable());
	// 
	// 	getIdentityError(Ra_mult, NR);
	// 
	// 	t2.start();
	// 	wsReal **Ra_mult2 = _mult2Matrix(mat, NR, NR, Ra_inv2, NR, NR);
	// 	printf("GJMatrix multiplication : %s\n", t2.getReadable());
	// 
	// 	getIdentityError(Ra_mult2, NR);
	//	t2.start();
	//	wsReal **Ra_mult = _mult3Matrix(mat, NR, NR, Ra_inv, NR, NR, Ra_mult2,
	//		NR, NR);
	//	printf("Matrix multiplication : %s\n", t2.getReadable());
	// 	printf("Output : \n");
	// 	for (int i=0 ; i<10 ; i++) {
	// 		for (int j=0 ; j<10 ; j++)
	// 			printf("%g      ", Ra_mult2[i][j]);
	// 		printf("\n");
	// 	}
	//  	printf("Output : \n");
	//  	for (int i=0 ; i<NR ; i++) {
	//  		for (int j=0 ; j<NR ; j++)
	//  			printf("%g	", Ra_inv[i][j]);
	//  		printf("\n");
	//  	}
	// 	printf("Mult output : \n");
	// 	for (int i=0 ; i<NR ; i++) {
	// 		for (int j=0 ; j<NR ; j++)
	// 			printf("%g	", Ra_mult[i][j]);
	// 		printf("\n");
	// 	}
	// 	printf("Output : \n");
	// 	for (int i=0 ; i<NR ; i++) {
	// 		for (int j=0 ; j<NR ; j++)
	// 			printf("%g	", Ra_inv[i][j]-Ra_inv2[i][j]);
	// 		printf("\n");
	// 	}

	t2.start();
	wsReal **m1 = new wsReal*[NR];
	wsReal **m2 = new wsReal*[NR];
	for (int i=0,k=1 ; i<NR ; i++) {
		sseMalloc(m1[i], wsReal, NC);
		sseMalloc(m2[i], wsReal, NC);
		for (int j=0 ; j<NC ; j++,k++) {
			m1[i][j] = (wsReal)k;
			m2[i][j] = (wsReal)k;
		}
	}
	printf("Array make : %s\n", t2.getReadable());

	t2.start();
	wsReal **res1 = krnkMM(m1, NR, NC, m2, NR, NC);
	printf("Kronecker product (np) : %s\n", t2.getReadable());

	t2.start();
	wsReal **res2 = krnkMM(m1, NR, NC, m2, NR, NC, 1);
	printf("Kronecker product (p) : %s\n", t2.getReadable());
	for (int i=0 ; i<(NR*NR) ; i++) {
		for (int j=0 ; j<(NC*NC) ; j++) {
			if (res1[i][j] == res2[i][j])
				continue;

			printf("Result is differ between res1 and res2 in %d[%g],%d[%g]\n",
				i, res1[i][j], j, res2[i][j]);
			exit(1);
		}
	}
	for (wsUint i=0 ; i<NR*NR ; i++)
		sseFree(res1[i]);
	DEALLOC(res1);
	exit(3);
#endif

	// 	cRconn r;
	// 	r.init();

	//	launchR("test.r");
	halt("Unreachble statement reached");
	return 0; /* Do not terminate */
	//return 1; /* Do terminate */
}

void test_matrix()
{
	cTimer xxx;
	
	wsUint N_sz[5] = {1000, 2000, 3000, 4000, 5000 };
	char *Sp_out = strdup(OPT_STRING(out));
	char *S_newOut = NULL;
	wsAlloc(S_newOut, char, 512);

	LOG("test_matrix() start\n");

	/* Sparse matrix test */ {
		wsUint sz = 1000;
		wsReal **Ra_m1 = sseMrand(sz, sz);
		cStdMatrix M(sz, sz, Ra_m1);
		cSpsMatrix S(M);
		cStdMatrix M2(S);
		cTimer q;
		q.start();
		M * S;
		LOG("Sparse(full) %s\n", q.getReadable());
		q.start();
		M * M2;
		LOG("Non-sparse(full) %s\n", q.getReadable());
		LOG(" Diff full sparse : %g\n",
			compMMfrobenius(M.get(), M.row(), M.col(), M2.get(), M2.row(), M2.col()));

		wsMat	Ra_m2 = sseMrand01(sz, sz, 0.05);
		cStdMatrix M3(sz, sz, Ra_m2);
		cSpsMatrix S2(M3);
		cStdMatrix M4(S2);
		q.start();
		M3 * S2;
		LOG("Sparse(spar) %s\n", q.getReadable());
		q.start();
		M3 * M4;
		LOG("Non-sparse(spar %s\n", q.getReadable());
		LOG(" Diff sparse : %g\n",
			compMMfrobenius(M3.get(), M3.row(), M3.col(), M4.get(), M4.row(), M4.col()));
	}

	if (!IS_ASSIGNED(nperm))
		OPT_NUMBER(nperm) = 1;

	wsUint P = OPT_NUMBER(nperm);

	/* For all test case */
	for (wsUint z=0 ; z<5 ; z++) {
		/* Set new output prefix */
		wsUint sz = N_sz[z];
		sprintf(S_newOut, "%s_%d", Sp_out, sz);
		OPTION().assign("out", S_newOut, 1);
		OPTION().FORCE_OPT_STRING(out);

		/* Make both of sym/non-sym matrix */
		LOG(" --- SIZE %d ---\n", sz);

		if (P == 1) {
			wsReal **Ra_ret = NULL;
			wsReal **Ra_ret2 = NULL;
			wsReal **Ra_m1 = sseMrand(sz, sz);
			wsReal **Ra_m2 = sseMrand(sz, sz);
			xxx.start();
			wsReal **Ra_m2t = transpose(Ra_m2, sz, sz);
			LOG("	time(M')           %s\n", xxx.getReadable());
			wsSym Ra_s1 = sseSrand(sz);
			wsSym Ra_s2 = sseSrand(sz);
			if (OPT_ENABLED(verbose)) {
				exportMatrix("m1", Ra_m1, sz, sz);
				exportMatrix("m2", Ra_m2, sz, sz);
			}

//			xxx.start();
//			Ra_ret = multMM(Ra_m1, sz, sz, Ra_m2, sz, sz);
//			LOG("	dMM	%s\n", xxx.getReadable());
//			deallocMatrix(Ra_ret, sz);

			xxx.start();
			Ra_ret = sseMpM(Ra_m1, sz, sz, Ra_m2, sz, sz);
			LOG("	time(M1��M2)       %s\n", xxx.getReadable());
			if (OPT_ENABLED(verbose))
				exportMatrix("MM", Ra_ret, sz, sz);
			sseUnmat(Ra_ret, sz);

			xxx.start();
			Ra_ret = sseMtS(Ra_m2, sz, sz, Ra_s1, sz);
			LOG("	time(M'��S)        %s\n", xxx.getReadable());
			xxx.start();
			wsMat Ra_ttt = transpose(Ra_m2, sz, sz);
			Ra_ret2 = sseMpS(Ra_ttt, sz, sz, Ra_s1, sz);
			sseUnmat(Ra_ttt, sz);
			LOG("	trans-MS %s\n", xxx.getReadable());
			LOG("	MtS-MS diff %g\n", compMMfrobenius(Ra_ret, sz, sz, Ra_ret2, sz, sz));
			sseUnmat(Ra_ret, sz);
			sseUnmat(Ra_ret2, sz);

			xxx.start();
			Ra_ret = sseMSMt(Ra_m1, sz, sz, Ra_m2, sz);
			wsReal R_retT = W0;
			for (wsUint i=0 ; i<sz ; i++)
				R_retT += Ra_ret[i][i];
			LOG("	trMSMt(indir)      %s\n", xxx.getReadable());
			sseUnmat(Ra_ret, sz);
			xxx.start();
			{
				cStdMatrix M(sz, sz, Ra_m1, MATDEL_NONE);
				cSymMatrix S(Ra_m2, sz, 1);
				wsReal R_retC = M.trMSMt(S);
				LOG("	trMSMt(dir)	   %s, diff %g\n", xxx.getReadable(),
					R_retC - R_retT);
				
				/* MSM type, to call MSMt */

			}

			xxx.start();
			wsSym Ra_des = sseMpMt(Ra_m2, sz, sz);
			LOG("	time(A��A')        %s\n", xxx.getReadable());
			xxx.start();
			wsSym Ra_sym = sym_sseMpM(Ra_m2, sz, sz, Ra_m2t, sz, sz);
			LOG("	time(sym(A��B'))   %s\n", xxx.getReadable());
			LOG("	diff(AA', AB')   = %g\n", compSSfrobenius(Ra_des, sz,
				Ra_sym, sz));
			sseUnmat(Ra_des, sz);
			sseUnmat(Ra_sym, sz);

			Ra_ret = sseMpMt(Ra_m1, sz, sz, Ra_m1, sz, sz);
			if (OPT_ENABLED(verbose))
				exportMatrix("evMMt", Ra_ret, sz, sz);

			//char b = 0;
//			MAT_t M_ev;
//			xxx.start();
//			wsReal *Ra_v1 = eigenDecomp(Ra_ret, sz, &b, &M_ev);
//			LOG("	e1	%s\n", xxx.getReadable());
//			exportMatrix("ev1", M_ev, sz, sz);
//			exportMatrix("ex1", &Ra_v1, 1, sz);
//			deallocMatrix(M_ev, sz, (void *)1);
// 			xxx.start();
// 			wsReal *Ra_v1 = eigenDecomp(Ra_ret, sz, &b, &M_ev);
// 			LOG("	e1	%s\n", xxx.getReadable());
// 			if (OPT_ENABLED(verbose)) {
// 				exportMatrix("ev1", M_ev, sz, sz);
// 				exportMatrix("ex1", &Ra_v1, 1, sz);
// 			}
// 			sseUnmat(M_ev, sz);
// 			xxx.start();
// 			wsReal *Ra_v2 = EIGENDECOMPOSITION(Ra_ret, sz, &M_ev, 1);
// 			LOG("	e2	%s\n", xxx.getReadable());
// 			if (OPT_ENABLED(verbose)) {
// 				exportMatrix("ev2", M_ev, sz, sz);
// 				exportMatrix("ex2", &Ra_v2, 1, sz);
// 			}

			xxx.start();
			cSymMatrix VVV(Ra_ret, sz, 1);
			cSymMatrix &RV = VVV.inv();
			RV.rem();
			LOG("	time(S^-1)         %s\n", xxx.getReadable());
			sseUnmat(Ra_ret, sz);

			xxx.start();
			Ra_ret = sseMpMt(Ra_m1, sz, sz, Ra_m2, sz, sz);
			LOG("	time(A��B')        %s\n", xxx.getReadable());
			if (OPT_ENABLED(verbose))
				exportMatrix("MMt", Ra_ret, sz, sz);

			Ra_ret2 = sseMpM(Ra_m1, sz, sz, Ra_m2t, sz, sz);
			/* Get difference between sseMM and sseMMt */
			LOG("	diff(MM,MMt)     = %g\n", compMMfrobenius(Ra_ret, sz,
				sz, Ra_ret2, sz, sz));
			sseUnmat(Ra_ret2, sz);
			sseUnmat(Ra_ret, sz);

			xxx.start();
			Ra_ret = sseSpM(Ra_s1, sz, Ra_m2, sz, sz);
			LOG("	time(S��M)         %s\n", xxx.getReadable());
			if (OPT_ENABLED(verbose))
				exportMatrix("SM", Ra_ret, sz, sz);

			/* Get difference between sseSM and sseSpMt */
			wsReal **Ra_sExp = sseSym2Mat(Ra_s1, sz);
			Ra_ret2 = sseMpM(Ra_sExp, sz, sz, Ra_m2, sz, sz);
			if (OPT_ENABLED(verbose))
				exportMatrix("SMtrue", Ra_ret2, sz, sz);
			LOG("	diff(MM,SM)      = %g\n", compMMfrobenius(Ra_ret, sz, sz,
				Ra_ret2, sz, sz));
			sseUnmat(Ra_ret2, sz);
			sseUnmat(Ra_sExp, sz);
			sseUnmat(Ra_ret, sz);

			xxx.start();
			Ra_ret = sseSpMt(Ra_s1, sz, Ra_m2, sz, sz);
			LOG("	time(S��M')        %s\n", xxx.getReadable());
			if (OPT_ENABLED(verbose))
				exportMatrix("SMt", Ra_ret, sz, sz);

			Ra_ret2 = sseSpM(Ra_s1, sz, Ra_m2t, sz, sz);
			/* Get difference between sseSM and sseSpMt */
			if (OPT_ENABLED(verbose))
				exportMatrix("SMttrue", Ra_ret2, sz, sz);
			LOG("	diff(SM,SMt)     = %g\n", compMMfrobenius(Ra_ret, sz, sz,
				Ra_ret2, sz, sz));
			sseUnmat(Ra_ret2, sz);
			sseUnmat(Ra_ret, sz);

			xxx.start();
			Ra_ret = sseMpS(Ra_m1, sz, sz, Ra_s2, sz);
			LOG("	time(M��S)         %s\n", xxx.getReadable());
			if (OPT_ENABLED(verbose))
				exportMatrix("MS", Ra_ret, sz, sz);

			Ra_sExp = sseSym2Mat(Ra_s2, sz);
			Ra_ret2 = sseMpMt(Ra_m1, sz, sz, Ra_sExp, sz, sz);
			/* Get difference between sseSM and sseSpMt */
			if (OPT_ENABLED(verbose))
				exportMatrix("MStrue", Ra_ret2, sz, sz);
			LOG("	diff(MS,MMt)     = %g\n", compMMfrobenius(Ra_ret, sz, sz,
				Ra_ret2, sz, sz));
			sseUnmat(Ra_ret2, sz);
			sseUnmat(Ra_ret, sz);

			xxx.start();
			Ra_ret = sseSpS(Ra_s1, sz, Ra_s2, sz);
			LOG("	time(S1��S2)       %s\n", xxx.getReadable());
			if (OPT_ENABLED(verbose))
				exportMatrix("SS", Ra_ret, sz, sz);
			//sseUnmat(Ra_ret, sz);

			wsReal **Ra_sExp2 = sseSym2Mat(Ra_s1, sz);
			Ra_ret2 = sseMpMt(Ra_sExp2, sz, sz, Ra_sExp, sz, sz);
			if (OPT_ENABLED(verbose))
				exportMatrix("SStrue", Ra_ret2, sz, sz);
			LOG("	diff(SS,MMt)     = %g\n", compMMfrobenius(Ra_ret, sz, sz,
				Ra_ret2, sz, sz));
			if (OPT_ENABLED(verbose)) {
				exportMatrix("s1", Ra_sExp2, sz, sz);
				exportMatrix("s2", Ra_sExp, sz, sz);
			}
			sseUnmat(Ra_sExp2, sz);
			sseUnmat(Ra_ret2, sz);
			sseUnmat(Ra_ret, sz);
			sseUnmat(Ra_sExp, sz);

			xxx.start();
			Ra_ret = sseSpS(Ra_s1, sz);
			LOG("	time(S��S')        %s\n", xxx.getReadable());
			Ra_ret2 = sseSpS(Ra_s1, sz, Ra_s1, sz);
			if (OPT_ENABLED(verbose))
				exportMatrix("SS11true", Ra_ret2, sz, sz);
			LOG("	diff(selfSS,SS)  = %g\n", compSSfrobenius(Ra_ret, sz,
				Ra_ret2, sz));
			if (OPT_ENABLED(verbose))
				exportSymMat("SS11", Ra_ret, sz);
			sseUnmat(Ra_ret2, sz);
			sseUnmat(Ra_ret, sz);

			sseUnmat(Ra_m1, sz);
			sseUnmat(Ra_m2, sz);
			sseUnmat(Ra_s1, sz);
			sseUnmat(Ra_s2, sz);
		} else {
			wsReal R_mm = W0;
			wsReal R_mmt = W0;
			wsReal R_sm = W0;
			wsReal R_smt = W0;
			wsReal R_ms = W0;
			wsReal R_ss = W0;

			for (wsUint p=0 ; p<P ; p++) {
				wsReal **Ra_ret = NULL;
				wsReal **Ra_m1 = sseMrand(sz, sz);
				wsReal **Ra_m2 = sseMrand(sz, sz);
				wsSym Ra_s1 = sseSrand(sz);
				wsSym Ra_s2 = sseSrand(sz);

				xxx.start();
				Ra_ret = sseMpM(Ra_m1, sz, sz, Ra_m2, sz, sz);
				R_mm += xxx.get();
				sseUnmat(Ra_ret, sz);

				xxx.start();
				Ra_ret = sseMpMt(Ra_m1, sz, sz, Ra_m2, sz, sz);
				R_mmt += xxx.get();
				sseUnmat(Ra_ret, sz);

				xxx.start();
				Ra_ret = sseSpM(Ra_s1, sz, Ra_m2, sz, sz);
				R_sm += xxx.get();
				sseUnmat(Ra_ret, sz);

				xxx.start();
				Ra_ret = sseSpMt(Ra_s1, sz, Ra_m2, sz, sz);
				R_smt += xxx.get();
				sseUnmat(Ra_ret, sz);

				xxx.start();
				Ra_ret = sseMpS(Ra_m1, sz, sz, Ra_s2, sz);
				R_ms += xxx.get();
				sseUnmat(Ra_ret, sz);

				xxx.start();
				Ra_ret = sseSpS(Ra_s1, sz, Ra_s2, sz);
				R_ss += xxx.get();
				sseUnmat(Ra_ret, sz);

				sseUnmat(Ra_m1, sz);
				sseUnmat(Ra_m2, sz);
				sseUnmat(Ra_s1, sz);
				sseUnmat(Ra_s2, sz);
			}
			LOG("	MM	" REAL_FMT " ms\n", R_mm/(wsReal)P);
			LOG("	MMt	" REAL_FMT " ms\n", R_mmt/(wsReal)P);
			LOG("	SM	" REAL_FMT " ms\n", R_sm/(wsReal)P);
			LOG("	SMt	" REAL_FMT " ms\n", R_smt/(wsReal)P);
			LOG("	MS	" REAL_FMT " ms\n", R_ms/(wsReal)P);
			LOG("	SS	" REAL_FMT " ms\n", R_ss/(wsReal)P);
		}
	}

	DEALLOC(S_newOut);
	LOG("test_matrix() end\n");
}

void test_matrix_function()
{
	double r = pT(2.0, 1.0);
	LOG("t-distribution x %g df %g : %g\n", 2.0, 1.0, r);
#ifndef USE_R
	halt("--testmatfunc cannot performed without R functionality");
#else
#endif
}

void test_matrix_class()
{
#ifndef USE_R
	halt("--testmatclass cannot performed without R functionality");
#else
	wsUint N_rsz[8] = { 1, 1, 1, 1, 3, 1, 2, 11 };
	wsUint N_sz[8] = { 1, 2, 3, 5, 9, 10, 20, 55 };
	for (wsUint i=0 ; i<8 ; i++) {
		char xx[1024];
		wsUint n = N_sz[i];
		wsUint b = N_rsz[i];
		wsUint x = n/b;

		cStdMatrix	M(n, n, NULL, MATDEL_NONE);
		cSymMatrix	S(n, 1);
		cDiagMatrix	D(n, 1);
		cBlkMatrix	B(n/b, b, 1, 1);
		cVector		Vf(n, 1);
		cVector		Vb(b, 1);
		cVector		Vf2(n, 1);
		cVector		Vb2(b, 1);
		cStdMatrix	Mb(b, b + 1, NULL, MATDEL_NONE);
		cStdMatrix	MbEq(n, b + 1, NULL, MATDEL_NONE);
		cStdMatrix	E;

		sprintf(xx, "m%d", n);
		M.file(xx, 20);
		sprintf(xx, "s%d", n);
		S.file(xx, 20);
		sprintf(xx, "b%d", n);
		B.file(xx, 20);
		sprintf(xx, "d%d", n);
		D.file(xx, 20);
		sprintf(xx, "d%d", n);
		D.file(xx, 20);
		sprintf(xx, "mb%d", n);
		Mb.file(xx, 20);
		Mb.setR("Mb");
		sprintf(xx, "mbe%d", n);
		MbEq.file(xx, 20);
		MbEq.setR("Mbe");

		R().sendVector("V", Vf.get(), Vf.size());
		R().sendVector("Vb", Vb.get(), Vb.size());
		R().sendVector("V2", Vf2.get(), Vf2.size());
		R().sendVector("Vb2", Vb2.get(), Vb2.size());

		sprintf(xx, "Vb <- rep(Vb, each=%d)", x);
		R().Rparse(xx);
		sprintf(xx, "Vb2 <- rep(Vb2, each=%d)", x);
		R().Rparse(xx);
		sprintf(xx, "write.table(Vb, \"%s.vb%d\", row.names=F, col.names=F, quote=F)",
			OPT_STRING(out), n);
		R().Rparse(xx);
		sprintf(xx, "write.table(Vb2, \"%s.vb2%d\", row.names=F, col.names=F, quote=F)",
			OPT_STRING(out), n);
		R().Rparse(xx);
		sprintf(xx, "write.table(V, \"%s.v%d\", row.names=F, col.names=F, quote=F)",
			OPT_STRING(out), n);
		R().Rparse(xx);
		sprintf(xx, "write.table(V2, \"%s.v2%d\", row.names=F, col.names=F, quote=F)",
			OPT_STRING(out), n);
		R().Rparse(xx);

		/* MM */
		cStdMatrix M2 = M*M;
		M.setR("M");
		R().Rparse("ret <- M%*%M");
		E.getR("ret");
		LOG("MM inaccuracy : %g\n", M2.normF(E));
		sprintf(xx, "m%d.MM", n);
		M2.file(xx, 20);
		M2.rem();

		/* MMt */
		cStdMatrix M3 = M.Mt(M);
		R().Rparse("ret <- M%*%t(M)");
		E.getR("ret");
		LOG("MMt inaccuracy : %g\n", M3.normF(E));
		sprintf(xx, "m%d.MMt", n);
		M3.file(xx, 20);
		M3.rem();

		/* MS */
		cStdMatrix M4 = M*S;
		S.setR("S");
		R().Rparse("ret <- M%*%S");
		E.getR("ret");
		LOG("MS inaccuracy : %g\n", M4.normF(E));
		sprintf(xx, "m%d.MS", n);
		M4.file(xx, 20);
		M4.rem();

		/* MD */
		cStdMatrix M5 = M*D;
		D.setR("D");
		R().Rparse("ret <- M%*%D");
		E.getR("ret");
		LOG("MD inaccuracy : %g\n", M5.normF(E));
		sprintf(xx, "m%d.MD", n);
		M5.file(xx, 20);
		M5.rem();

		/* MB */
		cStdMatrix M6 = M*B;
		B.setR("B");
		R().Rparse("ret <- M%*%B");
		E.getR("ret");
		LOG("MB inaccuracy : %g\n", M6.normF(E));
		sprintf(xx, "m%d.MB", n);
		M6.file(xx, 20);
		M6.rem();

		/* SM */
		cStdMatrix S2 = S*M;
		R().Rparse("ret <- S%*%M");
		E.getR("ret");
		LOG("SM inaccuracy : %g\n", S2.normF(E));
		sprintf(xx, "m%d.SM", n);
		S2.file(xx, 20);
		S2.rem();

		/* SS */
		cStdMatrix S3 = S*S;
		R().Rparse("ret <- S%*%S");
		E.getR("ret");
		LOG("SS inaccuracy : %g\n", S3.normF(E));
		sprintf(xx, "m%d.SS", n);
		S3.file(xx, 20);
		S3.rem();

		/* SD */
		cStdMatrix S4 = S*D;
		R().Rparse("ret <- S%*%D");
		E.getR("ret");
		LOG("SD inaccuracy : %g\n", S4.normF(E));
		sprintf(xx, "m%d.SD", n);
		S4.file(xx, 20);
		S4.rem();

		/* SB */
		cStdMatrix S5 = S*B;
		R().Rparse("ret <- S%*%B");
		E.getR("ret");
		LOG("SB inaccuracy : %g\n", S5.normF(E));
		sprintf(xx, "m%d.SB", n);
		S5.file(xx, 20);
		S5.rem();

		/* BM (rep) */
		cStdMatrix &B1 = B*Mb;
		R().Rparse("ret <- S%*%Mb[rep(1:dim(Mb)[1], each=%d),]", b);
		E.getR("ret");
		LOG("BM(rep) inaccuracy : %g\n", B1.normF(E));
		sprintf(xx, "m%d.BMrep", n);
		B1.file(xx, 20);
		B1.rem();

		/* BM (full) */
		cStdMatrix &B2 = B*Mb;
		R().Rparse("ret <- S%*%Mbe", b);
		E.getR("ret");
		LOG("BM(full) inaccuracy : %g\n", B2.normF(E));
		sprintf(xx, "m%d.BMrep", n);
		B2.file(xx, 20);
		B2.rem();

		/* qf */
		wsReal R_qf = Vf.qf(M);
		R().Rparse("ret <- t(as.matrix(V))%*%M%*%as.matrix(V)");
		E.getR("ret");
		LOG("qf aMa(1n nn(M) n1) : %g (my) and %g (R)\n", R_qf, E.get()[0][0]);
		/* qf */
		R_qf = Vf.qf(M, Vf2);
		R().Rparse("ret <- t(as.matrix(V))%*%M%*%as.matrix(V2)");
		E.getR("ret");
		LOG("qf aMx (1n nn(M) n1) : %g (my) and %g (R)\n", R_qf, E.get()[0][0]);
		/* qf II */
		R_qf = Vf.qf(B);
		R().Rparse("ret <- t(as.matrix(V))%*%B%*%as.matrix(V)");
		E.getR("ret");
		LOG("qf aBa (1n nn(B) n1) : %g (my) and %g (R)\n", R_qf, E.get()[0][0]);
		/* qf II */
		R_qf = Vf.qf(B, Vf2);
		R().Rparse("ret <- t(as.matrix(V))%*%B%*%as.matrix(V2)");
		E.getR("ret");
		LOG("qf aBx (1n nn(B) n1) : %g (my) and %g (R)\n", R_qf, E.get()[0][0]);
		/* qf III */
		R_qf = Vb.qf(B);
		R().Rparse("ret <- t(as.matrix(Vb))%*%B%*%as.matrix(Vb)");
		E.getR("ret");
		LOG("qf bBb (1k nn(B) k1) : %g (my) and %g (R)\n", R_qf, E.get()[0][0]);
		/* qf IV */
		R_qf = Vb.qf(B, Vb2);
		R().Rparse("ret <- t(as.matrix(Vb))%*%B%*%as.matrix(Vb2)");
		E.getR("ret");
		LOG("qf bBy (1k nn(B) k1) : %g (my) and %g (R)\n", R_qf, E.get()[0][0]);

		/* BS */
		cMatrix &B3 = B-S;
		R().Rparse("ret <- B-S");
		E.getR("ret");
		LOG("B-S inaccuracy : %g\n", ((cSymMatrix *)&B3)->normF(E));
		sprintf(xx, "m%d.BmS", n);
		B3.file(xx, 20);
		B3.rem();		
	}
#endif
}

/*

# 3*3 blk matrix
b <- matrix(rnorm(9), nrow=3)

# 2*3 mat
v3 <- matrix(rnorm(6), nrow=2)
# 2*6 mat
v6 <- matrix(rnorm(12), nrow=2)

xb <- function(x, b, m) {
	# Get K (block size of b)
	K <- dim(b)[1]
	# Get dim of x [p*m] or [p*Km]
	cx <- dim(x)[2]

	R <- matrix(0, nrow=dim(x)[1], ncol=dim(b)[1]*m)
	if (cx == K*m) { # [p*Km] %*% [Km*Km]
		for (i in 1:m) {
			s <- (i-1)*K + 1
			e <- i*K
			cat(s, "~", e, "\n")
			R[,s:e] <- x[,s:e] %*% b
		}
	} else if (cx == m) { # [p*m] %*% [Km*Km]
		cb <- colSums(b)
		for (i in 1:cx) {
			s <- (i-1)*K + 1
			e <- i*K
			# x[1,i]*col(b)[1] ... x[1,i]*col(b)[k]
			#   ...
			# x[p,i]*col(b)[1] ... x[p,i]*col(b)[k]
			R[,s:e] <- x[,i,drop=F] %x% t(cb)
		}
	} else stop("Non-compatible")
	return(R)
}

xb(v3, b, 1)
xb(v6, b, 2)

xb(v3, b, 3)

m <- function(fn) { as.matrix(read.table(fn)) }
F <- function(m,n) { sum((m-n)^2) }

mmt <- list()
mm <- ms <- md <- mb <- list()
sm <- ss <- sd <- sb <- list()

W.mmt <- list()
W.mm <- W.ms <- W.md <- W.mb <- list()
W.sm <- W.ss <- W.sd <- W.sb <- list()

j <- 1
M <- list()
B <- list()
D <- list()
S <- list()
V <- list()
V2 <- list()
Vb <- list()
Vb2 <- list()
for (i in c(1, 2, 3, 5, 9, 10, 20, 55)) {
	M[[j]] <- m(paste("res.m", i, sep=""))
	B[[j]] <- m(paste("res.b", i, sep=""))
	D[[j]] <- m(paste("res.d", i, sep=""))
	S[[j]] <- m(paste("res.s", i, sep=""))
	V[[j]] <- m(paste("res.v", i, sep=""))
	V2[[j]] <- m(paste("res.v2", i, sep=""))
	Vb[[j]] <- m(paste("res.vb", i, sep=""))
	Vb2[[j]] <- m(paste("res.vb2", i, sep=""))
d	mm[[j]] <- M[[j]]%*%M[[j]]
	mmt[[j]] <- M[[j]]%*%t(M[[j]])
	ms[[j]] <- M[[j]]%*%S[[j]]
	md[[j]] <- M[[j]]%*%D[[j]]
	mb[[j]] <- M[[j]]%*%B[[j]]

	ss[[j]] <- S[[j]]%*%S[[j]]
	sm[[j]] <- S[[j]]%*%M[[j]]
	sd[[j]] <- S[[j]]%*%D[[j]]
	sb[[j]] <- S[[j]]%*%B[[j]]

	W.mm[[j]] <- m(paste("res.m", i, ".MM", sep=""))
	W.mmt[[j]] <- m(paste("res.m", i, ".MMt", sep=""))
	W.ms[[j]] <- m(paste("res.m", i, ".MS", sep=""))
	W.md[[j]] <- m(paste("res.m", i, ".MD", sep=""))
	W.mb[[j]] <- m(paste("res.m", i, ".MB", sep=""))
	
	W.sm[[j]] <- m(paste("res.m", i, ".SM", sep=""))
	W.ss[[j]] <- m(paste("res.m", i, ".SS", sep=""))
	W.sd[[j]] <- m(paste("res.m", i, ".SD", sep=""))
	W.sb[[j]] <- m(paste("res.m", i, ".SB", sep=""))
	
	cat("Fmm[",i,"] <- ", F(mm[[j]], W.mm[[j]]), "\n")
	cat("Fmmt[",i,"] <- ", F(mmt[[j]], W.mmt[[j]]), "\n")
	cat("Fms[",i,"] <- ", F(ms[[j]], W.ms[[j]]), "\n")
	cat("Fmd[",i,"] <- ", F(md[[j]], W.md[[j]]), "\n")
	cat("Fmb[",i,"] <- ", F(mb[[j]], W.mb[[j]]), "\n")

	cat("Fsm[",i,"] <- ", F(sm[[j]], W.sm[[j]]), "\n")
	cat("Fss[",i,"] <- ", F(ss[[j]], W.ss[[j]]), "\n")
	cat("Fsd[",i,"] <- ", F(sd[[j]], W.sd[[j]]), "\n")
	cat("Fsb[",i,"] <- ", F(sb[[j]], W.sb[[j]]), "\n")

	j <- j+1
}

*/

#ifdef USE_SSE4

/* _mm_cvtepu8_epi16	SSE4
 * _mm_unpackhi_epi8	SSE2
 * _mm_mullo_epi16		SSE2
 * _mm_shuffle_epi8		SSSE3
 *
 * Requires SSE4
 */
inline __m128i sseGenoMul(__m128i aNew, __m128i bNew) {
#define MASK (char)0x80
		/* Separate LO part */
		__m128i Z = _mm_setzero_si128();
		__m128i Rlo = _mm_mullo_epi16(_mm_cvtepu8_epi16(aNew), _mm_cvtepu8_epi16(bNew));
		/* Separate HI part */
		__m128i Rhi = _mm_mullo_epi16(_mm_unpackhi_epi8(aNew,Z), _mm_unpackhi_epi8(bNew,Z));
		/* Combin & ret */
		__m128i maskLo = _mm_set_epi8(MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK, 14, 12, 10, 8, 6, 4, 2, 0);
		__m128i maskHi = _mm_set_epi8(14, 12, 10, 8, 6, 4, 2, 0, MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK);
		return _mm_or_si128(_mm_shuffle_epi8(Rlo, maskLo), _mm_shuffle_epi8(Rhi, maskHi));
#undef MASK
}

union char32
{
	WISARD_pc i[2];
	unsigned char v[16];
	__m128i s;
};

inline __m128i sseGenoMulCompl(__m128i a, __m128i b) {
		return sseGenoMul(a, b);
}

typedef union _w128 {
		__m128i s;
		unsigned int i[4];
} w128;

/* Reduce sixteen 8-bit integers into eight 16-bit integers by their adjacent sum */
#define _mm_rdc_8to16(s,Z)	_mm_add_epi16(_mm_cvtepu8_epi16(s), _mm_unpackhi_epi8((s), (Z)))
/* Reduce eight 16-bit integers into four 32-bit integers by their adjacent sum */
#define _mm_rdc_16to32(s,Z)	_mm_add_epi32(_mm_cvtepu16_epi32(s), _mm_unpackhi_epi16((s), (Z)))
#define _mm_part_ing_sum(part,zero,dest) { __m128i s16 = _mm_rdc_8to16(part, zero); \
	dest = _mm_add_epi32(dest, _mm_rdc_16to32(s16, zero)); \
	part = _mm_setzero_si128(); }
#define _mm_part_last_sum(part,zero,accu,dest) { __m128i s16 = _mm_rdc_8to16(part, zero); \
	dest = _mm_add_epi32(accu, _mm_rdc_16to32(s16, zero)); }

/* Equivalent code
wsUint R = 0;
for (wsUint i=0 ; i<N_col ; i++)
	R += Na_1[i] * Na_2[i];
return R;
 * Required compatibility : SSE4
 */
inline wsUint sseVVchar(char *Na_1, wsUint N_col, char *Na_2)
{
	wsUint	N_med = 0;
	wsUint	N_ret = 0;
#ifdef USE_SSE
	w128	W_sse;
	memset(&W_sse, 0x00, sizeof(w128));
	N_med = get16(N_col);
	__m128i sse_part	= _mm_setzero_si128();
	__m128i sse_XY		= _mm_setzero_si128();
	for (wsUint i=0 ; i<N_med ; i+=16) {
		__m128i *a = (__m128i *)(Na_1 + i);
		__m128i *b = (__m128i *)(Na_2 + i);
		sse_part = _mm_add_epi8(sse_part, sseGenoMulCompl(*a, *b)); //SSE2

		/* Avoid an overflow of partial sum */
		if ((i%512) == (512-16)) {
			/* Lower part */
			__m128i Z	= _mm_setzero_si128();
			/* As eight 16-bit */
			__m128i s16	= _mm_rdc_8to16(sse_part, Z);
			sse_XY		= _mm_add_epi32(sse_XY, _mm_rdc_16to32(s16, Z));
			sse_part	= _mm_setzero_si128();
		}
	}
	/* Lower part */
	__m128i Z	= _mm_setzero_si128();
	/* As eight 16-bit */
	__m128i s16	= _mm_rdc_8to16(sse_part, Z);
	W_sse.s		= _mm_add_epi32(sse_XY, _mm_rdc_16to32(s16, Z));
	N_ret += W_sse.i[0] + W_sse.i[1] + W_sse.i[2] + W_sse.i[3];
#endif
	/* Summation using CPU */
	for (wsUint i=N_med ; i<N_col ; i++)
		N_ret += Na_1[i] * Na_2[i];

	return N_ret;
}

/* Equivalent code
wsUint R = 0;
for (wsUint i=0 ; i<N_col ; i++)
	R += Na_1[i] * Na_2[i];
return R;
 * Required compatibility : SSSE3
 */
inline wsReal sseVVchar2(char *Na_1, wsUint N_col, char *Na_2)
{
	wsUint	N_med = 0;
	wsUint	N_ret = 0;
#ifdef USE_SSE
	w128	W_sse;
	memset(&W_sse, 0x00, sizeof(w128));
	N_med = get16(N_col);
	__m128i sse_XY		= _mm_setzero_si128();
	for (wsUint i=0,j=0 ; i<N_med ; i+=16) {
		__m128i *a = (__m128i *)(Na_1 + i);
		__m128i *b = (__m128i *)(Na_2 + i);
		sse_XY = _mm_add_epi16(sse_XY, _mm_maddubs_epi16(*a, *b)); //SSSE3
		if (++j == 256) {
			unsigned short Na_sseSumBuf[32], *Rp_sseStorePtr;
		        Rp_sseStorePtr = (unsigned short *)alignSSE(Na_sseSumBuf);
		        _mm_store_si128((__m128i *)Rp_sseStorePtr, sse_XY);
		        for (wsUint X=0 ; X<8 ; X++)
                		N_ret += Rp_sseStorePtr[X];
			sse_XY = _mm_setzero_si128(); //SSE2
			j = 0;
		}
	}
	unsigned short Na_sseSumBuf[32], *Rp_sseStorePtr;
	Rp_sseStorePtr = (unsigned short *)alignSSE(Na_sseSumBuf);
	_mm_store_si128((__m128i *)Rp_sseStorePtr, sse_XY); //SSE2
	for (wsUint X=0 ; X<8 ; X++)
		N_ret += Rp_sseStorePtr[X];
#endif
	/* Summation using CPU */
	for (wsUint i=N_med ; i<N_col ; i++)
		N_ret += Na_1[i] * Na_2[i];

	return (wsReal)N_ret;
}

inline wsReal sseVVchar2ic(char *Na_1, wsUint N_col, char *Na_2)
{
	wsUint	N_med = 0;
	wsUint	N_ret = 0;
#ifdef USE_SSE
	w128	W_sse;
	memset(&W_sse, 0x00, sizeof(w128));
	N_med = get16(N_col);
	__m128i sse_XY		= _mm_setzero_si128();
	__m128i sse_m1		= _mm_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1);
	for (wsUint i=0,j=0 ; i<N_med ; i+=16) {
		__m128i *a = (__m128i *)(Na_1 + i);
		__m128i *b = (__m128i *)(Na_2 + i);
		__m128i aa = _mm_cmpgt_epi8(*a, sse_m1);
		__m128i bb = _mm_cmpgt_epi8(*b, sse_m1);
		sse_XY = _mm_add_epi16(sse_XY, _mm_maddubs_epi16(
			_mm_and_si128(*a, aa),
			_mm_and_si128(*b, bb))); //SSSE3
		wsUint Q = 0;
		for (wsUint k=i,l=0 ; l<16 ; l++,k++) {
			if (isMissing(Na_1[k]) || isMissing(Na_2[k])) continue;
			Q += Na_1[k] * Na_2[k];
		}
		if (++j == 256) {
			unsigned short Na_sseSumBuf[32], *Rp_sseStorePtr;
			Rp_sseStorePtr = (unsigned short *)alignSSE(Na_sseSumBuf);
			_mm_store_si128((__m128i *)Rp_sseStorePtr, sse_XY);
			for (wsUint X=0 ; X<8 ; X++)
				N_ret += Rp_sseStorePtr[X];
			sse_XY = _mm_setzero_si128(); //SSE2
			j = 0;
		}
	}
	unsigned short Na_sseSumBuf[32], *Rp_sseStorePtr;
	Rp_sseStorePtr = (unsigned short *)alignSSE(Na_sseSumBuf);
	_mm_store_si128((__m128i *)Rp_sseStorePtr, sse_XY); //SSE2
	for (wsUint X=0 ; X<8 ; X++)
		N_ret += Rp_sseStorePtr[X];
#endif
	/* Summation using CPU */
	for (wsUint i=N_med ; i<N_col ; i++) {
		if (isMissing(Na_1[i]) || isMissing(Na_2[i])) continue;
		N_ret += Na_1[i] * Na_2[i];
	}

	return (wsReal)N_ret;
}

#endif

void char2double(char *Na_i, wsReal *Ra_o, wsUint N_sz)
{

}

//using namespace svg;
void test_function()
{
#ifdef USE_FTP
	cFtp F;
	F.DoOpen("ftp.ncbi.nlm.nih.gov");
	F.DoLogin("user anonymous:biznok@gmail.com");
	F.DoCD("cd snp/organisms/human_9606/rs_fasta");
//	F.DoList("ls");
	F.DoBinary();
	F.DoGet("rs_chAltOnly.fas.gz");
#endif
#ifdef USE_R
	R().Rparse("cat('hello\n')");
	R().sendOption(OPTION());
	R().Rparse("fnlist(OPTS)");
#endif
	xOperation*	X_op	= parse("3 != 3");
	xOperand*	X_res	= understand(X_op);
	LOG("%s\n", X_res->toString());

	{
		/* Build two 8 shorts from 16 chars */
		/* Build four 4 integers from two 8 shorts */
		/* Convert four 4 integers into eight 2 doubles */
		wsUint Na_szs[] = { 1, 3, 5, 7, 9, 10, 11, 12, 15, 20, 25, 30, 32,
			33, 35, 1000, 9999, 9876540 };
		cTimer v;
		wsUint N_sz = (wsUint)len(Na_szs, wsUint);
		for (wsUint i=0 ; i<N_sz ; i++) {
			wsUint	N_szcur = Na_szs[i];
			LOG("Testing size [%d]...\n", N_szcur);
			char* Na_vec = NULL;
			sseMalloc(Na_vec, char, N_szcur);
			wsFloat *Ra_resT = NULL, *Ra_resC = NULL;
			sseMalloc(Ra_resT, wsFloat, N_szcur);
			sseMalloc(Ra_resC, wsFloat, N_szcur);

			for (wsUint j=0 ; j<N_szcur ; j++) {
				Na_vec[j] = (char)wsAdvrand(128);
			}
			v.start();
			for (wsUint j=0 ; j<N_szcur ; j++) {
				Ra_resT[j] = (wsFloat)Na_vec[j];
			}
			LOG("non-SSE %s\n", v.getReadable());
			v.start();
			//char2double(Na_vec, Ra_resC, N_szcur);
			LOG("SSE %s\n", v.getReadable());
		}
	}

	{
		wsUint Na_szs[] = { 1, 3, 5, 7, 9, 10, 11, 12, 15, 20, 25, 30, 32,
			33, 35, 1000, 9999, 9876540 };
		cTimer v;
		wsUint N_sz = (wsUint)len(Na_szs, wsUint);
		for (wsUint i=0 ; i<N_sz ; i++) {
			wsUint N_szcur = Na_szs[i];
			LOG("Testing size [%d]...\n", N_szcur);
			wsReal *Ra_vec = sseVector(N_szcur);
			wsFloat *Ra_resT = NULL, *Ra_resC = NULL;
			sseMalloc(Ra_resT, wsFloat, N_szcur);
			sseMalloc(Ra_resC, wsFloat, N_szcur);

			for (wsUint j=0 ; j<N_szcur ; j++) {
				Ra_vec[j] = wsUnifrand();
			}
			v.start();
			for (wsUint j=0 ; j<N_szcur ; j++) {
				Ra_resT[j] = (wsFloat)Ra_vec[j];
			}
			LOG("non-SSE %s\n", v.getReadable());
			v.start();
			sse2corr(Ra_vec, Ra_resC, N_szcur);
			LOG("SSE %s\n", v.getReadable());

			for (wsUint j=0 ; j<N_szcur ; j++)
				if (Ra_resT[j] != Ra_resC[j]) {
					halt("Validation failed at size [%d]", N_szcur);
				}
			sseFree(Ra_vec);
			sseFree(Ra_resT);
			sseFree(Ra_resC);
		}
	}

	wsUint N_path = 0;
	char **Sa_paths = loadStringValues("~/a,~bin/aa,~wkekf", &N_path);
	for (wsUint i=0 ; i<N_path ; i++) {
		cStrFile x(Sa_paths[i], "xx", 1);
		printf("Open [%s] [%s]\n", Sa_paths[i], x.getPath());
	}

	// Demo page shows sample usage of the Simple SVG library.

// 	Dimensions dimensions(100, 100);
// 	Document doc("my_svg.svg", Layout(dimensions, Layout::BottomLeft));
// 
// 	// Red image border.
// 	svgPolygon border(Stroke(1, Color::Red));
// 	border << Point(0, 0) << Point(dimensions.width, 0)
// 		<< Point(dimensions.width, dimensions.height) << Point(0, dimensions.height);
// 	doc << border;
// 
// 	// Long notation.  Local variable is created, children are added to varaible.
// 	LineChart chart(5.0);
// 	svgPolyline polyline_a(Stroke(.5, Color::Blue));
// 	svgPolyline polyline_b(Stroke(.5, Color::Aqua));
// 	svgPolyline polyline_c(Stroke(.5, Color::Fuchsia));
// 	polyline_a << Point(0, 0) << Point(10, 30)
// 		<< Point(20, 40) << Point(30, 45) << Point(40, 44);
// 	polyline_b << Point(0, 10) << Point(10, 22)
// 		<< Point(20, 30) << Point(30, 32) << Point(40, 30);
// 	polyline_c << Point(0, 12) << Point(10, 15)
// 		<< Point(20, 14) << Point(30, 10) << Point(40, 2);
// 	chart << polyline_a << polyline_b << polyline_c;
// 	doc << chart;
// 
// 	// Condensed notation, parenthesis isolate temporaries that are inserted into parents.
// 	doc << (LineChart(Dimensions(65, 5))
// 		<< (svgPolyline(Stroke(.5, Color::Blue)) << Point(0, 0) << Point(10, 8) << Point(20, 13))
// 		<< (svgPolyline(Stroke(.5, Color::Orange)) << Point(0, 10) << Point(10, 16) << Point(20, 20))
// 		<< (svgPolyline(Stroke(.5, Color::Cyan)) << Point(0, 5) << Point(10, 13) << Point(20, 16)));
// 
// 	doc << Circle(Point(80, 80), 20, Fill(Color(100, 200, 120)), Stroke(1, Color(200, 250, 150)));
// 
// 	doc << Text(Point(5, 77), "Simple SVG", Color::Silver, Font(10, "Verdana"));
// 
// 	doc << (svgPolygon(Color(200, 160, 220), Stroke(.5, Color(150, 160, 200))) << Point(20, 70)
// 		<< Point(25, 72) << Point(33, 70) << Point(35, 60) << Point(25, 55) << Point(18, 63));
// 
// 	doc << svgRectangle(Point(70, 55), 20, 15, Color::Yellow);
// 
// 	doc.save();
	
	cTimer t;
	/* tau */ if (0) {
		for (wsUint Q=1 ; Q<9000 ; Q++) {
//			notice("Testing %d...", Q);
// 		Q = 38;
// 		char r1[38] = { 0,2,1,1,2,1,1,2,2,1,1,0,2,0,0,0,2,0,0,0,2,2,0,2,2,2,0,0,2,1,0,2,2,2,0,0,2,2 };
// 		char r2[38] = { 0,1,1,1,2,1,1,0,2,1,2,1,0,0,1,2,1,0,2,2,1,2,2,1,1,2,2,0,2,0,0,0,1,2,1,0,0,2 };
		char *v1 = NULL;
		char *v2 = NULL;
		MULTI_ALIGNED_MALLOC(v1, char, Q);
		MULTI_ALIGNED_MALLOC(v2, char, Q);
 		for (wsUint i=0 ; i<Q ; i++) {
// 			v1[i] = r1[i];
// 			v2[i] = r2[i];
//		}
			if (rand()%100 == 0)
				v1[i] = -9;
			else v1[i] = rand()%3;
		for (wsUint i=0 ; i<Q ; i++)
			if (rand()%100 == 0)
				v2[i] = -9;
			else v2[i] = rand()%3;
		}
		t.start();
//		wsReal x = _corTauComplete(v1, v2, Q, NULL, NULL, (char *)NULL);
		_corPearsonIncomplete(v1, v2, Q, NULL, NULL, (char *)NULL);
//		wsReal u1 = t.get();
//		noticenf("%s, ", t.getReadable());
#ifdef USE_R
// 		t.start();
// 		R().sendVector("v1", v1, Q);
// 		R().sendVector("v2", v2, Q);
// 		R().Rparse("v3 <- cor(v1, v2, use='pair')");
// 		wsReal qq;
// 		R().recvVector("v3", &qq, 1);
// 		wsReal u2 = t.get();
// //		noticenf("%s, (%g times faster)\n", t.getReadable(), u2/u1);
// 		noticenf("%d	%g\n", Q, u2/u1);
// 		//LOG("diff %g - %g = %.7g\n", x, qq, x - qq);
// 		if (fabs(x-qq) > 1e-5) {
// 			FILE *f1 = fopen("v1.txt", "w+");
// 			FILE *f2 = fopen("v2.txt", "w+");
// 
// 			for (wsUint i=0 ; i<Q ; i++) {
// 				isMissing(v1[i]) ? fprintf(f1, "NA	") : fprintf(f1, "%d	", v1[i]);
// 				isMissing(v2[i]) ? fprintf(f2, "NA	") : fprintf(f2, "%d	", v2[i]);
// 			}
// 			fclose(f1);
// 			fclose(f2);
// 
// 			halt("Failed at %d\n", Q);
// 		}
// 		R().sendVector("v1", v1, Q);
// 		R().sendVector("v2", v2, Q);
//		R().Rparse("cat(paste(cor(v1, v2, method='kendall', use='pair'), '\n'))");
#endif
//		LOG("Val %g\n", x);
		}
//		FILE *f1 = fopen("v1.txt", "w+");
//		FILE *f2 = fopen("v2.txt", "w+");
//		for (wsUint i=0 ; i<Q ; i++) {
//			isMissing(v1[i]) ? fprintf(f1, "NA	") : fprintf(f1, "%d	", v1[i]);
//			isMissing(v2[i]) ? fprintf(f2, "NA	") : fprintf(f2, "%d	", v2[i]);
//		}
//		fclose(f1);
//		fclose(f2);
	}
	/* SSE-char */ if (0) {
		cTimer t;
#ifdef USE_SSE4
#define NNN 100000003
		char *v1 = NULL;
		MULTI_ALIGNED_MALLOC(v1, char, NNN);
		char *v2 = NULL;
		MULTI_ALIGNED_MALLOC(v2, char, NNN);
		wsUint R = 0;
		for (wsUint i=0 ; i<NNN ; i++) {
			v1[i] = rand()%3;
			v2[i] = rand()%3;
		}
		t.start();
		for (wsUint i=0 ; i<NNN ; i++)
			R += v1[i] * v2[i];
		LOG("Non-SSE %d %s\n", R, t.getReadable());
		t.start();
		wsUint R2 = sseVVchar(v1, NNN, v2);
		LOG("SSE %d %s\n", R2, t.getReadable());
		t.start();
		wsUint R3 = sseVVchar2(v1, NNN, v2);
		LOG("SSE2 %d %s\n", R3, t.getReadable());
		//		sprintf("Rsse %d, Rnorm %d\n", R2, R);
		//		sprintf("Rsse %d, Rnorm %d\n", R2, R);
		DEALLOC_ALIGNED(v1);
		DEALLOC_ALIGNED(v2);
#endif
	}
	/* SSE-char2 */ {
		cTimer t;
#ifdef USE_SSE4
#define NNN 100000003
		char *v1 = NULL;
		MULTI_ALIGNED_MALLOC(v1, char, NNN);
		char *v2 = NULL;
		MULTI_ALIGNED_MALLOC(v2, char, NNN);
		wsUint R = 0;
		for (wsUint i=0 ; i<NNN ; i++) {
			if (rand()%10 == 0) {
				v1[i] = -9;
				v2[i] = rand()%3;
				continue;
			}
			if (rand()%10 == 0) {
				v1[i] = rand()%3;
				v2[i] = -9;
				continue;
			}
			v1[i] = rand()%3;
			v2[i] = rand()%3;
		}
		t.start();
		for (wsUint i=0 ; i<NNN ; i++) {
			if (isMissing(v1[i]) || isMissing(v2[i])) continue;
			R += v1[i] * v2[i];
		}
		LOG("Non-SSE %d %s\n", R, t.getReadable());
		t.start();
		wsUint R3 = sseVVchar2ic(v1, NNN, v2);
		LOG("SSE2 %d %s\n", R3, t.getReadable());
		//		sprintf("Rsse %d, Rnorm %d\n", R2, R);
		//		sprintf("Rsse %d, Rnorm %d\n", R2, R);
		DEALLOC_ALIGNED(v1);
		DEALLOC_ALIGNED(v2);
#endif
	}

	/* SSE */ {
		wsReal v = W0;
		wsReal v2 = v;
		wsMat b = sseMrand(1, 100);
		wsUint m = getMed(100);
		sse_t sv = sseSet(0.0);
		for (wsUint i=0 ; i<m ; i+=sseJmp)
			sv = sseAdd(sv, *((sse_t *)(b[0]+i)));
		sseSum(sv, v);
		for (wsUint i=m ; i<100 ; i++)
			v += b[0][i];
		for (wsUint i=0 ; i<100 ; i++)
			v2 += b[0][i];
		LOG("SSE diff = %g\n", fabs(v-v2));
	}
	/* WIS */ {
		cWisWriter a;
		wsMat b = sseMrand(100, 100);
		xWisData *c = a.root();
		xWisData *d = c->add(CT_ICOR, 0, 100, 100, ICOR_UNDEF);
		d->add(CT_MTRX, b, 100, 100, 4);
		a.write("hoho");

//		cWisReader q("res.hoho");
//		xWisData *x = q.fetch("ICOR/MTRX");
//		xAttrMTRX *xx = (xAttrMTRX *)x->Sp_attr;
	}
	/* loadStringValues */ {
		wsUint n = 0;
		char *a1 = Z("a,b,c,d,e");
		const char *t1[5] = { "a", "b", "c", "d", "e" };
		char **v1 = loadStringValues(a1, &n);
		for (wsUint i=0 ; i<5 ; i++)
			if (strcmp(t1[i], v1[i]))
				halt("loadStringValues() failed");

// 		char *a2 = Z("a[b,c],d,e");
// 		const char *t2[3] = { "a[b,c]", "d", "e" };
// 		char **v2 = loadStringValues(a2, &n);
// 		for (wsUint i=0 ; i<3 ; i++)
// 			if (strcmp(t2[i], v2[i]))
//  				halt("loadStringValues() failed");
	}
}

} // End namespace ONETOOL
