#include "analyses/tridge.h"
#include "utils/vector.h"
#include "utils/matrix.h"

namespace ONETOOL {

cStdMatrix funCvTrRidge(cVector& V_y, cStdMatrix& M_Xt, cVector& V_w,
	cVector& V_lam, cVector& V_tau, wsUintCst N_k, wsStrCst S_type="LIKE",
	wsRealCst R_eps=1e-10, wsUintCst N_maxIter=100);

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cTridgeAnalysis::cTridgeAnalysis(cIO *Cp_inpIO) : cAnalysis(Cp_inpIO)
{
}

cTridgeAnalysis::~cTridgeAnalysis()
{
}

void cTridgeAnalysis::run()
{
	// Assumes binary pheno
	if (getIO()->isContinuous()) halt_fmt(WISARD_CANT_OPT_W_SOMESTATE,
		"--tridge", "continuous phenotype");
	// Assumes single pheno
	if (getIO()->sizePheno() > 1) halt_fmt(WISARD_CANT_OPT_W_SOMESTATE,
		"--tridge", "multiple phenotype");
	wsUint N_anaSamp, N_samp = getIO()->sizeSample();
	wsRealCst* Ra_y = getIO()->getPheno();
	char* Ba_misPheno = getIO()->getPheCovMissing(&N_anaSamp);
	char** Na_geno = getIO()->getGenotype();
	wsUint N_vrt = getIO()->sizeVariant(), I = 0;
	cStdMatrix M_Xt(N_vrt+1, N_anaSamp);
	wsMat Ra_Xt = M_Xt.get();
	cVector V_anaY(N_anaSamp);
	wsVec Ra_anaY = V_anaY.get();
	M_Xt.r2v_ptr(0).set(W1);
	LOOP (i, N_samp) {
		if (Ba_misPheno[i]) continue;
		LOOP (j, N_vrt)
			Ra_Xt[j+1][I] = Na_geno[i][j];
		Ra_anaY[I++] = Ra_y[i];
	}
	if (I != N_anaSamp) halt("SYSERR: N_anaSamp does not match");
	// w.vec <- c(rep(0,p0.num-s.num),rep(1,s.num))
	cVector V_w(N_vrt+1);
	V_w.set(W1);
	V_w.get()[0] = W0;
	// lam.vec <- exp(seq(log(1e+3),log(1e-3),length.out=100))
	cVector V_lam = seq(1e+3, 1e-3, 100, SEQ_LOG, SEQ_EXP);
	// tau.vec <-exp(seq(log(1e+0), log(1e-3), length.out=100))
	cVector V_tau = seq(1e+0, 1e-3, 100, SEQ_LOG, SEQ_EXP);
	// FIXME: missing treatmen
	// cv.tr.ridge <-cv.tr.ridge.fun(y.vec, x.mat, w.vec, lam.vec, tau.vec, k.num=5, eps=1e-10, iter.max=1e+2)
	cStdMatrix M_res = funCvTrRidge(V_anaY, M_Xt, V_w, V_lam, V_tau, 5);
}

#endif

cVector seq(wsReal R_s, wsReal R_e, wsUint N_len, xSeqType X_in, xSeqType X_out)
{
	/* Input value processing */
	switch (X_in) {
	case SEQ_ASIS: break;
	case SEQ_LOG:
		R_s = log(R_s);
		R_e = log(R_e);
		break;
	case SEQ_EXP:
		R_s = exp(R_s);
		R_e = exp(R_e);
		break;
	}
	wsReal R_intv = (R_e - R_s) / ((wsReal)N_len - W1);
	wsVec Ra_out = sseVector(N_len);
	wsReal R_cur = R_s;
	for (wsUint i=0 ; i<N_len ; i++) {
		switch (X_out) {
		case SEQ_ASIS:
			Ra_out[i] = R_cur; break;
		case SEQ_LOG:
			Ra_out[i] = log(R_cur); break;
		case SEQ_EXP:
			Ra_out[i] = exp(R_cur); break;
		}
		R_cur += R_intv;
	}
	return cVector(Ra_out, N_len);
}

/*
  #############################################################################################################    
  #x.mat, w.vec, b.vec must include intercetpt!!!!!
  #ridge.fun returns a vector of parameter for a fixed lamda!!!! 
*/
void funRidge(cVector& V_y, cStdMatrix& M_Xt,cVector& V_w, cVector& V_b,
	wsVec Ra_newB, wsRealCst R_lam, wsUint* Na_sset, wsUintCst N_sset,
	wsRealCst R_eps=1e-10, wsUint N_maxIter=1e+2)
// ridge.fun <- function(y.vec,x.mat,w.vec,b.vec,lam,eps=1e-10,iter.max=1e+2){# func
{
#ifdef TRIDGEvalidate
	V_w.file("w");
	V_y.file("y");
	M_Xt.file("Xt");
#endif
	// N_vrt = p.num
	// N_anaSamp = n.num
	// p.num <- dim(x.mat)[2]; n.num <- dim(x.mat)[1]
	// N_maxIter = iter.max
	// M_Xt = t(x.mat)
	// #s = length(s.set)
    // n.set <- w.vec==0; s.set <- w.vec!=0
	cVector V_ws = V_w.divide(Na_sset, N_sset);
	wsUint N_vrt = M_Xt.row();
	wsUint N_anaSamp = M_Xt.col();
	
    // cost.vec <- rep(0,iter.max); #grad.vec <- rep(0,p.num)    
	wsVec Ra_cost = sseEmptyVec(N_maxIter);
    // new.b.vec <- rep(0,p.num)
	cVector V_nb(N_vrt);
    
    // for(iter in 1:iter.max){# iter 
	LOOP (i, N_maxIter) {
		cVector V_Xb = V_b * M_Xt; // [np][p-] = [n-]
		V_Xb.minSelf(1e+10);
		// xb.vec <- pmin(drop(x.mat%*%b.vec),1e+10) 
		
		cVector V_b2 = V_b.sq(V_w);
		wsReal R_cost2 = V_b2.sum() * R_lam;
		cVector V_XbLgt = V_Xb.log1exp();
		cVector V_yXb = V_Xb * V_y;
		cVector V_cost1 = V_XbLgt - V_yXb;
		Ra_cost[i] = V_cost1.sum() + R_cost2; // [p-]
		// cost.vec[iter] <- sum(-y.vec*xb.vec+log(1+exp(xb.vec)))+lam*sum(w.vec*b.vec^2) 
    
		cVector V_mu = V_Xb.logit(); // [n-]
		// mu.vec <- exp(xb.vec)/(1+exp(xb.vec))
		cVector V_d = V_mu.p1p(); // [n-]
		V_d.maxSelf(1e-10);
		// d.vec <- pmax(mu.vec*(1-mu.vec),1e-10)
		cDiagMatrix M_d(V_d); // [pp]
		// d.mat <- diag(d.vec)
		cDiagMatrix M_dInv = M_d.inv(); // [nn]
		// inv.d.mat <- diag(1/d.vec)
		cDiagMatrix M_dSqrt = M_d.sqrt(); // [nn]
		// sqr.d.mat <- diag(sqrt(d.vec))
		
		cStdMatrix M_nXt;
		cStdMatrix M_sXt = M_Xt.divide(Na_sset, N_sset, &M_nXt);
		//cStdMatrix M_nXt = M_Xt.subset(Na_idxN); // [nS]'
		//cStdMatrix M_sXt = M_Xt.subset(Na_idxS); // [ns]'
		
		cStdMatrix M_dnXt = M_nXt * M_dSqrt; // [nn][nS]' = [nS]
		cStdMatrix M_dsXt = M_sXt * M_dSqrt; // [nn][ns]' = [ns]
		// pverbose("sum(dnXt) %g sum(dsXt) %g\n", M_dnXt.sum(), M_dsXt.sum());
      
		// xn.mat <- sqr.d.mat%*%x.mat[,n.set] #tilde version of X_n
		// xs.mat <- sqr.d.mat%*%x.mat[,s.set] #tilde version of X_s   
		cVector V_ym = V_y - V_mu;
		cVector V_ym2 = M_dInv * V_ym;
		V_ym2 += V_Xb;
		cVector V_s = M_dSqrt * V_ym2; // [nn]([n-]([nn][n-]=[n-]))=[n-]
		// s.vec <- sqr.d.mat%*%(xb.vec+drop(inv.d.mat%*%(y.vec-mu.vec))) #tilde version of s_0 
      
		cSymMatrix M_dnXt2 = M_dnXt.Mt();
		cSymMatrix& M_t0 = M_dnXt2.inv(); // [Sn][nS] = [SS]
		// t0.mat <- solve(t(xn.mat)%*%xn.mat)
		// cIdtMatrix M_id(N_anaSamp);
		// id.mat <- diag(rep(1,n.num))
		cStdMatrix M_dnX = M_dnXt.transpose();
		cSymMatrix M_p = M_dnX.MMt(M_t0); // [nn]-[nS][SS][Sn] = [nn]
		M_p.subIself();
		// p.mat <- id.mat-xn.mat%*%t0.mat%*%t(xn.mat) #tilde version of pi
      
		cVector V_newB;
		if (N_vrt < N_anaSamp) {
			// if(p.num<n.num){
			cStdMatrix M_t1 = M_dsXt * M_p; // [sn][nn] = [sn]
			// t1.mat <- t(xs.mat)%*%p.mat 
		
			cStdMatrix M_b0 = M_t1.Mt(M_dsXt);
			cDiagMatrix M_s(V_ws);
			M_s *= W2 * R_lam;
			M_b0 += M_s;
			cStdMatrix& M_b0inv = M_b0.inv();

			cVector V_b1 = M_t1 * V_s;
			V_newB = M_b0inv * V_b1;
		
        
			// ([sn][ns] + [ss]=[ss]) [sn][n-]=[s-] = [s-]
			// new.b.vec[s.set] <- solve(t1.mat%*%xs.mat+2*lam*diag(w.vec[s.set]))%*%t1.mat%*%s.vec
		} else {
			cDiagMatrix M_s(V_ws);
			cDiagMatrix& M_sInv = M_s.inv();
			// [ss][sn][nn] = [sn]
			cStdMatrix M_t21 = M_sInv * M_dsXt;
			// [sn][nn] = [sn]
			cStdMatrix M_t2 = M_t21 * M_p;
			// t2.mat <- diag(1/w.vec[s.set])%*%t(xs.mat)%*%p.mat
			// [sn] ([nn][ns][sn]+[nn]=[nn]) [n-] = [s-]
			// [sn][nn] = [nn]
			cStdMatrix M_b1 = M_dsXt * M_p;
			// [nn][ns] = [ns]
			cStdMatrix M_b2 = M_t2.tM(M_b1);
			M_b2.addDiag(W2+R_lam);
			cStdMatrix& M_b2inv = M_b2.inv();
			cVector V_b3 = M_b2inv * V_s;
			V_newB = M_t2 * V_b3;
			delete &M_b2inv;
			// new.b.vec[s.set] <- t2.mat%*%solve(p.mat%*%xs.mat%*%t2.mat+2*lam*id.mat)%*%s.vec
		}

		// [SS][Sn] ([n-] - [ns][s-] = [n-]) = [S-]
		cVector V_subS2 = V_newB * M_dsXt;
		V_s -= V_subS2;
		cVector V_newB1 = M_dnXt * V_s;
		cVector V_newB2 = M_t0 * V_newB1;
	//	new.b.vec[n.set] <- t0.mat%*%t(xn.mat)%*%(s.vec-drop(xs.mat%*%new.b.vec[s.set]))
		V_nb.replace(Na_sset, N_sset, V_newB.get());
		V_nb.replace(Na_sset, N_sset, V_newB2.get(), 1);

		cVector V_diffB = V_nb - V_b;
		wsReal R_diff = V_diffB.asum();
//		pverbose("Diff %g\n", R_diff);
		if (R_diff < R_eps) {
			pverbose("Done %d \n", i);
			break;
		}

		memcpy(V_b.get(), V_nb.get(), sizeof(wsReal)*N_vrt);
		// b.vec <- new.b.vec         
    } // iter 

//    cost.vec <- cost.vec[1:iter]
//    return(list(coef=b.vec,cost=cost.vec))     
    memcpy(Ra_newB, V_b.get(), sizeof(wsReal)*N_vrt);
 } // func
 
 /*
  #############################################################################################################      
  #x.mat, w.vec must include intercetpt !!!!!
  #ridge.path.fun returns a matrix of parameters for a sequence of lambdas!!!!   
 */

 cStdMatrix funRidgePath(cVector& V_y, cStdMatrix& M_Xt, cVector& V_w,
	 cVector& V_lam, wsUint* Na_sset, wsUint N_sset,
	 wsRealCst R_eps=1e-10, wsUintCst N_maxIter=100)
//  ridge.path.fun <- function(y.vec,x.mat,w.vec,lam.vec,eps=1e-10,iter.max=1e+2){# func
{
//   V_y   = y.vec
//   M_Xt  = t(x.mat)
//   V_lan = lam.vec
//   N_vrt = p.num

	wsUint N_lam = V_lam.size();
	wsVec Ra_lam = V_lam.get();
//    lam.num <- length(lam.vec)
//    p.num <- dim(x.mat)[2]
	wsUintCst N_vrt = M_Xt.row();

	wsMat		Ra_b = sseMatrix(N_lam, N_vrt);
//    b.mat <- matrix(0,ncol=lam.num,nrow=p.num)
	cVector V_b(N_vrt);
//    b.vec <- rep(0,p.num)

	LOOP (i, N_lam) {
		pverbose("Do lambda [%g]\n", Ra_lam[i]);
//    for(lam.pos in 1:lam.num){# lam
		funRidge(V_y, M_Xt, V_w, V_b, Ra_b[i], Ra_lam[i], Na_sset,
			N_sset, R_eps, N_maxIter);
	} // lam
//      b.vec <- ridge.fun(y.vec,x.mat,w.vec,b.vec,lam.vec[lam.pos],eps,iter.max)$coef
//      b.mat[,lam.pos] <- b.vec 
//    }# lam
    
	// [lp]
	return cStdMatrix(N_lam, N_vrt, Ra_b);
//    return(list(coef=b.mat,lam.vec=lam.vec))    
} //# func

/*
  #############################################################################################################      
  #x.mat, w.vec must include intercetpt !!!!!
  #tr.ridge.path.fun returns a matrix of parameters for a sequence of lambdas and a fixed tau!!!!     
*/
cStdMatrix funTrRidgePath(cVector& V_y, cStdMatrix& M_Xt, cVector& V_w,
	cVector& V_lam, wsRealCst R_tau, wsUint* Na_sset, wsUint N_sset,
	wsRealCst R_eps=1e-10, wsUintCst N_maxIter=100)
{
//  tr.ridge.path.fun <- function(y.vec,x.mat,w.vec,lam.vec,tau,eps=1e-10,iter.max=1e+2){
	cStdMatrix M_coef = funRidgePath(V_y, M_Xt, V_w, V_lam, Na_sset, N_sset,
		R_eps, N_maxIter);
//    tr.ridge.path <- ridge.path.fun(y.vec,x.mat,w.vec,lam.vec,eps,iter.max)$coef
//    tr.ridge.path[abs(tr.ridge.path)<tau] <- 0
	M_coef.trcLess(R_tau);
	M_coef.setClean(MATDEL_NONE);

	// [lp]
	return cStdMatrix(M_coef.row(), M_coef.col(), M_coef.get());
//    return(list(coef=tr.ridge.path,lam.vec=lam.vec,tau=tau))
}

cStdMatrix funTrRidgePath(cStdMatrix& M_ridgePath, wsRealCst R_tau)
{
	// b.mat <- ridge.path
	// b.mat[-1,][abs(b.mat[-1,])<tau] <- 0 # truncation except for the intercept 
	wsMat Ra_ret = sseMatrixP(M_ridgePath.row(), M_ridgePath.col(),
		M_ridgePath.get());
	cStdMatrix M_coef(M_ridgePath.row(), M_ridgePath.col(), Ra_ret);
	M_coef.trcLess(R_tau, 1);
	M_coef.setClean(MATDEL_NONE);

	// [lp]
	return cStdMatrix(M_coef.row(), M_coef.col(), M_coef.get());
	//    return(list(coef=tr.ridge.path,lam.vec=lam.vec,tau=tau))
}

inline wsUint fround(wsReal v)
{
	return (wsUint)(v + REAL_CONST(0.5));
}

/*
  #############################################################################################################      
  #x.mat, w.vec must include intercetpt !!!!!
  #cv.tr.ridge.path.fun returns a vector of parameters for a sequence of lambdas and taus!!!!     
*/
cStdMatrix funCvTrRidge(cVector& V_y, cStdMatrix& M_Xt, cVector& V_w,
	cVector& V_lam, cVector& V_tau, wsUintCst N_k, wsStrCst S_type/*="LIKE"*/,
	wsRealCst R_eps/*=1e-10*/, wsUintCst N_maxIter/*=100*/)
//  cv.tr.ridge.fun <- function(y.vec,x.mat,w.vec,lam.vec,tau.vec,k.num=5,e.type="LIKE",eps=1e-10,iter.max=1e+2){# func 
{
	wsUintCst N_anaSamp = M_Xt.col();
	wsUint N_anaCase = 0, N_anaCtrl = 0;
	wsUint* Na_idxCase = NULL;
	wsUint* Na_idxCtrl = NULL;
	if (N_k > 1) {
		wsVec Ra_anaY = V_y.get();
		wsUint N_cv = N_k;
		/* FInd out # of case/ctrl & record their pos */
		N_anaCase = N_anaCtrl = 0;
		wsAlloc(Na_idxCase, wsUint, N_anaSamp);
		wsAlloc(Na_idxCtrl, wsUint, N_anaSamp);
		LOOP(i, N_anaSamp) {
			if (Ra_anaY[i] == WISARD_AFFECTED)
				Na_idxCase[N_anaCase++] = i;
			else if (Ra_anaY[i] == WISARD_UNAFFECTED)
				Na_idxCtrl[N_anaCtrl++] = i;
			else halt("Phenotype of [%d]th sample is not case nor control!",
				i);
		}
		/* #--cv should not exceed #case/#ctrl */
		if (N_cv>N_anaCase || N_cv>N_anaCtrl)
			halt("# of cross-validation [%d] cannot exceed the number of "
				"case[%d] or control[%d]!", N_cv, N_anaCase, N_anaCtrl);
	}
	wsUintCst N_lam = V_lam.size();
	wsUintCst N_tau = V_tau.size();

	wsUint N_sset = 0;
	wsVec Ra_w = V_w.get();
	LOOP (i, V_w.size()) if (Ra_w[i] != W0) N_sset++;
	wsUint* Na_sset = NULL;
	wsAlloc(Na_sset, wsUint, N_sset);
	N_sset = 0;
	LOOP (i, V_w.size()) if (Ra_w[i] != W0) Na_sset[N_sset++] = i;

//    n.num <- dim(x.mat)[1]
//    lam.num <- length(lam.vec)
//    tau.num <- length(tau.vec)
    
//    k0.list <- split((1:n.num)[y.vec==0],1:k.num) # warnings may happen 
//    k1.list <- split((1:n.num)[y.vec==1],1:k.num) # warnings may happen

	// mat <- matrix(0,nrow=p.num,ncol=lam.num)
	vector<cStdMatrix> Mv_ridgeList;
	// ridge.list <- list(mat) 
	// for (k.pos in 2:k.num){ ridge.list[[k.pos]] <-mat }

	// for (k.pos in 1:k.num){ #  fold
	LOOP (k, N_k) {
		LOG("Optimizing fold [%d]\n", k);
		vInt Na_testIdxCase, Na_testIdxCtrl;
		wsUint	N_sCase = fround((wsReal)N_anaCase / (wsReal)N_k * (wsReal)k);
		wsUint	N_eCase = fround((wsReal)N_anaCase / (wsReal)N_k * (wsReal)(k+1));
		wsUint	N_sCtrl = fround((wsReal)N_anaCtrl / (wsReal)N_k * (wsReal)k);
		wsUint	N_eCtrl = fround((wsReal)N_anaCtrl / (wsReal)N_k * (wsReal)(k+1));

		/* Insert [s,e) case/ctrl */
		for (wsUint x=N_sCase ; x<N_eCase ; x++) Na_testIdxCase.push_back(Na_idxCase[x]);
		for (wsUint x=N_sCtrl ; x<N_eCtrl ; x++) Na_testIdxCtrl.push_back(Na_idxCtrl[x]);

		// Build N_curTe
		wsUint N_te = (wsUint)(Na_testIdxCase.size() + Na_testIdxCtrl.size());
		wsUint* Na_curTe = NULL;
		wsAlloc(Na_curTe, wsUint, N_te);
		wsUint N_idxTe = 0;
		FOREACH(vInt_it, Na_testIdxCase, ii) Na_curTe[N_idxTe++] = *ii;
		FOREACH(vInt_it, Na_testIdxCtrl, ii) Na_curTe[N_idxTe++] = *ii;

		// for(k.pos in 1:k.num){# fold
		// ts.set <- c(k0.list[[k.pos]],k1.list[[k.pos]])
		cVector V_trY;
		cVector V_teY = V_y.divide(Na_curTe, N_te, &V_trY);
		cStdMatrix M_trXt;
		cStdMatrix M_teXt = M_Xt.divideCol(Na_curTe, N_te, &M_trXt);

		// [lp]
		// ridge.list[[k.pos]] <- ridge.path.fun(y.vec[-ts.set],x.mat[-ts.set,],w.vec,lam.vec,eps,iter.max)$coef
		Mv_ridgeList.push_back(funRidgePath(V_trY, M_trXt, V_w, V_lam,
			Na_sset, N_sset, R_eps, N_maxIter));
	}

	// NL.cv.mat <- matrix(0,lam.num,tau.num)
	// PE.cv.mat <-matrix(0, lam.num, tau.num)
	// AC.cv.mat <-matrix(0, lam.num, tau.num)
	cStdMatrix M_cvNL(N_tau, N_lam);
	cStdMatrix M_cvPE(N_tau, N_lam);
	cStdMatrix M_cvAC(N_tau, N_lam);

	// cat("truncating ridges...","\n")
	// for(tau.pos in 1:tau.num){# tau 
	LOOP (i, N_tau) {
		// NL.mat <-matrix(0, k.num, lam.num)
		// PE.mat <-matrix(0, k.num, lam.num)
		// AC.mat <-matrix(0, k.num, lam.num)
		cStdMatrix M_NL(N_k, N_lam);
		cStdMatrix M_PE(N_k, N_lam);
		cStdMatrix M_AC(N_k, N_lam);
		wsMat	Ra_NL = M_NL.get();
		wsMat	Ra_PE = M_PE.get();

		LOOP (k, N_k) {
			vInt Na_testIdxCase, Na_testIdxCtrl;
			wsUint	N_sCase = fround((wsReal)N_anaCase / (wsReal)N_k * (wsReal)k);
			wsUint	N_eCase = fround((wsReal)N_anaCase / (wsReal)N_k * (wsReal)(k+1));
			wsUint	N_sCtrl = fround((wsReal)N_anaCtrl / (wsReal)N_k * (wsReal)k);
			wsUint	N_eCtrl = fround((wsReal)N_anaCtrl / (wsReal)N_k * (wsReal)(k+1));

			/* Insert [s,e) case/ctrl */
			for (wsUint x=N_sCase ; x<N_eCase ; x++) Na_testIdxCase.push_back(Na_idxCase[x]);
			for (wsUint x=N_sCtrl ; x<N_eCtrl ; x++) Na_testIdxCtrl.push_back(Na_idxCtrl[x]);

			// Build N_curTe
			wsUint N_te = (wsUint)(Na_testIdxCase.size() + Na_testIdxCtrl.size());
			wsUint* Na_curTe = NULL;
			wsAlloc(Na_curTe, wsUint, N_te);
			wsUint N_idxTe = 0;
			FOREACH (vInt_it, Na_testIdxCase, ii) Na_curTe[N_idxTe++] = *ii;
			FOREACH (vInt_it, Na_testIdxCtrl, ii) Na_curTe[N_idxTe++] = *ii;

			// for(k.pos in 1:k.num){# fold
			// ts.set <- c(k0.list[[k.pos]],k1.list[[k.pos]])
			cVector V_trY;
			cVector V_teY = V_y.divide(Na_curTe, N_te, &V_trY);
			cStdMatrix M_trXt;
			cStdMatrix M_teXt = M_Xt.divideCol(Na_curTe, N_te, &M_trXt);

			// [lp]
			cStdMatrix M_tr = funTrRidgePath(Mv_ridgeList[k], V_tau.get()[i]);
			// tr.mat <- tr.ridge.path.fun(y.vec[-ts.set],x.mat[-ts.set,],w.vec,lam.vec,tau.vec[tau.pos],eps,iter.max)$coef

			// xb.mat <- x.mat[ts.set,]%*%tr.mat
			// [np][pl] = [nl]
			// [pn][lp] = [ln]
			cStdMatrix M_Xb = M_tr * M_teXt;
			{
//				if(e.type=="NL"){# if else 
				cStdMatrix	M_err2 = M_Xb.log1exp();
				// -y.vec[ts.set]*xb.mat (column-wise product)
				cStdMatrix	M_err1 = V_teY.colP(M_Xb); // [lN]
				cStdMatrix	M_err = M_err2 - M_err1; // [lN]
				cVector		V_csErr = M_err.sumR(); // [l]
				if (N_lam != V_csErr.size()) halt("SYSERR: Dim err");
				memcpy(Ra_NL[k], V_csErr.get(), sizeof(wsReal)*N_lam);

				// [Np][pl] = [Nl]
				//xb.mat <- x.mat[ts.set,]%*%tr.mat

				// [Nl] [Nl]
				//err.mat[k.pos,] <- colSums(-y.vec[ts.set]*xb.mat+log(1+exp(xb.mat)))# a vector times each column of a matrix 
			}

			{
				wsMat Ra_Xb = M_Xb.get();
				//wsVec Ra_teY = V_teY.get();
				memset(Ra_PE[k], 0x00, sizeof(wsReal)*M_Xb.row());
				LOOP (ix, M_Xb.row()) LOOP (iy, M_Xb.col())
					if (Ra_Xb[ix][iy] > W0) Ra_PE[k][ix] += W1; 
				// pr.mat <- x.mat[ts.set,]%*%tr.mat>0 # a logical matrix (T=1, F=0)
				// err.mat[k.pos,] <- colSums(pr.mat==y.vec[ts.set])# a vector is equal to each column of a matrix 
			}

//			LOOP (i, N_lam) {
				// AC.mat[k.pos,lam.pos] <- performance.auc(y.vec[ts.set],pr.mat[,lam.pos])
//			}
		} // fold

		// NL.cv.mat[, tau.pos] <-colMeans(NL.mat)
		// PE.cv.mat[, tau.pos] <-colMeans(PE.mat)
		// AC.cv.mat[, tau.pos] <-colMeans(AC.mat)
		cVector V_mNL = M_NL.meanC();
		cVector V_mPE = M_PE.meanC();
		cVector V_mAC = M_AC.meanC();
		memcpy(M_cvNL.get()[i], V_mNL.get(), sizeof(wsReal)*N_lam);
		memcpy(M_cvPE.get()[i], V_mPE.get(), sizeof(wsReal)*N_lam);
		memcpy(M_cvAC.get()[i], V_mAC.get(), sizeof(wsReal)*N_lam);
	} // tau

	// NL.opt.pos <- which(NL.cv.mat==min(NL.cv.mat),arr.ind=TRUE)[1,]# row=lam, col=tau
	wsUint N_idxOptTau;
	wsUint N_idxOptLam;
	M_cvNL.getMin(&N_idxOptTau, &N_idxOptLam);
	cVector V_subLam = V_lam.subset(0, N_idxOptLam);
	cStdMatrix M_bNLlam = funRidgePath(V_y, M_Xt, V_w, V_subLam,
		Na_sset, N_sset, R_eps, N_maxIter);
	cStdMatrix M_bNLtau = funTrRidgePath(M_bNLlam, V_tau.get()[N_idxOptTau]);
    /*
    plot(as.vector(cv.err.mat))
    opt.pos <- which(cv.err.mat==min(cv.err.mat),arr.ind=TRUE)[1,]# row=lam, col=tau
	*/

	/*
    lam.vec <- lam.vec[1:opt.pos[1]]
    tau <- tau.vec[opt.pos[2]]
    lam <- lam.vec[opt.pos[1]]
    */
	wsRealCst R_tauOpt = V_tau.get()[N_idxOptTau];
	//wsRealCst R_lamOpt = V_lam.get()[N_idxOptLam];
	cStdMatrix M_ret = funTrRidgePath(V_y, M_Xt, V_w, V_lam, R_tauOpt,
		Na_sset, N_sset, R_eps, N_maxIter);
	M_ret.setClean(MATDEL_NONE);
	return cStdMatrix(M_ret.row(), M_ret.col(), M_ret.get());
//    tr.mat <- tr.ridge.path.fun(y.vec,x.mat,w.vec,lam.vec,tau,eps,iter.max)$coef
//    b.vec <- tr.mat[,opt.pos[1]]
//    return(list(coef=b.vec,lam=lam,tau=tau))  
} // func 

/*
  #############################################################################################################      
  ### opttional 
  
  plot.path.fun <- function(B.mat,L.vec){
    plot(L.vec,B.mat[1,],ylim=c(min(B.mat),max(B.mat)),type="n",main="path graph",ylab="coefficietns",xlab="lambda")
    for(id in 1:dim(B.mat)[1]){lines(L.vec,B.mat[id,],col=id)}
  }
    
  #############################################################################################################      
  ### example 
  
  n.num <- 200
  p.num <- 20
  p0.num <- p.num+1
  s.num <- 8 
  
  x.mat <- matrix(rnorm(n.num*p.num),nrow=n.num)
  x.mat <- cbind(1,x.mat)
  
  b.vec <- 1/(1:p0.num)
  b.vec <- b.vec*(-1)^(1:p0.num)
  
  xb.vec <- pmin(drop(x.mat%*%b.vec),1e+10)
  mu.vec <- exp(xb.vec)/(1+exp(xb.vec))
  y.vec <- rbinom(n.num,1,prob=mu.vec)
  w.vec <- c(rep(0,5),rep(1,p.num+1-5))
  
  iter.max <- 100
  eps <- 1e-6

  par(mfrow=c(2,3))

  lam.vec <- exp(seq(log(1e+3),log(1e-3),length.out=50))
  ridge.path <- ridge.path.fun(y.vec,x.mat,w.vec=w.vec,lam.vec=lam.vec)$coef
  plot.path.fun(ridge.path,1:length(lam.vec))
  
  for(tau in c(1,1/2,1/4,1/8)){
    tr.ridge.path <- tr.ridge.path.fun(y.vec=y.vec,x.mat=x.mat,w.vec=w.vec,lam.vec=lam.vec,tau=tau)$coef  
    plot.path.fun(tr.ridge.path,1:length(lam.vec))  
  }
    
  tau.vec <- exp(seq(log(1),log(1e-3),length.out=20))
  cv.tr.ridge <- cv.tr.ridge.fun(y.vec,x.mat,w.vec,lam.vec,tau.vec,k.num=5,e.type="NL",eps=1e-10,iter.max=1e+2)
    
  tau.vec <- exp(seq(log(1),log(1e-3),length.out=20))
  cv.tr.ridge <- cv.tr.ridge.fun(y.vec,x.mat,w.vec,lam.vec,tau.vec,k.num=5,e.type="PE",eps=1e-10,iter.max=1e+2)
*/

} // End namespace ONETOOL
