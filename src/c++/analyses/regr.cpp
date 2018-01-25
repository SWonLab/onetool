#include "analyses/regr.h"
#include "analyses/emai.h"
#include "utils/dcdflib.h"
#include "utils/data.h"
#include "utils/matrix.h"
#include "analyses/annot.h"
#include "global/Rconn.h"
#include "analyses/setmgr.h"
#include "utils/vis.h"
#include "utils/util.h"
#include "utils/stat.h"

#define NVAR_DEF_REGR	5

namespace ONETOOL {

int thrRegression(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result);
void _regression(xAnaRegr *Xp_a);

/*
 *
 * PUBLIC FUNCTION DEFINITION
 *
 */

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

inline void NormX(int n,double* x){
	double w = W1 / sqrt(sseVV(x, n));
	sseVpC(x, w, x, n);
}

typedef enum _xGLMtype {
	GLM_BINOMIAL,
	GLM_GAUSSIAN
} xGLMtype;

void plsFitterwt(wsUint N_fac, cStdMatrix& M_Xt,
				 cVector& V_y, wsReal* wt,
				 wsReal* Ch,
				 wsMat Ra_Wh /* N_var * nf */,
				 wsMat Ra_Th_t, /* nf * N_samp */
				 wsMat Ra_Ph /* N_var * nf */,
				 wsMat Ra_Uh_t /* nf * N_samp */);

// n = nsamp
wsMat scalew(cStdMatrix& M_Xt, cVector& V_wei, char center = TRUE,
	char scale = TRUE)
{
	if (M_Xt.col() != V_wei.size()) halt("Invalid dimension in scalew");
	wsReal R_sum = V_wei.sum();
	if (V_wei.isAnyNeg() || R_sum == W0)
		halt("weights must be non-negative and not all zero");
	V_wei /= R_sum;
	cVector V_cen;
	if (center) {
		V_cen = M_Xt.sumC(NULL);
//		_center = apply(weights * X, 2, function(x) sum(x))
		V_cen *= V_wei;
	} else
		V_cen.init(M_Xt.row());
//		_center = 0

//	X <- sweep(X, 2, _center)
// 	norm <- apply(X * X * V_wei, 2, function(x) sum(x,na.rm=na.rm))
// 	norm[norm <= 1e-07 * max(norm)] <- 1
// 	if (scale)
// 		X <- sweep(X, 2, sqrt(norm), "/")
//	return(X)

	return NULL; /* FIXME : Temp result */
}

void cgplsone_fit(xGLMtype X_typeY, cVector &V_y, cStdMatrix &M_Xt, wsUint nf=2,
	wsReal *Ra_w=NULL, wsReal tol=1e-9, wsReal deps=1e-20,
	wsUint nitermax=100, wsUint scale=0)
{
/*
	varmu<- family$variance
	dlink <- family$mu.eta
	link <- family$linkfun
	linkinv <- family$linkinv
	

	binom�� ���
	link	= logit = log(p/(1-p))
	linkinv	= e(p) / (1+e(p))
	dlink	= exp(v) / (1+exp(v))^2
	varmu	= p * (1-p)

	gau�� ���
	link	= nothing
	linkinv	= nothing
*/

	wsUint	N_var	= M_Xt.row(); // N_var
	wsUint	N_samp	= M_Xt.col(); // N_samp
	wsReal*	Ra_y	= V_y.get();

	/*
	fk <- link(mustart)
	weight ������ 0 = 0.25, 1 = 0.75
	*/
	wsReal*	Ra_curY = sseVector(N_samp);
	wsReal* Ra_mu = NULL;
	switch (X_typeY) {
	case GLM_GAUSSIAN:
		memcpy(Ra_mu, Ra_y, sizeof(wsReal)*N_samp);
		memcpy(Ra_curY, Ra_y, sizeof(wsReal)*N_samp);
		break;
	case GLM_BINOMIAL:
		if (Ra_w == NULL) for (wsUint i=0 ; i<N_samp ; i++)
			Ra_mu[i] = Ra_y[i] == WISARD_AFFECTED ? REAL_CONST(0.75) : REAL_CONST(0.25);
		else for (wsUint i=0 ; i<N_samp ; i++)
			//mustart <- (weights * y + 0.5)/(weights + 1)
			Ra_mu[i] = (Ra_w[i] * Ra_y[i] + REAL_CONST(0.5)) / (Ra_w[i] + W1);

		for (wsUint i=0 ; i<N_samp ; i++)
			Ra_curY[i] = log(Ra_mu[i] / (W1 - Ra_mu[i]));
		break;
	}
	cVector V_curY(Ra_curY, N_samp, 1);
//	mu <- linkinv(fk)

	/*
	v <- (dlink(fk)^2/varmu(mu))
	*/
	wsReal*	Ra_curW = sseVector(N_samp);
	cVector	V_v(Ra_curW, N_samp);
	switch (X_typeY) {
	case GLM_BINOMIAL:
		for (wsUint i=0 ; i<N_samp ; i++) {
			wsReal R_ev		= exp(Ra_curY[i]);
			wsReal R_1ev	= W1 + R_ev;
			wsReal R_sqNum	= R_ev / SQR(R_1ev);

			Ra_curW[i] = SQR(R_sqNum) / (Ra_mu[i] * (W1 - Ra_mu[i]));
		}
		break;
	case GLM_GAUSSIAN:
		sseVinit(Ra_curW, N_samp, W1);
		break;
	}
//	MAT_t Ra_curXt = scalew(M_Xt, V_v, TRUE, scale);

	cStdMatrix M_Wh, M_Ph, M_Th_t, M_Uh_t;
// 	MAT_t Ra_Wh /* N_var * nf */,
// 	MAT_t Ra_Th_t, /* nf * N_samp */
// 	MAT_t Ra_Ph /* N_var * nf */,
// 	MAT_t Ra_Uh_t /* nf * N_samp */
	cVector		Qk(nf);
// 	Tk <- matrix(0,nrow=nr,ncol=nf)
// 	Wk <- matrix(0,nrow=nc,ncol=nf)
// 	Qk <- rep(0,nf)
// 	Pk <- matrix(0,nrow=nc,ncol=nf)
// 	Uk <- matrix(0,nrow=nr,ncol=nf)

	cVector		beta(W1, N_var);
	cVector		betaold	= beta/REAL_CONST(1000.0);
	wsReal		dbeta	= 1;
	wsUint		N_iter	= 1;
	char		B_stop	= 0;
	char		B_conv	= 0;

	while (dbeta>tol && N_iter<nitermax) {
		plsFitterwt(nf, M_Xt, V_curY, Ra_curW, Qk.get(),
			M_Wh.get(), M_Th_t.get(), M_Ph.get(), M_Uh_t.get());
// 		fit1 <- .C("plsFitterwt",
// 		as.integer(nc),as.integer(nr),as.integer(nf),as.double(Ek),as.double(fk),
// 			as.double(v),as.double(Qk),as.double(Wk),as.double(Tk),as.double(Pk),as.double(Uk),PACKAGE="plstools")
// 			Qk <- fit1[[7]]
// 		Wk <- matrix(fit1[[8]],nr=nc,nc=nf)
// 		Pk <- matrix(fit1[[10]],nr=nc,nc=nf)
// 		Tk <- matrix(fit1[[9]],nr=nr,nc=nf)
// 		Uk <- matrix(fit1[[11]],nr=nr,nc=nf)

		cVector V_eta = M_Th_t.tV(Qk);
//		eta <- as.numeric(weighted.mean(fk,v)+Tk%*%Qk) // [nr*nf] [nf*1]
		V_eta += V_curY.wMean(Ra_curW, N_samp);

		wsReal *Ra_eta = V_eta.get();
		switch (X_typeY) {
		case GLM_BINOMIAL:
		/*
		mu <- linkinv(eta)
		*/
			/* linkinv = e(p) / (1+e(p)) */
			/* dlink = exp(v) / (1+exp(v))^2 */
			for (wsUint i=0 ; i<N_samp ; i++) {
				wsReal	R_expP	= exp(Ra_eta[i]);
				wsReal	R_1expP	= W1 + R_expP;
				Ra_mu[i] = R_expP / R_1expP;

				/*
				v <- (dlink(eta)^2/varmu(mu))
				*/
				wsReal R_d1 = Ra_mu[i] / R_1expP;
				Ra_curW[i] = SQR(R_d1) / (Ra_mu[i] * (W1 - Ra_mu[i]));
				if (NA(Ra_curW[i])) {
					B_stop = 1;
					break;
				}

				/*
				fk <- eta + (1/dlink(eta))*(y-mu)
				*/
				Ra_curY[i] = (Ra_y[i] - Ra_mu[i]) / R_d1 + Ra_eta[i];
			}
			break;
		case GLM_GAUSSIAN:
		/*
		mu <- linkinv(eta)
		*/
			/* linkinv = e(p) / (1+e(p)) */
			/* dlink = exp(v) / (1+exp(v))^2 */
			for (wsUint i=0 ; i<N_samp ; i++) {
//				wsReal	R_expP	= exp(Ra_eta[i]);
				Ra_mu[i] = W1;

				/*
				v <- (dlink(eta)^2/varmu(mu))
				*/
				Ra_curW[i] = W1;
				if (NA(Ra_curW[i])) {
					B_stop = 1;
					break;
				}

				/*
				fk <- eta + (1/dlink(eta))*(y-mu)
				*/
				Ra_curY[i] = (Ra_y[i] - Ra_mu[i]) + Ra_eta[i];
			}
			break;
		}
		if (B_stop == 1) break;
//		if(sum(v=="NaN")>0)
//			break
		
//		Ra_curXt	= scalew(M_Xt, V_v, TRUE, scale);
		betaold		= beta;

		// Whx <- Wk %*% solve(crossprod(Pk, Wk))
		cStdMatrix	M_PW	= M_Ph.tM(M_Wh);
		cStdMatrix&	M_PWi	= M_PW.inv();
		cStdMatrix	M_Whx	= M_Wh * M_PWi;
		delete &M_PWi;

		// beta <- Whx %*% Qk
		beta = M_Whx * Qk;

		// dbeta <- max(abs(beta - betaold)/(abs(betaold) + deps))
		cVector V_betaDiff = beta - betaold;
		wsReal *Ra_betaDiff = V_betaDiff.get();
		dbeta = W0;
		for (wsUint i=0 ; i<N_var ; i++) {
			wsReal R_cur = fabs(Ra_betaDiff[i]) / (fabs(betaold.get()[i]) + deps);
			if (R_cur > dbeta) dbeta = R_cur;
		}

		N_iter++;
	}
	if(dbeta < tol)
		B_conv = TRUE;

	/* FIXME : Temp */
	if (!B_conv)
		LOGwarn("Failed to converge\n");

//	glm1 <- glm(y~as.matrix(Tk),family=family)
// 	Whx <- Wk %*% solve(crossprod(Pk, Wk))
// 	bh <- Whx %*% Qk
// 	mx <- apply(X, 2, function(x) mean(x, na.rm = TRUE))
// 	meta <- mean(eta,na.rm = TRUE)
// 	a <- meta - sum(mx*bh)
// 	coefx <- c(a,bh)
// 	names(coefx) <- c("(Intercept)",names(X))
// 	res <- list(bh = coefx,y = y,x = X,qh = Qk,th = as.data.frame(Tk),
// 				ph = Pk,wh = Wk,nf = nf,niter = niter,convergence = convergence,call = match.call())
}

/* plsFitter: last modified 06/14/09 15:56:49 by pbady */
// plsFitterwt 06/27/09 12:13:36
void plsFitterwt(wsUint N_fac, cStdMatrix& M_Xt,
	cVector& V_y, wsReal* wt,
	wsReal* Ch,
	wsMat Ra_Wh /* N_var * nf */,
	wsMat Ra_Th_t, /* nf * N_samp */
	wsMat Ra_Ph /* N_var * nf */,
	wsMat Ra_Uh_t /* nf * N_samp */)
{
	wsUint		N_var	= M_Xt.row();
	wsUint		N_samp	= M_Xt.col();
	wsMat		Ra_Xt	= M_Xt.get();
	if (V_y.size() != N_samp) halt("DIMERR");
	double		auxi,S;
	wsReal*		Ra_yh	= sseVector(N_samp);
	wsReal*		Ra_uh	= sseVector(N_samp);
	wsReal*		Ra_wh	= sseVector(N_var);
	cVector		V_th(N_samp), V_ph(N_var);
	wsMat		Ra_Xh	= sseMatrixP(N_var, N_samp, Ra_Xt);
	cStdMatrix	M_Xh(N_samp, N_var, Ra_Xh);
	wsReal*		Ra_xhi	= sseVector(N_var);
	memcpy(Ra_yh, V_y.get(), sizeof(wsReal)*N_samp);
	cVector		V_yh(Ra_yh, N_samp);

	for (wsUint h=0 ; h<N_fac ; h++){
		wsReal *Ra_th = V_th.get();
		wsReal *Ra_ph = V_ph.get();
		sseVinit(Ra_wh, N_var, W1/sqrt((double)N_var));
		sseVset0(Ra_xhi, N_var);

		for (wsUint i=0 ; i<N_var ; i++) {
			auxi		= V_yh.VV(Ra_Xh[i], wt, N_samp);
			S			= V_yh.V2(wt, N_samp);
			Ra_wh[i]	= auxi/S;
		}
		NormX(N_var, Ra_wh);
		for (wsUint i=0 ; i<N_samp ; i++){
			for (wsUint j=0 ; j<N_var ; j++)
				Ra_xhi[j] = Ra_Xh[j][i];
//			extractVec(N_samp,N_var,i,Ra_Xh,Ra_xhi,1);
			auxi		= sseVV(Ra_xhi, N_var, Ra_wh);
			S			= sseVV(Ra_wh, N_var, Ra_wh);
			Ra_th[i]	= auxi/S;
		}
		for (wsUint i=0 ; i<N_var ; i++) {
			auxi		= sseVVV(Ra_Xh[i], V_th.get(), wt, N_samp);
			S			= sseVVV(V_th.get(), V_th.get(), wt, N_samp);
			Ra_ph[i]	= auxi/S;
		}
		cStdMatrix M_Xt1 = V_th.tV(V_ph);
		M_Xh -= M_Xt1;
		auxi	= sseVVV(V_th.get(), Ra_yh, wt, N_samp);
		S		= sseVVV(V_th.get(), V_th.get(), wt, N_samp);
		wsReal ch = auxi/S;
		for (wsUint i=0 ; i<N_samp ; i++)
			Ra_uh[i] = Ra_yh[i]/ch;

		// results
		Ch[h] = ch;
		for (wsUint i=0 ; i<N_var ; i++) {
			Ra_Ph[i][h] = Ra_ph[i];
			Ra_Wh[i][h] = Ra_wh[i];
		}
		memcpy(Ra_Th_t[h], Ra_th, sizeof(wsReal)*N_samp);
		memcpy(Ra_Uh_t[h], Ra_uh, sizeof(wsReal)*N_samp);
	}
	sseFree(Ra_wh);
	sseFree(Ra_uh);
	sseFree(Ra_xhi);
}

class LassoFit {
	int numFeatures;
public:
	// Number of lambda values
	int numberOfLambdas;

	// Intercepts
	double* intercepts;

	// Compressed weights for each solution
	double** compressedWeights;

	// Pointers to compressed weights
	int* indices;

	// Number of weights for each solution
	int* numberOfWeights;

	// Number of non-zero weights for each solution
	int* nonZeroWeights;

	// The value of lambdas for each solution
	double* lambdas;

	// R^2 value for each solution
	double* rsquared;

	// Total number of passes over data
	int numberOfPasses;


	LassoFit(int numberOfLambdas, int maxAllowedFeaturesAlongPath, int numFeatures) {
		intercepts = sseVector(numberOfLambdas);
		compressedWeights = sseEmptyMatrix(numberOfLambdas, maxAllowedFeaturesAlongPath);
		indices = NULL;
		wsCalloc(indices, int, maxAllowedFeaturesAlongPath);
		numberOfWeights = NULL;
		wsCalloc(numberOfWeights, int, numberOfLambdas);
		lambdas = sseVector(numberOfLambdas);
		rsquared = sseVector(numberOfLambdas);
		nonZeroWeights = NULL;
		wsCalloc(nonZeroWeights, int, numberOfLambdas);
		this->numFeatures = numFeatures;
	}

	double* getWeights(int lambdaIdx) {
		double* weights = sseVector(numFeatures);
		for (int i = 0; i < numberOfWeights[lambdaIdx]; i++) {
			weights[indices[i]] = compressedWeights[lambdaIdx][i];
		}
		return weights;
	}

	string toString() {
		string sb;
		int numberOfSolutions = numberOfLambdas;
		sb.append("Compression R2 values:\n");
		for (int i = 0; i < numberOfSolutions; i++) {
			char S_buf[2048];
			sprintf(S_buf, "%d	%d	%g	%g\n", i+1, nonZeroWeights[i],
				rsquared[i], lambdas[i]);
			sb += S_buf;
//			sb.append((i + 1) + "\t" +  + "\t" + MathUtil.getFormattedDouble(rsquared[i], 4) + "\t"
//				+ MathUtil.getFormattedDouble(lambdas[i], 5) + "\n");
		}
		return sb;
	}
};

class LassoFitGenerator {
    // In order to speed up the compression, we limit the number of
    // observations,
    // but this limit is dependent on the number of features that we should
    // learn
    // their weights. In other words, for learning weights of more features, we
    // need more observations.
#define MAX_OBSERVATIONS_TO_FEATURES_RATIO 100

#define EPSILON 1.0e-6

    // The default number of lambda values to use
#define DEFAULT_NUMBER_OF_LAMBDAS 100

    // Convergence threshold for coordinate descent
    // Each inner coordination loop continues until the relative change
    // in any coefficient is less than this threshold
#define CONVERGENCE_THRESHOLD 1.0e-4

#define SMALL 1.0e-5
#define MIN_NUMBER_OF_LAMBDAS 5
#define MAX_RSQUARED 0.99999

    wsVec	Ra_targets;
    wsMat	Ra_obs;
    wsUint	N_feat;
    wsUint	N_obs;

public:
	wsUint getMaxAllowedObservations(int N_maxFeat) {
		return 100000;
	}

	void init(cVector &V_y, cStdMatrix &M_Xt) {
		wsUint N_inpObs = M_Xt.col();
		N_feat = M_Xt.row();
		N_feat--;
		if (N_inpObs > getMaxAllowedObservations(N_feat)) {
			halt("Number of observations (%d) exceeds the maximum allowed number: %d",
				N_inpObs, getMaxAllowedObservations(N_feat));
		}
		this->N_obs = N_inpObs;
		wsMat M_XtWo1 = &(M_Xt.get()[1]);
		Ra_obs = sseMatrixP(N_feat, N_obs, M_XtWo1);

		Ra_targets = sseVector(this->N_obs);
		memcpy(Ra_targets, V_y.get(), sizeof(wsReal)*N_obs);
	}

	void setNumberOfFeatures(int N_inpFeat) {
		N_feat = N_inpFeat;
	}

	void setFeatureValues(int N_idx, wsVec Ra_val, wsUint N_sz) {
		for (wsUint i=0 ; i<N_sz ; i++)
			Ra_obs[N_idx][i] = Ra_val[i];
	}

	wsVec getFeatureValues(int idx) {
		return Ra_obs[idx];
	}

	void setObservationValues(int N_idx, wsVec Ra_val) {
		for (wsUint f = 0; f < N_feat; f++) {
			Ra_obs[f][N_idx] = Ra_val[f];
		}
	}

	LassoFit* getLassoFit(int N_maxAllowedFeatPerModel) {
		cTimer	duration;
		wsReal	R_feat = (wsReal)N_feat;
		duration.start();

		if (N_maxAllowedFeatPerModel < 0) {
			N_maxAllowedFeatPerModel = N_feat;
		}
		int N_lambda = DEFAULT_NUMBER_OF_LAMBDAS;
		int maxAllowedFeaturesAlongPath = (int)MIN(N_maxAllowedFeatPerModel * (wsReal)1.2, R_feat);

		// lambdaMin = flmin * lambdaMax
		double flmin = (N_obs < N_feat ? 5e-2 : 1e-4);

		/********************************
		* Standardize features and target: Center the target and features
		* (mean 0) and normalize their vectors to have the same standard
		* deviation
		*/
		wsVec	Ra_meanFeat = sseVector(N_feat);
		wsVec	Ra_featureStds = sseVector(N_feat);
		wsVec	Ra_feat2residCorr = sseEmptyVec(N_feat);

		wsReal	R_factor = W1 / sqrt((double)N_obs);
		for (wsUint j=0 ; j<N_feat ; j++) {
			double mean = sseVsum(Ra_obs[j], N_obs) / (wsReal)N_obs;//MathUtil.getAvg(Ra_obs[j]);
			Ra_meanFeat[j] = mean;
			for (wsUint i=0 ; i<N_obs ; i++)
				Ra_obs[j][i] = R_factor * (Ra_obs[j][i] - mean);
			Ra_featureStds[j] = sqrt(sseVV(Ra_obs[j], N_obs));//sqrt(MathUtil.getDotProduct(Ra_obs[j], Ra_obs[j]));

			sseVpC(Ra_obs[j], W1/Ra_featureStds[j], Ra_obs[j], N_obs);
//			MathUtil.divideInPlace(Ra_obs[j], (float) featureStds[j]);
		}

		wsReal targetMean = sseVsum(Ra_targets, N_obs) / (wsReal)N_obs;//MathUtil.getAvg(Ra_targets);
		for (wsUint i=0 ; i < N_obs; i++) {
			Ra_targets[i] = R_factor * (Ra_targets[i] - targetMean);
		}
		float targetStd = (float)sqrt(sseVV(Ra_targets, N_obs));//sqrt(MathUtil.getDotProduct(Ra_targets, Ra_targets));
		sseVpC(Ra_targets, W1/targetStd, Ra_targets, N_obs);
//		MathUtil.divideInPlace(Ra_targets, targetStd);

		for (wsUint j=0 ; j<N_feat ; j++) {
			Ra_feat2residCorr[j] = sseVV(Ra_targets, N_obs, Ra_obs[j]);//MathUtil.getDotProduct(Ra_targets, Ra_obs[j]);
		}

		wsMat	Ra_feat2featCorr = sseMatrix(N_feat, maxAllowedFeaturesAlongPath);
		wsVec	Ra_activeWgt = sseEmptyVec(N_feat);
		int*	correlationCacheIndices = NULL;
		wsCalloc(correlationCacheIndices, int, N_feat);//new int[N_feat];
		wsVec	denseActiveSet = sseVector(N_feat);

		LassoFit *fit = new LassoFit(N_lambda, maxAllowedFeaturesAlongPath, N_feat);
		fit->numberOfLambdas = 0;

		double alf = pow(max(EPSILON, flmin), 1.0 / (N_lambda - 1));
		wsReal	R_r2 = 0.0;
		fit->numberOfPasses = 0;
		int numberOfInputs = 0;
		int minimumNumberOfLambdas = min(MIN_NUMBER_OF_LAMBDAS, N_lambda);

		double curLambda = 0;
		double maxDelta;
		for (int N_iter = 1; N_iter <= N_lambda; N_iter++) {
//			LOG("Starting iteration %d of Compression.\n", N_iter);

			/**********
			* Compute lambda for this round
			*/
			if (N_iter == 1) {
				curLambda = DBL_MAX; // first lambda is infinity
			} else if (N_iter == 2) {
				curLambda = W0;
				for (wsUint j=0 ; j<N_feat ; j++) {
					curLambda = max(curLambda, fabs(Ra_feat2residCorr[j]));
				}
				curLambda = alf * curLambda;
			} else {
				curLambda = curLambda * alf;
			}

			double prevRsq = R_r2;
			double v;
			while (true) {
				fit->numberOfPasses++;
				maxDelta = 0.0;
				for (wsUint k=0 ; k<N_feat ; k++) {
					double prevWeight = Ra_activeWgt[k];
					double u = Ra_feat2residCorr[k] + prevWeight;
					v = (u >= 0 ? u : -u) - curLambda;
					// Computes sign(u)(|u| - curLambda)+
					Ra_activeWgt[k] = (v > 0 ? (u >= 0 ? v : -v) : 0.0);

					// Is the weight of this variable changed?
					// If not, we go to the next one
					if (Ra_activeWgt[k] == prevWeight) {
						continue;
					}

					// If we have not computed the correlations of this
					// variable with other variables, we do this now and
					// cache the result
					if (correlationCacheIndices[k] == 0) {
						numberOfInputs++;
						if (numberOfInputs > maxAllowedFeaturesAlongPath) {
							// we have reached the maximum
							break;
						}
						for (wsUint j=0 ; j<N_feat ; j++) {
							// if we have already computed correlations for
							// the jth variable, we will reuse it here.
							if (correlationCacheIndices[j] != 0) {
								Ra_feat2featCorr[j][numberOfInputs - 1] = Ra_feat2featCorr[k][correlationCacheIndices[j] - 1];
							} else {
								// Correlation of variable with itself if one
								if (j == k) {
									Ra_feat2featCorr[j][numberOfInputs - 1] = 1.0;
								} else {
									Ra_feat2featCorr[j][numberOfInputs - 1] = sseVV(
										Ra_obs[j], N_obs, Ra_obs[k]);
								}
							}
						}
						correlationCacheIndices[k] = numberOfInputs;
						fit->indices[numberOfInputs - 1] = k;
					}

					// How much is the weight changed?
					double delta = Ra_activeWgt[k] - prevWeight;
					R_r2 += delta * (2.0 * Ra_feat2residCorr[k] - delta);
					maxDelta = max((delta >= 0 ? delta : -delta), maxDelta);

					for (wsUint j=0;  j<N_feat ; j++) {
						Ra_feat2residCorr[j] -= Ra_feat2featCorr[j][correlationCacheIndices[k] - 1]
						* delta;
					}
				}

				if (maxDelta < CONVERGENCE_THRESHOLD || numberOfInputs > maxAllowedFeaturesAlongPath) {
					break;
				}

				for (int ii = 0; ii < numberOfInputs; ii++) {
					denseActiveSet[ii] = Ra_activeWgt[fit->indices[ii]];
				}

				do {
					fit->numberOfPasses++;
					maxDelta = 0.0;
					for (int l = 0; l < numberOfInputs; l++) {
						int k = fit->indices[l];
						double prevWeight = Ra_activeWgt[k];
						double u = Ra_feat2residCorr[k] + prevWeight;
						v = (u >= 0 ? u : -u) - curLambda;
						Ra_activeWgt[k] = (v > 0 ? (u >= 0 ? v : -v) : 0.0);
						if (Ra_activeWgt[k] == prevWeight) {
							continue;
						}
						double delta = Ra_activeWgt[k] - prevWeight;
						R_r2 += delta * (2.0 * Ra_feat2residCorr[k] - delta);
						maxDelta = max((delta >= 0 ? delta : -delta), maxDelta);
						for (int j = 0; j < numberOfInputs; j++) {
							Ra_feat2residCorr[fit->indices[j]] -= Ra_feat2featCorr[fit->indices[j]][correlationCacheIndices[k] - 1]
							* delta;
						}
					}
				} while (maxDelta >= CONVERGENCE_THRESHOLD);

				for (int ii = 0; ii < numberOfInputs; ii++) {
					denseActiveSet[ii] = Ra_activeWgt[fit->indices[ii]] - denseActiveSet[ii];
				}
				for (wsUint j=0 ; j<N_feat ; j++) {
					if (correlationCacheIndices[j] == 0) {
						Ra_feat2residCorr[j] -= sseVV(denseActiveSet, numberOfInputs,
							Ra_feat2featCorr[j]);
					}
				}
			}

			if (numberOfInputs > maxAllowedFeaturesAlongPath) {
				break;
			}
			if (numberOfInputs > 0) {
				for (int ii = 0; ii < numberOfInputs; ii++) {
					fit->compressedWeights[N_iter - 1][ii] = Ra_activeWgt[fit->indices[ii]];
				}
			}
			fit->numberOfWeights[N_iter - 1] = numberOfInputs;
			fit->rsquared[N_iter - 1] = R_r2;
			fit->lambdas[N_iter - 1] = curLambda;
			fit->numberOfLambdas = N_iter;

			if (N_iter < minimumNumberOfLambdas) {
				continue;
			}

			int me = 0;
			for (int j = 0; j < numberOfInputs; j++) {
				if (fit->compressedWeights[N_iter - 1][j] != 0.0) {
					me++;
				}
			}
			if (me > N_maxAllowedFeatPerModel || ((R_r2 - prevRsq) < (SMALL * R_r2))
				|| R_r2 > MAX_RSQUARED) {
					break;
			}
		}

		for (int k = 0; k < fit->numberOfLambdas; k++) {
			fit->lambdas[k] = targetStd * fit->lambdas[k];
			int nk = fit->numberOfWeights[k];
			for (int l = 0; l < nk; l++) {
				fit->compressedWeights[k][l] = targetStd * fit->compressedWeights[k][l] / Ra_featureStds[fit->indices[l]];
				if (fit->compressedWeights[k][l] != 0) {
					fit->nonZeroWeights[k]++;
				}
			}
			double product = 0;
			for (int i = 0; i < nk; i++) {
				product += fit->compressedWeights[k][i] * Ra_meanFeat[fit->indices[i]];
			}
			fit->intercepts[k] = targetMean - product;
		}

		// First lambda was infinity; fixing it
		fit->lambdas[0] = exp(2 * log(fit->lambdas[1]) - log(fit->lambdas[2]));

//		LOG("Elapsed time for compression: %s\n", duration.getReadable());
		return fit;
	}

	void setTargets(wsVec Ra_inpTargets) {
		memcpy(Ra_targets, Ra_inpTargets, sizeof(wsReal)*N_obs);
	}

	void setTarget(int N_idx, wsReal R_target) {
		Ra_targets[N_idx] = R_target;
	}

	LassoFit* fit(int maxAllowedFeaturesPerModel) {
		LassoFit* fit = getLassoFit(maxAllowedFeaturesPerModel);
		return fit;
	}
};

void mmLasso(cVector& V_y, cStdMatrix& M_Xt, double R_lambda,
			 int N_maxIter, double R_thr, double *est, double *max_diff)
{
	int		N_iter		= 0;
	wsUint	N_samp		= M_Xt.col();
	wsUint	N_var		= M_Xt.row();

	double	R_diff		= 0;
	wsVec	Ra_param	= NULL;
	wsVec	Ra_oParam	= NULL;
	wsReal*	Ra_x2div;

	Ra_param		= sseVector(N_var);
	cVector	V_oldParam(N_var);
	Ra_oParam	= V_oldParam.get();
	Ra_x2div		= sseVector(N_var);

	// Initialization of parameters
	for (wsUint j=0 ; j<N_var ; j++) {
		Ra_param[j]		= W1 + wsUnifrand();
		Ra_oParam[j]	= Ra_param[j];
	}

	// Calculate alpha and sum of alpha for each row
	cVector		V_alpha	= M_Xt.asumC();
	cStdMatrix	M_alpha	= M_Xt.acDiv(V_alpha);

	wsVec		s_alpha	= sseVector(N_samp);
	wsMat		X		= M_Xt.get();
	wsMat		alpha	= sseMatrix(N_var, N_samp);
	memset(s_alpha, 0x00, sizeof(wsReal)*N_samp);
	for (wsUint i=0 ; i<N_samp ; i++) {
		for (wsUint j=0 ; j<N_var ; j++)
			s_alpha[i] += fabs(X[j][i]);

		for (wsUint j=0 ; j<N_var ; j++)
			alpha[j][i] = fabs(X[j][i])/s_alpha[i];
	}

	// Calculate sum x_ij^2 / alpha_ij
	wsMat Ra_Xt = M_Xt.get();
	for (wsUint j=0 ; j<N_var ; j++) {
		wsVec Ra_V = sseVector(N_samp);
		wsReal R_v = W0;
		sseVsquare(Ra_Xt[j], N_samp, Ra_V);
		for (wsUint i=0 ; i<N_samp ; i++) if (alpha[j][i])
			R_v += Ra_V[i] / alpha[j][i];
		Ra_x2div[j] = R_v;
	}
	cVector V_X2div(Ra_x2div, N_var, 1);
//	VEC_t Ra_X2div = V_X2div.get();

	wsVec xsq_div = sseVector(N_var);
	memset(xsq_div, 0x00, sizeof(wsReal)*N_var);
	for (wsUint j=0 ; j<N_var ; j++) {
		for (wsUint i=0 ; i<N_samp ; i++) if (alpha[j][i])
			xsq_div[j] += X[j][i]*X[j][i]/alpha[j][i];
	}

	// Start Iteration

	for (N_iter=0 ; N_iter<N_maxIter ; N_iter++) {
		cVector V_Xb = V_oldParam * M_Xt;
		cVector V_resid = V_y - V_Xb;
		// Need total param_old
// 		for(i=0;i<N_samp;i++) {
// 			resid[i] = Y[i];
// 			for(j=0;j<N_var;j++)
// 				resid[i] -= Ra_X[i][j]*param_old[j];
// 		}

		cVector V_temp = M_Xt * V_resid;
		wsReal*	temp = V_temp.get();
		for (wsUint j=0 ; j<N_var ; j++) {
// 			temp[j] = 0;
// 			for (i=0 ; i<N_samp ; i++)
// 				temp[j] += Ra_X[i][j]*resid[i];

			// Calculate parameters for each index
			Ra_param[j] = (Ra_x2div[j]*Ra_oParam[j] + temp[j]) /
				(Ra_x2div[j] + R_lambda/fabs(Ra_oParam[j]));
		}

		*max_diff = 0;
		for (wsUint j=0 ; j<N_var ; j++) {
			R_diff = fabs(Ra_param[j] -Ra_oParam[j]);
			if (R_diff >  *max_diff) *max_diff = R_diff;
		}

		if(*max_diff <= R_thr)
			break;

		// Update param_old
		for (wsUint j=0 ; j<N_var ; j++) {
			Ra_oParam[j] = Ra_param[j];
			if (Ra_oParam[j] == 0) Ra_oParam[j] += R_thr/2.0;
		}

		pverbose("Iteration = %d \n",N_iter+1);
	}

	// End Iteration

	for (wsUint j=0 ; j<N_var ; j++)
		est[j] = Ra_param[j];

	// free memory allocation
	sseFree(Ra_param);
//	sseFree(Ra_oParam);
	sseFree(Ra_x2div);
}

cRegrAnalysis::cRegrAnalysis(cIO* Cp_inpIO,
	cSetManagerAnalysis *Cp_inpAnaMgr) : cAnalysis(Cp_inpIO)
{
	wsUintCst	N_oCov		= getIO()->sizeCovar();
	wsUintCst	N_gxeCov	= getIO()->sizeGxEcovar();
	/* G2, G1*G2 */
	wsUint	N_gxg		= OPT_ENABLED(gxg) ? 2 : 0;
	N_pheno		= getIO()->sizePheno();
	B_isCont	= getIO()->isContinuous();
	Cp_anaSM	= Cp_inpAnaMgr;
	if (OPT_ENABLED(gxe))
		/* If --gxe, # of covariate will be doubled (covGxE) */
		N_col		= (N_oCov + N_gxeCov)+2+N_gxg;
	else
		N_col		= N_oCov+2+N_gxg; /* Intercept and genotype */

	/* If --gxg */
	if (OPT_ENABLED(gxg)) {
		/* Need --set or --setconsec */
		ASSERT_OPTS_AND(set, setconsec);
		/* --gxe not allowed */
		if (OPT_ENABLED(gxe)) halt_fmt(WISARD_CANT_EXCL_OPT, "--gxg", "--gxe");
	}

//	Ra_sampwgt		= NULL;
	/* Load sample weights if assigned */
	//wsUint		N_samp = getIO()->getSampleSize();
// 	if (IS_ASSIGNED(sampleweight)) {
// 		cStrFile	C_sw(OPT_STRING(sampleweight), "Sample weighting file");
// 		mSamp&		Xm_samp = getIO()->getSampleData();
// 		vSampPtr&	Xv_samp = getIO()->getSample();
// 		char*		Sp_buf = NULL;
// 		MULTI_MALLOC(Sp_buf, char, 4096);
// 
// 		Ra_sampwgt = sseVector(N_samp);
// 		sseVinit(Ra_sampwgt, N_samp, WISARD_NAN);
// 
// 		for (wsUint L=1 ; C_sw.gets(Sp_buf, 4096) ; L++) {
// 			char *a, *b, *c, *d;
// 			getString(&Sp_buf, &a);
// 			getString(&a, &b);
// 			getString(&b, &c);
// 			mSamp_it X_find = Xm_samp.find(a);
// 			if (X_find == Xm_samp.end()) continue;
// 			xSample&	X_samp = X_find->second;
// 
// 			/* Get weight */
// 			wsReal R_wgt = (wsReal)strtod(b, &d);
// 			if (d && d[0]) halt("Invalid weight value [%s] in line [%d]",
// 				c, L);
// 			/* Weight can't be zero since it takes an inverse */
// 			if (R_wgt == W0)
// 				halt("Sample weight of [%s::%s] cannot be zero",
// 					Sp_buf, a);
// 
// 			/* If the sample is not available, skip */
// 			if (X_samp.N_idx == SAMP_NODATA) continue;
// 			Ra_sampwgt[X_samp.N_idx] = W1/R_wgt;
// 		}
// 		DEALLOC(Sp_buf);
// 		/* Check all samples have their weight */
// 		for (wsUint i=0 ; i<N_samp ; i++)
// 			if (NA(Ra_sampwgt[i])) halt("Sample [%s::%s] does not have "
// 				"sample weight!", Xv_samp[i]->S_FID.c_str(),
// 				Xv_samp[i]->S_IID.c_str());
// 
// 		LOGnote("Sample weights successfully loaded\n");
// 	}

	/* Remove individuals having missing phenotypes or covariates */
	Ba_isExcl	= getIO()->getPheCovMissing();

	/* Resize phenotype and covariates */
//	N_anaSamp	= _makeY(Ba_isExcl, &Ra_anaY);//, &Ra_anaCov);
	N_anaSamp	= _makeYs(Ba_isExcl, &Ra_anaYs);//, &Ra_anaCov);

	Ra_anaXt	= NULL;

	LOG("%d samples will be included in regression analysis\n", N_anaSamp);
}

cRegrAnalysis::~cRegrAnalysis()
{
	DEALLOC(Ba_isExcl);
//	sseFree(Ra_anaY);
	sseUnmat(Ra_anaYs, N_pheno);
// 	if (Ra_anaCov)
// 		deallocMatrix(Ra_anaCov, N_col-2, (void *)1);
}

void cRegrAnalysis::run()
{
	/* Do not perform analysis when there is no sample to analyze */
	if (N_anaSamp == 0) {
		LOGwarn("No sample was remained for regression analysis\n");
		return;
	}

	/* LASSO */
	if (OPT_ENABLED(lasso)) {
		if (!IS_ASSIGNED(lassolambda)) {
			OPTION().assign("lassolambda");
			OPTION().FORCE_OPT_REAL(lassolambda);
		}
	}

	wsUint		i;
	vVariant&	Xv_vrt		= Cp_IO->getVariant();

/**/cExporter*	Cp_regr		= cExporter::summon(B_isCont?"linear.regr.res":"logistic.regr.res");
	LOGoutput("Result of regression analysis is exported to [%s.%s]\n",
		OPT_STRING(out), B_isCont?"linear.regr.res":"logistic.regr.res");
	vSampPtr	Xa_samp		= Cp_IO->getSample();
	vCovar&		Xa_ci		= Cp_IO->getCovInfo();
	wsUintCst	N_oCov		= Cp_IO->sizeCovar();
	wsUintCst	N_gxeCov	= Cp_IO->sizeGxEcovar();
	wsUintCst	N_vrt		= Cp_IO->sizeVariant();

	/* Make X/Xt from phenotype-available samples */
	//	wsReal	**Ra_anaX	= NULL; /* checked */
	//	wsReal	**Ra_anaXt	= NULL; /* checked */
	_makeX(Ba_isExcl, &Ra_anaXt);

	if (OPT_ENABLED(donull))
		_doNull(Ra_anaXt);

	/* Annotation buffer */
	char S_anno[256] = {0, };
	if (IS_ASSIGNED(annogene))
		strcpy(S_anno, "ANNOT");

	/* Notify */
	LOG("%s regression start...\n", B_isCont?"Linear":"Logistic");

	/* Print header */
	char S_bufImp[16]	= { 0, };
	strcpy(S_bufImp, OPT_ENABLED(avail) ? "NMISSING" : "NIMPUTED");
	headerVariant(Cp_regr);
	char S_bufPheno[32] = { 0, };
	if (N_pheno > 1)
		strcpy(S_bufPheno, "PHENO	");
	Cp_regr->put(S_bufPheno);
	Cp_regr->fmt("	%s	%s	BETA	STAT	P_REGR", S_bufImp,
		B_isCont?"GINV":"NITER");

	/* Print covariates */
	for (i=0 ; i<N_oCov ; i++)
		Cp_regr->fmt("	BETA_%s	STAT_%s	P_%s", Xa_ci[i].Sp_varName,
			Xa_ci[i].Sp_varName, Xa_ci[i].Sp_varName);

	/* Print interactions (covGxE) */
	vInt &Xv_covIdx = getIO()->getGxEcovariatesIndices();
	if (OPT_ENABLED(gxe)) for (wsUint i=0 ; i<N_gxeCov ; i++) {
		xCovar &X_cov = Xa_ci[Xv_covIdx[i]];
		Cp_regr->fmt("	BETA_%s*GENO	STAT_%s*GENO	P_%s*GENO",
			X_cov.Sp_varName, X_cov.Sp_varName, X_cov.Sp_varName);
	}
	Cp_regr->put("\n");

	/* Return array size for each SNP calculate for thread>1 (covGxE) */
	wsUint		N_szRetPerSNP = NVAR_DEF_REGR + (
		OPT_ENABLED(gxe) ? ((N_oCov + N_gxeCov)*3) : N_oCov*3 );
	wsMat		Ra_rets		= NULL;
	bool		B_adjust	= OPT_ENABLED(adjust);
	vPheno&		Xv_pheno	= getIO()->getPhenoInfo();

	/* [[R]] Output buffer */
	wsMat	Ra_pvals	= sseMatrix(N_pheno, (wsUint)Xv_vrt.size());
	int**	Na_chrss	= (int **)sseUmatrix(N_pheno,(wsUint)Xv_vrt.size());
	wsUint*	Na_pval		= NULL;
	wsCalloc(Na_pval, wsUint, Xv_vrt.size());

	if (B_isCont)
		X_type = REGR_LINEAR;
	else
		X_type = REGR_LOGISTIC;

	/* Perform regression */
	if (Cp_anaSM && (OPT_ENABLED(lasso) || OPT_ENABLED(pls))) {
		LOG("Set-wise regression will be performed...\n");
		mGeneDef&	Xm_gene = Cp_anaSM->getGeneDef();
		cExporter*	Cp_la	= cExporter::summon("lasso.res");
		headerVariant(Cp_la);
		Cp_la->put("	SET	LAMBDA	BETA\n");
		wsUint		N_lasso	= 0;

		FOREACH (mGeneDef_it, Xm_gene, it) {
			vInt Xv_gidx = it->second;

			/* Do lasso */
			wsUint N_imp;
			wsMat Ra_curXt = NULL;
			wsMat Ra_curYs = NULL;

			/* Rebuild Xt matrix */
			wsUint N_curSamp = _setXSs(Xv_gidx, Ra_anaXt, Ra_anaYs,
				&Ra_curXt, &Ra_curYs, &N_imp);
			cStdMatrix	M_Xt((wsUint)Xv_gidx.size() + 1 + N_oCov, N_curSamp,
				Ra_curXt, MATDEL_NONE);

			if (OPT_ENABLED(lasso) || OPT_ENABLED(pls)) {
				char		S_bufPhe[256] = { 0, };
				for (wsUint i=0 ; i<N_pheno ; i++) {
					cVector		V_y(Ra_curYs[i], N_curSamp, 1);
					wsReal*		Ra_est = sseVector(M_Xt.row());
//					wsReal		R_maxDiff = W0;
					vInt		Xv_fidx;

					if (OPT_ENABLED(lasso)) {
// 						mmLasso(V_y, M_Xt, OPT_REAL(lassolambda), 1000, 1e-9,
//							Ra_est, &R_maxDiff);
						/*
						 * LassoFitGenerator is initialized
						 */
						LassoFitGenerator fitGenerator;
						fitGenerator.init(V_y, M_Xt);
                
						/*
						 * Generate the Lasso fit. The -1 arguments means that
						 * there would be no limit on the maximum number of 
						 * features per model
						 */
						LassoFit *fit = fitGenerator.fit(-1);
//						LOG(fit->toString().c_str());

						/* Fill Ra_est */
						wsUint i = 0;
						if (OPT_ENABLED(lassoall)) for (i=0 ; i < Xv_gidx.size() ; i++) {
							N_lasso++;
							Xv_fidx.push_back(Xv_gidx[i]);
							Ra_est[i] = fit->lambdas[i];
						} else for (i=0 ; fit->lambdas[i]>=OPT_REAL(lassolambda) && i < Xv_gidx.size() ; i++) {
							N_lasso++;
							Xv_fidx.push_back(Xv_gidx[i]);
							Ra_est[i] = fit->lambdas[i];
						}
// 						if (i > 0)
// 							LOG("[%d] reported\n", i);
					}
					if (OPT_ENABLED(pls)) {
						cgplsone_fit(GLM_GAUSSIAN, V_y, M_Xt);
					}
//						plsFitterwt(2, M_Xt, V_y, NULL, NULL, NULL, NULL, NULL, NULL);

					/* Print phenotpe info if necessary */
					if (N_pheno > 1)
						sprintf(S_bufPhe, "	%s", Xv_pheno[i].S_name.c_str());

					/* Print output */
					wsUint k=0;
					FOREACHDO (vInt_it, Xv_fidx, j, k++) {
						entryVariant(Cp_la, Xv_vrt[*j]);
						Cp_la->fmt("	%s%s	%g\n", it->first.c_str(),
							S_bufPhe, Ra_est[k]);
					}
				}
			} /* END OF IF : lasso & PLS */

			/* Deallocate allocated buffers */
			wsUint N_xCol = N_oCov + 1 + (wsUint)Xv_gidx.size();
			for (wsUint i=N_oCov+1 ; i<N_xCol ; i++)
				sseFree(Ra_curXt[i]);
			DEALLOC(Ra_curXt);

			/* Multiple regression */
		} /* END OF FOREACH : geneDef */

		LOG("%d/%d sets tested...\n", Xm_gene.size(), Xm_gene.size());
		if (OPT_ENABLED(lasso)) {
			if (N_lasso)
				LOGnote("[%d] LASSO results were reported with lambda cut-off [%g]\n", N_lasso, OPT_REAL(lassolambda));
			else
				LOGwarn("No results were reported due to too high lambda cut-off [%g]\n", OPT_REAL(lassolambda));
		}

		delete Cp_la;
	} else if (OPT_ENABLED(gxg)) {
		/* For each chromosome, do regression */
		for (wsUint i=1 ; i<=NCHR_SPECIES ; i++)
			_regressionMain(*Cp_regr, Ra_pvals, Na_chrss, Na_pval, i);
	} else {
		if (OPT_NUMBER(thread) == 1)
			_regressionMain(*Cp_regr, Ra_pvals, Na_chrss, Na_pval);
		else {
			/* Need to do this first */
			getIO()->getMAF();

			/* Prepare buffer */
			Ra_rets = sseMatrix(N_pheno, (N_szRetPerSNP)*N_vrt*2);

			/* Do regression */
			xAnaThread X_at = { getIO(), this };
			WORKER().run(thrRegression, forVariant_equal, &X_at, Ra_rets,
				sizeof(int)*3);

			/* Print out results */
			vDbl*	Xv_ress = new vDbl[N_pheno];

			wsUint k=0;
			LOOP (u, N_pheno) {
				wsReal*	Rp_pval	= Ra_pvals[u];
				int*	Np_chrs = Na_chrss[u];
				wsReal*	Ra_ret	= Ra_rets[u];

				vPheno	&Xa_phe	= getIO()->getPhenoInfo();
				const char	*S_p	= Xa_phe[u].S_name.c_str();

				for (i=0 ; i<Xv_vrt.size() ; i++) {
					wsReal*	Rp_ret	= Ra_ret + k;
					wsReal	R_pval	= Rp_ret[NVAR_DEF_REGR-1];

					/* --remna */
					if (NA(R_pval) && OPT_ENABLED(remna))
						continue;

					/* [[R]] collect pval and chr*/
					if (!NA(R_pval)) {
						Rp_pval[Na_pval[u]]		= R_pval;
						Np_chrs[Na_pval[u]++]	= Xv_vrt[i].chr;
					}

					/* --pvalrange */
					if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), R_pval))
						continue;

					/* Print common part */
					if (N_pheno > 1)
						sprintf(S_bufPheno, "%s	", S_p);
					entryVariant(Cp_regr, Xv_vrt[i]);
					Cp_regr->put(S_bufPheno);

					/* Print context */
					for (wsUint j=0 ; j<N_szRetPerSNP ; j++,k++) {
						if (NA(Ra_ret[k])) Cp_regr->put("	NA");
						else Cp_regr->fmt("	%g", Ra_ret[k]);
					}
					/* Multiple testing */
					if (B_adjust) Xv_ress[u].push_back(Rp_ret[NVAR_DEF_REGR-2]);
					Cp_regr->put("\n");
				}
			}
//			sseFree(Ra_ret);
			sseUnmat(Ra_rets, N_pheno);
			delete [] Xv_ress;

			/* Multiple testing adjustment in multi-threading env. */
			if (B_adjust) for (wsUint u=0 ; u<N_pheno ; u++) {
				vPheno &Xa_phe = getIO()->getPhenoInfo();
				/* Make name */
				char S_buf[512] = { 0, };
				if (N_pheno == 1)
					sprintf(S_buf, "%s.regr.adjust.res", X_type==REGR_LOGISTIC?"logistic":"linear");
				else
					sprintf(S_buf, "%s.regr.%s.adjust.res",
						X_type==REGR_LOGISTIC?"logistic":"linear",
						Xa_phe[u].S_name.c_str());
				multcomp(getIO(), Xv_ress[u], S_buf);
			}
		}

		LOG("%d/%d variants tested...\n", Xv_vrt.size(), Xv_vrt.size());
	} /* End of set/single */

	/* [[R]] Draw plot */
	char S_title[1024];
	sprintf(S_title, "%s regression", B_isCont?"Linear":"Logistic");
	qqVariant(S_title, B_isCont?"linear":"logistic", Ra_pvals, getIO());
	mhtVariant(S_title, B_isCont ? "linear" : "logistic", Ra_pvals, Na_chrss, getIO());
	sseUnmat(Ra_pvals, N_pheno);

	LOOP (i, N_pheno) sseFree(Na_chrss[i]);
	DEALLOC(Na_chrss);

	DEALLOC(Na_pval);

	delete Cp_regr;
	sseUnmat(Ra_anaXt, N_col);
}

int cRegrAnalysis::fitLogisticByNR(wsUint N_curSamp, wsUint N_col,
	wsReal **Ra_curXt, wsReal *Ra_curY, wsReal *Ra_coef, wsReal *Ra_V,
	wsUint N_maxIter, bool &B_conv, wsReal **Rp_phiHat/*=NULL*/,
	wsReal *Rp_sigma)
{
	///////////////////////////////////////
	// Newton-Raphson to fit logistic model   

	wsUint	i, j, k;
	B_conv = false;

	wsReal *p = NULL; /* checked */
	sseMalloc(p, wsReal, N_curSamp);

	wsReal **T = sseMatrix(N_col, N_col); /* checked */

	int		N_iter = 0;
	while (!B_conv && N_iter<(int)N_maxIter) {
		// Determine p and V

		wsReal **Ra_t = sseMpM(&Ra_coef, 1, N_col, Ra_curXt, N_col, N_curSamp);
		for (i=0 ; i<N_curSamp ; i++) {
			//wsReal t = sseVV(Ra_coef, N_col, Ra_curX[i]);
			p[i] = W1/(W1+exp(-Ra_t[0][i]));
			Ra_V[i] = p[i] * (W1 - p[i]);
		}
		sseUnmat(Ra_t, 1);

		// Update coefficients
		// b <- b +  solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 
		for (j=0 ; j<N_col ; j++)
			for (k=j ; k<N_col ; k++) {
				T[j][k] = T[k][j] =
					sseVVV(Ra_curXt[j], Ra_V, Ra_curXt[k], N_curSamp);
			}

		// solve( t(X) %*% V %*% X )

		wsReal **Ra_invT = NULL; /* checked */
		if (!EDinverse(T, N_col, &Ra_invT))
			goto _fin;

		// solve( t(X) %*% V %*% X ) %*% t(X)
		wsReal **T2 = sseSpM(Ra_invT, N_col, Ra_curXt, N_col, N_curSamp); /* checked */
		sseUnmat(Ra_invT, N_col);

		// (y-p)
		wsReal *Ra_yp = NULL; /* checked */
		sseMalloc(Ra_yp, wsReal, N_curSamp);
		sseVsV(Ra_curY, p, Ra_yp, N_curSamp);

		// solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 
		wsReal *Ra_newCoef = sseMpV(T2, N_col, N_curSamp, Ra_yp);
		sseFree(Ra_yp);
		sseUnmat(T2, N_col);

		// Update coefficients, and check for 
		// convergence
		wsReal R_delta = W0;
		for (j=0 ; j<N_col ; j++) {
			R_delta += fabs(Ra_newCoef[j]);
			Ra_coef[j] += Ra_newCoef[j];
		}
		sseFree(Ra_newCoef);

		//		printf("Delta %f\n", delta);
		if (NA(R_delta)) {
			B_conv = false;
			break;
		} if (R_delta < 1e-6) {
			B_conv = true;
			break;
		}

		// Next iteration
		N_iter++;
	}
_fin:
	sseUnmat(T, N_col);
	/* If sigma */
	if (Rp_sigma) {
		wsReal R_ss = W0;
		for (wsUint j=0 ; j<N_curSamp ; j++) {
			if (Ra_curY[j] == WISARD_AFFECTED)
				R_ss += REAL_CONST(-2.0) * log(p[j]);
			else
				R_ss += REAL_CONST(-2.0) * log(W1 - p[j]);
		}
		*Rp_sigma = R_ss / (wsReal)(N_curSamp - N_col);
	}

	if (Rp_phiHat)
		*Rp_phiHat = p;
	else
		sseFree(p);

	return N_iter;
}

/* Return structure
 * 1 6  5+0*3 + 5+0*3+2
 * 2 9
 * 3 12
 * [0]   [1]  [2]   [3]  [4]  [5]  [6]  [7]... [NVAR_DEF_REGR+(nCov-1)*3] . [NVAR_DEF_REGR+(nCov-1)*3+2]
 * N_imp beta tStat pVal gInv b_c1 t_c1 p_c1   b_cn                         p_cn
 */
int thrRegression(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	xAnaThread*		Xp_at		= (xAnaThread *)Vp_shareData;
	cRegrAnalysis*	Cp_reg		= (cRegrAnalysis *)(Xp_at->Vp_data);
	xRegrT			X_type		= Cp_reg->getType();
	cIO*			Cp_IO		= Xp_at->Cp_IO;
	vVariant&		Xa_SNP		= Cp_IO->getVariant();
	xMaf*			Xp_maf		= Cp_IO->getMAF();
	wsUintCst			N_oCov		= Cp_IO->sizeCovar();
	wsUintCst			N_gxeCov	= Cp_IO->sizeGxEcovar();
	wsUintCst			N_col		= Cp_reg->getBetaSize();
	int*			Np_data		= (int *)Vp_data;
	wsUint			N_s			= (wsUint)Np_data[0];
	wsUint			N_e			= (wsUint)Np_data[1];
	wsUint			N_pheno		= Cp_IO->sizePheno();
	wsMat			Ra_anaYs	= Cp_reg->getAnaYs();
	wsUintCst			N_anaSamp	= Cp_reg->getAnaSampSize();
	/* Set default X matrix */
	wsMat			Ra_anaXt	= NULL;
	Cp_reg->_getX(&Ra_anaXt);

	/* Return pointer */
	wsMat			Ra_rets		= (wsMat)Vp_result;

	/* Return size PER SNP (covGxE) */
	wsUint	N_szRetPerSNP = NVAR_DEF_REGR + (
		OPT_ENABLED(gxe) ? ((N_oCov+N_gxeCov)*3) : N_oCov*3 );

	xAnaRegr X = { X_type, 0, };
	X.N_col		= N_col;
	X.N_gxe		= N_gxeCov;
	X.N_pheno	= N_pheno;
	wsAlloc(X.Ra_b, wsReal*, N_pheno);
	X.Ra_t		= sseMatrix(N_pheno, N_col);
	X.Ra_swgt	= Cp_IO->getSampleWeight();

	//LOG("Thread %d : Range %d~%d\n", N_idx, N_s, N_e);
	for (wsUint j=N_s ; j<N_e ; j++) {
		xVariant&	X_curSNP	= Xa_SNP[j];
		wsUint		N_posIns	= N_szRetPerSNP*j;
		wsMat		Ra_curXt	= NULL;
		wsMat		Ra_curYs	= NULL;
		wsUint		N_imp		= 0;
		wsUint		N_curSamp;

		/* If MAF == 0 */
		if (Xp_maf[j].R_maf == W0) {
			for (wsUint u=0 ; u<N_pheno ; u++)
				for (wsUint i=0 ; i<N_szRetPerSNP ; i++)
					Ra_rets[u][N_posIns+i] = WISARD_NAN;
			continue;
		}
		/* Set genotype */
		N_curSamp = Cp_reg->_setXs(X_curSNP, j, Ra_anaXt, Ra_anaYs,
			&Ra_curXt, &Ra_curYs, &N_imp);

		/* If N_curSamp == 0 */
		if (N_curSamp == 0) {
			for (wsUint u=0 ; u<N_pheno ; u++)
				for (wsUint i=0 ; i<N_szRetPerSNP ; i++)
					Ra_rets[u][N_posIns+i] = WISARD_NAN;
			continue;
		}
		X.N_samp	= N_curSamp;
		X.Ra_Xt		= Ra_curXt;
		X.Ra_Ys		= Ra_curYs;
		X.X_type	= X_type;

		wsReal R_tDF = (wsReal)(N_curSamp - N_col - 1); // REGR_LINEAR
		wsReal*	Ra_V	= NULL; // REGR_LOGISTIC
		switch (X_type) {
		case REGR_LINEAR2:
		case REGR_LINEAR:
			X.B_gInv	= -1;
			break;
		case REGR_LOGISTIC:
			Ra_V		= sseVector(N_anaSamp);
		}
		sseMinit(X.Ra_t, N_pheno, N_col, WISARD_NAN);
		_regression(&X);

		/* Process results */
		for (wsUint u=0 ; u<N_pheno ; u++) {
			wsReal*	Ra_ret	= Ra_rets[u] + N_posIns;
			wsReal*	Ra_t	= X.Ra_t[u];
			wsReal*	Ra_b	= X.Ra_b[u];

			/* # of imputed genotypes */
			Ra_ret[0]		= N_imp;

			/* Perform Wald test to each beta */
			if (NA(Ra_t[N_col-1])) {
				for (wsUint i=1 ; i<=NVAR_DEF_REGR ; i++)
					Ra_ret[i] = WISARD_NAN;
			} else {
				wsReal R_pVal = WISARD_NAN;
				switch (X_type) {
				case REGR_LINEAR2:
				case REGR_LINEAR:
					R_pVal = ptdist(&Ra_t[N_col-1], &R_tDF);
					Ra_ret[1] = X.B_gInv == -1 ? WISARD_NAN : (wsReal)X.B_gInv;
					break;
				case REGR_LOGISTIC:
					R_pVal = PVchisq(SQR(Ra_t[N_col-1]), 1.0);
					Ra_ret[1] = (wsReal)X.N_iter;
					break;
				}
				/* NVAR_DEF_REGR-1 variables should be filled in here */
				Ra_ret[2] = Ra_b[N_col-1];
				Ra_ret[3] = Ra_t[N_col-1];
				Ra_ret[4] = R_pVal;
			}

			/* Print covariate test */
			wsUint N_from = NVAR_DEF_REGR;

			for (wsUint i=1 ; i<=N_oCov ; i++) {
				if (NA(Ra_t[i])) {
					Ra_ret[N_from+3*i-3]	= WISARD_NAN;
					Ra_ret[N_from+3*i-2]	= WISARD_NAN;
					Ra_ret[N_from+3*i-1]	= WISARD_NAN;
				} else {
					wsReal R_pValCov = WISARD_NA;
					switch (X_type) {
					case REGR_LINEAR2:
					case REGR_LINEAR:
						R_pValCov = ptdist(Ra_t+i, &R_tDF); break;
					case REGR_LOGISTIC:
						R_pValCov = PVchisq(SQR(Ra_t[i]), 1.0); break;
					}
					R_pValCov = ptdist(Ra_t+i, &R_tDF);

					Ra_ret[N_from+3*i-3]	= Ra_b[i];
					Ra_ret[N_from+3*i-2]	= Ra_t[i];
					Ra_ret[N_from+3*i-1]	= R_pValCov;
				}
			}
			/* Print interaction test (covGxE) */
			if (OPT_ENABLED(gxe)) for (wsUint i=1 ; i<=N_gxeCov ; i++) {
				if (NA(Ra_t[i+N_oCov])) {
					Ra_ret[N_from+N_oCov*3 + 3*i-3]	= WISARD_NAN;
					Ra_ret[N_from+N_oCov*3 + 3*i-2]	= WISARD_NAN;
					Ra_ret[N_from+N_oCov*3 + 3*i-1]	= WISARD_NAN;
				} else {
					wsReal R_pValCov = WISARD_NA;
					switch (X_type) {
					case REGR_LINEAR2:
					case REGR_LINEAR:
						R_pValCov = ptdist(&Ra_t[i+N_oCov], &R_tDF); break;
					case REGR_LOGISTIC:
						R_pValCov = PVchisq(SQR(Ra_t[i+N_oCov]), 1.0); break;
					}
					R_pValCov = ptdist(&Ra_t[i+N_oCov], &R_tDF);

					Ra_ret[N_from+N_oCov*3 + 3*i-3]	= Ra_b[i+1];
					Ra_ret[N_from+N_oCov*3 + 3*i-2]	= Ra_t[i+N_oCov];
					Ra_ret[N_from+N_oCov*3 + 3*i-1]	= R_pValCov;
				}
			}
			sseFree(Ra_b);
			X.Ra_b[u] = NULL;
		}

		if (Ra_curXt != Ra_anaXt) {
			sseUnmat(Ra_curXt, N_col);
			sseUnmat(Ra_curYs, N_pheno);
		}

		/* Now cleanup */
		switch (X_type) {
		case REGR_LOGISTIC:
			sseFree(Ra_V); break;
		default:
			break;
		}

		/* Count done */
		Np_data[2]++;

		/* Print progress */
		if (N_idx==0 && (j%100)==0) {
			wsUint N_sum = 0;
			for (int x=2,y=0 ; y<OPT_NUMBER(thread) ; x+=3,y++)
				N_sum += Np_data[x];
			notice("[%d/%d] variants tested...\r", N_sum, Xa_SNP.size());
		}
	}
	sseUnmat(Ra_anaXt, N_col);
	sseUnmat(X.Ra_t, N_pheno);
	DEALLOC(X.Ra_b);

	return 0;
}

wsUint cRegrAnalysis::_makeY(const char *Ba_filter, wsReal **Ra_resY)//,
//	wsReal ***Ra_resCov)
{
	wsUint	i, j, N_resSamp;
	wsUint	N_origSamp	= Cp_IO->sizeSample();
	wsRealCst	*Ra_origY	= Cp_IO->getPheno();
// 	wsReal	**Ra_origC	= Cp_IO->getCovariates();
// 	wsUint	N_cov		= Cp_IO->getCovariateSize();

	/* Count the size of result sample */
	for (N_resSamp=i=0 ; i<N_origSamp ; i++)
		if (Ba_filter[i] == 0) N_resSamp++;

	sseMalloc(*Ra_resY, wsReal, N_resSamp);

	/* Resize phenotype and covariates */
	if (Cp_IO->isContinuous()) {
		for (i=j=0 ; i<N_origSamp ; i++) {
			if (Ba_filter[i]) continue;

			(*Ra_resY)[j++] = Ra_origY[i];
		}
	} else {
		for (i=j=0 ; i<N_origSamp ; i++) {
			if (Ba_filter[i]) continue;

			(*Ra_resY)[j++] = Ra_origY[i];
		}
	}
// 	if (Ra_origC) {
// 		*Ra_resCov	= sseMatrix(N_cov, N_anaSamp);
// 		for (i=0 ; i<N_cov ; i++) for (j=k=0 ; j<N_anaSamp ; j++) {
// 			if (Ba_isExcl[j]) continue;
// 			(*Ra_resCov)[i][k++] = Ra_origC[i][j];
// 		}
// 	} else
// 		*Ra_resCov = NULL;

	return N_resSamp;
}

wsUint cRegrAnalysis::_makeYs(const char *Ba_filter, wsMat *Ra_resY)
{
	wsUint	i, J, j, N_resSamp;
	wsUint	N_origSamp	= Cp_IO->sizeSample();
	wsUint	N_pheno		= Cp_IO->sizePheno();
	wsMat	Ra_origY	= Cp_IO->getPhenos();

	/* Count the size of result sample */
	for (N_resSamp=i=0 ; i<N_origSamp ; i++)
		if (Ba_filter[i] == 0) N_resSamp++;
	*Ra_resY = sseMatrix(N_pheno, N_resSamp);

	/* Resize phenotype and covariates */
	for (i=0 ; i<N_pheno ; i++)
		for (j=J=0 ; j<N_origSamp ; j++) {
			if (Ba_filter[j]) continue;
			(*Ra_resY)[i][J++] = Ra_origY[i][j];
		}

	return N_resSamp;
}

wsUint cRegrAnalysis::_makeX(const char *Ba_filter, wsMat *Ra_resXt)
{
	wsUint	i, j, k;
	wsUint	N_origSamp	= Cp_IO->sizeSample();
	wsUint	N_origCov	= Cp_IO->sizeCovar();
	wsReal	**Ra_cov	= Cp_IO->getCovariates();
	*Ra_resXt			= sseMatrix(N_col, N_anaSamp);

	/* First column to be intercept
	 * 2~(n-1) column to be covariates
	 * n column to be X */
	for (i=0 ; i<N_anaSamp ; i++) {
		(*Ra_resXt)[0][i]	= W1;
	}
	for (i=1 ; i<=N_origCov ; i++) {
		for (j=k=0 ; j<N_origSamp ; j++) {
			if (Ba_isExcl[j]) continue;
			(*Ra_resXt)[i][k] = Ra_cov[i-1][j];
			k++;
		}
	}

	/* Currently meaningless */
	return 0;
}

void cRegrAnalysis::_getX(wsReal ***Ra_retAnaXt)
{
	if (!Ra_anaXt)
		halt_fmt(WISARD_SYST_NULL_INPUT, "Genotype matrix");
//		halt("SYSERR: getX() failed");

	*Ra_retAnaXt	= sseMatrixP(N_col, N_anaSamp, Ra_anaXt);
}

wsUintCst cRegrAnalysis::_setX(xVariant& X_snp, wsUint N_idxSNP,
	wsMat Ra_anaXt, wsReal *Ra_anaY, wsReal ***Rp_curXt, wsReal **Rp_curY,
	wsUint *Np_imp)
{
	wsUint		j, k;
	wsUint		N_curSamp;
	wsUint		N_origSamp	= Cp_IO->sizeSample();

	wsUintCst		N_origCov	= Cp_IO->sizeCovar();
	wsUintCst		N_gxeCov	= Cp_IO->sizeGxEcovar();
	//MAT_t		Ra_cov		= Cp_IO->getCovariates();
	wsMat		Ra_covGxE	= Cp_IO->getGxEcovariates();

	char		**Na_data	= Cp_IO->getGenotype();
	vSampPtr&	Xa_samp		= Cp_IO->getSample();
	vVariant&	Xa_snp		= Cp_IO->getVariant();

	(*Np_imp) = 0;
	if (OPT_ENABLED(avail)) {
		char *Ba_curExcl = NULL;
		wsAlloc(Ba_curExcl, char, N_origSamp);
		memcpy(Ba_curExcl, Ba_isExcl, sizeof(char)*N_origSamp);

		/* Get the number */
		N_curSamp = 0;
		for (j=0 ; j<N_origSamp ; j++) {
			if (Ba_isExcl[j]) continue;
			if (isMissing(Na_data[j][N_idxSNP])) {
				(*Np_imp)++;
				Ba_curExcl[j] = 1;
			} else
				N_curSamp++;
		}

		if (N_curSamp == N_anaSamp) {
			/* If there is no missing genotype, use as is */
			N_curSamp	= N_anaSamp;
			*Rp_curXt	= Ra_anaXt;
			*Rp_curY	= Ra_anaY;

			/* But we need to add cov*geno term when --gxe (covGxE) */
			if (OPT_ENABLED(gxe)) for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_curExcl[j]) continue;
				char N_geno = Na_data[j][N_idxSNP];

				if (isMissing(N_geno))
					halt_fmt(WISARD_SYST_CANT_MISGENO, Xa_samp[j]->S_IID.c_str(),
						Xa_snp[N_idxSNP].name);
//					halt("SYSERR: Genotype should be exists!");
				for (wsUint i=0 ; i<N_gxeCov ; i++)
						(*Rp_curXt)[N_origCov+1+i][k] =
							Na_data[j][N_idxSNP]*Ra_covGxE[i][j];
				(*Rp_curXt)[N_col-1][k] =
					Na_data[j][N_idxSNP];
				k++;
			} else for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_curExcl[j]) continue;
				char N_geno = Na_data[j][N_idxSNP];

				if (isMissing(N_geno))
					halt_fmt(WISARD_SYST_CANT_MISGENO, Xa_samp[j]->S_IID.c_str(),
						Xa_snp[N_idxSNP].name);
				(*Rp_curXt)[N_col-1][k] =
					Na_data[j][N_idxSNP];
				k++;
			}
		} else {
			/* Rebuild X, Xt, Y */
			_makeY(Ba_curExcl, Rp_curY);//, &Ra_curCov);
			//N_curSamp =
			_makeX(Ba_curExcl, Rp_curXt);

			/* If --gxe, add interaction (covGxE) */
			if (OPT_ENABLED(gxe)) {
				/* Fill out genotype part */
				for (j=k=0 ; j<N_origSamp ; j++) {
					if (Ba_curExcl[j]) continue;
					char N_geno = Na_data[j][N_idxSNP];

					if (isMissing(N_geno))
						halt_fmt(WISARD_SYST_NULL_REGRGENO, Xa_snp[N_idxSNP].name,
							Xa_samp[j]->S_FID.c_str(), Xa_samp[j]->S_IID.c_str());
					for (wsUint i=0 ; i<N_gxeCov ; i++)
							(*Rp_curXt)[N_origCov+1+i][k] =
								N_geno*Ra_covGxE[i][j];
					(*Rp_curXt)[N_col-1][k++] = N_geno;
				}
			} else {
				/* Fill out genotype part */
				for (j=k=0 ; j<N_origSamp ; j++) {
					if (Ba_curExcl[j]) continue;
					char N_geno = Na_data[j][N_idxSNP];

					if (isMissing(N_geno))
						halt_fmt(WISARD_SYST_NULL_REGRGENO, Xa_snp[N_idxSNP].name,
							Xa_samp[j]->S_FID.c_str(), Xa_samp[j]->S_IID.c_str());
					(*Rp_curXt)[N_col-1][k++] = N_geno;
				}
			}
		}

		DEALLOC(Ba_curExcl);
	} else {
		vSampPtr&	Xa_samp	= Cp_IO->getSample();
		xMaf&		X_maf	= Cp_IO->getMAF()[N_idxSNP];

		/* Build X (covGxE) */
		if (OPT_ENABLED(gxe)) {
			for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_isExcl[j]) continue;
				char N_geno = Na_data[j][N_idxSNP];
				wsReal R_geno;

				if (isMissing(N_geno)) {
					R_geno = _imputeGeno(X_snp, Xa_samp[j], X_maf.R_maf);
					(*Np_imp)++;
				} else
					R_geno = (wsReal)Na_data[j][N_idxSNP];

				for (wsUint i=0 ; i<N_gxeCov ; i++)
						Ra_anaXt[N_origCov+1+i][k] =
							R_geno*Ra_covGxE[i][j];
				Ra_anaXt[N_col-1][k] = R_geno;

				k++;
			}
		} else {
			wsReal R_Q = W0;
			for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_isExcl[j]) continue;
				char N_geno = Na_data[j][N_idxSNP];
				wsReal R_geno;

				if (isMissing(N_geno)) {
					R_geno = _imputeGeno(X_snp, Xa_samp[j], X_maf.R_maf);
					(*Np_imp)++;
				} else
					R_geno = (wsReal)Na_data[j][N_idxSNP];

				Ra_anaXt[N_col-1][k] = R_geno;
				R_Q += R_geno;

				k++;
			}
			/* DO NOT ALLOW SUM == 0 */
			if (R_Q == W0)
				return 0;
		}
		N_curSamp	= N_anaSamp;
		*Rp_curXt	= Ra_anaXt;
		*Rp_curY	= Ra_anaY;
	}

	return N_curSamp;
}

wsUintCst cRegrAnalysis::_setXs(xVariant& X_snp, wsUint N_idxSNP,
	wsMat Ra_anaXt, wsMat Ra_anaYs, wsMat *Rp_curXt, wsMat *Rp_curYs,
	wsUint *Np_imp, wsUint N_idxGxG/*=0xffffffff*/)
{
	wsUint		j, k;
	wsUint		N_curSamp;
	wsUint		N_origSamp	= Cp_IO->sizeSample();

	wsUintCst		N_origCov	= Cp_IO->sizeCovar();
	/* --gxg then ignore --gxe, so N_gxeCov == 2, for G2 and G1xG2 */
	wsUintCst		N_gxeCov	= N_idxGxG != 0xffffffff ?
		2 : Cp_IO->sizeGxEcovar();
	//MAT_t		Ra_cov		= Cp_IO->getCovariates();
	wsMat		Ra_covGxE	= Cp_IO->getGxEcovariates();

	char**		Na_data		= Cp_IO->getGenotype();
	vSampPtr&	Xa_samp		= Cp_IO->getSample();
	vVariant&	Xa_snp		= Cp_IO->getVariant();

	(*Np_imp) = 0;
	if (OPT_ENABLED(avail)) {
		char *Ba_curExcl = NULL;
		wsAlloc(Ba_curExcl, char, N_origSamp);
		memcpy(Ba_curExcl, Ba_isExcl, sizeof(char)*N_origSamp);

		/* Get the number */
		N_curSamp = 0;
		for (j=0 ; j<N_origSamp ; j++) {
			if (Ba_isExcl[j]) continue;
			if (isMissing(Na_data[j][N_idxSNP])) {
				(*Np_imp)++;
				Ba_curExcl[j] = 1;
			} else
				N_curSamp++;
		}

		if (N_curSamp == N_anaSamp) {
			/* If there is no missing genotype, use as is */
			N_curSamp	= N_anaSamp;
			*Rp_curXt	= Ra_anaXt;
			*Rp_curYs	= Ra_anaYs;

			/* But we need to add cov*geno term when --gxe (covGxE) */
			if (OPT_ENABLED(gxg))
				setGxGeno(Cp_IO, N_idxSNP, N_idxGxG, (*Rp_curXt)[N_col-3],
					(*Rp_curXt)[N_col-2], (*Rp_curXt)[N_col-1],
					Ba_curExcl, 1);
			else if (OPT_ENABLED(gxe)) for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_curExcl[j]) continue;
				char N_geno = Na_data[j][N_idxSNP];

				if (isMissing(N_geno))
					halt_fmt(WISARD_SYST_CANT_MISGENO, Xa_samp[j]->S_IID.c_str(),
					Xa_snp[N_idxSNP].name);
				//					halt("SYSERR: Genotype should be exists!");
				for (wsUint i=0 ; i<N_gxeCov ; i++)
					(*Rp_curXt)[N_origCov+1+i][k] =
						N_geno*Ra_covGxE[i][j];
				(*Rp_curXt)[N_col-1][k] = N_geno;
				k++;
			} else
				setGeno(Cp_IO, N_idxSNP, (*Rp_curXt)[N_col-1], Ba_curExcl, 1);
		} else {
			/* Rebuild X, Xt, Y */
			_makeYs(Ba_curExcl, Rp_curYs);//, &Ra_curCov);
			//N_curSamp =
			_makeX(Ba_curExcl, Rp_curXt);

			/* If --gxe, add interaction (covGxE) */
			if (OPT_ENABLED(gxg))
				setGxGeno(Cp_IO, N_idxSNP, N_idxGxG, (*Rp_curXt)[N_col-3],
					(*Rp_curXt)[N_col-2], (*Rp_curXt)[N_col-1],
					Ba_curExcl, 1);
			else if (OPT_ENABLED(gxe)) {
				/* Fill out genotype part */
				for (j=k=0 ; j<N_origSamp ; j++) {
					if (Ba_curExcl[j]) continue;
					char N_geno = Na_data[j][N_idxSNP];

					if (isMissing(N_geno))
						halt_fmt(WISARD_SYST_NULL_REGRGENO, Xa_snp[N_idxSNP].name,
						Xa_samp[j]->S_FID.c_str(), Xa_samp[j]->S_IID.c_str());
					for (wsUint i=0 ; i<N_gxeCov ; i++)
						(*Rp_curXt)[N_origCov+1+i][k] =
							N_geno*Ra_covGxE[i][j];
					(*Rp_curXt)[N_col-1][k++] = N_geno;
				}
			} else {
				/* Fill out genotype part */
				setGeno(Cp_IO, N_idxSNP, (*Rp_curXt)[N_col-1], Ba_curExcl, 1);
// 				for (j=k=0 ; j<N_origSamp ; j++) {
// 					if (Ba_curExcl[j]) continue;
// 					char N_geno = Na_data[j][N_idxSNP];
// 
// 					if (isMissing(N_geno))
// 						halt_fmt(WISARD_SYST_NULL_REGRGENO, Xa_snp[N_idxSNP].name,
// 						Xa_samp[j]->S_FID.c_str(), Xa_samp[j]->S_IID.c_str());
// 					(*Rp_curXt)[N_col-1][k++] = N_geno;
// 				}
			}
		}

		DEALLOC(Ba_curExcl);

		/* End of --avail case */
	} else {
		vSampPtr&	Xa_samp = Cp_IO->getSample();
		xMaf&		X_maf	= Cp_IO->getMAF()[N_idxSNP];

		/* Build X (covGxE) */
		if (OPT_ENABLED(gxg)) {
			wsReal R_Q1 = W0;
			wsReal R_Q2 = W0;
			wsReal R_Q3 = W0;
			for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_isExcl[j]) continue;
				char N_g1 = Na_data[j][N_idxSNP];
				char N_g2 = Na_data[j][N_idxGxG];
//				char N_g3 = N_g1 * N_g2;
				wsReal R_g1 = WISARD_NAN, R_g2 = WISARD_NAN, R_g3 = WISARD_NAN;

				if (isMissing(N_g1)) {
					R_g1 = _imputeGeno(X_snp, Xa_samp[j], X_maf.R_maf);
					(*Np_imp)++;
				} else R_g1 = (wsReal)N_g1;
				if (isMissing(N_g2)) {
					R_g2 = _imputeGeno(X_snp, Xa_samp[j], X_maf.R_maf);
					(*Np_imp)++;
				} else R_g2 = (wsReal)N_g2;

				R_g3 = R_g1 * R_g2;

				Ra_anaXt[N_col-3][k] = R_g1;
				Ra_anaXt[N_col-2][k] = R_g2;
				Ra_anaXt[N_col-1][k] = R_g3;
				R_Q1 += R_g1;
				R_Q2 += R_g2;
				R_Q3 += R_g3;

				k++;
			}
			/* DO NOT ALLOW SUM == 0 */
			if (R_Q1 == W0 || R_Q2 == W0 || R_Q3 == W0)
				return 0;
		} else if (OPT_ENABLED(gxe)) {
			for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_isExcl[j]) continue;
				char N_geno = Na_data[j][N_idxSNP];
				wsReal R_geno;

				if (isMissing(N_geno)) {
					R_geno = _imputeGeno(X_snp, Xa_samp[j], X_maf.R_maf);
					(*Np_imp)++;
				} else R_geno = (wsReal)N_geno;

				for (wsUint i=0 ; i<N_gxeCov ; i++)
					Ra_anaXt[N_origCov+1+i][k] =
						R_geno*Ra_covGxE[i][j];
				Ra_anaXt[N_col-1][k] = R_geno;

				k++;
			}
		} else {
			if (IS_ASSIGNED(dosage)) {
				setGeno(Cp_IO, N_idxSNP, Ra_anaXt[N_col-1], Ba_isExcl, 1);
			} else {
			wsReal R_Q = W0;
			for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_isExcl[j]) continue;
				char N_geno = Na_data[j][N_idxSNP];
				wsReal R_geno;

				if (isMissing(N_geno)) {
					R_geno = _imputeGeno(X_snp, Xa_samp[j], X_maf.R_maf);
					(*Np_imp)++;
				} else R_geno = (wsReal)N_geno;

				Ra_anaXt[N_col-1][k] = R_geno;
				R_Q += R_geno;

				k++;
			}
			/* DO NOT ALLOW SUM == 0 */
			if (R_Q == W0)
				return 0;
			}
		}
		N_curSamp	= N_anaSamp;
		*Rp_curXt	= Ra_anaXt;
		*Rp_curYs	= Ra_anaYs;
	}

	return N_curSamp;
}

wsUintCst cRegrAnalysis::_setXSs(vInt& Xv_idx, wsMat Ra_anaXt, wsMat Ra_anaYs,
	wsMat *Rp_curXt, wsMat *Rp_curYs, wsUint *Np_imp)
{
	wsUint		j, k;
	wsUint		N_curSamp;
	wsUint		N_origSamp	= Cp_IO->sizeSample();

	wsUintCst		N_origCov	= Cp_IO->sizeCovar();
	wsUintCst		N_col		= (wsUint)Xv_idx.size() + N_origCov + 1;
//	MAT_t		Ra_cov		= Cp_IO->getCovariates();

	char**		Na_data		= Cp_IO->getGenotype();
	vSampPtr&	Xa_samp		= Cp_IO->getSample();
	vVariant&	Xa_snp		= Cp_IO->getVariant();
	xMaf*		Xp_maf		= Cp_IO->getMAF();

	/* Allocate Rp_curXt */
	wsAlloc(*Rp_curXt, wsVec, N_col);
	for (j=N_origCov+1 ; j<N_col ; j++)
		(*Rp_curXt)[j] = sseVector(N_anaSamp);
	/* Set common part */
	for (j=0 ; j<=N_origCov ; j++)
		(*Rp_curXt)[j] = Ra_anaXt[j];

	(*Np_imp) = 0;

	wsUint I = 1;
	FOREACHDO (vInt_it, Xv_idx, i, I++) {
		wsReal		R_Q			= W0;
		wsUint		N_idxSNP	= *i;
		xVariant&	X_snp		= Xa_snp[N_idxSNP];
		xMaf&		X_maf		= Xp_maf[N_idxSNP];
		for (j=k=0 ; j<N_origSamp ; j++) {
			if (Ba_isExcl[j]) continue;
			char N_geno = Na_data[j][N_idxSNP];
			wsReal R_geno;

			if (isMissing(N_geno)) {
				R_geno = _imputeGeno(X_snp, Xa_samp[j], X_maf.R_maf);
				(*Np_imp)++;
			} else R_geno = (wsReal)N_geno;

			(*Rp_curXt)[N_origCov+I][k] = R_geno;
			R_Q += R_geno;

			k++;
		}
		/* DO NOT ALLOW SUM == 0 */
		if (R_Q == W0)
			return 0;
	}
	N_curSamp	= N_anaSamp;
	*Rp_curYs	= Ra_anaYs;

	return N_curSamp;
}

wsUintCst cRegrAnalysis::_setXXs(xVariant& X_s1, xVariant &X_s2,
	wsUint N_x1, wsUint N_x2, wsMat Ra_anaXt, wsMat Ra_anaYs,
	wsMat *Rp_curXt, wsMat *Rp_curYs, wsUint *Np_imp)
{
	wsUint		j, k;
	wsUint		N_curSamp;
	wsUint		N_origSamp	= Cp_IO->sizeSample();

//	wsUintCst		N_origCov	= Cp_IO->getCovariateSize();
//	wsUintCst		N_gxeCov	= Cp_IO->getGxEcovariateSize();
	//MAT_t		Ra_cov		= Cp_IO->getCovariates();
//	MAT_t		Ra_covGxE	= Cp_IO->getGxEcovariates();

	char**		Na_data		= Cp_IO->getGenotype();
	vSampPtr&	Xa_samp		= Cp_IO->getSample();
	vVariant&	Xa_snp		= Cp_IO->getVariant();

	(*Np_imp) = 0;
	if (OPT_ENABLED(avail)) {
		char *Ba_curExcl = NULL;
		wsAlloc(Ba_curExcl, char, N_origSamp);
		memcpy(Ba_curExcl, Ba_isExcl, sizeof(char)*N_origSamp);

		/* Get the number */
		N_curSamp = 0;
		for (j=0 ; j<N_origSamp ; j++) {
			if (Ba_isExcl[j]) continue;
			if (isMissing(Na_data[j][N_x1]) || isMissing(Na_data[j][N_x2])) {
				(*Np_imp)++;
				Ba_curExcl[j] = 1;
			} else
				N_curSamp++;
		}

		if (N_curSamp == N_anaSamp) {
			/* If there is no missing genotype, use as is */
			N_curSamp	= N_anaSamp;
			*Rp_curXt	= Ra_anaXt;
			*Rp_curYs	= Ra_anaYs;

			/* But we need to add cov*geno term when --gxe (covGxE) */
			if (OPT_ENABLED(gxe)) for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_curExcl[j]) continue;
				char N_g1 = Na_data[j][N_x1];
				char N_g2 = Na_data[j][N_x2];

				if (isMissing(N_g1) || isMissing(N_g2))
					halt_fmt(WISARD_SYST_CANT_MISGENO, Xa_samp[j]->S_IID.c_str(),
						Xa_snp[N_x1].name);
				//					halt("SYSERR: Genotype should be exists!");
// 				for (wsUint i=0 ; i<N_gxeCov ; i++)
// 					(*Rp_curXt)[N_origCov+1+i][k] =
// 						N_geno*Ra_covGxE[i][j];
				(*Rp_curXt)[N_col-3][k] = N_g1;
				(*Rp_curXt)[N_col-2][k] = N_g2;
				(*Rp_curXt)[N_col-1][k] = N_g1*N_g2;
				k++;
			} else // for (j=k=0 ; j<N_origSamp ; j++) {
				// 				if (Ba_curExcl[j]) continue;
				// 				char N_geno = Na_data[j][N_idxSNP];
				// 
				// 				if (isMissing(N_geno))
				// 					halt_fmt(WISARD_SYST_CANT_MISGENO, Xa_samp[j]->S_IID.c_str(),
				// 					Xa_snp[N_idxSNP].name);
				// 				(*Rp_curXt)[N_col-1][k] =
				// 					Na_data[j][N_idxSNP];
				// 				k++;
				// 			}
				setGxGeno(Cp_IO, N_x1, N_x2, (*Rp_curXt)[N_col-3],
					(*Rp_curXt)[N_col-2], (*Rp_curXt)[N_col-1],
					Ba_curExcl, 1);
		} else {
			/* Rebuild X, Xt, Y */
			_makeYs(Ba_curExcl, Rp_curYs);//, &Ra_curCov);
			//N_curSamp =
			_makeX(Ba_curExcl, Rp_curXt);

			/* If --gxe, add interaction (covGxE) */
			if (OPT_ENABLED(gxe)) {
				/* Fill out genotype part */
				for (j=k=0 ; j<N_origSamp ; j++) {
					if (Ba_curExcl[j]) continue;
					char N_g1 = Na_data[j][N_x1];
					char N_g2 = Na_data[j][N_x2];

					if (isMissing(N_g1) || isMissing(N_g2))
						halt_fmt(WISARD_SYST_NULL_REGRGENO, Xa_snp[N_x1].name,
							Xa_samp[j]->S_FID.c_str(), Xa_samp[j]->S_IID.c_str());
// 					for (wsUint i=0 ; i<N_gxeCov ; i++)
// 						(*Rp_curXt)[N_origCov+1+i][k] =
// 							N_geno*Ra_covGxE[i][j];
					(*Rp_curXt)[N_col-3][k] = N_g1;
					(*Rp_curXt)[N_col-2][k] = N_g2;
					(*Rp_curXt)[N_col-1][k] = N_g1*N_g2;
					k++;
				}
			} else {
				/* Fill out genotype part */
				setGxGeno(Cp_IO, N_x1, N_x2, (*Rp_curXt)[N_col-3],
					(*Rp_curXt)[N_col-2], (*Rp_curXt)[N_col-1],
					Ba_curExcl, 1);
// 				for (j=k=0 ; j<N_origSamp ; j++) {
// 					if (Ba_curExcl[j]) continue;
// 					char N_geno = Na_data[j][N_idxSNP];
// 
// 					if (isMissing(N_geno))
// 						halt_fmt(WISARD_SYST_NULL_REGRGENO, Xa_snp[N_idxSNP].name,
// 						Xa_samp[j]->S_FID.c_str(), Xa_samp[j]->S_IID.c_str());
// 					(*Rp_curXt)[N_col-1][k++] = N_geno;
// 				}
			}
		}

		DEALLOC(Ba_curExcl);
	} else {
		vSampPtr&	Xa_samp	= Cp_IO->getSample();
		xMaf&		X_maf1	= Cp_IO->getMAF()[N_x1];
		xMaf&		X_maf2	= Cp_IO->getMAF()[N_x2];

		/* Build X (covGxE) */
		if (OPT_ENABLED(gxe)) {
			for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_isExcl[j]) continue;
				char N_g1 = Na_data[j][N_x1];
				char N_g2 = Na_data[j][N_x2];
				wsReal R_g1, R_g2;

				if (isMissing(N_g1) || isMissing(N_g2)) {
					R_g1 = _imputeGeno(X_s1, Xa_samp[j], X_maf1.R_maf);
					R_g2 = _imputeGeno(X_s2, Xa_samp[j], X_maf2.R_maf);
					(*Np_imp)++;
				} else {
					R_g1 = (wsReal)N_x1;
					R_g2 = (wsReal)N_x2;
				}

// 				for (wsUint i=0 ; i<N_gxeCov ; i++)
// 					Ra_anaXt[N_origCov+1+i][k] =
// 						R_geno*Ra_covGxE[i][j];
				Ra_anaXt[N_col-3][k] = R_g1;
				Ra_anaXt[N_col-2][k] = R_g2;
				Ra_anaXt[N_col-1][k] = R_g1*R_g2;
				k++;
			}
		} else {
			wsReal R_Q = W0;
			for (j=k=0 ; j<N_origSamp ; j++) {
				if (Ba_isExcl[j]) continue;
				char N_g1 = Na_data[j][N_x1];
				char N_g2 = Na_data[j][N_x2];
				wsReal R_g1, R_g2;

				if (isMissing(N_g1)) {
					R_g1 = _imputeGeno(X_s1, Xa_samp[j], X_maf1.R_maf);
					(*Np_imp)++;
				} else R_g1 = (wsReal)N_g1;

				if (isMissing(N_g2)) {
					R_g2 = _imputeGeno(X_s2, Xa_samp[j], X_maf2.R_maf);
					(*Np_imp)++;
				} else R_g2 = (wsReal)N_g2;

				Ra_anaXt[N_col-3][k] = R_g1;
				Ra_anaXt[N_col-2][k] = R_g2;
				Ra_anaXt[N_col-1][k] = R_g1 * R_g2;
				R_Q += R_g1 + R_g2;

				k++;
			}
			/* DO NOT ALLOW SUM == 0 */
			if (R_Q == W0)
				return 0;
		}
		N_curSamp	= N_anaSamp;
		*Rp_curXt	= Ra_anaXt;
		*Rp_curYs	= Ra_anaYs;
	}

	return N_curSamp;
}

/*
 * Main regression function
 */
void _regression(xAnaRegr *Xp_a)
{
	wsUint		N_gxe	= Xp_a->N_gxe;
	wsUint		N_pheno	= Xp_a->N_pheno;
	wsUint		N_samp	= Xp_a->N_samp;
	wsUint		N_col	= Xp_a->N_col;
	wsMat		Ra_Ys	= Xp_a->Ra_Ys;
	wsUint		N_cov	= N_col - 2 - N_gxe;
	if (N_samp == 0) return;

	wsReal R_tDF	= (wsReal)(N_samp - N_col - 1);

	cStdMatrix	M_Xt(N_col, N_samp, Xp_a->Ra_Xt, MATDEL_NONE);
	wsSym		Ra_XXi	= NULL;	// REGR_LINEAR
	wsSym		Ra_Si	= NULL;	// REGR_LOGISTIC
	/* if linreg do XXi */
#ifdef USE_LINEAR2
	if (Xp_a->X_type == REGR_LINEAR) {
#else
	if (Xp_a->X_type == REGR_LINEAR || Xp_a->X_type == REGR_LINEAR2) {
#endif
		cSymMatrix	M_XX	= M_Xt.Mt();
		cSymMatrix&	M_XXi	= M_XX.inv(NULL, &(Xp_a->B_gInv));
		if (!M_XXi.sane()) return;
		M_XXi.setDontDealloc();

		Ra_XXi	= M_XXi.get();
		delete &M_XXi;
	}
#ifdef USE_LINEAR2
	else if (Xp_a->X_type == REGR_LINEAR2) {
		for (wsUint i=0 ; i<N_col ; i++)
			Xp_a->Ra_XX[N_col-1][i] = sseVV(Xp_a->Ra_Xt[i], N_samp, Xp_a->Ra_Xt[N_col-1]);
		if (OPT_ENABLED(gxe)) {
			for (wsUint i=1 ; i<=N_gxe ; i++)
				for (wsUint j=0; j<=(i+N_gxe) ; j++)
					Xp_a->Ra_XX[i+N_cov][j] = sseVV(Xp_a->Ra_Xt[i+N_cov], N_samp, Xp_a->Ra_Xt[j]);
		}
		cSymMatrix	M_XX(Xp_a->Ra_XX, N_col, 1);
		cSymMatrix&	M_XXi	= M_XX.inv(NULL, &(Xp_a->B_gInv));
		if (!M_XXi.sane()) return;
		M_XXi.setDontDealloc();

		Ra_XXi	= M_XXi.get();
		delete &M_XXi;
	}
#endif
	for (wsUint k=0 ; k<N_pheno ; k++) {
		cVector	V_y(Ra_Ys[k], N_samp, 1);
		wsReal*	Ra_t	= Xp_a->Ra_t[k];

		switch (Xp_a->X_type) {
		case REGR_LINEAR2:
		case REGR_LINEAR: {
				cSymMatrix	M_XXi(Ra_XXi, N_col, 1);
				cVector		V_Xy;

				if (Xp_a->Ra_swgt) {
					cDiagMatrix	M_swgt(N_samp, Xp_a->Ra_swgt, 1);
					cStdMatrix	M_XW = M_Xt * M_swgt;
					V_Xy	= M_XW * V_y;
					cSymMatrix	M_XX	= M_Xt.MMt(M_swgt);
					cSymMatrix&	M_XXi_	= M_XX.inv();
					wsSym		Ra_XXi_	= M_XXi_.get();
					for (wsUint i=0 ; i<M_XXi.row() ; i++)
						memcpy(Ra_XXi[i], Ra_XXi_[i], sizeof(wsReal)*(i+1));
					delete &M_XXi_;
				} else {
#ifdef USE_LINEAR2
					if (Xp_a->X_type == REGR_LINEAR2) {
						Xp_a->Ra_Xy[k][N_col-1] = sseVV(Ra_Ys[k], N_samp, Xp_a->Ra_Xt[N_col-1]);
						if (OPT_ENABLED(gxe)) for (wsUint i=1 ; i<=N_gxe ; i++)
							Xp_a->Ra_Xy[k][N_cov + i] = sseVV(Ra_Ys[k], N_samp, Xp_a->Ra_Xt[i+N_cov]);
						V_Xy.init(N_col, Xp_a->Ra_Xy[k], NULL, 0, 1);
					} else
#endif
						V_Xy	= M_Xt * V_y;
				}
				cVector		V_b		= M_XXi * V_Xy;
				wsReal*		Ra_b	= V_b.get();
				V_b.setDontDealloc();
				Xp_a->Ra_b[k]		= Ra_b;
				cVector		V_res	= M_Xt.tV(V_b);
				V_res.subFrom(V_y);
				if (Xp_a->Ra_swgt) {
					wsVec Ra_res = V_res.get();
					sseVpVsqrt(Ra_res, Xp_a->Ra_swgt, Ra_res, V_res.size());
				}
				wsReal		R_sig2	= V_res.ss() / R_tDF;

				/* Perform Wald test to each beta */
				for (wsUint i=1 ; i<N_col ; i++)
					Ra_t[i] = Ra_b[i] / sqrt(R_sig2*Ra_XXi[i][i]);
			} break;
		case REGR_LOGISTIC: {
				wsReal*	Ra_V	= sseVector(N_samp);
				Xp_a->Ra_b[k]	= sseEmptyVec(N_col);
				wsReal*	Ra_b	= Xp_a->Ra_b[k];
				bool	B_conv	= false;
				Xp_a->N_iter = (short)cRegrAnalysis::fitLogisticByNR(N_samp, N_col,
					Xp_a->Ra_Xt, Xp_a->Ra_Ys[k], Xp_a->Ra_b[k], Ra_V, 20, B_conv);
				cDiagMatrix	M_V(N_samp, Ra_V);

				// S <- solve( t(X) %*% V %*% X )
				cSymMatrix	M_S		= M_Xt.MMt(M_V);
				cSymMatrix&	M_Si	= M_S.inv();
				if (!M_Si.sane()) return;
				Ra_Si				= M_Si.get();

				/* Perform Wald test to each beta */
				for (wsUint i=1 ; i<N_col ; i++)
					Ra_t[i] = Ra_b[i] / sqrt(Ra_Si[i][i]);
				delete &M_Si;
			} break;
		} /* END OF switch, regression type */
	}

	sseUnmat(Ra_XXi, N_col);
}

void cRegrAnalysis::_regressionMain(cExporter &C_regr, wsMat Ra_pvals,
	int** Na_chrss,
	wsUint *Na_pval, int N_chr/*=-1*/)
{
	bool		B_adjust	= OPT_ENABLED(adjust);
	vVariant&	Xa_SNP		= Cp_IO->getVariant();
	xMaf*		Xp_maf		= Cp_IO->getMAF();
	wsUintCst	N_oCov		= Cp_IO->sizeCovar();
	/* --gxe and --gxg are M.E.
	 * so N_gxeCov == 1 when --gxg (== N_chr!=-1), for GxG */
	wsUintCst	N_gxeCov	= N_chr != -1 ? 1 : Cp_IO->sizeGxEcovar();
	vPheno&		Xa_phe		= getIO()->getPhenoInfo();
	wsUint		N_idxSNP	= 0;
	wsUint		i;

	xAnaRegr X = { X_type, 0 };
	X.N_col		= N_col;
	X.N_gxe		= N_gxeCov;
	X.N_pheno	= N_pheno;
	wsAlloc(X.Ra_b, wsReal*, N_pheno);
	X.Ra_t		= sseMatrix(N_pheno, N_col);
	X.Ra_swgt	= getIO()->getSampleWeight();

#ifdef USE_LINEAR2
	X.N_idxSNP	= N_idxSNP;	// 140512 _linear2
	wsAlloc(X.Ra_Xy, wsReal*, N_pheno);
	for (wsUint i=0 ; i<N_pheno ; i++)
		X.Ra_Xy[i] = sseMpV(Ra_anaXt, N_col, N_anaSamp, Ra_anaYs[i]);
	X.Ra_XX = sseMpMt(Ra_anaXt, N_col, N_anaSamp);
#endif

	vector<double>	*Xv_chis = new vector<double>[N_pheno];
	/* Annotation */
	char*	S_bufAnno		= NULL;
	wsAlloc(S_bufAnno, char, 65536);
	FOREACHDO (vVariant_it, Xa_SNP, it, N_idxSNP++) {
		/* If --gxg then N_chr is NOT 1, so filtering it */
		if (N_chr != -1 && it->chr != N_chr) continue;

		wsUint	N_curSamp	= 0;
		wsMat	Ra_curXt	= NULL;
		wsMat	Ra_curYs	= NULL;
		wsUint	N_imp		= 0xffffffff;
		sseMinit(X.Ra_t, N_pheno, N_col, WISARD_NAN);

		wsReal	R_tDF		= WISARD_NAN; // REGR_LINEAR
		memset(X.Ra_b, 0x00, sizeof(wsReal*)*N_pheno);

		/* Annotation */
		if (IS_ASSIGNED(annogene))
			sprintf(S_bufAnno,"	%s", it->anno);

		X.N_iter	= -1;
		X.B_gInv	= -1;
		if (IS_ASSIGNED(dosage) || (!IS_ASSIGNED(dosage) && Xp_maf[N_idxSNP].R_maf != W0)) {
			N_curSamp	= _setXs(*it, N_idxSNP, Ra_anaXt, Ra_anaYs,
				&Ra_curXt, &Ra_curYs, &N_imp);
			if (N_curSamp) {
				X.N_samp	= N_curSamp;
				X.Ra_Xt		= Ra_curXt;
				X.Ra_Ys		= Ra_curYs;
				X.X_type	= X_type;
#ifdef USE_LINEAR2
				if (N_curSamp == N_anaSamp && X_type == REGR_LINEAR)
					X.X_type	= REGR_LINEAR2;
#endif
				_regression(&X);

				/* Set DF */
				R_tDF		= (wsReal)(N_curSamp - N_col - 1);
			}
		}

		for (wsUint k=0 ; k<N_pheno ; k++) {
			wsReal*	Ra_b	= X.Ra_b[k];
			wsReal*	Ra_t	= X.Ra_t[k];
			wsReal	R_pval	= Ra_t[N_col-1];
			char	B_na	= NA(R_pval);

			/* --remna */
			if (OPT_ENABLED(remna) && B_na) continue;

			/* --pvalrange */
			if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), R_pval))
				continue;

			char	S_bufImp[256]		= { 0, };
			char	S_bufAdd[256]		= { 0, };
			char	S_bufRes[256]		= { 0, };

			/* Imputation buffer */
			if (N_imp == 0xffffffff) strcpy(S_bufImp, "<NA>");
			else sprintf(S_bufImp, "%d", N_imp);

			/* Iter buffer */
			if ((X_type == REGR_LOGISTIC && X.N_iter == -1) ||
				(X_type == REGR_LINEAR && X.B_gInv == -1))
				strcpy(S_bufAdd, "<NA>");
			else sprintf(S_bufAdd, "%d", X_type==REGR_LOGISTIC ? X.N_iter : X.B_gInv);

			/* Result buffer */
			if (B_na) strcpy(S_bufRes, "NA	NA	NA");
			else {
				wsReal R_Plogistic = PVchisq(SQR(Ra_t[N_col-1]), 1.0);
				switch (X_type) {
				case REGR_LOGISTIC:
					R_Plogistic = PVchisq(SQR(Ra_t[N_col-1]), 1.0); break;
				case REGR_LINEAR2:
				case REGR_LINEAR:
					R_Plogistic = ptdist(Ra_t+(N_col-1), &R_tDF); break;
				}

				/* [[R]] Insert */
				Ra_pvals[k][Na_pval[k]]		= R_Plogistic;
				Na_chrss[k][Na_pval[k]++]	= it->chr;

				sprintf(S_bufRes, "%g	%g	%g", Ra_b[N_col-1],
					Ra_t[N_col-1], R_Plogistic);
			}
			/* For multiple comparsion */
			if (B_adjust) Xv_chis[k].push_back(SQR(Ra_t[N_col-1]));

			entryVariant(&C_regr, *it);
			if (N_pheno > 1) {
				const char *S_p = Xa_phe[k].S_name.c_str();
				C_regr.fmt("%s	", S_p);
			}
			C_regr.fmt("	%s	%s	%s", S_bufImp, S_bufAdd, S_bufRes);

			/* Print covariate test */
			for (i=1 ; i<=N_oCov ; i++) {
				if (NA(Ra_t[i])) C_regr.put("	NA	NA	NA");
				else {
					wsReal R_pValCov = WISARD_NAN;
					switch (X_type) {
					case REGR_LOGISTIC:
						R_pValCov = PVchisq(SQR(Ra_t[i]), 1.0); break;
					case REGR_LINEAR2:
					case REGR_LINEAR:
						R_pValCov = ptdist(Ra_t+i, &R_tDF); break;
					}
					C_regr.fmt("	%g	%g	%g", Ra_b[i], Ra_t[i],
						R_pValCov);
				}
			}

			/* Print interaction test (covGxE) */
			if (OPT_ENABLED(gxe)) for (i=1 ; i<=N_gxeCov ; i++) {
				if (NA(Ra_t[i+N_oCov]))
					C_regr.put("	NA	NA	NA");
				else {
					wsReal R_pValCov = WISARD_NA;
					switch (X_type) {
					case REGR_LOGISTIC:
						R_pValCov = PVchisq(SQR(Ra_t[i+N_oCov]), 1.0); break;
					case REGR_LINEAR2:
					case REGR_LINEAR:
						R_pValCov = ptdist(Ra_t+(i+N_oCov), &R_tDF); break;
					}
					C_regr.fmt("	%g	%g	%g", Ra_b[i+N_oCov],
						Ra_t[i+N_oCov], R_pValCov);
				}
			}

			sseFree(Ra_b);

			/* Supply line feed */
			C_regr.put("\n");
		}

		/* If this is newly allocated one */
		if (N_curSamp != N_anaSamp) {
			sseUnmat(Ra_curYs, N_pheno);
			sseUnmat(Ra_curXt, N_col);
		}

		if ((N_idxSNP%100) == 0)
			notice("%d/%d variants tested...\r", N_idxSNP, Xa_SNP.size());
	}
	DEALLOC(S_bufAnno);
	DEALLOC(X.Ra_b);
	sseUnmat(X.Ra_t, N_pheno);
	//X.Ra_t		= sseMatrix(N_pheno, N_col);
#ifdef USE_LINEAR2
	sseUnmat(X.Ra_XX, N_col);
	sseUnmat(X.Ra_Xy, N_pheno);
#endif

	/* Perform multiple comparison */
	if (B_adjust) for (wsUint i=0 ; i<N_pheno ; i++) {
		char	S_type[32];
		vPheno&	Xa_phe = getIO()->getPhenoInfo();

		/* Set type */
		switch (X_type) {
		case REGR_LOGISTIC:
			strcpy(S_type, "logistic"); break;
		case REGR_LINEAR2:
		case REGR_LINEAR:
			strcpy(S_type, "linear"); break;
		}

		/* Make name */
		char S_buf[512] = { 0, };
		if (N_pheno == 1)
			sprintf(S_buf, "%s.regr.adjust.res", S_type);
		else
			sprintf(S_buf, "%s.regr.%s.adjust.res", S_type, Xa_phe[i].S_name.c_str());
		multcomp(getIO(), Xv_chis[i], S_buf);
	}
	delete [] Xv_chis;

	LOG("%d/%d variants tested...\n", N_idxSNP, Xa_SNP.size());
}

void cRegrAnalysis::_doNull(wsMat Ra_curXt)
{
	wsUint N_covCol = N_col - 1;
//	wsReal	*Ra_betaHat	= NULL;
	wsMat	Ra_betaHats	= NULL;
	wsAlloc(Ra_betaHats, wsReal*, N_pheno);

	wsMat	Ra_resids	= NULL;
	if (B_isCont) {
		/* (X'X)^-1 */
		wsReal **Ra_XtX = sseMpMt(Ra_curXt, N_covCol, N_anaSamp,
			Ra_curXt, N_covCol, N_anaSamp); /* checked */
		wsReal **Ra_XtX_inv = SO_invMatrix(Ra_XtX, N_covCol); /* checked */
		if (Ra_XtX_inv == NULL) {
			LOGwarn("Null test failed\n");
			sseUnmat(Ra_XtX, N_covCol);
			sseUnmat(Ra_XtX_inv, N_covCol);
			return;
		}
		sseUnmat(Ra_XtX, N_covCol);

		/* X'y */
//		wsReal *Ra_Xty = sseMpV(Ra_curXt, N_covCol, N_anaSamp, Ra_anaY); /* checked */

		for (wsUint i=0 ; i<N_pheno ; i++) {
			wsReal *Ra_Xty = sseMpV(Ra_curXt, N_covCol, N_anaSamp, Ra_anaYs[i]);
			Ra_betaHats[i] = sseMpV(Ra_XtX_inv, N_covCol, N_covCol, Ra_Xty);
			sseFree(Ra_Xty);
		}

		/* Get beta estimate ^beta = (X'X)^-1X'y */
// 		Ra_betaHat = sseMpV(Ra_XtX_inv, N_covCol, N_covCol,
// 			Ra_Xty); /* checked */
//		sseFree(Ra_Xty);

		/* Get residual e = y- X %*% ^beta */
		// 	wsReal *Ra_resid = sseMtV(Ra_curXt, N_covCol, N_anaSamp, Ra_betaHat, N_covCol); /* checked */
		// 	sseVsV(Ra_anaY, Ra_resid, Ra_resid, N_anaSamp);
		wsAlloc(Ra_resids, wsReal*, N_pheno);
		for (wsUint i=0 ; i<N_pheno ; i++) {
			Ra_resids[i] = sseMtV((wsMatCst)Ra_curXt, N_covCol, N_anaSamp, Ra_betaHats[i], N_covCol); /* checked */
			sseVsV(Ra_anaYs[i], Ra_resids[i], Ra_resids[i], N_anaSamp);
		}
	} else {
		wsReal	*Ra_V		= sseVector(N_anaSamp); /* checked */
//		Ra_betaHat = sseVector(N_covCol);
		Ra_betaHats = sseEmptyMatrix(N_pheno, N_covCol);

		/* Initialize */
// 		for (wsUint i=0 ; i<N_col ; i++)
// 			Ra_betaHat[i] = W0;

		bool B_conv = false;
		for (wsUint i=0 ; i<N_pheno ; i++)
			fitLogisticByNR(N_anaSamp, N_covCol, Ra_curXt, Ra_anaYs[i], Ra_betaHats[i],
				Ra_V, 20, B_conv);

		/* Get residual e = y- X %*% ^beta */
		/* Get deviance residual */
		// 	wsReal *Ra_resid = sseMtV(Ra_curXt, N_covCol, N_anaSamp, Ra_betaHat, N_covCol); /* checked */
		// 	sseVsV(Ra_anaY, Ra_resid, Ra_resid, N_anaSamp);
		Ra_resids	= sseMatrix(N_pheno, N_anaSamp);
		for (wsUint i=0 ; i<N_pheno ; i++) {
			cVector		V_beta(Ra_betaHats[i], N_covCol, 1);
			cStdMatrix	M_Xt(N_covCol, N_anaSamp, Ra_curXt, MATDEL_NONE);
			cVector		V_predict = V_beta * M_Xt;
			wsReal*		Ra_pred =V_predict.get();

			/* Convert to pi */
			for (wsUint j=0 ; j<N_anaSamp ; j++)
				Ra_pred[j] = W1 / (W1 + exp(-Ra_pred[j]));
			/* Get deviance residual */
			for (wsUint j=0 ; j<N_anaSamp ; j++) {
				if (Ra_anaYs[i][j] == WISARD_AFFECTED)
					Ra_resids[i][j] = sqrt(REAL_CONST(-2.0) * log(Ra_pred[j]));
				else
					Ra_resids[i][j] = -sqrt(REAL_CONST(-2.0) * log(W1 - Ra_pred[j]));
			}
		}

		sseFree(Ra_V);
	}

	/* Export bhat */
	vCovar&		Xa_cov	= getIO()->getCovInfo();
	cExporter*	Cp_bHat	= cExporter::summon(B_isCont?"beta0.linear.regr.res":
			"beta0.logistic.regr.res");
	vPheno&		Xa_phe	= getIO()->getPhenoInfo();
	for (wsUint k=0 ; k<N_pheno ; k++)
		for (wsUint i=0 ; i<N_covCol ; i++) {
			char S_varN[512];
			if (i == 0) sprintf(S_varN, "%s.INTERCEPT", Xa_phe[k].S_name.c_str());
			else strcpy(S_varN, Xa_cov[i-1].Sp_varName);
			Cp_bHat->fmt("%s.%s	%g\n", Xa_phe[k].S_name.c_str(), S_varN,
				Ra_betaHats[k][i]);
		}
	delete Cp_bHat;

	/* Export residual */
//	wsUint N_samp = getIO()->getSampleSize();
	//char *S_varName = "RESIDUAL";
	char **Sa_varNames = NULL;
	wsAlloc(Sa_varNames, char *, N_pheno);
	for (wsUint i=0 ; i<N_pheno ; i++) {
		Sa_varNames[i] = new char[512];
		sprintf(Sa_varNames[i], "%s.RESIDUAL", Xa_phe[i].S_name.c_str());
	}
// 	exportSampleMatrix(B_isCont?"res0.linear.regr.res":"res0.logistic.regr.res",
// 		getIO(), Ba_isExcl, (char)1, Na_sampIdx, 1, &S_varName, &Ra_resid);
	exportSampleMatrix(B_isCont?"res0.linear.regr.res":"res0.logistic.regr.res",
		getIO(), Ba_isExcl, (char)1, N_pheno, Sa_varNames, Ra_resids);
//	sseFree(Ra_resid);
// 	sseFree(Ra_betaHat);
	sseUnmat(Ra_resids, N_pheno);
	sseUnmat(Ra_betaHats, N_pheno);
}

#endif

} // End namespace ONETOOL
