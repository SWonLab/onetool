#include "analyses/gemma.h"
#include "utils/matrix.h"
#include "analyses/annot.h"
#include "utils/dcdflib.h"
#include "input/impute.h"
#include "utils/stat.h"

namespace ONETOOL {

class FUNC_PARAM
{

public:
	bool calc_null;
	wsUint ni_test;
	wsUint n_cvt;
	cVector& eval;
	cStdMatrix& Uab;
	cVector& ab;
	size_t e_mode;
	wsMat Pab, Iab;
	wsMat PPab, PPPab;
};

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

#define R_2 W2

struct _xGemma;
typedef wsReal(*hRcsvFun)(wsUint i, wsUint j, wsUint k, wsStrCst S_lb,
	wsUint N_level, struct _xGemma &X_data);
typedef wsReal(*hTermFun)(wsUint i, wsUint k, struct _xGemma &X_data);
typedef pair<hTermFun, hRcsvFun>	xPairFun;
typedef map<string,xPairFun>		mRcsvFun;
typedef mRcsvFun::iterator			mRcsvFun_it;

wsUint GetabIndex(const wsUint a, const wsUint b, const wsUint n_cvt);

inline void CalcPab(const wsUint n_cvt, const size_t e_mode, cVector& Hi_eval, cStdMatrix& Uab_t, cVector& _ab, wsMat Pab)
{
	wsUint index_ab, index_aw, index_bw, index_ww;
	double p_ab;
	double ps_ab, ps_aw, ps_bw, ps_ww;
	wsVec ab = _ab.get();

	for (wsUint p=0; p<=n_cvt+1; ++p) {
		for (wsUint a=p+1; a<=n_cvt+2; ++a) {
			for (wsUint b=a; b<=n_cvt+2; ++b) {
				index_ab=GetabIndex(a, b, n_cvt);
				if (p==0) {
					cVector Uab_t_ = Uab_t.r2v_ptr(index_ab);
					p_ab = Hi_eval.sum(Uab_t_);
					if (e_mode!=0) { p_ab=ab[index_ab]-p_ab; }
					Pab[0][index_ab] = p_ab;
				}
				else {
					index_aw=GetabIndex(a, p, n_cvt);
					index_bw=GetabIndex(b, p, n_cvt);
					index_ww=GetabIndex(p, p, n_cvt);

					ps_ab=Pab[p-1][index_ab];
					ps_aw=Pab[p-1][index_aw];
					ps_bw=Pab[p-1][index_bw];
					ps_ww=Pab[p-1][index_ww];

					p_ab=ps_ab-ps_aw*ps_bw/ps_ww;
					Pab[p][index_ab] = p_ab;
				}
			}
		}
	}
	return;
}

inline void CalcPPab(const wsUint n_cvt, const size_t e_mode, cVector& HiHi_eval, cStdMatrix& Uab_t, cVector& _ab, wsMat Pab, wsMat PPab)
{
	wsUint index_ab, index_aw, index_bw, index_ww;
	double p2_ab;
	double ps2_ab, ps_aw, ps_bw, ps_ww, ps2_aw, ps2_bw, ps2_ww;
	wsVec ab = _ab.get();

	for (wsUint p=0; p<=n_cvt+1; ++p) {
		for (wsUint a=p+1; a<=n_cvt+2; ++a) {
			for (wsUint b=a; b<=n_cvt+2; ++b) {
				index_ab=GetabIndex(a, b, n_cvt);
				if (p==0) {
					cVector tmp = Uab_t.r2v_ptr(index_ab);
					p2_ab = HiHi_eval.sum(tmp);
					if (e_mode!=0) { p2_ab=p2_ab-ab[index_ab]+2.0*Pab[0][index_ab]; }
					PPab[0][index_ab] = p2_ab;
				}
				else {
					index_aw=GetabIndex(a, p, n_cvt);
					index_bw=GetabIndex(b, p, n_cvt);
					index_ww=GetabIndex(p, p, n_cvt);

					ps2_ab=PPab[p-1][index_ab];
					ps_aw=Pab[p-1][index_aw];
					ps_bw=Pab[p-1][index_bw];
					ps_ww=Pab[p-1][index_ww];
					ps2_aw=PPab[p-1][index_aw];
					ps2_bw=PPab[p-1][index_bw];
					ps2_ww=PPab[p-1][index_ww];

					p2_ab=ps2_ab+ps_aw*ps_bw*ps2_ww/(ps_ww*ps_ww);
					p2_ab-=(ps_aw*ps2_bw+ps_bw*ps2_aw)/ps_ww;
					PPab[p][index_ab] = p2_ab;

				}
			}
		}
	}
	return;
}


inline void CalcPPPab(const wsUint n_cvt, const size_t e_mode, cVector& HiHiHi_eval, cStdMatrix& Uab_t, cVector& _ab,
	wsMat Pab, wsMat PPab, wsMat PPPab)
{
	wsUint index_ab, index_aw, index_bw, index_ww;
	double p3_ab;
	double ps3_ab, ps_aw, ps_bw, ps_ww, ps2_aw, ps2_bw, ps2_ww, ps3_aw, ps3_bw, ps3_ww;
	wsVec ab = _ab.get();

	for (wsUint p=0; p<=n_cvt+1; ++p) {
		for (wsUint a=p+1; a<=n_cvt+2; ++a) {
			for (wsUint b=a; b<=n_cvt+2; ++b) {
				index_ab=GetabIndex(a, b, n_cvt);
				if (p==0) {
					cVector tmp = Uab_t.r2v_ptr(index_ab);
					p3_ab = HiHiHi_eval.sum(tmp);
					if (e_mode!=0) { p3_ab=ab[index_ab]-p3_ab+3.0*PPab[0][index_ab]-3.0*Pab[0][index_ab]; }
					PPPab[0][index_ab] = p3_ab;
				}
				else {
					index_aw=GetabIndex(a, p, n_cvt);
					index_bw=GetabIndex(b, p, n_cvt);
					index_ww=GetabIndex(p, p, n_cvt);

					ps3_ab=PPPab[p-1][index_ab];
					ps_aw=Pab[p-1][index_aw];
					ps_bw=Pab[p-1][index_bw];
					ps_ww=Pab[p-1][index_ww];
					ps2_aw=PPab[p-1][index_aw];
					ps2_bw=PPab[p-1][index_bw];
					ps2_ww=PPab[p-1][index_ww];
					ps3_aw=PPPab[p-1][index_aw];
					ps3_bw=PPPab[p-1][index_bw];
					ps3_ww=PPPab[p-1][index_ww];

					p3_ab=ps3_ab-ps_aw*ps_bw*ps2_ww*ps2_ww/(ps_ww*ps_ww*ps_ww);
					p3_ab-=(ps_aw*ps3_bw+ps_bw*ps3_aw+ps2_aw*ps2_bw)/ps_ww;
					p3_ab+=(ps_aw*ps2_bw*ps2_ww+ps_bw*ps2_aw*ps2_ww+ps_aw*ps_bw*ps3_ww)/(ps_ww*ps_ww);

					PPPab[p][index_ab] = p3_ab;
				}
			}
		}
	}
	return;
}

inline wsUint GetabIndex(const wsUint a, const wsUint b, const wsUint n_cvt) {
	if (a>n_cvt+2 || b>n_cvt+2 || a<=0 || b<=0) {
		halt("error in GetabIndex.");
		return 0;
	}
	wsUint index;
	wsUint l, h;
	if (b>a) { l=a; h=b; }
	else { l=b; h=a; }

	wsUint n=n_cvt+2;
	index=(2*n-l+2)*(l-1)/2+h-l;

	return index;
}


void CalcRLWald(const double &l, const FUNC_PARAM &params, double &beta, double &se,
	double &p_wald, double &t_wald)
{
	wsUint n_cvt=params.n_cvt;
	wsUint n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	wsUint sz = params.eval.size();

	int df=(int)sz-(int)n_cvt-1;

	cStdMatrix Pab(n_cvt+2, n_index);
	cVector Hi_eval(sz);
	cVector v_temp = params.eval * l;

	if (params.e_mode==0) { Hi_eval.set(1.0); }
	else { memcpy(Hi_eval.get(), v_temp.get(), sizeof(wsReal)*sz); }
	v_temp += 1.0;
	Hi_eval /= v_temp;

	CalcPab(n_cvt, params.e_mode, Hi_eval, params.Uab, params.ab, Pab.get());

	wsUint index_yy=GetabIndex(n_cvt+2, n_cvt+2, n_cvt);
	wsUint index_xx=GetabIndex(n_cvt+1, n_cvt+1, n_cvt);
	wsUint index_xy=GetabIndex(n_cvt+2, n_cvt+1, n_cvt);
	double P_yy=Pab[n_cvt][index_yy];
	double P_xx=Pab[n_cvt][index_xx];
	double P_xy=Pab[n_cvt][index_xy];
	double Px_yy=Pab[n_cvt+1][index_yy];

	beta=P_xy/P_xx;
	double tau=(double)df/Px_yy;
	se=sqrt(1.0/(tau*P_xx));

	t_wald=(P_yy-Px_yy)*tau;
	p_wald=pf(t_wald, 1.0, df, 0);
}


void CalcRLScore(const double &l, const FUNC_PARAM &params,
	double &beta, double &se, double &p_score, double &t_score)
{
	wsUint n_cvt=params.n_cvt;
	wsUint n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	wsUint sz = params.eval.size();

	int df=(int)sz-(int)n_cvt-1;

	cStdMatrix Pab(n_cvt+2, n_index);
	cVector Hi_eval(sz);
	cVector v_temp = params.eval * l;

	if (params.e_mode==0) { Hi_eval.set(1.0); }
	else { memcpy(Hi_eval.get(), v_temp.get(), sizeof(wsReal)*sz); }
	v_temp += 1.0;
	Hi_eval /= v_temp;

	CalcPab(n_cvt, params.e_mode, Hi_eval, params.Uab, params.ab, Pab.get());

	wsUint index_yy=GetabIndex(n_cvt+2, n_cvt+2, n_cvt);
	wsUint index_xx=GetabIndex(n_cvt+1, n_cvt+1, n_cvt);
	wsUint index_xy=GetabIndex(n_cvt+2, n_cvt+1, n_cvt);
	double P_yy=Pab[n_cvt][index_yy];
	double P_xx=Pab[n_cvt][index_xx];
	double P_xy=Pab[n_cvt][index_xy];
	double Px_yy=Pab[n_cvt+1][index_yy];

	beta=P_xy/P_xx;
	double tau=(double)df/Px_yy;
	se=sqrt(1.0/(tau*P_xx));

	t_score=(double)sz*P_xy*P_xy/(P_yy*P_xx);
	p_score=pf(t_score, 1.0, df, 0);
}

inline double LogL_dev1(double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *)params;
	wsUint n_cvt=p->n_cvt;
	wsUint ni_test=p->ni_test;
	wsUint sz = p->eval.size();

	wsUint nc_total;
	if (p->calc_null==true) { nc_total=n_cvt; }
	else { nc_total=n_cvt+1; }

	double dev1=0.0, trace_Hi=0.0;
	wsUint index_yy;

	cVector Hi_eval(sz);
	cVector HiHi_eval(sz);
	cVector v_temp = p->eval * l;

	if (p->e_mode==0) { Hi_eval.set(1.0); }
	else { memcpy(Hi_eval.get(), v_temp.get(), sizeof(wsReal)*sz); }
	v_temp += 1.0;
	Hi_eval /= v_temp;

	memcpy(HiHi_eval.get(), Hi_eval.get(), sizeof(wsReal)*sz);
	HiHi_eval *= Hi_eval;

	trace_Hi = Hi_eval.sum();

	if (p->e_mode!=0) { trace_Hi=(double)ni_test-trace_Hi; }

	CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, p->Pab);
	CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, p->Pab, p->PPab);

	double trace_HiK=((double)ni_test-trace_Hi)/l;

	index_yy=GetabIndex(n_cvt+2, n_cvt+2, n_cvt);

	double P_yy=p->Pab[nc_total][index_yy];
	double PP_yy=p->PPab[nc_total][index_yy];
	double yPKPy=(P_yy-PP_yy)/l;
	dev1=-0.5*trace_HiK+0.5*(double)ni_test*yPKPy/P_yy;

	return dev1;
}

inline double LogRL_dev1(double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *)params;
	wsUint n_cvt=p->n_cvt;
	wsUint ni_test=p->ni_test;
	wsUint sz = p->eval.size();

	double df;
	wsUint nc_total;
	if (p->calc_null==true) { nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else { nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0; }

	double dev1=0.0, trace_Hi=0.0;
	wsUint index_ww;

	cVector Hi_eval(sz);
	cVector HiHi_eval(sz);
	cVector v_temp = p->eval * l;

	if (p->e_mode==0) { Hi_eval.set(1.0); }
	else { memcpy(Hi_eval.get(), v_temp.get(), sizeof(wsReal)*sz); }
	v_temp += 1.0;
	Hi_eval /= v_temp;

	memcpy(HiHi_eval.get(), Hi_eval.get(), sizeof(wsReal)*sz);
	HiHi_eval *= Hi_eval;

	v_temp.set(1.0);
	trace_Hi = Hi_eval.sum(v_temp);

	if (p->e_mode!=0) {
		trace_Hi=(double)ni_test-trace_Hi;
	}

	CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, p->Pab);
	CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, p->Pab, p->PPab);

	//calculate tracePK and trace PKPK
	double trace_P=trace_Hi;
	double ps_ww, ps2_ww;
	LOOP(i, nc_total) {
		index_ww=GetabIndex(i+1, i+1, n_cvt);
		ps_ww=p->Pab[i][index_ww];
		ps2_ww=p->PPab[i][index_ww];
		trace_P-=ps2_ww/ps_ww;
	}
	double trace_PK=(df-trace_P)/l;

	//calculate yPKPy, yPKPKPy
	index_ww=GetabIndex(n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=p->Pab[nc_total][index_ww];
	double PP_yy=p->PPab[nc_total][index_ww];
	double yPKPy=(P_yy-PP_yy)/l;

	dev1=-0.5*trace_PK+0.5*df*yPKPy/P_yy;

	return dev1;
}

double LogRL_f(double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *)params;
	wsUint n_cvt=p->n_cvt;
	wsUint ni_test=p->ni_test;
	wsUint sz = p->eval.size();

	double df;
	wsUint nc_total;
	if (p->calc_null==true) { nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else { nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0; }

	double f=0.0, logdet_h=0.0, logdet_hiw=0.0, d;
	wsUint index_ww;

	cVector Hi_eval(sz);
	cVector v_temp = p->eval * l;

	if (p->e_mode==0) { Hi_eval.set(1.0); }
	else { memcpy(Hi_eval.get(), v_temp.get(), sizeof(wsReal)*sz); }
	v_temp += 1.0;
	Hi_eval /= v_temp;

	wsVec _v_temp = v_temp.get();
	for (size_t i=0; i<sz; ++i) {
		d=_v_temp[i];
		logdet_h+=log(fabs(d));
	}

	CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, p->Pab);
	v_temp.set(1.0);
	CalcPab(n_cvt, p->e_mode, v_temp, p->Uab, p->ab, p->Iab);

	//calculate |WHiW|-|WW|
	logdet_hiw=0.0;
	for (wsUint i=0; i<nc_total; ++i) {
		index_ww=GetabIndex(i+1, i+1, n_cvt);
		d=p->Pab[i][index_ww];
		logdet_hiw+=log(d);
		d=p->Iab[i][index_ww];
		logdet_hiw-=log(d);
	}
	index_ww=GetabIndex(n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=p->Pab[nc_total][index_ww];

	double c=0.5*df*(log(df)-log(2*WISARD_PI)-1.0);
	f=c-0.5*logdet_h-0.5*logdet_hiw-0.5*df*log(P_yy);

	return f;
}

double LogL_f(double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *)params;
	wsUint n_cvt=p->n_cvt;
	wsUint ni_test=p->ni_test;
	wsUint sz = p->eval.size();

	wsUint nc_total;
	if (p->calc_null==true) { nc_total=n_cvt; }
	else { nc_total=n_cvt+1; }

	double f=0.0, logdet_h=0.0, d;
	wsUint index_yy;

	cVector Hi_eval(sz);
	cVector v_temp = p->eval * l;

	if (p->e_mode==0) { Hi_eval.set(1.0); }
	else { memcpy(Hi_eval.get(), v_temp.get(), sizeof(wsReal)*sz); }
	v_temp += 1.0;
	Hi_eval /= v_temp;

	wsVec _v_temp = v_temp.get(); 
	for (size_t i=0; i<sz; ++i) {
		d=_v_temp[i];
		logdet_h+=log(fabs(d));
	}

	CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, p->Pab);

	double c=0.5*(double)ni_test*(log((double)ni_test)-log(2*WISARD_PI)-1.0);

	index_yy=GetabIndex(n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=p->Pab[nc_total][index_yy];
	f=c-0.5*logdet_h-0.5*(double)ni_test*log(P_yy);

	return f;
}

inline double LogRL_dev2(double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *)params;
	wsUint n_cvt=p->n_cvt;
	wsUint ni_test=p->ni_test;
	wsUint sz = p->eval.size();

	double df;
	wsUint nc_total;
	if (p->calc_null==true) { nc_total=n_cvt; df=(double)ni_test-(double)n_cvt; }
	else { nc_total=n_cvt+1; df=(double)ni_test-(double)n_cvt-1.0; }

	double dev2=0.0, trace_Hi=0.0, trace_HiHi=0.0;
	wsUint index_ww;

	cVector Hi_eval(sz);
	cVector HiHi_eval(sz);
	cVector HiHiHi_eval(sz);
	cVector v_temp = p->eval * l;

	if (p->e_mode==0) { Hi_eval.set(1.0); }
	else { memcpy(Hi_eval.get(), v_temp.get(), sizeof(wsReal)*sz); }
	v_temp += 1.0;
	Hi_eval /= v_temp;

	memcpy(HiHi_eval.get(), Hi_eval.get(), sizeof(wsReal)*sz);
	HiHi_eval *= Hi_eval;
	memcpy(HiHiHi_eval.get(), HiHi_eval.get(), sizeof(wsReal)*sz);
	HiHiHi_eval *= Hi_eval;

	v_temp.set(1.0);
	trace_Hi = Hi_eval.sum(v_temp);
	trace_HiHi = HiHi_eval.sum(v_temp);

	if (p->e_mode!=0) {
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}

	CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, p->Pab);
	CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, p->Pab, p->PPab);
	CalcPPPab(n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, p->Pab, p->PPab, p->PPPab);

	//calculate tracePK and trace PKPK
	double trace_P=trace_Hi, trace_PP=trace_HiHi;
	double ps_ww, ps2_ww, ps3_ww;
	LOOP (i, nc_total) {
		index_ww=GetabIndex(i+1, i+1, n_cvt);
		ps_ww=p->Pab[i][index_ww];
		ps2_ww=p->PPab[i][index_ww];
		ps3_ww=p->PPPab[i][index_ww];
		trace_P-=ps2_ww/ps_ww;
		trace_PP+=ps2_ww*ps2_ww/(ps_ww*ps_ww)-2.0*ps3_ww/ps_ww;
	}
	double trace_PKPK=(df+trace_PP-2.0*trace_P)/(l*l);

	//calculate yPKPy, yPKPKPy
	index_ww=GetabIndex(n_cvt+2, n_cvt+2, n_cvt);
	double P_yy=p->Pab[nc_total][index_ww];
	double PP_yy=p->PPab[nc_total][index_ww];
	double PPP_yy=p->PPPab[nc_total][index_ww];
	double yPKPy=(P_yy-PP_yy)/l;
	double yPKPKPy=(P_yy+PPP_yy-2.0*PP_yy)/(l*l);

	dev2=0.5*trace_PKPK-0.5*df*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);

	return dev2;
}

inline double LogL_dev2(double l, void *params)
{
	FUNC_PARAM *p=(FUNC_PARAM *)params;
	wsUint n_cvt=p->n_cvt;
	wsUint ni_test=p->ni_test;
	wsUint sz = p->eval.size();

	wsUint nc_total;
	if (p->calc_null==true) { nc_total=n_cvt; }
	else { nc_total=n_cvt+1; }

	double dev2=0.0, trace_Hi=0.0, trace_HiHi=0.0;
	wsUint index_yy;

	cVector Hi_eval(sz);
	cVector HiHi_eval(sz);
	cVector HiHiHi_eval(sz);
	cVector v_temp = p->eval * l;

	if (p->e_mode==0) { Hi_eval.set(1.0); }
	else { memcpy(Hi_eval.get(), v_temp.get(), sizeof(wsReal)*sz); }
	v_temp += 1.0;
	Hi_eval /= v_temp;

	memcpy(HiHi_eval.get(), Hi_eval.get(), sizeof(wsReal)*sz);
	HiHi_eval *= Hi_eval;
	memcpy(HiHiHi_eval.get(), HiHi_eval.get(), sizeof(wsReal)*sz);
	HiHiHi_eval *= Hi_eval;

	v_temp.set(1.0);
	trace_Hi = Hi_eval.sum(v_temp);
	trace_HiHi = HiHi_eval.sum(v_temp);

	if (p->e_mode!=0) {
		trace_Hi=(double)ni_test-trace_Hi;
		trace_HiHi=2*trace_Hi+trace_HiHi-(double)ni_test;
	}

	CalcPab(n_cvt, p->e_mode, Hi_eval, p->Uab, p->ab, p->Pab);
	CalcPPab(n_cvt, p->e_mode, HiHi_eval, p->Uab, p->ab, p->Pab, p->PPab);
	CalcPPPab(n_cvt, p->e_mode, HiHiHi_eval, p->Uab, p->ab, p->Pab, p->PPab, p->PPPab);

	double trace_HiKHiK=((double)ni_test+trace_HiHi-2*trace_Hi)/(l*l);

	index_yy=GetabIndex(n_cvt+2, n_cvt+2, n_cvt);
	double P_yy		= p->Pab[nc_total][index_yy];
	double PP_yy	= p->PPab[nc_total][index_yy];
	double PPP_yy	= p->PPPab[nc_total][index_yy];

	double yPKPy	= (P_yy-PP_yy)/l;
	double yPKPKPy	= (P_yy+PPP_yy-2.0*PP_yy)/(l*l);

	dev2=0.5*trace_HiKHiK-0.5*(double)ni_test*(2.0*yPKPKPy*P_yy-yPKPy*yPKPy)/(P_yy*P_yy);

	return dev2;
}





void CalcUab(cStdMatrix& UtW, cVector& Uty, cVector& Utx, cStdMatrix& Uab_t)
{
	wsUint index_ab;
	wsUint n_cvt=UtW.row();
	wsUint N_samp = Uty.size();

	for (wsUint b=1; b<=n_cvt+2; ++b) {
		index_ab=GetabIndex(n_cvt+1, b, n_cvt);
		wsVec Uab_col= Uab_t[index_ab];

		if (b==n_cvt+2) { memcpy(Uab_col, Uty.get(), sizeof(wsReal)*N_samp); }
		else if (b==n_cvt+1) { memcpy(Uab_col, Utx.get(), sizeof(wsReal)*N_samp); }
		else {
			wsVec UtW_col=UtW[b-1];
			memcpy(Uab_col, UtW_col, sizeof(wsReal)*N_samp);
		}

		sseVpV(Uab_col, Utx.get(), Uab_col, N_samp);
	}
}

void CalcUab(cStdMatrix& UtW, cVector& Uty, cStdMatrix& Uab_t)
{
	wsUint index_ab;
	wsUint n_cvt=UtW.row();
	wsUint N_samp = Uty.size();
	cVector u_a(Uty.size());

	for (wsUint a=1; a<=n_cvt+2; ++a) {
		if (a==n_cvt+1) { continue; }

		if (a==n_cvt+2) { memcpy(u_a.get(), Uty.get(), sizeof(wsReal)*N_samp); }
		else
			memcpy(u_a.get(), UtW.get()[a-1], sizeof(wsReal)*N_samp);

		for (wsUint b=a; b>=1; --b) {
			if (b==n_cvt+1) { continue; }

			index_ab=GetabIndex(a, b, n_cvt);
			wsVec Uab_col = Uab_t.get()[index_ab];

			if (b==n_cvt+2) { memcpy(Uab_col, Uty.get(), sizeof(wsReal)*N_samp); }
			else
				memcpy(Uab_col, UtW.get()[b-1], sizeof(wsReal)*N_samp);

			sseVpV(Uab_col, u_a.get(), Uab_col, N_samp);
		}
	}
}


void CalcLambda(const char func_name, FUNC_PARAM &params, const double l_min, const double l_max,
	const size_t n_region, double &lambda, double &logf)
{
	if (func_name!='R' && func_name!='L' && func_name!='r' && func_name!='l') {
		halt("func_name only takes 'R' or 'L': 'R' for log-restricted likelihood, 'L' for log-likelihood.");
		return;
	}

	vector<pair<double, double> > lambda_lh;

	//evaluate first order derivates in different intervals
	double lambda_l, lambda_h, lambda_interval=log(l_max/l_min)/(double)n_region;
	double dev1_l, dev1_h, logf_l, logf_h;

	for (wsUint i=0; i<n_region; ++i) {
		lambda_l=l_min*exp(lambda_interval*i);
		lambda_h=l_min*exp(lambda_interval*(i+1.0));

		if (func_name=='R' || func_name=='r') {
			dev1_l=LogRL_dev1(lambda_l, &params);
			dev1_h=LogRL_dev1(lambda_h, &params);
		}
		else {
			dev1_l=LogL_dev1(lambda_l, &params);
			dev1_h=LogL_dev1(lambda_h, &params);
		}

		if (dev1_l*dev1_h<=0) {
			lambda_lh.push_back(make_pair(lambda_l, lambda_h));
		}
	}

	//if derivates do not change signs in any interval
	if (lambda_lh.empty()) {
		if (func_name=='R' || func_name=='r') {
			logf_l=LogRL_f(l_min, &params);
			logf_h=LogRL_f(l_max, &params);
		}
		else {
			logf_l=LogL_f(l_min, &params);
			logf_h=LogL_f(l_max, &params);
		}

		if (logf_l>=logf_h) { lambda=l_min; logf=logf_l; }
		else { lambda=l_max; logf=logf_h; }
	}
	else {
		//if derivates change signs
		//int status;
		int iter=0, max_iter=100;
		double l;

		double (*f)(double, void *) = func_name=='R' || func_name=='r' ? LogRL_dev1 : LogL_dev1;
		double (*df)(double, void *) = func_name=='R' || func_name=='r' ? LogRL_dev2 : LogL_dev2;

		for (vector<double>::size_type i=0; i<lambda_lh.size(); ++i) {
			lambda_l=lambda_lh[i].first; lambda_h=lambda_lh[i].second;

			// Try BRENT first
			l = Brent_fmin(lambda_l, lambda_h, LogL_dev1, &params, 1e-1);

			// Apply NR next
			double h;
			do {
				iter++;
				h = f(l, &params) / df(l, &params);
				l -= h;
			} while (fabs(h)>=1e-5 && iter<max_iter && l>l_min && l<l_max);

			if (l<l_min) { l=l_min; }
			if (l>l_max) { l=l_max; }
			if (func_name=='R' || func_name=='r') { logf_l=LogRL_f(l, &params); }
			else { logf_l=LogL_f(l, &params); }

			if (i==0) { logf=logf_l; lambda=l; }
			else if (logf<logf_l) { logf=logf_l; lambda=l; }
			else {}
		}

		if (func_name=='R' || func_name=='r') {
			logf_l=LogRL_f(l_min, &params);
			logf_h=LogRL_f(l_max, &params);
		}
		else {
			logf_l=LogL_f(l_min, &params);
			logf_h=LogL_f(l_max, &params);
		}

		if (logf_l>logf) { lambda=l_min; logf=logf_l; }
		if (logf_h>logf) { lambda=l_max; logf=logf_h; }
	}

	return;
}

void CalcLambda(const char func_name, cVector& eval, cStdMatrix& UtW, cVector& Uty, double l_min, const double l_max, const size_t n_region, double &lambda, double &logl_H0)
{
	if (func_name!='R' && func_name!='L' && func_name!='r' && func_name!='l') {
		halt("func_name only takes 'R' or 'L': 'R' for log-restricted likelihood, 'L' for log-likelihood.");
		return;
	}

	wsUint n_cvt=UtW.row(), ni_test=UtW.col();
	wsUint n_index=(n_cvt+2+1)*(n_cvt+2)/2;

	cStdMatrix Uab_t(n_index, ni_test);
	cVector ab(n_index);
	CalcUab(UtW, Uty, Uab_t);

	//	if (e_mode!=0) {
	//		gsl_vector_set_zero (ab);
	//		Calcab (W, y, ab);
	//	}

	wsMat Pab			= sseMatrix(n_cvt+2, n_index);
	wsMat PPab			= sseMatrix(n_cvt+2, n_index);
	wsMat PPPab			= sseMatrix(n_cvt+2, n_index);
	wsMat Iab			= sseMatrix(n_cvt+2, n_index);
	FUNC_PARAM param0={ true, ni_test, n_cvt, eval, Uab_t, ab, 0, Pab, Iab, PPab, PPPab };

	CalcLambda(func_name, param0, l_min, l_max, n_region, lambda, logl_H0);

	sseUnmat(Pab, n_cvt+2);
	sseUnmat(PPab, n_cvt+2);
	sseUnmat(PPPab, n_cvt+2);
	sseUnmat(Iab, n_cvt+2);
}

typedef struct _xGemma
{
	mStrMat		*Xp_mat;
	mStrReal	*Xp_vec;
	mRcsvFun	Xm_funs;
	cVector		V_DLinv;
	cVector		*Vp_y;
	cVector		*Vp_y2;
	cVector		*Vp_D;
	cStdMatrix	*Mp_XtP;
	wsUint		N_silent;
	wsUint		N_x;		/* Size of X, including both case : have genotype or not */

	/* Updated for every loop */
	wsReal		R_dL;
	wsReal		R_delta;
	wsReal		R_yPy, R_yPPy, R_trP, R_lXPX;
	wsReal		R_xRx, R_xRy, R_yRy, R_LdetWHW;
} xGemma;

wsReal getRcsv(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	mRcsvFun	&Xm_funs	= X_d.Xm_funs;
	mStrMat		*Xp_mat		= X_d.Xp_mat;
	mStrReal	*Xp_vec		= X_d.Xp_vec;

	/* Print calling status if required */
	if (X_d.N_silent == 0) {
		for (wsUint I=0 ; I<N_level ; I++) LOGnf(" ");
		LOG("%s (%d,%d,%d)", S_lb, i, j, k);
	}

	/* If j == 0, call terminal function */
	if (j == 0xffffffff) {
		char S_idx[16];
		sprintf(S_idx, "%d,%d", i, k);
		string S_lbIK = S_lb;
		S_lbIK.append(S_idx);
		mStrReal_it X_find = Xp_vec->find(S_lbIK);

		if (X_find == Xp_vec->end()) {
			wsReal R_v0 = Xm_funs[S_lb].first(i, k, X_d);
			(*Xp_vec)[S_lbIK] = R_v0;
			if (X_d.N_silent == 0) LOGnf(" = %g (calc)\n", R_v0);
			return R_v0;
		} else {
			if (X_d.N_silent == 0) LOGnf(" = %g (cached)\n", X_find->second);
			return X_find->second;
		}
	}
	if (X_d.N_silent == 0) LOGnf("\n");

	/* Check is it exists */
	string S_lbs = (string)S_lb;
	mStrMat_it X_find = Xp_mat->find(S_lbs);
	if (X_find == Xp_mat->end())
		halt("SYSERR : Cascading matrix [%s] does not exists", S_lb);
	vMatrix *Xp_lbMat = X_find->second;
	if (Xp_lbMat->size() <= j) {
		/* Fill matrix */
		for (wsUint I=(wsUint)Xp_lbMat->size() ; I<=j ; I++)
			Xp_lbMat->push_back(new cStdMatrix(X_d.N_x, X_d.N_x,
				WISARD_NAN));
	}

	/* Check value is available */
	cMatrix	*Cp_mat	= (*Xp_lbMat)[j];
	wsReal **Ra_mat = Cp_mat->get();
	wsReal R_val = Ra_mat[i][k];

	/* Need to fill */
	if (R_val != R_val) {
		if (X_d.N_silent != 2) {
			for (wsUint I=0 ; I<N_level ; I++) LOGnf(" ");
			LOG("get(%s,%d,%d,%d)\n", S_lb, i, j, k);
		}
		/* Call recursive function */
		R_val = Ra_mat[i][k] =
			Xm_funs[S_lb].second(i, j, k, S_lb, N_level+1, X_d);
	}
	return R_val;
}

/*
 *
 * W by W part
 *
 */

wsReal ww0(wsUint i, wsUint k, xGemma &X_d)
{
	/* Both are NOT INTERCEPT */
	cVector V1 = X_d.V_DLinv * X_d.Mp_XtP->r2v(i);
	V1 *= X_d.Mp_XtP->r2v(k);
	wsReal R_ret = V1.sum();
	if (X_d.N_silent == 0) LOG("ww0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal wwN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal F_jk = getRcsv(i, j-1, k, S_lb, N_level, X_d);
	wsReal F_ji = getRcsv(i, j-1, j, S_lb, N_level, X_d);
	wsReal F_ki = getRcsv(k, j-1, j, S_lb, N_level, X_d);
	wsReal F_ii = getRcsv(j, j-1, j, S_lb, N_level, X_d);

	wsReal R_ret = F_jk - F_ji*F_ki / F_ii;
	if (X_d.N_silent == 0) LOG("wwN[%d,%d,%d] = %g\n", i, j, k, R_ret);
	return R_ret;
}

wsReal w2w0(wsUint i, wsUint k, xGemma &X_d)
{
	cVector V2 = X_d.V_DLinv.sq();
	/* Both are NOT INTERCEPT */
	cVector V1 = V2 * X_d.Mp_XtP->r2v(i);
	cVector Vk = X_d.Mp_XtP->r2v(k);
	V1 *= Vk;
	wsReal R_ret = V1.sum();
	if (X_d.N_silent == 0) LOG("w2w0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal w2wN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal G_jk = getRcsv(i, j-1, k, S_lb, N_level, X_d);

	wsReal F_ji = getRcsv(i, j-1, j, "wPw", N_level, X_d);
	wsReal F_ki = getRcsv(k, j-1, j, "wPw", N_level, X_d);
	wsReal F_ii = getRcsv(j, j-1, j, "wPw", N_level, X_d);

	wsReal G_ki = getRcsv(k, j-1, j, S_lb, N_level, X_d);
	wsReal G_ji = getRcsv(i, j-1, j, S_lb, N_level, X_d);
	wsReal G_ii = getRcsv(j, j-1, j, S_lb, N_level, X_d);

	wsReal R_ret = G_jk + F_ji*F_ki*G_ii/SQR(F_ii) - (F_ji*G_ki + F_ki*G_ji)/F_ii;
	if (X_d.N_silent == 0) LOG("w2wN[%d,%d,%d] = %g\n", i, j, k, R_ret);
	return R_ret;
}

wsReal w3w0(wsUint i, wsUint k, xGemma &X_d)
{
	cVector V3 = X_d.V_DLinv.cb();
	/* Both are NOT INTERCEPT */
	cVector V1 = V3 * X_d.Mp_XtP->r2v(i);
	V1 *= X_d.Mp_XtP->r2v(k);
	wsReal R_ret = V1.sum();
	if (X_d.N_silent == 0) LOG("w3w0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal w3wN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal H_jk = getRcsv(i, j-1, k, S_lb, N_level, X_d);

	wsReal F_ji = getRcsv(i, j-1, j, "wPw", N_level, X_d);
	wsReal F_ki = getRcsv(k, j-1, j, "wPw", N_level, X_d);
	wsReal G_ii = getRcsv(j, j-1, j, "wPPw", N_level, X_d);
	wsReal F_ii = getRcsv(j, j-1, j, "wPw", N_level, X_d);

	wsReal H_ki = getRcsv(k, j-1, j, "wPPPw", N_level, X_d);
	wsReal H_ji = getRcsv(i, j-1, j, "wPPPw", N_level, X_d);
	wsReal G_ji = getRcsv(i, j-1, j, "wPPw", N_level, X_d);
	wsReal G_ki = getRcsv(k, j-1, j, "wPPw", N_level, X_d);

	wsReal H_ii = getRcsv(j, j-1, j, "wPPPw", N_level, X_d);

	wsReal R_ret = H_jk - F_ji*F_ki*SQR(G_ii) / CUBE(F_ii)
		- (F_ji*H_ki + F_ki*H_ji + G_ji*G_ki) / F_ii
		+ (G_ii*(F_ji*G_ki + F_ki*G_ji) + F_ji*F_ki*H_ii) / SQR(F_ii);
	if (X_d.N_silent == 0) LOG("w3wN[%d,%d,%d] = %g\n", i, j, k, R_ret);
	return R_ret;
}

/*
 *
 * Y by W part
 *
 */

wsReal yw0(wsUint i, wsUint k, xGemma &X_d)
{
	/* Both are NOT INTERCEPT */
	cVector V1 = X_d.V_DLinv * *(X_d.Vp_y);
	cVector Vk = X_d.Mp_XtP->r2v(k);
	V1 *= Vk;
	wsReal R_ret = V1.sum();
	if (X_d.N_silent == 0) LOG("yw0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal ywN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal F_yk = getRcsv(i, j-1, k, S_lb, N_level, X_d);
	wsReal F_yi = getRcsv(i, j-1, j, S_lb, N_level, X_d);
	wsReal F_ki = getRcsv(k, j-1, j, "wPw", N_level, X_d);
	wsReal F_ii = getRcsv(j, j-1, j, "wPw", N_level, X_d);

	wsReal R_ret = F_yk - F_yi*F_ki / F_ii;
	if (X_d.N_silent == 0) LOG("ywN[%d,%d,%d] = %g\n", i, j, k, R_ret);
	return R_ret;
}

wsReal y2w0(wsUint i, wsUint k, xGemma &X_d)
{
	cVector V2 = X_d.V_DLinv.sq();
	/* Both are NOT INTERCEPT */
	cVector V1 = V2 * *(X_d.Vp_y);
	cVector Vk = X_d.Mp_XtP->r2v(k);
	V1 *= Vk;
	wsReal R_ret = V1.sum();
	if (X_d.N_silent == 0) LOG("y2w0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal y2wN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal G_yk = getRcsv(i, j-1, k, S_lb, N_level, X_d);

	wsReal F_yi = getRcsv(i, j-1, j, "yPw", N_level, X_d);
	wsReal F_ki = getRcsv(k, j-1, j, "wPw", N_level, X_d);
	wsReal F_ii = getRcsv(j, j-1, j, "wPw", N_level, X_d);

	wsReal G_ki = getRcsv(k, j-1, j, "wPPw", N_level, X_d);
	wsReal G_yi = getRcsv(i, j-1, j, S_lb, N_level, X_d);
	wsReal G_ii = getRcsv(j, j-1, j, "wPPw", N_level, X_d);

	wsReal R_ret = G_yk + F_yi*F_ki*G_ii/SQR(F_ii) - (F_yi*G_ki + F_ki*G_yi)/F_ii;
	if (X_d.N_silent == 0) LOG("y2wN[%d,%d,%d] = %g\n", i, j, k, R_ret);
	return R_ret;
}

wsReal y3w0(wsUint i, wsUint k, xGemma &X_d)
{
	cVector V3 = X_d.V_DLinv.cb();
	/* Both are NOT INTERCEPT */
	cVector V1 = V3 * *(X_d.Vp_y);
	cVector Xk = X_d.Mp_XtP->r2v(k);
	V1 *= Xk;

	wsReal R_ret = V1.sum();
	if (X_d.N_silent == 0) LOG("y3w0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal y3wN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal H_yk = getRcsv(i, j-1, k, S_lb, N_level, X_d);

	wsReal F_yi = getRcsv(i, j-1, j, "yPw", N_level, X_d);
	wsReal F_ki = getRcsv(k, j-1, j, "wPw", N_level, X_d);
	wsReal G_ii = getRcsv(j, j-1, j, "wPPw", N_level, X_d);
	wsReal F_ii = getRcsv(j, j-1, j, "wPw", N_level, X_d);

	wsReal H_ki = getRcsv(k, j-1, j, "wPPPw", N_level, X_d);
	wsReal H_yi = getRcsv(i, j-1, j, S_lb, N_level, X_d);
	wsReal G_yi = getRcsv(i, j-1, j, "yPPw", N_level, X_d);
	wsReal G_ki = getRcsv(k, j-1, j, "wPPw", N_level, X_d);

	wsReal H_ii = getRcsv(j, j-1, j, "wPPPw", N_level, X_d);

	wsReal R_ret = H_yk - F_yi*F_ki*SQR(G_ii) / CUBE(F_ii)
		- (F_yi*H_ki + F_ki*H_yi + G_yi*G_ki) / F_ii
		+ (G_ii*(F_yi*G_ki + F_ki*G_yi) + F_yi*F_ki*H_ii) / SQR(F_ii);
	if (X_d.N_silent == 0) LOG("y3wN[%d,%d,%d] = %g\n", i, j, k, R_ret);
	return R_ret;
}

/*
 *
 * Y by Y part
 *
 */

wsReal yy0(wsUint i, wsUint k, xGemma &X_d)
{
	cVector V1 = X_d.V_DLinv * *(X_d.Vp_y2);
	/* Both are ALWAYS not INTERCEPT */
	wsReal R_ret = V1.sum();
	if (X_d.N_silent == 0) LOG("yy0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal yyN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal F_yy = getRcsv(i, j-1, k, S_lb, N_level, X_d);
	wsReal F_yi = getRcsv(i, j-1, j, "yPw", N_level, X_d);
	wsReal F_ii = getRcsv(j, j-1, j, "wPw", N_level, X_d);

	wsReal R_ret = F_yy - SQR(F_yi) / F_ii;
	if (X_d.N_silent == 0) LOG("yyN[%d,%d,%d] = %g\n", i, j, k, R_ret);
	return R_ret;
}

wsReal y2y0(wsUint i, wsUint k, xGemma &X_d)
{
	cVector V2 = X_d.V_DLinv.sq();
	cVector V3 = V2 * *(X_d.Vp_y2);
	/* Both are NOT INTERCEPT */
	wsReal R_ret = V3.sum();
	if (X_d.N_silent == 0) LOG("y2y0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal y2yN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal G_yy = getRcsv(i, j-1, k, S_lb, N_level, X_d);

	wsReal F_yi = getRcsv(i, j-1, j, "yPw", N_level, X_d);
	wsReal G_ii = getRcsv(j, j-1, j, "wPPw", N_level, X_d);
	wsReal F_ii = getRcsv(j, j-1, j, "wPw", N_level, X_d);

	wsReal G_yi = getRcsv(k, j-1, j, "yPPw", N_level, X_d);

	wsReal R_ret = G_yy + SQR(F_yi)*G_ii/SQR(F_ii) - W2*F_yi*G_yi/F_ii;
	if (X_d.N_silent == 0) LOG("y2yN[%d,%d,%d] = %g\n", i, j, k, R_ret);
	return R_ret;
}

wsReal y3y0(wsUint i, wsUint k, xGemma &X_d)
{
	cVector V3 = X_d.V_DLinv.cb();
	cVector Vr = V3 * *(X_d.Vp_y2);
	/* Both are NOT INTERCEPT */
	wsReal R_ret = Vr.sum();
	if (X_d.N_silent == 0) LOG("y3y0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal y3yN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal H_yy = getRcsv(i, j-1, k, S_lb, N_level, X_d);

	wsReal F_yi = getRcsv(i, j-1, j, "yPw", N_level, X_d);
	wsReal G_ii = getRcsv(j, j-1, j, "wPPw", N_level, X_d);
	wsReal F_ii = getRcsv(j, j-1, j, "wPw", N_level, X_d);

	wsReal H_yi = getRcsv(i, j-1, j, "yPPPw", N_level, X_d);
	wsReal G_yi = getRcsv(i, j-1, j, "yPPw", N_level, X_d);

	wsReal H_ii = getRcsv(j, j-1, j, "wPPPw", N_level, X_d);

	wsReal R_ret = H_yy - SQR(F_yi)*SQR(G_ii) / CUBE(F_ii)
		- (R_2*F_yi*H_yi + SQR(G_yi)) / F_ii
		+ (R_2*F_yi*G_yi*G_ii + SQR(F_yi)*H_ii) / SQR(F_ii);
	if (X_d.N_silent == 0) LOG("y3yN[%d,%d,%d] = %g\n", i, j, k, R_ret);
	return R_ret;
}

/*
 *
 * tr(A) part
 *
 */

wsReal trP0(wsUint i, wsUint k, xGemma &X_d)
{
	/* Return SQsum of DLinv */
	wsReal R_ret = X_d.V_DLinv.sum();
	if (X_d.N_silent == 0) LOG("trP0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal trPN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal T_ii = getRcsv(i, j-1, k, S_lb, N_level, X_d);
	wsReal F_ii = getRcsv(i, j-1, k, "wPw", N_level, X_d);
	wsReal G_ii = getRcsv(i, j-1, k, "wPPw", N_level, X_d);

	return T_ii - G_ii / F_ii;
}

/*
 *
 * tr(A%*%A) part
 *
 */

wsReal trPP0(wsUint i, wsUint k, xGemma &X_d)
{
	/* Return SQsum of DLinv */
	wsReal R_ret = X_d.V_DLinv.ss();
	if (X_d.N_silent == 0) LOG("trPP0[%d,%d] = %g\n", i, k, R_ret);
	return R_ret;
}

wsReal trPPN(wsUint i, wsUint j, wsUint k, wsStrCst S_lb, wsUint N_level,
	xGemma &X_d)
{
	wsReal T_ii = getRcsv(i, j-1, k, S_lb, N_level, X_d);
	wsReal F_ii = getRcsv(i, j-1, k, "wPw", N_level, X_d);
	wsReal G_ii = getRcsv(i, j-1, k, "wPPw", N_level, X_d);
	wsReal H_ii = getRcsv(i, j-1, k, "wPPPw", N_level, X_d);

	wsReal R_GF = G_ii/F_ii;
	return T_ii - SQR(R_GF) - R_2*H_ii/F_ii;
}

/*
 *
 * cFemmaAnalysis definition
 *
 */

cFemmaAnalysis::cFemmaAnalysis(cIO *Cp_IO, cAnalysis *Cp_inpAnaPhi) :
	cAnalysis(Cp_IO)
{
	Mp_Wp_t = NULL;
	Mp_Wp = NULL;
	Vp_yp_t = NULL;
	//wsReal **Ra_ocv	= Cp_IO->getCovariates();
	//REAL_c *Ra_y	= Cp_IO->getPheno();
	wsUint N_oSamp	= Cp_IO->sizeSample();
	Cp_anaPhi		= Cp_inpAnaPhi;
	N_anaRow		= Cp_IO->sizeCovar()+2;

	/* Note : FEMMA does not supports --avail,--time */
	if (OPT_ENABLED(avail))
		LOG("[NOTE] FEMMA does not support --avail, it will be ignored\n");
	if (OPT_ENABLED(time))
		LOG("[NOTE] FEMMA does not support --time, it will be ignored\n");

	/* Get phenotype-covariates availability */
	Ba_isInc	= Cp_IO->getAvailPheno(&N_anaSamp, &Ra_yGS);
	LOG("[%d] samples used in GEMMA\n", N_anaSamp);

	/* Get reduced covariates matrix */
	Ra_DM		= getIO()->getFiltNullDesignMat(Ba_isInc, 0, N_anaSamp);

	/* Get phi matrix */
	wsReal **Ra_origPhi = getFullCorMat(Cp_anaPhi);
	Ra_phi = sseMsubsetRect(Ra_origPhi, N_oSamp, Ba_isInc, 1, N_anaSamp);

	/* Perform eigendecomposition */
#ifdef USE_ED2
	wsReal *Ra_D	= NULL;
	getEigenPhi(Ra_phi, N_anaSamp, &Ra_Pt, &Ra_D);
#else
	char B_stat = 0;
	wsReal *Ra_D	= eigenDecomp(Ra_phi, N_anaSamp, &B_stat, &Ra_Pt, 1);
#endif
	V_D.init(N_anaSamp, Ra_D);

	cStdMatrix M_U(N_anaSamp, N_anaSamp, Ra_Pt, MATDEL_NONE);
	cStdMatrix M_W(N_anaRow-1, N_anaSamp, Ra_DM, MATDEL_NONE);
	cVector V_y(Ra_yGS, N_anaSamp, 1);
	UtW = M_W.Mt(M_U);
	Uty = V_y.Mt(M_U);

	CalcLambda('L', V_D, UtW, Uty, 1e-5, 10e5, 10, l_mle_null, X_0.R_logL);
	LOG("Null lambda (%s) [%g], logL [%g]\n", OPT_ENABLED(ml) ? "ML" : "REML",
		X_0.R_lambda, X_0.R_logL);
	//		if (X_0.B_failed)
	//			halt_fmt(WISARD_FAIL_ESTIMATE_PARAM, "ratio of variance(lambda)");
}

cFemmaAnalysis::~cFemmaAnalysis()
{
	if (Ra_Pt) {
		sseUnmat(Ra_Pt, N_anaSamp);
		Ra_Pt = NULL;
	}
	if (Mp_Wp_t) delete Mp_Wp_t;
	if (Mp_Wp) delete Mp_Wp;
	if (Vp_yp_t) delete Vp_yp_t;
	sseFree(Ra_yGS);
	sseUnmat(Ra_phi, N_anaSamp);
	deallocMatrix(Ra_DM, N_anaRow-1, (void *)1);
	DEALLOC(Ba_isInc);
}

wsReal cFemmaAnalysis::NR(wsReal R_param, wsUint N_idxSNP, cStdMatrix &M_XtP,
	wsReal *Rp_tau, wsReal ***Rp_ab, cSymMatrix *Mp_WHWinv,
	wsReal *Rp_LdetH, wsReal *Rp_yPy, wsReal *Rp_LdetXHX/*=NULL*/)
{
	cDiagMatrix M_D(V_D);
	M_D.setDontDealloc();
	wsUint N_loop	= 0;
	wsReal R_lambda	= W1;
	wsReal R_sqLam	= W0;

	for ( ; N_loop<30 ; N_loop++) {
		/* Check lambda */
		if (R_lambda != R_lambda)
			return R_lambda;

		wsReal	R_trQ		= W0;
		wsReal	R_trQQ		= W0;

		/* Calculate D_\lambda^-1 = \lambdaD + I */
		cDiagMatrix&	M_D_Linv= M_D * R_lambda + 1;
		*Rp_LdetH = log(M_D_Linv.tr());
		M_D_Linv.inv();
		if (OPT_ENABLED(ml))
			R_trQ	= sseVsum(M_D_Linv.get()[0], N_anaSamp, &R_trQQ);

		/* Calculate (AA)^-1 = W_p'D_\lambda^-1W_p*/
		//		wsReal R_LdetWHW	= W0;
		cSymMatrix	M_AA	= Mp_Wp_t->MMt(M_D_Linv);
		*Mp_WHWinv			= M_AA.inv();//&R_LdetWHW);

		/* Calculate BB = W_p (AA)^-1 W_p' */
		cSymMatrix	M_BB	= Mp_Wp->MMt(*Mp_WHWinv);

		/* Calculate P'P_1P = D_\lambda^-1 - D_\lambda^-1(BB)D_\lambda^-1 */
		cSymMatrix	M_DBD	= M_D_Linv.MMt(M_BB);
		cSymMatrix	M_PP1P	= M_D_Linv - M_DBD;
//		wsReal R_trP1sP1s = sseTrSS(M_PP1P.get(), M_PP1P.row()); /* REML */
//		wsReal R_trP1s = M_PP1P.tr(); /* REML */
		M_DBD.rem();

/**/	cStdMatrix	M_PXt	= M_XtP.transpose();

		/* PP1P - PP1P %*% Xp %*% (Xp' %*% PP1P %*% Xp)^-1 %*% Xp' %*% PP1P */

		/* Calculate (Xp' %*% PP1P %*% Xp)^-1 */
//		wsReal R_LdetXPX	= W0;
/**/	cSymMatrix	M_XPX	= M_XtP.MMt(M_PP1P);
		M_XtP.rem();
/**/	cSymMatrix&	M_XPXi	= M_XPX.inv(Rp_LdetXHX);
		M_XPX.rem();
//		if (Rp_LdetXHX)
//			*Rp_LdetXHX = R_LdetWHW + R_LdetXPX;

/**/	cStdMatrix	M_PXp	= M_PP1P * M_PXt;
		M_PXt.rem();

		cSymMatrix	M_PPxPm	= M_PXp.MMt(M_XPXi);
		M_PXp.rem();
		//wsReal R_trPms = M_PPxPm.tr(); /* REML */
		cSymMatrix	M_PPxP	= M_PP1P - M_PPxPm;
		/* trPx, trPxPx - REML */
		if (!OPT_ENABLED(ml)) {
			R_trQ		= M_PPxP.tr();
			R_trQQ		= sseTrSS(M_PPxP.get(), M_PPxP.row());
		}

		/* Calculate K */
		cVector V_K		= *Vp_yp_t * M_PPxP;
//		wsReal *Ra_y	= Vp_yp_t->get();
//		wsReal **Ra_K	= sseMS(&Ra_y, 1, Vp_yp_t->size(),
//			M_PPxP.get(), M_PPxP.row());
		*Rp_yPy			= V_K.sum(Ra_yGS, N_anaSamp);
		wsReal R_yPPy	= V_K.ss();
		wsReal R_yPPPy	= V_K.qf(M_PPxP);
//		*Rp_yPy			= sseVV(Ra_K[0], N_samp, Ra_y);
//		wsReal R_yPPy	= sseVV(Ra_K[0], N_samp, Ra_K[0]);
//		wsReal R_yPPPy	= sseVpSpV(Ra_K[0], M_PPxP.get(), N_samp);

		/* */
		wsReal R_yPpPy		= (*Rp_yPy - R_yPPy) / R_lambda;
		wsReal R_yPpPpPy	= (*Rp_yPy + R_yPPPy - R_2*R_yPPy)
			/ R_sqLam;

		/* Calculate first derivate */
		wsReal R_trHinvPi	= (R_param - R_trQ) / R_lambda;
		wsReal R_d1			= REAL_CONST(-0.5) * R_trHinvPi +
			R_param/R_2 * R_yPpPy / *Rp_yPy;

		/* Set the value of tau */
		*Rp_tau = R_param / *Rp_yPy;

		/* Calculate second derivate */
		wsReal R_trQPQP	= (R_param + R_trQQ - R_2*R_trQ)
			/ R_sqLam;
		wsReal R_d2 = -R_param/R_2*
			(R_2*(R_yPpPpPy-*Rp_yPy)-SQR(R_yPpPy))/SQR(*Rp_yPy) +
			REAL_CONST(0.5) * R_trQPQP;

		wsReal R_dec = R_d1 / R_d2;
		/* Converged */
		if (fabs(R_dec) < 1e-8) {
			/* (W_p' %*% D_\lambda^-1 %*% W_p)^-1 %*% W_p' %*% D_\lambda^-1 %*% y_p */
			cStdMatrix	M_XPtDinv	= M_XtP * M_D_Linv;
			cStdMatrix	M_AH		= M_XPXi * M_XPtDinv;
			wsReal		*Ra_Vp		= Vp_yp_t->get();
			*Rp_ab					= sseMpMt(&Ra_Vp, 1, Vp_yp_t->size(),
				M_AH.get(), M_AH.row(), M_AH.col());

			break;
		}
		M_XPXi.rem();

		/* Go */
		R_lambda -= R_dec;
		R_sqLam = SQR(R_lambda);
	}

	return R_lambda;
}

double test(double R_lambda, void *Vp_dt)
{
	xGemma		*Xp_dt = (xGemma *)Vp_dt;
	mStrMat		Xm_mat;
	mStrReal	Xm_vec;
	vMatrix		X_wPw, X_yPw, X_yPy;
	vMatrix		X_wPPw, X_yPPw, X_yPPy;
	vMatrix		X_wPPPw, X_yPPPw, X_yPPPy, X_trP, X_trPP;

	/* Set map to NULL */
	Xm_mat.insert(make_pair((string)"wPw", &X_wPw));
	Xm_mat.insert(make_pair((string)"yPw", &X_yPw));
	Xm_mat.insert(make_pair((string)"yPy", &X_yPy));
	Xm_mat.insert(make_pair((string)"wPPw", &X_wPPw));
	Xm_mat.insert(make_pair((string)"yPPw", &X_yPPw));
	Xm_mat.insert(make_pair((string)"yPPy", &X_yPPy));
	Xm_mat.insert(make_pair((string)"wPPPw", &X_wPPPw));
	Xm_mat.insert(make_pair((string)"yPPPw", &X_yPPPw));
	Xm_mat.insert(make_pair((string)"yPPPy", &X_yPPPy));
	Xm_mat.insert(make_pair((string)"trP", &X_trP));
	Xm_mat.insert(make_pair((string)"trPP", &X_trPP));

	Xp_dt->V_DLinv.init(Xp_dt->Vp_D->size());
	wsReal*	Ra_DLinv	= Xp_dt->V_DLinv.get();
	wsReal*	Ra_D		= Xp_dt->Vp_D->get();
	wsUint	N_dim		= Xp_dt->Vp_D->size();
	wsUint	N_Med		= getMed(N_dim);
#ifdef USE_SSE
	sse_t	sse_1		= sseSet(1.0);
	sse_t	sse_L		= sseSet(R_lambda);
	for (wsUint i=0 ; i<N_Med ; i+=sseJmp) {
		sse_t *sse_D = (sse_t *)(Ra_D + i);
		sse_t *sse_V = (sse_t *)(Ra_DLinv + i);

		*sse_V = sseDiv(sse_1, sseAdd(sseMul(*sse_D, sse_L), sse_1));
	}
#endif
	for (wsUint i=N_Med ; i<N_dim ; i++)
		Ra_DLinv[i] = W1 / ((R_lambda*Ra_D[i]) + W1);

	/* Adapt to xGemma */
	Xp_dt->Xp_mat	= &Xm_mat;
	Xp_dt->Xp_vec	= &Xm_vec;

	/* Get yPy, yPPy, yPPPy */
	wsUint	k = Xp_dt->Mp_XtP->row()-1;
	wsReal	R_yPy	= getRcsv(k, k, k, "yPy", 0, *Xp_dt);
	wsReal	R_yRy	= getRcsv(k, k-1, k, "yPy", 0, *Xp_dt);
	wsReal	R_yPPy	= getRcsv(k, k, k, "yPPy", 0, *Xp_dt);
	wsReal	R_yPPPy	= getRcsv(k, k, k, "yPPPy", 0, *Xp_dt);
	wsReal	R_trP, R_trPP;
//	wsReal	R_xPx	= getRcsv(k+1, k, k+1, "wPw", 0, *Xp_dt);

	if (OPT_ENABLED(ml)) {
		R_trP	= Xp_dt->V_DLinv.sum();
		R_trPP	= Xp_dt->V_DLinv.ss();
	} else {
		R_trP	= getRcsv(k, k, k, "trP", 0, *Xp_dt);
		R_trPP	= getRcsv(k, k, k, "trPP", 0, *Xp_dt);
	}

	/* Calculate required terms for log-likelihood derivatives */
	wsReal R_yPpPy		= (R_yPy - R_yPPy) / R_lambda;
	wsReal R_yPpPpPy	= (R_yPy + R_yPPPy - (R_yPPy+R_yPPy)) / SQR(R_lambda);
	Xp_dt->R_xRx	= getRcsv(k, k-1, k, "wPw", 0, *Xp_dt);
	Xp_dt->R_xRy	= getRcsv(k, k-1, k, "yPw", 0, *Xp_dt);
	Xp_dt->R_yRy	= R_yRy;
	Xp_dt->R_LdetWHW = W0;
	for (wsUint i=0 ; i<=k ; i++)
		Xp_dt->R_LdetWHW += log(fabs(getRcsv(i, i-1, i, "wPw", 0, *Xp_dt)));

	/* Set m terms */
	wsReal R_m = OPT_ENABLED(ml) ? Xp_dt->Mp_XtP->col() :
		Xp_dt->Mp_XtP->col() - Xp_dt->Mp_XtP->row();

	wsReal R_trQp = (R_m - R_trP) / R_lambda;
	wsReal R_trQpQp = (R_m + R_trPP - (R_trP+R_trP)) / SQR(R_lambda);

	/* Calculate log-likelihood derivatives */
	wsReal R_ddL = (R_trQpQp - R_m*(R_2 * R_yPpPpPy * R_yPy
		- SQR(R_yPpPy))/SQR(R_yPy)) / R_2;
	wsReal R_dL = (-R_trQp + R_m*R_yPpPy/R_yPy) / R_2;

	/* Calculate delta... */
	Xp_dt->R_delta	= R_dL / R_ddL;
	Xp_dt->R_dL		= R_dL;
	Xp_dt->R_yPy	= R_yPy;
	Xp_dt->R_yPPy	= R_yPPy;
	Xp_dt->R_trP	= R_trP;
// 	if (!k) {
// 		Xp_dt->R_lXPX = log(fabs(Xp_dt->V_DLinv.sum()));
// 	} else {
// 		Xp_dt->R_lXPX	= W0;
// 		for (wsUint i=k ; i>0 ; i--) {
// 			wsReal Q = getRcsv(i, i-1, i, "wPw", 0, *Xp_dt);
// 			Xp_dt->R_lXPX += Q;
// 		}
//
// 		Xp_dt->R_lXPX *= fabs(Xp_dt->V_DLinv.sum());
// 		Xp_dt->R_lXPX = log(Xp_dt->R_lXPX);
// 	}

	/* Clear matrices */
	FOREACH (mStrMat_it, Xm_mat, i) {
		vMatrix &v = *(i->second);
		FOREACH (vMatrix_it, v, j)
			delete *j;
	}

	return R_dL;
}

typedef struct _xBrentLogL {
	wsReal		R_LdetH;
	cDiagMatrix	*Mp_D;
	wsUint		N_samp;
	cStdMatrix	*Mp_Wp_t;
	wsReal		R_param;
	wsReal		R_yPy;
	wsReal		R_LdetXHX;
	cStdMatrix	&M_XtP;
	wsReal		R_tau;
	wsReal		*Ra_yGS;
	cStdMatrix	*Mp_Wp;
	cVector		*Vp_yp_t;
	wsReal		**Ra_ab;
	cSymMatrix	*Mp_WHWinv;
} xBrentLogL;

#define R_2PI	(WISARD_PI * R_2)

void cFemmaAnalysis::_doTest(cStdMatrix &M_Xt, wsUint N_anaSamp,
	wsUint N_cova, wsMat Ra_Pt, cDiagMatrix *Mp_D, wsReal R_l0,
	int N_idxSNP, cVector *Vp_yp_t, cVector &V_y2p_t, xFemmaRes *Xp_ret)
{
	wsUint		N_add		= N_idxSNP==-1 ? 0 : 1;
	wsUint		N_cov		= N_cova + N_add;
	wsReal		R_param		= (wsReal)(N_anaSamp - (OPT_ENABLED(ml) ? 0 : N_cov));
	Xp_ret->R_lambda		= REAL_CONST(0.5);

	/* Set initial */
	cStdMatrix	M_P(N_anaSamp, N_anaSamp, Ra_Pt, MATDEL_NONE);
	wsReal		R_LdetXX = W0;

	cVector		V_x		= M_Xt.r2v_ptr(-1);
#ifdef USE_ED2
/**/cStdMatrix	M_XtP	= M_Xt.Mt(M_P);
#else
/**/cStdMatrix	M_XtP	= M_Xt * M_P;
#endif
	cSymMatrix	M_XtX	= M_Xt.Mt(); /* REML */
	//wsReal R_LdetXtX = W0;
	cSymMatrix	&M_XtXi	= M_XtX.inv(&R_LdetXX);
	M_XtX.rem();
	delete &M_XtXi;
	cVector		V_xP	= M_XtP.r2v_ptr(-1);

	/* Set cascading data */
	cVector		V_D(Mp_D->get()[0], Mp_D->row(), 1);
	xGemma		X_dt;
	if (OPT_ENABLED(verbose))
		X_dt.N_silent	= 0;
	else
		X_dt.N_silent	= 2;
	X_dt.N_x		= M_Xt.row();
	X_dt.Mp_XtP	= &M_XtP;
	X_dt.Vp_D		= &V_D;
	X_dt.Vp_y		= Vp_yp_t;
	X_dt.Vp_y2		= &V_y2p_t;
	/* Set function */
	X_dt.Xm_funs.insert(make_pair("wPw", xPairFun(ww0, wwN)));
	X_dt.Xm_funs.insert(make_pair("wPPw", xPairFun(w2w0, w2wN)));
	X_dt.Xm_funs.insert(make_pair("wPPPw", xPairFun(w3w0, w3wN)));
	X_dt.Xm_funs.insert(make_pair("yPw", xPairFun(yw0, ywN)));
	X_dt.Xm_funs.insert(make_pair("yPPw", xPairFun(y2w0, y2wN)));
	X_dt.Xm_funs.insert(make_pair("yPPPw", xPairFun(y3w0, y3wN)));
	X_dt.Xm_funs.insert(make_pair("yPy", xPairFun(yy0, yyN)));
	X_dt.Xm_funs.insert(make_pair("yPPy", xPairFun(y2y0, y2yN)));
	X_dt.Xm_funs.insert(make_pair("yPPPy", xPairFun(y3y0, y3yN)));
	X_dt.Xm_funs.insert(make_pair("trP", xPairFun(trP0, trPN)));
	X_dt.Xm_funs.insert(make_pair("trPP", xPairFun(trPP0, trPPN)));

	wsReal R_lBound = 10e-5;
	wsReal R_uBound = 10e5;

	/* Set initial lambda and perform NR */
	wsReal	R_lambda	= REAL_CONST(0.5);
	wsReal	R_dLprev	= WISARD_NAN;
	wsReal	R_tau		= WISARD_NAN;
	//cVector	V_est;
//	wsReal	R_est		= WISARD_NAN;
	for (wsUint i=0 ; i<30 ; i++) {
		/* Do loop */
		test(R_lambda, &X_dt);
		if (R_lBound > R_lambda || R_uBound < R_lambda) {
			Xp_ret->R_lambda = R_lambda;
			break;
		}

		/* Update lambda */
		R_lambda -= X_dt.R_delta;
		/* i != 0 then compare */
		if (i && fabs(R_dLprev-X_dt.R_dL)<1e-8) {
			wsReal R_m = (OPT_ENABLED(ml) ? M_Xt.col() : M_Xt.col()-M_Xt.row());
			/* Set tau */
			R_tau = R_m / X_dt.R_yPy;

			/* Get the estimates of coefficients */
			//wsReal		R_LdetXtPX = WISARD_NAN;
// 			cDiagMatrix	M_DLinv(Mp_D->row(), X_dt.V_DLinv.get(), 1);
// 			cSymMatrix  M_oo = M_XtP.MMt(M_DLinv);
// 			wsReal R_w2 = M_oo.detL();
// 			M_oo.inv(&R_w2);
// 			cSymMatrix	M_XtPX		= M_XtP.MMt(M_DLinv);
// 			cSymMatrix	&M_XtPXi	= M_XtPX.inv(&R_LdetXtPX);
// 			M_XtPX.rem();
//			cVector		V_XtPy		= M_XtP.MV(M_DLinv, *Vp_yp_t);
			//V_est = M_XtPXi * V_XtPy;
//			R_est = X_dt.R_xRy / X_dt.R_xRx;
//			delete &M_XtPXi;
//			V_XtPy.rem();

			/* Get eta term */
			wsReal R_eta = OPT_ENABLED(ml) ? W0 :
				- R_LdetXX + X_dt.R_LdetWHW;
			cDiagMatrix& M_nD = *Mp_D * R_lambda;
			M_nD += W1;

			/* Get final log-likelihood */
			Xp_ret->R_logL = (R_m * (log(R_m) - log(R_2PI) - W1 - log(X_dt.R_yPy))
				- M_nD.Lsum() + R_eta) / R_2;
			delete &M_nD;
			Xp_ret->R_lambda = R_lambda;
			break;
		}

		R_dLprev = X_dt.R_dL;
	}

	/* NR failed, then try Brent */
	if (NA(R_tau) || R_lBound > Xp_ret->R_lambda || R_uBound < Xp_ret->R_lambda) {
	/* Try Brent's method to get an estimate */
// 	xBrentLogL X_bll = {
// 		W0,
// 		Mp_D,
// 		N_samp, Mp_Wp_t, R_param,
// 		W0,
// 		W0,
// 		M_XtP,
// 		W0,
// 		Ra_yGS,
// 		Mp_Wp,
// 		Vp_yp_t,
// 		NULL,
// 		NULL
// 	};
 		wsReal R_lambda = Brent_fmin(R_lBound, R_uBound, test, &X_dt, 1e-5);
		wsReal R_m = (OPT_ENABLED(ml) ? M_Xt.col() : M_Xt.col()-M_Xt.row());
		/* Set tau */
		R_tau = R_m / X_dt.R_yPy;

		/* Get the estimates of coefficients */
//		wsReal		R_LdetXtPX = WISARD_NAN;
//		cDiagMatrix	M_DLinv(Mp_D->row(), X_dt.V_DLinv.get(), 1);
// 		cSymMatrix  M_oo = M_XtP.MMt(M_DLinv);
// 		wsReal R_w2 = M_oo.detL();
// 		M_oo.inv(&R_w2);
// 		cSymMatrix	M_XtPX		= M_XtP.MMt(M_DLinv);
// 		cSymMatrix	&M_XtPXi	= M_XtPX.inv(&R_LdetXtPX);
// 		M_XtPX.rem();
//		cVector		V_XtPy		= M_XtP.MV(M_DLinv, *Vp_yp_t);
		//V_est					= M_XtPXi * V_XtPy;
//		R_est = X_dt.R_xRy / X_dt.R_xRx;

//		delete &M_XtPXi;
//		V_XtPy.rem();

		/* Get eta term */
		cDiagMatrix& M_nD = *Mp_D * R_lambda;
		M_nD += W1;
// 		if (!OPT_ENABLED(ml)) {
// 			LOG("From %g->", X_dt.R_LdetWHW);
// 			cDiagMatrix& M_nDinv = M_nD.inv();
// 			cSymMatrix R = M_XtP.MMt(M_nDinv);
// 			cSymMatrix& Rinv = R.inv(&X_dt.R_LdetWHW);
// 			delete &Rinv;
// 			LOGnf("%g\n", X_dt.R_LdetWHW);
// 		}
		wsReal R_eta = OPT_ENABLED(ml) ? W0 :
			R_LdetXX - X_dt.R_LdetWHW;

		/* Get final log-likelihood */
		Xp_ret->R_logL = (R_m * (log(R_m) - log(R_2PI) - W1 - log(X_dt.R_yPy))
			- M_nD.Lsum() + R_eta) / R_2;
		delete &M_nD;
		Xp_ret->R_lambda = R_lambda;
	}

	/* Failed to get lambda */
	if (NA(Xp_ret->R_lambda)) {
		Xp_ret->B_failed = 1;
		return;
	}

	wsReal R_m = (OPT_ENABLED(ml) ? M_Xt.col() : M_Xt.col()-M_Xt.row());
//	Xp_ret->R_logL = X_dt.R_dL;

	/* Do test only if N_idxSNP is not -1 */
	if (N_idxSNP >= 0) {
		//wsReal *Ra_est = V_est.get();
		cDiagMatrix	M_DLinv(Mp_D->row(), X_dt.V_DLinv.get(), 1);

//		wsReal	R_varBeta	= Vp_yp_t->qf(M_DLinv)/(R_m * V_xP.qf(M_DLinv));
//		wsReal	R_varBeta	= Vp_yp_t->qf(M_DLinv)/(R_m * X_dt.R_xRx);


		/* Perform Wald test */
//		Xp_ret->R_Twald = SQR(Ra_est[M_Xt.row()-1]) / R_varBeta;
		Xp_ret->R_Twald = (X_dt.R_yRy-X_dt.R_yPy)*R_m/X_dt.R_yPy;//SQR(R_est) / R_varBeta;	
		Xp_ret->R_Pwald = pf(Xp_ret->R_Twald, 1.0, R_param, 0);

		/* Perform LRT */
		Xp_ret->R_Tlrt	= W2 * (Xp_ret->R_logL - R_l0);
		Xp_ret->R_Plrt	= PVchisq(Xp_ret->R_Tlrt, 1.0);

		/* Perform score test */

		Xp_ret->R_Tscore	= (wsReal)M_Xt.col() * SQR(X_dt.R_xRy) / (X_dt.R_xRx * X_dt.R_yRy);
		Xp_ret->R_Pscore	=  pf(Xp_ret->R_Tscore, 1.0, R_param, 0);
// 		cVector		V_yHat		= Mp_Wp_t->tV(Ra_est, M_Xt.row()-1);
// 		cVector		V_e			= *Vp_yp_t - V_yHat;
// 		wsReal		R_eVx		= V_e.qf(M_DLinv, V_x);
//
// 		cSymMatrix	M_WHW		= Mp_Wp_t->MMt(M_DLinv);
// 		cSymMatrix	&M_WHWinv	= M_WHW.inv();
// 		M_WHW.rem();
// 		cVector		V_xV		= V_x * M_DLinv;
// 		cVector		V_xVZ		= V_xV * *Mp_Wp;
// 		wsReal		R_denom		= V_x.qf(M_DLinv) - V_xVZ.qf(M_WHWinv);
// 		delete &M_WHWinv;
//
// 		Xp_ret->R_Tscore = SQR(R_eVx) / R_denom;
// 		Xp_ret->R_Pscore = chiprobP(Xp_ret->R_Tscore, 1.0);
	}
	M_XtP.rem();

	Xp_ret->B_failed = 0;
}

int femma(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	if (N_idx < 0 || N_idx >= OPT_NUMBER(thread))
		halt("Error index [%d]", N_idx);
	char			S_bufAnno[512];
	xAnaThread*		Xp_at		= (xAnaThread *)Vp_shareData;
	cFemmaAnalysis*	Cp_fa		= (cFemmaAnalysis *)(Xp_at->Vp_data);
	xFemmaRes&		Xp_0		= Cp_fa->get0();
	xFemmaRes*		Xa_res		= (xFemmaRes *)Vp_result;
	int*			Np_data		= (int *)Vp_data;
	cIO*			Cp_IO		= Cp_fa->getIO();
	vVariant&		Xa_snp		= Cp_IO->getVariant();
	xMaf*			Xp_maf		= Cp_IO->getMAF();
	wsMat			Ra_DM		= Cp_fa->get0DM();
	wsUint			N_anaSamp	= Cp_fa->getAnaSampSize();
	wsUint			N_anaRow	= Cp_IO->sizeCovar()+2;
	const char		*Ba_isInc	= Cp_fa->getIsInc();

	wsUint			N_s			= (wsUint)Np_data[0];
	wsUint			N_e			= (wsUint)Np_data[1];
	wsMat			Ra_Pt		= Cp_fa->getPt();
	cVector&		V_D			= Cp_fa->getD();
	cStdMatrix&		UtW			= Cp_fa->getUtW();
	cVector&		Uty			= Cp_fa->getUty();
	wsRealCst		l_mle_null	= Cp_fa->getLmleNull();

	wsMat			Ra_fDM		= sseMatrix(N_anaRow, N_anaSamp);
	for (wsUint k=0 ; k<(N_anaRow-1) ; k++)
		memcpy(Ra_fDM[k], Ra_DM[k], sizeof(wsReal)*N_anaSamp);

	wsUint n_cvt=UtW.row();
	wsUint n_index=(n_cvt+2+1)*(n_cvt+2)/2;
	cVector ab(n_index);
	cStdMatrix Uab_t(n_index, N_anaSamp);
	CalcUab(UtW, Uty, Uab_t);
	wsMat Pab			= sseMatrix(n_cvt+2, n_index);
	wsMat PPab			= sseMatrix(n_cvt+2, n_index);
	wsMat PPPab			= sseMatrix(n_cvt+2, n_index);
	wsMat Iab			= sseMatrix(n_cvt+2, n_index);

	for (wsUint Q=N_s ; Q<N_e ; Q++) {
		cTimer		t;
		xFemmaRes	*Xp_r	= Xa_res + Q;
		Xp_r->B_failed = 1;
		/* Calculate X_p = P^t %*% X ==> X_p' = X' %*% P */
		if (Xp_maf[Q].R_allMaf > W0) {
			//setGeno(Cp_IO, Ba_isInc);
			/* Fill X */
			setGeno(Cp_IO, Q, Ra_fDM[N_anaRow-1], Ba_isInc);
			// 			for (wsUint j=0,J=0 ; j<N_samp ; j++) {
			// 				if (!Ba_isInc[j]) continue;
			// 				char N_geno = Na_data[j][Q];
			// 				if (isMissing(N_geno))
			// 					Ra_fDM[N_anaRow-1][J++] = X_snp.maf * W2;
			// 				else
			// 					Ra_fDM[N_anaRow-1][J++] = (wsReal)N_geno;
			// 			}
			t.start();
			wsVec _Utx = sseMpV(Ra_Pt, N_anaSamp, N_anaSamp, Ra_fDM[N_anaRow-1]);
			cVector Utx(_Utx, N_anaSamp);
			CalcUab(UtW, Uty, Utx, Uab_t);

			FUNC_PARAM param1={ false, N_anaSamp, n_cvt, V_D, Uab_t, ab, 0, Pab, Iab, PPab, PPPab };

			//3 is before 1, for beta
			double beta, se;
			CalcRLScore(l_mle_null, param1, beta, se, Xp_r->R_Pscore, Xp_r->R_Tscore);

			double lambda_mle, lambda_remle, logl_H1;
			CalcLambda('R', param1, 1e-5, 10e5, 10, lambda_remle, logl_H1);
			CalcRLWald(lambda_remle, param1, beta, se, Xp_r->R_Pwald, Xp_r->R_Twald);

			CalcLambda('L', param1, 1e-5, 10e5, 10, lambda_mle, logl_H1);
			Xp_r->R_Tlrt = 2.0*(logl_H1-Xp_0.R_logL);
			Xp_r->R_Plrt=PVchisq(Xp_r->R_Tlrt, W1);
			Xp_r->B_failed = 0;
			Xp_r->R_lambda = lambda_remle;
			Xp_r->R_time = t.get();

			if (IS_ASSIGNED(annogene))
				strcpy(S_bufAnno, Xa_snp[Q].anno);
		}
		Np_data[2]++;

		/* Print progress */
		if (N_idx==0 && (Q%100)==0) {
			wsUint N_sum = 0;
			for (int x=2 ; x<OPT_NUMBER(thread) ; x+=3)
				N_sum += Np_data[x];
			notice("%d/%d variants tested...\r", N_sum, Xa_snp.size());
		}
	}
	sseUnmat(Ra_fDM, N_anaRow);
	sseUnmat(Pab, n_cvt+2);
	sseUnmat(PPab, n_cvt+2);
	sseUnmat(PPPab, n_cvt+2);
	sseUnmat(Iab, n_cvt+2);
	return 0;
}

void cFemmaAnalysis::run()
{
	vVariant&	Xa_snp	= Cp_IO->getVariant();
	xMaf*		Xp_maf	= Cp_IO->getMAF();
	/**/cExporter*	Cp_gma	= cExporter::summon("gemma.res");

	/* Export header */
	char		S_time[512]		={ 0, };
	if (OPT_ENABLED(time))
		strcpy(S_time, "	TIME");
	headerVariant(Cp_gma);
	Cp_gma->fmt("%s	Lambda	STAT_WALD	P_WALD	STAT_LRT	P_LRT	STAT_SCORE	P_SCORE\n",
		S_time);

	/* For all SNPs */
	wsUint	i			= 0;
	// 	char	**Na_data	= getIO()->getSNPData();
	// 	wsUint	N_samp		= getIO()->getSampleSize();

	if (OPT_NUMBER(thread) == 1) {
		wsMat	Ra_fDM	= sseMatrix(N_anaRow, N_anaSamp);
		for (wsUint j=0 ; j<(N_anaRow-1) ; j++)
			memcpy(Ra_fDM[j], Ra_DM[j], sizeof(wsReal)*N_anaSamp);

		wsUint n_cvt=UtW.row();
		wsUint n_index=(n_cvt+2+1)*(n_cvt+2)/2;
		cVector ab(n_index);
		cStdMatrix Uab_t(n_index, N_anaSamp);
		CalcUab(UtW, Uty, Uab_t);

		wsMat Pab			= sseMatrix(n_cvt+2, n_index);
		wsMat PPab			= sseMatrix(n_cvt+2, n_index);
		wsMat PPPab			= sseMatrix(n_cvt+2, n_index);
		wsMat Iab			= sseMatrix(n_cvt+2, n_index);

		FOREACHDO(vVariant_it, Xa_snp, X_snp, i++) {
			/* Set to initially 'FAILED' */
			xFemmaRes	Xr;
			Xr.B_failed = 1;

			/* Test only maf > 0 */
			if (Xp_maf[i].R_allMaf > W0) {
				/* Calculate X_p = P^t %*% X ==> X_p' = X' %*% P */

				/* Fill X */
				setGeno(Cp_IO, i, Ra_fDM[N_anaRow-1], Ba_isInc);
				/* Calc */
				wsVec _Utx = sseMpV(Ra_Pt, N_anaSamp, N_anaSamp, Ra_fDM[N_anaRow-1]);
				cVector Utx(_Utx, N_anaSamp);
				CalcUab(UtW, Uty, Utx, Uab_t);

				FUNC_PARAM param1={ false, N_anaSamp, n_cvt, V_D, Uab_t, ab, 0, Pab, Iab, PPab, PPPab };

				//3 is before 1, for beta
				double beta, se;
				CalcRLScore(l_mle_null, param1, beta, se, Xr.R_Pscore, Xr.R_Tscore);

				double lambda_mle, lambda_remle, logl_H1;
				CalcLambda('R', param1, 1e-5, 10e5, 10, lambda_remle, logl_H1);
				CalcRLWald(lambda_remle, param1, beta, se, Xr.R_Pwald, Xr.R_Twald);

				CalcLambda('L', param1, 1e-5, 10e5, 10, lambda_mle, logl_H1);
				Xr.R_Tlrt	= 2.0*(logl_H1-X_0.R_logL);
				Xr.R_Plrt	= PVchisq(Xr.R_Tlrt, W1);

				Xr.B_failed	= 0;
				Xr.R_lambda	= lambda_remle;
			}

			if (OPT_ENABLED(time))
				sprintf(S_time, "	%s", cTimer::fmt(Xr.R_time));

			/* --remna */
			if (((NA(Xr.R_Pwald) && NA(Xr.R_Tlrt) && NA(Xr.R_Pscore))
				|| Xr.B_failed) && OPT_ENABLED(remna))
				continue;

			/* Print result */
			entryVariant(Cp_gma, *X_snp);
			if (Xr.B_failed)
				Cp_gma->fmt("%s	NA	NA	NA	NA	NA	NA	NA\n", S_time);
			else
				Cp_gma->fmt("%s	%g	%g	%g	%g	%g	%g	%g\n", S_time,
				Xr.R_lambda,
				Xr.R_Twald, Xr.R_Pwald, Xr.R_Tlrt, Xr.R_Plrt,
				Xr.R_Tscore, Xr.R_Pscore);

			if ((i%100) == 0)
				notice("%d/%d variants tested...\r", i, Xa_snp.size());
		}

		sseUnmat(Ra_fDM, N_anaRow);
		sseUnmat(Pab, n_cvt+2);
		sseUnmat(PPab, n_cvt+2);
		sseUnmat(PPPab, n_cvt+2);
		sseUnmat(Iab, n_cvt+2);
	}
	else {
		/* Prepare full DMs */
		// 		MULTI_MALLOC(Ra_fDMs, MAT_t, OPT_NUMBER(thread));
		// 		for (int j=0 ; j<OPT_NUMBER(thread) ; j++) {
		// 			Ra_fDMs[j] = sseMatrix(N_anaRow, N_anaSamp);
		// 			for (wsUint k=0 ; k<(N_anaRow-1) ; k++)
		// 				memcpy(Ra_fDMs[j][k], Ra_DM[k], sizeof(wsReal)*N_anaSamp);
		// 		}
		/* Prepare thread info */
		xAnaThread X_th ={
			getIO(),
			this
		};
		xFemmaRes *Xa_r = NULL;
		wsAlloc(Xa_r, xFemmaRes, Xa_snp.size());
		WORKER().run(femma, forVariant_equal, &X_th, Xa_r, sizeof(int)*3);

		FOREACHDO(vVariant_it, Xa_snp, X_snp, i++) {
			xFemmaRes &Xr = Xa_r[i];

			if (OPT_ENABLED(time))
				sprintf(S_time, "	%s", cTimer::fmt(Xr.R_time));

			/* --remna */
			if (((NA(Xr.R_Pwald) && NA(Xr.R_Tlrt) && NA(Xr.R_Pscore))
				|| Xr.B_failed) && OPT_ENABLED(remna))
				continue;

			entryVariant(Cp_gma, *X_snp);
			if (Xr.B_failed)
				Cp_gma->fmt("%s	NA	NA	NA	NA	NA	NA	NA\n", S_time);
			else
				Cp_gma->fmt("%s	%g	%g	%g	%g	%g	%g	%g\n", S_time,
				Xr.R_lambda,
				Xr.R_Twald, Xr.R_Pwald, Xr.R_Tlrt, Xr.R_Plrt,
				Xr.R_Tscore, Xr.R_Pscore);
		}

		/* Dealloc DMs */
		// 		for (int j=0 ; j<OPT_NUMBER(thread) ; j++)
		// 			deallocMatrix(Ra_fDMs[j], N_anaRow, (void *)1);
		// 		DEALLOC(Ra_fDMs);
		DEALLOC(Xa_r);
	}

	LOG("%d/%d variants tested...\n", Xa_snp.size(), Xa_snp.size());
	delete Cp_gma;

	if (Ra_Pt) {
		sseUnmat(Ra_Pt, N_anaSamp);
		Ra_Pt = NULL;
	}
}

#endif

} // End namespace ONETOOL
