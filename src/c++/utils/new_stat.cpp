#include <math.h>
#include "global/common.h"
#include "utils/util.h"
#include "utils/new_stat.h"

namespace ONETOOL {

#define M_LN_SQRT_2PI 0.918938533204672741780329736406 

	double	x_alngam(double *x)
		/*
		**********************************************************************

		double alngam(double *x)
		double precision LN of the GAMma function


		Function


		Returns the natural logarithm of GAMMA(X).


		Arguments


		X --> value at which scaled log gamma is to be returned
		X is DOUBLE PRECISION


		Method


		If X .le. 6.0, then use recursion to get X below 3
		then apply rational approximation number 5236 of
		Hart et al, Computer Approximations, John Wiley and
		Sons, NY, 1968.

		If X .gt. 6.0, then use recursion to get X to at least 12 and
		then use formula 5423 of the same source.

		**********************************************************************
		*/
	{
#define hln2pi 0.91893853320467274178e0
		static double coef[5] = {
			0.83333333333333023564e-1,-0.27777777768818808e-2,0.79365006754279e-3,
			-0.594997310889e-3,0.8065880899e-3
		};
		static double scoefd[4] = {
			0.62003838007126989331e2,0.9822521104713994894e1,-0.8906016659497461257e1,
			0.1000000000000000000e1
		};
		static double scoefn[9] = {
			0.62003838007127258804e2,0.36036772530024836321e2,0.20782472531792126786e2,
			0.6338067999387272343e1,0.215994312846059073e1,0.3980671310203570498e0,
			0.1093115956710439502e0,0.92381945590275995e-2,0.29737866448101651e-2
		};
		static int K1 = 9;
		static int K3 = 4;
		static int K5 = 5;
		double _alngam,offset,prod,xx;
		int i,n;
		/*
		..
		.. Executable Statements ..
		*/
		if(!(*x <= 6.0e0)) goto S70;
		prod = 1.0e0;
		xx = *x;
		if(!(*x > 3.0e0)) goto S30;
	S10:
		if(!(xx > 3.0e0)) goto S20;
		xx -= 1.0e0;
		prod *= xx;
		goto S10;
	S30:
	S20:
		if(!(*x < 2.0e0)) goto S60;
	S40:
		if(!(xx < 2.0e0)) goto S50;
		prod /= xx;
		xx += 1.0e0;
		goto S40;
	S60:
	S50:
		_alngam = x_dblEvalPolyX(scoefn,9,xx-2.0)/x_dblEvalPolyX(scoefd,4,xx-2.0);
		/*
		COMPUTE RATIONAL APPROXIMATION TO GAMMA(X)
		*/
		_alngam *= prod;
		_alngam = log(_alngam);
		goto S110;
	S70:
		offset = hln2pi;
		/*
		IF NECESSARY MAKE X AT LEAST 12 AND CARRY CORRECTION IN OFFSET
		*/
		n = x_fifidint(12.0e0-*x);
		if(!(n > 0)) goto S90;
		prod = 1.0e0;
		for(i=1; i<=n; i++) prod *= (*x+(double)(i-1));
		offset -= log(prod);
		xx = *x+(double)n;
		goto S100;
	S90:
		xx = *x;
	S100:
		/*
		COMPUTE POWER SERIES
		*/
		_alngam = x_dblEvalPolyX(coef,5,1.0/SQR(xx))/xx;
		_alngam += (offset+(xx-0.5e0)*log(xx)-xx);
	S110:
		return _alngam;
#undef hln2pi
	}
	void	x_bratio(double a,double b,double x,double y,double *w,
		double *w1,int *ierr)
		/*
		-----------------------------------------------------------------------

		EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B)

		--------------------

		IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1
		AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES

		W  = IX(A,B)
		W1 = 1 - IX(A,B)

		IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
		IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
		W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
		THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
		ONE OF THE FOLLOWING VALUES ...

		IERR = 1  IF A OR B IS NEGATIVE
		IERR = 2  IF A = B = 0
		IERR = 3  IF X .LT. 0 OR X .GT. 1
		IERR = 4  IF Y .LT. 0 OR Y .GT. 1
		IERR = 5  IF X + Y .NE. 1
		IERR = 6  IF X = A = 0
		IERR = 7  IF Y = B = 0

		--------------------
		WRITTEN BY ALFRED H. MORRIS, JR.
		NAVAL SURFACE WARFARE CENTER
		DAHLGREN, VIRGINIA
		REVISED ... NOV 1991
		-----------------------------------------------------------------------
		*/
	{
		double a0,b0,eps,lambda,t,x0,y0,z;
		int ierr1,ind,n;
		double T2,T3,T4,T5;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
		FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
		*/
		eps = x_spmpar(1);
		*w = *w1 = 0.0;
		if(a < 0.0 || b < 0.0) goto S270;
		if(a == 0.0 && b == 0.0) goto S280;
		if(x < 0.0 || x > 1.0) goto S290;
		if(y < 0.0 || y > 1.0) goto S300;
		z = x+y-1.0;
		if(fabs(z) > 3.0*eps) goto S310;
		*ierr = 0;
		if(x == 0.0) goto S210;
		if(y == 0.0) goto S230;
		if(a == 0.0) goto S240;
		if(b == 0.0) goto S220;
		eps = x_fifdmax1(eps,1.e-15);
		if(x_fifdmax1(a,b) < 1.e-3*eps) goto S260;
		ind = 0;
		a0 = a;
		b0 = b;
		x0 = x;
		y0 = y;
		if(x_fifdmin1(a0,b0) > 1.0) goto S40;
		/*
		PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
		*/
		if(x <= 0.5) goto S10;
		ind = 1;
		a0 = b;
		b0 = a;
		x0 = y;
		y0 = x;
	S10:
		if(b0 < x_fifdmin1(eps,eps* a0)) goto S90;
		if(a0 < x_fifdmin1(eps,eps* b0) && b0* x0 <= 1.0) goto S100;
		if(x_fifdmax1(a0,b0) > 1.0) goto S20;
		if(a0 >= x_fifdmin1(0.2,b0)) goto S110;
		if(pow(x0,a0) <= 0.9) goto S110;
		if(x0 >= 0.3) goto S120;
		n = 20;
		goto S140;
	S20:
		if(b0 <= 1.0) goto S110;
		if(x0 >= 0.3) goto S120;
		if(x0 >= 0.1) goto S30;
		if(pow(x0* b0,a0) <= 0.7) goto S110;
	S30:
		if(b0 > 15.0) goto S150;
		n = 20;
		goto S140;
	S40:
		/*
		PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
		*/
		if(a > b) goto S50;
		lambda = a-(a+b)*x;
		goto S60;
	S50:
		lambda = (a+b)*y-b;
	S60:
		if(lambda >= 0.0) goto S70;
		ind = 1;
		a0 = b;
		b0 = a;
		x0 = y;
		y0 = x;
		lambda = fabs(lambda);
	S70:
		if(b0 < 40.0 && b0* x0 <= 0.7) goto S110;
		if(b0 < 40.0) goto S160;
		if(a0 > b0) goto S80;
		if(a0 <= 100.0) goto S130;
		if(lambda > 0.03* a0) goto S130;
		goto S200;
	S80:
		if(b0 <= 100.0) goto S130;
		if(lambda > 0.03* b0) goto S130;
		goto S200;
	S90:
		/*
		EVALUATION OF THE APPROPRIATE ALGORITHM
		*/
		*w = x_fpser(a0,b0,x0,eps);
		*w1 = 0.5+(0.5-*w);
		goto S250;
	S100:
		*w1 = x_apser(a0,b0,x0,eps);
		*w = 0.5+(0.5-*w1);
		goto S250;
	S110:
		*w = x_bpser(a0,b0,x0,eps);
		*w1 = 0.5+(0.5-*w);
		goto S250;
	S120:
		*w1 = x_bpser(b0,a0,y0,eps);
		*w = 0.5+(0.5-*w1);
		goto S250;
	S130:
		T2 = 15.0*eps;
		*w = x_bfrac(&a0,&b0,&x0,&y0,&lambda,&T2);
		*w1 = 0.5+(0.5-*w);
		goto S250;
	S140:
		*w1 = x_bup(&b0,&a0,&y0,&x0,&n,&eps);
		b0 += (double)n;
	S150:
		T3 = 15.0*eps;
		x_bgrat(&b0,&a0,&y0,&x0,w1,&T3,&ierr1);
		*w = 0.5+(0.5-*w1);
		goto S250;
	S160:
		n = (int)b0;
		b0 -= (double)n;
		if(b0 != 0.0) goto S170;
		n -= 1;
		b0 = 1.0;
	S170:
		*w = x_bup(&b0,&a0,&y0,&x0,&n,&eps);
		if(x0 > 0.7) goto S180;
		*w += x_bpser(a0,b0,x0,eps);
		*w1 = 0.5+(0.5-*w);
		goto S250;
	S180:
		if(a0 > 15.0) goto S190;
		n = 20;
		*w += x_bup(&a0,&b0,&x0,&y0,&n,&eps);
		a0 += (double)n;
	S190:
		T4 = 15.0*eps;
		x_bgrat(&a0,&b0,&x0,&y0,w,&T4,&ierr1);
		*w1 = 0.5+(0.5-*w);
		goto S250;
	S200:
		T5 = 100.0*eps;
		*w = x_basym(&a0,&b0,&lambda,&T5);
		*w1 = 0.5+(0.5-*w);
		goto S250;
	S210:
		/*
		TERMINATION OF THE PROCEDURE
		*/
		if(a == 0.0) goto S320;
	S220:
		*w = 0.0;
		*w1 = 1.0;
		return;
	S230:
		if(b == 0.0) goto S330;
	S240:
		*w = 1.0;
		*w1 = 0.0;
		return;
	S250:
		if(ind == 0) return;
		t = *w;
		*w = *w1;
		*w1 = t;
		return;
	S260:
		/*
		PROCEDURE FOR A AND B .LT. 1.E-3*EPS
		*/
		*w = b/(a+b);
		*w1 = a/(a+b);
		return;
	S270:
		/*
		ERROR RETURN
		*/
		*ierr = 1;
		return;
	S280:
		*ierr = 2;
		return;
	S290:
		*ierr = 3;
		return;
	S300:
		*ierr = 4;
		return;
	S310:
		*ierr = 5;
		return;
	S320:
		*ierr = 6;
		return;
	S330:
		*ierr = 7;
		return;
	}
	/*
	-----------------------------------------------------------------------

	COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B .GE. 8

	--------

	IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
	LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).

	-----------------------------------------------------------------------
	*/
	double	x_algdiv(double a,double b)
	{
		static double c0 = .833333333333333e-01;
		static double c1 = -.277777777760991e-02;
		static double c2 = .793650666825390e-03;
		static double c3 = -.595202931351870e-03;
		static double c4 = .837308034031215e-03;
		static double c5 = -.165322962780713e-02;
		double c,d,h,s11,s3,s5,s7,s9,t,u,v,w,x,x2;
		/*
		..
		.. Executable Statements ..
		*/
		if(a <= b) goto S10;
		h = b/ a;
		c = 1.0/(1.0+h);
		x = h/(1.0+h);
		d = a+(b-0.5);
		goto S20;
	S10:
		h = a/ b;
		c = h/(1.0+h);
		x = 1.0/(1.0+h);
		d = b+(a-0.5);
	S20:
		/*
		SET SN = (1 - X**N)/(1 - X)
		*/
		x2 = x*x;
		s3 = 1.0+(x+x2);
		s5 = 1.0+(x+x2*s3);
		s7 = 1.0+(x+x2*s5);
		s9 = 1.0+(x+x2*s7);
		s11 = 1.0+(x+x2*s9);
		/*
		SET W = DEL(B) - DEL(A + B)
		*/
		t = SQR(1.0/ b);
		w = ((((c5*s11*t+c4*s9)*t+c3*s7)*t+c2*s5)*t+c1*s3)*t+c0;
		w *= (c/ b);
		/*
		COMBINE THE RESULTS
		*/
		u = d * x_alnrel(a/b);
		v = a * (log(b)-1.0);
		if(u <= v) goto S30;
		return w-v-u;
	S30:
		return w-u-v;
	}
	double	x_alnrel(double a)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF THE FUNCTION LN(1 + A)
		-----------------------------------------------------------------------
		*/
	{
		static double p1 = -1.29418923021993;
		static double p2 = .405303492862024;
		static double p3 = -.0178874546012214;
		static double q1 = -1.62752256355323;
		static double q2 = .747811014037616;
		static double q3 = -.0845104217945565;
		double t,t2,w;
		/*
		..
		.. Executable Statements ..
		*/
		if(fabs(a) > 0.375) goto S10;
		t = a/(a+2.0);
		t2 = t*t;
		w = (((p3*t2+p2)*t2+p1)*t2+1.0)/(((q3*t2+q2)*t2+q1)*t2+1.0);
		return 2.0*t*w;
	S10:
		return log(1.0+a);
	}
	double	x_apser(double a,double b,double x,double eps)
		/*
		-----------------------------------------------------------------------
		APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
		A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
		A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
		-----------------------------------------------------------------------
		*/
	{
		static double g = .577215664901533e0;
		double apser,aj,bx,c,j,s,t,tol;
		/*
		..
		.. Executable Statements ..
		*/
		bx = b*x;
		t = x-bx;
		if(b*eps > 2.e-2) goto S10;
		c = log(x)+x_psi(b)+g+t;
		goto S20;
	S10:
		c = log(bx)+g+t;
	S20:
		tol = 5.0e0*eps*fabs(c);
		j = 1.0e0;
		s = 0.0e0;
	S30:
		j += 1.0e0;
		t *= (x-bx/j);
		aj = t/j;
		s += aj;
		if(fabs(aj) > tol) goto S30;
		apser = -(a*(c+s));
		return apser;
	}
	double	x_basym(double *a,double *b,double *lambda,double *eps)
		/*
		-----------------------------------------------------------------------
		ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
		LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
		IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
		A AND B ARE GREATER THAN OR EQUAL TO 15.
		-----------------------------------------------------------------------
		*/
	{
		static double e0 = 1.12837916709551;
		static double e1 = .353553390593274;
		int num = 20;
		/*
		------------------------
		****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
		ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
		THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
		------------------------
		E0 = 2/SQRT(PI)
		E1 = 2**(-3/2)
		------------------------
		*/
		int K3 = 1;
		double bsum,dsum,f,h,h2,hn,j0,j1,r,r0,r1,s,sum,t,t0,t1,u,w,w0,z,z0,
			z2,zn,znm1;
		int i,im1,imj,j,m,mm1,mmj,n,np1;
		double a0[21],b0[21],c[21],d[21],T1,T2;
		/*
		..
		.. Executable Statements ..
		*/
		if(*a >= *b) goto S10;
		h = *a/ *b;
		r0 = 1.0/(1.0+h);
		r1 = (*b-*a)/ *b;
		w0 = 1.0/sqrt(*a*(1.0+h));
		goto S20;
	S10:
		h = *b/ *a;
		r0 = 1.0/(1.0+h);
		r1 = (*b-*a)/ *a;
		w0 = 1.0/sqrt(*b*(1.0+h));
	S20:
		T1 = -(*lambda/ *a);
		T2 = *lambda/ *b;
		f = *a*x_rlog1(&T1)+*b*x_rlog1(&T2);
		t = exp(-f);
		if(t == 0.0) return 0.0;
		z0 = sqrt(f);
		z = 0.5*(z0/e1);
		z2 = f+f;
		a0[0] = 2.0/3.0*r1;
		c[0] = -(0.5*a0[0]);
		d[0] = -c[0];
		j0 = 0.5/e0*x_erfc1(&K3,&z0);
		j1 = e1;
		sum = j0+d[0]*w0*j1;
		s = 1.0;
		h2 = h*h;
		hn = 1.0;
		w = w0;
		znm1 = z;
		zn = z2;
		for(n=2; n<=num; n+=2) {
			hn = h2*hn;
			a0[n-1] = 2.0*r0*(1.0+h*hn)/((double)n+2.0);
			np1 = n+1;
			s += hn;
			a0[np1-1] = 2.0*r1*s/((double)n+3.0);
			for(i=n; i<=np1; i++) {
				r = -(0.5*((double)i+1.0));
				b0[0] = r*a0[0];
				for(m=2; m<=i; m++) {
					bsum = 0.0;
					mm1 = m-1;
					for(j=1; j<=mm1; j++) {
						mmj = m-j;
						bsum += (((double)j*r-(double)mmj)*a0[j-1]*b0[mmj-1]);
					}
					b0[m-1] = r*a0[m-1]+bsum/(double)m;
				}
				c[i-1] = b0[i-1]/((double)i+1.0);
				dsum = 0.0;
				im1 = i-1;
				for(j=1; j<=im1; j++) {
					imj = i-j;
					dsum += (d[imj-1]*c[j-1]);
				}
				d[i-1] = -(dsum+c[i-1]);
			}
			j0 = e1*znm1+((double)n-1.0)*j0;
			j1 = e1*zn+(double)n*j1;
			znm1 = z2*znm1;
			zn = z2*zn;
			w = w0*w;
			t0 = d[n-1]*w*j0;
			w = w0*w;
			t1 = d[np1-1]*w*j1;
			sum += (t0+t1);
			if(fabs(t0)+fabs(t1) <= *eps*sum) goto S80;
		}
	S80:
		u = exp(-x_bcorr(*a,*b));
		return e0*t*u*sum;
	}
	double	x_bcorr(double a0,double b0)
		/*
		-----------------------------------------------------------------------

		EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
		LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
		IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8.

		-----------------------------------------------------------------------
		*/
	{
		static double c0 = .833333333333333e-01;
		static double c1 = -.277777777760991e-02;
		static double c2 = .793650666825390e-03;
		static double c3 = -.595202931351870e-03;
		static double c4 = .837308034031215e-03;
		static double c5 = -.165322962780713e-02;
		double a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2;
		/*
		..
		.. Executable Statements ..
		*/
		a = x_fifdmin1(a0,b0);
		b = x_fifdmax1(a0,b0);
		h = a/b;
		c = h/(1.0+h);
		x = 1.0/(1.0+h);
		x2 = x*x;
		/*
		SET SN = (1 - X**N)/(1 - X)
		*/
		s3 = 1.0+(x+x2);
		s5 = 1.0+(x+x2*s3);
		s7 = 1.0+(x+x2*s5);
		s9 = 1.0+(x+x2*s7);
		s11 = 1.0+(x+x2*s9);
		/*
		SET W = DEL(B) - DEL(A + B)
		*/
		t = SQR(1.0/b);
		w = ((((c5*s11*t+c4*s9)*t+c3*s7)*t+c2*s5)*t+c1*s3)*t+c0;
		w *= (c/b);
		/*
		COMPUTE  DEL(A) + W
		*/
		t = SQR(1.0/a);
		return (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/a+w;
	}
	double	x_betaln(double a0,double b0)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
		-----------------------------------------------------------------------
		E = 0.5*LN(2*PI)
		--------------------------
		*/
	{
		static double e = .918938533204673;
		static double a,b,c,h,u,v,w,z;
		int i,n;
		/*
		..
		.. Executable Statements ..
		*/
		a = x_fifdmin1(a0,b0);
		b = x_fifdmax1(a0,b0);
		if(a >= 8.0) goto S100;
		if(a >= 1.0) goto S20;
		/*
		-----------------------------------------------------------------------
		PROCEDURE WHEN A .LT. 1
		-----------------------------------------------------------------------
		*/
		if(b >= 8.0) goto S10;
		return x_gamln(a)+(x_gamln(b)-x_gamln(a+b));
	S10:
		return x_gamln(a)+x_algdiv(a,b);
	S20:
		/*
		-----------------------------------------------------------------------
		PROCEDURE WHEN 1 .LE. A .LT. 8
		-----------------------------------------------------------------------
		*/
		if(a > 2.0) goto S40;
		if(b > 2.0) goto S30;
		return x_gamln(a)+x_gamln(b)-x_gsumln(&a,&b);
	S30:
		w = 0.0;
		if(b < 8.0) goto S60;
		return x_gamln(a)+x_algdiv(a,b);
	S40:
		/*
		REDUCTION OF A WHEN B .LE. 1000
		*/
		if(b > 1000.0) goto S80;
		n = int(a-1.0);
		w = 1.0;
		for(i=1; i<=n; i++) {
			a -= 1.0;
			h = a/b;
			w *= (h/(1.0+h));
		}
		w = log(w);
		if(b < 8.0) goto S60;
		return w+x_gamln(a)+x_algdiv(a,b);
	S60:
		/*
		REDUCTION OF B WHEN B .LT. 8
		*/
		n = int(b-1.0);
		z = 1.0;
		for(i=1; i<=n; i++) {
			b -= 1.0;
			z *= (b/(a+b));
		}
		return w+log(z)+(x_gamln(a)+(x_gamln(b)-x_gsumln(&a,&b)));
	S80:
		/*
		REDUCTION OF A WHEN B .GT. 1000
		*/
		n = int(a-1.0);
		w = 1.0;
		for(i=1; i<=n; i++) {
			a -= 1.0;
			w *= (a/(1.0+a/b));
		}
		return log(w)-(double)n*log(b)+(x_gamln(a)+x_algdiv(a,b));
	S100:
		/*
		-----------------------------------------------------------------------
		PROCEDURE WHEN A .GE. 8
		-----------------------------------------------------------------------
		*/
		w = x_bcorr(a,b);
		h = a/b;
		c = h/(1.0+h);
		u = -((a-0.5)*log(c));
		v = b*x_alnrel(h);
		if(u <= v) goto S110;
		return -(0.5*log(b))+e+w-v-u;
	S110:
		return -(0.5*log(b))+e+w-u-v;
	}
	double	x_bfrac(double *a,double *b,double *x,double *y,double *lambda,
		double *eps)
		/*
		-----------------------------------------------------------------------
		CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
		IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
		-----------------------------------------------------------------------
		*/
	{
		double bfrac,alpha,an,anp1,beta,bn,bnp1,c,c0,c1,e,n,p,r,r0,s,t,w,yp1;
		/*
		..
		.. Executable Statements ..
		*/
		bfrac = x_brcomp(a,b,x,y);
		if(bfrac == 0.0) return bfrac;
		c = 1.0+*lambda;
		c0 = *b/ *a;
		c1 = 1.0+1.0/ *a;
		yp1 = *y+1.0;
		n = 0.0;
		p = 1.0;
		s = *a+1.0;
		an = 0.0;
		bn = anp1 = 1.0;
		bnp1 = c/c1;
		r = c1/c;
	S10:
		/*
		CONTINUED FRACTION CALCULATION
		*/
		n += 1.0;
		t = n/ *a;
		w = n*(*b-n)**x;
		e = *a/s;
		alpha = p*(p+c0)*e*e*(w**x);
		e = (1.0+t)/(c1+t+t);
		beta = n+w/s+e*(c+n*yp1);
		p = 1.0+t;
		s += 2.0;
		/*
		UPDATE AN, BN, ANP1, AND BNP1
		*/
		t = alpha*an+beta*anp1;
		an = anp1;
		anp1 = t;
		t = alpha*bn+beta*bnp1;
		bn = bnp1;
		bnp1 = t;
		r0 = r;
		r = anp1/bnp1;
		if(fabs(r-r0) <= *eps*r) goto S20;
		/*
		RESCALE AN, BN, ANP1, AND BNP1
		*/
		an /= bnp1;
		bn /= bnp1;
		anp1 = r;
		bnp1 = 1.0;
		goto S10;
	S20:
		/*
		TERMINATION
		*/
		bfrac *= r;
		return bfrac;
	}
	void	x_bgrat(double *a,double *b,double *x,double *y,double *w,
		double *eps,int *ierr)
		/*
		-----------------------------------------------------------------------
		ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
		THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
		THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED.
		IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
		-----------------------------------------------------------------------
		*/
	{
		double bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum,t,t2,u,v,z;
		int i,n,nm1;
		double c[30],d[30],T1;
		/*
		..
		.. Executable Statements ..
		*/
		bm1 = *b-0.5e0-0.5e0;
		nu = *a+0.5e0*bm1;
		if(*y > 0.375e0) goto S10;
		T1 = -*y;
		lnx = x_alnrel(T1);
		goto S20;
	S10:
		lnx = log(*x);
	S20:
		z = -(nu*lnx);
		if(*b*z == 0.0e0) goto S70;
		/*
		COMPUTATION OF THE EXPANSION
		SET R = EXP(-Z)*Z**B/GAMMA(B)
		*/
		r = *b*(1.0e0+x_gam1(*b))*exp(*b*log(z));
		r *= (exp(*a*lnx)*exp(0.5e0*bm1*lnx));
		u = x_algdiv(*b,*a)+*b*log(nu);
		u = r*exp(-u);
		if(u == 0.0e0) goto S70;
		x_grat1(b,&z,&r,&p,&q,eps);
		v = 0.25e0*SQR(1.0e0/nu);
		t2 = 0.25e0*lnx*lnx;
		l = *w/u;
		j = q/r;
		sum = j;
		t = cn = 1.0e0;
		n2 = 0.0e0;
		for(n=1; n<=30; n++) {
			bp2n = *b+n2;
			j = (bp2n*(bp2n+1.0e0)*j+(z+bp2n+1.0e0)*t)*v;
			n2 += 2.0e0;
			t *= t2;
			cn /= (n2*(n2+1.0e0));
			c[n-1] = cn;
			s = 0.0e0;
			if(n == 1) goto S40;
			nm1 = n-1;
			coef = *b-(double)n;
			for(i=1; i<=nm1; i++) {
				s += (coef*c[i-1]*d[n-i-1]);
				coef += *b;
			}
		S40:
			d[n-1] = bm1*cn+s/(double)n;
			dj = d[n-1]*j;
			sum += dj;
			if(sum <= 0.0e0) goto S70;
			if(fabs(dj) <= *eps*(sum+l)) goto S60;
		}
	S60:
		/*
		ADD THE RESULTS TO W
		*/
		*ierr = 0;
		*w += (u*sum);
		return;
	S70:
		/*
		THE EXPANSION CANNOT BE COMPUTED
		*/
		*ierr = 1;
		return;
	}
	double	x_bpser(double a,double b,double x,double eps)
		/*
		-----------------------------------------------------------------------
		POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1
		OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
		-----------------------------------------------------------------------
		*/
	{
		double bpser,a0,apb,b0,c,n,sum,t,tol,u,w,z;
		int i,m;
		/*
		..
		.. Executable Statements ..
		*/
		bpser = 0.0;
		if(x == 0.0) return bpser;
		/*
		-----------------------------------------------------------------------
		COMPUTE THE FACTOR X**A/(A*BETA(A,B))
		-----------------------------------------------------------------------
		*/
		a0 = x_fifdmin1(a,b);
		if(a0 < 1.0) goto S10;
		z = a*log(x)-x_betaln(a,b);
		bpser = exp(z)/ a;
		goto S100;
	S10:
		b0 = x_fifdmax1(a,b);
		if(b0 >= 8.0) goto S90;
		if(b0 > 1.0) goto S40;
		/*
		PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
		*/
		bpser = pow(x,a);
		if(bpser == 0.0) return bpser;
		apb = a+b;
		if(apb > 1.0) goto S20;
		z = 1.0+x_gam1(apb);
		goto S30;
	S20:
		u = a+b-1.;
		z = (1.0+x_gam1(u))/apb;
	S30:
		c = (1.0+x_gam1(a))*(1.0+x_gam1(b))/z;
		bpser *= (c*(b/apb));
		goto S100;
	S40:
		/*
		PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
		*/
		u = x_gamln1(a0);
		m = int(b0-1.0);
		if(m < 1) goto S60;
		c = 1.0;
		for(i=1; i<=m; i++) {
			b0 -= 1.0;
			c *= (b0/(a0+b0));
		}
		u = log(c)+u;
	S60:
		z = a*log(x)-u;
		b0 -= 1.0;
		apb = a0+b0;
		if(apb > 1.0) goto S70;
		t = 1.0+x_gam1(apb);
		goto S80;
	S70:
		u = a0+b0-1.;
		t = (1.0+x_gam1(u))/apb;
	S80:
		bpser = exp(z)*(a0/ a)*(1.0+x_gam1(b0))/t;
		goto S100;
	S90:
		/*
		PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
		*/
		u = x_gamln1(a0)+x_algdiv(a0,b0);
		z = a*log(x)-u;
		bpser = a0/ a*exp(z);
	S100:
		if(bpser == 0.0 || a <= 0.1*eps) return bpser;
		/*
		-----------------------------------------------------------------------
		COMPUTE THE SERIES
		-----------------------------------------------------------------------
		*/
		sum = n = 0.0;
		c = 1.0;
		tol = eps/ a;
	S110:
		n += 1.0;
		c *= ((0.5+(0.5-b/n))*x);
		w = c/(a+n);
		sum += w;
		if(fabs(w) > tol) goto S110;
		bpser *= (1.0+a*sum);
		return bpser;
	}
	double	x_brcmp1(int *mu,double *a,double *b,double *x,double *y)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
		-----------------------------------------------------------------------
		*/
	{
		static double Const = .398942280401433e0;
		double brcmp1,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
		int i,n;
		/*
		-----------------
		CONST = 1/SQRT(2*PI)
		-----------------
		*/
		double T1,T2,T3,T4;
		/*
		..
		.. Executable Statements ..
		*/
		a0 = x_fifdmin1(*a,*b);
		if(a0 >= 8.0e0) goto S130;
		if(*x > 0.375e0) goto S10;
		lnx = log(*x);
		T1 = -*x;
		lny = x_alnrel(T1);
		goto S30;
	S10:
		if(*y > 0.375e0) goto S20;
		T2 = -*y;
		lnx = x_alnrel(T2);
		lny = log(*y);
		goto S30;
	S20:
		lnx = log(*x);
		lny = log(*y);
	S30:
		z = *a*lnx+*b*lny;
		if(a0 < 1.0e0) goto S40;
		z -= x_betaln(*a,*b);
		return x_esum(*mu,z);
	S40:
		/*
		-----------------------------------------------------------------------
		PROCEDURE FOR A .LT. 1 OR B .LT. 1
		-----------------------------------------------------------------------
		*/
		b0 = x_fifdmax1(*a,*b);
		if(b0 >= 8.0e0) goto S120;
		if(b0 > 1.0e0) goto S70;
		/*
		ALGORITHM FOR B0 .LE. 1
		*/
		brcmp1 = x_esum(*mu,z);
		if(brcmp1 == 0.0e0) return brcmp1;
		apb = *a+*b;
		if(apb > 1.0e0) goto S50;
		z = 1.0e0+x_gam1(apb);
		goto S60;
	S50:
		u = *a+*b-1.e0;
		z = (1.0e0+x_gam1(u))/apb;
	S60:
		c = (1.0e0+x_gam1(*a))*(1.0e0+x_gam1(*b))/z;
		return brcmp1*(a0*c)/(1.0e0+a0/b0);
	S70:
		/*
		ALGORITHM FOR 1 .LT. B0 .LT. 8
		*/
		u = x_gamln1(a0);
		n = int(b0-1.0e0);
		if(n < 1) goto S90;
		c = 1.0e0;
		for(i=1; i<=n; i++) {
			b0 -= 1.0e0;
			c *= (b0/(a0+b0));
		}
		u = log(c)+u;
	S90:
		z -= u;
		b0 -= 1.0e0;
		apb = a0+b0;
		if(apb > 1.0e0) goto S100;
		t = 1.0e0+x_gam1(apb);
		goto S110;
	S100:
		u = a0+b0-1.e0;
		t = (1.0e0+x_gam1(u))/apb;
	S110:
		return a0*x_esum(*mu,z)*(1.0e0+x_gam1(b0))/t;
	S120:
		/*
		ALGORITHM FOR B0 .GE. 8
		*/
		u = x_gamln1(a0)+x_algdiv(a0,b0);
		T3 = z-u;
		return a0*x_esum(*mu,T3);
	S130:
		/*
		-----------------------------------------------------------------------
		PROCEDURE FOR A .GE. 8 AND B .GE. 8
		-----------------------------------------------------------------------
		*/
		if(*a > *b) goto S140;
		h = *a/ *b;
		x0 = h/(1.0e0+h);
		y0 = 1.0e0/(1.0e0+h);
		lambda = *a-(*a+*b)**x;
		goto S150;
	S140:
		h = *b/ *a;
		x0 = 1.0e0/(1.0e0+h);
		y0 = h/(1.0e0+h);
		lambda = (*a+*b)**y-*b;
	S150:
		e = -(lambda/ *a);
		if(fabs(e) > 0.6e0) goto S160;
		u = x_rlog1(&e);
		goto S170;
	S160:
		u = e-log(*x/x0);
	S170:
		e = lambda/ *b;
		if(fabs(e) > 0.6e0) goto S180;
		v = x_rlog1(&e);
		goto S190;
	S180:
		v = e-log(*y/y0);
	S190:
		T4 = -(*a*u+*b*v);
		z = x_esum(*mu,T4);
		return Const*sqrt(*b*x0)*z*exp(-x_bcorr(*a,*b));
	}
	double	x_brcomp(double *a,double *b,double *x,double *y)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF X**A*Y**B/BETA(A,B)
		-----------------------------------------------------------------------
		*/
	{
		static double Const = .398942280401433;
		double brcomp,a0,apb,b0,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z;
		int i,n;
		/*
		-----------------
		CONST = 1/SQRT(2*PI)
		-----------------
		*/
		/*
		..
		.. Executable Statements ..
		*/
		brcomp = 0.0;
		if(*x == 0.0 || *y == 0.0) return brcomp;
		a0 = x_fifdmin1(*a,*b);
		if(a0 >= 8.0) goto S130;
		if(*x > 0.375) goto S10;
		lnx = log(*x);
		lny = x_alnrel(-*x);
		goto S30;
	S10:
		if(*y > 0.375) goto S20;
		lnx = x_alnrel(-*y);
		lny = log(*y);
		goto S30;
	S20:
		lnx = log(*x);
		lny = log(*y);
	S30:
		z = *a*lnx+*b*lny;
		if(a0 < 1.0) goto S40;
		z -= x_betaln(*a,*b);
		return exp(z);
	S40:
		/*
		-----------------------------------------------------------------------
		PROCEDURE FOR A .LT. 1 OR B .LT. 1
		-----------------------------------------------------------------------
		*/
		b0 = x_fifdmax1(*a,*b);
		if(b0 >= 8.0) goto S120;
		if(b0 > 1.0) goto S70;
		/*
		ALGORITHM FOR B0 .LE. 1
		*/
		brcomp = exp(z);
		if(brcomp == 0.0) return brcomp;
		apb = *a+*b;
		if(apb > 1.0) goto S50;
		z = 1.0+x_gam1(apb);
		goto S60;
	S50:
		u = *a+*b-1.;
		z = (1.0+x_gam1(u))/apb;
	S60:
		c = (1.0+x_gam1(*a))*(1.0+x_gam1(*b))/z;
		brcomp = brcomp*(a0*c)/(1.0+a0/b0);
		return brcomp;
	S70:
		/*
		ALGORITHM FOR 1 .LT. B0 .LT. 8
		*/
		u = x_gamln1(a0);
		n = int(b0-1.0);
		if(n < 1) goto S90;
		c = 1.0;
		for(i=1; i<=n; i++) {
			b0 -= 1.0;
			c *= (b0/(a0+b0));
		}
		u = log(c)+u;
	S90:
		z -= u;
		b0 -= 1.0;
		apb = a0+b0;
		if(apb > 1.0) goto S100;
		t = 1.0+x_gam1(apb);
		goto S110;
	S100:
		u = a0+b0-1.;
		t = (1.0+x_gam1(u))/apb;
	S110:
		return a0*exp(z)*(1.0+x_gam1(b0))/t;
	S120:
		/*
		ALGORITHM FOR B0 .GE. 8
		*/
		u = x_gamln1(a0)+x_algdiv(a0,b0);
		return a0*exp(z-u);
	S130:
		/*
		-----------------------------------------------------------------------
		PROCEDURE FOR A .GE. 8 AND B .GE. 8
		-----------------------------------------------------------------------
		*/
		if(*a > *b) goto S140;
		h = *a/ *b;
		x0 = h/(1.0+h);
		y0 = 1.0/(1.0+h);
		lambda = *a-(*a+*b)**x;
		goto S150;
	S140:
		h = *b/ *a;
		x0 = 1.0/(1.0+h);
		y0 = h/(1.0+h);
		lambda = (*a+*b)**y-*b;
	S150:
		e = -(lambda/ *a);
		if(fabs(e) > 0.6) goto S160;
		u = x_rlog1(&e);
		goto S170;
	S160:
		u = e-log(*x/x0);
	S170:
		e = lambda/ *b;
		if(fabs(e) > 0.6) goto S180;
		v = x_rlog1(&e);
		goto S190;
	S180:
		v = e-log(*y/y0);
	S190:
		z = exp(-(*a*u+*b*v));
		brcomp = Const*sqrt(*b*x0)*z*exp(-x_bcorr(*a,*b));
		return brcomp;
	}
	double	x_bup(double *a,double *b,double *x,double *y,int *n,double *eps)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
		EPS IS THE TOLERANCE USED.
		-----------------------------------------------------------------------
		*/
	{
		int K1 = 1;
		int K2 = 0;
		double bup,ap1,apb,d,l,r,t,w;
		int i,k,kp1,mu,nm1;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		OBTAIN THE SCALING FACTOR EXP(-MU) AND
		EXP(MU)*(X**A*Y**B/BETA(A,B))/A
		*/
		apb = *a+*b;
		ap1 = *a+1.0e0;
		mu = 0;
		d = 1.0e0;
		if(*n == 1 || *a < 1.0e0) goto S10;
		if(apb < 1.1e0*ap1) goto S10;
		mu = int(fabs(x_exparg(&K1)));
		k = int(x_exparg(&K2));
		if(k < mu) mu = k;
		t = mu;
		d = exp(-t);
	S10:
		bup = x_brcmp1(&mu,a,b,x,y)/ *a;
		if(*n == 1 || bup == 0.0e0) return bup;
		nm1 = *n-1;
		w = d;
		/*
		LET K BE THE INDEX OF THE MAXIMUM TERM
		*/
		k = 0;
		if(*b <= 1.0e0) goto S50;
		if(*y > 1.e-4) goto S20;
		k = nm1;
		goto S30;
	S20:
		r = (*b-1.0e0)**x/ *y-*a;
		if(r < 1.0e0) goto S50;
		t = nm1; k = nm1;
		if(r < t) k = (int)r;
	S30:
		/*
		ADD THE INCREASING TERMS OF THE SERIES
		*/
		for(i=1; i<=k; i++) {
			l = i-1;
			d = (apb+l)/(ap1+l)**x*d;
			w += d;
		}
		if(k == nm1) goto S70;
	S50:
		/*
		ADD THE REMAINING TERMS OF THE SERIES
		*/
		kp1 = k+1;
		for(i=kp1; i<=nm1; i++) {
			l = i-1;
			d = (apb+l)/(ap1+l)**x*d;
			w += d;
			if(d <= *eps*w) goto S70;
		}
	S70:
		/*
		TERMINATE THE PROCEDURE
		*/
		return bup*w;
	}

	double	x_dlanor(double *x)
		/*
		**********************************************************************

		double dlanor(double *x)
		Double precision Logarithm of the Asymptotic Normal


		Function


		Computes the logarithm of the cumulative normal distribution
		from abs( x ) to infinity for abs( x ) >= 5.


		Arguments


		X --> Value at which cumulative normal to be evaluated
		DOUBLE PRECISION X


		Method


		23 term expansion of formula 26.2.12 of Abramowitz and Stegun.
		The relative error at X = 5 is about 0.5E-5.


		Note


		ABS(X) must be >= 5 else there is an error stop.

		**********************************************************************
		*/
	{
#define dlsqpi 0.91893853320467274177e0
		static double coef[12] = {
			-1.0e0,3.0e0,-15.0e0,105.0e0,-945.0e0,10395.0e0,-135135.0e0,2027025.0e0,
			-34459425.0e0,654729075.0e0,-13749310575.e0,316234143225.0e0
		};
		static int K1 = 12;
		double approx,correc,xx,xx2;
		/*
		..
		.. Executable Statements ..
		*/
		xx = fabs(*x);
		if(xx < 5.0e0) halt(" Argument too small in DLANOR");
		approx = -dlsqpi-0.5e0*xx*xx-log(xx);
		xx2 = xx*xx;
		correc = x_dblEvalPolyX(coef, 12, 1.0/xx2)/xx2;
		return approx+x_dln1px(correc);
#undef dlsqpi
	}
	double	x_dln1mx(double x)
		/*
		**********************************************************************

		double dln1mx(double *x)
		Double precision LN(1-X)


		Function


		Returns ln(1-x) for small x (good accuracy if x .le. 0.1).
		Note that the obvious code of
		LOG(1.0-X)
		won't work for small X because 1.0-X loses accuracy


		Arguments


		X --> Value for which ln(1-x) is desired.
		X is DOUBLE PRECISION


		Method


		If X > 0.1, the obvious code above is used ELSE
		The Taylor series for 1-x is expanded to 20 terms.

		**********************************************************************
		*/
	{
		return x_dln1px(-x);
	}
	double	x_dln1px(double a)
		/*
		**********************************************************************

		double dln1px(double *a)
		Double precision LN(1+X)


		Function


		Returns ln(1+x)
		Note that the obvious code of
		LOG(1.0+X)
		won't work for small X because 1.0+X loses accuracy


		Arguments


		X --> Value for which ln(1-x) is desired.
		X is DOUBLE PRECISION


		Method


		Renames ALNREL from:
		DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
		Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
		Trans. Math.  Softw. 18 (1993), 360-373.

		**********************************************************************
		-----------------------------------------------------------------------
		EVALUATION OF THE FUNCTION LN(1 + A)
		-----------------------------------------------------------------------
		*/
	{
		static double p1 = -.129418923021993e+01;
		static double p2 = .405303492862024e+00;
		static double p3 = -.178874546012214e-01;
		static double q1 = -.162752256355323e+01;
		static double q2 = .747811014037616e+00;
		static double q3 = -.845104217945565e-01;
		double t,t2,w;
		/*
		..
		.. Executable Statements ..
		*/
		if(fabs(a) > 0.375) return log(1. + a);
		t = a/(a+2.0e0);
		t2 = t*t;
		w = (((p3*t2+p2)*t2+p1)*t2+1.0e0)/(((q3*t2+q2)*t2+q1)*t2+1.0e0);
		return 2.0e0*t*w;
	}
	double	x_dlnbet(double a0,double b0)
		/*
		**********************************************************************

		double dlnbet(a0,b0)
		Double precision LN of the complete BETa


		Function


		Returns the natural log of the complete beta function,
		i.e.,

		ln( Gamma(a)*Gamma(b) / Gamma(a+b)


		Arguments


		A,B --> The (symmetric) arguments to the complete beta
		DOUBLE PRECISION A, B


		Method


		Renames BETALN from:
		DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
		Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
		Trans. Math.  Softw. 18 (1993), 360-373.

		**********************************************************************
		-----------------------------------------------------------------------
		EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
		-----------------------------------------------------------------------
		E = 0.5*LN(2*PI)
		--------------------------
		*/
	{
		static double e = .918938533204673;
		double a,b,c,h,u,v,w,z;
		int i,n;
		/*
		..
		.. Executable Statements ..
		*/
		a = x_fifdmin1(a0,b0);
		b = x_fifdmax1(a0,b0);
		if(a >= 8.0) goto S100;
		if(a >= 1.0) goto S20;
		/*
		-----------------------------------------------------------------------
		PROCEDURE WHEN A .LT. 1
		-----------------------------------------------------------------------
		*/
		if(b >= 8.0) goto S10;
		return x_gamln(a)+(x_gamln(b)-x_gamln(a+b));
	S10:
		return x_gamln(a)+x_algdiv(a,b);
	S20:
		/*
		-----------------------------------------------------------------------
		PROCEDURE WHEN 1 .LE. A .LT. 8
		-----------------------------------------------------------------------
		*/
		if(a > 2.0) goto S40;
		if(b > 2.0) goto S30;
		return x_gamln(a)+x_gamln(b)-x_gsumln(&a,&b);
	S30:
		w = 0.0;
		if(b < 8.0) goto S60;
		return x_gamln(a)+x_algdiv(a,b);
	S40:
		/*
		REDUCTION OF A WHEN B .LE. 1000
		*/
		if(b > 1000.0) goto S80;
		n = int(a-1.0);
		w = 1.0;
		for(i=1; i<=n; i++) {
			a -= 1.0;
			h = a/b;
			w *= (h/(1.0+h));
		}
		w = log(w);
		if(b < 8.0) goto S60;
		return w+x_gamln(a)+x_algdiv(a,b);
	S60:
		/*
		REDUCTION OF B WHEN B .LT. 8
		*/
		n = int(b-1.0);
		z = 1.0;
		for(i=1; i<=n; i++) {
			b -= 1.0;
			z *= (b/(a+b));
		}
		return w+log(z)+(x_gamln(a)+(x_gamln(b)-x_gsumln(&a,&b)));
	S80:
		/*
		REDUCTION OF A WHEN B .GT. 1000
		*/
		n = int(a-1.0);
		w = 1.0;
		for(i=1; i<=n; i++) {
			a -= 1.0;
			w *= (a/(1.0+a/b));
		}
		return log(w)-(double)n*log(b)+(x_gamln(a)+x_algdiv(a,b));
	S100:
		/*
		-----------------------------------------------------------------------
		PROCEDURE WHEN A .GE. 8
		-----------------------------------------------------------------------
		*/
		w = x_bcorr(a,b);
		h = a/b;
		c = h/(1.0+h);
		u = -((a-0.5)*log(c));
		v = b*x_alnrel(h);
		if(u <= v) goto S110;
		return -(0.5*log(b))+e+w-v-u;
	S110:
		return -(0.5*log(b))+e+w-u-v;
	}
	double	x_dlngam(double a)
		/*
		**********************************************************************

		double dlngam(double *a)
		Double precision LN of the GAMma function


		Function


		Returns the natural logarithm of GAMMA(X).


		Arguments


		X --> value at which scaled log gamma is to be returned
		X is DOUBLE PRECISION


		Method


		Renames GAMLN from:
		DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
		Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
		Trans. Math.  Softw. 18 (1993), 360-373.

		**********************************************************************
		-----------------------------------------------------------------------
		EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
		-----------------------------------------------------------------------
		WRITTEN BY ALFRED H. MORRIS
		NAVAL SURFACE WARFARE CENTER
		DAHLGREN, VIRGINIA
		--------------------------
		D = 0.5*(LN(2*PI) - 1)
		--------------------------
		*/
	{
		static double c0 = .833333333333333e-01;
		static double c1 = -.277777777760991e-02;
		static double c2 = .793650666825390e-03;
		static double c3 = -.595202931351870e-03;
		static double c4 = .837308034031215e-03;
		static double c5 = -.165322962780713e-02;
		static double d = .418938533204673e0;
		double t,w;
		int i,n;
		double T1;
		/*
		..
		.. Executable Statements ..
		*/
		if(a > 0.8e0) goto S10;
		return x_gamln1(a)-log(a);
	S10:
		if(a > 2.25e0) goto S20;
		t = a-1.0;
		return x_gamln1(t);
	S20:
		if(a >= 10.0e0) goto S40;
		n = int(a-1.25e0);
		t = a;
		w = 1.0e0;
		for(i=1; i<=n; i++) {
			t -= 1.0e0;
			w = t*w;
		}
		T1 = t-1.0;
		return x_gamln1(T1)+log(w);
	S40:
		t = SQR(1.0e0/ a);
		w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/ a;
		return d+w+(a-0.5e0)*(log(a)-1.0e0);
	}
	double	x_dstrem(double z)
	{
		/*
		**********************************************************************
		double dstrem(double *z)
		Double precision Sterling Remainder
		Function
		Returns   Log(Gamma(Z))  -  Sterling(Z)  where   Sterling(Z)  is
		Sterling's Approximation to Log(Gamma(Z))
		Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5 ) * LOG( Z ) - Z
		Arguments
		Z --> Value at which Sterling remainder calculated
		Must be positive.
		DOUBLE PRECISION Z
		Method
		If Z >= 6 uses 9 terms of series in Bernoulli numbers
		(Values calculated using Maple)
		Otherwise computes difference explicitly
		**********************************************************************
		*/
#define hln2pi 0.91893853320467274178e0
#define ncoef 10
		static double coef[ncoef] = {
			0.0e0,0.0833333333333333333333333333333e0,
			-0.00277777777777777777777777777778e0,0.000793650793650793650793650793651e0,
			-0.000595238095238095238095238095238e0,
			0.000841750841750841750841750841751e0,-0.00191752691752691752691752691753e0,
			0.00641025641025641025641025641026e0,-0.0295506535947712418300653594771e0,
			0.179644372368830573164938490016e0
		};
		static int K1 = 10;
		double dstrem,sterl;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		For information, here are the next 11 coefficients of the
		remainder term in Sterling's formula
		-1.39243221690590111642743221691
		13.4028640441683919944789510007
		-156.848284626002017306365132452
		2193.10333333333333333333333333
		-36108.7712537249893571732652192
		691472.268851313067108395250776
		-0.152382215394074161922833649589D8
		0.382900751391414141414141414141D9
		-0.108822660357843910890151491655D11
		0.347320283765002252252252252252D12
		-0.123696021422692744542517103493D14
		*/
		if(z <= 0.0e0) halt("Zero or negative argument in DSTREM");
		if(!(z > 6.0e0)) goto S10;
		dstrem = x_dblEvalPolyX(coef,10,1.0/SQR(z))*z;
		goto S20;
	S10:
		sterl = hln2pi+(z-0.5e0)*log(z)-z;
		dstrem = x_dlngam(z)-sterl;
	S20:
		return dstrem;
#undef hln2pi
#undef ncoef
	}
	double	x_dt1(double p,double q,double df)
		/*
		**********************************************************************

		double dt1(double *p,double *q,double *df)
		Double precision Initialize Approximation to
		INVerse of the cumulative T distribution


		Function


		Returns  the  inverse   of  the T   distribution   function, i.e.,
		the integral from 0 to INVT of the T density is P. This is an
		initial approximation


		Arguments


		P --> The p-value whose inverse from the T distribution is
		desired.
		P is DOUBLE PRECISION

		Q --> 1-P.
		Q is DOUBLE PRECISION

		DF --> Degrees of freedom of the T distribution.
		DF is DOUBLE PRECISION

		**********************************************************************
		*/
	{
		static double coef[4][5] = {
			{    1.0e0,     1.0e0,    0.0e0,   0.0e0,  0.0e0 },
			{    3.0e0,    16.0e0,    5.0e0,   0.0e0,  0.0e0 },
			{  -15.0e0,    17.0e0,   19.0e0,   3.0e0,  0.0e0 },
			{ -945.0e0, -1920.0e0, 1482.0e0, 776.0e0, 79.0e0 }
		};
		static double denom[4] = {
			4.0e0,96.0e0,384.0e0,92160.0e0
		};
		static int ideg[4] = {
			2,3,4,5
		};
		double denpow,sum,term,x,xx;
		/*
		..
		.. Executable Statements ..
		*/
		x = fabs(x_dblNormalDistInv(p,q));
		xx = x*x;
		sum = x;
		denpow = 1.0;
		for (wsUint i=0; i<4; i++) {
			term = x_dblEvalPolyX(&coef[i][0],ideg[i],xx)*x;
			denpow *= df;
			sum += (term/(denpow*denom[i]));
		}
		return (p >= 0.5) ? sum : -sum;
	}
	double	x_erf1(double *x)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF THE REAL ERROR FUNCTION
		-----------------------------------------------------------------------
		*/
	{
		static double c = .564189583547756;
		static double a[5] = {
			.771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
			.479137145607681e-01,.128379167095513e+00
		};
		static double b[3] = {
			.301048631703895e-02,.538971687740286e-01,.375795757275549e+00
		};
		static double p[8] = {
			-1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
			4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
			4.51918953711873e+02,3.00459261020162e+02
		};
		static double q[8] = {
			1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
			2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
			7.90950925327898e+02,3.00459260956983e+02
		};
		static double r[5] = {
			2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
			4.65807828718470e+00,2.82094791773523e-01
		};
		static double s[4] = {
			9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
			1.80124575948747e+01
		};
		double erf1,ax,bot,t,top,x2;
		/*
		..
		.. Executable Statements ..
		*/
		ax = fabs(*x);
		if(ax > 0.5) goto S10;
		t = *x**x;
		top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0;
		bot = ((b[0]*t+b[1])*t+b[2])*t+1.0;
		return *x*(top/bot);
	S10:
		if(ax > 4.0) goto S20;
		top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[
			7];
		bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[
			7];
		erf1 = 0.5+(0.5-exp(-(*x**x))*top/bot);
		if(*x < 0.0) erf1 = -erf1;
		return erf1;
	S20:
		if(ax >= 5.8) goto S30;
		x2 = *x**x;
		t = 1.0/x2;
		top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
		bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0;
		erf1 = (c-top/(x2*bot))/ax;
		erf1 = 0.5+(0.5-exp(-x2)*erf1);
		if(*x < 0.0) erf1 = -erf1;
		return erf1;
	S30:
		return x_fifdsign(1.0,*x);
	}
	double	x_erfc1(int *ind,double *x)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

		ERFC1(IND,X) = ERFC(X)            IF IND = 0
		ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
		-----------------------------------------------------------------------
		*/
	{
		static double c = .564189583547756e0;
		static double a[5] = {
			.771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
			.479137145607681e-01,.128379167095513e+00
		};
		static double b[3] = {
			.301048631703895e-02,.538971687740286e-01,.375795757275549e+00
		};
		static double p[8] = {
			-1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
			4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
			4.51918953711873e+02,3.00459261020162e+02
		};
		static double q[8] = {
			1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
			2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
			7.90950925327898e+02,3.00459260956983e+02
		};
		static double r[5] = {
			2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
			4.65807828718470e+00,2.82094791773523e-01
		};
		static double s[4] = {
			9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
			1.80124575948747e+01
		};
		static int K1 = 1;
		double erfc1,ax,bot,e,t,top,w;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		ABS(X) .LE. 0.5
		*/
		ax = fabs(*x);
		if(ax > 0.5e0) goto S10;
		t = *x**x;
		top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
		bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
		erfc1 = 0.5e0+(0.5e0-*x*(top/bot));
		if(*ind != 0) erfc1 = exp(t)*erfc1;
		return erfc1;
	S10:
		/*
		0.5 .LT. ABS(X) .LE. 4
		*/
		if(ax > 4.0e0) goto S20;
		top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[
			7];
		bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[
			7];
		erfc1 = top/bot;
		goto S40;
	S20:
		/*
		ABS(X) .GT. 4
		*/
		if(*x <= -5.6e0) goto S60;
		if(*ind != 0) goto S30;
		if(*x > 100.0e0) goto S70;
		if(*x**x > -x_exparg(&K1)) goto S70;
	S30:
		t = SQR(1.0e0/ *x);
		top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
		bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
		erfc1 = (c-t*top/bot)/ax;
	S40:
		/*
		FINAL ASSEMBLY
		*/
		if(*ind == 0) goto S50;
		if(*x < 0.0e0) erfc1 = 2.0e0*exp(*x**x)-erfc1;
		return erfc1;
	S50:
		w = *x**x;
		t = w;
		e = w-t;
		erfc1 = (0.5e0+(0.5e0-e))*exp(-t)*erfc1;
		if(*x < 0.0e0) erfc1 = 2.0e0-erfc1;
		return erfc1;
	S60:
		/*
		LIMIT VALUE FOR LARGE NEGATIVE X
		*/
		erfc1 = 2.0e0;
		if(*ind != 0) erfc1 = 2.0e0*exp(*x**x);
		return erfc1;
	S70:
		/*
		LIMIT VALUE FOR LARGE POSITIVE X
		WHEN IND = 0
		*/
		return 0.0;
	}
	double	x_esum(int mu,double x)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF EXP(MU + X)
		-----------------------------------------------------------------------
		*/
	{
		double w;
		/*
		..
		.. Executable Statements ..
		*/
		if(x > 0.0e0) goto S10;
		if(mu < 0) goto S20;
		w = (double)mu+x;
		if(w > 0.0e0) goto S20;
		return exp(w);
	S10:
		if(mu > 0) goto S20;
		w = (double)mu+x;
		if(w < 0.0e0) goto S20;
		return exp(w);
	S20:
		return exp(mu)*exp(x);
	}
	double	x_exparg(int *l)
		/*
		--------------------------------------------------------------------
		IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
		EXP(W) CAN BE COMPUTED.

		IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
		WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.

		NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
		--------------------------------------------------------------------
		*/
	{
		/*
		..
		.. Executable Statements ..
		*/
		if(*l == 0) goto S50;
		return -708.389334568083577311;//0.99999*-708.3964185322689;//((double)m*lnb);
	S50:
		return 709.775615066259866112;//0.99999*709.7827128933888;//((double)m*lnb);
	}
	double	x_fpser(double a,double b,double x,double eps)
		/*
		-----------------------------------------------------------------------

		EVALUATION OF I (A,B)
		X

		FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.

		-----------------------------------------------------------------------

		SET  FPSER = X**A
		*/
	{
		static int K1 = 1;
		double fpser,an,c,s,t,tol;
		/*
		..
		.. Executable Statements ..
		*/
		fpser = 1.0;
		if(a <= 1.e-3*eps) goto S10;
		fpser = 0.0;
		t = a*log(x);
		if(t < x_exparg(&K1)) return fpser;
		fpser = exp(t);
	S10:
		/*
		NOTE THAT 1/B(A,B) = B
		*/
		fpser = b/ a*fpser;
		tol = eps/ a;
		an = a+1.0;
		t = x;
		s = t/an;
	S20:
		an += 1.0;
		t = x*t;
		c = t/an;
		s += c;
		if(fabs(c) > tol) goto S20;
		fpser *= (1.0+a*s);
		return fpser;
	}
	double	x_gam1(double a)
		/*
		------------------------------------------------------------------
		COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5
		------------------------------------------------------------------
		*/
	{
		static double s1 = .273076135303957e+00;
		static double s2 = .559398236957378e-01;
		static double p[7] = {
			.577215664901533e+00,-.409078193005776e+00,-.230975380857675e+00,
			.597275330452234e-01,.766968181649490e-02,-.514889771323592e-02,
			.589597428611429e-03
		};
		static double q[5] = {
			.100000000000000e+01,.427569613095214e+00,.158451672430138e+00,
			.261132021441447e-01,.423244297896961e-02
		};
		static double r[9] = {
			-.422784335098468e+00,-.771330383816272e+00,-.244757765222226e+00,
			.118378989872749e+00,.930357293360349e-03,-.118290993445146e-01,
			.223047661158249e-02,.266505979058923e-03,-.132674909766242e-03
		};
		double bot,d,t,top,w,T1;
		/*
		..
		.. Executable Statements ..
		*/
		t = a;
		d = a-0.5;
		if(d > 0.0) t = d-0.5;
		T1 = t;
		if(T1 < 0) goto S40;
		else if(T1 == 0) goto S10;
		else  goto S20;
	S10:
		return 0.0;
	S20:
		top = (((((p[6]*t+p[5])*t+p[4])*t+p[3])*t+p[2])*t+p[1])*t+p[0];
		bot = (((q[4]*t+q[3])*t+q[2])*t+q[1])*t+1.0e0;
		w = top/bot;
		if(d > 0.0) goto S30;
		return a*w;
	S30:
		return t/ a*(w-1.0);
	S40:
		top = (((((((r[8]*t+r[7])*t+r[6])*t+r[5])*t+r[4])*t+r[3])*t+r[2])*t+r[1])*t+
			r[0];
		bot = (s2*t+s1)*t+1.0;
		w = top/bot;
		if(d > 0.0) goto S50;
		return a*(w+1.0);
	S50:
		return t*w/ a;
	}
	void	x_gaminv(double *a,double *x,double *x0,double *p,double *q,
		int *ierr)
		/*
		----------------------------------------------------------------------
		INVERSE INCOMPLETE GAMMA RATIO FUNCTION

		GIVEN POSITIVE A, AND NONEGATIVE P AND Q WHERE P + Q = 1.
		THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER
		ITERATION IS EMPLOYED. THE ROUTINE ATTEMPTS TO COMPUTE X
		TO 10 SIGNIFICANT DIGITS IF THIS IS POSSIBLE FOR THE
		PARTICULAR COMPUTER ARITHMETIC BEING USED.

		------------

		X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
		AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT
		NUMBER AVAILABLE. OTHERWISE, GAMINV ATTEMPTS TO OBTAIN
		A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE
		IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.

		X0 IS AN OPTIONAL INITIAL APPROXIMATION FOR X. IF THE USER
		DOES NOT WISH TO SUPPLY AN INITIAL APPROXIMATION, THEN SET
		X0 .LE. 0.

		IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
		WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING
		VALUES ...

		IERR =  0    THE SOLUTION WAS OBTAINED. ITERATION WAS
		NOT USED.
		IERR.GT.0    THE SOLUTION WAS OBTAINED. IERR ITERATIONS
		WERE PERFORMED.
		IERR = -2    (INPUT ERROR) A .LE. 0
		IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A
		IS TOO LARGE.
		IERR = -4    (INPUT ERROR) P + Q .NE. 1
		IERR = -6    20 ITERATIONS WERE PERFORMED. THE MOST
		RECENT VALUE OBTAINED FOR X IS GIVEN.
		THIS CANNOT OCCUR IF X0 .LE. 0.
		IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X.
		THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
		IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE
		ROUTINE IS NOT CERTAIN OF ITS ACCURACY.
		ITERATION CANNOT BE PERFORMED IN THIS
		CASE. IF X0 .LE. 0, THIS CAN OCCUR ONLY
		WHEN P OR Q IS APPROXIMATELY 0. IF X0 IS
		POSITIVE THEN THIS CAN OCCUR WHEN A IS
		EXCEEDINGLY CLOSE TO X AND A IS EXTREMELY
		LARGE (SAY A .GE. 1.E20).
		----------------------------------------------------------------------
		WRITTEN BY ALFRED H. MORRIS, JR.
		NAVAL SURFACE WEAPONS CENTER
		DAHLGREN, VIRGINIA
		-------------------
		*/
	{
		static double a0 = 3.31125922108741e0;
		static double a1 = 11.6616720288968e0;
		static double a2 = 4.28342155967104e0;
		static double a3 = .213623493715853e0;
		static double b1 = 6.61053765625462e0;
		static double b2 = 6.40691597760039e0;
		static double b3 = 1.27364489782223e0;
		static double b4 = .036117081018842e0;
		static double c = .577215664901533e0;
		static double ln10 = 2.302585e0;
		static double tol = 1.e-5;
		static double amin[2] = {
			500.0e0,100.0e0
		};
		static double bmin[2] = {
			1.e-28,1.e-13
		};
		static double dmin[2] = {
			1.e-06,1.e-04
		};
		static double emin[2] = {
			2.e-03,6.e-03
		};
		static double eps0[2] = {
			1.e-10,1.e-08
		};
		static int K1 = 1;
		static int K2 = 2;
		static int K3 = 3;
		static int K8 = 0;
		double am1,amax,ap1,ap2,ap3,apn,b,c1,c2,c3,c4,c5,d,e,e2,eps,g,h,pn,qg,qn,
			r,rta,s,s2,sum,t,u,w,_xmax,_xmin,xn,y,z;
		int iop;
		double T4,T5,T6,T7,T9;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		****** E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
		E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
		XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
		LARGEST POSITIVE NUMBER.
		*/
		e = x_spmpar(K1);
		_xmin = x_spmpar(K2);
		_xmax = x_spmpar(K3);
		*x = 0.0e0;
		if(*a <= 0.0e0) goto S300;
		t = *p+*q-1.e0;
		if(fabs(t) > e) goto S320;
		*ierr = 0;
		if(*p == 0.0e0) return;
		if(*q == 0.0e0) goto S270;
		if(*a == 1.0e0) goto S280;
		e2 = 2.0e0*e;
		amax = 0.4e-10/(e*e);
		iop = 1;
		if(e > 1.e-10) iop = 2;
		eps = eps0[iop-1];
		xn = *x0;
		if(*x0 > 0.0e0) goto S160;
		/*
		SELECTION OF THE INITIAL APPROXIMATION XN OF X
		WHEN A .LT. 1
		*/
		if(*a > 1.0e0) goto S80;
		T4 = *a+1.0e0;
		g = x_Xgamm(&T4);
		qg = *q*g;
		if(qg == 0.0e0) goto S360;
		b = qg/ *a;
		if(qg > 0.6e0**a) goto S40;
		if(*a >= 0.30e0 || b < 0.35e0) goto S10;
		t = exp(-(b+c));
		u = t*exp(t);
		xn = t*exp(u);
		goto S160;
	S10:
		if(b >= 0.45e0) goto S40;
		if(b == 0.0e0) goto S360;
		y = -log(b);
		s = 0.5e0+(0.5e0-*a);
		z = log(y);
		t = y-s*z;
		if(b < 0.15e0) goto S20;
		xn = y-s*log(t)-log(1.0e0+s/(t+1.0e0));
		goto S220;
	S20:
		if(b <= 0.01e0) goto S30;
		u = ((t+2.0e0*(3.0e0-*a))*t+(2.0e0-*a)*(3.0e0-*a))/((t+(5.0e0-*a))*t+2.0e0);
		xn = y-s*log(t)-log(u);
		goto S220;
	S30:
		c1 = -(s*z);
		c2 = -(s*(1.0e0+c1));
		c3 = s*((0.5e0*c1+(2.0e0-*a))*c1+(2.5e0-1.5e0**a));
		c4 = -(s*(((c1/3.0e0+(2.5e0-1.5e0**a))*c1+((*a-6.0e0)**a+7.0e0))*c1+(
			(11.0e0**a-46.0)**a+47.0e0)/6.0e0));
		c5 = -(s*((((-(c1/4.0e0)+(11.0e0**a-17.0e0)/6.0e0)*c1+((-(3.0e0**a)+13.0e0)*
			*a-13.0e0))*c1+0.5e0*(((2.0e0**a-25.0e0)**a+72.0e0)**a-61.0e0))*c1+((
			(25.0e0**a-195.0e0)**a+477.0e0)**a-379.0e0)/12.0e0));
		xn = (((c5/y+c4)/y+c3)/y+c2)/y+c1+y;
		if(*a > 1.0e0) goto S220;
		if(b > bmin[iop-1]) goto S220;
		*x = xn;
		return;
	S40:
		if(b**q > 1.e-8) goto S50;
		xn = exp(-(*q/ *a+c));
		goto S70;
	S50:
		if(*p <= 0.9e0) goto S60;
		T5 = -*q;
		xn = exp((x_alnrel(T5)+x_gamln1(*a))/ *a);
		goto S70;
	S60:
		xn = exp(log(*p*g)/ *a);
	S70:
		if(xn == 0.0e0) goto S310;
		t = 0.5e0+(0.5e0-xn/(*a+1.0e0));
		xn /= t;
		goto S160;
	S80:
		/*
		SELECTION OF THE INITIAL APPROXIMATION XN OF X
		WHEN A .GT. 1
		*/
		if(*q <= 0.5e0) goto S90;
		w = log(*p);
		goto S100;
	S90:
		w = log(*q);
	S100:
		t = sqrt(-(2.0e0*w));
		s = t-(((a3*t+a2)*t+a1)*t+a0)/((((b4*t+b3)*t+b2)*t+b1)*t+1.0e0);
		if(*q > 0.5e0) s = -s;
		rta = sqrt(*a);
		s2 = s*s;
		xn = *a+s*rta+(s2-1.0e0)/3.0e0+s*(s2-7.0e0)/(36.0e0*rta)-((3.0e0*s2+7.0e0)*
			s2-16.0e0)/(810.0e0**a)+s*((9.0e0*s2+256.0e0)*s2-433.0e0)/(38880.0e0**a*
				rta);
		xn = x_fifdmax1(xn,0.0e0);
		if(*a < amin[iop-1]) goto S110;
		*x = xn;
		d = 0.5e0+(0.5e0-*x/ *a);
		if(fabs(d) <= dmin[iop-1]) return;
	S110:
		if(*p <= 0.5e0) goto S130;
		if(xn < 3.0e0**a) goto S220;
		y = -(w+x_gamln(*a));
		d = x_fifdmax1(2.0e0,*a*(*a-1.0e0));
		if(y < ln10*d) goto S120;
		s = 1.0e0-*a;
		z = log(y);
		goto S30;
	S120:
		t = *a-1.0e0;
		T6 = -(t/(xn+1.0e0));
		xn = y+t*log(xn)-x_alnrel(T6);
		T7 = -(t/(xn+1.0e0));
		xn = y+t*log(xn)-x_alnrel(T7);
		goto S220;
	S130:
		ap1 = *a+1.0e0;
		if(xn > 0.70e0*ap1) goto S170;
		w += x_gamln(ap1);
		if(xn > 0.15e0*ap1) goto S140;
		ap2 = *a+2.0e0;
		ap3 = *a+3.0e0;
		*x = exp((w+*x)/ *a);
		*x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2)))/ *a);
		*x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2)))/ *a);
		*x = exp((w+*x-log(1.0e0+*x/ap1*(1.0e0+*x/ap2*(1.0e0+*x/ap3))))/ *a);
		xn = *x;
		if(xn > 1.e-2*ap1) goto S140;
		if(xn <= emin[iop-1]*ap1) return;
		goto S170;
	S140:
		apn = ap1;
		t = xn/apn;
		sum = 1.0e0+t;
	S150:
		apn += 1.0e0;
		t *= (xn/apn);
		sum += t;
		if(t > 1.e-4) goto S150;
		t = w-log(sum);
		xn = exp((xn+t)/ *a);
		xn *= (1.0e0-(*a*log(xn)-xn-t)/(*a-xn));
		goto S170;
	S160:
		/*
		SCHRODER ITERATION USING P
		*/
		if(*p > 0.5e0) goto S220;
	S170:
		if(*p <= 1.e10*_xmin) goto S350;
		am1 = *a-0.5e0-0.5e0;
	S180:
		if(*a <= amax) goto S190;
		d = 0.5e0+(0.5e0-xn/ *a);
		if(fabs(d) <= e2) goto S350;
	S190:
		if(*ierr >= 20) goto S330;
		*ierr += 1;
		x_gratio(a,&xn,&pn,&qn,&K8);
		if(pn == 0.0e0 || qn == 0.0e0) goto S350;
		r = x_rcomp(a,&xn);
		if(r == 0.0e0) goto S350;
		t = (pn-*p)/r;
		w = 0.5e0*(am1-xn);
		if(fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S200;
		*x = xn*(1.0e0-t);
		if(*x <= 0.0e0) goto S340;
		d = fabs(t);
		goto S210;
	S200:
		h = t*(1.0e0+w*t);
		*x = xn*(1.0e0-h);
		if(*x <= 0.0e0) goto S340;
		if(fabs(w) >= 1.0e0 && fabs(w)*t*t <= eps) return;
		d = fabs(h);
	S210:
		xn = *x;
		if(d > tol) goto S180;
		if(d <= eps) return;
		if(fabs(*p-pn) <= tol**p) return;
		goto S180;
	S220:
		/*
		SCHRODER ITERATION USING Q
		*/
		if(*q <= 1.e10*_xmin) goto S350;
		am1 = *a-0.5e0-0.5e0;
	S230:
		if(*a <= amax) goto S240;
		d = 0.5e0+(0.5e0-xn/ *a);
		if(fabs(d) <= e2) goto S350;
	S240:
		if(*ierr >= 20) goto S330;
		*ierr += 1;
		x_gratio(a,&xn,&pn,&qn,&K8);
		if(pn == 0.0e0 || qn == 0.0e0) goto S350;
		r = x_rcomp(a,&xn);
		if(r == 0.0e0) goto S350;
		t = (*q-qn)/r;
		w = 0.5e0*(am1-xn);
		if(fabs(t) <= 0.1e0 && fabs(w*t) <= 0.1e0) goto S250;
		*x = xn*(1.0e0-t);
		if(*x <= 0.0e0) goto S340;
		d = fabs(t);
		goto S260;
	S250:
		h = t*(1.0e0+w*t);
		*x = xn*(1.0e0-h);
		if(*x <= 0.0e0) goto S340;
		if(fabs(w) >= 1.0e0 && fabs(w)*t*t <= eps) return;
		d = fabs(h);
	S260:
		xn = *x;
		if(d > tol) goto S230;
		if(d <= eps) return;
		if(fabs(*q-qn) <= tol**q) return;
		goto S230;
	S270:
		/*
		SPECIAL CASES
		*/
		*x = _xmax;
		return;
	S280:
		if(*q < 0.9e0) goto S290;
		T9 = -*p;
		*x = -x_alnrel(T9);
		return;
	S290:
		*x = -log(*q);
		return;
	S300:
		/*
		ERROR RETURN
		*/
		*ierr = -2;
		return;
	S310:
		*ierr = -3;
		return;
	S320:
		*ierr = -4;
		return;
	S330:
		*ierr = -6;
		return;
	S340:
		*ierr = -7;
		return;
	S350:
		*x = xn;
		*ierr = -8;
		return;
	S360:
		*x = _xmax;
		*ierr = -8;
		return;
	}
	double	x_gamln(double a)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
		-----------------------------------------------------------------------
		WRITTEN BY ALFRED H. MORRIS
		NAVAL SURFACE WARFARE CENTER
		DAHLGREN, VIRGINIA
		--------------------------
		D = 0.5*(LN(2*PI) - 1)
		--------------------------
		*/
	{
		static double c0 = .833333333333333e-01;
		static double c1 = -.277777777760991e-02;
		static double c2 = .793650666825390e-03;
		static double c3 = -.595202931351870e-03;
		static double c4 = .837308034031215e-03;
		static double c5 = -.165322962780713e-02;
		static double d = .418938533204673;
		double t,w;
		int i,n;
		double T1;
		/*
		..
		.. Executable Statements ..
		*/
		if(a > 0.8) goto S10;
		return x_gamln1(a)-log(a);
	S10:
		if(a > 2.25) goto S20;
		return x_gamln1(a-1.0);
	S20:
		if(a >= 10.0e0) goto S40;
		n = int(a-1.25);
		t = a;
		w = 1.0;
		for (i=1; i<=n; i++) {
			t -= 1.0;
			w = t*w;
		}
		T1 = t-1.0;
		return x_gamln1(T1)+log(w);
	S40:
		t = SQR(1.0 / a);
		w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/ a;
		return d+w+(a-0.5)*(log(a)-1.0);
	}
	double	x_gamln1(double a)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
		-----------------------------------------------------------------------
		*/
	{
		static double p0 = .577215664901533e+00;
		static double p1 = .844203922187225e+00;
		static double p2 = -.168860593646662e+00;
		static double p3 = -.780427615533591e+00;
		static double p4 = -.402055799310489e+00;
		static double p5 = -.673562214325671e-01;
		static double p6 = -.271935708322958e-02;
		static double q1 = .288743195473681e+01;
		static double q2 = .312755088914843e+01;
		static double q3 = .156875193295039e+01;
		static double q4 = .361951990101499e+00;
		static double q5 = .325038868253937e-01;
		static double q6 = .667465618796164e-03;
		static double r0 = .422784335098467e+00;
		static double r1 = .848044614534529e+00;
		static double r2 = .565221050691933e+00;
		static double r3 = .156513060486551e+00;
		static double r4 = .170502484022650e-01;
		static double r5 = .497958207639485e-03;
		static double s1 = .124313399877507e+01;
		static double s2 = .548042109832463e+00;
		static double s3 = .101552187439830e+00;
		static double s4 = .713309612391000e-02;
		static double s5 = .116165475989616e-03;
		double w,x;
		/*
		..
		.. Executable Statements ..
		*/
		if(a >= 0.6) goto S10;
		w = ((((((p6*a+p5)*a+p4)*a+p3)*a+p2)*a+p1)*a+p0)/((((((q6*a+q5)*a+
			q4)*a+q3)*a+q2)*a+q1)*a+1.0);
		return -(a*w);
	S10:
		x = a-1.0;
		w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0)/(((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x+1.0);
		return x*w;
	}
	double	x_Xgamm(double *a)
		/*
		-----------------------------------------------------------------------

		EVALUATION OF THE GAMMA FUNCTION FOR REAL ARGUMENTS

		-----------

		GAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
		BE COMPUTED.

		-----------------------------------------------------------------------
		WRITTEN BY ALFRED H. MORRIS, JR.
		NAVAL SURFACE WEAPONS CENTER
		DAHLGREN, VIRGINIA
		-----------------------------------------------------------------------
		*/
	{
		static const double d = .41893853320467274178e0;
		static const double pi = 3.1415926535898e0;
		static const double r1 = .820756370353826e-03;
		static const double r2 = -.595156336428591e-03;
		static const double r3 = .793650663183693e-03;
		static const double r4 = -.277777777770481e-02;
		static const double r5 = .833333333333333e-01;
		static const double p[7] = {
			.539637273585445e-03,.261939260042690e-02,.204493667594920e-01,
			.730981088720487e-01,.279648642639792e+00,.553413866010467e+00,1.0e0
		};
		static const double q[7] = {
			-.832979206704073e-03,.470059485860584e-02,.225211131035340e-01,
			-.170458969313360e+00,-.567902761974940e-01,.113062953091122e+01,1.0e0
		};
		static const int K2 = 3;
		static const int K3 = 0;
		double Xgamm,bot,g,lnx,t,top,w,x,z;
		double s = numeric_limits<double>::quiet_NaN();
		int i,j,m,n,T1;
		/*
		..
		.. Executable Statements ..
		*/
		Xgamm = 0.0e0;
		x = *a;
		if(fabs(*a) >= 15.0e0) {
			/*
			-----------------------------------------------------------------------
			EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
			-----------------------------------------------------------------------
			*/
			if(fabs(*a) >= 1.e3) return Xgamm;
			if(*a <= 0.0e0) {
				x = -*a;
				n = (int)x;
				t = x-(double)n;
				if(t > 0.9e0) t = 1.0e0-t;
				s = sin(pi*t)/pi;
				if(x_fifmod(n,2) == 0) s = -s;
				if(s == 0.0e0) return Xgamm;
			}
			/*
			COMPUTE THE MODIFIED ASYMPTOTIC SUM
			*/
			t = 1.0e0/(x*x);
			g = ((((r1*t+r2)*t+r3)*t+r4)*t+r5)/x;
			/*
			ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
			BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.
			*/
			lnx = log(x);
			/*
			FINAL ASSEMBLY
			*/
			z = x;
			g = d+g+(z-0.5e0)*(lnx-1.e0);
			w = g;
			t = g-w;
			if(w > 0.99999e0*x_exparg(0)) return Xgamm;
			Xgamm = exp(w)*(1.0e0+t);
			if(*a < 0.0e0) {
				if (s != s)
					halt("Invalid value of s : Should not be NaN but it is");
				Xgamm = 1.0e0/(Xgamm*s)/x;
			}
			return Xgamm;
		}
		/*
		-----------------------------------------------------------------------
		EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15
		-----------------------------------------------------------------------
		*/
		t = 1.0e0;
		m = x_fifidint(*a)-1;
		/*
		LET T BE THE PRODUCT OF A-J WHEN A .GE. 2
		*/
		T1 = m;
		if(T1 < 0) {
			/*
			LET T BE THE PRODUCT OF A+J WHEN A .LT. 1
			*/
			t = *a;
			if(*a <= 0.0e0) {
				m = -m-1;
				if (m != 0) {
					for(j=1; j<=m; j++) {
						x += 1.0e0;
						t = x*t;
					}
				}
				x += (0.5e0+0.5e0);
				t = x*t;
				if(t == 0.0e0) return Xgamm;
			}
			/*
			THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
			CODE MAY BE OMITTED IF DESIRED.
			*/
			if(fabs(t) < 1.e-30) {
				if(fabs(t)*x_spmpar(K2) <= 1.0001e0)
					return Xgamm;
				return 1.0/t;
			}
		} else {
			if(T1 != 0) {
				for(j=1; j<=m; j++) {
					x -= 1.0;
					t = x*t;
				}
			}
			x -= 1.0e0;
		}
		/*
		COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1
		*/
		top = p[0];
		bot = q[0];
		for(i=1; i<7; i++) {
			top = p[i]+x*top;
			bot = q[i]+x*bot;
		}
		Xgamm = top/bot;
		/*
		TERMINATION
		*/
		if(*a < 1.0e0)
			Xgamm /= t;
		else
			Xgamm *= t;
		return Xgamm;
	}
	void	x_grat1(double *a,double *x,double *r,double *p,double *q,
		double *eps)
	{
		static int K2 = 0;
		double a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma,g,h,j,l,sum,t,tol,w,z,T1,T3;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		-----------------------------------------------------------------------
		EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
		P(A,X) AND Q(A,X)
		IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED.
		THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A).
		-----------------------------------------------------------------------
		*/
		if(*a**x == 0.0e0) goto S120;
		if(*a == 0.5e0) goto S100;
		if(*x < 1.1e0) goto S10;
		goto S60;
	S10:
		/*
		TAYLOR SERIES FOR P(A,X)/X**A
		*/
		an = 3.0e0;
		c = *x;
		sum = *x/(*a+3.0e0);
		tol = 0.1e0**eps/(*a+1.0e0);
	S20:
		an += 1.0e0;
		c = -(c*(*x/an));
		t = c/(*a+an);
		sum += t;
		if(fabs(t) > tol) goto S20;
		j = *a**x*((sum/6.0e0-0.5e0/(*a+2.0e0))**x+1.0e0/(*a+1.0e0));
		z = *a*log(*x);
		h = x_gam1(*a);
		g = 1.0e0+h;
		if(*x < 0.25e0) goto S30;
		if(*a < *x/2.59e0) goto S50;
		goto S40;
	S30:
		if(z > -.13394e0) goto S50;
	S40:
		w = exp(z);
		*p = w*g*(0.5e0+(0.5e0-j));
		*q = 0.5e0+(0.5e0-*p);
		return;
	S50:
		l = x_rexp(&z);
		w = 0.5e0+(0.5e0+l);
		*q = (w*j-l)*g-h;
		if(*q < 0.0e0) goto S90;
		*p = 0.5e0+(0.5e0-*q);
		return;
	S60:
		/*
		CONTINUED FRACTION EXPANSION
		*/
		a2nm1 = a2n = 1.0e0;
		b2nm1 = *x;
		b2n = *x+(1.0e0-*a);
		c = 1.0e0;
	S70:
		a2nm1 = *x*a2n+c*a2nm1;
		b2nm1 = *x*b2n+c*b2nm1;
		am0 = a2nm1/b2nm1;
		c += 1.0e0;
		cma = c-*a;
		a2n = a2nm1+cma*a2n;
		b2n = b2nm1+cma*b2n;
		an0 = a2n/b2n;
		if(fabs(an0-am0) >= *eps*an0) goto S70;
		*q = *r*an0;
		*p = 0.5e0+(0.5e0-*q);
		return;
	S80:
		/*
		SPECIAL CASES
		*/
		*p = 0.0e0;
		*q = 1.0e0;
		return;
	S90:
		*p = 1.0e0;
		*q = 0.0e0;
		return;
	S100:
		if(*x >= 0.25e0) goto S110;
		T1 = sqrt(*x);
		*p = x_erf1(&T1);
		*q = 0.5e0+(0.5e0-*p);
		return;
	S110:
		T3 = sqrt(*x);
		*q = x_erfc1(&K2,&T3);
		*p = 0.5e0+(0.5e0-*q);
		return;
	S120:
		if(*x <= *a) goto S80;
		goto S90;
	}
	void	x_gratio(double *a,double *x,double *ans,double *qans,int *ind)
		/*
		----------------------------------------------------------------------
		EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
		P(A,X) AND Q(A,X)

		----------

		IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X
		ARE NOT BOTH 0.

		ANS AND QANS ARE VARIABLES. GRATIO ASSIGNS ANS THE VALUE
		P(A,X) AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER.
		IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS
		POSSIBLE (UP TO 14 SIGNIFICANT DIGITS). OTHERWISE, IF
		IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE
		6-TH SIGNIFICANT DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY
		IS REQUESTED TO WITHIN 1 UNIT OF THE 3RD SIGNIFICANT DIGIT.

		ERROR RETURN ...
		ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE,
		WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT.
		P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN
		X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE.
		----------------------------------------------------------------------
		WRITTEN BY ALFRED H. MORRIS, JR.
		NAVAL SURFACE WEAPONS CENTER
		DAHLGREN, VIRGINIA
		--------------------
		*/
	{
		static double alog10 = 2.30258509299405e0;
		static double d10 = -.185185185185185e-02;
		static double d20 = .413359788359788e-02;
		static double d30 = .649434156378601e-03;
		static double d40 = -.861888290916712e-03;
		static double d50 = -.336798553366358e-03;
		static double d60 = .531307936463992e-03;
		static double d70 = .344367606892378e-03;
		static double rt2pin = .398942280401433e0;
		static double rtpi = 1.77245385090552e0;
		static double third = .333333333333333e0;
		static double acc0[3] = {
			5.e-15,5.e-7,5.e-4
		};
		static double big[3] = {
			20.0e0,14.0e0,10.0e0
		};
		static double d0[13] = {
			.833333333333333e-01,-.148148148148148e-01,.115740740740741e-02,
			.352733686067019e-03,-.178755144032922e-03,.391926317852244e-04,
			-.218544851067999e-05,-.185406221071516e-05,.829671134095309e-06,
			-.176659527368261e-06,.670785354340150e-08,.102618097842403e-07,
			-.438203601845335e-08
		};
		static double d1[12] = {
			-.347222222222222e-02,.264550264550265e-02,-.990226337448560e-03,
			.205761316872428e-03,-.401877572016461e-06,-.180985503344900e-04,
			.764916091608111e-05,-.161209008945634e-05,.464712780280743e-08,
			.137863344691572e-06,-.575254560351770e-07,.119516285997781e-07
		};
		static double d2[10] = {
			-.268132716049383e-02,.771604938271605e-03,.200938786008230e-05,
			-.107366532263652e-03,.529234488291201e-04,-.127606351886187e-04,
			.342357873409614e-07,.137219573090629e-05,-.629899213838006e-06,
			.142806142060642e-06
		};
		static double d3[8] = {
			.229472093621399e-03,-.469189494395256e-03,.267720632062839e-03,
			-.756180167188398e-04,-.239650511386730e-06,.110826541153473e-04,
			-.567495282699160e-05,.142309007324359e-05
		};
		static double d4[6] = {
			.784039221720067e-03,-.299072480303190e-03,-.146384525788434e-05,
			.664149821546512e-04,-.396836504717943e-04,.113757269706784e-04
		};
		static double d5[4] = {
			-.697281375836586e-04,.277275324495939e-03,-.199325705161888e-03,
			.679778047793721e-04
		};
		static double d6[2] = {
			-.592166437353694e-03,.270878209671804e-03
		};
		static double e00[3] = {
			.25e-3,.25e-1,.14e0
		};
		static double x00[3] = {
			31.0e0,17.0e0,9.7e0
		};
		int K1 = 1;
		int K2 = 0;
		double a2n,a2nm1,acc,am0,amn,an,an0,apn,b2n,b2nm1,c,c0,c1,c2,c3,c4,c5,c6,
			cma,e,e0,g,h,j,l,r,rta,rtx,s,sum,t,t1,tol,twoa,u,w,x0,y,z;
		int i,iop,m,max,n;
		double wk[20],T3;
		int T4,T5;
		double T6,T7;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		--------------------
		****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
		FLOATING POINT NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
		*/
		e = x_spmpar(K1);
		if(*a < 0.0e0 || *x < 0.0e0) goto S430;
		if(*a == 0.0e0 && *x == 0.0e0) goto S430;
		if(*a**x == 0.0e0) goto S420;
		iop = *ind+1;
		if(iop != 1 && iop != 2) iop = 3;
		acc = x_fifdmax1(acc0[iop-1],e);
		e0 = e00[iop-1];
		x0 = x00[iop-1];
		/*
		SELECT THE APPROPRIATE ALGORITHM
		*/
		if(*a >= 1.0e0) goto S10;
		if(*a == 0.5e0) goto S390;
		if(*x < 1.1e0) goto S160;
		t1 = *a*log(*x)-*x;
		u = *a*exp(t1);
		if(u == 0.0e0) goto S380;
		r = u*(1.0e0+x_gam1(*a));
		goto S250;
	S10:
		if(*a >= big[iop-1]) goto S30;
		if(*a > *x || *x >= x0) goto S20;
		twoa = *a+*a;
		m = x_fifidint(twoa);
		if(twoa != (double)m) goto S20;
		i = m/2;
		if(*a == (double)i) goto S210;
		goto S220;
	S20:
		t1 = *a*log(*x)-*x;
		r = exp(t1)/x_Xgamm(a);
		goto S40;
	S30:
		l = *x/ *a;
		if(l == 0.0e0) goto S370;
		s = 0.5e0+(0.5e0-l);
		z = x_rlog(&l);
		if(z >= 700.0e0/ *a) goto S410;
		y = *a*z;
		rta = sqrt(*a);
		if(fabs(s) <= e0/rta) goto S330;
		if(fabs(s) <= 0.4e0) goto S270;
		t = pow(1.0e0/ *a,2.0);
		t1 = (((0.75e0*t-1.0e0)*t+3.5e0)*t-105.0e0)/(*a*1260.0e0);
		t1 -= y;
		r = rt2pin*rta*exp(t1);
	S40:
		if(r == 0.0e0) goto S420;
		if(*x <= x_fifdmax1(*a,alog10)) goto S50;
		if(*x < x0) goto S250;
		goto S100;
	S50:
		/*
		TAYLOR SERIES FOR P/R
		*/
		apn = *a+1.0e0;
		t = *x/apn;
		wk[0] = t;
		for(n=2; n<=20; n++) {
			apn += 1.0e0;
			t *= (*x/apn);
			if(t <= 1.e-3) goto S70;
			wk[n-1] = t;
		}
		n = 20;
	S70:
		sum = t;
		tol = 0.5e0*acc;
	S80:
		apn += 1.0e0;
		t *= (*x/apn);
		sum += t;
		if(t > tol) goto S80;
		max = n-1;
		for(m=1; m<=max; m++) {
			n -= 1;
			sum += wk[n-1];
		}
		*ans = r/ *a*(1.0e0+sum);
		*qans = 0.5e0+(0.5e0-*ans);
		return;
	S100:
		/*
		ASYMPTOTIC EXPANSION
		*/
		amn = *a-1.0e0;
		t = amn/ *x;
		wk[0] = t;
		for(n=2; n<=20; n++) {
			amn -= 1.0e0;
			t *= (amn/ *x);
			if(fabs(t) <= 1.e-3) goto S120;
			wk[n-1] = t;
		}
		n = 20;
	S120:
		sum = t;
	S130:
		if(fabs(t) <= acc) goto S140;
		amn -= 1.0e0;
		t *= (amn/ *x);
		sum += t;
		goto S130;
	S140:
		max = n-1;
		for(m=1; m<=max; m++) {
			n -= 1;
			sum += wk[n-1];
		}
		*qans = r/ *x*(1.0e0+sum);
		*ans = 0.5e0+(0.5e0-*qans);
		return;
	S160:
		/*
		TAYLOR SERIES FOR P(A,X)/X**A
		*/
		an = 3.0e0;
		c = *x;
		sum = *x/(*a+3.0e0);
		tol = 3.0e0*acc/(*a+1.0e0);
	S170:
		an += 1.0e0;
		c = -(c*(*x/an));
		t = c/(*a+an);
		sum += t;
		if(fabs(t) > tol) goto S170;
		j = *a**x*((sum/6.0e0-0.5e0/(*a+2.0e0))**x+1.0e0/(*a+1.0e0));
		z = *a*log(*x);
		h = x_gam1(*a);
		g = 1.0e0+h;
		if(*x < 0.25e0) goto S180;
		if(*a < *x/2.59e0) goto S200;
		goto S190;
	S180:
		if(z > -.13394e0) goto S200;
	S190:
		w = exp(z);
		*ans = w*g*(0.5e0+(0.5e0-j));
		*qans = 0.5e0+(0.5e0-*ans);
		return;
	S200:
		l = x_rexp(&z);
		w = 0.5e0+(0.5e0+l);
		*qans = (w*j-l)*g-h;
		if(*qans < 0.0e0) goto S380;
		*ans = 0.5e0+(0.5e0-*qans);
		return;
	S210:
		/*
		FINITE SUMS FOR Q WHEN A .GE. 1
		AND 2*A IS AN INTEGER
		*/
		sum = exp(-*x);
		t = sum;
		n = 1;
		c = 0.0e0;
		goto S230;
	S220:
		rtx = sqrt(*x);
		sum = x_erfc1(&K2,&rtx);
		t = exp(-*x)/(rtpi*rtx);
		n = 0;
		c = -0.5e0;
	S230:
		if(n == i) goto S240;
		n += 1;
		c += 1.0e0;
		t = *x*t/c;
		sum += t;
		goto S230;
	S240:
		*qans = sum;
		*ans = 0.5e0+(0.5e0-*qans);
		return;
	S250:
		/*
		CONTINUED FRACTION EXPANSION
		*/
		tol = x_fifdmax1(5.0e0*e,acc);
		a2nm1 = a2n = 1.0e0;
		b2nm1 = *x;
		b2n = *x+(1.0e0-*a);
		c = 1.0e0;
	S260:
		a2nm1 = *x*a2n+c*a2nm1;
		b2nm1 = *x*b2n+c*b2nm1;
		am0 = a2nm1/b2nm1;
		c += 1.0e0;
		cma = c-*a;
		a2n = a2nm1+cma*a2n;
		b2n = b2nm1+cma*b2n;
		an0 = a2n/b2n;
		if(fabs(an0-am0) >= tol*an0) goto S260;
		*qans = r*an0;
		*ans = 0.5e0+(0.5e0-*qans);
		return;
	S270:
		/*
		GENERAL TEMME EXPANSION
		*/
		if(fabs(s) <= 2.0e0*e && *a*e*e > 3.28e-3) goto S430;
		c = exp(-y);
		T3 = sqrt(y);
		w = 0.5e0*x_erfc1(&K1,&T3);
		u = 1.0e0/ *a;
		z = sqrt(z+z);
		if(l < 1.0e0) z = -z;
		T4 = iop-2;
		if(T4 < 0) goto S280;
		else if(T4 == 0) goto S290;
		else  goto S300;
	S280:
		if(fabs(s) <= 1.e-3) goto S340;
		c0 = ((((((((((((d0[12]*z+d0[11])*z+d0[10])*z+d0[9])*z+d0[8])*z+d0[7])*z+d0[
			6])*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
		c1 = (((((((((((d1[11]*z+d1[10])*z+d1[9])*z+d1[8])*z+d1[7])*z+d1[6])*z+d1[5]
			)*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
		c2 = (((((((((d2[9]*z+d2[8])*z+d2[7])*z+d2[6])*z+d2[5])*z+d2[4])*z+d2[3])*z+
			d2[2])*z+d2[1])*z+d2[0])*z+d20;
		c3 = (((((((d3[7]*z+d3[6])*z+d3[5])*z+d3[4])*z+d3[3])*z+d3[2])*z+d3[1])*z+
			d3[0])*z+d30;
		c4 = (((((d4[5]*z+d4[4])*z+d4[3])*z+d4[2])*z+d4[1])*z+d4[0])*z+d40;
		c5 = (((d5[3]*z+d5[2])*z+d5[1])*z+d5[0])*z+d50;
		c6 = (d6[1]*z+d6[0])*z+d60;
		t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
		goto S310;
	S290:
		c0 = (((((d0[5]*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-third;
		c1 = (((d1[3]*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
		c2 = d2[0]*z+d20;
		t = (c2*u+c1)*u+c0;
		goto S310;
	S300:
		t = ((d0[2]*z+d0[1])*z+d0[0])*z-third;
	S310:
		if(l < 1.0e0) goto S320;
		*qans = c*(w+rt2pin*t/rta);
		*ans = 0.5e0+(0.5e0-*qans);
		return;
	S320:
		*ans = c*(w-rt2pin*t/rta);
		*qans = 0.5e0+(0.5e0-*ans);
		return;
	S330:
		/*
		TEMME EXPANSION FOR L = 1
		*/
		if(*a*e*e > 3.28e-3) goto S430;
		c = 0.5e0+(0.5e0-y);
		w = (0.5e0-sqrt(y)*(0.5e0+(0.5e0-y/3.0e0))/rtpi)/c;
		u = 1.0e0/ *a;
		z = sqrt(z+z);
		if(l < 1.0e0) z = -z;
		T5 = iop-2;
		if(T5 < 0) goto S340;
		else if(T5 == 0) goto S350;
		else  goto S360;
	S340:
		c0 = ((((((d0[6]*z+d0[5])*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z-
			third;
		c1 = (((((d1[5]*z+d1[4])*z+d1[3])*z+d1[2])*z+d1[1])*z+d1[0])*z+d10;
		c2 = ((((d2[4]*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z+d20;
		c3 = (((d3[3]*z+d3[2])*z+d3[1])*z+d3[0])*z+d30;
		c4 = (d4[1]*z+d4[0])*z+d40;
		c5 = (d5[1]*z+d5[0])*z+d50;
		c6 = d6[0]*z+d60;
		t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u+c0;
		goto S310;
	S350:
		c0 = (d0[1]*z+d0[0])*z-third;
		c1 = d1[0]*z+d10;
		t = (d20*u+c1)*u+c0;
		goto S310;
	S360:
		t = d0[0]*z-third;
		goto S310;
	S370:
		/*
		SPECIAL CASES
		*/
		*ans = 0.0e0;
		*qans = 1.0e0;
		return;
	S380:
		*ans = 1.0e0;
		*qans = 0.0e0;
		return;
	S390:
		if(*x >= 0.25e0) goto S400;
		T6 = sqrt(*x);
		*ans = x_erf1(&T6);
		*qans = 0.5e0+(0.5e0-*ans);
		return;
	S400:
		T7 = sqrt(*x);
		*qans = x_erfc1(&K2,&T7);
		*ans = 0.5e0+(0.5e0-*qans);
		return;
	S410:
		if(fabs(s) <= 2.0e0*e) goto S430;
	S420:
		if(*x <= *a) goto S370;
		goto S380;
	S430:
		/*
		ERROR RETURN
		*/
		*ans = 2.0e0;
		return;
	}
	double	x_gsumln(double *a,double *b)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
		FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
		-----------------------------------------------------------------------
		*/
	{
		double x;
		/*
		..
		.. Executable Statements ..
		*/
		x = *a+*b-2.;
		if(x > 0.25) goto S10;
		return x_gamln1(1.0+x);
	S10:
		if(x > 1.25) goto S20;
		return x_gamln1(x)+x_alnrel(x);
	S20:
		return x_gamln1(x-1.0)+log(x*(1.0+x));
	}
	double	x_psi(double xx)
		/*
		---------------------------------------------------------------------

		EVALUATION OF THE DIGAMMA FUNCTION

		-----------

		PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
		BE COMPUTED.

		THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
		APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY
		CODY, STRECOK AND THACHER.

		---------------------------------------------------------------------
		PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK
		PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY
		A.H. MORRIS (NSWC).
		---------------------------------------------------------------------
		*/
	{
		static double dx0 = 1.461632144968362341262659542325721325;
		static double piov4 = .785398163397448;
		static double p1[7] = {
			.895385022981970e-02,.477762828042627e+01,.142441585084029e+03,
			.118645200713425e+04,.363351846806499e+04,.413810161269013e+04,
			.130560269827897e+04
		};
		static double p2[4] = {
			-.212940445131011e+01,-.701677227766759e+01,-.448616543918019e+01,
			-.648157123766197e+00
		};
		static double q1[6] = {
			.448452573429826e+02,.520752771467162e+03,.221000799247830e+04,
			.364127349079381e+04,.190831076596300e+04,.691091682714533e-05
		};
		static double q2[4] = {
			.322703493791143e+02,.892920700481861e+02,.546117738103215e+02,
			.777788548522962e+01
		};
		static int K1 = 3;
		static int K2 = 1;
		double psi,aug,den,sgn,upper,w,x,xmax1,xmx0,xsmall,z;
		int i,m,n,nq;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		---------------------------------------------------------------------
		MACHINE DEPENDENT CONSTANTS ...
		XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
		WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
		AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
		ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
		PSI MAY BE REPRESENTED AS ALOG(X).
		XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
		MAY BE REPRESENTED BY 1/X.
		---------------------------------------------------------------------
		*/
		xmax1 = 2147483647;
		xmax1 = x_fifdmin1(xmax1,1.0/x_spmpar(K2));
		xsmall = 1.e-9;
		x = xx;
		aug = 0.0;
		if(x >= 0.5) goto S50;
		/*
		---------------------------------------------------------------------
		X .LT. 0.5,  USE REFLECTION FORMULA
		PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
		---------------------------------------------------------------------
		*/
		if(fabs(x) > xsmall) goto S10;
		if(x == 0.0) goto S100;
		/*
		---------------------------------------------------------------------
		0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
		FOR  PI*COTAN(PI*X)
		---------------------------------------------------------------------
		*/
		aug = -(1.0/x);
		goto S40;
	S10:
		/*
		---------------------------------------------------------------------
		REDUCTION OF ARGUMENT FOR COTAN
		---------------------------------------------------------------------
		*/
		w = -x;
		sgn = piov4;
		if(w > 0.0) goto S20;
		w = -w;
		sgn = -sgn;
	S20:
		/*
		---------------------------------------------------------------------
		MAKE AN ERROR EXIT IF X .LE. -XMAX1
		---------------------------------------------------------------------
		*/
		if(w >= xmax1) goto S100;
		nq = x_fifidint(w);
		w -= (double)nq;
		nq = x_fifidint(w*4.0);
		w = 4.0*(w-(double)nq*.25);
		/*
		---------------------------------------------------------------------
		W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
		ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
		QUADRANT AND DETERMINE SIGN
		---------------------------------------------------------------------
		*/
		n = nq/2;
		if(n+n != nq) w = 1.0-w;
		z = piov4*w;
		m = n/2;
		if(m+m != n) sgn = -sgn;
		/*
		---------------------------------------------------------------------
		DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
		---------------------------------------------------------------------
		*/
		n = (nq+1)/2;
		m = n/2;
		m += m;
		if(m != n) goto S30;
		/*
		---------------------------------------------------------------------
		CHECK FOR SINGULARITY
		---------------------------------------------------------------------
		*/
		if(z == 0.0) goto S100;
		/*
		---------------------------------------------------------------------
		USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
		SIN/COS AS A SUBSTITUTE FOR TAN
		---------------------------------------------------------------------
		*/
		aug = sgn*(cos(z)/sin(z)*4.0);
		goto S40;
	S30:
		aug = sgn*(sin(z)/cos(z)*4.0);
	S40:
		x = 1.0-x;
	S50:
		if(x > 3.0) goto S70;
		/*
		---------------------------------------------------------------------
		0.5 .LE. X .LE. 3.0
		---------------------------------------------------------------------
		*/
		den = x;
		upper = p1[0]*x;
		for(i=1; i<=5; i++) {
			den = (den+q1[i-1])*x;
			upper = (upper+p1[i+1-1])*x;
		}
		den = (upper+p1[6])/(den+q1[5]);
		xmx0 = x-dx0;
		psi = den*xmx0+aug;
		return psi;
	S70:
		/*
		---------------------------------------------------------------------
		IF X .GE. XMAX1, PSI = LN(X)
		---------------------------------------------------------------------
		*/
		if(x >= xmax1) goto S90;
		/*
		---------------------------------------------------------------------
		3.0 .LT. X .LT. XMAX1
		---------------------------------------------------------------------
		*/
		w = 1.0/(x*x);
		den = w;
		upper = p2[0]*w;
		for(i=1; i<=3; i++) {
			den = (den+q2[i-1])*w;
			upper = (upper+p2[i+1-1])*w;
		}
		aug = upper/(den+q2[3])-0.5/x+aug;
	S90:
		psi = aug+log(x);
		return psi;
	S100:
		/*
		---------------------------------------------------------------------
		ERROR RETURN
		---------------------------------------------------------------------
		*/
		return 0.0;
	}
	double	x_rcomp(double *a,double *x)
		/*
		-------------------
		EVALUATION OF EXP(-X)*X**A/GAMMA(A)
		-------------------
		RT2PIN = 1/SQRT(2*PI)
		-------------------
		*/
	{
		static double rt2pin = .398942280401433e0;
		double _rcomp,t,t1,u;
		/*
		..
		.. Executable Statements ..
		*/
		_rcomp = 0.0e0;
		if(*a >= 20.0e0) goto S20;
		t = *a*log(*x)-*x;
		if(*a >= 1.0e0) goto S10;
		_rcomp = *a*exp(t)*(1.0e0+x_gam1(*a));
		return _rcomp;
	S10:
		_rcomp = exp(t)/x_Xgamm(a);
		return _rcomp;
	S20:
		u = *x/ *a;
		if(u == 0.0e0) return _rcomp;
		t = SQR(1.0e0/ *a);
		t1 = (((0.75e0*t-1.0e0)*t+3.5e0)*t-105.0e0)/(*a*1260.0e0);
		t1 -= (*a*x_rlog(&u));
		_rcomp = rt2pin*sqrt(*a)*exp(t1);
		return _rcomp;
	}
	double	x_rexp(double *x)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF THE FUNCTION EXP(X) - 1
		-----------------------------------------------------------------------
		*/
	{
		static double p1 = .914041914819518e-09;
		static double p2 = .238082361044469e-01;
		static double q1 = -.499999999085958e+00;
		static double q2 = .107141568980644e+00;
		static double q3 = -.119041179760821e-01;
		static double q4 = .595130811860248e-03;
		double rexp,w;
		/*
		..
		.. Executable Statements ..
		*/
		if(fabs(*x) > 0.15e0) goto S10;
		rexp = *x*(((p2**x+p1)**x+1.0e0)/((((q4**x+q3)**x+q2)**x+q1)**x+1.0e0));
		return rexp;
	S10:
		w = exp(*x);
		if(*x > 0.0e0) goto S20;
		rexp = w-1.0;
		return rexp;
	S20:
		rexp = w*(0.5e0+(0.5e0-1.0e0/w));
		return rexp;
	}
	double	x_rlog(double *x)
		/*
		-------------------
		COMPUTATION OF  X - 1 - LN(X)
		-------------------
		*/
	{
		static double a = .566749439387324e-01;
		static double b = .456512608815524e-01;
		static double p0 = .333333333333333e+00;
		static double p1 = -.224696413112536e+00;
		static double p2 = .620886815375787e-02;
		static double q1 = -.127408923933623e+01;
		static double q2 = .354508718369557e+00;
		double rlog,r,t,u,w,w1;
		/*
		..
		.. Executable Statements ..
		*/
		if(*x < 0.61e0 || *x > 1.57e0) goto S40;
		if(*x < 0.82e0) goto S10;
		if(*x > 1.18e0) goto S20;
		/*
		ARGUMENT REDUCTION
		*/
		u = *x-0.5e0-0.5e0;
		w1 = 0.0e0;
		goto S30;
	S10:
		u = *x-0.7e0;
		u /= 0.7e0;
		w1 = a-u*0.3e0;
		goto S30;
	S20:
		u = 0.75e0**x-1.e0;
		w1 = b+u/3.0e0;
	S30:
		/*
		SERIES EXPANSION
		*/
		r = u/(u+2.0e0);
		t = r*r;
		w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0e0);
		rlog = 2.0e0*t*(1.0e0/(1.0e0-r)-r*w)+w1;
		return rlog;
	S40:
		r = *x-0.5e0-0.5e0;
		rlog = r-log(*x);
		return rlog;
	}
	double	x_rlog1(double *x)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF THE FUNCTION X - LN(1 + X)
		-----------------------------------------------------------------------
		*/
	{
		static double a = .566749439387324e-01;
		static double b = .456512608815524e-01;
		static double p0 = .333333333333333e+00;
		static double p1 = -.224696413112536e+00;
		static double p2 = .620886815375787e-02;
		static double q1 = -.127408923933623e+01;
		static double q2 = .354508718369557e+00;
		double rlog1,h,r,t,w,w1;
		/*
		..
		.. Executable Statements ..
		*/
		if(*x < -0.39e0 || *x > 0.57e0) goto S40;
		if(*x < -0.18e0) goto S10;
		if(*x > 0.18e0) goto S20;
		/*
		ARGUMENT REDUCTION
		*/
		h = *x;
		w1 = 0.0e0;
		goto S30;
	S10:
		h = *x+0.3e0;
		h /= 0.7e0;
		w1 = a-h*0.3e0;
		goto S30;
	S20:
		h = 0.75e0**x-0.25e0;
		w1 = b+h/3.0e0;
	S30:
		/*
		SERIES EXPANSION
		*/
		r = h/(h+2.0e0);
		t = r*r;
		w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0e0);
		rlog1 = 2.0e0*t*(1.0e0/(1.0e0-r)-r*w)+w1;
		return rlog1;
	S40:
		w = *x+0.5e0+0.5e0;
		rlog1 = *x-log(w);
		return rlog1;
	}
	/*
	SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
	THE COMPUTER BEING USED.
	*/
	double x_spmpar(int i)
	{
		if (i > 1) goto S10;
		return 0.0000000000000002220446049250313080847263336181640625; // 2^(1-53)
	S10:
		if (i > 2) goto S20;
		return 2.225073858507201383090232717332404064219215980462331e-308; // 2^(-1021+2) * 0.5 * 0.5 * 0.5
	S20:
		// z = 4503599627370496.0; // 2^(53-1)
		// w = ((z-1)*2+1)/(2*z);
		// w = 0.99999999999999988897769753748434595763683319091796875;
		// z = 4.4942328371557897693232629769725618340449424473557664e307 // 2^(1024-2)
		// spmpar = w*z*2*2;
		return 1.7976931348623157081452742373170435679807056752584500e308;
	}
	double x_stvaln(double *p)
		/*
		**********************************************************************

		double stvaln(double *p)
		STarting VALue for Neton-Raphon
		calculation of Normal distribution Inverse


		Function


		Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
		infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P


		Arguments


		P --> The probability whose normal deviate is sought.
		P is DOUBLE PRECISION


		Method


		The  rational   function   on  page 95    of Kennedy  and  Gentle,
		Statistical Computing, Marcel Dekker, NY , 1980.

		**********************************************************************
		*/
	{
		static double xden[5] = {
			0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,
			0.38560700634e-2
		};
		static double xnum[5] = {
			-0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,
			-0.453642210148e-4
		};
		double sign,y,z;
		/*
		..
		.. Executable Statements ..
		*/
		if(!(*p <= 0.5e0)) goto S10;
		sign = -1.0e0;
		z = *p;
		goto S20;
	S10:
		sign = 1.0e0;
		z = 1.0e0-*p;
	S20:
		y = sqrt(-(2.0e0*log(z)));
		return sign * (y+x_dblEvalPolyX(xnum,5,y)/x_dblEvalPolyX(xden,5,y));
	}
	/************************************************************************
	FIFDINT:
	Truncates a double precision number to an integer and returns the
	value in a double.4
	***********************************************************************/
	double x_fifdint(double a)
		// a     -     number to be truncated
	{
		return (double)((int)a);
	}
	/************************************************************************
	FIFDMAX1:
	returns the maximum of two numbers a and b
	************************************************************************/
	double x_fifdmax1(double a,double b)
		/* a     -      first number */
		/* b     -      second number */
	{
		if (a < b) return b;
		else return a;
	}
	/************************************************************************
	FIFDMIN1:
	returns the minimum of two numbers a and b
	************************************************************************/
	double x_fifdmin1(double a,double b)
		/* a     -     first number */
		/* b     -     second number */
	{
		if (a < b) return a;
		else return b;
	}
	/************************************************************************
	FIFDSIGN:
	transfers the sign of the variable "sign" to the variable "mag"
	************************************************************************/
	double x_fifdsign(double mag,double sign)
		/* mag     -     magnitude */
		/* sign    -     sign to be transfered */
	{
		if (mag < 0) mag = -mag;
		if (sign < 0) mag = -mag;
		return mag;

	}
	/************************************************************************
	FIFIDINT:
	Truncates a double precision number to a long integer
	************************************************************************/
	long	x_fifidint(double a)
		/* a - number to be truncated */
	{
		if (a < 1.0) return (long) 0;
		else return (long) a;
	}
	/************************************************************************
	FIFMOD:
	returns the modulo of a and b
	************************************************************************/
	long	x_fifmod(long a,long b)
		/* a - numerator */
		/* b - denominator */
	{
		return a % b;
	}
	double	x_stirlerr(double n)
	{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
		/*
		error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
		*/
		const static double sferr_halves[31] = {
			0.0, /* n=0 - wrong, place holder only */
			0.1534264097200273452913848,  /* 0.5 */
			0.0810614667953272582196702,  /* 1.0 */
			0.0548141210519176538961390,  /* 1.5 */
			0.0413406959554092940938221,  /* 2.0 */
			0.03316287351993628748511048, /* 2.5 */
			0.02767792568499833914878929, /* 3.0 */
			0.02374616365629749597132920, /* 3.5 */
			0.02079067210376509311152277, /* 4.0 */
			0.01848845053267318523077934, /* 4.5 */
			0.01664469118982119216319487, /* 5.0 */
			0.01513497322191737887351255, /* 5.5 */
			0.01387612882307074799874573, /* 6.0 */
			0.01281046524292022692424986, /* 6.5 */
			0.01189670994589177009505572, /* 7.0 */
			0.01110455975820691732662991, /* 7.5 */
			0.010411265261972096497478567, /* 8.0 */
			0.009799416126158803298389475, /* 8.5 */
			0.009255462182712732917728637, /* 9.0 */
			0.008768700134139385462952823, /* 9.5 */
			0.008330563433362871256469318, /* 10.0 */
			0.007934114564314020547248100, /* 10.5 */
			0.007573675487951840794972024, /* 11.0 */
			0.007244554301320383179543912, /* 11.5 */
			0.006942840107209529865664152, /* 12.0 */
			0.006665247032707682442354394, /* 12.5 */
			0.006408994188004207068439631, /* 13.0 */
			0.006171712263039457647532867, /* 13.5 */
			0.005951370112758847735624416, /* 14.0 */
			0.005746216513010115682023589, /* 14.5 */
			0.005554733551962801371038690  /* 15.0 */
		};
		double nn;

		if (n <= 15.0) {
			nn = n + n;
			if (nn == (int)nn) return sferr_halves[(int)nn];
			double np1 = n + 1.;
			return x_alngam(&np1) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI;
		}

		nn = n*n;
		if (n>500) return((S0-S1/nn)/n);
		if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
		if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
		/* 15 < n <= 35 : */
		return (S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n;
	}
	double	x_truncate(double x)
	{
		if (x >= 0) return floor(x);
		else return ceil(x);
	}
	double x_betaPraw(double x, double pin, double qin, int lower_tail, int log_p)
	{
		if (log_p) halt("log_p is not allowed");

		double w, wc;
		int ierr;
		x_bratio(pin, qin, x, 0.5 - x + 0.5, &w, &wc, &ierr); /* -> ./toms708.c */
															  /* ierr = 8 is about inaccuracy in extreme cases */
		if (ierr && (ierr != 8 || log_p))
			LOG("pbeta_raw() -> bratio() gave error code %d\n", ierr);
		return lower_tail ? w : wc;
	} /* pbeta_raw() */

	double lbeta(double a, double b)
	{
		double corr, p, q;

#ifdef IEEE_754
		if (ISNAN(a) || ISNAN(b))
			return a + b;
#endif
		p = q = a;
		if (b < p) p = b;/* := min(a,b) */
		if (b > q) q = b;/* := max(a,b) */

						 /* both arguments must be >= 0 */
		if (p < 0) numeric_limits<double>::quiet_NaN();
		else if (p == 0) {
			return ML_POSINF;
		}
		else if (q == numeric_limits<double>::infinity() ||
			q == -numeric_limits<double>::infinity()) { /* q == +Inf */
			return ML_NEGINF;
		}

		if (p >= 10) {
			/* p and q are big. */
			corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
			return log(q) * -0.5 + M_LN_SQRT_2PI + corr
				+ (p - 0.5) * log(p / (p + q)) + q * util_log1p(-p / (p + q));
		}
		else if (q >= 10) {
			/* p is small, but q is big. */
			corr = lgammacor(q) - lgammacor(p + q);
			return lgammafn(p) + corr + p - p * log(p + q)
				+ (q - 0.5) * util_log1p(-p / (p + q));
		}
		else {
			/* p and q are small: p <= q < 10. */
			/* R change for very small args */
			double pq = p + q;
			if (p < 1e-306) return x_alngam(&p) + (x_alngam(&q) - x_alngam(&pq));
		}
		return log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
	}

	double gammafn(double x)
	{
		const static double gamcs[42] = {
			+.8571195590989331421920062399942e-2,
			+.4415381324841006757191315771652e-2,
			+.5685043681599363378632664588789e-1,
			-.4219835396418560501012500186624e-2,
			+.1326808181212460220584006796352e-2,
			-.1893024529798880432523947023886e-3,
			+.3606925327441245256578082217225e-4,
			-.6056761904460864218485548290365e-5,
			+.1055829546302283344731823509093e-5,
			-.1811967365542384048291855891166e-6,
			+.3117724964715322277790254593169e-7,
			-.5354219639019687140874081024347e-8,
			+.9193275519859588946887786825940e-9,
			-.1577941280288339761767423273953e-9,
			+.2707980622934954543266540433089e-10,
			-.4646818653825730144081661058933e-11,
			+.7973350192007419656460767175359e-12,
			-.1368078209830916025799499172309e-12,
			+.2347319486563800657233471771688e-13,
			-.4027432614949066932766570534699e-14,
			+.6910051747372100912138336975257e-15,
			-.1185584500221992907052387126192e-15,
			+.2034148542496373955201026051932e-16,
			-.3490054341717405849274012949108e-17,
			+.5987993856485305567135051066026e-18,
			-.1027378057872228074490069778431e-18,
			+.1762702816060529824942759660748e-19,
			-.3024320653735306260958772112042e-20,
			+.5188914660218397839717833550506e-21,
			-.8902770842456576692449251601066e-22,
			+.1527474068493342602274596891306e-22,
			-.2620731256187362900257328332799e-23,
			+.4496464047830538670331046570666e-24,
			-.7714712731336877911703901525333e-25,
			+.1323635453126044036486572714666e-25,
			-.2270999412942928816702313813333e-26,
			+.3896418998003991449320816639999e-27,
			-.6685198115125953327792127999999e-28,
			+.1146998663140024384347613866666e-28,
			-.1967938586345134677295103999999e-29,
			+.3376448816585338090334890666666e-30,
			-.5793070335782135784625493333333e-31
		};

		int i, n;
		double y;
		double sinpiy, value;

#ifdef NOMORE_FOR_THREADS
		static int ngam = 0;
		static double xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;

		/* Initialize machine dependent constants, the first time gamma() is called.
		FIXME for threads ! */
		if (ngam == 0) {
			ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);/*was .1*d1mach(3)*/
			gammalims(&xmin, &xmax);/*-> ./gammalims.c */
			xsml = exp(fmax2(log(DBL_MIN), -log(DBL_MAX)) + 0.01);
			/*   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
			dxrel = sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)) */
		}
#else
		/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
		* (xmin, xmax) are non-trivial, see ./gammalims.c
		* xsml = exp(.01)*DBL_MIN
		* dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
		*/
# define ngam 22
# define xmin -170.5674972726612
# define xmax  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel 1.490116119384765696e-8
#endif

		if (x != x) return x;

		/* If the argument is exactly zero or a negative integer
		* then return NaN. */
		if (x == 0 || (x < 0 && x == (long)x)) {
			//ML_ERROR(ME_DOMAIN, "gammafn");
			return numeric_limits<double>::quiet_NaN();
		}

		y = fabs(x);

		if (y <= 10) {

			/* Compute gamma(x) for -10 <= x <= 10
			* Reduce the interval and find gamma(1 + y) for 0 <= y < 1
			* first of all. */

			n = (int) x;
			if(x < 0) --n;
			y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
			--n;
			value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
			if (n == 0)
				return value;/* x = 1.dddd = 1+y */

			if (n < 0) {
				/* compute gamma(x) for -10 <= x < 1 */

				/* exact 0 or "-n" checked already above */

				/* The answer is less than half precision */
				/* because x too near a negative integer. */
				if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
					//			ML_ERROR(ME_PRECISION, "gammafn");
				}

				/* The argument is so close to 0 that the result would overflow. */
				if (y < xsml) {
					//			ML_ERROR(ME_RANGE, "gammafn");
					if(x > 0) return ML_POSINF;
					else return ML_NEGINF;
				}

				n = -n;

				for (i = 0; i < n; i++) {
					value /= (x + i);
				}
				return value;
			}
			else {
				/* gamma(x) for 2 <= x <= 10 */

				for (i = 1; i <= n; i++) {
					value *= (y + i);
				}
				return value;
			}
		}
		else {
			/* gamma(x) for	 y = |x| > 10. */

			if (x > xmax) {			/* Overflow */
									// 		ML_ERROR(ME_RANGE, "gammafn");
				return ML_POSINF;
			}

			if (x < xmin) {			/* Underflow */
									//	    ML_ERROR(ME_UNDERFLOW, "gammafn");
				return 0.;
			}

			if(y <= 50 && y == (int)y) { /* compute (n - 1)! */
				value = 1.;
				for (i = 2; i < y; i++) value *= i;
			}
			else { /* normal case */
				value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
					((2*y == (int)2*y)? x_stirlerr(y) : lgammacor(y)));
			}
			if (x > 0)
				return value;

			if (fabs((x - (int)(x - 0.5))/x) < dxrel){

				/* The answer is less than half precision because */
				/* the argument is too near a negative integer. */

				//    ML_ERROR(ME_PRECISION, "gammafn");
			}

			sinpiy = sin(WISARD_PI * y);
			if (sinpiy == 0) {		/* Negative integer arg - overflow */
									//	    ML_ERROR(ME_RANGE, "gammafn");
				return ML_POSINF;
			}

			return -WISARD_PI / (y * sinpiy * value);
		}
#ifdef xmax
#	undef xmax
#endif
	}

#define fpu 3e-308
	/* acu_min:  Minimal value for accuracy 'acu' which will depend on (a,p);
	acu_min >= fpu ! */
#define acu_min 1e-300
#define lower fpu
#define upper 1-2.22e-16

#define const1 2.30753
#define const2 0.27061
#define const3 0.99229
#define const4 0.04481
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2)) */

	double lgammafn_sign(double x, int *sgn)
	{
		double ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
		static double xmax = 0.;
		static double dxrel = 0.;

		if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
			xmax = d1mach(2)/log(d1mach(2));/* = 2.533 e305	 for IEEE double */
			dxrel = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
		}
#else
		/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
		xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
		dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
		*/
#define xmax  2.5327372760800758e+305
#define dxrel 1.490116119384765696e-8
#endif

		if (sgn != NULL) *sgn = 1;

#ifdef IEEE_754
		if(ISNAN(x)) return x;
#endif

		if (x < 0 && fmod(floor(-x), 2.) == 0)
			if (sgn != NULL) *sgn = -1;

		if (x <= 0 && x == x_truncate(x)) { /* Negative integer argument */
											//		ML_ERROR(ME_RANGE, "lgamma");
			return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
		}

		y = fabs(x);

		if (y < 1e-306) return -log(x); // denormalized range, R change
		if (y <= 10) return log(fabs(gammafn(x)));
		/*
		ELSE  y = |x| > 10 ---------------------- */

		if (y > xmax) {
			//		ML_ERROR(ME_RANGE, "lgamma");
			return ML_POSINF;
		}

		if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
			if(x > 1e17)
				return(x*(log(x) - 1.));
			else if(x > 4934720.)
				return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
			else
#endif
				return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
		}
		/* else: x < -10; y = -x */
		sinpiy = fabs(sin(WISARD_PI * y));

		if (sinpiy == 0) { /* Negative integer argument ===
						   Now UNNECESSARY: caught above */
			halt(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
		}

		ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

		if(fabs((x - x_truncate(x - 0.5)) * ans / x) < dxrel) {
			/* The answer is less than half precision because
			* the argument is too near a negative integer. */

			//		ML_ERROR(ME_PRECISION, "lgamma");
		}

		return ans;
#ifdef xmax
#	undef xmax
#endif
	}

	double lgammafn(double x)
	{
		return lgammafn_sign(x, NULL);
	}

	double lgammacor(double x)
	{
		const static double algmcs[15] = {
			+.1666389480451863247205729650822e+0,
			-.1384948176067563840732986059135e-4,
			+.9810825646924729426157171547487e-8,
			-.1809129475572494194263306266719e-10,
			+.6221098041892605227126015543416e-13,
			-.3399615005417721944303330599666e-15,
			+.2683181998482698748957538846666e-17,
			-.2868042435334643284144622399999e-19,
			+.3962837061046434803679306666666e-21,
			-.6831888753985766870111999999999e-23,
			+.1429227355942498147573333333333e-24,
			-.3547598158101070547199999999999e-26,
			+.1025680058010470912000000000000e-27,
			-.3401102254316748799999999999999e-29,
			+.1276642195630062933333333333333e-30
		};

		double tmp;

		/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
		*   xbig = 2 ^ 26.5
		*   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig  94906265.62425156
#define xmax  3.745194030963158e306

		if (x < 10) numeric_limits<double>::quiet_NaN();
		else if (x >= xmax) {
			//		ML_ERROR(ME_UNDERFLOW, "lgammacor");
			/* allow to underflow below */
		}
		else if (x < xbig) {
			tmp = 10 / x;
			return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
		}
		return 1 / (x * 12);
	}

	double qbeta(double alpha, double p, double q, int lower_tail, int log_p)
	{
		int B_swapTail, i_pb, i_inn;
		double a, adj, logbeta, g, h, pp, p_, prev, qq, r, s, t, tx, w, y, yprev;
		double acu;
		volatile double xinbta;

		/* test for admissibility of parameters */

#ifdef IEEE_754
		if (ISNAN(p) || ISNAN(q) || ISNAN(alpha))
			return p + q + alpha;
#endif
		if(p < 0. || q < 0.) return numeric_limits<double>::quiet_NaN();

		R_Q_P01_boundaries(alpha, 0, 1);

		p_ = R_DT_qIv(alpha);/* lower_tail prob (in any case) */

		if(log_p && (p_ == 0. || p_ == 1.))
			return p_; /* better than NaN or infinite loop;
					   FIXME: suboptimal, since -Inf < alpha ! */

					   /* initialize */
		logbeta = lbeta(p, q);

		/* change tail if necessary;  afterwards   0 < a <= 1/2	 */
		if (p_ <= 0.5) {
			a = p_;	pp = p; qq = q; B_swapTail = 0;
		} else { /* change tail, swap  p <-> q :*/
			a = (!lower_tail && !log_p)? alpha : 1 - p_;
			pp = q; qq = p; B_swapTail = 1;
		}

		/* calculate the initial approximation */

		/* y := {fast approximation of} qnorm(1 - a) :*/
		r = sqrt(-2 * log(a));
		y = r - (const1 + const2 * r) / (1. + (const3 + const4 * r) * r);
		if (pp > 1 && qq > 1) {
			r = (y * y - 3.) / 6.;
			s = 1. / (pp + pp - 1.);
			t = 1. / (qq + qq - 1.);
			h = 2. / (s + t);
			w = y * sqrt(h + r) / h - (t - s) * (r + 5. / 6. - 2. / (3. * h));
			xinbta = pp / (pp + qq * exp(w + w));
		} else {
			r = qq + qq;
			t = 1. / (9. * qq);
			t = r * pow(1. - t + y * sqrt(t), 3.0);
			if (t <= 0.)
				xinbta = 1. - exp((util_log1p(-a)+ log(qq) + logbeta) / qq);
			else {
				t = (4. * pp + r - 2.) / t;
				if (t <= 1.)
					xinbta = exp((log(a * pp) + logbeta) / pp);
				else
					xinbta = 1. - 2. / (t + 1.);
			}
		}

		/* solve for x by a modified newton-raphson method, */
		/* using the function pbeta_raw */

		r = 1 - pp;
		t = 1 - qq;
		yprev = 0.;
		adj = 1;
		/* Sometimes the approximation is negative! */
		if (xinbta < lower)
			xinbta = 0.5;
		else if (xinbta > upper)
			xinbta = 0.5;

		/* Desired accuracy should depend on  (a,p)
		* This is from Remark .. on AS 109, adapted.
		* However, it's not clear if this is "optimal" for IEEE double prec.

		* acu = fmax2(acu_min, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));

		* NEW: 'acu' accuracy NOT for squared adjustment, but simple;
		* ---- i.e.,  "new acu" = sqrt(old acu)

		*/
		acu = fmax2(acu_min, pow(10., -13 - 2.5/(pp * pp) - 0.5/(a * a)));
		tx = prev = 0.;	/* keep -Wall happy */

		for (i_pb=0; i_pb < 1000; i_pb++) {
			y = x_betaPraw(xinbta, pp, qq, /*lower_tail = */1, 0);
			if (y != y) return numeric_limits<double>::quiet_NaN();

			y = (y - a) *
				exp(logbeta + r * log(xinbta) + t * util_log1p(-xinbta));
			if (y * yprev <= 0.)
				prev = fmax2(fabs(adj),fpu);
			g = 1;
			for (i_inn=0; i_inn < 1000;i_inn++) {
				adj = g * y;
				if (fabs(adj) < prev) {
					tx = xinbta - adj; /* trial new x */
					if (tx >= 0. && tx <= 1) {
						if (prev <= acu)	goto L_converged;
						if (fabs(y) <= acu) goto L_converged;
						if (tx != 0. && tx != 1)
							break;
					}
				}
				g /= 3;
			}
			if (fabs(tx - xinbta) < 1e-15*xinbta) goto L_converged;
			xinbta = tx;
			yprev = y;
		}
		/*-- NOT converged: Iteration count --*/
		//    ML_ERROR(ME_PRECISION, "qbeta");

	L_converged:
		return B_swapTail ? 1 - xinbta : xinbta;
	}

} // END NAMESPACE onetool
