/*  Copyright 2012 Daniel Wilson.
 *
 *  random.h
 *  Part of the myutils library.
 *
 *  The myutils library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  The myutils library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the myutils library. If not, see <http://www.gnu.org/licenses/>.
 */
/********************************************/
/*	random.h 23rd February 2005				*/
/*	(c) Danny Wilson and Numerical Recipes	*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <cmath>
#include <time.h>
#include <vector>
#include "myutils/vector.h"
#include "myutils/matrix.h"
#include "myutils/lotri_matrix.h"

#include "myutils/myerror.h"

namespace myutils {
class Random {
protected:
	/* protected member variables */
	int seed;
	/* protected member variables used by ran2() */
	int idum;
	int idum2,iy;
	int *iv;
	const int NTAB;
	int protected_ncalls;
	/* protected member variables used by binomial() */
	int nold;
	double pold,pc,plog,pclog,en,oldg;
	/* protected member variables used by poisson() */
	double sq,alxm,g,oldm;
	/* protected member variables used by Z() */
	int iset;
	double gset;
	
protected:
	int autosetseed(void)
	{
		time_t lt;
		lt=time(NULL);

		return (int)lt;
	}
	/* uniform random number generation */
	inline double ran2(void)
	{
		++protected_ncalls;
		const int IM1=2147483563,IM2=2147483399;
		const int IA1=40014,IA2=40692,IQ1=53668,IQ2=52774;
		const int IR1=12211,IR2=3791,IMM1=IM1-1;
		const int NDIV=1+IMM1/NTAB;
		const double EPS=3.0e-16,RNMX=1.0-EPS,AM=1.0/double(IM1);
		int j,k;
		double temp;

		if (idum <= 0) {
			idum=(idum==0 ? 1 : -idum);
			idum2=idum;
			for (j=NTAB+7;j>=0;j--) {
				k=idum/IQ1;
				idum=IA1*(idum-k*IQ1)-k*IR1;
				if (idum < 0) idum += IM1;
				if (j < NTAB) iv[j] = idum;
			}
			iy=iv[0];
		}
		k=idum/IQ1;
		idum=IA1*(idum-k*IQ1)-k*IR1;
		if (idum < 0) idum += IM1;
		k=idum2/IQ2;
		idum2=IA2*(idum2-k*IQ2)-k*IR2;
		if (idum2 < 0) idum2 += IM2;
		j=iy/NDIV;
		iy=iv[j]-idum2;
		iv[j] = idum;
		if (iy < 1) iy += IMM1;
		if ((temp=AM*iy) > RNMX) {
			return RNMX;
		}
		else {
			return temp;
		}
	}
	void rerror(const char* error_text)
	// Standard error handler
	{
		printf("Random Package run-time error...\n");
		printf("%s\n", error_text);
		printf("...now exiting to system...\n");
		exit(13);
	}
	/* 0 < a <= 1. From Devroye (1986) p. 425 */
	double ahrens_dieter74_gamma(const double a)
	{
		double b = (exp(1.)+a)/exp(1.);
		double c = 1./a;
		double U,V,W,X;
		while(true) {
			U = ran2();
			W = ran2();
			V = b * U;
			if(V<=1) {
				X = pow(V,c);
				if(W<=exp(-X)) break;
			}
			else {
				X = -log(c*(b-V));
				if(W<=pow(X,a-1.)) break;
			}
		}
		return X;
	}
	/* a > 1. From Devroye (1986) p. 410 */
	double best78_gamma(const double a)
	{
		double b = a - 1.;
		double c = 3.*a - 0.75;
		double U,V,W,X,Y,Z;
		while(true) {
			U = ran2();
			V = ran2();
			W = U*(1.-U);
			Y = sqrt(c/W)*(U-0.5);
			X = b + Y;
			if(X>=0) {
				Z = 64. * pow(W,3.) * pow(V,2.);
				if(Z <= 1.0 - 2.0*pow(Y,2.)/X) break;
				if(log(Z) <= 2.*(b * log(X/b) - Y)) break;
			}
		}
		return X;
	}
	/* positive integers for ia only. From Numerical Recipes */
	double gamdev(const int ia)
	{
		int j;
		double am,e,s,v1,v2,x,y;

		if (ia < 1) error("Error in routine gamma");
		if (ia < 6) {
			x=1.0;
			for (j=1;j<=ia;j++) x *= ran2();
			x = -log(x);
		} else {
			do {
				do {
					do {
						v1=ran2();
						v2=2.0*ran2()-1.0;
					} while (v1*v1+v2*v2 > 1.0);
					y=v2/v1;
					am=ia-1;
					s=sqrt(2.0*am+1.0);
					x=s*y+am;
				} while (x <= 0.0);
				e=(1.0+y*y)*exp(am*log(x/am)-s*y);
			} while (ran2() > e);
		}
		return x;
	}
	double gammln(const double xx)
	{
		int j;
		double x,y,tmp,ser;
		static const double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
			-0.5395239384953e-5};

		y=x=xx;
		tmp=x+5.5;
		tmp -= (x+0.5)*log(tmp);
		ser=1.000000000190015;
		for (j=0;j<6;j++) ser += cof[j]/++y;
		return -tmp+log(2.5066282746310005*ser/x);
	}

public:
	/* Default constructor */
	Random() : NTAB(32) {
		iv = new int[NTAB];
		setseed(-autosetseed());
		nold = -1;
		pold = -1.0;
		oldm = -1.0;
		iset = 0;
	}
	/* Copy constructor */
	Random(const Random &ran) : NTAB(32) {
		seed = ran.seed;
		iv = new int[NTAB];
		int i;
		for(i=0;i<NTAB;i++) iv[i] = ran.iv[i];
		idum = ran.idum;
		idum2 = ran.idum2;
		iy = ran.iy;
		protected_ncalls = ran.protected_ncalls;
		/* protected member variables used by binomial() */
		nold = ran.nold;
		pold = ran.pold;
		pc = ran.pc;
		plog = ran.plog;
		pclog = ran.pclog;
		en = ran.en;
		oldg = ran.oldg;
		/* protected member variables used by poisson() */
		sq = ran.sq;
		alxm = ran.alxm;
		g = ran.g;
		oldm = ran.oldm;
		/* protected member variables used by Z() */
		iset = ran.iset;
		gset = ran.gset;
	}
	/* Assignment operator */
	Random& operator=(const Random &ran) {
		seed = ran.seed;
		int i;
		for(i=0;i<NTAB;i++) iv[i] = ran.iv[i];
		idum = ran.idum;
		idum2 = ran.idum2;
		iy = ran.iy;
		protected_ncalls = ran.protected_ncalls;
		/* protected member variables used by binomial() */
		nold = ran.nold;
		pold = ran.pold;
		pc = ran.pc;
		plog = ran.plog;
		pclog = ran.pclog;
		en = ran.en;
		oldg = ran.oldg;
		/* protected member variables used by poisson() */
		sq = ran.sq;
		alxm = ran.alxm;
		g = ran.g;
		oldm = ran.oldm;
		/* protected member variables used by Z() */
		iset = ran.iset;
		gset = ran.gset;
		return *this;
	}
	/* Destructor */
	~Random() {
		delete[] iv;
	}
	/* Equality operator */
	bool operator==(const Random &ran) const {
		int i;
		for(i=0;i<NTAB;i++) if(iv[i] != ran.iv[i]) return false;
		if(idum != ran.idum) return false;
		if(idum2 != ran.idum2) return false;
		if(iy != ran.iy) return false;
		if(protected_ncalls != ran.protected_ncalls) return false;
		/* protected member variables used by binomial() */
		if(nold != ran.nold) return false;
		if(pold != ran.pold) return false;
		if(pc != ran.pc) return false;
		if(plog != ran.plog) return false;
		if(pclog != ran.pclog) return false;
		if(en != ran.en) return false;
		if(oldg != ran.oldg) return false;
		/* protected member variables used by poisson() */
		if(sq != ran.sq) return false;
		if(alxm != ran.alxm) return false;
		if(g != ran.g) return false;
		if(oldm != ran.oldm) return false;
		/* protected member variables used by Z() */
		if(iset != ran.iset) return false;
		if(gset != ran.gset) return false;
		return true;		
	}
	/* Inqquality operator */
	bool operator!=(const Random &ran) const {
		return !operator==(ran);
	}
	/* seed_in must be a negative integer */
	Random& setseed(const int seed_in) {
		if(seed_in>0) error("Random must be seeded with a negative integer");
		seed=seed_in;
		idum=seed;
		idum2=123456789;
		iy=0;
		protected_ncalls=0;
		return *this;
	}
	/* seed_in must be a negative integer. set_ncalls is # calls to ran2() */
	Random& setseed(const int seed_in, const int set_ncalls) {
		if(seed_in>0) error("Random must be seeded with a negative integer");
		if(set_ncalls<0) error("ncalls must be non-negative");
		if(seed!=seed_in || protected_ncalls>set_ncalls) setseed(seed_in);
		while(protected_ncalls<set_ncalls) ran2();
		return *this;
	}
	int getseed() {
		return seed;
	}
	Random& setidum(const int idum_in, const int idum2_in, const int iy_in, const std::vector<int> &iv_in) {
		if(iv_in.size()!=NTAB) error("Random::setidum(): iv must have size NTAB");
		seed=1;		/* positive seed indicates it was not properly set */
		idum=idum_in;
		idum2=idum2_in;
		iy=iy_in;
		int i;
		for(i=0;i<NTAB;i++) iv[i] = iv_in[i];
		return *this;
	}
	Random& getidum(int &idum_out, int &idum2_out, int &iy_out, std::vector<int> &iv_out) {
		idum_out = idum;
		idum2_out = idum2;
		iy_out = iy;
		iv_out = std::vector<int>(NTAB);
		int i;
		for(i=0;i<NTAB;i++) iv_out[i] = iv[i];
		return *this;
	}
	Random& printidum() {
		printf("idum = %d\nidum2 = %d\niy = %d\niv = %d",idum,idum2,iy,iv[0]);
		int i;
		for(i=1;i<NTAB;i++) printf(", %d",iv[i]);
		printf("\n");
		return *this;
	}
	int ncalls() {
		return protected_ncalls;
	}
	int bernoulli(const double p)
	{
		double rnumber = ran2();
		if (p<=rnumber) return 0;
		else return 1;
	}
	bool bernoulliTF(const double p)
	{
		double rnumber = ran2();
		if (p<=rnumber) return false;
		else return true;
	}
	double beta(const double a, const double b)
	{
		if(a<=0.0 || b<=0.0) error("Error in beta: a and b parameters must be >0");
		double gam1,gam2;

		if(a == 1.0) gam1 = exponential(1.0);
		else if(a == (double)((int) a)) gam1 = gamdev((int)a);
		else if(a < 1.0) gam1 = ahrens_dieter74_gamma(a);
		else gam1 = best78_gamma(a);

		if(b == 1.0) gam2 = exponential(1.0);
		else if(b == (double)((int) b)) gam2 = gamdev((int)b);
		else if(b < 1.0) gam2 = ahrens_dieter74_gamma(b);
		else gam2 = best78_gamma(b);

		return gam1/(gam1+gam2);
	}
	double binomial(const int n, const double pp)
	{
#ifndef PI
		const double PI=3.141592653589793238;
#endif
		int j;
		// Static members made class members 13/04/09
		//static int nold=(-1);
		double am,em,g,angle,p,bnl,sq,t,y;
		//static double pold=(-1.0),pc,plog,pclog,en,oldg;

		p=(pp <= 0.5 ? pp : 1.0-pp);
		am=n*p;
		if (n < 25) {
			bnl=0.0;
			for (j=0;j<n;j++)
				if (ran2() < p) ++bnl;
		} else if (am < 1.0) {
			g=exp(-am);
			t=1.0;
			for (j=0;j<=n;j++) {
				t *= ran2();
				if (t < g) break;
			}
			bnl=(j <= n ? j : n);
		} else {
			if (n != nold) {
				en=n;
				oldg=gammln(en+1.0);
				nold=n;
			} if (p != pold) {
				pc=1.0-p;
				plog=log(p);
				pclog=log(pc);
				pold=p;
			}
			sq=sqrt(2.0*am*pc);
			do {
				do {
					angle=PI*ran2();
					y=tan(angle);
					em=sq*y+am;
				} while (em < 0.0 || em >= (en+1.0));
				em=floor(em);
				t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
					-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
			} while (ran2() > t);
			bnl=em;
		}
		if (p != pp) bnl=n-bnl;
		return bnl;
	}
	//double *dirichlet(const int *p, const int k)
	//{
	//	double *result;
	//	result=(double *)malloc((unsigned) k*sizeof(double));
	//	if (!result) error("Allocation failure in dirichlet");

	//	double total=0.0;
	//	int i;
	//	for (i=0;i<k;i++)
	//	{
	//		result[i] = gamdev(p[i]);
	//		total += result[i];
	//	}
	//	for (i=0;i<k;i++)
	//		result[i] /= total;

	//	return result;
	//}
	//std::vector<double> dirichlet(const std::vector<int> &p, const int k)
	//{
	//	std::vector<double> result(k);
	//	
	//	double total=0.0;
	//	int i;
	//	for (i=0;i<k;i++)
	//	{
	//		result[i] = gamdev(p[i]);
	//		total += result[i];
	//	}
	//	for (i=0;i<k;i++)
	//		result[i] /= total;

	//	return result;
	//}
	/* random variates returned in double* r */
	void dirichlet(const double *a, const int k, double *r)
	{
		double total=0.0;
		int i;
		for(i=0;i<k;i++) {
			if(a[i] == 1.0) r[i] = exponential(1.0);
			else if(a[i] == (double)((int) a[i])) r[i] = gamdev((int)a[i]);
			else if(a[i] < 1.0) r[i] = ahrens_dieter74_gamma(a[i]);
			else r[i] = best78_gamma(a[i]);
			total += r[i];
		}
		for(i=0;i<k;i++) r[i] /= total;
	}
	/* random variates returned in std::vector<double> r */
	void dirichlet(const std::vector<double> &a, std::vector<double> &r)
	{
		double total=0.0;
		int i;
		int k = (int) a.size();
		if(r.size()!=k) r.resize(k);
		for(i=0;i<k;i++) {
			if(a[i] == 1.0) r[i] = exponential(1.0);
			else if(a[i] == (double)((int) a[i])) r[i] = gamdev((int)a[i]);
			else if(a[i] < 1.0) r[i] = ahrens_dieter74_gamma(a[i]);
			else r[i] = best78_gamma(a[i]);
			total += r[i];
		}
		for(i=0;i<k;i++) r[i] /= total;
	}
	int discrete(const int a, const int b)
	{
		double rnumber = ran2();					// uniform continuous [0,1]
		rnumber *= (b-a+1);							// uniform continuous [0,b-a+1]
		int result = static_cast<int>(rnumber);		// uniform discrete [0,b-a]
		return result + a;							// uniform discrete [a,b]
	}
	double exponential(const double mean)
	{
		double dum;

		do
			dum=ran2();
		while (dum == 0.0);
		return -log(dum)*mean;
	}
	double exponential_ratio() {
		double dum1,dum2;
		do dum1 = ran2(); while(dum1 == 0.0);
		do dum2 = ran2(); while(dum2 == 0.0);
		return log(dum1)/log(dum2);
	}
	/* b is the scale parameter, c the shape parameter. mean = bc, variance = bbc */
	double gamma(const double b, const double c)
	{
		if (b<=0) error("Error in gamma: 1st parameter should be >0");
		if (c<=0) error("Error in gamma: 2nd parameter should be >0");
		if (c == 1.0) return exponential(b);
		int cint = (int) c;
		if (c == (double) cint) return b*gamdev(cint);
		if (c<1.0) return b*ahrens_dieter74_gamma(c);
		return b*best78_gamma(c);
	}
	/* If X ~ geometric(p) then E(X) = (1-p)/p and E(X+1) = 1/p */
	int geometric(const double p)
	{
		return (int)ceil(log(U())/log(1.-p)-1.);
	}
	double inverse(const double a, const double b) {
		if(a<=0.0) error("Lower bound for inverse distribution must be positive");
		if(b<=a) error("Upper bound must be greater than lower bound for inverse distribution");
		return a*pow(b/a,U());
	}
	/* Returns X where Y=log(X) ~ Normal(mu,sigma) */
	double log_normal(const double mu, const double sigma) {
		return exp(normal(mu,sigma));
	}
	/* Returns the minimum of n uniform(0,1) random deviates */
	double minU(const int n)
	{
		return 1.-pow(1.-ran2(),1.0/(double)n);
	}
	int *multinomial(const double* p, const int n, const int k)
	{
		int *result;
		result=(int *)malloc((unsigned) k*sizeof(int));
		if (!result) error("Allocation failure in multinomial");
		int i;
		for (i=0;i<k;i++) result[i]=0;

		double pmax=p[0], pnow;
		for (i=1;i<k;i++)
		{
			pnow=p[i];
			if (pnow>pmax) pmax=pnow;
		}

		int j=n, rnum2; 
		double rnum1,ratio;
		do
		{
			rnum1 = ran2();
			rnum2 = discrete(0,k-1);
			ratio = p[rnum2]/pmax;
			if (rnum1 <= ratio)
			{
				++result[rnum2];
				--j;
			}
		} while (j>0);

		return result;
	}
	int *multinomial(const double* p, const double pmax, const int n, const int k)
	{
		int *result;
		result=(int *)malloc((unsigned) k*sizeof(int));
		if (!result) error("Allocation failure in multinomial");
		for (int i=0;i<k;i++) result[i]=0;

		int j=n, rnum2; 
		double rnum1,ratio;
		do
		{
			rnum1 = ran2();
			rnum2 = discrete(0,k-1);
			ratio = p[rnum2]/pmax;
			if (rnum1 <= ratio)
			{
				++result[rnum2];
				--j;
			}
		} while (j>0);

		return result;
	}
	std::vector<int> multinomial(const std::vector<double> &p, const int n, const int k)
	{
		std::vector<int> result(k);
		int i;
		for (i=0;i<k;i++) result[i]=0;

		double pmax=p[0], pnow;
		for (i=1;i<k;i++)
		{
			pnow=p[i];
			if (pnow>pmax) pmax=pnow;
		}

		int j=n, rnum2; 
		double rnum1,ratio;
		do
		{
			rnum1 = ran2();
			rnum2 = discrete(0,k-1);
			ratio = p[rnum2]/pmax;
			if (rnum1 <= ratio)
			{
				++result[rnum2];
				--j;
			}
		} while (j>0);

		return result;
	}
	std::vector<int> multinomial(const std::vector<double> &p, const double pmax, const int n, const int k)
	{
		std::vector<int> result(k);
		int i;
		for (i=0;i<k;i++) result[i]=0;

		int j=n, rnum2; 
		double rnum1,ratio;
		do
		{
			rnum1 = ran2();
			rnum2 = discrete(0,k-1);
			ratio = p[rnum2]/pmax;
			if (rnum1 <= ratio)
			{
				++result[rnum2];
				--j;
			}
		} while (j>0);

		return result;
	}
	/* p and result have length k. Sum of result equals n */
	void multinomial(const double* p, const int k, int* result, const int n)
	{
		int i;
		for (i=0;i<k;i++) result[i]=0;

		double pmax=p[0], pnow;
		for (i=1;i<k;i++)
		{
			pnow=p[i];
			if (pnow>pmax) pmax=pnow;
		}

		int j=n, rnum2; 
		double rnum1,ratio;
		do
		{
			rnum1 = ran2();
			rnum2 = discrete(0,k-1);
			ratio = p[rnum2]/pmax;
			if (rnum1 <= ratio)
			{
				++result[rnum2];
				--j;
			}
		} while (j>0);
	}
	/* Returns the random variates in the Vector MN */
	void multivariate_normal(Vector<double> &mu, Matrix<double> &Sigma, Vector<double> &MN) {
		Matrix<double> temp;
		Vector<double> z;
		return multivariate_normal(mu,Sigma,MN,temp,z);
	}
	/* Returns the random variates in the Vector MN */
	void multivariate_normal(Vector<double> &mu, Matrix<double> &Sigma, Vector<double> &MN, Matrix<double> &temp, Vector<double> &z, bool *cholesky_fail=0) {
		/*	Cholesky decomposition from Numerical Recipies in C++ */
		/*	Note that eigen decomposition is stabler, and might better pick up
			non-positive definite Sigma. If not picked up, the empirical
			variance-covariance matrix for the simulations will not equal Sigma. */
		int i,j,k;
		double sum;

		int n = Sigma.nrows();
		if(n!=Sigma.ncols()) error("multivariate_normal(): Sigma is not a square matrix");
		if(n!=mu.size()) error("multivariate_normal(): mu and Sigma have incompatible sizes");
		if(cholesky_fail!=0) *cholesky_fail = false;
		temp.resize(n,n);
		z.resize(n);
		MN.resize(n);
		for(i=0;i<n;i++)
			for(j=0;j<=i;j++)
				temp[i][j] = temp[j][i] = (double)Sigma[i][j];

		for(i=0;i<n;i++) {
			for(j=i;j<n;j++) {
				for(sum=temp[i][j],k=i-1;k>=0;k--) sum -= temp[i][k]*temp[j][k];
				if(i==j) {
					if(sum <= 0.0) {/* Sigma, with rounding errors, is not positive definite */
						if(cholesky_fail!=0) {
							*cholesky_fail = true;
							return;
						}
						printf("\nSigma = \n");
						int ii,jj;
						for(ii=0;ii<n;ii++) {
							for(jj=0;jj<n;jj++) printf("%.3g\t",Sigma[ii][jj]);
							printf("\n");
						}
						error("multivariate_normal(): Cholesky decomposition failed, not positive definite");
					}
					/* Temporarily use z to store the diagonal */
					z[i] = sqrt(sum);
				}
				else temp[j][i] = sum/z[i];
			}
		}
		for(i=0;i<n;i++) temp[i][i] = z[i];
		/*	Simulate MultiNormal(mu, Sigma), where Sigma is the variance-covariance matrix.
			Compute the Cholesky decomposition Sigma = L . L', where ' denotes the transpose.
			Generate a vector of i.i.d. standard normal variates Z. Then
					M = L' . Z + mu
			has the desired distribution.*/
		for(i=0;i<n;i++) {
			MN[i] = mu[i];
			z[i] = Z();
		}
		for(i=0;i<n;i++) {
			for(k=0;k<=i;k++) MN[i] += temp[i][k] * z[k];
		}

	}
	/* Returns the random variates in the Vector MN */
	void multivariate_normal(Vector<double> &mu, LowerTriangularMatrix<double> &Cholesky, Vector<double> &MN, Vector<double> &z) {
		/*	Cholesky decomposition from Numerical Recipies in C++ */
		/*	Note that eigen decomposition is stabler, and might better pick up
			non-positive definite Sigma. If not picked up, the empirical
			variance-covariance matrix for the simulations will not equal Sigma. */
		int i,k;
		int n = Cholesky.n();
		if(n!=mu.size()) error("multivariate_normal(): mu and Sigma have incompatible sizes");
		z.resize(n);
		MN.resize(n);

		/*	Simulate MultiNormal(mu, Sigma), where Sigma is the variance-covariance matrix.
			Compute the Cholesky decomposition Sigma = L . L', where ' denotes the transpose.
			Generate a vector of i.i.d. standard normal variates Z. Then
					M = L' . Z + mu
			has the desired distribution.*/
		for(i=0;i<n;i++) {
			MN[i] = mu[i];
			z[i] = Z();
		}
		for(i=0;i<n;i++) {
			for(k=0;k<=i;k++) MN[i] += Cholesky[i][k] * z[k];
		}

	}
	// NB sigma is the standard deviation
	double normal(const double mu, const double sigma)
	{
		double X = Z();		// X ~ N(0,1)
		X *= sigma;				// X ~ N(0,sigma)
		X += mu;				// X ~ N(mu,sigma)
		return X;
	}
	double poisson(const double xm)
	{
#ifndef PI
		const double PI=3.141592653589793238;
#endif
		// Static members made class members 13/04/09
		//static double sq,alxm,g,oldm=(-1.0);
		double em,t,y;

		if (xm < 12.0) {
			if (xm != oldm) {
				oldm=xm;
				g=exp(-xm);
			}
			em = -1;
			t=1.0;
			do {
				++em;
				t *= ran2();
			} while (t > g);
		} else {
			if (xm != oldm) {
				oldm=xm;
				sq=sqrt(2.0*xm);
				alxm=log(xm);
				g=xm*alxm-gammln(xm+1.0);
			}
			do {
				do {
					y=tan(PI*ran2());
					em=sq*y+xm;
				} while (em < 0.0);
					em=floor(em);
					t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
			} while (ran2() > t);
		}
		return em;
	}
	/*b=mean of full distribution t=cutoff point*/
	double trunc_exponential(const double b, const double t)
	{
		return -b*log(1.0-(1.0-exp(-t/b))*ran2());
	}
	/* truncated geometric with range 1..t. Mean of non-truncated distn would be 1/p. */
	int trunc_geometric(const double p, const int t)
	{
		const double a = pow(1.-p,(double)t);
		return (int)ceil(log(a-(a-1.)*ran2())/log(1.-p));
	}
	inline double U(void){return ran2();}
	double uniform(const double a, const double b)
	{
		double rnumber = ran2();			// continuous uniform [0,1]
		rnumber *= (b-a);					// continuous uniform [0,b-a]
		rnumber += a;						// continuous uniform [a,b]
		return rnumber;
	}
	double Z(void)
	{
		// Static members made class members 13/04/09
		//static int iset=0;
		//static double gset;
		double fac,rsq,v1,v2;

		if (idum < 0) iset=0;
		if (iset == 0) {
			do {
				v1=2.0*ran2()-1.0;
				v2=2.0*ran2()-1.0;
				rsq=v1*v1+v2*v2;
			} while (rsq >= 1.0 || rsq == 0.0);
			fac=sqrt(-2.0*log(rsq)/rsq);
			gset=v1*fac;
			iset=1;
			return v2*fac;
		} else {
			iset=0;
			return gset;
		}
	}
};
};

#endif
