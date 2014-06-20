/* amoeba.h */
#ifndef _AMOEBA_H_
#define _AMOEBA_H_

#include <vector>
#include <cmath>
#include <myutils.h>

using namespace std;
using namespace myutils;

class AmoebaFunction {
public:
	// Pure virtual function for the minimand: to be defined in the derived class
	virtual double f(const vector<double>& x) = 0;
	// Constructor to set default value of EPS
	AmoebaFunction() {}
};

class Amoeba {
public:
	AmoebaFunction &AmoebaFunc;
	bool coutput;
	
	int n_iterations;			// number of iterations taken to find function_minimum
	double function_minimum;	// value of BFGSFunc.f() at its minimum

	vector< vector<double> > P;	// storage for the vertices of the simplex to be explored during minimization
	vector<double> Y;			// storage for the function evaluations to be tried during minimization

	const double FTOL;			// fractional tolerance to be achieved in function evaluation
	
public:
	Amoeba(AmoebaFunction &AmoebaFunc_in) : AmoebaFunc(AmoebaFunc_in), coutput(false), n_iterations(0), FTOL(1.0e-10) {}
	const vector<double> minimize(const vector<double>& parameters, const double lambda = 1.0) {
		// The value of parameters is taken as the initial point and lambda is used to generate an initial simplex
		const int n = parameters.size();
		P = vector< vector<double> >(n+1,vector<double>(n,0.0));
		Y = vector<double>(n+1,0.0);
		int i,j;
		for(i=0;i<n+1;i++) {
			for(j=0;j<n;j++) {
				P[i][j] = parameters[j];
				if(j==i-1) P[i][j] += lambda;
			}
			// Initial function evaluation
			Y[i] = AmoebaFunc.f(P[i]);
		}
		// Minimize
		amoeba(P,Y,FTOL,n_iterations);
		if(coutput) {
			if(false) cout << "Minimization failed" << endl;
			else {
				for(i=0;i<n+1;i++) {
					cout << "Function is minimized at f(";
					for(j=0;j<n;j++) cout << P[i][j] << " ";
					cout << "\b) = " << Y[i] << endl;
				}
			}
		}
		// The choice of which vertex to return is arbitrary, but should not matter if convergence is reached
		function_minimum = Y[0];
		return P[0];
	}
protected:
	/*
	 Multidimensional minimization of the function funk(x) where x[1..ndim] is a vector in ndim dimensions, 
	 by the downhill simplex method of Nelder and Mead. The matrix p[1..ndim+1] [1..ndim] is input. Its 
	 ndim+1 rows are ndim-dimensional vectors which are the vertices of the starting simplex. Also input is 
	 the vector y[1..ndim+1], whose components must be pre- initialized to the values of funk evaluated at 
	 the ndim+1 vertices (rows) of p; and ftol the fractional convergence tolerance to be achieved in the 
	 function value (n.b.!). On output, p and y will have been reset to ndim+1 new points all within ftol 
	 of a minimum function value, and nfunk gives the number of function evaluations taken.
	 */
	void amoeba(vector< vector<double> > &p, vector<double> &y, const double ftol, int &nfunk) {
		const int NMAX=5000;
		const double TINY=1.0e-10;
		int i,ihi,ilo,inhi,j;
		double rtol,ysave,ytry;
		
		int mpts=p.size();
		int ndim=p[0].size();
		vector<double> psum(ndim);
		nfunk=0;
		get_psum(p,psum);
		for (;;) {
			ilo=0;
			ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
			for (i=0;i<mpts;i++) {
				if (y[i] <= y[ilo]) ilo=i;
				if (y[i] > y[ihi]) {
					inhi=ihi;
					ihi=i;
				} else if (y[i] > y[inhi] && i != ihi) inhi=i;
			}
			rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
			if (rtol < ftol) {
				SWAP(y[0],y[ilo]);
				for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i]);
				break;
			}
			if (nfunk >= NMAX) error("amoeba(): NMAX exceeded");
			nfunk += 2;
			ytry=amotry(p,y,psum,ihi,-1.0);
			if (ytry <= y[ilo])
				ytry=amotry(p,y,psum,ihi,2.0);
			else if (ytry >= y[inhi]) {
				ysave=y[ihi];
				ytry=amotry(p,y,psum,ihi,0.5);
				if (ytry >= ysave) {
					for (i=0;i<mpts;i++) {
						if (i != ilo) {
							for (j=0;j<ndim;j++)
								p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
							y[i]=AmoebaFunc.f(psum);
						}
					}
					nfunk += ndim;
					get_psum(p,psum);
				}
			} else --nfunk;
		}
	}
	/*
	 Extrapolates by a factor fac through the face of the simplex across from the high point, tries it, 
	 and replaces the high point if the new point is better.
	 */
	double amotry(vector< vector<double> > &p, vector<double> &y, vector<double> &psum, const int ihi, const double fac) {
		int j;
		double fac1,fac2,ytry;
		
		int ndim=p[0].size();
		vector<double> ptry(ndim);
		fac1=(1.0-fac)/ndim;
		fac2=fac1-fac;
		for (j=0;j<ndim;j++)
			ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
		ytry=AmoebaFunc.f(ptry);
		if (ytry < y[ihi]) {
			y[ihi]=ytry;
			for (j=0;j<ndim;j++) {
				psum[j] += ptry[j]-p[ihi][j];
				p[ihi][j]=ptry[j];
			}
		}
		return ytry;
	}
	void get_psum(vector< vector<double> > &p, vector<double> &psum) {
		int i,j;
		double sum;
		
		int mpts=p.size();
		int ndim=p[0].size();
		for (j=0;j<ndim;j++) {
			for (sum=0.0,i=0;i<mpts;i++)
				sum += p[i][j];
			psum[j]=sum;
		}
	}	
};

#endif // _AMOEBA_H_