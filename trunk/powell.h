#ifndef _POWELL_MINIMISATION_H_
#define _POWELL_MINIMISATION_H_

#include <vector>
#include <math.h>
#include <limits>
#include <myutils.h>
#include "brent.h"

#pragma warning( disable : 4355 )

using namespace std;
using namespace myutils;

class PowellFunction {
public:
	virtual double f(const vector<double>& x) = 0;
};

class Powell : public BrentFunction {
public:
	PowellFunction &PowFunc;
	Brent brent;

	bool coutput;
	int ITMAX;					// maximum number of iterations
	double TINY;				// a small number
	double TOL;					// tolerance

	int N;						// number of dimensions [= p.size()]
	vector<double> p;			// parameter vector for minimum of PowFunc.f()
	Matrix<double> xi;			// Matrix of vector directions
	double function_minimum;	// value of PowFunc.f() at its minimum
	int n_iterations;			// number of iterations taken to find function_minimum

//	int BrentFunc_i;			// the column in xi that is being minimized one-dimensionally
	vector<double> BrentFunc_xt;// parameters to be fed into one-dimensional minimization
	vector<double> BrentFunc_xi;

	bool fail;

public:
	Powell(PowellFunction &PowFunc_in) : PowFunc(PowFunc_in), ITMAX(200), TINY(1.0e-25), TOL(1.0e-8), coutput(false), brent(*this) {}

	const vector<double>& minimize(const vector<double>& parameters, const double tol) {
		fail = false;
		p = parameters;
		n_iterations = 0;
		N = (int)parameters.size();
		xi = Matrix<double>(N,N,0.0);
		int i;
		for(i=0;i<N;i++) xi[i][i] = 1.;
		powell(tol, n_iterations, function_minimum);
		if(coutput) {
			if(fail) cout << "Minimization failed" << endl;
			else {
				cout << "Function is minimized at f(";
				for(i=0;i<N;i++) cout << p[i] << " ";
				cout << "\b) = " << function_minimum << endl;
			}
		}
		return p;
	}
	double f(const double x) {
		for(int j=0;j<N;j++)
			BrentFunc_xt[j] = p[j] + x * BrentFunc_xi[j];
		return PowFunc.f(BrentFunc_xt);
	}

protected:
	void powell(const double ftol, int &iter, double &fret)
	{
		fail = false;
		int i,j,ibig;
		double del,fp,fptt,t;

		BrentFunc_xt = vector<double>(N);
		BrentFunc_xi = vector<double>(N);

		vector<double> pt = p;
		vector<double> ptt(N);//,xit(N);
		fret = PowFunc.f(p);
		
		for (iter=0;;++iter) {
			fp=fret;
			ibig=0;
			del=0.0;
			for (i=0;i<N;i++) {
				for (j=0;j<N;j++) BrentFunc_xi[j]=xi[j][i]; /*copying is so we can minimize along this direction*/
				fptt=fret;
				fret = linmin();
				if (fptt-fret > del) {
					del=fptt-fret;
					ibig=i+1;
				}
			}
			if (2.0*(fp-fret) <= ftol*(FABS(fp)+FABS(fret))+TINY) {
				return;
			}
			if (iter == ITMAX) {
				//error("Powell: Too many iterations");
				fail = true;
				return;
			}
			for (j=0;j<N;j++) {
				ptt[j]=2.0*p[j]-pt[j];
				BrentFunc_xi[j]=p[j]-pt[j];
				pt[j]=p[j];
			}
			fptt=PowFunc.f(ptt);
			if (fptt < fp) {
				t=2.0*(fp-2.0*fret+fptt)*SQR(fp-fret-del)-del*SQR(fp-fptt);
				if (t < 0.0) {
					fret = linmin();
					for (j=0;j<N;j++) {
						xi[j][ibig-1]=xi[j][N-1];
						xi[j][N-1]=BrentFunc_xi[j];
					}
				}
			}
		}
	}
	inline double linmin() {
		double xmin,fret;
		xmin = brent.minimize(0.0,1.0,TOL);
		fret = brent.function_minimum;
		for(int j=0;j<N;j++) {
			BrentFunc_xi[j] *= xmin;
			p[j] += BrentFunc_xi[j];
		}
		if(brent.fail) fail = true;
		return fret;
	}
	inline double FABS(const double &a) {
		return a < 0.0 ? -a : a;
	}
	inline double SQR(const double a) {
		return a*a;
	}
};

#pragma warning( default : 4355 )

#endif // _POWELL_MINIMISATION_H_
