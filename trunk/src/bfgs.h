/*
 *
 *  bfgs.h
 *
 */
#ifndef _BFGS_H_
#define _BFGS_H_

#include <vector>
#include <math.h>
#include <limits>
#include <myutils.h>

using namespace std;
using namespace myutils;

class BFGSFunction {
public:
	// Pure virtual function for the minimand: to be defined in the derived class
	virtual double f(const vector<double>& x) = 0;

	// Concrete function that ensures the last evaluation of f is stored in last_f
	double record_f(const vector<double>& x) {
		last_f = f(x);
		return last_f;
	}

	// Virtual function for the gradient with default implementation. The gradient is output in g
	virtual void df(const vector<double> &x, vector<double> &g) {
		const int n = x.size();
		xh = x;
		double fold = last_f;
		int j;
		for(j=0;j<n;j++) {
			double temp = x[j];
			double h = EPS*fabs(temp);
			if(h==0.0) h = EPS;
			xh[j] = temp + h;										// Trick to reduce finite-precision error
			h = xh[j] - temp;
			double fh = record_f(xh);
			xh[j] = temp;
			g[j] = (fh-fold)/h;
		}
	}
	
	// Approximate square root of machine precision
	double EPS;
	// Last evaluation of the minimand
	double last_f;
	// Storage for the perturbed parameter used in calculating the gradient
	vector<double> xh;
	// Constructor to set default value of EPS
	BFGSFunction() : EPS(1.0e-8) {}
};

class BFGS /*: public BrentFunction*/ {
public:
	BFGSFunction &BFGSFunc;
	bool coutput;
	
	double gtol;				// tolerance
	int n_iterations;			// number of iterations taken to find function_minimum
	double function_minimum;	// value of BFGSFunc.f() at its minimum
	double STPMX;				// limits the maximum step size to avoid bad areas
	
//	vector<double> dg,g,hdg,pnew,xi;
								// storage for temporary vectors: could save some time on memory allocation
	vector<double> p;			// storage for the parameter values to be manipulated during minimization

	Matrix<double> hessin;		// Approximate inverse Hessian matrix
	
	bool fail;
	
public:
	BFGS(BFGSFunction &BFGSFunc_in) : BFGSFunc(BFGSFunc_in), coutput(false), n_iterations(0), gtol(1.0e-8), STPMX(100.0) {}
	
	const vector<double>& minimize(const vector<double>& parameters, const double tol) {
		fail = false;
		p = parameters;

		dfpmin(p, gtol, n_iterations, function_minimum);
		if(coutput) {
			if(fail) cout << "Minimization failed" << endl;
			else {
				cout << "Function is minimized at f(";
				int i;
				for(i=0;i<p.size();i++) cout << p[i] << " ";
				cout << "\b) = " << function_minimum << endl;
			}
		}
		return p;
	}
protected:
	/*
	 Given a starting point p[0..n-1] the Broyden-Fletcher-Goldfarb-Shanno variant of Davidon-
	 Fletcher-Powell minimization is performed on a function whose value and gradient are provided
	 by the member variable BFGSFunction. The convergence requirement on zeroing the gradient
	 is input as gtol. Returned quantities are p[0..n-1] (the location of the minimum), iter (the
	 number of iterations that were performed), and fret (the minimum value of the function). The
	 routine lnsrch is called to perform approximate line minimizations.
	*/
	void dfpmin(vector<double> &p, const double gtol, int &iter, double &fret) {
		const int ITMAX = 200;
		const double EPS = numeric_limits<double>::epsilon();
		const double TOLX = 4.0*EPS;
		// Here ITMAX is the maximum allowed number of iterations, EPS is the machine precision,
		// TOLX is the convergence criterion on x values, and STPMX is the scaled maximum step length
		// allowed in line searches.
		bool check;
		double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
		int n = p.size();
		vector<double> dg(n),g(n),hdg(n),pnew(n),xi(n);
		hessin = Matrix<double>(n,n);
		fp = BFGSFunc.record_f(p);									// Calculate starting function value and gradient
		BFGSFunc.df(p,g);
		int i,j;
		for(i=0;i<n;i++) {											// Initialize the inverse Hessian to the unit matrix
			for(j=0;j<n;j++) {
				hessin[i][j] = 0.0;
			}
			hessin[i][i] = 1.0;
			xi[i] = -g[i];											// Initial line direction
			sum += p[i]*p[i];
		}
		stpmax = STPMX*MAX(sqrt(sum),(double)n);
		int its;
		for(its=0;its<ITMAX;its++) {								// Main loop over the iterations
			iter = its;
			lnsrch(p,fp,g,xi,pnew,fret,stpmax,check);
			// The new function evaluation occurs in lnsrch, save the function value in fp for the
			// next line search. It is usually safe to ignore the value of check
			fp = fret;
			for(i=0;i<n;i++) {
				xi[i] = pnew[i]-p[i];								// Update the line direction
				p[i] = pnew[i];										// and the current point
			}
			test = 0.0;												// Test for convergence on delta_x
			for(i=0;i<n;i++) {
				temp = fabs(xi[i])/MAX(fabs(p[i]),1.0);
				if(temp>test) test = temp;
			}
			if(test<TOLX) {
				return;
			}
			for(i=0;i<n;i++) {
				dg[i] = g[i];										// Save the old gradient
			}
			BFGSFunc.df(p,g);										// and get the new gradient
			test = 0.0;												// Test for convergence on zero gradient
			den = MAX(fabs(fret),1.0);
			for(i=0;i<n;i++) {
				temp = fabs(g[i])*MAX(fabs(p[i]),1.0)/den;
				if(temp>test) test = temp;
			}
			if(test<gtol) {
				return;
			}
			for(i=0;i<n;i++) {										// Compute difference of gradients
				dg[i] = g[i]-dg[i];
			}
			for(i=0;i<n;i++) {										// and difference times current matrix
				hdg[i] = 0.0;
				for(j=0;j<n;j++) {
					hdg[i] += hessin[i][j]*dg[j];
				}
			}
			fac = fae = sumdg = sumxi = 0.0;						// Calculate dot products for the denominators
			for(i=0;i<n;i++) {
				fac += dg[i]*xi[i];
				fae += dg[i]*hdg[i];
				sumdg += dg[i]*dg[i];
				sumxi += xi[i]*xi[i];
			}
			if(fac>sqrt(EPS*sumdg*sumxi)) {							// Skip update if fac not sufficiently positive
				fac = 1.0/fac;
				fad = 1.0/fae;
				// The vector that makes BFGS different from DFP:
				for(i=0;i<n;i++) dg[i] = fac*xi[i]-fad*hdg[i];
				for(i=0;i<n;i++) {
					for(j=i;j<n;j++) {
						hessin[i][j] += fac*xi[i]*xi[j]
						-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
						hessin[j][i] = hessin[i][j];
					}
				}
			}
			for(i=0;i<n;i++) {										// Now calculate the next direction to go
				xi[i] = 0.0;
				for(j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
			}
		}
		error("dfpmin(): too many iterations");			
	}
	/*
	 Given an n-dimensional point xold[0..n-1] the value of the function and gradient there, fold
	 and g[0..n-1], and a direction p[0..n-1], finds a new point x[0..n-1] along the direction
	 p from xold where the function or functor func has decreased "sufficiently". The new function
	 value is returned in f. stpmax is an input quantity that limits the length of the steps so that you
	 do not try to evaluate the function in regions where it is undefined or subject to overflow. p is
	 usually the Newton direction. The output quantity check is false on a normal exit. It is true
	 when x is too close to xold. In a minimization algorithm, this usually signals convergence and
	 can be ignored. However, in a zero-finding algorithm the calling program should check whether
	 the convergence is spurious.
	*/
	void lnsrch(const vector<double> &xold, const double fold, const vector<double> &g, vector<double> &p, vector<double> &x, double &f, const double stpmax, bool &check) {
		const double ALF = 1.0e-4, TOLX = numeric_limits<double>::epsilon();
		// ALF ensures sufficient decrease in function value, TOLX is the convergence criterion on delta_x
		
		double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
		double rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
		int i,n=xold.size();
		check=false;
		for(i=0;i<n;i++) sum += p[i]*p[i];
		sum = sqrt(sum);
		if(sum>stpmax) {
			for(i=0;i<n;i++) {
				p[i] *= stpmax/sum;									// Scale if attempted step is too big
			}
		}
		for(i=0;i<n;i++) {
			slope += g[i]*p[i];
		}
		if(slope>=0.0) error("lnsrch: Roundoff problem");
		test = 0.0;													// Compute lambda_min
		for(i=0;i<n;i++) {
			temp = fabs(p[i])/MAX(fabs(xold[i]),1.0);
			if(temp>test) test = temp;
		}
		alamin = TOLX/test;
		alam = 1.0;													// Always try full Newton step first
		for(;;) {													// Start of iteration loop
			for(i=0;i<n;i++) x[i] = xold[i] + alam*p[i];
			f = BFGSFunc.record_f(x);
			if(alam<alamin) {										// Convergence on delta_x. For zero finding
				for(i=0;i<n;i++) x[i] = xold[i];					// the calling program should verify convergence
				check = true;
				return;
			} else if(f <= fold+ALF*alam*slope) {					// Sufficient function decrease
				return;
			} else {												// Backtrack
				if(alam==1.0) {
					tmplam = -slope/(2.0*(f-fold-slope));			// First time
				} else {											// Subsequent backtracks
					rhs1 = f-fold-alam*slope;
					rhs2 = f2-fold-alam2*slope;
					a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
					b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
					if(a==0.0) {
						tmplam = -slope/(2.0*b);
					} else {
						disc = b*b-3.0*a*slope;
						if(disc<0.0) {
							tmplam = 0.5*alam;
						} else if(b<=0.0) {
							tmplam = (-b+sqrt(disc))/(3.0*a);
						} else {
							tmplam = -slope/(b+sqrt(disc));
						}
					}
					if(tmplam>0.5*alam) {
						tmplam = 0.5*alam;							// lambda <= 0.5 lambda_1
					}
				}
			}
			alam2 = alam;
			f2 = f;
			alam = MAX(tmplam,0.1*alam);							// lambda >= 0.1 lambda_1
		}
	}
};

#endif // _BFGS_H_
