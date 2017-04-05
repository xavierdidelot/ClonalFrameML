/*  Copyright 2013 Daniel Wilson.
 *
 *  brent.h
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
 *
 *  Parts of this code are based on code in Numerical Recipes in C++
 *  WH Press, SA Teukolsky, WT Vetterling, BP Flannery (2002).
 *
 */
#ifndef _BRENT_MINIMISATION_
#define _BRENT_MINIMISATION_

#include <limits>
#include "myutils/myerror.h"

using namespace std;

/*	Class Brent performs parabolic interpolation and Brent's method on a one-
	dimensional member function, BrentFunc.f(x). BrentFunc must be an instance of a class
	derived from the abstract class BrentFunction. Its member function f(x) takes only a
	single parameter, but using a derived class allows for it to be controlled by other
	member variables and/or call other member functions, enabling a neater alternative to
	using function pointers and global variables.

	See Numerical Recipes in C++ [Press et al 2002] for details of the algorithm.
*/

class BrentFunction {
public:
	virtual double f(const double x) = 0;
};

/*	An example derived class might look like MyFunction below. By passing an instance of
	MyFunction to an instance of Brent in its constructor, the function MyFunction::f(x)
	can be minimized with respect to x, whilst having an auxilliary variable y, which is not
	minimized.

class MyFunction : public BrentFunction {
	double y;

public:
	MyFunction(const double y_in) : y(y_in) {}
	double f(const double x) {
		return (x+y)*(x+y);
	}
};
*/

class Brent {
public:
	BrentFunction & BrentFunc;

	bool coutput;
	double evala_BrentFunc, evalb_BrentFunc, evalc_BrentFunc;
	double pointa,pointb,pointc;
	double GLIMIT, TINY, tolerance;
	int ITMAX;
	double ZEPS,EPS;
	double function_minimum;
	bool bracketed;

	bool fail;

public:
	Brent(BrentFunction &BrentFunc_in) : BrentFunc(BrentFunc_in), GLIMIT(100.), TINY(1.e-20), ITMAX(100), coutput(false), EPS(3.0e-8) {}
	double minimize(const double pointa_in, const double pointb_in, const double tol) {
		fail = false;
		ZEPS=numeric_limits<double>::epsilon()*1.0e-3;
		pointa = pointa_in;
		pointb = pointb_in;
		pointc = 0.0;
		tolerance = tol;
		mnbrak(pointa, pointb, pointc, evala_BrentFunc, evalb_BrentFunc, evalc_BrentFunc);
		if(coutput) {
			cout << "Function is bracketed by:" << endl;
			cout << "f(" << pointa << ") = " << evala_BrentFunc << endl;
			cout << "f(" << pointb << ") = " << evalb_BrentFunc << endl;
			cout << "f(" << pointc << ") = " << evalc_BrentFunc << endl;
		}
		double result = 0.0;
		function_minimum = brent(pointa, pointb, pointc, result);
		if(coutput)
			cout << "Function is minimized at f(" << result << ") = " << function_minimum << endl;
		return result;
	};

	double rootfind(double x1, double x2, double tol) {
	  //Using Brent�s method, find the root of a function func known to lie between x1 and x2. The
	  //root, returned as zbrent, will be refined until its accuracy is tol.
	  int iter;
	  double a=x1,b=x2,c=x2,d,e,min1,min2;
	  double fa=BrentFunc.f(a),fb=BrentFunc.f(b),fc,p,q,r,s,tol1,xm;
	  bracketed = true;
	  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
	    if(coutput)
	      cout << "f(" << x1 << ") = " << fa << "\tf(" << x2 << ") = " << fb << endl;
	    //myutils::warning("Root must be bracketed in rootfind");
	    bracketed = false;
	    return 0.0;
	  }
	  fc=fb;
	  for (iter=1;iter<=ITMAX;iter++) {
	    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
	      c=a; //Rename a, b, c and adjust bounding interval d.
	      fc=fa;
	      e=d=b-a;
	    }
	    if (fabs(fc) < fabs(fb)) {
	      a=b;
	      b=c;
	      c=a;
	      fa=fb;
	      fb=fc;
	      fc=fa;
	    }
	    tol1=2.0*EPS*fabs(b)+0.5*tol; //Convergence check.
	    xm=0.5*(c-b);
	    if (fabs(xm) <= tol1 || fb == 0.0) return b;
	    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
	      s=fb/fa; //Attempt inverse quadratic interpolation.
	      if (a == c) {
		p=2.0*xm*s;
		q=1.0-s;
	      } else {
		q=fa/fc;
		r=fb/fc;
		p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
		q=(q-1.0)*(r-1.0)*(s-1.0);
	      }
	      if (p > 0.0) q = -q; //Check whether in bounds.
	      p=fabs(p);
	      min1=3.0*xm*q-fabs(tol1*q);
	      min2=fabs(e*q);
	      if (2.0*p < (min1 < min2 ? min1 : min2)) {
		e=d; //Accept interpolation.
		d=p/q;
	      } else {
		d=xm; //Interpolation failed, use bisection.
		e=d;
	      }
	    } else { //Bounds decreasing too slowly, use bisection.
	      d=xm;
	      e=d;
	    }
	    a=b; //Move last best guess to a.
	    fa=fb;
	    if (fabs(d) > tol1) //Evaluate new trial root.
	      b += d;
	    else
	      b += SIGN(tol1,xm);
	    fb=BrentFunc.f(b);
	  }
	  myutils::warning("Maximum number of iterations exceeded in zbrent");
	  return 0.0; //Never get here.
	}
protected:
	/*	The hard work is done by algorithms modified from
		Numerical Recipes in C++ [Press et al 2002]		*/
	inline void shft3(double &a, double &b, double &c, const double d) {
		a=b;
		b=c;
		c=d;
	}
	inline void shft2(double &a, double &b, const double c) {
		a=b;
		b=c;
	}
	void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc) {
		const double GOLD=1.618034;
		double ulim,u,r,q,fu;

		fa = BrentFunc.f(ax);
		fb = BrentFunc.f(bx);
		if (fb > fa) {
			SWAP(ax,bx);
			SWAP(fb,fa);
		}
		cx=bx+GOLD*(bx-ax);
		fc=BrentFunc.f(cx);
		while (fb > fc) {
			r=(bx-ax)*(fb-fc);
			q=(bx-cx)*(fb-fa);
			u=bx-((bx-cx)*q-(bx-ax)*r)/
				(2.0*SIGN(MAX(FABS(q-r),TINY),q-r));
			ulim=bx+GLIMIT*(cx-bx);
			if ((bx-u)*(u-cx) > 0.0) {
				fu=BrentFunc.f(u);
				if (fu < fc) {
					ax=bx;
					bx=u;
					fa=fb;
					fb=fu;
					return;
				} else if (fu > fb) {
					cx=u;
					fc=fu;
					return;
				}
				u=cx+GOLD*(cx-bx);
				fu=BrentFunc.f(u);
			} else if ((cx-u)*(u-ulim) > 0.0) {
				fu=BrentFunc.f(u);
				if (fu < fc) {
					shft3(bx,cx,u,cx+GOLD*(cx-bx));
					shft3(fb,fc,fu,BrentFunc.f(u));
				}
			} else if ((u-ulim)*(ulim-cx) >= 0.0) {
				u=ulim;
				fu=BrentFunc.f(u);
			} else {
				u=cx+GOLD*(cx-bx);
				fu=BrentFunc.f(u);
			}
			shft3(ax,bx,cx,u);
			shft3(fa,fb,fc,fu);
		}
	}
	inline void SWAP(double &a, double &b) {
		double dum=a;a=b;b=dum;
	}
	inline double SIGN(const double &a, const double &b) {
		return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
	}
	inline double MAX(const double &a, const double &b) {
		return b > a ? (b) : (a);
	}
	inline double FABS(const double &a) {
		return a < 0.0 ? -a : a;
	}
	double brent(const double ax, const double bx, const double cx, double &xmin)
	{
		const double CGOLD=0.3819660;
		int iter;
		double a,b,d=0.0,etemp,fu,fv,fw,fx;
		double p,q,r,tol1,tol2,u,v,w,x,xm;
		double e=0.0;

		a=(ax < cx ? ax : cx);
		b=(ax > cx ? ax : cx);
		x=w=v=bx;
		fw=fv=fx=BrentFunc.f(x);
		for (iter=0;iter<ITMAX;iter++) {
			xm=0.5*(a+b);
			tol2=2.0*(tol1=tolerance*FABS(x)+ZEPS);
			if (FABS(x-xm) <= (tol2-0.5*(b-a))) {
				xmin=x;
				return fx;
			}
			if (FABS(e) > tol1) {
				r=(x-w)*(fx-fv);
				q=(x-v)*(fx-fw);
				p=(x-v)*q-(x-w)*r;
				q=2.0*(q-r);
				if (q > 0.0) p = -p;
				q=FABS(q);
				etemp=e;
				e=d;
				if (FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
					d=CGOLD*(e=(x >= xm ? a-x : b-x));
				else {
					d=p/q;
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				}
			} else {
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			}
			u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
			fu=BrentFunc.f(u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				shft3(v,w,x,u);
				shft3(fv,fw,fx,fu);
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					v=w;
					w=u;
					fv=fw;
					fw=fu;
				} else if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
			}
		}
		//myutils::error("Brent: Too many iterations");
		fail = true;
		xmin=x;
		return fx;
	}


};

class ConstrainedBrent {
public:
	BrentFunction & BrentFunc;

	bool coutput;
	double evala_BrentFunc, evalb_BrentFunc, evalc_BrentFunc;
	double pointa,pointb,pointc;
	double GLIMIT, TINY, tolerance;
	int ITMAX;
	double ZEPS;
	double function_minimum;
	double min_x,max_x;

public:
	ConstrainedBrent(BrentFunction &BrentFunc_in) : BrentFunc(BrentFunc_in), GLIMIT(100.), TINY(1.e-20), ITMAX(100), coutput(false) {}
	double minimize(const double pointa_in, const double pointb_in, const double tol, const double min_x_in, const double max_x_in) {
		min_x = min_x_in;
		max_x = max_x_in;
		ZEPS=numeric_limits<double>::epsilon()*1.0e-3;
		pointa = pointa_in;
		pointb = pointb_in;
		pointc = min_x;

		if(pointa<min_x || pointa>max_x) error("ConstrainedBrent::minimize(): point a falls outside range");
		if(pointb<min_x || pointb>max_x) error("ConstrainedBrent::minimize(): point b falls outside range");

		tolerance = tol;
		mnbrak(pointa, pointb, pointc, evala_BrentFunc, evalb_BrentFunc, evalc_BrentFunc);
		if(coutput) {
			cout << "Function is bracketed by:" << endl;
			cout << "f(" << pointa << ") = " << evala_BrentFunc << endl;
			cout << "f(" << pointb << ") = " << evalb_BrentFunc << endl;
			cout << "f(" << pointc << ") = " << evalc_BrentFunc << endl;
		}
		double result = 0.0;
		function_minimum = brent(pointa, pointb, pointc, result);
		if(coutput)
			cout << "Function is minimized at f(" << result << ") = " << function_minimum << endl;
		return result;
	};
protected:
	/*	The hard work is done by algorithms modified from
		Numerical Recipes in C++ [Press et al 2002]		*/
	inline void shft3(double &a, double &b, double &c, const double d) {
		a=b;
		b=c;
		c=d;
	}
	inline void shft2(double &a, double &b, const double c) {
		a=b;
		b=c;
	}
	void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc) {
		const double GOLD=1.618034;
		double ulim,u,r,q,fu;

		fa = BrentFunc.f(ax);
		fb = BrentFunc.f(bx);
		if (fb > fa) {
			SWAP(ax,bx);
			SWAP(fb,fa);
		}
		cx=bx+GOLD*(bx-ax);
		if(cx<min_x) cx = min_x;
		if(cx>max_x) cx = max_x;

		fc=BrentFunc.f(cx);
		while (fb > fc) {
			r=(bx-ax)*(fb-fc);
			q=(bx-cx)*(fb-fa);
			u=bx-((bx-cx)*q-(bx-ax)*r)/
				(2.0*SIGN(MAX(FABS(q-r),TINY),q-r));
			if(u<min_x) u = min_x;
			if(u>max_x) u = max_x;
			ulim=bx+GLIMIT*(cx-bx);
			if ((bx-u)*(u-cx) > 0.0) {
				fu=BrentFunc.f(u);
				if (fu < fc) {
					ax=bx;
					bx=u;
					fa=fb;
					fb=fu;
					return;
				} else if (fu > fb) {
					cx=u;
					fc=fu;
					return;
				}
				u=cx+GOLD*(cx-bx);
				if(u<min_x) u = min_x;
				if(u>max_x) u = max_x;
				fu=BrentFunc.f(u);
			} else if ((cx-u)*(u-ulim) > 0.0) {
				fu=BrentFunc.f(u);
				if (fu < fc) {
					shft3(bx,cx,u,cx+GOLD*(cx-bx));
					if(u<min_x) u = min_x;
					if(u>max_x) u = max_x;
					shft3(fb,fc,fu,BrentFunc.f(u));
				}
			} else if ((u-ulim)*(ulim-cx) >= 0.0) {
				u=ulim;
				if(u<min_x) u = min_x;
				if(u>max_x) u = max_x;
				fu=BrentFunc.f(u);
			} else {
				u=cx+GOLD*(cx-bx);
				if(u<min_x) u = min_x;
				if(u>max_x) u = max_x;
				fu=BrentFunc.f(u);
			}
			shft3(ax,bx,cx,u);
			shft3(fa,fb,fc,fu);
		}
	}
	inline void SWAP(double &a, double &b) {
		double dum=a;a=b;b=dum;
	}
	inline double SIGN(const double &a, const double &b) {
		return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
	}
	inline double MAX(const double &a, const double &b) {
		return b > a ? (b) : (a);
	}
	inline double FABS(const double &a) {
		return a < 0.0 ? -a : a;
	}
	double brent(const double ax, const double bx, const double cx, double &xmin)
	{
		const double CGOLD=0.3819660;
		int iter;
		double a,b,d=0.0,etemp,fu,fv,fw,fx;
		double p,q,r,tol1,tol2,u,v,w,x,xm;
		double e=0.0;

		a=(ax < cx ? ax : cx);
		b=(ax > cx ? ax : cx);
		x=w=v=bx;

		fw=fv=fx=BrentFunc.f(x);
		for (iter=0;iter<ITMAX;iter++) {
			xm=0.5*(a+b);
			tol2=2.0*(tol1=tolerance*FABS(x)+ZEPS);
			if (FABS(x-xm) <= (tol2-0.5*(b-a))) {
				xmin=x;
				return fx;
			}
			if (FABS(e) > tol1) {
				r=(x-w)*(fx-fv);
				q=(x-v)*(fx-fw);
				p=(x-v)*q-(x-w)*r;
				q=2.0*(q-r);
				if (q > 0.0) p = -p;
				q=FABS(q);
				etemp=e;
				e=d;
				if (FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
					d=CGOLD*(e=(x >= xm ? a-x : b-x));
				else {
					d=p/q;
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				}
			} else {
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			}
			u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
			fu=BrentFunc.f(u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				shft3(v,w,x,u);
				shft3(fv,fw,fx,fu);
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					v=w;
					w=u;
					fv=fw;
					fw=fu;
				} else if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
			}
		}
		myutils::error("Brent: Too many iterations");
		xmin=x;
		return fx;
	}


};

#endif // _BRENT_MINIMISATION_
