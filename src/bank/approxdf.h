NOTREACHED
/*  Copyright 2012 Daniel Wilson.
 *
 *  approxdf.h
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
#ifndef _APPROXDF_H_
#define _APPROXDF_H_

#include "myutils/vector.h"
#include "myutils/myerror.h"
#include <math.h>

namespace myutils {

class approxdf {
public:
	int n;
	Vector<double> CDF,G,EV,PR;
public:
	approxdf() {
		n = 0;
		CDF = G = EV = PR = Vector<double>(0);
	}
	approxdf(Vector<double> &EV_in, Vector<double> &PR_in) {
		initialize(EV_in,PR_in);
	}
	void initialize(Vector<double> &EV_in, Vector<double> &PR_in) {
		n = EV_in.size();
		if(PR_in.size()!=n) error("approxdf(): EV and PR must have same length");
		EV = Vector<double>(n);
		PR = Vector<double>(n);
		int i;
		for(i=0;i<n;i++) {
			EV[i] = EV_in[i];
			PR[i] = PR_in[i];
			if(i>0 && EV[i]<EV[i-1]) error("approxdf(): EV must be increasing");
		}
		//# EV is a list (n) of evaluation points
		//# PR is a list (n) of (unnormalized) density estimates
		//# PR <- predict(fit[[fitname]],EV)
		//# PDF is a list (n-1) of (unnormalized) density estimates in
		//# the intervals defined by EV, using a piecewise linear
		//# approximation to the p.d.f. Note that this implies that the
		//# approximation to the c.d.f. is piecewise quadratic.
		Vector<double> PDF(n-1);
		for(i=1;i<n;i++) PDF[i-1] = .5*(EV[i]-EV[i-1])*(PR[i]+PR[i-1]);
		//# CDF is a list (n) of the partial sum of PDF
		CDF = Vector<double>(n); CDF[0] = 0;
		for(i=1;i<n;i++) CDF[i] = CDF[i-1]+PDF[i-1];
		//# PR, PDF and CDF are then normalized
		double t = CDF[n-1];
		for(i=0;i<n;i++) PR[i]/=t;
		for(i=0;i<n;i++) CDF[i]/=t;
		for(i=0;i<n-1;i++) PDF[i]/=t;
		//# G is list (n-1) of estimates of the gradient between each
		//# pair of evaluation points defined by EV
		G = Vector<double>(n-1);
		for(i=1;i<n;i++) G[i-1] = (PR[i]-PR[i-1])/(EV[i]-EV[i-1]);
	}
	double cdf(const double x) {
		int wh;
		for(wh=0;wh<n;wh++) {
			if(EV[wh]>x) {
				--wh;
				break;
			}
		}
		if(wh==-1 || wh==n) error("cdf(): x lies outside original range");
		//# A piecewise quadratic approximation to the c.d.f.
		return CDF[wh]+(x-EV[wh])*(PR[wh]+G[wh]/2*(x-EV[wh]));
	}
	double icdf(const double U) {
		int wh;
		for(wh=0;wh<n;wh++) {
			if(CDF[wh]>U) {
				--wh;
				break;
			}
		}
		if(wh==-1) error("icdf(): U is less than 0");
		if(wh==n) error("icdf(): U is greater than 1");
		//# A piecewise inverse-quadratic approximation to the i.c.d.f.
		return ((G[wh]*EV[wh]-PR[wh]+sqrt(PR[wh]*PR[wh]+2*G[wh]*(U-CDF[wh])))/G[wh]);
	}
	double pdf(const double x) {
		int wh;
		for(wh=0;wh<n;wh++) {
			if(EV[wh]>x) {
				--wh;
				break;
			}
		}
		if(wh==-1 || wh==n) error("cdf(): x lies outside original range");
		//# A piecewise linear approximation to the p.d.f.
	    return PR[wh]+(x-EV[wh])*G[wh];
	}
};

};	//namespace myutils

#endif//_APPROXDF_H_
