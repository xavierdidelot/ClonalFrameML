/*  Copyright 2012 Daniel Wilson.
 *
 *  lgamma.h
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
#ifndef _LGAMMA_H_
#define _LGAMMA_H_

#include <limits>

namespace myutils {
/* From Numerical Recipes in C++ */
/*double lgamma(const double xx) {
	if(xx==1.0) return 0.0;
	int j;
	double x,y,tmp,ser;
	static const double cof[6] = {76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
		-0.5395239384953e-5};

	y = x = xx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;
	for(j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp + log(2.5066282746310005*ser/x);
}*/
//double lbeta(const double x, const double y) {
//	return lgamma(x)+lgamma(y)-lgamma(x+y);
//}
// gcf, gser and gammp from Numerical Recipes in C++ (Press et al)
void gcf(double &gammcf, const double a, const double x, double &gln)
{
	const int ITMAX=100;
	const double EPS=numeric_limits<double>::epsilon();
	const double FPMIN=numeric_limits<double>::min()/EPS;
	int i;
	double an,b,c,d,del,h;

	gln=lgamma(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS) break;
	}
	if (i > ITMAX) error("myutils::gcf(): a too large, ITMAX too small");
	gammcf=exp(-x+a*log(x)-gln)*h;
}
void gser(double &gamser, const double a, const double x, double &gln)
{
	const int ITMAX=100;
	const double EPS=numeric_limits<double>::epsilon();
	int n;
	double sum,del,ap;

	gln=lgamma(a);
	if (x <= 0.0) {
		if (x < 0.0) error("myutils::gser(): x less than 0");
		gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=0;n<ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				gamser=sum*exp(-x+a*log(x)-gln);
				return;
			}
		}
		error("myutils::gser(): a too large, ITMAX too small");
		return;
	}
}
double gammp(const double a, const double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0)
		error("myutils::gammp(): Invalid arguments");
	if (x < a+1.0) {
		gser(gamser,a,x,gln);
		return gamser;
	} else {
		gcf(gammcf,a,x,gln);
		return 1.0-gammcf;
	}
}
// Overloaded lgamma with 2 arguments becomes the incomplete gamma function
double lgamma(const double a, const double x) {
	return log(gammp(a,x));
}
};

#endif//_LGAMMA_H_