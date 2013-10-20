/*  Copyright 2012 Daniel Wilson.
 *
 *  ess.h
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
#ifndef _EFFECTIVE_SAMPLE_SIZE_H_
#define _EFFECTIVE_SAMPLE_SIZE_H_

#include <math.h>
#include <vector.h>

namespace myutils {
	inline double effectiveSampleSize(double* statistic, const int samples) {
		//int maxLag = samples;
		int maxLag = 1000;

		Vector<double> gammaStat(maxLag,0.0);
		Vector<double> varGammaStat(maxLag,0.0);
		double meanStat = 0.0;
  		double varStat,varVarStat,assVarCor,del1, del2;

		int i,j,lag;

  		for(i=0; i<samples; i++) meanStat += statistic[i];
		meanStat /= (double)samples;

		for(lag=0; lag<maxLag; lag++) {
			for(j=0; j<samples-lag; j++) {
    			del1=statistic[j] - meanStat;
				del2=statistic[j+lag] - meanStat;
				gammaStat[lag] += ( del1*del2 );
				varGammaStat[lag] += (del1*del1*del2*del2);
			}

			gammaStat[lag] /= (double)samples;
			varGammaStat[lag] /= (double)(samples-lag);
			varGammaStat[lag] -= (gammaStat[0] * gammaStat[0]);
		}

		varStat = gammaStat[0];
		varVarStat = varGammaStat[0];
		assVarCor = 1.0;

		lag=1;
		while ((lag < maxLag-3) && (gammaStat[lag] + gammaStat[lag+1] > 0)) {
			varStat += (2.0*(gammaStat[lag]+gammaStat[lag+1]));
			varVarStat += (2.0*(varGammaStat[lag] + varGammaStat[lag+1]));
			assVarCor += (2.0*((gammaStat[lag] * gammaStat[lag]) + (gammaStat[lag+1] * gammaStat[lag+1])) / (gammaStat[0] * gammaStat[0]));
			if (gammaStat[lag]+gammaStat[lag+1] < gammaStat[lag+2]+gammaStat[lag+3] ) break;
			lag += 2;
		}

		// standard error of mean
		double stdErrorOfMean = sqrt(varStat/samples);
		// variance of statistic
		double variance = gammaStat[0];
		// standard error of variance
		double stdErrorOfVariance = sqrt(varVarStat/samples);
		// effective sample size
		double ESS = gammaStat[0] * samples / varStat;
		// M
		int M = lag;
		return ESS;
	}
};

#endif//_EFFECTIVE_SAMPLE_SIZE_H_