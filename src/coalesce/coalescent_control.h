/*  Copyright 2013 Daniel Wilson.
 *
 *  coalescent_control.h
 *  Part of the coalesce library.
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
#ifndef _CONTROL_H_
#define _CONTROL_H_

#pragma warning(disable: 4786)

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

using namespace std;

#include "myutils/matrix.h"
#include "myutils/controlwizard.h"
using myutils::Matrix;
using myutils::ControlWizard;
using myutils::TP_UNRECOGNISED;
using myutils::TP_INT;
using myutils::TP_DOUBLE;
using myutils::TP_STRING;
using myutils::TP_VEC_INT;
using myutils::TP_VEC_DOUBLE;
using myutils::TP_EXT_VEC_DOUBLE;
//using myutils::DATA_TYPE;

class Control
{
public:
	int nsamp;							// Sample size
	vector<double> ntimes;				// Times (in Ne gens) of samples, ordering unimportant
	double Negens;						// Expresses 1 unit of Ne gens in the same units as ntimes
	int loci;							// Number of independent loci simulated
	int seq_len;						// Total length of sequences simulated
	vector<int> len;					// Lengths for each locus
	double r;							// Per site rate of crossing-over per Negens (standard model, lambda = 0)
										// or TWICE the per site rate of initiation of recombination per Negens
										// (bacterial model, lambda > 0)
	vector<double> rmap;				// Map for heterogeneous recombination rates
	double lambda;						// 1/mean tract length

	double M;							// headline per site mutation rate
	int n_states;						// number of states (e.g. 4 nucleotides, 64 codons)
	vector<double> state_freq;			// initial state frequencies
	vector<double> state_rel_mut_rate;	// mutation rates relative to headline mutation rate
	vector<double> state_M;				// 1/(state-specific per site mutation rate)
										// (to be calculated)
	Matrix<double> mut_matrix;			// transition matrix for the states
	vector<char> state_name;			// letters for the states

	int nruns;							// store this in control class also
	int update_interval;

	bool coutput;

	/* Variables for structured coalescent */
	int ndemes;
	vector<int> deme_config;			// for each sample member, their starting deme (0..ndemes-1)
	Matrix<double> mig;					// ndemes * ndemes matrix: double the backwards in time migration rate from i to j
	vector<double> N_deme_over_D;		// for each deme, the pop size relative to the total
	
public:
	Control()
	{
		coutput=true;
		/* Set defaults */
		nsamp = 2;
		ntimes = vector<double>(2,0.0);
		Negens = 1.0;
		loci = 1;
		seq_len = 1;
		len = vector<int>(1,1);
		r = lambda = M = 0.0;
		rmap = vector<double>(0);
		ndemes = 0;
	}
	Control& read_input(char* filename)
	{
		seq_len=-1;
		Negens=1.0; //By default
		
		ControlWizard control_file;
		control_file.coutput=coutput;
		control_file.add_ITEM("n",TP_INT,&nsamp);
		control_file.add_ITEM("ntimes",TP_VEC_DOUBLE,&ntimes);
		control_file.add_item("Negens",TP_DOUBLE,&Negens);
		control_file.add_ITEM("loci",TP_INT,&loci);
		control_file.add_ITEM("len",TP_VEC_INT,&len);
		control_file.add_ITEM("lambda",TP_DOUBLE,&lambda);
		control_file.add_ITEM("n_states",TP_INT,&n_states);
		control_file.add_ITEM("mu",TP_DOUBLE,&M);
		control_file.add_ITEM("r",TP_DOUBLE,&r);
		vector<double> temp_mut_matrix;
		control_file.add_ITEM("mut_matrix",TP_EXT_VEC_DOUBLE,&temp_mut_matrix);
		control_file.add_item("state_rel_mut_rate",TP_EXT_VEC_DOUBLE,&state_rel_mut_rate);
		control_file.add_item("state_M",TP_EXT_VEC_DOUBLE,&state_M);
		control_file.add_ITEM("state_freq",TP_EXT_VEC_DOUBLE,&state_freq);
		control_file.add_item("nruns",TP_INT,&nruns);
		control_file.add_item("seq_len",TP_INT,&seq_len);
		control_file.add_ITEM("update_interval",TP_INT,&update_interval);
		control_file.read_input(filename);
		if(coutput)control_file.check_required();
		else
			if(!control_file.got_required)error("Not all necessary items found in control file");

		/*Check for necessary parameters*/
		if(!control_file.got_required)error("read_input(): necessary parameters not found");
		
		vector_to_Matrix(&temp_mut_matrix,&mut_matrix,n_states,n_states);

		if(ntimes.size()!=nsamp)error("ntimes inconsistent in size with n");
		sort(ntimes.begin(),ntimes.end());
		if(coutput) {
			int o;
			for(o=0;o<(int)ntimes.size();o++)printf("%g ",ntimes[o]);
		}

		if(len.size()!=loci)error("len inconsistent in size with loci");
		if(state_freq.size()!=n_states)error("state_freq inconsistent in size with n_states");

		if(update_interval==0)update_interval=1;

		int ind_seq_len=0;
		int i;
		for(i=0;i<loci;i++)ind_seq_len+=len[i];
		if(seq_len==-1)seq_len=ind_seq_len;
		else if(seq_len!=ind_seq_len) error("seq_len and the sum of len do not equate");

		if(coutput)for(i=0;i<nsamp;i++){ntimes[i]/=Negens;printf("%g ",ntimes[i]);}
		M*=Negens;
		r*=Negens;

		bool one_or_the_other=false;
		bool got_state_M=false;
		if(state_M.size()>0)
		{
			one_or_the_other=true;
			got_state_M=true;
		}
		if(state_rel_mut_rate.size()>0)
			one_or_the_other=true;
		if(!one_or_the_other)error("read_input(): neither state_M or state_rel_mut_rate received");
		if(got_state_M)
		{
			if(state_M.size()!=n_states)error("state_M inconsistent in size with n_states");
			state_rel_mut_rate.resize(n_states);
			int i;
			for(i=0;i<n_states;i++)
				state_rel_mut_rate[i]=1.0/(state_M[i]*M);
		}
		else
		{
			if(state_rel_mut_rate.size()!=n_states)error("state_rel_mut_rate inconsistent in size with n_states");
			state_M.resize(n_states);
			int i;
			for(i=0;i<n_states;i++)
				state_M[i]=1.0/(state_rel_mut_rate[i]*M);
		}

		return *this;
	}
	void error(char* error_text)
	{
		printf("\nRun-time error in Control::");
		printf("%s%\n", error_text);
		printf("Exiting to system...\n");
		exit(13);
	}
	Control& display_params()
	{
		printf("\nParameters\n\n");
		printf("Sample size, n = %d\n",nsamp);
		printf("Number of loci = %d\n",loci);
		printf("Locus lengths  = "); int i; for(i=0;i<(int)len.size();i++)printf("%d ",len[i]); printf("\n");
		printf("Sequence length = %d\n",seq_len);
		printf("Number of states = %d\n",n_states);
		printf("Per site recombination rate (scaled time) = %g\n",r);
		printf("Mean tract length (per codon) = %g\n",1/lambda);
		printf("Headline per site mutation rate (scaled time) = %g\n",M);
		printf("Number of runs = %d\n",nruns);
		return *this;
	}
	Control& vector_to_Matrix(vector<double> *vec, Matrix<double> *mat, int rows, int cols)
	{
		if((int)vec->size()<(rows*cols))error("vector_to_Matrix(): vector too small to fill matrix");
		if((int)vec->size()>(rows*cols))error("vector_to_Matrix(): vector too large to fit matrix");

		mat->resize(rows,cols);

		int current_row=0;
		int i,j;
		for(i=0;i<(int)vec->size();i+=cols)
		{
			for(j=0;j<cols;j++)
			{
				mat->element[current_row][j]=vec->at(i+j);
			}
			++current_row;
		}

		return *this;
	}
	Control& read_mut_matrix(Matrix<double> &G, vector<double> &pi)
	{
		/*Check it is a rate matrix*/
		int i,j;
		if(G.nrows()!=G.ncols())error("read_mut_matrix(): not a square matrix");
		for(i=0;i<G.nrows();i++)
		{
			double rowsum=0.0;
			for(j=0;j<G.ncols();j++)
				if(i!=j)
				{
					if(G[i][j]<0.0)error("read_mut_matrix(): negative rates in matrix");
					rowsum+=G[i][j];
				}
			if(rowsum!=-G[i][i]) error("read_mut_matrix(): not a rate matrix");
		}
		n_states=G.nrows();
		mut_matrix.resize(n_states,n_states);
		/*Check that the state frequencies are okay*/
		if(pi.size()!=n_states)error("read_mut_matrix(): pi inconsistent in size with G");
		double pisum=0.0;
		for(i=0;i<n_states;i++)
		{
			if(pi[i]<0.0)error("read_mut_matrix(): negative equilibrium frequencies");
			pisum+=pi[i];
		}
		if(pisum!=1.0)error("read_mut_matrix(): pi does not sum to one");
		/*Perform the copying*/
		state_freq=pi;
		state_rel_mut_rate.resize(n_states,0.0);
		state_M.resize(n_states,0.0);
		M=0.0;
		for(i=0;i<n_states;i++)M-=G[i][i]*state_freq[i];
		for(i=0;i<n_states;i++)
		{
			for(j=0;j<n_states;j++)
				mut_matrix[i][j]=-G[i][j]/G[i][i];
			mut_matrix[i][i]=0.0;
			state_rel_mut_rate[i]=-G[i][i]/M;
			state_M[i]=-1.0/G[i][i];
		}
		return *this;
	}
};


#endif
