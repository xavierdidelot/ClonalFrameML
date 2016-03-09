/*  Copyright 2013 Daniel Wilson.
 *
 *  mutation.h
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
#ifndef _MUTATION_H_
#define _MUTATION_H_

#include <vector>
#include "myutils/myutils.h"

using std::vector;
using namespace myutils;

class Mutation_Matrix {
public:
	int n_states;
	vector<double> state_freq;
	vector<string> state_char;
	Random *ran;

	Matrix<double> C;			/*continuous-time rate matrix*/
	Matrix<double> D;			/*discrete-time transition matrix*/
	vector<double> mutation_rate;
	vector<double> mutation_mean;
	
protected:
	/*assumes C_in is a valid rate matrix of size n_states_in*/
	void initialize(const int n_states_in, Matrix<double> *C_in) {
		n_states = n_states_in;
		C = *C_in;
		D.resize(n_states,n_states);
		mutation_rate.resize(n_states);
		mutation_mean.resize(n_states);
		int i,j;
		for(i=0;i<n_states;i++) {
			for(j=0;j<n_states;j++)
				D[i][j] = -C[i][j]/C[i][i];
			D[i][i] = 0.0;
			mutation_rate[i] = -C[i][i];
			mutation_mean[i] = -1./C[i][i];
		}
	}

public:
	Mutation_Matrix() {};
	Mutation_Matrix(const Mutation_Matrix &M) {
		n_states = M.n_states;
		state_freq = M.state_freq;
		state_char = M.state_char;
		ran = M.ran;
		C = M.C;
		D = M.D;
		mutation_rate = M.mutation_rate;
		mutation_mean = M.mutation_mean;
	}
	Mutation_Matrix& operator=(const Mutation_Matrix& M) {
		//if(this==&M) return *this;
		n_states = M.n_states;
		state_freq = M.state_freq;
		state_char = M.state_char;
		ran = M.ran;
		C = M.C;
		D = M.D;
		mutation_rate = M.mutation_rate;
		mutation_mean = M.mutation_mean;
		return *this;
	}
	void set_state_freq(vector<double> state_freq_in) {
		if(state_freq_in.size()!=n_states) error("Mutation_Matrix::set_state_freq(state_freq_in): state_freq_in inconsistent in size with n_states");
		double tot = 0.0;
		int i;
		for(i=0;i<n_states;i++) tot += state_freq_in[i];
		state_freq.resize(n_states);
		for(i=0;i<n_states;i++) state_freq[i]=state_freq_in[i]/tot;
	}
	void set_state_char(vector<string> state_char_in) {
		if(state_char_in.size()!=n_states) error("Mutation_Matrix::set_state_char(state_char_in): state_char_in inconsistent in size with n_states");
		state_char = state_char_in;
	}
	void set_ran(Random *ran_in) {
		ran = ran_in;
	}
	inline double get_rate(const int state) {
		if(state<0||state>=n_states) error("Mutation_Matrix::get_rate(state): unknown state");
		return mutation_rate[state];
	}
	int draw() {
		int state;
		double rp;
		double U = ran->U();
		for(state=0;state<n_states-1;state++) {
			rp = state_freq[state];
			if(U < rp) break;
			U -= rp;
		}
		return state;
	}
	int mutate(const int state) {
		if(state<0||state>=n_states) error("Mutation_Matrix::mutate(state): unknown state");
		int new_state = state;
		double U = ran->U();
		double rp;
		for(new_state=0;new_state<n_states-1;new_state++) {
			rp = D[state][new_state];
			if(U < rp) break;
			U -= rp;
		}
		return new_state;
	};
	int mutate_edge(const int state, const double time) {
		int old_state = state;
		if(old_state<0||old_state>=n_states) error("Mutation_Matrix::mutate_edge(state,time): unknown state");
		if(time<0.0) error("Mutation_Matrix::mutate_edge(state,time): time must be non-negative");
		double time_remaining = time;
		double next_mutation = ran->exponential(mutation_mean[old_state]);
		double U,rp;
		int new_state = state;
		while(next_mutation<time_remaining) {
			U = ran->U();
			for(new_state=0;new_state<n_states-1;new_state++) {
				rp = D[old_state][new_state];
				if(U < rp) break;
				U -= rp;
			}
			time_remaining -= next_mutation;
			next_mutation = ran->exponential(mutation_mean[old_state]);
		}
		return new_state;
	}
	int mutate_edge(const int state, const double time, int &nmut) {
		int old_state = state;
		if(old_state<0||old_state>=n_states) error("Mutation_Matrix::mutate_edge(state,time): unknown state");
		if(time<0.0) error("Mutation_Matrix::mutate_edge(state,time): time must be non-negative");
		double time_remaining = time;
		double next_mutation = ran->exponential(mutation_mean[old_state]);
		double U,rp;
		int new_state = state;
		while(next_mutation<time_remaining) {
			U = ran->U();
			for(new_state=0;new_state<n_states-1;new_state++) {
				rp = D[old_state][new_state];
				if(U < rp) break;
				U -= rp;
			}
			time_remaining -= next_mutation;
			next_mutation = ran->exponential(mutation_mean[old_state]);
			++nmut;
		}
		return new_state;
	}
	double expected_rate() {
		double result = 0.0;
		int i;
		for(i=0;i<n_states;i++)
			result -= state_freq[i] * C[i][i];
		return result;
	}
};

class Nucleotide_Mutation_Matrix : public Mutation_Matrix {
public:
	void set_defaults() {
		n_states = 4;
		state_freq.resize(n_states,0.25);
		state_char.resize(n_states);
		state_char[0] = string(1,'A');
		state_char[1] = string(1,'G');
		state_char[2] = string(1,'C');
		state_char[3] = string(1,'T');
	}
};

class JC69 : public Nucleotide_Mutation_Matrix {
public:
	JC69(const double lambda, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		C.resize(n_states,n_states);
		int i,j;
		for(i=0;i<4;i++) {
			C[i][i] = 0.0;
			for(j=0;j<4;j++)
				if(i!=j) {
					C[i][j] = lambda/4.;
					C[i][i] -= C[i][j];
				}
		}
		initialize(n_states,&C);
	}
	int fast_mutate(const int state, const double time) {
		double lambda = 4.*C[0][1];
		double p0 = .25 + .75*exp(-lambda*time);
		double p1 = .25 - .25*exp(-lambda*time);
		double p2 = p1;
		p1 += p0;
		p2 += p1;
		double U = ran->U();
		switch(state) {
		case 0:
			if(U<p0) return 0;
			if(U<p1) return 2;
			if(U<p2) return 3;
			return 1;
		case 1:
			if(U<p0) return 1;
			if(U<p1) return 2;
			if(U<p2) return 3;
			return 0;
		case 2:
			if(U<p0) return 2;
			if(U<p1) return 0;
			if(U<p2) return 1;
			return 3;
		case 3:
			if(U<p0) return 3;
			if(U<p1) return 0;
			if(U<p2) return 1;
			return 2;
		}
		error("JC69::fast_mutate(): state not recognised");
		return -1;
	}
};

class F81 : public Nucleotide_Mutation_Matrix {
public:
	F81(const double lambda, const vector<double> state_freq_in, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		set_state_freq(state_freq_in);
		C.resize(n_states,n_states);
		int i,j;
		for(i=0;i<4;i++) {
			C[i][i] = 0.0;
			for(j=0;j<4;j++)
				if(i!=j) {
					C[i][j] = lambda*state_freq[j];
					C[i][i] -= C[i][j];
				}
		}
		initialize(n_states,&C);
	}
};

class K80 : public Nucleotide_Mutation_Matrix {
public:
	K80(const double lambda, const double kappa, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		C.resize(n_states,n_states);
		update(lambda,kappa);
	}
	K80& update(const double lambda, const double kappa) {
		int i,j;
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				C[i][j] = (i==j) ? 0.0 : lambda/4.;
		C[0][1]*=kappa;C[1][0]*=kappa;
		C[2][3]*=kappa;C[3][2]*=kappa;
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				if(i!=j)
					C[i][i] -= C[i][j];
		initialize(n_states,&C);
		return *this;
	}
	int fast_mutate(const int state, const double time) {
		double lambda = 4.*C[0][2];
		double kappa = C[0][1]/C[0][2];
		double p0 = .25 + .25*exp(-2.*lambda*time) + .5*exp(-lambda*time*(1.+kappa));
		double p1 = p0 + .25 - .25*exp(-2.*lambda*time);
		double p2 = p1 + .25 - .25*exp(-2.*lambda*time);
		double U = ran->U();
		switch(state) {
		case 0:
			if(U<p0) return 0;
			if(U<p1) return 2;
			if(U<p2) return 3;
			return 1;
		case 1:
			if(U<p0) return 1;
			if(U<p1) return 2;
			if(U<p2) return 3;
			return 0;
		case 2:
			if(U<p0) return 2;
			if(U<p1) return 0;
			if(U<p2) return 1;
			return 3;
		case 3:
			if(U<p0) return 3;
			if(U<p1) return 0;
			if(U<p2) return 1;
			return 2;
		}
		error("K80::fast_mutate(): state not recognised");
		return -1;
	}
};

class HKY85 : public Nucleotide_Mutation_Matrix {
public:
	HKY85(const double lambda, const double kappa, const vector<double> state_freq_in, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		set_state_freq(state_freq_in);
		C.resize(n_states,n_states);
		int i,j;
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				C[i][j] = (i==j) ? 0.0 : lambda*state_freq[j];
		C[0][1]*=kappa;C[1][0]*=kappa;
		C[2][3]*=kappa;C[3][2]*=kappa;
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				if(i!=j)
					C[i][i] -= C[i][j];
		initialize(n_states,&C);
	}
};

class TN93 : public Nucleotide_Mutation_Matrix {
public:
	TN93(const double lambda, const double kappa_R, const double kappa_Y, const vector<double> state_freq_in, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		set_state_freq(state_freq_in);
		C.resize(n_states,n_states);
		int i,j;
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				C[i][j] = (i==j) ? 0.0 : lambda*state_freq[j];
		C[0][1]*=kappa_R;C[1][0]*=kappa_R;
		C[2][3]*=kappa_Y;C[3][2]*=kappa_Y;
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				if(i!=j)
					C[i][i] -= C[i][j];
		initialize(n_states,&C);
	}
};

class Codon_Mutation_Matrix : public Mutation_Matrix {
public:
	void set_defaults() {
		n_states = 64;
		state_freq.resize(n_states,1./61.);
		state_freq[10]=state_freq[11]=state_freq[14]=0.0;
		string default_char = string(3,'-');
		state_char.resize(n_states,default_char);
		int i,j,k,l;
		vector<char> base(4,'-');
		base[0] = 'U'; base[1] = 'C'; base[2] = 'A'; base[3] = 'G';
		for(i=0,l=0;i<4;i++)
			for(j=0;j<4;j++)
				for(k=0;k<4;k++,l++) {
					state_char[l][0] = base[i];
					state_char[l][1] = base[j];
					state_char[l][2] = base[k];
				}
		//for(i=0;i<64;i++) state_char[i] = i+1;
	}
	virtual Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega) = 0;
	virtual Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega, const vector<double> &pi) = 0;
	virtual Codon_Mutation_Matrix& build_C(const double mu, const double kappa, const double omega, const vector<double> &pi) = 0;

};

class NY98 : public Codon_Mutation_Matrix {
protected:
	NY98() {}
public:
	NY98(Random *ran_in) {
		set_defaults();
		ran = ran_in;
		C.resize(n_states,n_states);
		D.resize(n_states,n_states);
	}
	NY98(const double mu, const double kappa, const double omega, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		C.resize(n_states,n_states);
		D.resize(n_states,n_states);
		update(mu,kappa,omega,state_freq);
	}
	NY98(const double mu, const double kappa, const double omega, const vector<double> &pi, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		set_state_freq(pi);
		C.resize(n_states,n_states);
		D.resize(n_states,n_states);
		update(mu,kappa,omega,state_freq);
	}
	Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega) {
		build_C(mu,kappa,omega,state_freq);
		initialize(n_states,&C);
		return *this;
	}
	Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega, const vector<double> &pi) {
		set_state_freq(pi);
		build_C(mu,kappa,omega,state_freq);
		initialize(n_states,&C);
		return *this;
	}
	Codon_Mutation_Matrix& build_C(const double mu, const double kappa, const double omega, const vector<double> &pi)
	{
		int i,j;
		/*Initialize to zero*/
		for(i=0;i<64;i++){for(j=0;j<64;j++)C[i][j]=0.0;}
		
		C[0][1]=kappa*mu;
		C[0][2]=omega*mu;
		C[0][3]=omega*mu;
		C[0][4]=kappa*omega*mu;
		C[0][8]=omega*mu;
		C[0][12]=omega*mu;
		C[0][16]=kappa*omega*mu;
		C[0][32]=omega*mu;
		C[0][48]=omega*mu;
		C[1][2]=omega*mu;
		C[1][3]=omega*mu;
		C[1][5]=kappa*omega*mu;
		C[1][9]=omega*mu;
		C[1][13]=omega*mu;
		C[1][17]=kappa*omega*mu;
		C[1][33]=omega*mu;
		C[1][49]=omega*mu;
		C[2][3]=kappa*mu;
		C[2][6]=kappa*omega*mu;
		C[2][10]=omega*mu;
		C[2][14]=omega*mu;
		C[2][18]=kappa*mu;	/*Synonymous!*/
		C[2][34]=omega*mu;
		C[2][50]=omega*mu;
		C[3][7]=kappa*omega*mu;
		C[3][11]=omega*mu;
		C[3][15]=omega*mu;
		C[3][19]=kappa*mu;	/*Synonymous!*/
		C[3][35]=omega*mu;
		C[3][51]=omega*mu;
		C[4][5]=kappa*mu;
		C[4][6]=mu;
		C[4][7]=mu;
		C[4][8]=omega*mu;
		C[4][12]=omega*mu;
		C[4][20]=kappa*omega*mu;
		C[4][36]=omega*mu;
		C[4][52]=omega*mu;
		C[5][6]=mu;
		C[5][7]=mu;
		C[5][9]=omega*mu;
		C[5][13]=omega*mu;
		C[5][21]=kappa*omega*mu;
		C[5][37]=omega*mu;
		C[5][53]=omega*mu;
		C[6][7]=kappa*mu;
		C[6][10]=omega*mu;
		C[6][14]=omega*mu;
		C[6][22]=kappa*omega*mu;
		C[6][38]=omega*mu;
		C[6][54]=omega*mu;
		C[7][11]=omega*mu;
		C[7][15]=omega*mu;
		C[7][23]=kappa*omega*mu;
		C[7][39]=omega*mu;
		C[7][55]=omega*mu;
		C[8][9]=kappa*mu;
		C[8][10]=omega*mu;
		C[8][11]=omega*mu;
		C[8][12]=kappa*omega*mu;
		C[8][24]=kappa*omega*mu;
		C[8][40]=omega*mu;
		C[8][56]=omega*mu;
		C[9][10]=omega*mu;
		C[9][11]=omega*mu;
		C[9][13]=kappa*omega*mu;
		C[9][25]=kappa*omega*mu;
		C[9][41]=omega*mu;
		C[9][57]=omega*mu;
		C[10][11]=kappa*mu;
		C[10][14]=kappa*mu;
		C[10][26]=kappa*omega*mu;
		C[10][42]=omega*mu;
		C[10][58]=omega*mu;
		C[11][15]=kappa*omega*mu;
		C[11][27]=kappa*omega*mu;
		C[11][43]=omega*mu;
		C[11][59]=omega*mu;
		C[12][13]=kappa*mu;
		C[12][14]=omega*mu;
		C[12][15]=omega*mu;
		C[12][28]=kappa*omega*mu;
		C[12][44]=omega*mu;
		C[12][60]=omega*mu;
		C[13][14]=omega*mu;
		C[13][15]=omega*mu;
		C[13][29]=kappa*omega*mu;
		C[13][45]=omega*mu;
		C[13][61]=omega*mu;
		C[14][15]=kappa*omega*mu;
		C[14][30]=kappa*omega*mu;
		C[14][46]=omega*mu;
		C[14][62]=omega*mu;
		C[15][31]=kappa*omega*mu;
		C[15][47]=omega*mu;
		C[15][63]=omega*mu;
		C[16][17]=kappa*mu;
		C[16][18]=mu;
		C[16][19]=mu;
		C[16][20]=kappa*omega*mu;
		C[16][24]=omega*mu;
		C[16][28]=omega*mu;
		C[16][32]=omega*mu;
		C[16][48]=omega*mu;
		C[17][18]=mu;
		C[17][19]=mu;
		C[17][21]=kappa*omega*mu;
		C[17][25]=omega*mu;
		C[17][29]=omega*mu;
		C[17][33]=omega*mu;
		C[17][49]=omega*mu;
		C[18][19]=kappa*mu;
		C[18][22]=kappa*omega*mu;
		C[18][26]=omega*mu;
		C[18][30]=omega*mu;
		C[18][34]=omega*mu;
		C[18][50]=omega*mu;
		C[19][23]=kappa*omega*mu;
		C[19][27]=omega*mu;
		C[19][31]=omega*mu;
		C[19][35]=omega*mu;
		C[19][51]=omega*mu;
		C[20][21]=kappa*mu;
		C[20][22]=mu;
		C[20][23]=mu;
		C[20][24]=omega*mu;
		C[20][28]=omega*mu;
		C[20][36]=omega*mu;
		C[20][52]=omega*mu;
		C[21][22]=mu;
		C[21][23]=mu;
		C[21][25]=omega*mu;
		C[21][29]=omega*mu;
		C[21][37]=omega*mu;
		C[21][53]=omega*mu;
		C[22][23]=kappa*mu;
		C[22][26]=omega*mu;
		C[22][30]=omega*mu;
		C[22][38]=omega*mu;
		C[22][54]=omega*mu;
		C[23][27]=omega*mu;
		C[23][31]=omega*mu;
		C[23][39]=omega*mu;
		C[23][55]=omega*mu;
		C[24][25]=kappa*mu;
		C[24][26]=omega*mu;
		C[24][27]=omega*mu;
		C[24][28]=kappa*omega*mu;
		C[24][40]=omega*mu;
		C[24][56]=omega*mu;
		C[25][26]=omega*mu;
		C[25][27]=omega*mu;
		C[25][29]=kappa*omega*mu;
		C[25][41]=omega*mu;
		C[25][57]=omega*mu;
		C[26][27]=kappa*mu;
		C[26][30]=kappa*omega*mu;
		C[26][42]=omega*mu;
		C[26][58]=omega*mu;
		C[27][31]=kappa*omega*mu;
		C[27][43]=omega*mu;
		C[27][59]=omega*mu;
		C[28][29]=kappa*mu;
		C[28][30]=mu;
		C[28][31]=mu;
		C[28][44]=omega*mu;
		C[28][60]=omega*mu;
		C[29][30]=mu;
		C[29][31]=mu;
		C[29][45]=omega*mu;
		C[29][61]=omega*mu;
		C[30][31]=kappa*mu;
		C[30][46]=mu;
		C[30][62]=omega*mu;
		C[31][47]=mu;
		C[31][63]=omega*mu;
		C[32][33]=kappa*mu;
		C[32][34]=mu;
		C[32][35]=omega*mu;
		C[32][36]=kappa*omega*mu;
		C[32][40]=omega*mu;
		C[32][44]=omega*mu;
		C[32][48]=kappa*omega*mu;
		C[33][34]=mu;
		C[33][35]=omega*mu;
		C[33][37]=kappa*omega*mu;
		C[33][41]=omega*mu;
		C[33][45]=omega*mu;
		C[33][49]=kappa*omega*mu;
		C[34][35]=kappa*omega*mu;
		C[34][38]=kappa*omega*mu;
		C[34][42]=omega*mu;
		C[34][46]=omega*mu;
		C[34][50]=kappa*omega*mu;
		C[35][39]=kappa*omega*mu;
		C[35][43]=omega*mu;
		C[35][47]=omega*mu;
		C[35][51]=kappa*omega*mu;
		C[36][37]=kappa*mu;
		C[36][38]=mu;
		C[36][39]=mu;
		C[36][40]=omega*mu;
		C[36][44]=omega*mu;
		C[36][52]=kappa*omega*mu;
		C[37][38]=mu;
		C[37][39]=mu;
		C[37][41]=omega*mu;
		C[37][45]=omega*mu;
		C[37][53]=kappa*omega*mu;
		C[38][39]=kappa*mu;
		C[38][42]=omega*mu;
		C[38][46]=omega*mu;
		C[38][54]=kappa*omega*mu;
		C[39][43]=omega*mu;
		C[39][47]=omega*mu;
		C[39][55]=kappa*omega*mu;
		C[40][41]=kappa*mu;
		C[40][42]=omega*mu;
		C[40][43]=omega*mu;
		C[40][44]=kappa*omega*mu;
		C[40][56]=kappa*omega*mu;
		C[41][42]=omega*mu;
		C[41][43]=omega*mu;
		C[41][45]=kappa*omega*mu;
		C[41][57]=kappa*omega*mu;
		C[42][43]=kappa*mu;
		C[42][46]=kappa*omega*mu;
		C[42][58]=kappa*omega*mu;
		C[43][47]=kappa*omega*mu;
		C[43][59]=kappa*omega*mu;
		C[44][45]=kappa*mu;
		C[44][46]=omega*mu;
		C[44][47]=omega*mu;
		C[44][60]=kappa*omega*mu;
		C[45][46]=omega*mu;
		C[45][47]=omega*mu;
		C[45][61]=kappa*omega*mu;
		C[46][47]=kappa*mu;
		C[46][62]=kappa*omega*mu;
		C[47][63]=kappa*omega*mu;
		C[48][49]=kappa*mu;
		C[48][50]=mu;
		C[48][51]=mu;
		C[48][52]=kappa*omega*mu;
		C[48][56]=omega*mu;
		C[48][60]=omega*mu;
		C[49][50]=mu;
		C[49][51]=mu;
		C[49][53]=kappa*omega*mu;
		C[49][57]=omega*mu;
		C[49][61]=omega*mu;
		C[50][51]=kappa*mu;
		C[50][54]=kappa*omega*mu;
		C[50][58]=omega*mu;
		C[50][62]=omega*mu;
		C[51][55]=kappa*omega*mu;
		C[51][59]=omega*mu;
		C[51][63]=omega*mu;
		C[52][53]=kappa*mu;
		C[52][54]=mu;
		C[52][55]=mu;
		C[52][56]=omega*mu;
		C[52][60]=omega*mu;
		C[53][54]=mu;
		C[53][55]=mu;
		C[53][57]=omega*mu;
		C[53][61]=omega*mu;
		C[54][55]=kappa*mu;
		C[54][58]=omega*mu;
		C[54][62]=omega*mu;
		C[55][59]=omega*mu;
		C[55][63]=omega*mu;
		C[56][57]=kappa*mu;
		C[56][58]=omega*mu;
		C[56][59]=omega*mu;
		C[56][60]=kappa*omega*mu;
		C[57][58]=omega*mu;
		C[57][59]=omega*mu;
		C[57][61]=kappa*omega*mu;
		C[58][59]=kappa*mu;
		C[58][62]=kappa*omega*mu;
		C[59][63]=kappa*omega*mu;
		C[60][61]=kappa*mu;
		C[60][62]=mu;
		C[60][63]=mu;
		C[61][62]=mu;
		C[61][63]=mu;
		C[62][63]=kappa*mu;

		/*Remove the STOP codons from the scheme*/
		for(i=0;i<64;i++){
			C[10][i]=0.0;	C[i][10]=0.0;
			C[11][i]=0.0;	C[i][11]=0.0;
			C[14][i]=0.0;	C[i][14]=0.0;
		}

		/*Fill in the lower triangle*/
		for(i=0;i<64;i++){
			for(j=i+1;j<64;j++)C[j][i]=C[i][j];}

		/*Apply the equilibrium frequencies*/
		for(i=0;i<64;i++)
			for(j=0;j<64;j++)
				C[i][j]*=pi[j];
		
		/*Compute the diagonal*/
		for(i=0;i<64;i++)
		{
			double rowsum=0.0;
			for(j=0;j<64;j++) rowsum+=C[i][j];
			C[i][i]=-rowsum;
		}
		return *this;
	}
};

class NY98_61 : public NY98 {
public:
	void set_defaults() {
		n_states = 61;
		state_freq.resize(n_states,1./61.);
		string default_char = string(3,'-');
		state_char.resize(n_states,default_char);
		int i,j,k,l,m;
		vector<char> base(4,'-');
		base[0] = 'U'; base[1] = 'C'; base[2] = 'A'; base[3] = 'G';
		for(i=0,l=0,m=0;i<4;i++)
			for(j=0;j<4;j++)
				for(k=0;k<4;k++,l++,m++) {
					state_char[m][0] = base[i];
					state_char[m][1] = base[j];
					state_char[m][2] = base[k];
					if(l==10 || l==11 || l==14) --m;
				}
		//for(i=0;i<61;i++) state_char[i] = i+1;
	}
	NY98_61(Random *ran_in) {
		set_defaults();
		ran = ran_in;
		C.resize(n_states,n_states);
		D.resize(n_states,n_states);
	}
	NY98_61(const double mu, const double kappa, const double omega, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		C.resize(n_states,n_states);
		D.resize(n_states,n_states);
		update(mu,kappa,omega,state_freq);
	}
	NY98_61(const double mu, const double kappa, const double omega, const vector<double> &pi, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		set_state_freq(pi);
		C.resize(n_states,n_states);
		D.resize(n_states,n_states);
		update(mu,kappa,omega,state_freq);
	}
	Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega) {
		build_C(mu,kappa,omega,state_freq);
		initialize(n_states,&C);
		return *this;
	}
	Codon_Mutation_Matrix& update(const double mu, const double kappa, const double omega, const vector<double> &pi) {
		set_state_freq(pi);
		build_C(mu,kappa,omega,state_freq);
		initialize(n_states,&C);
		return *this;
	}
	Codon_Mutation_Matrix& build_C(const double mu, const double kappa, const double omega, const vector<double> &pi)
	{
		int i,j;
		/*Initialize to zero*/
		for(i=0;i<61;i++){for(j=0;j<61;j++)C[i][j]=0.0;}
		
		C[0][1]=kappa*mu;
		C[0][2]=omega*mu;
		C[0][3]=omega*mu;
		C[0][4]=kappa*omega*mu;
		C[0][8]=omega*mu;
		C[0][10]=omega*mu;
		C[0][13]=kappa*omega*mu;
		C[0][29]=omega*mu;
		C[0][45]=omega*mu;
		C[1][2]=omega*mu;
		C[1][3]=omega*mu;
		C[1][5]=kappa*omega*mu;
		C[1][9]=omega*mu;
		C[1][11]=omega*mu;
		C[1][14]=kappa*omega*mu;
		C[1][30]=omega*mu;
		C[1][46]=omega*mu;
		C[2][3]=kappa*mu;
		C[2][6]=kappa*omega*mu;
		C[2][15]=kappa*mu;
		C[2][31]=omega*mu;
		C[2][47]=omega*mu;
		C[3][7]=kappa*omega*mu;
		C[3][12]=omega*mu;
		C[3][16]=kappa*mu;
		C[3][32]=omega*mu;
		C[3][48]=omega*mu;
		C[4][5]=kappa*mu;
		C[4][6]=mu;
		C[4][7]=mu;
		C[4][8]=omega*mu;
		C[4][10]=omega*mu;
		C[4][17]=kappa*omega*mu;
		C[4][33]=omega*mu;
		C[4][49]=omega*mu;
		C[5][6]=mu;
		C[5][7]=mu;
		C[5][9]=omega*mu;
		C[5][11]=omega*mu;
		C[5][18]=kappa*omega*mu;
		C[5][34]=omega*mu;
		C[5][50]=omega*mu;
		C[6][7]=kappa*mu;
		C[6][19]=kappa*omega*mu;
		C[6][35]=omega*mu;
		C[6][51]=omega*mu;
		C[7][12]=omega*mu;
		C[7][20]=kappa*omega*mu;
		C[7][36]=omega*mu;
		C[7][52]=omega*mu;
		C[8][9]=kappa*mu;
		C[8][10]=kappa*omega*mu;
		C[8][21]=kappa*omega*mu;
		C[8][37]=omega*mu;
		C[8][53]=omega*mu;
		C[9][11]=kappa*omega*mu;
		C[9][22]=kappa*omega*mu;
		C[9][38]=omega*mu;
		C[9][54]=omega*mu;
		C[10][11]=kappa*mu;
		C[10][12]=omega*mu;
		C[10][25]=kappa*omega*mu;
		C[10][41]=omega*mu;
		C[10][57]=omega*mu;
		C[11][12]=omega*mu;
		C[11][26]=kappa*omega*mu;
		C[11][42]=omega*mu;
		C[11][58]=omega*mu;
		C[12][28]=kappa*omega*mu;
		C[12][44]=omega*mu;
		C[12][60]=omega*mu;
		C[13][14]=kappa*mu;
		C[13][15]=mu;
		C[13][16]=mu;
		C[13][17]=kappa*omega*mu;
		C[13][21]=omega*mu;
		C[13][25]=omega*mu;
		C[13][29]=omega*mu;
		C[13][45]=omega*mu;
		C[14][15]=mu;
		C[14][16]=mu;
		C[14][18]=kappa*omega*mu;
		C[14][22]=omega*mu;
		C[14][26]=omega*mu;
		C[14][30]=omega*mu;
		C[14][46]=omega*mu;
		C[15][16]=kappa*mu;
		C[15][19]=kappa*omega*mu;
		C[15][23]=omega*mu;
		C[15][27]=omega*mu;
		C[15][31]=omega*mu;
		C[15][47]=omega*mu;
		C[16][20]=kappa*omega*mu;
		C[16][24]=omega*mu;
		C[16][28]=omega*mu;
		C[16][32]=omega*mu;
		C[16][48]=omega*mu;
		C[17][18]=kappa*mu;
		C[17][19]=mu;
		C[17][20]=mu;
		C[17][21]=omega*mu;
		C[17][25]=omega*mu;
		C[17][33]=omega*mu;
		C[17][49]=omega*mu;
		C[18][19]=mu;
		C[18][20]=mu;
		C[18][22]=omega*mu;
		C[18][26]=omega*mu;
		C[18][34]=omega*mu;
		C[18][50]=omega*mu;
		C[19][20]=kappa*mu;
		C[19][23]=omega*mu;
		C[19][27]=omega*mu;
		C[19][35]=omega*mu;
		C[19][51]=omega*mu;
		C[20][24]=omega*mu;
		C[20][28]=omega*mu;
		C[20][36]=omega*mu;
		C[20][52]=omega*mu;
		C[21][22]=kappa*mu;
		C[21][23]=omega*mu;
		C[21][24]=omega*mu;
		C[21][25]=kappa*omega*mu;
		C[21][37]=omega*mu;
		C[21][53]=omega*mu;
		C[22][23]=omega*mu;
		C[22][24]=omega*mu;
		C[22][26]=kappa*omega*mu;
		C[22][38]=omega*mu;
		C[22][54]=omega*mu;
		C[23][24]=kappa*mu;
		C[23][27]=kappa*omega*mu;
		C[23][39]=omega*mu;
		C[23][55]=omega*mu;
		C[24][28]=kappa*omega*mu;
		C[24][40]=omega*mu;
		C[24][56]=omega*mu;
		C[25][26]=kappa*mu;
		C[25][27]=mu;
		C[25][28]=mu;
		C[25][41]=omega*mu;
		C[25][57]=omega*mu;
		C[26][27]=mu;
		C[26][28]=mu;
		C[26][42]=omega*mu;
		C[26][58]=omega*mu;
		C[27][28]=kappa*mu;
		C[27][43]=mu;
		C[27][59]=omega*mu;
		C[28][44]=mu;
		C[28][60]=omega*mu;
		C[29][30]=kappa*mu;
		C[29][31]=mu;
		C[29][32]=omega*mu;
		C[29][33]=kappa*omega*mu;
		C[29][37]=omega*mu;
		C[29][41]=omega*mu;
		C[29][45]=kappa*omega*mu;
		C[30][31]=mu;
		C[30][32]=omega*mu;
		C[30][34]=kappa*omega*mu;
		C[30][38]=omega*mu;
		C[30][42]=omega*mu;
		C[30][46]=kappa*omega*mu;
		C[31][32]=kappa*omega*mu;
		C[31][35]=kappa*omega*mu;
		C[31][39]=omega*mu;
		C[31][43]=omega*mu;
		C[31][47]=kappa*omega*mu;
		C[32][36]=kappa*omega*mu;
		C[32][40]=omega*mu;
		C[32][44]=omega*mu;
		C[32][48]=kappa*omega*mu;
		C[33][34]=kappa*mu;
		C[33][35]=mu;
		C[33][36]=mu;
		C[33][37]=omega*mu;
		C[33][41]=omega*mu;
		C[33][49]=kappa*omega*mu;
		C[34][35]=mu;
		C[34][36]=mu;
		C[34][38]=omega*mu;
		C[34][42]=omega*mu;
		C[34][50]=kappa*omega*mu;
		C[35][36]=kappa*mu;
		C[35][39]=omega*mu;
		C[35][43]=omega*mu;
		C[35][51]=kappa*omega*mu;
		C[36][40]=omega*mu;
		C[36][44]=omega*mu;
		C[36][52]=kappa*omega*mu;
		C[37][38]=kappa*mu;
		C[37][39]=omega*mu;
		C[37][40]=omega*mu;
		C[37][41]=kappa*omega*mu;
		C[37][53]=kappa*omega*mu;
		C[38][39]=omega*mu;
		C[38][40]=omega*mu;
		C[38][42]=kappa*omega*mu;
		C[38][54]=kappa*omega*mu;
		C[39][40]=kappa*mu;
		C[39][43]=kappa*omega*mu;
		C[39][55]=kappa*omega*mu;
		C[40][44]=kappa*omega*mu;
		C[40][56]=kappa*omega*mu;
		C[41][42]=kappa*mu;
		C[41][43]=omega*mu;
		C[41][44]=omega*mu;
		C[41][57]=kappa*omega*mu;
		C[42][43]=omega*mu;
		C[42][44]=omega*mu;
		C[42][58]=kappa*omega*mu;
		C[43][44]=kappa*mu;
		C[43][59]=kappa*omega*mu;
		C[44][60]=kappa*omega*mu;
		C[45][46]=kappa*mu;
		C[45][47]=mu;
		C[45][48]=mu;
		C[45][49]=kappa*omega*mu;
		C[45][53]=omega*mu;
		C[45][57]=omega*mu;
		C[46][47]=mu;
		C[46][48]=mu;
		C[46][50]=kappa*omega*mu;
		C[46][54]=omega*mu;
		C[46][58]=omega*mu;
		C[47][48]=kappa*mu;
		C[47][51]=kappa*omega*mu;
		C[47][55]=omega*mu;
		C[47][59]=omega*mu;
		C[48][52]=kappa*omega*mu;
		C[48][56]=omega*mu;
		C[48][60]=omega*mu;
		C[49][50]=kappa*mu;
		C[49][51]=mu;
		C[49][52]=mu;
		C[49][53]=omega*mu;
		C[49][57]=omega*mu;
		C[50][51]=mu;
		C[50][52]=mu;
		C[50][54]=omega*mu;
		C[50][58]=omega*mu;
		C[51][52]=kappa*mu;
		C[51][55]=omega*mu;
		C[51][59]=omega*mu;
		C[52][56]=omega*mu;
		C[52][60]=omega*mu;
		C[53][54]=kappa*mu;
		C[53][55]=omega*mu;
		C[53][56]=omega*mu;
		C[53][57]=kappa*omega*mu;
		C[54][55]=omega*mu;
		C[54][56]=omega*mu;
		C[54][58]=kappa*omega*mu;
		C[55][56]=kappa*mu;
		C[55][59]=kappa*omega*mu;
		C[56][60]=kappa*omega*mu;
		C[57][58]=kappa*mu;
		C[57][59]=mu;
		C[57][60]=mu;
		C[58][59]=mu;
		C[58][60]=mu;
		C[59][60]=kappa*mu;

		/*Fill in the lower triangle*/
		for(i=0;i<61;i++){
			for(j=i+1;j<61;j++)C[j][i]=C[i][j];}

		/*Apply the equilibrium frequencies*/
		for(i=0;i<61;i++)
			for(j=0;j<61;j++)
				C[i][j]*=pi[j];
		
		/*Compute the diagonal*/
		for(i=0;i<61;i++)
		{
			double rowsum=0.0;
			for(j=0;j<61;j++) rowsum+=C[i][j];
			C[i][i]=-rowsum;
		}
		return *this;
	}
};

class FSM_Binary : public Mutation_Matrix {
/****************************************************************/
/*	Mutations occur at rate lambda/2 per unit time.				*/
/*																*/
/*	Transition probability matrix, given time t is				*/
/*																*/
/*	P[0,0] = P[1,1] = 1/2 + 1/2*exp(-lambda*t)					*/
/*	P[0,1] = P[1,0] = 1/2 - 1/2*exp(-lambda*t)					*/
/*																*/
/*	Reversible model, so Pr(observing unordered pair ab) = 		*/
/*			(2-delta[a,b])*pi[a] P[a,b]^(2t)					*/
/*	where delta is the Kronecker delta and pi the equilibrium	*/
/*	frequency which is 1/2.										*/
/*																*/
/*	So Pr(observing unordered pair ab|mrca at t)				*/
/*			=	1/4 + 1/4*exp(-lambda*t)	if a=b				*/
/*			or	1/2 - 1/2*exp(-lambda*t)	otherwise			*/
/*																*/
/*	Expected pairwise diversity in a coalescent model, where	*/
/*	time is measured in units of PNe generations (P is ploidy	*/
/*	Ne is effective population size), is lambda/(1+2*lambda)	*/
/*																*/
/****************************************************************/
public:
	void set_defaults() {
		n_states = 2;
		state_freq.resize(n_states,0.5);
		state_char.resize(n_states);
		state_char[0] = string(1,'0');
		state_char[1] = string(1,'1');
	}
	FSM_Binary(Random *ran_in) {
		set_defaults();
		ran = ran_in;
		C.resize(n_states,n_states);
		initialize(2,&C);
	}
	FSM_Binary(const double lambda, Random *ran_in) {
		set_defaults();
		ran = ran_in;
		C.resize(n_states,n_states);
		update(lambda);
	}
	FSM_Binary& update(const double lambda) {
		C[0][0] = -lambda/2.; C[0][1] = lambda/2.;
		C[1][0] = lambda/2.; C[1][1] = -lambda/2.;
		initialize(n_states,&C);
		return *this;
	}
	int fast_mutate(const int state, const double time) {
		return (ran->bernoulliTF(0.5+0.5*exp(-2.0*C[0][1]*time))) ? state : !state;
	}
};

/*multinomial sampler*/
#endif // _MUTATION_H_
