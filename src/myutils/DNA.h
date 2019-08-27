/*  Copyright 2012 Daniel Wilson.
 *
 *  DNA.h
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
/*	DNA.h 23rd February 2005				*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _DNA_H_
#define _DNA_H_

#pragma warning(disable: 4786)

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "myutils.h"
#include <map>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace myutils;

class DNA {
public:
	vector<string> label;
	vector<string> sequence;
	int nseq;
	int lseq;
	vector<double> ntimes;
	bool coutput;
	map<char,int> baseToInt;		// converts TUCAG- to 112345
	map<int,char> intToBase;		// converts 012345 to NTCAG-

protected:
	vector<int> _uniqueHaps;
	vector<int> _sites;
	LowerTriangularMatrix<int> __B;
	vector<int> _M;
	vector<double> _F;
	vector<double> _four;
	LowerTriangularMatrix< vector<double> > _G;
	LowerTriangularMatrix<double> _A;
	LowerTriangularMatrix<double> _B;
	LowerTriangularMatrix<double> _CC;
	Matrix<double> _D;

public:
	DNA() {
		coutput = false;
		init();
	}
	DNA(const char* filename) {
		coutput = false;
		readFASTA_1pass(filename);
		init();
	}
	DNA& init() {
		baseToInt['T'] = 1;
		baseToInt['U'] = baseToInt['T'];
		baseToInt['C'] = 2;
		baseToInt['A'] = 3;
		baseToInt['G'] = 4;
		baseToInt['-'] = 5;
		intToBase[0] = 'N';
		intToBase[1] = 'T';
		intToBase[2] = 'C';
		intToBase[3] = 'A';
		intToBase[4] = 'G';
		intToBase[5] = '-';
		return *this;
	}
	/*DNA& readFASTA(const char* filename) {
		ifstream in1(filename);
		if(!in1.is_open()) {
			string errmsg = "DNA::readFASTA(): File ";
			errmsg += string(filename);
			errmsg += " not found";
			error(errmsg.c_str());
		}

		int str;
		nseq = 0;
		while(!in1.eof()) {
			str = in1.get();
			if((char)str=='>') {
				++nseq;
			}
		}
		in1.close();
		if(coutput) cout << "Read in " << nseq << " sequence" << endl;
		if(nseq==0) {
			lseq = 0;
			return *this;
		}

		ifstream in2(filename);
		if(!in2.is_open())error("File not found second time");

		lseq = 0;
		string junk;
		while(!in2.eof()) {
			str = in2.get();
			if((char)str=='>') {
				getline(in2,junk);
				if (!junk.empty()&&*junk.rbegin()=='\r') junk.erase(junk.length()-1,1);
				while(!in2.eof()) {
					str = in2.get();
					if((char)str=='>') break;
					if(str!=-1 && (char)str!='\n' && (char)str!='\r')
						++lseq;
				}
				if(coutput) cout << "Sequences are " << lseq << " long" << endl;
				break;
			}
		}
		in2.close();

		string blank(lseq,' ');
		sequence.resize(nseq,blank);
		label.resize(nseq);
		ntimes.resize(nseq,0.0);

		ifstream in3(filename);
		if(!in3.is_open())error("File not found third time");
		int NSEQ = 0; int LSEQ = 0;
		while(true) {
			str = in3.get();
			if(in3.eof()) error("Cannot find sequences!");
			if((char)str=='>') {
				getline(in3,label[NSEQ]);
				if (!label[NSEQ].empty()&&*label[NSEQ].rbegin()=='\r') label[NSEQ].erase(label[NSEQ].length()-1,1);
				break;
			}
		}
		while(true) {
			str = in3.get();
			if(in3.eof()) break;
			if(LSEQ<lseq)
				sequence[NSEQ][LSEQ] = (char)str;
			if(str!=-1 && (char)str!='\n' && (char)str!='\r')
				++LSEQ;
			if((char)str=='>') {
				++NSEQ;
				getline(in3,label[NSEQ]);
				if (!label[NSEQ].empty()&&*label[NSEQ].rbegin()=='\r') label[NSEQ].erase(label[NSEQ].length()-1,1);
				LSEQ=0;
			}
		}
		in3.close();

		if(coutput) for(NSEQ=0;NSEQ<nseq;NSEQ++) {
			cout << label[NSEQ] << endl;
			cout << sequence[NSEQ] << endl;
		}
		
		return *this;
	}*/
	DNA& readFASTA_1pass(const char* filename) {
		ifstream in1(filename);
		if(!in1.is_open()) {
			string errmsg = "DNA::readFASTA_1pass(): File ";
			errmsg += string(filename);
			errmsg += " not found";
			error(errmsg.c_str());
		}
		
		nseq = 0;
		lseq = -1;
		string s;
		getline(in1,s);
		if (!s.empty()&&*s.rbegin()=='\r') s.erase(s.length()-1,1);
		s.erase(remove(s.begin(),s.end(),' '),s.end());
		if(s.length()>0 && s[0]!='>') {
			string errmsg = "DNA::readFASTA_1pass(): File ";
			errmsg += string(filename);
			errmsg += " did not begin with '>'";
			error(errmsg.c_str());
		}
		label.push_back(s.substr(1));
		string newseq = "";
		while(!in1.eof()) {
			getline(in1,s);
			if (!s.empty()&&*s.rbegin()=='\r') s.erase(s.length()-1,1);
			s.erase(remove(s.begin(),s.end(),' '),s.end());
			if(s.length()>0 && s[0]=='>') {
				if(lseq==-1) lseq = newseq.length();
				if(newseq.length()!=lseq) {
					string errmsg = "DNA::readFASTA_1pass(): File ";
					errmsg += string(filename);
					errmsg += " sequences had different lengths";
					error(errmsg.c_str());
				}
				sequence.push_back(newseq);
				newseq = "";
				++nseq;
				label.push_back(s.substr(1));
			} else {
				newseq += s;
			}
		}
		if(lseq==-1) lseq = newseq.length();
		if(newseq.length()!=lseq) {
			string errmsg = "DNA::readFASTA_1pass(): File ";
			errmsg += string(filename);
			errmsg += " sequences had different lengths";
			error(errmsg.c_str());
		}
		sequence.push_back(newseq);
		newseq = "";
		++nseq;
		ntimes = vector<double>(nseq,0.0);
		in1.close();
		
		if(sequence.size()!=label.size()) {
			string errmsg = "DNA::readFASTA_1pass(): File ";
			errmsg += string(filename);
			errmsg += " different number of sequences and labels";
			error(errmsg.c_str());
		}
		
		if(coutput) for(int NSEQ=0;NSEQ<nseq;NSEQ++) {
			cout << label[NSEQ] << endl;
			cout << sequence[NSEQ] << endl;
		}
		
		return *this;
	}
	DNA& writeFASTA(const char* filename) {
		ofstream fout(filename);
		int n,pos;
		for(n=0;n<nseq;n++)
		{
			fout << ">" << label[n] << endl;
			for(pos=0;pos<lseq;pos++)
				fout << sequence[n][pos];
			fout << endl;
		}
		fout.close();
		return *this;
	}
	DNA& writeFASTA(vector<char> &code, const char* filename) {
		ofstream fout(filename);
		int n,pos;
		for(n=0;n<nseq;n++)
		{
			fout << ">" << label[n] << endl;
			for(pos=0;pos<lseq;pos++)
				fout << code[sequence[n][pos]];
			fout << endl;
		}
		fout.close();
		return *this;
	}
	DNA& writeFASTA(vector<string> &code, const char* filename) {
		ofstream fout(filename);
		int n,pos;
		for(n=0;n<nseq;n++)
		{
			fout << ">" << label[n] << endl;
			for(pos=0;pos<lseq;pos++)
				fout << code[sequence[n][pos]];
			fout << endl;
		}
		fout.close();
		return *this;
	}
	/*Subscript operator*/
	inline string& operator[](int pos) {
		return sequence[pos];
	}
	/*Resize wipes current entries*/
	DNA& resize(const int NSEQ, const int LSEQ) {
		if(NSEQ<0) error("DNA: NSEQ must be non-negative");
		if(LSEQ<0) error("DNA: LSEQ must be non-negative");
		nseq = NSEQ;
		lseq = LSEQ;

		string blank(lseq,' ');
		string empty = "";
		
		sequence.resize(nseq,blank);
		label.resize(nseq,empty);
		ntimes.resize(nseq);
		int i;
		for(i=0;i<nseq;i++) {
			sequence[i] = blank;
			label[i] = empty;
			ntimes[i] = 0.0;
		}
		return *this;
	}
	DNA& resize(const int NSEQ, const string& str) {
		if(NSEQ<0) error("DNA: NSEQ must be non-negative");
		nseq = NSEQ;
		sequence.resize(nseq,str);
		lseq = (int)sequence[0].size();
		string empty = "";
		label.resize(nseq,empty);
		ntimes.resize(nseq,0.0);
		return *this;
	}
	DNA& clear() {
		string blank(lseq,' ');
		sequence.resize(nseq,blank);
		string empty = "";
		label.resize(nseq,empty);
		ntimes.resize(nseq,0.0);
		return *this;
	}
public:
	/* Number of segregating sites */
	double S() {
		double result = 0.0;
		if(nseq==0) return 0.0;

		int i,j;
		for(j=0;j<lseq;j++) {
			char hap = sequence[0][j];
			for(i=1;i<nseq;i++)
				if(sequence[i][j]!=hap) {
					++result;
					break;
				}
		}

		return result;
	}

	/* Number of unique haplotypes */
	double H() {
		int result = 1;

		if(nseq==0) return 0.0;
		vector<int> uniqueHaps(nseq,-1);
		uniqueHaps[0] = 0;

		int i,ii,j;
		bool unique;
		for(i=1;i<nseq;i++) {
			unique = true;
			for(ii=0;ii<result;ii++) {
				for(j=0;j<lseq;j++)
					if(sequence[i][j]!=sequence[uniqueHaps[ii]][j]) break;
				if(j==lseq) unique = false;
			}
			if(unique==true) {
				uniqueHaps[result] = i;
				++result;
			}
		}

		return (double)result;
	}

	/* Average number of pairwise differences */
	double pi() {
		double result = 0.0;

		int i,j,k;
		for(i=0;i<nseq;i++)
			for(j=0;j<i;j++)
				for(k=0;k<lseq;k++)
					result += (sequence[i][k]==sequence[j][k]) ? 0.0 : 1.0;
		result *= 2.0/(double)(nseq)/(double)(nseq-1);

		return result;
	}

	/* Variance in number of pairwise differences */
	double Varpi() {
		double E,EE,pi;
		int i,j,k;

		E = EE = 0.0;
		for(i=0;i<nseq;i++)
			for(j=0;j<i;j++) {
				pi = 0.0;
				for(k=0;k<lseq;k++)
					pi += (sequence[i][k]==sequence[j][k]) ? 0.0 : 1.0;
				E += pi;
				EE += pi*pi;
			}
		E *= 2.0/(double)(nseq)/(double)(nseq-1);
		EE *= 2.0/(double)(nseq)/(double)(nseq-1);

		double result = EE - E*E;
		return result;
	}

	double Tajima() {
		double D = 0.0;
		int i,j,k,n,L;
		n = nseq;
		L = lseq;
		double a1,a2,b1,b2,c1,c2,e1,e2,khat,S;
		bool segregating;
		khat = S = 0.0;
		for(k=0;k<L;k++) {
			segregating = false;
			for(i=0;i<n;i++)
				for(j=0;j<i;j++)
					if(sequence[i][k]!=sequence[j][k]) {
						++khat;
						segregating = true;
					}
			if(segregating) ++S;
		}
		if(S==0) return 0.0;
		khat /= (double)(n*(n-1)/2);
		a1 = a2 = 0.0;
		for(i=1;i<=n-1;i++) {
			a1 += 1./(double)i;
			a2 += 1./(double)(i*i);
		}
		b1 = (double)(n+1)/(double)(3*(n-1));
		b2 = (double)(2*(n*n+n+3))/(double)(9*n*(n-1));
		c1 = b1 - 1./a1;
		c2 = b2 - (double)(n+2)/a1/(double)(n) + a2/a1/a1;
		e1 = c1/a1;
		e2 = c2/(a1*a1+a2);
		D = (khat - S/a1)/sqrt(e1*S+e2*S*(S-1.));
		return D;
	}
	/* This function counts the average number of pairwise differences, where the matrix
	   diff defines those differences using 0 = different or 1 = identical.
	   If diff is the identity matrix then this function is equivalent to double pi(). */
	double pi(Matrix<int> &diff, map<char,int> &chmap) {
		double result = 0.0;

		int i,j,k;
		for(i=0;i<nseq;i++)
			for(j=0;j<i;j++)
				for(k=0;k<lseq;k++)
					result += (diff[chmap[(char)sequence[i][k]]][chmap[(char)sequence[j][k]]]==0);
		result *= 2.0/(double)(nseq)/(double)(nseq);

		return result;
	}
	/* Hudson and Kaplan's Rm, the minimum # recombinations. See Myers and Griffiths(2003)*/
	double Rm() {
		if(nseq==0) return 0.0;
		if(lseq==0) return 0.0;

		/* Determine which sites are biallelic segregating */
		_sites = vector<int>(lseq,0);
		int i,j,k;
		int S = 0;
		char hap0,hap1;
		bool segregating;
		for(j=0;j<lseq;j++) {
			segregating = false;
			hap0 = sequence[0][j];
			for(i=1;i<nseq;i++) {
				if(!segregating && sequence[i][j]!=hap0) {
					segregating = true;
					hap1 = sequence[i][j];
				}
				else if(segregating && sequence[i][j]!=hap0 && sequence[i][j]!=hap1) {
					segregating = false;	// define segregating only for biallelic sites
					break;
				}
			}
			if(segregating) {
				_sites[S] = j;
				++S;
			}
		}
		if(S<2) return 0.0;

		/* Calculate the compatibility matrix */
		__B = LowerTriangularMatrix<int>(S,0);	// so j>=k always
		// __B[j][k] = 0 for compatible, 1 for incompatible
		bool comb[3];
		for(j=0;j<S;j++)
			for(k=0;k<j;k++)
			{
				hap0 = sequence[0][_sites[j]];
				hap1 = sequence[0][_sites[k]];
				comb[0] = false;				// hap0  hap1'
				comb[1] = false;				// hap0' hap1
				comb[2] = false;				// hap0' hap1'
				for(i=1;i<nseq;i++) {
					if(sequence[i][_sites[j]]==hap0 && sequence[i][_sites[k]]!=hap1) comb[0] = true;
					if(sequence[i][_sites[j]]!=hap0 && sequence[i][_sites[k]]==hap1) comb[1] = true;
					if(sequence[i][_sites[j]]!=hap0 && sequence[i][_sites[k]]!=hap1) comb[2] = true;
					if(comb[0] && comb[1] && comb[2]) break;
				}
				__B[j][k] = (comb[0] && comb[1] && comb[2]) ? 1 : 0;			
			}

		/* Calculate the dynamic programming partition matrix */
		_M = vector<int>(S,0);
//		int maxM = 0;
		_M[S-1] = 0;
		_M[S-2] = __B[S-1][S-2];
		for(i=S-3;i>=0;i--) {
			_M[i] = __B[i+1][i] + _M[i+1];
			for(k=i+2;k<S;k++) if(__B[k][i]+_M[k]>_M[i]) _M[i] = __B[k][i]+_M[k];
		}

		return (double)_M[0];
	}
	void RecCorrelations(double* result) {
		RecCorrelations(result,true);
	}
	void RecCovariances(double* result) {
		RecCorrelations(result,false);
	}
	void RecCorrelations(double* result, bool normalize) {
		result[0] = result[1] = result[2] = 0.0;

		if(nseq==0) return;
		if(lseq==0) return;

		/* Determine which sites are biallelic segregating */
		_sites = vector<int>(lseq,0);
		int i,j,k;
		int S = 0;
		char hap0,hap1;
		bool segregating;
		for(j=0;j<lseq;j++) {
			segregating = false;
			hap0 = sequence[0][j];
			for(i=1;i<nseq;i++) {
				if(!segregating && sequence[i][j]!=hap0) {
					segregating = true;
					hap1 = sequence[i][j];
				}
				else if(segregating && sequence[i][j]!=hap0 && sequence[i][j]!=hap1) {
					segregating = false;	// define segregating only for biallelic sites
					break;
				}
			}
			if(segregating) {
				_sites[S] = j;
				++S;
			}
		}
		if(S<3) return;
		
		/* Calculate frequency statistics */
		_F = vector<double>(S,1.0);							/* _F is the marginal frequency of hap0 at site j */
		for(j=0;j<S;j++) {
			hap0 = sequence[0][_sites[j]];
			for(i=1;i<nseq;i++)
				if(sequence[i][_sites[j]]==hap0) _F[j]++;
			_F[j] /= (double)nseq;
		}

		_four = vector<double>(4,0.0);							/* _G[j][k] is the frequency of AB (_G[j][k][0]),	*/
		_G = LowerTriangularMatrix< vector<double> >(S,_four);	/* Ab (1), aB (2), ab (3) for sites j and k		*/
		for(j=0;j<S;j++)
		  for(k=0;k<j;k++) {
			  hap0 = sequence[0][_sites[j]];
			  hap1 = sequence[0][_sites[k]];
			  for(i=0;i<nseq;i++) {
				  if(sequence[i][_sites[j]]==hap0 && sequence[i][_sites[k]]==hap1) ++_G[j][k][0];
				  else if(sequence[i][_sites[j]]==hap0 && sequence[i][_sites[k]]!=hap1) ++_G[j][k][1];
				  else if(sequence[i][_sites[j]]!=hap0 && sequence[i][_sites[k]]==hap1) ++_G[j][k][2];
				  else if(sequence[i][_sites[j]]!=hap0 && sequence[i][_sites[k]]!=hap1) ++_G[j][k][3];
				  else warning("Unexpected choice");
			  }
			  for(i=0;i<4;i++) _G[j][k][i] /= (double)nseq;
		  }
  
		/* Calculate LD statistics for pairs of sites */
		_A = LowerTriangularMatrix<double>(S,0.0);			//	rsq
		_B = LowerTriangularMatrix<double>(S,0.0);			//	Dprime
		_CC = LowerTriangularMatrix<double>(S,0.0);			//	G4
		_D = Matrix<double>(S,S,0.0);

		double temp;
		for(i=0;i<S;i++) {
			for(j=0;j<i;j++) {
				temp = _G[i][j][0] - _F[i]*_F[j];
				_A[i][j] = pow(temp,2.0)/(_F[i]*(1.-_F[i])*_F[j]*(1.-_F[j]));
				_B[i][j] = (temp < 0.0) ? -temp/MIN(_F[i]*_F[j],(1.-_F[i])*(1.-_F[j])) : temp/MIN(_F[i]*(1.-_F[j]),(1.-_F[i])*_F[j]);
				_CC[i][j] = (_G[i][j][0]>0.0 && _G[i][j][1]>0.0 && _G[i][j][2]>0.0 && _G[i][j][3]>0.0) ? 1.0 : 0.0;
				_D[i][j] = _D[j][i] = _sites[i] - _sites[j];
			}
		}

		double  E[4] = {0.0,0.0,0.0,0.0};
		double EE[4] = {0.0,0.0,0.0,0.0};
		double ED[3] = {0.0,0.0,0.0};
		int ctr;
		for(i=0,ctr=0;i<S;i++)
			for(j=0;j<i;j++,ctr++) {
				E[0] += _A[i][j]; E[1] += _B[i][j]; E[2] += _CC[i][j]; E[3] += _D[i][j];
				EE[0] += _A[i][j]*_A[i][j]; EE[1] += _B[i][j]*_B[i][j]; EE[2] += _CC[i][j]*_CC[i][j]; EE[3] += _D[i][j]*_D[i][j];
				ED[0] += _A[i][j]*_D[i][j]; ED[1] += _B[i][j]*_D[i][j]; ED[2] += _CC[i][j]*_D[i][j];
			}

		if(normalize)	// Calculate correlation
			for(k=0;k<3;k++)
				result[k] = (ED[k]-E[k]*E[3]/(double)(ctr))/sqrt((EE[k]-E[k]*E[k]/(double)(ctr)))/sqrt((EE[3]-E[3]*E[3]/(double)(ctr)));
		else			// Calculate covariance
			for(k=0;k<3;k++)
				result[k] = (ED[k]-E[k]*E[3]/(double)(ctr))/(double)(ctr);

		return;
	}
	/* Convert the DNA sequence into an amino acid sequence */
	DNA& translate(const int offset,vector<string> &polypeptide) {
		if(offset<0) error("DNA::transcribe(): cannot have negative offset");
		if((lseq-offset)%3!=0) error("DNA::transcribe(): DNA length minus offset isn't a multiple of 3");
		const int tlen = (lseq-offset)/3;
		string blank(" ",tlen);
		polypeptide = vector<string>(nseq,blank);
		int i,j,ctr;
		for(i=offset,ctr=0;i<lseq;i+=3,ctr+=1) {
			for(j=0;j<nseq;j++)	{
				string triplet = sequence[j].substr(i,3);
				polypeptide[j][ctr] = codonToPeptide(tripletToCodon(triplet));
			}
		}
		return *this;
	}
	/* Convert the DNA sequence into a sequence of codons (characters 0-63, 64 for indel, -1 for unknown) */
	DNA& tocodon(const int offset,vector<string> &codonsequence) {
		if(offset<0) error("DNA::tocodon(): cannot have negative offset");
		if((lseq-offset)%3!=0) error("DNA::tocodon(): DNA length minus offset isn't a multiple of 3");
		const int tlen = (lseq-offset)/3;
		string blank(" ",tlen);
		codonsequence = vector<string>(nseq,blank);
		int i,j,ctr;
		for(i=offset,ctr=0;i<lseq;i+=3,ctr+=1) {
			for(j=0;j<nseq;j++)	{
				string triplet = sequence[j].substr(i,3);
				codonsequence[j][ctr] = (char)tripletToCodon(triplet);
			}
		}
		return *this;
	}
	/* Convert the DNA sequence into a sequence of nucleotide numbers (TCAG=0123, 4 for indel, -1 for unknown) */
	DNA& tonucleotide(const int offset,Matrix<int> &ntsequence) {
		if(offset<0) error("DNA::tonucleotide(): cannot have negative offset");
		if(offset>=lseq) error("DNA::tonucleotide(): cannot offset the whole sequence");
		const int tlen = lseq-offset;
		ntsequence = Matrix<int>(nseq,tlen);
		int i,j,ctr;
		for(i=offset,ctr=0;i<lseq;i++,ctr++) {
			for(j=0;j<nseq;j++)	{
				ntsequence[j][ctr] = baseToInt[sequence[j][i]]-1;
			}
		}
		return *this;
	}
	/* Convert the DNA sequence into a sequence of codon numbers (0-63, 64 for indel, -1 for unknown) */
	DNA& tocodon(const int offset,Matrix<int> &codonsequence) {
		if(offset<0) error("DNA::tocodon(): cannot have negative offset");
		if((lseq-offset)%3!=0) error("DNA::tocodon(): DNA length minus offset isn't a multiple of 3");
		const int tlen = (lseq-offset)/3;
		codonsequence = Matrix<int>(nseq,tlen);
		int i,j,ctr;
		for(i=offset,ctr=0;i<lseq;i+=3,ctr+=1) {
			for(j=0;j<nseq;j++)	{
				string triplet = sequence[j].substr(i,3);
				codonsequence[j][ctr] = tripletToCodon(triplet);
			}
		}
		return *this;
	}
	/* Convert the DNA sequence into a sequence of codon numbers (0-60, 61 for indel, -1 for unknown) */
	DNA& tocodon61(const int offset,Matrix<int> &codonsequence) {
		if(offset<0) error("DNA::tocodon(): cannot have negative offset");
		if((lseq-offset)%3!=0) error("DNA::tocodon(): DNA length minus offset isn't a multiple of 3");
		const int tlen = (lseq-offset)/3;
		codonsequence = Matrix<int>(nseq,tlen);
		int i,j,ctr;
		for(i=offset,ctr=0;i<lseq;i+=3,ctr+=1) {
			for(j=0;j<nseq;j++)	{
				string triplet = sequence[j].substr(i,3);
				codonsequence[j][ctr] = tripletToCodon61(triplet);
			}
		}
		return *this;
	}
	/* Convert the DNA sequence into a sequence of codon numbers (0-60, 61 for indel, -1 for unknown) */
	DNA& tocodon61_noerror(const int offset,Matrix<int> &codonsequence) {
		if(offset<0) error("DNA::tocodon(): cannot have negative offset");
		if((lseq-offset)%3!=0) error("DNA::tocodon(): DNA length minus offset isn't a multiple of 3");
		const int tlen = (lseq-offset)/3;
		codonsequence = Matrix<int>(nseq,tlen);
		int i,j,ctr;
		for(i=offset,ctr=0;i<lseq;i+=3,ctr+=1) {
			for(j=0;j<nseq;j++)	{
				string triplet = sequence[j].substr(i,3);
				codonsequence[j][ctr] = tripletToCodon61_noerror(triplet);
			}
		}
		return *this;
	}
	/* Convert the DNA sequence into a sequence of codon numbers (0-60, 61 for indel, -1 for unknown) */
	DNA& tocodon61_warning(const int offset,Matrix<int> &codonsequence) {
		if(offset<0) error("DNA::tocodon(): cannot have negative offset");
		if((lseq-offset)%3!=0) error("DNA::tocodon(): DNA length minus offset isn't a multiple of 3");
		const int tlen = (lseq-offset)/3;
		codonsequence = Matrix<int>(nseq,tlen);
		int i,j,ctr;
		for(i=offset,ctr=0;i<lseq;i+=3,ctr+=1) {
			for(j=0;j<nseq;j++)	{
				string triplet = sequence[j].substr(i,3);
				codonsequence[j][ctr] = tripletToCodon61_noerror(triplet);
				if(codonsequence[j][ctr]==-2) {
					stringstream wrn_txt;
					wrn_txt << "Found a STOP codon in sequence " << j << " position " << i << " (codon " << ctr << ")";
					warning(wrn_txt.str().c_str());
				}
			}
		}
		return *this;
	}
	/* Returns 0-63 for codons, 64 for indels and -1 for unknown */
	int tripletToCodon(string &tri) {
		const int a = baseToInt[tri[0]];
		const int b = baseToInt[tri[1]];
		const int c = baseToInt[tri[2]];
		bool indel = false;
		if(a==5) indel = true;
		if(b==5) indel = true;
		if(c==5) indel = true;
		if(indel==true) {
			if(a==5 && b==5 && c==5) return 64;
			else return -1;
		}
		/* return a value from 0 to 63 */
		return (a-1)*16 + (b-1)*4 + c - 1;
	}
	/* Returns 0-60 for non-STOP codons, 61 for indels and -1 for unknown */
	int tripletToCodon61(string &tri) {
		const int a = baseToInt[tri[0]];
		const int b = baseToInt[tri[1]];
		const int c = baseToInt[tri[2]];
		bool indel = false;
		if(a==5) indel = true;
		if(b==5) indel = true;
		if(c==5) indel = true;
		if(indel==true) {
			if(a==5 && b==5 && c==5) return 61;
			else return -1;
		}
		/* return a value from 0 to 63 */
		int ret = (a-1)*16 + (b-1)*4 + c - 1;
		/* remove STOP codons so value ranges from 0 to 60 */
		if(ret==10 || ret==11 || ret==14) error("DNA::tripletToCodon61(): STOP codon found");
		if(ret>=14) --ret; /* (shouldn't ever be equal to because of previous line) */
		if(ret>=11) --ret;
		if(ret>=10) --ret;
		return ret;
	}
	/* Returns 0-60 for non-STOP codons, 61 for indels and -1 for unknown */
	int tripletToCodon61_noerror(string &tri) {
		const int a = baseToInt[tri[0]];
		const int b = baseToInt[tri[1]];
		const int c = baseToInt[tri[2]];
		bool indel = false;
		if(a==5) indel = true;
		if(b==5) indel = true;
		if(c==5) indel = true;
		if(indel==true) {
			if(a==5 && b==5 && c==5) return 64;
			else return -1;
		}
		/* return a value from 0 to 63 */
		int ret = (a-1)*16 + (b-1)*4 + c - 1;
		/* remove STOP codons so value ranges from 0 to 60 */
		if(ret==10 || ret==11 || ret==14) return -2; // WARNING value instead of ERROR
		if(ret>=14) --ret; /* (shouldn't ever be equal to because of previous line) */
		if(ret>=11) --ret;
		if(ret>=10) --ret;
		return ret;
	}
	char codonToPeptide(const int codon) {
		switch(codon) {
		case 0:		return 'F';
		case 1:		return 'F';
		case 2:		return 'L';
		case 3:		return 'L';
		case 4:		return 'S';
		case 5:		return 'S';
		case 6:		return 'S';
		case 7:		return 'S';
		case 8:		return 'Y';
		case 9:		return 'Y';
		case 10:		return 'X';
		case 11:		return 'X';
		case 12:		return 'C';
		case 13:		return 'C';
		case 14:		return 'X';
		case 15:		return 'W';
		case 16:		return 'L';
		case 17:		return 'L';
		case 18:		return 'L';
		case 19:		return 'L';
		case 20:		return 'P';
		case 21:		return 'P';
		case 22:		return 'P';
		case 23:		return 'P';
		case 24:		return 'H';
		case 25:		return 'H';
		case 26:		return 'Q';
		case 27:		return 'Q';
		case 28:		return 'R';
		case 29:		return 'R';
		case 30:		return 'R';
		case 31:		return 'R';
		case 32:		return 'I';
		case 33:		return 'I';
		case 34:		return 'I';
		case 35:		return 'M';
		case 36:		return 'T';
		case 37:		return 'T';
		case 38:		return 'T';
		case 39:		return 'T';
		case 40:		return 'N';
		case 41:		return 'N';
		case 42:		return 'K';
		case 43:		return 'K';
		case 44:		return 'S';
		case 45:		return 'S';
		case 46:		return 'R';
		case 47:		return 'R';
		case 48:		return 'V';
		case 49:		return 'V';
		case 50:		return 'V';
		case 51:		return 'V';
		case 52:		return 'A';
		case 53:		return 'A';
		case 54:		return 'A';
		case 55:		return 'A';
		case 56:		return 'D';
		case 57:		return 'D';
		case 58:		return 'E';
		case 59:		return 'E';
		case 60:		return 'G';
		case 61:		return 'G';
		case 62:		return 'G';
		case 63:		return 'G';
		}
		return '?';
	}
};



#endif // _DNA_H_
