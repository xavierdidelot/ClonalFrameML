/*  Copyright 2012 Daniel Wilson.
 *
 *  MLST.h
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
#ifndef _MLST_H_
#define _MLST_H_

#pragma warning(disable: 4786)

#include <myerror.h>
#include <vector.h>
#include <matrix.h>
#include <DNA.h>

using namespace myutils;

class MLST {
 public:
  int n;                  // number of sequences
  int nloc;               // number of loci
  Vector<int> nhap;       // nhap[l] (l=0..nloc-1) gives the number of unique alleles at locus l
  Vector<DNA> allele;     // allele[l] (l=0..nloc-1) stores the DNA sequences of the nhap[l] unique alleles at locus l
  Matrix<int> count;      // count[l][i] (l=0..nloc-1,i=0..nhap[l]-1) is the count of unique allele i at locus l
  Matrix<int> haplotype;  // haplotype[i] (i=0..n-1) gives the allelic profile for sequence i, so that
                          // haplotype[i][l] (l=0..nloc-1) is allele number at locus l, so that the DNA sequence
                          // is accessed using allele[l][haplotype[i][l]]. However, a short-cut would be, rather than
                          // using MLST.allele[l][haplotype[i][l]], to use MLST.seq(i,l).
 public:
  string& seq(const int i, const int l) {
    return allele[l][haplotype[i][l]];
  }
  MLST() {};
  MLST(const int nloc_in, const char* filename[]) {
    nloc = nloc_in;
    Vector<DNA*> temp(nloc);
    int l;
    for(l=0;l<nloc;l++) temp[l] = new DNA(filename[l]);
    initialize(temp);
    for(l=0;l<nloc;l++) delete temp[l];
  }
  MLST(Vector<DNA*> &temp) {
    initialize(temp);
  }
  void initialize(Vector<DNA*> &temp) {
    nloc = temp.size();
    if(nloc<1) myutils::error("MLST::initialize(): must be at least one locus");
    int l;
    n = temp[0]->nseq;
    for(l=1;l<nloc;l++) if(temp[l]->nseq!=n) myutils::error("MLST(): all loci should have the same number of sequences");
    nhap.resize(nloc);
    allele.resize(nloc);
    haplotype = Matrix<int>(n,nloc,-1);
    count = Matrix<int>(nloc,n,0);
    Vector<int> convert(n);
    int i,j;
    for(l=0;l<nloc;l++) {
      nhap[l] = 0;
      for(i=0;i<n;i++)
	for(j=0;j<=i;j++)
	  if(temp[l]->sequence[i]==temp[l]->sequence[j]) {
	    ++count[l][j];
	    haplotype[i][l] = j;
	    break;
	  }
      int check_total = 0;
      for(i=0;i<n;i++) {
	nhap[l] += (count[l][i]>0) ? 1 : 0;
	check_total += count[l][i];
      }
      if(check_total!=n) myutils::error("MLST(): problem in counting haplotypes");
      allele[l].resize(nhap[l],temp[l]->lseq);
      int hap = 0;
      for(i=0;i<n;i++) {
	if(count[l][i]>0) {
	  allele[l][hap] = temp[l]->sequence[i];
	  count[l][hap] = count[l][i];
	  convert[i] = hap;
	  ++hap;
	}
      }
      if(hap!=nhap[l]) myutils::error("MLST(): hap and nhap disagree");
      for(;hap<n;hap++) count[l][hap] = 0;
      for(i=0;i<n;i++)
	haplotype[i][l] = convert[haplotype[i][l]];
	  /*if(nhap.element[l]>1) {
		double pi = allele[l].pi();
		double H = allele[l].H();
	  }*/
    }
	/*cout << "Allelic profiles of the " << n << " haplotypes" << endl;
	for(i=0;i<n;i++){
		cout << "Hap" << i << ":";
		for(j=0;j<nloc;j++) cout << " " << haplotype[i][j];
		cout << endl;
	}
	cout << endl;
	cout << "Frequency of the haplotypes at each locus" << endl;
	for(l=0;l<nloc;l++) {
		cout << "Loc" << l << ":";
		for(i=0;i<nhap[l];i++) cout << " " << count[l][i];
		cout << endl;
	}
	cout << endl;*/

  }
};

#endif//_MLST_H_
