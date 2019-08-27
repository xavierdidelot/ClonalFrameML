/*  Copyright 2012 Daniel Wilson.
 *
 *  matrix.h
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
/*	matrix.h 23rd February 2005				*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdlib.h>
#include <stdio.h>
#include "vector.h"
#include "utils.h"

/****************************************************************/
/*						myutils::Matrix							*/
/*																*/
/*	Matrix is a C++ style container whose memory storage is		*/
/*	designed so that elements can easily be viewed at debug		*/
/*	time in MSVC++ and to be compatible with some C code in		*/
/*	which matrices are stored as one-dimensional arrays, where	*/
/*	element (i,j) would be accessed as M[i*n+j].				*/
/*																*/
/*	Element (i,j) can be accessed in one of three ways:			*/
/*		M[i][j]				clearest syntax						*/
/*		M.element[i][j]		useful for viewing during debug		*/
/*		M.array[i*n+j]		compatible with C arrays			*/
/*																*/
/****************************************************************/

namespace myutils
{
template <typename T>
class safeArray {
public:
	T *element;
	int lo,hi;
public:
	safeArray(T *set_element, const int set_lo, const int set_hi) : element(set_element), lo(set_lo), hi(set_hi) {};
	inline T& operator[](int pos){
		if(pos<lo) error("safeArray::operator[](int pos): pos<min");
		if(pos>=hi) error("safeArray::operator[](int pos): pos>=max");
		return element[pos];
	}
	inline const T& operator[](int pos) const {
		if(pos<lo) error("safeArray::operator[](int pos): pos<min");
		if(pos>=hi) error("safeArray::operator[](int pos): pos>=max");
		return element[pos];
	}
};

template <typename T>
class Matrix
{
public:
	/*Preserve public access for back-compatibility*/
	T *array;
	T **element;

protected:
	unsigned long int protected_nrows;
	unsigned long int protected_ncols;
	int initialized;

public:
	/*Default constructor*/	Matrix()
	{
		initialized=0;
		initialize(0,0);
	}
	/*Constructor*/			Matrix(int nrows, int ncols)
	{
		initialize(nrows,ncols);
	}
	/*Constructor*/			Matrix(int nrows, int ncols, T value)
	{
		initialize(nrows,ncols);
		unsigned long int i,j;
		for(i=0;i<nrows;i++)
			for(j=0;j<ncols;j++)
				element[i][j]=value;
	}
	/*Destructor*/			~Matrix()
	{
		delete[] array;
		delete[] element;
	}
	Matrix<T>& initialize(int nrows, int ncols)
	{
		unsigned long int i;
		const unsigned long int newsize = (unsigned long int)(nrows)*(unsigned long int)(ncols);
		array = new T[newsize];
		if (!array) error("array allocation failure in Matrix::initialize()");

		element = new T*[(unsigned long int)nrows];
		if (!element) error("element allocation failure in Matrix::initialize()");
		for(i=0;i<nrows;i++) element[i] = &(array[i*(unsigned long int)ncols+0]);

		protected_nrows=(unsigned long int)nrows;
		protected_ncols=(unsigned long int)ncols;
		initialized=1;
		return *this;
	}
	/*All current data is lost when the Matrix is resized*/
	Matrix<T>& resize(int nrows, int ncols)
	{
		unsigned long int i;
		if (!initialized) return initialize(nrows,ncols);
		if((nrows==protected_nrows)&&(ncols==protected_ncols))return *this;
		
		delete[] array;
		delete[] element;

		const unsigned long int newsize = (unsigned long int)(nrows)*(unsigned long int)(ncols);
		array = new T[newsize];
		if (!array) error("array allocation failure in Matrix::resize()");
		
		element = new T*[(unsigned long int)nrows];
		if (!element) error("element allocation failure in Matrix::resize()");
		for(i=0;i<nrows;i++) element[i] = &(array[i*(unsigned long int)ncols+0]);

		protected_nrows=(unsigned long int)nrows;
		protected_ncols=(unsigned long int)ncols;
		return *this;
	}
	int nrows(){return (int)protected_nrows;}
	int ncols(){return (int)protected_ncols;}
	int nrows() const {return (int)protected_nrows;}
	int ncols() const {return (int)protected_ncols;}
/*	void error(char* error_text)
	{
		printf("Run-time error in Matrix::");
		printf("%s%\n", error_text);
		printf("Exiting to system...\n");
		exit(13);
	}*/
	/*Copy constructor*/	Matrix(const Matrix<T> &mat)
	/*	Copy constructor for the following cases:
			Matrix mat2(mat);
			Matrix mat2=mat;
		and when Matrix is returned from a function	*/
	{
	  initialize((int)mat.protected_nrows,(int)mat.protected_ncols);
		int i;
		for(i=0;i<protected_nrows*protected_ncols;i++)
			array[i] = mat.array[i];
	}
	/*Assignment operator*/	Matrix<T>& operator=(const Matrix<T>& mat)
	{
		//if(this==mat)return *this;
		resize(mat.nrows(),mat.ncols());
		int i;
		for(i=0;i<protected_nrows*protected_ncols;i++)
			array[i] = mat.array[i];
		return *this;
	}
#ifdef _MYUTILS_DEBUG
	/*DEBUG Subscript operator*/inline safeArray< T > operator[](unsigned long int pos){
		if(pos<0) error("Matrix::operator[](int row): row<0");
		if(pos>=protected_nrows) error("Matrix::operator[](int row): row>=nrows()");
		//return element[pos];
		return safeArray< T >(element[pos],0,protected_ncols);
	};
	/*DEBUG Subscript operator*/inline const safeArray< T > operator[](unsigned long int pos) const {
		if(pos<0) error("Matrix::operator[](int row): row<0");
		if(pos>=protected_nrows) error("Matrix::operator[](int row): row>=nrows()");
		//return element[pos];
		return const safeArray< T >(element[pos],0,protected_ncols);
	};
#else
	/*Subscript operator*/inline T* operator[](unsigned long int pos){return element[pos];};
	/*Subscript operator*/inline const T* operator[](unsigned long int pos) const {return element[pos];};
#endif

	/*Matrix multiplication*/
	Matrix<T> operator*(const Matrix<T>& mat)
	{
		if(ncols()!=mat.nrows()) error("Matrix multiplication: matrices are not conformable");
		Matrix<T> result(nrows(),mat.ncols(),0.0);
		int i,j,k;
		for(i=0;i<nrows();i++)
			for(j=0;j<mat.ncols();j++)
				for(k=0;k<ncols();k++)
					result[i][j] += element[i][k] * mat.element[k][j];
		return result;
	}
	/* Multiply the operands and store the result in this matrix */
	Matrix<T>& multiply(const Matrix<T>& op1, const Matrix<T>& op2) {
		if(op1.ncols()!=op2.nrows()) error("Matrix multiplication: matrices are not conformable");
		resize(op1.nrows(),op2.ncols());
		int i,j,k;
		for(i=0;i<op1.nrows();i++)
			for(j=0;j<op2.ncols();j++)
				for(k=0;k<op1.ncols();k++)
					element[i][j] += op1.element[i][k] * op2.element[k][j];
		return *this;
	}
	/*apply a function to every element of the matrix*/
	Matrix<T> map(T (* f)(T))
	{
	  Matrix<T> result((int)protected_nrows,(int)protected_ncols);
		int i,j;
		for(i=0;i<(int)protected_nrows;i++)
		  for(j=0;j<(int)protected_ncols;j++)
				result[i][j] = f(element[i][j]);
		return result;
	}

	/* Numerical Recipes in C++ routine for inverting a square real matrix */
	Matrix<T> invert() {
		if(protected_nrows!=protected_ncols) error("Matrix inversion: must be a symmetric matrix");
		Matrix<T> a = *this;
		Matrix<T> b(protected_nrows,protected_ncols,0);
		int i;
		for(i=0;i<protected_nrows;i++) b[i][i] = 1;

		int icol,irow,j,k,l,ll;
		double big,dum,pivinv;

		int n=a.nrows();
		int m=b.ncols();
		myutils::Vector<int> indxc(n);
		myutils::Vector<int> indxr(n);
		myutils::Vector<int> ipiv(n);
		for(j=0;j<n;j++) ipiv[j]=0;
		for(i=0;i<n;i++)
		{
			big=0.0;
			for(j=0;j<n;j++)
				if(ipiv[j]!=1)
					for(k=0;k<n;k++)
					{
						if(ipiv[k]==0)
						{
							if(fabs(a[j][k])>=big)
							{
								big=fabs(a[j][k]);
								irow=j;
								icol=k;
							}
						}
					}
			++(ipiv[icol]);
			if(irow!=icol)
			{
				for(l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
				for(l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
			}
			indxr[i]=irow;
			indxc[i]=icol;
			if(a[icol][icol] == 0.0) error("Matrix inversion: Singular Matrix");
			pivinv=1.0/a[icol][icol];
			a[icol][icol]=1.0;
			for(l=0;l<n;l++) a[icol][l] *= pivinv;
			for(l=0;l<m;l++) b[icol][l] *= pivinv;
			for(ll=0;ll<n;ll++)
				if (ll != icol)
				{
					dum=a[ll][icol];
					a[ll][icol]=0.0;
					for(l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
					for(l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
				}
		}
		for(l=n-1;l>=0;l--)
		{
			if(indxr[l]!=indxc[l])
				for(k=0;k<n;k++)
					SWAP(a[k][indxr[l]],a[k][indxc[l]]);
		}
		return a;
	}
};

template <typename T>
inline Matrix<T> IdentityMatrix(const int n) {
	Matrix<T> m(n,n,(T)0);
	int i;
	for(i=0;i<n;i++) m[i][i] = (T)1;
	return m;
};

};		// namespace myutils

#endif	// _MATRIX_H_
