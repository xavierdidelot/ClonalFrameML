/*  Copyright 2012 Daniel Wilson.
 *
 *  lotri_matrix.h
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
/*	lotri_matrix.h 23rd February 2005		*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _LOWER_TRIANGULAR_MATRIX_H_
#define _LOWER_TRIANGULAR_MATRIX_H_

#include <stdlib.h>
#include <stdio.h>

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
class LowerTriangularMatrix
{
public:
	/*Preserve public access for back-compatibility*/
	T *array;
	T **element;

protected:
	int _n;		/* dimension of the lower triangular square matrix	*/
	int _size;	/* number of elements of the matrix					*/
//	int protected_ncols;
	int initialized;

public:
	/*Default constructor*/	LowerTriangularMatrix()
	{
		initialized=0;
		initialize(0);
	}
	/*Constructor*/			LowerTriangularMatrix(int n)
	{
		initialize(n);
	}
	/*Constructor*/			LowerTriangularMatrix(int n, T value)
	{
		initialize(n);
		int i,j;
		for(i=0;i<n;i++)
			for(j=0;j<=i;j++)
				element[i][j]=value;
	}
	/*Destructor*/			~LowerTriangularMatrix()
	{
		delete[] array;
		delete[] element;
	}
	LowerTriangularMatrix<T>& initialize(int n)
	{
		int i;
		int size = n*(n+1)/2;
		array = new T[size];
		if (!array) error("array allocation failure in LowerTriangularMatrix::initialize()");

		element = new T*[n];
		if (!element) error("element allocation failure in LowerTriangularMatrix::initialize()");
		for(i=0;i<n;i++) element[i] = &(array[i*(i+1)/2+0]);

		_n = n;
		_size = size;
		initialized=1;
		return *this;
	}
	/*All current data is lost when the LowerTriangularMatrix is resized*/
	LowerTriangularMatrix<T>& resize(int n)
	{
		int i;
		int size = n*(n+1)/2;
		if (!initialized) return initialize(n);
		if(n==_n)return *this;
		
		delete[] array;
		delete[] element;

		array = new T[size];
		if (!array) error("array allocation failure in LowerTriangularMatrix::resize()");
		
		element = new T*[n];
		if (!element) error("element allocation failure in LowerTriangularMatrix::resize()");
		for(i=0;i<n;i++) element[i] = &(array[i*(i+1)/2+0]);

		_n = n;
		_size = size;
		return *this;
	}
	int n(){return _n;}
	int size(){return _size;}
	int n() const {return _n;}
	int size() const {return _size;}
	void error(const char* error_text)
	{
		printf("Run-time error in LowerTriangularMatrix::");
		printf("%s\n", error_text);
		printf("Exiting to system...\n");
		exit(13);
	}
	/*Copy constructor*/	LowerTriangularMatrix(const LowerTriangularMatrix<T> &mat)
	/*	Copy constructor for the following cases:
			LowerTriangularMatrix mat2(mat);
			LowerTriangularMatrix mat2=mat;
		and when LowerTriangularMatrix is returned from a function	*/
	{
		initialize(mat._n);
		int i;
		for(i=0;i<_size;i++)
			array[i] = mat.array[i];
	}
	/*Assignment operator*/	LowerTriangularMatrix<T>& operator=(const LowerTriangularMatrix<T>& mat)
	{
		//if(this==mat)return *this;
		resize(mat._n);
		int i;
		for(i=0;i<_size;i++)
			array[i] = mat.array[i];
		return *this;
	}
	/*Subscript operator*/inline T* operator[](int pos){return element[pos];};
	inline T& safe(int i, int j) {
		return (j<=i) ? element[i][j] : element[j][i];
	}
};
};

#endif // _LOWER_TRIANGULAR_MATRIX_H_