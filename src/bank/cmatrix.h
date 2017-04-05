/*  Copyright 2012 Daniel Wilson.
 *
 *  cmatrix.h
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
/*	cmatrix.h 23rd February 2005			*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _CMATRIX_H_
#define _CMATRIX_H_

#include <stdlib.h>
#include <stdio.h>

namespace myutils
{
/*Cannot accept objects of type class*/
template <typename T>
class CMatrix
{
public:
	/*Preserve public access for back-compatibility*/
	T **element;

protected:
	int protected_nrows;
	int protected_ncols;
	int initialized;

public:
	/*Default constructor*/	CMatrix()
	{
		initialized=0;
		initialize(0,0);
	}
	/*Constructor*/			CMatrix(int nrows, int ncols)
	{
		initialize(nrows,ncols);
	}
	/*Constructor*/			CMatrix(int nrows, int ncols, T value)
	{
		initialize(nrows,ncols);
		int i,j;
		for(i=0;i<nrows;i++)
			for(j=0;j<ncols;j++)
				element[i][j]=value;
	}
	/*Destructor*/			~CMatrix()
	{
		int i;
		for(i=protected_nrows-1;i>=0;i--) free((T*) element[i]);
		free((T**) element);
	}
	CMatrix<T>& initialize(int nrows, int ncols)
	{
		element=(T **) malloc((unsigned) nrows*sizeof(T*));
		if (!element) error("row allocation failure in Matrix::initialize()");

		int i;
		for(i=0;i<nrows;i++)
		{
			element[i]=(T *) malloc((unsigned) ncols*sizeof(T));
			if (!element[i]) error("column allocation failure in Matrix::initialize()");
		}
		protected_nrows=nrows;
		protected_ncols=ncols;
		initialized=1;
		return *this;
	}
	CMatrix<T>& resize(int nrows, int ncols)
	{
		int i;
		if (!initialized) initialize(nrows,ncols);
		else
		{
			if(nrows!=protected_nrows)
			{
				element=(T **) realloc(element,(unsigned) nrows*sizeof(T*));
				if (!element) error("row allocation failure in Matrix::resize()");

				if(nrows<protected_nrows)
				{
					for(i=protected_nrows-1;i>=nrows;i--)
						free ((T*) element[i]);
				}
				if(nrows>protected_nrows)
				{
					for(i=protected_nrows;i<nrows;i++)
					{
						element[i]=(T *) malloc((unsigned) protected_ncols*sizeof(T));
						if (!element[i]) error("column allocation failure 1 in Matrix::resize()");
					}
				}
				protected_nrows=nrows;
			}
			if(ncols!=protected_ncols)
			{
				for(i=0;i<nrows;i++)
				{
					element[i]=(T *) realloc(element[i],(unsigned) ncols*sizeof(T));
					if (!element[i]) error("column allocation failure 2 in Matrix::resize()");
				}
				protected_ncols=ncols;
			}

		}
		return *this;
	}
	int nrows(){return protected_nrows;}
	int ncols(){return protected_ncols;}
	void error(const char* error_text)
	{
		printf("Run-time error in Matrix::");
		printf("%s%\n", error_text);
		printf("Exiting to system...\n");
		exit(13);
	}
	/*Copy constructor*/	CMatrix(CMatrix<T>& mat)
	/*	Copy constructor for the following cases:
			Matrix mat2(mat);
			Matrix mat2=mat;
		and when Matrix is returned from a function	*/
	{
		initialize(mat.nrows(),mat.ncols());
		int i,j;
		for(i=0;i<protected_nrows;i++)
		{
			for(j=0;j<protected_ncols;j++)
			{
				element[i][j]=mat.element[i][j];
			}
		}
	}
	/*Assignment operator*/	CMatrix<T>& operator=(CMatrix<T>& mat)
	{
		if(this==&mat)return *this;
		resize(mat.nrows(),mat.ncols());
		int i,j;
		for(i=0;i<protected_nrows;i++)
		{
			for(j=0;j<protected_ncols;j++)
			{
				element[i][j]=mat.element[i][j];
			}
		}
		return *this;
	}
	/*Subscript operator*/inline T* operator[](int pos){return element[pos];};
};
};

#endif