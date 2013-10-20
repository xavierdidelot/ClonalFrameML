/*  Copyright 2012 Daniel Wilson.
 *
 *  vector.h
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
/*	vector.h 23rd February 2005				*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_VECTOR_H_
#define _MYUTILS_VECTOR_H_

#include <myerror.h>
#include <stdlib.h>
#include <stdio.h>
//#include <myutils.h>

namespace myutils
{
template <typename T>
class Vector
{
public:
	/*Preserve public access for back-compatibility*/
	T *element;

protected:
	int protected_size;
	int initialized;

public:
	/*Default constructor*/	Vector()
	{
		initialized=0;
		initialize(0);
	}
	/*Constructor*/			Vector(int size)
	{
		initialize(size);
	}
	/*Constructor*/			Vector(int size, T value)
	{
		initialize(size);
		int i;
		for(i=0;i<size;i++)
			element[i]=value;
	}
	/*Destructor*/			~Vector()
	{
		if(protected_size>0)
			delete[] element;
	}
	Vector<T>& initialize(int size)
	{
		element=new T[size];
		if (!element) error("allocation failure in Vector::initialize()");
		protected_size=size;
		initialized=1;
		return *this;
	}
	/*All current data is lost when the Matrix is resized*/
	Vector<T>& resize(int size)
	{
		if (!initialized) return initialize(size);
		if(size==protected_size)return *this;
		delete[] element;

		element=new T[size];
		if (!element) error("allocation failure in Vector::resize()");

		protected_size=size;
		return *this;
	}
	int size(){return protected_size;}
	int size() const {return protected_size;}
/*	void error(char* error_text)
	{
		printf("Run-time error in Vector::");
		printf("%s%\n", error_text);
		printf("Exiting to system...\n");
		exit(13);
	}*/
	/*Copy constructor*/	Vector(const Vector<T> &vec)
	/*	Copy constructor for the following cases:
			Vector vec2(vec);
			Vector vec2=vec;
		and when Vector is returned from a function	*/
	{
		initialize(vec.protected_size);
		int i;
		for(i=0;i<protected_size;i++)
			element[i]=vec.element[i];
	}
	/*Assignment operator*/	Vector<T>& operator=(const Vector<T>& vec)
	{
		resize(vec.size());
		int i;
		for(i=0;i<protected_size;i++)
			element[i]=vec.element[i];
		return *this;
	}
#ifdef _VECTOR_
	/*Copy constructor*/	Vector(const std::vector<T> &vec)
	/*	Copy constructor for the following cases:
			Vector vec2(vec);
			Vector vec2=vec;
		and when Vector is returned from a function	*/
	{
		initialize(vec.size());
		int i;
		for(i=0;i<protected_size;i++)
			element[i]=vec[i];
	}
	/*Assignment operator*/	Vector<T>& operator=(const std::vector<T>& vec)
	{
		resize(vec.size());
		int i;
		for(i=0;i<protected_size;i++)
			element[i]=vec[i];
		return *this;
	}
#endif
#ifdef _MYUTILS_DEBUG
	/*Subscript operator*/inline T& operator[](int pos){
		if(pos<0) error("Vector::operator[](int pos): pos<0");
		if(pos>=protected_size) error("Vector::operator[](int pos): pos>=size()");
		return element[pos];
	};
	/*Subscript operator*/inline const T& operator[](int pos) const {
		if(pos<0) error("Vector::operator[](int pos): pos<0");
		if(pos>=protected_size) error("Vector::operator[](int pos): pos>=size()");
		return element[pos];
	};
#else
	/*Subscript operator*/inline T& operator[](int pos){return element[pos];};
	/*Subscript operator*/inline const T& operator[](int pos) const {return element[pos];};
#endif
};
};

#endif // _MYUTILS_VECTOR_H_
