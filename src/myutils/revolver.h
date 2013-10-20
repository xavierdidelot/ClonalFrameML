/*  Copyright 2012 Daniel Wilson.
 *
 *  revolver.h
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
/*	revolver.h 28th August 2009				*/
/*	The revolver container has a fixed		*/
/*	number of elements that it releases		*/
/*	and takes back as required.	Its purpose	*/
/*	is to avoid unnecessary memory alloc-	*/
/*	ation and freeing.						*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_REVOLVER_H_
#define _MYUTILS_REVOLVER_H_

#include <myerror.h>
#include <stdlib.h>
#include <stdio.h>
//#include <myutils.h>

namespace myutils
{
template <typename T>
class Revolver
{
public:
	/*Preserve public access for back-compatibility*/
	T **element;

protected:
	int protected_size;
	int protected_available;

public:
	/*Default constructor*/	Revolver()
	{
		initialize(0);
	}
	/*Constructor*/			Revolver(int size)
	{
		initialize(size);
	}
	/*Constructor*/			Revolver(int size, T &value)
	{
		initialize(size,value);
	}
	/*Destructor*/			~Revolver()
	{
		// Do not delete items in use!!
		if(!full()) error("Revolver::~Revolver(): not full");
		for(i=0;i<protected_size;i++) delete element[i];
		if(protected_size>0)
			delete[] element;
	}
	bool full() const { return protected_available==protected_size; }
	bool full() { return protected_available==protected_size; }
	bool empty() const { return protected_available==0; }
	bool empty() { return protected_available==0; }
	int size(){return protected_size;}
	int size() const {return protected_size;}
	int navail(){return protected_available;}
	int navail() const {return protected_available;}
	Revolver<T>& initialize(int size)
	{
		element=new T*[size];
		if(!element) error("Revolver::initialize() allocation failure");
		int i;
		for(i=0;i<size;i++) {
			element[i] = new T;
			if(!element[i]) error("Revolver::initialize() allocation failure");
		}
		protected_size=size;
		protected_available=size;
		return *this;
	}
	Revolver<T>& initialize(int size, T &value)
	{
		element=new T*[size];
		if(!element) error("Revolver::initialize() allocation failure");
		int i;
		for(i=0;i<size;i++) {
			element[i] = new T(value);
			if(!element[i]) error("Revolver::initialize() allocation failure");
		}
		protected_size=size;
		protected_available=size;
		return *this;
	}
#ifdef _MYUTILS_DEBUG
	/* NB:- order is not stable in Revolver */
	/*Subscript operator*/inline T* operator[](int pos){
		if(pos<0) error("Revolver::operator[](int pos): pos<0");
		if(pos>=protected_size) error("Revolver::operator[](int pos): pos>=size()");
		return element[pos];
	};
#else
	/* NB:- order is not stable in Revolver */
	/*Subscript operator*/inline T* operator[](int pos){return element[pos];};
#endif
	/* Release an element for use */
	T* pop() {
		if(empty()) {
			if(size()==0) error("Revolver::pop(): zero-sized container");
			error("Revolver::pop(): empty container");
		}
		--protected_available;
		return element[protected_available];
	}
#ifdef _MYUTILS_DEBUG
	/* Return an element to the container, checking that it belongs to the container */
	Revolver& push(T* val) {
		if(full()) error("Revolver::push(): full container");
		int i;
		for(i=protected_available;i<protected_size;i++) {
			if(element[i]==val) {
				if(i!=protected_available)
					SWAP(element[protected_available],element[i]);
				++protected_available;
				return *this;
			}
		}
		error("Revolver::push(): element does not belong to container");
		return *this;
	}
#else
	/* Return an element to the container */
	Revolver& push(T* val) {
		if(full()) error("Revolver::push(): full container");
		element[protected_available] = val;
		++protected_available;
		return *this;
	}
#endif
};
};

#endif // _MYUTILS_REVOLVER_H_
