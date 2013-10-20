/*  Copyright 2012 Daniel Wilson.
 *
 *  census.h
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
/*	census.h 28th August 2009				*/
/*											*/
/*	Keeps track of the membership of a		*/
/*	finite number of individuals among a	*/
/*	finite number of populations.			*/
/*											*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_CENSUS_H_
#define _MYUTILS_CENSUS_H_

#include <myerror.h>
#include <vector.h>
#include <utils.h>
#include <iostream>
//#include <myutils.h>

using std::cout;
using std::endl;

namespace myutils {
class Census {
public:
	/*Default constructor*/
	Census() {
		Vector<int> where(0);
		initialize(0,0,where);
	}
	/*Constructor*/
	Census(const unsigned int npop, const unsigned int nind) {
		Vector<int> where(nind,0);
		initialize(npop,nind,where);
	}
	/*Constructor*/
	Census(const unsigned int npop, const unsigned int nind, Vector<int> &where) {
		initialize(npop,nind,where);
	}
	/*Copy constructor*/
	Census(const Census& cen) {
		_npop = cen._npop;
		_nind = cen._nind;
		_where = cen._where;
		_who = cen._who;
		_index = cen._index;
		_mind = cen._mind;
		_cind = cen._cind;
	}
	/*Assignment operator*/
	Census& operator=(const Census& cen) {
		_npop = cen._npop;
		_nind = cen._nind;
		_where = cen._where;
		_who = cen._who;
		_index = cen._index;
		_mind = cen._mind;
		_cind = cen._cind;
		return *this;
	}
	Census& initialize(const unsigned int npop, const unsigned int nind, Vector<int> &where) {
		if(npop<0) error("Census::Census(): number of populations must be non-negative");
		if(nind<0) error("Census::Census(): number of individuals must be non-negative");
		/* Accept the arguments */
		_npop = npop;
		_nind = nind;
		/* Initialize the membership lists */
		_where = Vector<int>(_nind);
		_mind = Vector<int>(_npop,0);
		int i;
		for(i=0;i<_nind;i++) {
			if(where[i]<0) error("Census::Census(): population cannot be negative");
			if(where[i]>=_npop) error("Census::Census(): population number exceeds maximum");
			_where[i] = where[i];
			_mind[where[i]]++;
		}
		/* Calculate the cumulative membership numbers */
		_cind = Vector<int>(_npop,0);
		int p;
		for(p=1;p<_npop;p++) {
			_cind[p] = _cind[p-1] + _mind[p-1];
		}
		if(_npop>0 && _cind[_npop-1]+_mind[_npop-1]!=_nind) error("Census::Census(): number of individuals doesn't match");
		/* Initialize the who list */
		_who = Vector<int>(_nind);
		_index = Vector<int>(_nind);
		Vector<int> _tind(_npop,0);
		for(i=0;i<_nind;i++) {
			const int pop = _where[i];
			const int ix = _cind[pop]+_tind[pop];
			_who[ix] = i;
			_index[i] = ix;
			++_tind[pop];
		}
		return *this;
	}
	/*Destructor*/
	~Census() {}
	/*Simple functions*/
	int npop() {return _npop;}
	int nind() {return _nind;}
	int nind(const int p) {return _mind[p];}
	Vector<int> where() {return _where;}
	int where(const int i) {return _where[i];}
	Vector<int> who(const int p) {
		if(p<0 || p>=_npop)  error("Census::who(): Population p out of range");
		Vector<int> ret(_mind[p]);
		int i;
		for(i=0;i<_mind[p];i++) {
			ret[i] = _who[_cind[p]+i];
		}
		return ret;
	}
	int who(const int p, const int i) {return _who[_cind[p]+i];}
	int ferocious_who(const int p, const int i) {
		if(p<0 || p>=_npop) error("Census::who(): Population p out of range");
		if(i<0 || i>=_mind[p]) error("Census::who(): Index i out of range for population p");
		return _who[_cind[p]+i];
	}
	int meek_who(const int p, const int i) {
		if(p<0 || p>=_npop) return -1;
		if(i<0 || i>=_mind[p]) return -1;
		return _who[_cind[p]+i];
	}
	/* Last individual in the population */
	int last(const int p) {
		if(p<0 || p>=_npop) error("Census::last(): population out of range");
		if(_mind[p]==0) error("Census::last(): population is empty");
		return _who[_cind[p]+_mind[p]-1];
	}
	/*Not-so simple functions*/
	int migrate(const int from, const int to) {
		const int ind = last(from);
		migrate(ind,from,to);
		return ind;
	}
	/* ind is the absolute index of the individual */
	Census& migrate(const int ind, const int from, const int to) {
		if(from==to) return *this;
		if(from<0 || from>=_npop) error("Census::migrate(): donor population out of range");
		if(_where[ind]!=from) error("Census::migrate(): individual is not member of donor population");
		if(to<0 || to>=_npop) error("Census::migrate(): recipient population out of range");
		if(from<to) {
			/* Swap into last position */ {
				const int ifrom = ind;
				const int ix_from = _index[ifrom];
				const int ix_to = _cind[from]+_mind[from]-1;
				const int ito = _who[ix_to];
				SWAP(_who[ix_from],_who[ix_to]);
				SWAP(_index[ifrom],_index[ito]);				
			}
			int p;
			/* Swap into last position of successive populations */
			for(p=from;p<to;p++) {
				/* 1. Add to new pop */
				--_mind[p];
				++_mind[p+1];
				--_cind[p+1];
				/* 2. Swap from first to last position */
				const int ix_from = _cind[p+1];
				const int ix_to = _cind[p+1]+_mind[p+1]-1;
				const int ifrom = _who[ix_from];
				const int ito = _who[ix_to];
				SWAP(_who[ix_from],_who[ix_to]);
				SWAP(_index[ifrom],_index[ito]);
			}
			/* Update _where */
			_where[ind] = to;
		}
		else {
			/* Swap into first position */ {
				const int ifrom = ind;
				const int ix_from = _index[ifrom];
				const int ix_to = _cind[from];
				const int ito = _who[ix_to];
				SWAP(_who[ix_from],_who[ix_to]);
				SWAP(_index[ifrom],_index[ito]);				
			}
			int p;
			/* Swap into first position of successive populations */
			for(p=from;p>to;p--) {
				/* 1. Add to new pop */
				--_mind[p];
				++_mind[p-1];
				++_cind[p];
				/* 2. Swap from last to first position */
				const int ix_from = _cind[p-1]+_mind[p-1]-1;
				const int ix_to = _cind[p-1];
				const int ifrom = _who[ix_from];
				const int ito = _who[ix_to];
				SWAP(_who[ix_from],_who[ix_to]);
				SWAP(_index[ifrom],_index[ito]);				
			}
			/* Update _where */
			_where[ind] = to;
		}
		return *this;
	}
	Census& inspect() {
		int i;
		cout << "_where = {" << _where[0];
		for(i=1;i<_nind;i++) cout << " " << _where[i];
		cout << "}" << endl;
		cout << "_who   = {" << _who[0];
		for(i=1;i<_nind;i++) cout << " " << _who[i];
		cout << "}" << endl;
		cout << "_index = {" << _index[0];
		for(i=1;i<_nind;i++) cout << " " << _index[i];
		cout << "}" << endl;
		cout << "_mind  = {" << _mind[0];
		for(i=1;i<_npop;i++) cout << " " << _mind[i];
		cout << "}" << endl;
		cout << "_cind  = {" << _cind[0];
		for(i=1;i<_npop;i++) cout << " " << _cind[i];
		cout << "}" << endl;
		
		return *this;
	}
	
protected:
	/* Number of populations */
	int _npop;
	/* Number of individuals */
	int _nind;
	/*	_where[i], i=0.._nind-1, has value [0,_npop-1],
		Population to which individual i belongs */
	Vector<int> _where;
	/*	_who[i], i=0.._nind-1, has value [0,_nind-1],
		Collapsed unordered list of individuals belonging to the
		population to which i corresponds */
	Vector<int> _who;
	/*	_index[i], i=0.._nind-1, has value [0,_nind-1],
		Position of individual i in vector _who */
	Vector<int> _index;
	/*	_mind[p], p=0.._npop-1, has value [0,_nind],
		Number of members of population p */
	Vector<int> _mind;
	/*	_cind[p], p=0.._npop-1, has value [0,_nind],
		Cumulative number of members of population p */
	Vector<int> _cind;
};
};

#endif // _MYUTILS_CENSUS_H_
