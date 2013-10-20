/*  Copyright 2012 Daniel Wilson.
 *
 *  sort.h
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
#ifndef _MYUTILS_SORT_H_
#define _MYUTILS_SORT_H_

#include <vector>
#include <functional>

namespace myutils {

/*	WARNING: this class has very limited utility.
	
	Syntax:
	sort(sortme.begin(),sortme.end(),sort_by_vector<T>(sortby));

	where sortby is the vector of interest, if sortme is a vector
	that starts of as the indeces of sortby, i.e. 0,1,2,...,size()-1
	then following the sort, it will be reordered according to sortby.
	*/
template<typename T>
class sort_by_vector : public std::binary_function<int,int,bool>
{
	const vector<T> &sort_by;
public:
	sort_by_vector(const vector<T> &sort_by_in) : sort_by(sort_by_in) {}

	bool operator()(int a, int b) const
	{
		return (sort_by.at(a)<sort_by.at(b));
	}
};

};		// namespace myutils

#endif	// _MYUTILS_SORT_H_