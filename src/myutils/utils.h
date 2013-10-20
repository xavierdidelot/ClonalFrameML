/*  Copyright 2012 Daniel Wilson.
 *
 *  utils.h
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
/*	utils.h 23rd February 2005				*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_UTILS_H_
#define _MYUTILS_UTILS_H_

namespace myutils {
template<typename T>
void SWAP(T &a, T &b) {
	T c = a;
	a = b;
	b = c;
}

template<typename T>
T MIN(T a, T b) {
	return (a<b) ? a : b;
}

template<typename T>
T MAX(T a, T b) {
	return (a<b) ? b : a;
}
};

#endif // _MYUTILS_UTILS_H_