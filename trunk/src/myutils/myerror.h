/*  Copyright 2012 Daniel Wilson.
 *
 *  myerror.h
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
/*	myerror.h 23rd February 2005			*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_ERROR_H
#define _MYUTILS_ERROR_H

#include <stdio.h>
#include <stdlib.h>
// For use with MPI programs
#ifdef _MYUTILS_MPI_ABORT_ON_EXIT
#include <mpi.h>
#endif

namespace myutils
{
	inline void error(const char* error_text)
	{
		printf("ERROR: ");
		printf("%s\n", error_text);
#ifdef _MYUTILS_MPI_ABORT_ON_EXIT
		MPI_Abort(MPI_COMM_WORLD,13);
#endif
		exit(13);
	}

	inline void warning(const char* warning_text)
	{
		printf("WARNING: ");
		printf("%s\n", warning_text);
		return;
	}

};

#endif