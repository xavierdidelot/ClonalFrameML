/*  Copyright 2012 Daniel Wilson.
 *
 *  myutils.h
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
/*	myutils.h 23rd February 2005			*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_H_
#define _MYUTILS_H_

#pragma warning(disable: 4786)

/*Includes all header files in the myutils directory*/
/*#include "cmatrix.h"
#include "matrix.h"
#include "random.h"
#include "error.h"
#include "DNA.h"
#include "vector.h"*/

#include <myerror.h>
#include <utils.h>
#include <cmatrix.h>
#include <vector.h>
#include <matrix.h>
#include <lotri_matrix.h>
#include <random.h>
#include <DNA.h>
#include <pause.h>
#include <sort.h>

/*#include "controlwizard.h" /* has problems in Linux with pointers */
/*#include "pause.h"	/* removed because conio.h is not standard */

#endif
