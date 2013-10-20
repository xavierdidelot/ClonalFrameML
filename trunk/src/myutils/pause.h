/*  Copyright 2012 Daniel Wilson.
 *
 *  pause.h
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
/*	pause.h 23rd February 2005				*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_PAUSE_H_
#define _MYUTILS_PAUSE_H_

#ifdef _WIN32

	#include <conio.h>
	#include <stdio.h>

	namespace myutils
	{
		inline void pause()
		{
			printf("\nPress any key\n");
			int ch=-99;
			while (ch==-99)
				ch=_getch();
		}
		inline void silent_pause()
		{
			int ch=-99;
			while (ch==-99)
				ch=_getch();
		}
	};

#else

	namespace myutils
	{
		inline void pause() {}
	};

#endif

#endif