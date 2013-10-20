/*  Copyright 2012 Daniel Wilson.
 *
 *  tsv.h
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
#ifndef _TSV_H_
#define _TSV_H_

#pragma warning(disable: 4786)
#include <myerror.h>
using myutils::error;
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <matrix.h>
using myutils::Matrix;
using namespace std;

class tsv {
public:
	bool coutput;
	Matrix<string> data;
	vector<string> fieldname;
	map<string, int> fieldnum;
	vector< vector<string> > fieldvalue;
	const int DUMPMAX;
	
	tsv(const int dumpmax = 1000) : DUMPMAX(dumpmax) {
		coutput=false;
	}
	tsv& read(const char* infilename)
	{
		ifstream infile(infilename);
		if(!infile.is_open()) error("Could not open file");
		int nfields = 0;
		int character = 0;
		string this_fieldname = "";
		fieldname.resize(0);
		if(coutput) cout << "Fields found: ";
		while(!infile.eof())
		{
//			if(!infile.good()) error("Problem reading file - is buffer too small?");
			character = infile.get();
			if(character=='\t')
			{
				if(coutput) cout << this_fieldname << " ";
				fieldname.push_back(this_fieldname);
				fieldnum[this_fieldname]=nfields;
				this_fieldname = "";
				++nfields;
			}
			else if(character=='\r'||character=='\n'||character==-1)
			{
				if(coutput) cout << this_fieldname << " ";
				fieldname.push_back(this_fieldname);
				fieldnum[this_fieldname]=nfields;
				this_fieldname = "";
				++nfields;
				character = infile.peek();
				if(character=='\r'||character=='\n'||character==-1)
					infile.get();
				break;
			}
			else
				this_fieldname += (char)character;
		}
		if(coutput) cout << "(" << nfields << " fields in total)" << endl << flush;

		int nrows = 0, ntries = 0;
		char* dump = new char[DUMPMAX];
		while(!infile.eof())
		{
			infile.getline(dump,DUMPMAX);
			if(dump[0]!='\0')/*check the line isn't blank using the end-of-string character*/
				++nrows;
			++ntries;
			if(coutput && ntries%1000==0) cout << "\r" << ntries << " attempts, " << nrows << " rows so far" << flush;
		}	cout << endl;
		if(coutput) cout << "Found " << nrows << " rows of data" << endl << flush;
		data.resize(nrows,nfields);

		infile.close();
		ifstream infile2(infilename);

		infile2.getline(dump,DUMPMAX);
		delete[] dump;

		int row = 0;
		int col = 0;
		string value = "";
		while(!infile2.eof())
		{
			character = infile2.get();
			if(character=='\t')
			{
				if(row<nrows && col<nfields)
					data[row][col] = value;
				value = "";
				++col;
			}
			else if(character=='\r'||character=='\n'||character==-1)
			{
				if(row<nrows && col<nfields)
					data[row][col] = value;
				value = "";
				++row;
				col = 0;
				character = infile2.peek();
				if(character=='\r' || character=='\n' || character==-1)
					infile2.get();
			}
			else
				value += (char)character;
		}
		record_values();
		return *this;
	};

	tsv& record_values()
	{
		fieldvalue.resize(fieldname.size());
		unsigned int f;
		for(f=0;f<fieldname.size();f++)
		{
			fieldvalue[f].resize(0);
			fieldvalue[f].push_back(data[0][f]);
			for(int i=0;i<data.nrows();i++)
			{
				string val = data[i][f];
				bool acopy = false;
				unsigned int k;
				for(k=0;k<fieldvalue[f].size();k++)
					if(val==fieldvalue[f][k])
					{
						acopy = true;
						break;
					}
				if(!acopy)
				{
					fieldvalue[f].push_back(val);
				}
			}
		}
		return *this;
	};

	vector<int> n_values()
	{
		vector<int> result(data.ncols(),0);
		int f;
		for(f=0;f<data.ncols();f++)
			result[f]=(int)fieldvalue[f].size();
		return result;
	};
	bool value_exist(string field, string value)
	{
		int f=fieldnum[field];
		unsigned int i;
		for(i=0;i<fieldvalue[f].size();i++)
			if(fieldvalue[f][i]==value)
				break;
		if(i>=fieldvalue[f].size()) return false;
		return true;
	};
	bool field_exist(string f)
	{
		bool exists = false;
		unsigned int i;
		for(i=0;i<fieldname.size();i++)
			if(f==fieldname[i]) {
				exists = true;
				break;
			}
		return exists;
	};
	int which(string f) {
		map<string,int>::iterator m = fieldnum.find(f);
		if(m==fieldnum.end()) return -1;
		return m->second;
	}
};

#endif //_TSV_H_