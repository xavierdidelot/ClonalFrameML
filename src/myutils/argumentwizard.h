/*  Copyright 2012 Daniel Wilson.
 *
 *  argumentwizard.h
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
/*	argumentwizard.h 23rd February 2005		*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _ARGUMENT_WIZARD_H_
#define _ARGUMENT_WIZARD_H_

#pragma warning(disable: 4786)

#include <string>
#include <list>
#include <map>
#include <vector>
#include <iostream>
#include <ctype.h>
#include "myerror.h"
#include <sstream>

namespace myutils
{
#ifndef _CONTROL_AND_ARGUMENT_WIZARD_TYPES_
#define _CONTROL_AND_ARGUMENT_WIZARD_TYPES_
typedef void RTRV;//functions that retrieve the data
typedef void GENERIC;//for the generic pointers
enum DATA_TYPE {TP_UNRECOGNISED,TP_INT,TP_DOUBLE,TP_STRING,TP_VEC_INT,TP_VEC_DOUBLE,TP_EXT_VEC_DOUBLE};
#endif // _CONTROL_AND_ARGUMENT_WIZARD_TYPES_

class ArgumentWizard {
/*MEMBER VARIABLES*/
public:
	std::list<std::string> required;
	bool coutput;
	bool unrecognised;
	bool got_required;
	bool case_sensitive;
	bool fail_noprefix;

	/* used to avoid function pointers in selecting data-read function */
	DATA_TYPE switcher;

protected:
	std::map<std::string,DATA_TYPE> label_map;
	std::map<std::string,GENERIC*> data_map;
	int argc,argn;
	std::vector<std::string> argv;

/*MEMBER FUNCTIONS*/
public:
	ArgumentWizard(){set_defaults();}

	void read_input(const int argc_in, const char* argv_in[]) {
		argc = argc_in;
		argv = std::vector<std::string>(argc);
		int i;
		for(i=0;i<argc;i++) argv[i] = argv_in[i];

		required.unique();

		std::string label;
		argn = 1;
		while(argn<argc) {
			if(read_label(label)) {
				data_format(label);
				read_data(label);
				required.remove(label);
			}
		}
		if(coutput) std::cout << "Finished reading in control file." << std::endl << std::endl;
		got_required = auto_check_required();
	}

	// GENERIC* lets any type of pointer be passed to the function
	void add_item(std::string label, const DATA_TYPE data_type, GENERIC* location) {
		if(!case_sensitive) remove_case(label);
		label_map[label] = data_type;
		data_map[label] = location;
	}

	// These are essential items
	void add_ITEM(std::string label, const DATA_TYPE data_type, GENERIC* location) {
		if(!case_sensitive) remove_case(label);
		label_map[label] = data_type;
		data_map[label] = location;
		required.push_back(label);
	}

	// Returns false if some required items are not found
	bool check_required() {
		if(!got_required) {
			std::cout << "The following required items have not been found: ";
			std::list<std::string>::iterator i;
			for(i=required.begin();i!=required.end();i++) std::cout << *i << " ";
			std::cout << std::endl;
		}
		else std::cout << "All required items were found" << std::endl;
		return got_required;
	}

protected:
	void set_defaults() {
		coutput = true;
		unrecognised = true;
		case_sensitive = false;
		fail_noprefix = true;
	}

	void remove_case(std::string &s) {
		int i;
		for(i=0;i<(int)s.length();i++) s[i] = tolower(s[i]);
	}

	/* Returns true if a label is found */
	bool read_label(std::string &word) {
		if(argn>=argc) error("Syntax error in ArgumentWizard::read_label: exceeded number of arguments");
		word = argv[argn];
		if(word[0]!='-') {
			if(fail_noprefix) error("Syntax error in ArgumentWizard::read_label: option must be prefixed with a \'-\'");
			++argn;
			return false;
		}
		
		std::string word2 = std::string(word.length()-1,' ');
		int i;
		for(i=1;i<(int)word.length();i++) word2[i-1] = word[i];
		word = word2;

		if(!case_sensitive) remove_case(word);
		++argn;
		return true;
	}

	void data_format(std::string &label) {
		label_map[label];
		switcher = label_map[label];
	}

	// Returns false if some required items are not found
	bool auto_check_required() {
		return (required.size()==0);
	}

protected:
	RTRV function_get_unrecognised(std::string &label) {
		if((label.size()>0)&&(coutput || unrecognised))
			printf("Label \"%s\" not recognised.\n",label.c_str());
	}
	template<typename T>
	RTRV function_get_single(T dummy, std::string &label) {
		if(argn>=argc) error("Syntax error in ArgumentWizard::function_get_single(): exceeded number of arguments");
		std::string word = argv[argn];
		//if(word[0]=='-') error("Syntax error in ArgumentWizard::function_get_single(): expecting a value but got an option");

		std::stringstream s;
		s << word;
		T value;
		s >> value;
		GENERIC* ptr = data_map[label];
		(*(static_cast<T*>(ptr))) = value;
		if(coutput) std::cout << label << " = " << value << std::endl;
		++argn;
	}
	template<typename T>
	RTRV function_get_vector(T dummy, std::string &label) {
		if(argn>=argc) error("Syntax error in ArgumentWizard::function_get_vector(): exceeded number of arguments");
		std::string word;
		GENERIC* g_ptr = data_map[label];
		std::vector<T>* ptr = static_cast<std::vector<T>*>(g_ptr);
		ptr->clear();

		if(coutput) std::cout << label << " = ";
		while(true) {
			if(argn==argc) break;
			word = argv[argn];
			if(word[0]=='-') break;
			stringstream s;
			s << word;
			T value;
			s >> value;
			ptr->push_back(value);
			if(coutput) std::cout << word << " ";
			++argn;
		}
		if(coutput) std::cout << std::endl;
	}
	RTRV function_get_string(std::string &label) {
		if(argn>=argc) error("Syntax error in ArgumentWizard::function_get_single(): exceeded number of arguments");
		std::string word = argv[argn];
		//if(word[0]=='-') error("Syntax error in ArgumentWizard::function_get_single(): expecting a value but got an option");

		GENERIC* ptr = data_map[label];
		(*(static_cast<std::string*>(ptr))) = word;
		if(coutput) std::cout << label << " = " << word << std::endl;
		++argn;
	}
	RTRV function_get_external_vector_double(std::string &label) {
		error("ArgumentWizard:: TP_EXT_VEC_DOUBLE not available");
	}
	void read_data(std::string &label) {
		switch(switcher) {
		case TP_UNRECOGNISED:
			function_get_unrecognised(label);			break;
		case TP_INT:
			function_get_single((int)0,label);			break;
		case TP_DOUBLE:
			function_get_single((double)0,label);		break;
		case TP_STRING:
			function_get_string(label);					break;
		case TP_VEC_INT:
			function_get_vector((int)0,label);			break;
		case TP_VEC_DOUBLE:
			function_get_vector((double)0,label);		break;
		case TP_EXT_VEC_DOUBLE:
			function_get_external_vector_double(label);	break;
		default:
			function_get_unrecognised(label);			break;
		}
	}
};		// class ArgumentWizard
};		// namespace myutils

#endif	// _ARGUMENT_WIZARD_H_
