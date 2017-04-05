/*  Copyright 2012 Daniel Wilson.
 *
 *  controlwizard.h
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
/************************************************/
/*	control_wizard.h 23rd February 2005		*/
/*	(c) Danny Wilson.							*/
/*	www.danielwilson.me.uk						*/
/************************************************/

#ifndef _CONTROL_WIZARD_H_
#define _CONTROL_WIZARD_H_

#pragma warning(disable: 4786)

#include <vector>
#include <map>
#include <list>
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <ctype.h>
#include <limits>

//using namespace std;

namespace myutils
{
#ifndef _CONTROL_AND_ARGUMENT_WIZARD_TYPES_
#define _CONTROL_AND_ARGUMENT_WIZARD_TYPES_
typedef void RTRV;//functions that retrieve the data
typedef void GENERIC;//for the generic pointers
enum DATA_TYPE {TP_UNRECOGNISED,TP_INT,TP_DOUBLE,TP_STRING,TP_VEC_INT,TP_VEC_DOUBLE,TP_EXT_VEC_DOUBLE,TP_VEC_STRING};
#endif // _CONTROL_AND_ARGUMENT_WIZARD_TYPES_

class ControlWizard
{
/*MEMBER VARIABLES*/
public:
	std::vector<int> line_delimiters;
	std::vector<int> label_delimiters;
	std::vector<int> white_space;
	std::vector<int> elem_delimiters;
	std::vector<int> rem_delimiters;
	std::vector<int> eof_delimiters;
	std::list<std::string> required;
	bool coutput;			// Set to true to print out all comments
	bool unrecognised;		// Set to true to print out unrecognised options
	bool got_required;
	bool case_sensitive;
	bool _EOF;

	/* used to avoid function pointers in selecting data-read function */
	DATA_TYPE switcher;

protected:
	std::map<std::string,DATA_TYPE> label_map;
	std::map<std::string,GENERIC*> data_map;

/*MEMBER FUNCTIONS*/
public:
	ControlWizard(){set_defaults();}

	void read_input(const char* filename)
	{
		std::ifstream infile(filename);
		if(infile.is_open()==false) {
			string errTxt(filename);
			errTxt += " not found";
			error(errTxt.c_str());
		}

		required.unique();

		std::string label;
		while(eof(infile)==false)
		{
			if(read_label(infile,label))
			{
				data_format(label);
				//(*this.*read_data)(infile,label);
				read_data(infile,label);
				required.remove(label);
			}
		}
		if(coutput)printf("Finished reading in control file.\n\n");
		got_required=auto_check_required();
		infile.close();
	}

	std::ifstream& read_input(std::ifstream &infile)
	{
		if(infile.is_open()==false)error("File not found");

		required.unique();

		std::string label;
		while(eof(infile)==false)
		{
			if(read_label(infile,label))
			{
				data_format(label);
				//(*this.*read_data)(infile,label);
				read_data(infile,label);
				required.remove(label);
			}
		}
		if(coutput)printf("Finished reading in control file.\n\n");
		got_required=auto_check_required();
		return infile;
	}

	void add_item(std::string label,const DATA_TYPE data_type,GENERIC* location)
	//GENERIC* lets any type of pointer be passed to the function
	{
		if(!case_sensitive) {
			//std::transform(label.begin(),label.end(),label.begin(),tolower);
			int i;
			for(i=0;i<(int)label.length();i++) label[i] = tolower(label[i]);
		}
		label_map[label]=data_type;
		//printf("Assigned label \"%s\"\n",label.c_str());
		data_map[label]=location;
	}

	void add_ITEM(std::string label,const DATA_TYPE data_type,GENERIC* location)
	//These are essential items
	{
		if(!case_sensitive) {
			//std::transform(label.begin(),label.end(),label.begin(),tolower);
			int i;
			for(i=0;i<(int)label.length();i++) label[i] = tolower(label[i]);
		}
		label_map[label]=data_type;
		data_map[label]=location;
		required.push_back(label);
	}

	bool check_required()
	//Returns false if some required items are not found
	{
		if(!got_required)
		{
			printf("The following required items have not been found: ");
			//std::copy(required.begin(),required.end(),std::ostream_iterator<std::string>(std::cout," "));
			std::list<std::string>::iterator i;
			for(i=required.begin();i!=required.end();i++) cout << *i << " ";
			cout << endl;
		}
		else printf("All required items were found\n");
		return got_required;
	}

	bool eof(std::ifstream &infile)
	{
		return (infile.eof() || _EOF);
	}

	char get(std::ifstream &infile)
	{
		char ch = infile.get();
		int i;
		if(!_EOF)
			for(i=0;i<(int)eof_delimiters.size();i++)
				if((int)ch==eof_delimiters[i]) {
					_EOF = true;
					break;
				}
		return ch;
	}
	
protected:
	void set_defaults()
	{
		coutput=true;
		unrecognised=true;
		case_sensitive=false;
		white_space.push_back(' ');
		white_space.push_back(-1);
		white_space.push_back(10);
		white_space.push_back(13);
		label_delimiters.push_back('=');
		line_delimiters.push_back(10);
		line_delimiters.push_back(13);
		elem_delimiters.push_back(',');
		//elem_delimiters.push_back('\t');
		rem_delimiters.push_back('#');
		_EOF = false;
	}

	void error(const char* error_text)
	{
		printf("Run-time error in ControlWizard::");
		printf("%s\n", error_text);
		printf("Exiting to system...\n");
		exit(13);
	}

	bool read_label(std::ifstream &infile, std::string &word)
	/*Returns true if a label is found*/
	{
		int character;
		word="";
		bool label_delim_found=false;
		bool line_delim_found=false;
		bool include_char=true;

		while(eof(infile)==false)
		{
			character=get(infile);

			label_delim_found=false;
			include_char=true;
			int i;
			for(i=0;i<(int)white_space.size();i++)
			{
				if(character==white_space[i])include_char=false;
			}
			for(i=0;i<(int)line_delimiters.size();i++)
			{
				if(character==line_delimiters[i])
				{
					include_char=false;
					line_delim_found=true;
				}
			}
			for(i=0;i<(int)label_delimiters.size();i++)
			{
				if(character==label_delimiters[i])
				{
					include_char=false;
					label_delim_found=true;
				}
			}
			for(i=0;i<(int)rem_delimiters.size();i++)
			{
				if(character==rem_delimiters[i])
				{
					snail(infile);
					return false;
				}
			}

			if(include_char==true)word += static_cast<char>(character);
			if(line_delim_found==true)
			{
				if(word.size()>0)printf("Incomplete line \"%s\"\n",word.c_str());
				break;
			}
			if(label_delim_found==true)break;
		}
		
		if(!case_sensitive) {
			//std::transform(word.begin(),word.end(),word.begin(),tolower);
			int ii;
			for(ii=0;ii<(int)word.length();ii++) word[ii] = tolower(word[ii]);
		}
		//cout << "Returning string: (" << word << ")" << endl;
		return label_delim_found;
	}

	void data_format(std::string& label)
	{
		label_map[label];
		//DATA_TYPE switcher=label_map[label];
		switcher=label_map[label];

		/*switch(switcher)
		{
		case TP_UNRECOGNISED:
			read_data=function_get_unrecognised;break;
		case TP_INT:
			read_data=function_get_int;break;
		case TP_DOUBLE:
			read_data=function_get_double;break;
		case TP_STRING:
			read_data=function_get_string;break;
		case TP_VEC_INT:
			read_data=function_get_vector_int;break;
		case TP_VEC_DOUBLE:
			read_data=function_get_vector_double;break;
		case TP_EXT_VEC_DOUBLE:
			read_data=function_get_external_vector_double;break;
		default:
			read_data=function_get_unrecognised;break;
		}*/
	}

	bool auto_check_required()
	//Returns false if some required items are not found
	{
		if(required.size()>0)return false;
		return true;
	}

	void get_single(std::ifstream &infile,std::string &word)
	{
		int character;
		word="";
		bool line_delim_found=false;
		bool include_char=true;

		while(eof(infile)==false)
		{
			character=get(infile);

			include_char=true;
			int i;
			for(i=0;i<(int)white_space.size();i++)
			{
				if(character==white_space[i])include_char=false;
			}
			for(i=0;i<(int)line_delimiters.size();i++)
			{
				if(character==line_delimiters[i])
				{
					include_char=false;
					line_delim_found=true;
				}
			}

			if(include_char==true)word += static_cast<char>(character);
			if(line_delim_found==true)break;
		}
	}

	bool get_element(std::ifstream &infile,std::string &word)
	/*Returns false when line delimiter or EOF is reached*/
	{
		int character;
		word="";
		bool line_delim_found=false;
		bool include_char=true;
		bool elem_delim_found=false;

		while(eof(infile)==false)
		{
			character=get(infile);

			include_char=true;
			int i;
			for(i=0;i<(int)white_space.size();i++)
			{
				if(character==white_space[i])include_char=false;
			}
			for(i=0;i<(int)line_delimiters.size();i++)
			{
				if(character==line_delimiters[i])
				{
					include_char=false;
					line_delim_found=true;
				}
			}
			for(i=0;i<(int)elem_delimiters.size();i++)
			{
				if(character==elem_delimiters[i])
				{
					include_char=false;
					elem_delim_found=true;
				}
			}

			if(include_char==true)word += static_cast<char>(character);
			if(line_delim_found==true)return false;
			if(elem_delim_found==true)return true;
		}
		return false;
	}

	void get_multiple(std::ifstream &infile,std::vector<std::string> &words)
	{
		bool loop=true;
		int elem=(int)words.size()-1;

		while(loop==true)
		{
			words.push_back("");
			++elem;
			loop=get_element(infile,words[elem]);
		}
	}

	void snail(std::ifstream &infile)
	/*Proceeds to next line*/
	{
		int character;
		bool line_delim_found=false;
		
		while(eof(infile)==false)
		{
			character=get(infile);

			int i;
			for(i=0;i<(int)line_delimiters.size();i++)
			{
				if(character==line_delimiters[i])line_delim_found=true;
			}

			if(line_delim_found==true)break;
		}
	}
	
protected:
//	RTRV (ControlWizard::*read_data)(std::ifstream &infile, std::string &label);
	void read_data(std::ifstream &infile, std::string &label) {
		switch(switcher) {
		case TP_UNRECOGNISED:
			function_get_unrecognised(infile,label);			break;
		case TP_INT:
			function_get_int(infile,label);						break;
		case TP_DOUBLE:
			function_get_double(infile,label);					break;
		case TP_STRING:
			function_get_string(infile,label);					break;
		case TP_VEC_INT:
			function_get_vector_int(infile,label);				break;
		case TP_VEC_DOUBLE:
			function_get_vector_double(infile,label);			break;
		case TP_EXT_VEC_DOUBLE:
			function_get_external_vector_double(infile,label);	break;
		case TP_VEC_STRING:
			function_get_vector_string(infile,label);			break;
		default:
			function_get_unrecognised(infile,label);			break;
		}
	}
	RTRV function_get_unrecognised(std::ifstream &infile, std::string &label)
	{
		if((label.size()>0)&&(coutput || unrecognised))
			printf("Label \"%s\" not recognised.\n",label.c_str());
		snail(infile);
	}
	RTRV function_get_int(std::ifstream &infile, std::string &label)
	{
		std::string word="";
		get_single(infile,word);
		int value=atoi(word.c_str());
		GENERIC* ptr=data_map[label];
		(*(static_cast<int*>(ptr)))=value;
		if(coutput)printf("%s = %d\n",label.c_str(),value);
	}
	RTRV function_get_double(std::ifstream &infile, std::string &label)
	{
		std::string word="";
		get_single(infile,word);
		double value;
		if(word=="1.#INF") value = numeric_limits<double>::infinity();
		else if(word=="-1.#IND") value = numeric_limits<double>::quiet_NaN();
		else if(word=="-1.#INF") value = numeric_limits<double>::signaling_NaN();
		else value=atof(word.c_str());
		GENERIC* ptr=data_map[label];
		(*(static_cast<double*>(ptr)))=value;
		if(coutput)printf("%s = %g\n",label.c_str(),value);
	}
	RTRV function_get_string(std::ifstream &infile, std::string &label)
	/* TP_STRING must be enclosed by "double quotes" and must fit on one line */
	{
		int character;
		std::string word="";
		bool line_delim_found=false;
		bool include_char=true;
		bool string_terminated=false;

		while(eof(infile)==false)
		{
			character=get(infile);

			include_char=true;
			int i;
			for(i=0;i<(int)white_space.size();i++)
			{
				if(character==white_space[i])include_char=false;
			}
			for(i=0;i<(int)line_delimiters.size();i++)
			{
				if(character==line_delimiters[i])
				{
					include_char=false;
					line_delim_found=true;
				}
			}
			if(character=='\"') {
				include_char = false;
				break;
			}

			if(include_char==true)
			{
				word += static_cast<char>(character);
				break;
			}
			if(line_delim_found==true)break;
		}

		if(line_delim_found==false)
		{
			while(eof(infile)==false)
			{
				character=get(infile);

				include_char=true;
				int i;
				if(character=='\"') {
					include_char = false;
					string_terminated = true;
				}
				for(i=0;i<(int)line_delimiters.size();i++)
				{
					if(character==line_delimiters[i])
					{
						include_char=false;
						line_delim_found=true;
					}
				}

				if(include_char==true)word += static_cast<char>(character);
				if(line_delim_found==true)break;
				if(string_terminated==true) {
					snail(infile);
					break;
				}
			}
		}
		
		GENERIC* ptr=data_map[label];
		(*(static_cast<std::string*>(ptr)))=word;
		if(coutput)printf("%s = %s\n",label.c_str(),word.c_str());
	}
	/*	RTRV function_get_string(std::ifstream &infile, std::string &label)
	{
		int character;
		std::string word="";
		bool line_delim_found=false;
		bool include_char=true;
		bool string_started=false;

		while(eof(infile)==false)
		{
			character=get(infile);

			include_char=true;
			int i;
			for(i=0;i<white_space.size();i++)
			{
				if(character==white_space[i])include_char=false;
			}
			for(i=0;i<line_delimiters.size();i++)
			{
				if(character==line_delimiters[i])
				{
					include_char=false;
					line_delim_found=true;
				}
			}

			if(include_char==true)
			{
				word += static_cast<char>(character);
				break;
			}
			if(line_delim_found==true)break;
		}

		if(line_delim_found==false)
		{
			while(eof(infile)==false)
			{
				character=get(infile);

				include_char=true;
				int i;
				for(i=0;i<line_delimiters.size();i++)
				{
					if(character==line_delimiters[i])
					{
						include_char=false;
						line_delim_found=true;
					}
				}

				if(include_char==true)word += static_cast<char>(character);
				if(line_delim_found==true)break;
			}
		}
		
		GENERIC* ptr=data_map[label];
		(*(static_cast<std::string*>(ptr)))=word;
		if(coutput)printf("%s = %s\n",label.c_str(),word.c_str());
	}*/
	RTRV function_get_vector_int(std::ifstream &infile, std::string &label)
	{
		std::vector<std::string> words;
		get_multiple(infile,words);

		GENERIC* g_ptr=data_map[label];
		std::vector<int>* ptr=static_cast<std::vector<int>*>(g_ptr);
		ptr->clear();
		int i;
		for(i=0;i<(int)words.size();i++)
			ptr->push_back(atoi(words[i].c_str()));
		if(coutput)printf("%s read in %lu elements\n",label.c_str(),ptr->size());
	}
	RTRV function_get_vector_double(std::ifstream &infile, std::string &label)
	{
		std::vector<std::string> words;
		get_multiple(infile,words);

		GENERIC* g_ptr=data_map[label];
		std::vector<double>* ptr=static_cast<std::vector<double>*>(g_ptr);
		ptr->clear();
		int i;
		for(i=0;i<(int)words.size();i++)
			ptr->push_back(atof(words[i].c_str()));
		if(coutput)printf("%s read in %lu elements\n",label.c_str(),ptr->size());
	}
	RTRV function_get_external_vector_double(std::ifstream &infile, std::string &label)
	{
		std::string filename="";
		std::string internal_call="Opening external file";
		data_map[internal_call]=&filename;
		function_get_string(infile,internal_call);

		std::ifstream extfile(filename.c_str());
		if(extfile.is_open()==false)error("function_get_external_vector_double():External file not found");

		if(coutput)printf("\t");
		function_get_vector_double(extfile,label);	
	}
	RTRV function_get_vector_string(std::ifstream &infile, std::string &label)
	{
		std::vector<std::string> words;
		get_multiple(infile,words);

		GENERIC* g_ptr=data_map[label];
		std::vector<string>* ptr=static_cast<std::vector<string>*>(g_ptr);
		ptr->clear();
		int i;
		for(i=0;i<(int)words.size();i++)
			ptr->push_back(words[i]);
		if(coutput)printf("%s read in %lu elements\n",label.c_str(),ptr->size());
	}
};
};

#endif
