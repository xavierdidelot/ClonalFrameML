/*
 *  newick.h
 *  newick
 *
 *  Created by Daniel Wilson on 05/03/2013.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _NEWICK_H_
#define _NEWICK_H_
#include <vector>
#include <string>
#include <myerror.h>
#include <sstream>
#include <iostream>

using std::vector;
using std::string;
using myutils::error;
using std::stringstream;
using std::endl;
using std::cout;
using myutils::warning;

namespace myutils
{
class NewickNode {
public:
	// Member variables
	NewickNode *anc;				// Ancestral node
	vector<NewickNode*> dec;		// Descendant nodes (any number)
	double len;						// Length
	string str;						// Name
	vector<NewickNode*> *allnodes;	// Pointer to all nodes in the tree
	
	// Member functions
	NewickNode() {
		initialize();
	}
	NewickNode(string token, NewickNode *anc_in) {
		initialize();
		anc = anc_in;
		if(!anc==0) {
			// Get pointer to allnodes
			allnodes = anc->allnodes;
			// Add self to list of descendants
			allnodes->push_back(this);
		}
		// Remember when a node is created to add it to one's descendants
		process_token(token);
	}
	void initialize(){
		anc = 0;
		dec = vector<NewickNode*>(0);
		len = 0.0;
		str = "";
		allnodes = 0;
	}
	void process_token(string token){
		// If this is part of a nexus file, assume all comments enclosed by square brackets have been removed
		// If this is the outermost node, assume the trailing semi-colon has already been removed
		// Locate left-most open bracket
		size_t lbrkt = token.find('(');
		// Locate right-most close bracket
		size_t rbrkt = token.rfind(')');
		// Locate right-most colon
		size_t rcoln = token.rfind(':');
		// Some checks
		if(lbrkt!=token.npos && rbrkt!=token.npos && lbrkt>rbrkt) {
			stringstream errTxt;
			errTxt << "Token: " << token << endl;
			errTxt << "Left bracket to right of right bracket: " << lbrkt << ", " << rbrkt;
			error(errTxt.str().c_str());
		}
		if(lbrkt==token.npos && rbrkt!=token.npos) {
			stringstream errTxt;
			errTxt << "Token: " << token << endl;
			errTxt << "Found right bracket but no left bracket";
			error(errTxt.str().c_str());
		}
		if(lbrkt!=token.npos && rbrkt==token.npos) {
			stringstream errTxt;
			errTxt << "Token: " << token << endl;
			errTxt << "Found left bracket but no right bracket";
			error(errTxt.str().c_str());
		}
		if(rbrkt==lbrkt+1) {
			stringstream errTxt;
			errTxt << "Token: " << token << endl;
			errTxt << "Empty brackets";
			error(errTxt.str().c_str());
		}
		// Some indicator variables
		// Has descendants within brackets
		bool has_brkt = (lbrkt!=token.npos);
		// Has a colon
		bool has_coln = (rcoln!=token.npos && (!has_brkt || rcoln>rbrkt));
		if(has_coln && has_brkt) {
			// Name the node
			if(rcoln>rbrkt+1) {
				str = token.substr(rbrkt+1,rcoln-rbrkt-1);
			} else {
				str = "";
			}
			// Get the length
			if(rcoln<token.length()-1) {
				len = atof(token.substr(rcoln+1,token.length()-1-rcoln).c_str());
			} else {
				len = 0.0;
			}
		} else if(has_coln && !has_brkt) {
			// Name the node
			if(rcoln>0) {
				str = token.substr(0,rcoln);
			} else {
				str = "";
			}
			// Get the length
			if(rcoln<token.length()-1) {
				len = atof(token.substr(rcoln+1,token.length()-1-rcoln).c_str());
			} else {
				len = 0.0;
			}
		} else if(!has_coln && has_brkt) {
			// Name the node
			if(rbrkt<token.length()-1) {
				str = token.substr(rbrkt+1,token.length()-1-rbrkt);
			} else {
				str = "";
			}
			// There is no length
			len = 0.0;
		} else if(!has_coln && !has_brkt) {
			// Name the node
			str = token;
			// There is no length
			len = 0.0;
		}
		// Deal with descendant nodes
		if(has_brkt) {
			string desc = token.substr(lbrkt+1,rbrkt-lbrkt-1);
			// Find all top-level commas
			vector<size_t> poscomma(0);
			size_t pos;
			// Keep track of the opening and closing of brackets within the string
			int nlbrkt = 0;
			int nrbrkt = 0;
			for(pos=0;pos<desc.length();pos++) {
				if(desc[pos]==',') {
					if(nlbrkt==nrbrkt) {
						poscomma.push_back(pos);
					}
				} else if(desc[pos]=='(') {
					++nlbrkt;
				} else if(desc[pos]==')') {
					++nrbrkt;
					if(nrbrkt>nlbrkt) {
						stringstream errTxt;
						errTxt << "Token: " << desc << endl;
						errTxt << "Found right bracket before left bracket";
						error(errTxt.str().c_str());
					}
				}
			}
			if(nlbrkt!=nrbrkt) {
				stringstream errTxt;
				errTxt << "Token: " << desc << endl;
				errTxt << "Too few right brackets";
				error(errTxt.str().c_str());
			}
			// For each descendant separated by commas, start a new node
			if(poscomma.size()==0) {
				stringstream errTxt;
				errTxt << "Token: " << desc << endl;
				errTxt << "Single descendant found";
				warning(errTxt.str().c_str());
				dec.push_back(new NewickNode(desc.substr(0,desc.length()),this));
			} else{
				dec.push_back(new NewickNode(desc.substr(0,poscomma[0]),this));
				int i;
				for(i=1;i<poscomma.size();i++) {
					dec.push_back(new NewickNode(desc.substr(poscomma[i-1]+1,poscomma[i]-1-poscomma[i-1]),this));
				}
				dec.push_back(new NewickNode(desc.substr(poscomma[poscomma.size()-1]+1,desc.length()-1-poscomma[poscomma.size()-1]),this));
			}
		}
		if(false) {
			cout << "Node processed:" << endl;
			cout << "address: " << this << endl;
			cout << "anc: " << anc << endl;
			cout << "str: " << str << endl;
			cout << "len: " << len << endl;
			cout << "descendant addresses:";
			int i;
			for(i=0;i<dec.size();i++) cout << " " << dec[i];
			cout << endl;
			cout << "allnodes address: " << allnodes << endl;
		}
	}
};

class NewickTree {
public:
	// Member variables
	NewickNode root;
	vector<NewickNode*> allnodes;  // Pointer to all nodes in the tree
	
	// Member functions
	NewickTree() {
	}
	NewickTree(string token) {
		process_token(token);
	}
	void process_token(string token) {
		// Check for a trailing semi-colon
		if(token[token.length()-1]!=';') {
			stringstream errTxt;
			errTxt << "Token: " << token << endl;
			errTxt << "Expected trailing semi-colon but none found";
			error(errTxt.str().c_str());
		}
		// Set member variables
		allnodes = vector< NewickNode* >(1,&root);
		root.allnodes = &allnodes;
		// Start from the root node, having removed the trailing semi-colon
		root.process_token(token.substr(0,token.length()-1));
	}
};
	
}; // namespace myutils

#endif // _NEWICK_H_
