/*  Copyright 2013 Daniel Wilson and Xavier Didelot.
 *
 *  xmfa.h
 *  Part of ClonalFrameML
 *
 *  ClonalFrameML is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  ClonalFrameML is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with ClonalFrameML. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <DNA.h>

void readXMFA(const char *filename,DNA * dna) {
		string unlink="";
		ifstream in(filename);
		if(!in.is_open()) {
			string errmsg = "readXMFA(): File "+string(filename)+" not found";
			error(errmsg.c_str());
		}
		
		dna->nseq = 0;
		int block=0;
		string s;
		getline(in,s);
		if(s.length()>0 && s[0]!='>') {
			string errmsg = "readXMFA(): File "+string(filename)+" did not begin with '>'";
			error(errmsg.c_str());
		}
		dna->label.push_back(s.substr(1));
		string newseq = "";
		while(!in.eof()) {
			getline(in,s);
			if(s.length()>0 && (s[0]=='>'||s[0]=='=')) {
				if (block==0) dna->sequence.push_back("");
				if (dna->nseq>=0) dna->sequence[dna->nseq]+=unlink+newseq;
				newseq = "";
				if(s[0]=='>') {dna->nseq++;if (block==0) dna->label.push_back(s.substr(1));} 
				else {block++;dna->nseq=-1;if (block==1) unlink=string(1000,'N');}
				} else newseq += s;
		}
		dna->nseq=dna->sequence.size();
		dna->lseq=dna->sequence[0].length();
		in.close();
}

