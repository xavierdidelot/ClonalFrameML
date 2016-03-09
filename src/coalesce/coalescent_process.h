/*  Copyright 2013 Daniel Wilson.
 *
 *  coalescent_process.h
 *  Part of the coalesce library.
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
#ifndef _COALESCENT_PROCESS_H_
#define _COALESCENT_PROCESS_H_


#include <sstream>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;

#include "myutils/matrix.h"
using myutils::Matrix;
#include "myutils/random.h"
using myutils::Random;
#include "myutils/myerror.h"
using myutils::error;
#include "coalesce/coalescent_record.h"
#include "coalesce/coalescent_control.h"
#include "coalesce/mutation.h"

class ptr_vector
{
	int size;
public:
	class mt_node **ptr;					//vec of generic ptrs

public:
	ptr_vector() {};
	ptr_vector& initialize(const int size_in)
	{
		size=size_in;
		ptr=(class mt_node**) malloc((size_t) size*sizeof(class mt_node*));
		int i;
		for(i=0;i<size;i++)
			ptr[i]=NULL;
		return *this;
	}
	ptr_vector& nullify(const int from, const int to)
	{
		int i;
		for(i=from;i<=to;i++)
			ptr[i]=NULL;
		return *this;
	}
	ptr_vector& nullify(const int position)
	{
		ptr[position]=NULL;
		return *this;
	}
	ptr_vector& nullify()
	{
		int i;
		for(i=0;i<size;i++)
			ptr[i]=NULL;
		return *this;
	}
	ptr_vector& copy(ptr_vector *donor, const int from, const int to)
	{
		int i;
		for(i=from;i<=to;i++)
			ptr[i]=donor->ptr[i];
		return *this;
	}
	ptr_vector& copy(ptr_vector *donor, const int position)
	{
		ptr[position]=donor->ptr[position];
		return *this;
	}
	ptr_vector& copy(ptr_vector *donor)
	{
		int i;
		for(i=0;i<size;i++)
			ptr[i]=donor->ptr[i];
		return *this;
	}
	inline ptr_vector& assign(class mt_node *target, const int position)
	{
		ptr[position]=target;
		return *this;
	}
	inline ptr_vector& assign(class mt_node *target, const int from, const int to)
	{
		int i;
		for(i=from;i<=to;i++)
			ptr[i]=target;
		return *this;
	}

	int get_size() {return size;};

	~ptr_vector()
	{
		//nullify();
		free((class mt_node**) ptr);
	}
};

class ap_node
{
public:
	/*Fixed once*/
	int id;

	/*Recyclable*/
	enum {NOT_IN_USE,IN_USE,FIXED_NODE} flag;
	int deme;			// records the deme the node belongs to
	double time;
	class ptr_vector AMP;
	double rlen;
	double L;
	int ltr;
	int rtr;
	int active_id;		// records the position in active_node
	int ctree_id;		// records the position in the conditional marginal tree
public:
	ap_node() {};
	ap_node& initialize(const int id_in, const int size)
	{
		id=id_in;
		active_id = ctree_id = -1;
		AMP.initialize(size);
		recycle();
		return *this;
	}
	ap_node& recycle()
	{
		flag = NOT_IN_USE;
		active_id = ctree_id = deme = -1;
		//time=0.0;
		//AMP.nullify();
		return *this;
	}
	ap_node& activate(double *time_in)
	{
		flag = IN_USE;
		time=*time_in;
		return *this;
	}
	
	~ap_node() {};
};

class eventChain {
public:
	enum eventType {NONE,COALESCENCE,RECOMBINATION,ADD_LINEAGE,END};
protected:
	class eventChainEvent {
	public:
		double k,rlen,time;
		eventType type;
		eventChainEvent() {
			k = rlen = time = 0.0;
			type = NONE;
		}
	};
	vector<eventChainEvent> ev;
public:
	double rho;
public:
	eventChain() {
		ev = vector<eventChainEvent>(0);
	}
	eventChain(const int size_in) {
		ev = vector<eventChainEvent>(size_in);
	}
	const int size() const {
		return (int)ev.size();
	}
	void resize(const int size_in) {
//		cout << "resize to " << size_in << endl;
		if(size_in<0) myutils::error("eventChain::resize(): cannot have a negative size");
		ev.resize(size_in);
//		cout << "done resizing" << endl;
	}
	eventChainEvent& operator[](const int pos) {
		return ev[pos];
	}
	double loglikelihood(const double rh) {
		if(ev.size()<=0) myutils::error("eventChain::loglikelihood(): chain has zero length");
		if(ev[0].type==NONE) myutils::error("eventChain::loglikelihood(): no chain exists");
		int e = 0;
		if(ev[e].type==END) myutils::error("eventChain::loglikelihood(): first event is the last");
		double rate = (ev[e].k*ev[e].rlen*rh + ev[e].k*(ev[e].k-1.0))/2.0;
		double prec = ev[e].k*ev[e].rlen*rh/2.0/rate;
		double L = 0.0;
		if(ev[e].type==ADD_LINEAGE) {
	       		L += - rate * (ev[e].time);
		}
		else if(ev[e].type==RECOMBINATION) {
			L += log(rate) - rate * (ev[e].time) + log(prec);
		}
		else if(ev[e].type==COALESCENCE) {
			L += log(rate) - rate * (ev[e].time) + log(1.-prec);
		}
		for(e=1;e<(int)ev.size();e++) {
			if(ev[e].type==END) break;
			rate = (ev[e].k*ev[e].rlen*rh + ev[e].k*(ev[e].k-1.0))/2.0;
			prec = ev[e].k*ev[e].rlen*rh/2.0/rate;
			if(ev[e].type==ADD_LINEAGE) {
				L += - rate * (ev[e].time - ev[e-1].time);
			}
			else if(ev[e].type==RECOMBINATION) {
				L += log(rate) - rate * (ev[e].time - ev[e-1].time) + log(prec);
			}
			else if(ev[e].type==COALESCENCE) {
				L += log(rate) - rate * (ev[e].time - ev[e-1].time) + log(1.-prec);
			}
			
		}
		if(e==ev.size()) error("eventChain::loglikelihood() chain has no end");
		return L;
	}
};

class coalescent
{
public:
	/*Fixed*/
	class Control *con;			//ptr to con
	class Random *ran;			//ptr to ran
	class marginal_tree *tree;	//vec of tree's
	class ap_node **node;		//vec of ptrs to ap_node's
	class ap_node **active_node;//vec of ptrs to ap_node's
	class ap_node **inactive_node;//vec of ptrs to ap_node's
	int nodes_reserved;
	int L;
	int **segregating_tree;
	int **internal_seg_tree;
	int seg_tree_id;
	Matrix<double> genotype;
	bool no_gene_conversion;

	/*Recyclable*/
	int n_inactive;
	int ARG_k;					//#lineages in ARG
	int gen;					//number of events
	double total_rlen;
	double rho;					//total_rlen/ARG_k
	int n_segregating;
		
	bool samples_waiting;
	double time_next_sample;			// contemporaneous samples
	vector<double>::iterator ntimes_itr;// iterator for rifling through con->ntimes
	int next_waiting_sample;

	int nrecTypeI;
	int nrecTypeII;
	int nrecTypeIII;
	
	int nco,nrec,naddbase,nmut;
	vector<int> nrecWatt;

	int ncoI,ncoIIa,ncoIIb,ncoIII;
	
	/* Variables for conditional simulation */
	int ARG_k_fixed;
	vector<double> ftimes;	// times for fixed events
	vector<ap_node*> fnode;

	/* Variables for structured coalescent */
	ap_node ***ptr_deme;				// con->ndemes * con->nsamp matrix: members of each deme
	Vector<int> k_deme;					// number of ancestral lineages in each deme
	Vector<double> rho_deme;			// effective recn rate of each deme
	Vector<double> sum_mig;				// total backwards-in-time mig for each deme i
	double rate_coal, rate_recn, rate_mign;
	Vector<double> coal_deme;

protected:
	vector<int> _uniqueHaps;
	vector<int> _sites;
	LowerTriangularMatrix<int> ____B;
	vector<int> _M;
	vector<double> _F;
	vector<double> _four;
	LowerTriangularMatrix< vector<double> > _G;
	LowerTriangularMatrix<double> _A,___B,___C;
	Matrix<double> _D;

public:
	coalescent() {};
	coalescent& initialize(class Control *con_in, class Random *ran_in)
	{
		con=con_in;
		ran=ran_in;
		if(con->nsamp<0) con->nsamp = 0;
		if(con->ntimes.size()==0) con->ntimes = vector<double>(con->nsamp,0.0);
		if(con->seq_len<0) con->seq_len = 0;
		if(con->Negens<0) con->Negens = 1;
		if(con->len.size()==0) con->len = vector<int>(1,con->seq_len);
		if(con->r<0.0) con->r = 0.0;
		if(con->lambda<0.0) con->lambda = 0.0;

		if(con->lambda==0.0)no_gene_conversion=true;
		else no_gene_conversion=false;

		L=con->seq_len;
		tree=(class marginal_tree*) malloc((size_t) L*sizeof(marginal_tree));
		int i;
		for(i=0;i<L;i++)
			tree[i].initialize(i,con->nsamp);

		internal_seg_tree=(int**) malloc((size_t) 2*sizeof(int*));
		internal_seg_tree[0]=(int*) malloc((size_t) L*sizeof(int));
		internal_seg_tree[1]=(int*) malloc((size_t) L*sizeof(int));
		seg_tree_id=0;
		segregating_tree=&(internal_seg_tree[seg_tree_id]);
		
		ARG_k=ARG_k_fixed=0;
		node=(class ap_node**) malloc((size_t) ARG_k*sizeof(ap_node*));
		active_node=(class ap_node**) malloc((size_t) ARG_k*sizeof(ap_node*));
		inactive_node=(class ap_node**) malloc((size_t) ARG_k*sizeof(ap_node*));
		if(con->ndemes>0) {
			ptr_deme = (ap_node***) malloc((size_t) con->ndemes*sizeof(ap_node**));
			for(i=0;i<con->ndemes;i++) ptr_deme[i] = (ap_node**) malloc((size_t) ARG_k*sizeof(ap_node*));
		}
		n_inactive=0;

		nodes_reserved=0;
		reserve_nodes(10*con->nsamp);

		genotype.initialize(con->nsamp,L);

		return *this;
	}
	coalescent& go()
	{
		recycle();
		if(add_next_sample()!=0.0)error("Most recent node does not occur at time zero");
		double current_time=0.0;
		gen=0;

		while((ARG_k>1)||(samples_waiting))
		{
			//event(&current_time);
			double denom = 1.0/((double)ARG_k*rho+(double)ARG_k*((double)ARG_k-1.0));
			current_time += constant_size_model(2.0*denom);

			if((samples_waiting)&&(current_time>=time_next_sample))
			{
				current_time = add_next_sample();
			}
			else
			{
				/*2nd, choose type of event*/
				double rnum1 = ran->U();
				double pr_recom=(double)ARG_k*rho*denom;
				if (rnum1 <= pr_recom) recombine(&current_time);
				else coalesce(&current_time);
			}

			++gen;
		}

		return *this;
	}
	coalescent& go(eventChain& e)
	{
		recycle();
		if(add_next_sample()!=0.0)error("Most recent node does not occur at time zero");
		double current_time=0.0;
		gen=0;
		e.rho = 2.0*con->r;
		if(e.size()<1000) e.resize(1000);

		while((ARG_k>1)||(samples_waiting))
		{
			if(e.size()<=gen) e.resize(2*e.size());
			double denom = 1.0/((double)ARG_k*rho+(double)ARG_k*((double)ARG_k-1.0));
			current_time += constant_size_model(2.0*denom);
			e[gen].k = (double)ARG_k;
			e[gen].rlen = rho/e.rho;
			e[gen].time = current_time;

			if((samples_waiting)&&(current_time>=time_next_sample))
			{
				current_time = add_next_sample();
				e[gen].time = current_time;
				e[gen].type = eventChain::ADD_LINEAGE;
			}
			else
			{
				/*2nd, choose type of event*/
				double rnum1 = ran->U();
				double pr_recom=(double)ARG_k*rho*denom;
				if (rnum1 <= pr_recom) {
					recombine(&current_time);
					e[gen].type = eventChain::RECOMBINATION;
				}
				else {
					coalesce(&current_time);
					e[gen].type = eventChain::COALESCENCE;
				}
			}

			++gen;
		}

		e[gen-1].type = eventChain::END;
		return *this;
	}
	coalescent& migrate()
	{
		recycle();
		if(add_next_sample()!=0.0)error("Most recent node does not occur at time zero");
		double current_time=0.0;
		gen=0;
		int i,j;

		//const char tab = '\t';
		//for(i=0;i<con->ndemes;i++) cout << tab << "k" << i << tab << "co" << tab << "mig";
		//cout << endl;
		//cout << setprecision(3);
		//Vector<int> k_avg(con->ndemes,0);
		while((ARG_k>1)||(samples_waiting))
		{
			//if(gen%2==0) {
				rate_coal = rate_recn = rate_mign = 0.0;
				for(i=0;i<con->ndemes;i++) {
					coal_deme[i] = (double)k_deme[i] * (double)(k_deme[i]-1) / con->N_deme_over_D[i];
					rate_coal += coal_deme[i];
					rate_recn += rho_deme[i];
					rate_mign += (double)k_deme[i] * sum_mig[i];
				}
				rate_coal /= 2.0;
				rate_recn /= 2.0;
				rate_mign /= 2.0;
			//}
			//cout << current_time;
			//for(i=0;i<con->ndemes;i++) {
			//	cout << tab << k_deme[i] << tab << coal_deme[i]/2. << tab << (double)k_deme[i] * sum_mig[i]/2.;
			//	k_avg[i] += k_deme[i];
			//}
			//cout << endl;
			double denom = (rate_coal + rate_recn + rate_mign);
			current_time += constant_size_model(1.0/denom);

			if((samples_waiting)&&(current_time>=time_next_sample))
			{
				current_time = add_next_sample();
			}
			else
			{
				/*2nd, choose type of event*/
				double rnum1 = ran->U() * denom;
				if(rnum1 <= rate_coal) migrate_coalesce(&current_time);
				else if(rnum1 <= rate_coal+rate_recn) migrate_recombine(&current_time);
				else migrate_migrate(&current_time);
			}

			for(i=0;i<con->ndemes;i++)
				for(j=0;j<k_deme[i];j++)
					if(ptr_deme[i][j]->flag==ap_node::NOT_IN_USE) error("migrate(): pointer problem");
			++gen;
		}

		//cout << current_time;
		//for(i=0;i<con->ndemes;i++) cout << tab << (double)k_avg[i]/(double)gen;
		//cout << endl;
		return *this;
	}
	coalescent& conditional(class marginal_tree &ctree)
	{
		recycle();
		int i;
		fnode = vector<ap_node*>(ctree.size,(ap_node*)NULL);
		ftimes = vector<double>(ctree.size,0.0);
		for(i=0;i<ctree.size;i++) ftimes[i] = ctree.node[i].time;
		ntimes_itr = ftimes.begin();
		if(add_conditional_event(ctree)!=0.0)error("Most recent node does not occur at time zero");
		double current_time=0.0;
		gen=0;

		while((ARG_k>1)||(samples_waiting))
		{
			//event(&current_time);
			//double denom = (double)ARG_k*rho+2.*(double)(ARG_k_fixed)*(double)(ARG_k-ARG_k_fixed);
			//if(ARG_k>ARG_k_fixed) denom += (double)(ARG_k-ARG_k_fixed)*((double)(ARG_k-ARG_k_fixed)-1.0);
			double denom = (double)ARG_k*rho+(double)(ARG_k)*(double)(ARG_k-1);
			if(ARG_k_fixed>1) denom -= (double)(ARG_k_fixed)*(double)(ARG_k_fixed-1);
			if(denom == 0.0) {
				if(samples_waiting) current_time = add_conditional_event(ctree);
				else error("conditional(): infinite time until next event");
			}
			else {
				denom = 1.0/denom;
				current_time += constant_size_model(2.0*denom);

				if((samples_waiting)&&(current_time>=time_next_sample))
				{
					current_time = add_conditional_event(ctree);
				}
				else
				{
					/*2nd, choose type of event*/
					double rnum1 = ran->U();
					double pr_recom=(double)ARG_k*rho*denom;
					if (rnum1 <= pr_recom) {
						conditionally_recombine(&current_time);
						/*if(tree[0].node[2].time!=0) {
							warning("MRCA is not supposed to be found during recombination");
						}*/
					}
					else {
						conditionally_coalesce(&current_time);
					}
				}
			}

			if(ARG_k_fixed==1) {
				warning("Single fixed lineage left");
			}
			if(ARG_k_fixed>0)
				for(i=0;i<(int)fnode.size();i++) if(fnode[i]!=NULL) if(fnode[i]->active_id==-1) {
					warning("fnode inconsistency");
				}
			if(ARG_k_fixed>ARG_k) error("conditional(): ARG_k_fixed > ARG_k");
			for(i=0;i<ARG_k;i++) if(active_node[i]->active_id!=i) error("conditional(): active_id's incorrect");
			int ctr_fnode = 0;
			for(i=0;i<ctree.size;i++) if(fnode[i]!=NULL) ++ctr_fnode;
			//if(ctr_fnode!=ARG_k_fixed) error("conditional(): ctr_fnode <> ARG_k_fixed");
			++gen;
		}
		return *this;
	}
	coalescent& mutate()
	{
		int i;
		for(i=0;i<L;i++)
		{
			mutate_tree(i,tree[i].size-1,mutate_mrca());
		}
		return *this;
	}
	coalescent& mutate(Mutation_Matrix *M) {
		int site;
		for(site=0;site<L;site++) mutate_tree(site,tree[site].size-1,M->draw(),M);
		return *this;
	}
	coalescent& mutate(const int site, Mutation_Matrix *M) {
		mutate_tree(site,tree[site].size-1,M->draw(),M);
		return *this;
	}
	coalescent& mutate(const int site, Mutation_Matrix *M, vector<int>& mutLog) {
		mutLog.clear();
		mutate_tree_and_record(site,tree[site].size-1,M->draw(),M,mutLog);
		return *this;
	}
	coalescent& output_FASTA(vector<char> &code, const char* filename)
	{
		FILE* fout=fopen(filename,"w");
//		fprintf(fout,"%d %d\n\n",con->nsamp,con->seq_len);

		int n;
		for(n=0;n<con->nsamp;n++)
		{
			fprintf(fout,">seq%d_%g\n",n,con->ntimes[n]*con->Negens);
			int pos;
			for(pos=0;pos<con->seq_len;pos++)
				fprintf(fout,"%c",code[(int)genotype[n][pos]]);
			fprintf(fout,"\n");
		}
		fclose(fout);
		return *this;
	}
	coalescent& output_FASTA(vector<string> code, const char* filename)
	{
//		FILE* fout=fopen(filename,"w");
//		fprintf(fout,"%d %d\n\n",con->nsamp,con->seq_len);
		ofstream fout(filename);
//		fout << con->nsamp << " " << con->seq_len << endl << endl;

		int n;
		for(n=0;n<con->nsamp;n++)
		{
			fout << ">seq" << n << "_" << con->ntimes[n]*con->Negens << endl;
			//fprintf(fout,">seq%d_%g\n",n,con->ntimes[n]*con->Negens);
			int pos;
			for(pos=0;pos<con->seq_len;pos++)
				fout << code[(int)genotype[n][pos]];
				//fprintf(fout,"%s",code[genotype[n][pos]].c_str());
			fout << endl;
			//fprintf(fout,"\n");
		}
		//fclose(fout);
		fout.close();
		return *this;
	}
	/* which is a vector true or false whether to include each sequence */
	coalescent& output_FASTA(vector<char> &code, const char* filename, vector<bool> &which)
	{
		if(which.size()!=con->nsamp) error("coalescent::output_FASTA(): which must have length nsamp");
		FILE* fout=fopen(filename,"w");

		int n;
		for(n=0;n<con->nsamp;n++)
		{
			if(which[n]) {
				fprintf(fout,">seq%d_%g\n",n,con->ntimes[n]*con->Negens);
				int pos;
				for(pos=0;pos<con->seq_len;pos++)
					fprintf(fout,"%c",code[(int)genotype[n][pos]]);
				fprintf(fout,"\n");
			}
		}
		fclose(fout);
		return *this;
	}
	/* which is a vector true or false whether to include each sequence */
	coalescent& output_FASTA(vector<string> &code, const char* filename, vector<bool> &which)
	{
		if(which.size()!=con->nsamp) error("coalescent::output_FASTA(): which must have length nsamp");
		ofstream fout(filename);

		int n;
		for(n=0;n<con->nsamp;n++)
		{
			if(which[n]) {
				fout << ">seq" << n << "_" << con->ntimes[n]*con->Negens << endl;
				int pos;
				for(pos=0;pos<con->seq_len;pos++)
					fout << code[(int)genotype[n][pos]];
				fout << endl;
			}
		}
		fout.close();
		return *this;
	}
	coalescent& output_MEP(const char* filename)
	{
		FILE* fout=fopen(filename,"w");
		fprintf(fout,"n=%d, mu=%g, r=%g, Negens=%g\n",con->nsamp,con->M/con->Negens,con->r/con->Negens,con->Negens);
		fprintf(fout,"Time points = ");
		int i;
		for(i=0;i<(int)con->ntimes.size();i++)fprintf(fout,"%g ",con->ntimes[i]*con->Negens);
		fprintf(fout,"\n\n");
		int mrca=2*con->nsamp-2;
		double t_height=0.0;
		for(i=0;i<con->seq_len;i++)
		{
			double temp=tree[i].node[mrca].time;
			if(t_height!=tree[i].node[mrca].time)
			{
				t_height=tree[i].node[mrca].time;
				fprintf(fout,"Position %d\tHeight %g\t\t%g\n",i,t_height,t_height*con->Negens);
			}
		}

		fclose(fout);
		return *this;
	}
	coalescent& output_tree(const int site)
	{
		/*This node is always the mrca*/
		int mrca=2*(con->nsamp-1);
		int i=site;

		/*Create names for the files*/
		stringstream ageout_file;
		ageout_file << "age" << i << ".dat";
		stringstream treeout_file;
		treeout_file << "tree" << i << ".dat";
		/*Open them for writing*/
		FILE *ageout = fopen(ageout_file.str().c_str(), "w");
		FILE *treeout = fopen(treeout_file.str().c_str(), "w");

		int tree_id=site;
		/*ageout contains ages for each of the nodes*/
		int j;
		for(j=0;j<=mrca;j++)
		{
			fprintf(ageout,"%d %g\n",mrca-j,10.0*tree[tree_id].node[j].time);
		}
		/*treeout contains the labels for each of the base nodes*/
		/*followed by a colon then a list of all nodes ancestral to it*/
		for(j=0;j<con->nsamp;j++)
		{
			int gt=(int)genotype[j][tree_id];
			fprintf(treeout,"%2d : %d ",gt+1,mrca-j);
			class mt_node* anc;
			class mt_node* nextanc=tree[tree_id].node[j].ancestor;
			do
			{
				anc=nextanc;
				fprintf(treeout,"%d ",mrca-anc->id);
				nextanc=anc->ancestor;
			}while(nextanc!=NULL);
			fprintf(treeout,"\n");
		}

		fclose(ageout);
		fclose(treeout);


		FILE *tpicout = fopen("tpic.bat","w");
		int number=1;
		for(i=0;i<number;i++)
			fprintf(tpicout,"treepic tree%d.dat tree%d.ps -a age%d.dat -nv -u %g 0.5 %%%%.2f\n",i,i,i,tree[i].node[mrca].time);
		fclose(tpicout);
		FILE *viewtpicout = fopen("viewtpic.bat","w");
		for(i=0;i<number;i++)
			fprintf(viewtpicout,"tree%d.ps\n",i);
		fclose(viewtpicout);

		return *this;
	}
	~coalescent()
	{
		free((class marginal_tree*) tree);
		segregating_tree=0;
		free((int*) internal_seg_tree[1]);
		free((int*) internal_seg_tree[0]);
		free((int**) internal_seg_tree);
		free((class ap_node**) active_node);
		free((class ap_node**) inactive_node);
		int i;
		if(con->ndemes>0) {
			for(i=con->ndemes-1;i>=0;i--) free((ap_node**) ptr_deme[i]);
			free((ap_node***) ptr_deme);
		}
		for(i=0;i<nodes_reserved;i++)
		{
			delete node[i];
			node[i]=0;
		}
		free((class ap_node**) node);
	}
	inline double* operator[](int pos){return genotype.element[pos];}

protected:
	coalescent& reserve_nodes(const int number)
	{
		if(number<=nodes_reserved) return *this;
		node=(class ap_node**) realloc(node,(size_t) number*sizeof(ap_node*));
		active_node=(class ap_node**) realloc(active_node,(size_t) number*sizeof(ap_node*));
		inactive_node=(class ap_node**) realloc(inactive_node,(size_t) number*sizeof(ap_node*));
		int i;
		if(con->ndemes>0) {
			for(i=0;i<con->ndemes;i++)
				ptr_deme[i] = (ap_node**) realloc(ptr_deme[i],(size_t) number*sizeof(ap_node*));
		}
		
		for(i=nodes_reserved;i<number;i++)
		{
			node[i]=new ap_node;
			node[i]->initialize(i,L);
			inactive_node[n_inactive]=node[i];
			++n_inactive;
		}

		nodes_reserved=number;
		return *this;
	}
	coalescent& recycle()
	{
		int i,j;
		for(i=0;i<ARG_k;i++)
			deactivate_node(i);
		for(i=0;i<L;i++)
		{
			tree[i].recycle();
			internal_seg_tree[0][i]=i;
			internal_seg_tree[1][i]=-1;
		}
		seg_tree_id=0;
		segregating_tree=&(internal_seg_tree[seg_tree_id]);
		n_segregating=L;
		ARG_k=ARG_k_fixed=0;
		n_inactive=nodes_reserved;

		samples_waiting=true;
		time_next_sample=0.0;
		ntimes_itr=con->ntimes.begin();
		next_waiting_sample=0;

		nrecTypeI=0;
		nrecTypeII=0;
		nrecTypeIII=0;

		nco=nrec=naddbase=nmut=0;
		nrecWatt = vector<int>(con->seq_len,0);
		ncoI=ncoIIa=ncoIIb=ncoIII=0;

		if(con->ndemes>0) {
			for(i=0;i<con->ndemes;i++)
				for(j=0;j<nodes_reserved;j++)
					ptr_deme[i][j] = NULL;
			k_deme = Vector<int>(con->ndemes,0);
			rho_deme = Vector<double>(con->ndemes,0.0);
			if(con->N_deme_over_D.size()!=con->ndemes) error("recycle(): con->N_deme_over_D wrong # demes");
			if(con->mig.nrows()!=con->ndemes) error("recycle(): con->mig wrong number of rows");
			if(con->mig.ncols()!=con->ndemes) error("recycle(): con->mig wrong number of columns");
			sum_mig = Vector<double>(con->ndemes,0.0);
			coal_deme = Vector<double>(con->ndemes,0.0);
			if(con->deme_config.size()!=con->nsamp) error("recycle(): con->deme_config wrong sample size");

			for(i=0;i<con->ndemes;i++)
				for(j=0;j<con->ndemes;j++)
					sum_mig[i] += con->mig[i][j];
		}

		if(!no_gene_conversion && con->lambda<=0.0) error("coalescent::recycle(): lambda<=0.0 in gene conversion model");

		return *this;
	}
	coalescent& event(double *time)
	{
		double denom = 1.0/((double)ARG_k*rho+(double)ARG_k*((double)ARG_k-1.0));
		(*time)+=constant_size_model(2.0*denom);

		if((samples_waiting)&&((*time)>=time_next_sample))
		{
			(*time) = add_next_sample();
		}
		else
		{
			/*2nd, choose type of event*/
			double rnum1 = ran->U();
			double pr_recom=(double)ARG_k*rho*denom;
			if (rnum1 <= pr_recom) recombine(time);
			else coalesce(time);
		}
			
		return *this;
	}
	coalescent& coalesce(double *time)
	{
		++nco;
		/*Create new lineage in the ARG*/
		class ap_node *new_node=create_node(time);
		//printf("           : ");int o;for(o=0;o<ARG_k;o++)printf("%d ",active_node[o]->id);printf("\n");
		/*NB If you deactivate the nodes before  */
		/*creating the new one then you overwrite*/
		/*memory you want to read from!          */
		/*However, must not allow the new node to*/
		/*be chosen as one of the coalescing     */
		/*nodes. Since it is added at the end of */
		/*the active_node vector, simply restrict*/
		/*the maximum node that can be chosen.   */

		/*Choose 1st lineage to coalesce*/
		int lin1=ran->discrete(0,ARG_k-2);			//New node is in position ARG_k-1, so do not
		class ap_node *ap_node1=active_node[lin1];	//allow this position to be chosen
		deactivate_node(lin1);
		//printf("           : ");for( o=0;o<ARG_k;o++)printf("%d ",active_node[o]->id);printf("\n");

		/*Unfortunately, now the new node has    */
		/*been switched into position lin1 in the*/
		/*vector. Therefore, restrict the maximum*/
		/*again and if lin1 is chosen, override  */
		/*it and choose the end lineage.         */

		/*Choose 2nd lineage to coalesce*/			//New node is now in position lin1, so if it
		int lin2=ran->discrete(0,ARG_k-2);			//is chosen, force it to choose the last old node
		if(lin2==lin1)lin2=ARG_k-1;					//instead. As a result, do not let the last old node
		class ap_node *ap_node2=active_node[lin2];	//(in position ARG_k-1) be chosen initially.
		deactivate_node(lin2);
		//printf("           : ");for( o=0;o<ARG_k;o++)printf("%d ",active_node[o]->id);printf("\n");

		if((new_node->id==ap_node1->id)||(new_node->id==ap_node2->id)
			||(ap_node1->id==ap_node2->id))error("coalesce(): nodes not chosen correctly");
	//	printf("Identity of New node: %d Coalescing nodes: %d and %d\n",new_node->id,ap_node1->id,ap_node2->id);
	//	printf("                           which point to: %d and %d\n",ap_node1->AMP.ptr[0]->id,ap_node2->AMP.ptr[0]->id);

		/*Perform copying and coalescing*/
		int imax=n_segregating;
		int i;
		for(i=0;i<imax;i++)
		{
			int tree_id=(*segregating_tree)[i];
			if(ap_node1->AMP.ptr[tree_id]==NULL)
			{
				/*Rule 2.ii  */
				if(ap_node2->AMP.ptr[tree_id]==NULL)
				{	++ncoI;
				new_node->AMP.assign(NULL,tree_id);}
				/*Rule 2.i   */
				else
				{	++ncoIIa;
				new_node->AMP.assign(ap_node2->AMP.ptr[tree_id],tree_id);}
			}
			else
			{
				/*Rule 2.ii  */
				if(ap_node2->AMP.ptr[tree_id]==NULL)
				{	++ncoIIb;
				new_node->AMP.assign(ap_node1->AMP.ptr[tree_id],tree_id);}
				/*Rule 2.iii */
				else
				{
					++ncoIII;
					new_node->AMP.assign(tree[tree_id].coalesce(*time,ap_node1->AMP.ptr[tree_id]->id,ap_node2->AMP.ptr[tree_id]->id),tree_id);
				}
			}
		}
		if(!samples_waiting)
		{
			deactivate_trees();
			/*i=0;
			while(i<n_segregating)
			{
				int tree_id=(*segregating_tree)[i];
				if(tree[tree_id].get_k()==1)deactivate_tree(tree_id);
				else i++;
			}*/
		}
		/*Recalculate rlen*/
		calc_rlen(new_node);

	//	printf("                            coalescing to: %d\n",new_node->AMP.ptr[0]->id);		
		return *this;
	}
	virtual coalescent& recombine(double *time)
	{
		++nrec;
		int rtype=-1;

		/*Create new lineages in the ARG*/
		class ap_node *new_node1=create_node(time);
		class ap_node *new_node2=create_node(time);

		/*First choose the lineage*/
		double rnum1=ran->U()*total_rlen;
		int lin;
		for(lin=0;lin<ARG_k;lin++)
		{
			if (rnum1<=active_node[lin]->rlen) break;
			rnum1 -= active_node[lin]->rlen;
		}
		if (lin>=ARG_k) error("recombine(): lineage not chosen correctly");
		class ap_node *old_node=active_node[lin];
		deactivate_node(lin);

		//active_node[lin]->edge_time=(*time)-active_node[lin]->time;
		if(new_node1==old_node)error("Aah!");
		if(new_node2==old_node)error("Aah!");

		/*Determine the number of breakpoints*/
		if (old_node->L==0.0)error("recombine(): recombination at an empty locus");

		int ltr=old_node->ltr;
		int rtr=old_node->rtr;

		/*Is it a swap?*/
		if(no_gene_conversion)
		{
			perform_single_crossover(&ltr,&rtr);
			if(tree[ltr-1].get_k()>1 && tree[ltr].get_k()>1 && ltr!=old_node->ltr && old_node->AMP.ptr[ltr-1]!=NULL && old_node->AMP.ptr[ltr]!=NULL)
				++nrecWatt[ltr-1];
			else if(tree[rtr-1].get_k()>1 && tree[rtr].get_k()>1 && rtr!=old_node->rtr && old_node->AMP.ptr[rtr-1]!=NULL && old_node->AMP.ptr[rtr]!=NULL)
				++nrecWatt[rtr-1];
			rtype=2;
			++nrecTypeII;
		}
/*		else
		{
			double a,b,c,swap_yn,rnum3,single_yn;
			a = con->lambda*old_node->L;
			b = exp(-a);
			c = a+b;
			swap_yn = b/c;	//=1 if L=0 so "swap", but has no effect
							//this error shouldve been caught anyway
			rnum3 = ran->U();

			if (rnum3<=swap_yn)
			{
				//It's a swap!
				//In which case all of the recipient's genome is
				//ancestral except for the locus of interest
				//ltr and rtr do not need modifying
				rtype=1;
				++nrecTypeI;
			}*/
			else
			{
				/*Is it a single cross-over?*/
				//single_yn = (1-a)/c;
				//rnum3 -= swap_yn;
				double rnum3 = ran->U() * old_node->rlen;
				/*	Before 11.08.06 the next line was the same, but changes to calc_node_rlen
					imply that the relative rate of single to double xovers is altered.		*/
				double single_yn = con->r/con->lambda*(1.-pow(1.-con->lambda,(double)(old_node->L-1)));
				if (rnum3<=single_yn)
				{
					/*It's a single cross-over!*/
					perform_single_crossover(&ltr,&rtr);
					rtype=2;
					++nrecTypeII;
					if(tree[ltr-1].get_k()>1 && tree[ltr].get_k()>1 && ltr!=old_node->ltr && old_node->AMP.ptr[ltr-1]!=NULL && old_node->AMP.ptr[ltr]!=NULL)
						++nrecWatt[ltr-1];
					else if(tree[rtr-1].get_k()>1 && tree[rtr].get_k()>1 && rtr!=old_node->rtr && old_node->AMP.ptr[rtr-1]!=NULL && old_node->AMP.ptr[rtr]!=NULL)
						++nrecWatt[rtr-1];
				}
				else
				{
					/*It's a double cross-over!*/
					perform_double_crossover(&ltr,&rtr);
					rtype=3;
					++nrecTypeIII;
					//++nrecWatt[ltr-1];
					//++nrecWatt[rtr-1];
				}
			}
		//}

		/*Copy the relevant parts of AMP*/
		int i,pos;
		for(i=0,pos=(*segregating_tree)[0];(pos<ltr)&&(i<n_segregating);pos=(*segregating_tree)[i+1],i++)
		{
			new_node1->AMP.assign(old_node->AMP.ptr[pos],pos);
			new_node2->AMP.assign(NULL,pos);
		}
		for(;(pos<rtr)&&(i<n_segregating);pos=(*segregating_tree)[i+1],i++)
		{
			new_node2->AMP.assign(old_node->AMP.ptr[pos],pos);
			new_node1->AMP.assign(NULL,pos);
		}
		for(;i<n_segregating;pos=(*segregating_tree)[i+1],i++)
		{
			new_node1->AMP.assign(old_node->AMP.ptr[pos],pos);
			new_node2->AMP.assign(NULL,pos);
		}

		/*Recalculate rlen*/
		calc_rlen(new_node1,new_node2);

		return *this;
	}
	coalescent& conditionally_coalesce(double *time)
	{
		++nco;
		/*Create new lineage in the ARG*/
		class ap_node *new_node=create_node(time);
		int lin1,lin2;
		class ap_node *ap_node1,*ap_node2;
		if(true){//ARG_k_fixed>=ARG_k) {
			while(true) {
				lin1=ran->discrete(0,ARG_k-2);
				ap_node1=active_node[lin1];
				/* Always accept if not a FIXED_NODE */
				if(ap_node1->flag!=ap_node::FIXED_NODE) break;
				/* otherwise accept with probability */
				else if(ran->U()<1.-(ARG_k_fixed-1.)/(ARG_k-2.)) break;
			}
			/* If ap_node1 is a FIXED_NODE make new_node a FIXED_NODE */
			if(ap_node1->flag==ap_node::FIXED_NODE) {
				new_node->flag = ap_node::FIXED_NODE;
				new_node->ctree_id = ap_node1->ctree_id;
			}
			deactivate_node(lin1);
			while(true) {
				/* Choose a different lineage */
				lin2=ran->discrete(0,ARG_k-2);
				if(lin2==lin1)lin2=ARG_k-1;
				ap_node2=active_node[lin2];
				/* Don't accept if both nodes have FIXED_NODE status */
				if(!(new_node->flag==ap_node::FIXED_NODE && ap_node2->flag==ap_node::FIXED_NODE)) break;
			}
			/* If ap_node2 is a FIXED_NODE make new_node a FIXED_NODE */
			if(ap_node2->flag==ap_node::FIXED_NODE) {
				new_node->flag = ap_node::FIXED_NODE;
				new_node->ctree_id = ap_node2->ctree_id;
			}
			//else new_node->ctree_id = -1;
			if(new_node->ctree_id>-1) fnode[new_node->ctree_id] = new_node;
			/*Finally, deactivate*/
			deactivate_node(lin2);
		}
		else {
		}
		/*Choose 1st lineage to coalesce*
		int lin1=ran->discrete(0,ARG_k-2);
		/*Do not choose FIXED_NODEs *
		while(active_node[lin1]->flag==ap_node::FIXED_NODE) lin1 = ran->discrete(0,ARG_k-2);
		class ap_node *ap_node1=active_node[lin1];
		deactivate_node(lin1);
		/*Choose 2nd lineage to coalesce*
		int lin2=ran->discrete(0,ARG_k-2);
		if(lin2==lin1)lin2=ARG_k-1;
		class ap_node *ap_node2=active_node[lin2];
		/*Sort out fnode stuff*
		new_node->flag = ap_node2->flag;
		new_node->ctree_id = ap_node2->ctree_id;
		if(new_node->ctree_id>-1) fnode[new_node->ctree_id] = new_node;
		/*Finally, deactivate*
		deactivate_node(lin2);*/

		if((new_node->id==ap_node1->id)||(new_node->id==ap_node2->id)
			||(ap_node1->id==ap_node2->id))error("coalesce(): nodes not chosen correctly");

		/*Give new_node the flag of ap_node2, which might be a FIXED_NODE
		int found_fnode = 0;
		int a;
		for(a=0;a<tree[0].size;a++) 
			if(fnode[a]==ap_node2) {
				fnode[a] = new_node;
				++found_fnode;
			}
		if(new_node->flag==ap_node::FIXED_NODE && found_fnode!=1) error("coalescent::conditionally_coalesce(): problem finding FIXED_NODE");*/

		/*Perform copying and coalescing*/
		int imax=n_segregating;
		int i;
		for(i=0;i<imax;i++)
		{
			int tree_id=(*segregating_tree)[i];
			if(ap_node1->AMP.ptr[tree_id]==NULL)
			{
				/*Rule 2.ii  */
				if(ap_node2->AMP.ptr[tree_id]==NULL)
				{	++ncoI;
				new_node->AMP.assign(NULL,tree_id);}
				/*Rule 2.i   */
				else
				{	++ncoIIa;
				new_node->AMP.assign(ap_node2->AMP.ptr[tree_id],tree_id);}
			}
			else
			{
				/*Rule 2.ii  */
				if(ap_node2->AMP.ptr[tree_id]==NULL)
				{	++ncoIIb;
				new_node->AMP.assign(ap_node1->AMP.ptr[tree_id],tree_id);}
				/*Rule 2.iii */
				else
				{
					++ncoIII;
					new_node->AMP.assign(tree[tree_id].coalesce(*time,ap_node1->AMP.ptr[tree_id]->id,ap_node2->AMP.ptr[tree_id]->id),tree_id);
				}
			}
		}
		if(!samples_waiting)
		{
			//deactivate_trees2();
			i=0;
			while(i<n_segregating)
			{
				int tree_id=(*segregating_tree)[i];
				if(tree[tree_id].get_k()==1)deactivate_tree(tree_id);
				else i++;
			}
		}
		/*Recalculate rlen*/
		calc_rlen(new_node);
		
		return *this;
	}
	coalescent& conditionally_recombine(double *time)
	{
		++nrec;
		int rtype=-1;

		/*Create new lineages in the ARG*/
		class ap_node *new_node1=create_node(time);
		class ap_node *new_node2=create_node(time);

		/*First choose the lineage*/
		double rnum1=ran->U()*total_rlen;
		int lin;
		for(lin=0;lin<ARG_k;lin++)
		{
			if (rnum1<=active_node[lin]->rlen) break;
			rnum1 -= active_node[lin]->rlen;
		}
		if (lin>=ARG_k) error("recombine(): lineage not chosen correctly");
		class ap_node *old_node=active_node[lin];
		/*new_node1 is always the recipient*/
		new_node1->flag = old_node->flag;
		new_node1->ctree_id = old_node->ctree_id;
		if(new_node1->ctree_id!=-1) fnode[new_node1->ctree_id] = new_node1;
		int old_node_flag = (int)old_node->flag;
		deactivate_node(lin);

		//active_node[lin]->edge_time=(*time)-active_node[lin]->time;
		if(new_node1==old_node)error("Aah!");
		if(new_node2==old_node)error("Aah!");

		/*Determine the number of breakpoints*/
		if (old_node->L==0.0)error("recombine(): recombination at an empty locus");

		int ltr=old_node->ltr;
		int rtr=old_node->rtr;

		/*Is it a swap?*/
		if(no_gene_conversion)
		{
			//error("coalescent::conditionally_recombine(): only donor-recipient style rec defined");
			perform_single_crossover(&ltr,&rtr);
			if(tree[ltr-1].get_k()>1 && tree[ltr].get_k()>1 && ltr!=old_node->ltr && old_node->AMP.ptr[ltr-1]!=NULL && old_node->AMP.ptr[ltr]!=NULL)
				++nrecWatt[ltr-1];
			else if(tree[rtr-1].get_k()>1 && tree[rtr].get_k()>1 && rtr!=old_node->rtr && old_node->AMP.ptr[rtr-1]!=NULL && old_node->AMP.ptr[rtr]!=NULL)
				++nrecWatt[rtr-1];
			rtype=2;
			++nrecTypeII;
		}
		else
		{
			double swap_yn,rnum3,single_yn;
			rnum3 = ran->U() * old_node->rlen;
			swap_yn = (old_node_flag == (int)ap_node::FIXED_NODE) ? 0.5 * con->r/con->lambda*pow(1.-con->lambda,(double)(old_node->L-1)) : 0.0;
			/* Before 11.08.06 swap_yn = (old_node_flag == (int)ap_node::FIXED_NODE) ? con->r/con->lambda*pow(1.-con->lambda,(double)(old_node->L-1)) : 0.0;*/

			if (rnum3<=swap_yn)
			{
				//It's a swap!
				//In which case all of the recipient's genome is
				//ancestral except for the locus of interest
				//ltr and rtr do not need modifying
				rtype=1;
				++nrecTypeI;
			}
			else
			{
				rnum3 -= swap_yn;
				/*Is it a single cross-over?*/
				//single_yn = (1-a)/c;
				//rnum3 -= swap_yn;
				/*	The following line was the same before 11.08.06, but the changes in calc_node_rlen
					imply that the relative probability of double crossovers are altered as a result.   */
				single_yn = con->r/con->lambda*(1.-pow(1.-con->lambda,(double)(old_node->L-1)));
				if (rnum3<=single_yn)
				{
					/*It's a single cross-over!*/
					perform_single_crossover(&ltr,&rtr);
					rtype=2;
					++nrecTypeII;
					if(tree[ltr-1].get_k()>1 && tree[ltr].get_k()>1 && ltr!=old_node->ltr && old_node->AMP.ptr[ltr-1]!=NULL && old_node->AMP.ptr[ltr]!=NULL)
						++nrecWatt[ltr-1];
					else if(tree[rtr-1].get_k()>1 && tree[rtr].get_k()>1 && rtr!=old_node->rtr && old_node->AMP.ptr[rtr-1]!=NULL && old_node->AMP.ptr[rtr]!=NULL)
						++nrecWatt[rtr-1];
				}
				else
				{
					/*It's a double cross-over!*/
					perform_double_crossover(&ltr,&rtr);
					rtype=3;
					++nrecTypeIII;
					//++nrecWatt[ltr-1];
					//++nrecWatt[rtr-1];
				}
			}
		}
		/*new_node1 is always the recipient
		int found_fnode = 0;
		int a;
		for(a=0;a<tree[0].size;a++)
			if(fnode[a]==old_node) {
				fnode[a] = new_node1;
				++found_fnode;
			}
		if(new_node1->flag==ap_node::FIXED_NODE && found_fnode!=1) error("coalescent::conditionally_recombine(): problem finding FIXED_NODE");*/

		/*Copy the relevant parts of AMP*/
		int i,pos;
		for(i=0,pos=(*segregating_tree)[0];(pos<ltr)&&(i<n_segregating);pos=(*segregating_tree)[i+1],i++)
		{
			new_node1->AMP.assign(old_node->AMP.ptr[pos],pos);
			new_node2->AMP.assign(NULL,pos);
		}
		for(;(pos<rtr)&&(i<n_segregating);pos=(*segregating_tree)[i+1],i++)
		{
			new_node2->AMP.assign(old_node->AMP.ptr[pos],pos);
			new_node1->AMP.assign(NULL,pos);
		}
		for(;i<n_segregating;pos=(*segregating_tree)[i+1],i++)
		{
			new_node1->AMP.assign(old_node->AMP.ptr[pos],pos);
			new_node2->AMP.assign(NULL,pos);
		}

		/*Recalculate rlen*/
		calc_rlen(new_node1,new_node2);

		return *this;
	}
	coalescent& migrate_coalesce(double *time)
	{
		++nco;
		/*First choose deme for coalescence*/
		double rdeme = ran->U() * rate_coal;
		int deme;
		for(deme=0;deme<con->ndemes;deme++) {
			if(rdeme <= coal_deme[deme]/2.) break;
			else rdeme -= coal_deme[deme]/2.;
		}
		if(deme==con->ndemes) error("migrate_coalesce(): deme chosen incorrectly");

		/*Create new lineage in the ARG*/
		class ap_node *new_node=create_node(time);
		int lin1,lin2;
		class ap_node *ap_node1,*ap_node2;
		lin1 = ran->discrete(0,k_deme[deme]-1);
		ap_node1 = ptr_deme[deme][lin1];
		if(ap_node1->deme!=deme) error("migrate_coalesce(): node 1 deme not right deme");
		SWAP(ptr_deme[deme][lin1],ptr_deme[deme][k_deme[deme]-1]);
		--k_deme[deme];
		deactivate_node(ap_node1->active_id);
		lin2 = ran->discrete(0,k_deme[deme]-1);
		ap_node2 = ptr_deme[deme][lin2];
		if(ap_node2->deme!=deme) error("migrate_coalesce(): node 2 deme not right deme");
		ptr_deme[deme][lin2] = new_node;
		deactivate_node(ap_node2->active_id);
		new_node->deme = deme;

		if((new_node->id==ap_node1->id)||(new_node->id==ap_node2->id)
			||(ap_node1->id==ap_node2->id))error("migrate_coalesce(): nodes not chosen correctly");

		/*Perform copying and coalescing*/
		int imax=n_segregating;
		int i;
		mt_node *ptr; double last_update;
		for(i=0;i<imax;i++)
		{
			int tree_id=(*segregating_tree)[i];
			if(ap_node1->AMP.ptr[tree_id]==NULL)
			{
				/*Rule 2.ii  */
				if(ap_node2->AMP.ptr[tree_id]==NULL)
				{	++ncoI;
				new_node->AMP.assign(NULL,tree_id);}
				/*Rule 2.i   */
				else
				{	++ncoIIa;
				new_node->AMP.assign(ap_node2->AMP.ptr[tree_id],tree_id);}
			}
			else
			{
				/*Rule 2.ii  */
				if(ap_node2->AMP.ptr[tree_id]==NULL)
				{	++ncoIIb;
				new_node->AMP.assign(ap_node1->AMP.ptr[tree_id],tree_id);}
				/*Rule 2.iii */
				else
				{
					++ncoIII;
					/*** Structured coalescent stuff ***/
					ptr = ap_node1->AMP.ptr[tree_id];
					last_update = ptr->last_update;
					ptr->edge_time += (*time-last_update) * con->N_deme_over_D[deme];
					ptr->last_update = *time;
					ptr = ap_node2->AMP.ptr[tree_id];
					last_update = ptr->last_update;
					ptr->edge_time += (*time-last_update) * con->N_deme_over_D[deme];
					ptr->last_update = *time;
					/***********************************/
					new_node->AMP.assign(tree[tree_id].migrate_coalesce(*time,ap_node1->AMP.ptr[tree_id]->id,ap_node2->AMP.ptr[tree_id]->id),tree_id);
				}
			}
		}
		if(!samples_waiting)
		{
			//deactivate_trees2();
			i=0;
			while(i<n_segregating)
			{
				int tree_id=(*segregating_tree)[i];
				if(tree[tree_id].get_k()==1)deactivate_tree(tree_id);
				else i++;
			}
		}
		/*Recalculate rlen*/
		migrate_calc_rlen(new_node);

		/*Recalculate rates*
		rate_recn -= rho_deme[deme];
		rho_deme[deme] -= ap_node1->rlen / con->N_deme_over_D[deme];
		rho_deme[deme] -= ap_node2->rlen / con->N_deme_over_D[deme];
		calc_node_rlen(new_node);
		rho_deme[deme] += new_node->rlen / con->N_deme_over_D[deme];
		rate_recn += rho_deme[deme];

		rate_coal -= k_deme[deme] / con->N_deme_over_D[deme];
		coal_deme[deme] -= k_deme[deme] / con->N_deme_over_D[deme];

		rate_mign -= sum_mig[deme];*/
		
		return *this;
	}
	coalescent& migrate_recombine(double *time)
	{
		++nrec;
		int rtype=-1;

		/*First choose deme*/
		double rdeme = ran->U() * rate_recn;
		int deme;
		for(deme=0;deme<con->ndemes;deme++) {
			if(rdeme <= rho_deme[deme]/2.) break;
			else rdeme -= rho_deme[deme]/2.;
		}
		if(deme==con->ndemes) error("migrate_recombine(): deme not chosen correctly");

		/*Create new lineages in the ARG*/
		class ap_node *new_node1=create_node(time);
		class ap_node *new_node2=create_node(time);

		/*First choose the lineage*/
		double rnum1=ran->U() * rho_deme[deme] / 2.0 * con->N_deme_over_D[deme];
		class ap_node *old_node;
		int lin;
		for(lin=0;lin<k_deme[deme];lin++)
		{
			old_node = ptr_deme[deme][lin];
			if(rnum1<=old_node->rlen) break;
			else rnum1 -= old_node->rlen;
		}
		if(lin>=k_deme[deme]) error("migrate_recombine(): lineage not chosen correctly");
		if(old_node->deme!=deme) error("migrate_recombine(): lineage deme not right deme");
		new_node1->deme = new_node2->deme = deme;
		ptr_deme[deme][lin] = new_node1;
		++k_deme[deme];
		ptr_deme[deme][k_deme[deme]-1] = new_node2;
		/*new_node1 is always the recipient*/
		deactivate_node(old_node->active_id);
		if(new_node1==old_node)error("Aah!");
		if(new_node2==old_node)error("Aah!");

		/*Determine the number of breakpoints*/
		if(old_node->L==0.0)error("recombine(): recombination at an empty locus");

		int ltr=old_node->ltr;
		int rtr=old_node->rtr;

		/*Is it a swap?*/
		if(no_gene_conversion)
		{
			//error("coalescent::conditionally_recombine(): only donor-recipient style rec defined");
			perform_single_crossover(&ltr,&rtr);
			if(tree[ltr-1].get_k()>1 && tree[ltr].get_k()>1 && ltr!=old_node->ltr && old_node->AMP.ptr[ltr-1]!=NULL && old_node->AMP.ptr[ltr]!=NULL)
				++nrecWatt[ltr-1];
			else if(tree[rtr-1].get_k()>1 && tree[rtr].get_k()>1 && rtr!=old_node->rtr && old_node->AMP.ptr[rtr-1]!=NULL && old_node->AMP.ptr[rtr]!=NULL)
				++nrecWatt[rtr-1];
			rtype=2;
			++nrecTypeII;
		}
		else
		{
			double swap_yn,rnum3,single_yn;
			rnum3 = ran->U() * old_node->rlen;
			swap_yn = con->r/con->lambda*pow(1.-con->lambda,(double)(old_node->L-1));

			if (rnum3<=swap_yn)
			{
				//It's a swap!
				//In which case all of the recipient's genome is
				//ancestral except for the locus of interest
				//ltr and rtr do not need modifying
				rtype=1;
				++nrecTypeI;
			}
			else
			{
				rnum3 -= swap_yn;
				/*Is it a single cross-over?*/
				//single_yn = (1-a)/c;
				//rnum3 -= swap_yn;
				single_yn = con->r/con->lambda*(1.-pow(1.-con->lambda,(double)(old_node->L-1)));
				if (rnum3<=single_yn)
				{
					/*It's a single cross-over!*/
					perform_single_crossover(&ltr,&rtr);
					rtype=2;
					++nrecTypeII;
					if(tree[ltr-1].get_k()>1 && tree[ltr].get_k()>1 && ltr!=old_node->ltr && old_node->AMP.ptr[ltr-1]!=NULL && old_node->AMP.ptr[ltr]!=NULL)
						++nrecWatt[ltr-1];
					else if(tree[rtr-1].get_k()>1 && tree[rtr].get_k()>1 && rtr!=old_node->rtr && old_node->AMP.ptr[rtr-1]!=NULL && old_node->AMP.ptr[rtr]!=NULL)
						++nrecWatt[rtr-1];
				}
				else
				{
					/*It's a double cross-over!*/
					perform_double_crossover(&ltr,&rtr);
					rtype=3;
					++nrecTypeIII;
					//++nrecWatt[ltr-1];
					//++nrecWatt[rtr-1];
				}
			}
		}

		/*Copy the relevant parts of AMP*/
		int i,pos;
		for(i=0,pos=(*segregating_tree)[0];(pos<ltr)&&(i<n_segregating);pos=(*segregating_tree)[i+1],i++)
		{
			new_node1->AMP.assign(old_node->AMP.ptr[pos],pos);
			new_node2->AMP.assign(NULL,pos);
		}
		for(;(pos<rtr)&&(i<n_segregating);pos=(*segregating_tree)[i+1],i++)
		{
			new_node2->AMP.assign(old_node->AMP.ptr[pos],pos);
			new_node1->AMP.assign(NULL,pos);
		}
		for(;i<n_segregating;pos=(*segregating_tree)[i+1],i++)
		{
			new_node1->AMP.assign(old_node->AMP.ptr[pos],pos);
			new_node2->AMP.assign(NULL,pos);
		}

		/*Recalculate rlen*/
		migrate_calc_rlen(new_node1,new_node2);

		/*Recalculate rates*
		rate_recn -= rho_deme[deme];
		rho_deme[deme] -= old_node->rlen / con->N_deme_over_D[deme];
		calc_node_rlen(new_node1);
		calc_node_rlen(new_node2);
		rho_deme[deme] += new_node1->rlen / con->N_deme_over_D[deme];
		rho_deme[deme] += new_node2->rlen / con->N_deme_over_D[deme];
		rate_recn += rho_deme[deme];

		rate_coal += (double)(k_deme[deme]-1) / con->N_deme_over_D[deme];
		coal_deme[deme] += (double)(k_deme[deme]-1) / con->N_deme_over_D[deme];

		rate_mign += sum_mig[deme];*/

		return *this;
	}
	coalescent& migrate_migrate(double *time) {
		/*Choose source deme*/
		int source;
		double rdeme = ran->U() * rate_mign;
		for(source=0;source<con->ndemes;source++) {
			if(rdeme <= k_deme[source]*sum_mig[source]/2.) break;
			else rdeme -= k_deme[source]*sum_mig[source]/2.;
		}
		if(source>=con->ndemes) error("migrate_migrate(): source deme not chosen correctly");

		/*Choose target deme*/
		int target;
		rdeme = ran->U() * sum_mig[source];
		for(target=0;target<con->ndemes;target++) {
			if(rdeme <= con->mig[source][target]) break;
			else rdeme -= con->mig[source][target];
		}
		if(target>=con->ndemes) error("migrate_migrate(): target deme not chosen correctly");

		/*Choose lineage*/
		int lin = ran->discrete(0,k_deme[source]-1);
		ap_node *old_node = ptr_deme[source][lin];
		if(old_node->deme!=source) error("migrate_migrate(): lineage belongs to wrong deme");

		/*Perform migration*/
		if(old_node->deme!=target) {
			SWAP(ptr_deme[source][lin],ptr_deme[source][k_deme[source]-1]);
			--k_deme[source];
			++k_deme[target];
			ptr_deme[target][k_deme[target]-1] = old_node;
			old_node->deme = target;
		}

		/*Update edge_time for migrating node*/
		int i,pos;
		mt_node *ptr;
		double last_update;
		for(i=0;i<n_segregating;i++) {
			pos = (*segregating_tree)[i];
			ptr = old_node->AMP.ptr[pos];
			if(ptr!=NULL) {
				last_update = ptr->last_update;
				ptr->edge_time += (*time-last_update) * con->N_deme_over_D[source];
				ptr->last_update = *time;
			}
		}

		/*Recalculate rlen*/
		migrate_calc_rlen();

		/*Recalculate rates*
		rate_recn -= rho_deme[source] + rho_deme[target];
		rho_deme[source] -= old_node->rlen / con->N_deme_over_D[source];
		rho_deme[target] += old_node->rlen / con->N_deme_over_D[target];
		rate_recn += rho_deme[source] + rho_deme[target];
		
		rate_coal -= coal_deme[source] + coal_deme[target];
		coal_deme[source] -= (double)k_deme[source] / con->N_deme_over_D[source];
		coal_deme[target] += (double)(k_deme[target]-1) / con->N_deme_over_D[target];
		rate_coal += coal_deme[source] + coal_deme[target];
		
		rate_mign += sum_mig[target] - sum_mig[source];*/

		return *this;
	}
	double add_next_sample()
	{
		double samptime=*ntimes_itr;
		double currenttime=samptime;
		class ap_node* new_node;
		class mt_node* new_tree_node;
		
		while((currenttime==*ntimes_itr)&&(ntimes_itr!=con->ntimes.end())) //relies on ordering of ntimes
		{
			//add the new node
			++naddbase;
			new_node=create_node(&(*ntimes_itr));
			int i;
			for(i=0;i<L;i++)
			{
				new_tree_node=tree[i].add_base_node(&(*ntimes_itr),next_waiting_sample);
				new_node->AMP.assign(new_tree_node,i);
			}
			calc_node_rlen(new_node);
			if(con->ndemes>0) {
				int deme = con->deme_config[next_waiting_sample];
				new_node->deme = deme;
				++k_deme[deme];
				ptr_deme[deme][k_deme[deme]-1] = new_node;

				/*rho_deme[deme] += new_node->rlen / con->N_deme_over_D[deme];
				rate_recn += new_node->rlen / con->N_deme_over_D[deme];

				rate_coal += (double)(k_deme[deme]-1) / con->N_deme_over_D[deme];
				coal_deme[deme] += (double)(k_deme[deme]-1) / con->N_deme_over_D[deme];

				rate_mign += sum_mig[deme];*/
			}
			++next_waiting_sample;
			//++ARG_k;

			++ntimes_itr;
		}
		
		if(ntimes_itr==con->ntimes.end()) samples_waiting=false;
		else time_next_sample=*ntimes_itr;

		if(con->ndemes>0) {
			migrate_calc_rlen();
		}
		else calc_rlen();
		return currenttime;
	}
	double add_conditional_event(class marginal_tree &ctree)
	{
		double samptime=*ntimes_itr;
		double currenttime=samptime;
		class ap_node* new_node;
		class mt_node* new_tree_node;

		int i;
		class mt_node* conditional_event;
		
		while((currenttime==*ntimes_itr)&&(ntimes_itr!=con->ntimes.end())) //relies on ordering of ntimes
		{
			conditional_event = &(ctree.node[next_waiting_sample]);
			if(conditional_event->descendant[0]==NULL) {	// add a base node
				++naddbase;
				new_node = create_node(&(*ntimes_itr));
				for(i=0;i<L;i++)
				{
					new_tree_node=tree[i].add_base_node(&(*ntimes_itr),next_waiting_sample);
					new_node->AMP.assign(new_tree_node,i);
				}
				fnode[next_waiting_sample] = new_node;
				new_node->flag = ap_node::FIXED_NODE;
				new_node->ctree_id = next_waiting_sample;
				calc_node_rlen(new_node);
				++ARG_k_fixed;
			}
			else {	// add a coalescence
				//fnode[conditional_event->id] = coalesce(fnode[conditional_event->descendant[0]->id],fnode[conditional_event->descendant[1]->id]);
				class ap_node *new_node = create_node(&(*ntimes_itr));
				class ap_node *ap_node1 = fnode[conditional_event->descendant[0]->id];
				deactivate_node(ap_node1->active_id);
				class ap_node *ap_node2 = fnode[conditional_event->descendant[1]->id];
				if(ap_node1==ap_node2) error("add_conditional_event(): lineage cannot coalesce with itself");
				deactivate_node(ap_node2->active_id);
				if((new_node->id==ap_node1->id)||(new_node->id==ap_node2->id)
					||(ap_node1->id==ap_node2->id))error("add_conditional_event(): nodes not chosen correctly");
				fnode[conditional_event->id] = new_node;
				fnode[conditional_event->descendant[0]->id] = NULL;
				fnode[conditional_event->descendant[1]->id] = NULL;
				/* give node FIXED_NODE status only if it is not the mrca */
				if(conditional_event->id<ctree.size-1) {
					new_node->flag = ap_node::FIXED_NODE;
					--ARG_k_fixed;
				}
				else ARG_k_fixed -= 2;
				new_node->ctree_id = conditional_event->id;
				/*Perform copying and coalescing*/
				int imax=n_segregating;
				for(i=0;i<imax;i++)
				{
					int tree_id=(*segregating_tree)[i];
					if(ap_node1->AMP.ptr[tree_id]==NULL)
					{
						/*Rule 2.ii  */
						if(ap_node2->AMP.ptr[tree_id]==NULL)
						{	++ncoI;
						new_node->AMP.assign(NULL,tree_id);}
						/*Rule 2.i   */
						else
						{	++ncoIIa;
						new_node->AMP.assign(ap_node2->AMP.ptr[tree_id],tree_id);}
					}
					else
					{
						/*Rule 2.ii  */
						if(ap_node2->AMP.ptr[tree_id]==NULL)
						{	++ncoIIb;
						new_node->AMP.assign(ap_node1->AMP.ptr[tree_id],tree_id);}
						/*Rule 2.iii */
						else
						{
							++ncoIII;
							new_node->AMP.assign(tree[tree_id].coalesce(*ntimes_itr,ap_node1->AMP.ptr[tree_id]->id,ap_node2->AMP.ptr[tree_id]->id),tree_id);
						}
					}
				}
				if(!samples_waiting)
				{
					//deactivate_trees2();
					i=0;
					while(i<n_segregating)
					{
						int tree_id=(*segregating_tree)[i];
						if(tree[tree_id].get_k()==1)deactivate_tree(tree_id);
						else i++;
					}
				}
				calc_node_rlen(new_node);
			}
			++next_waiting_sample;
			++ntimes_itr;
		}
		
		/* Can deactivate some trees when last conditional event occurs */
		if(ntimes_itr==ftimes.end()) {
			samples_waiting=false;
			deactivate_trees();
			for(i=0;i<ARG_k;i++) calc_node_rlen(active_node[i]);
		}
		else time_next_sample=*ntimes_itr;
		if(ARG_k_fixed>1 && samples_waiting==false) error("add_conditional_event(): not all fixed events completed");

		calc_rlen();
		return currenttime;
	}
	double constant_size_model(const double mean)
	{
		double time=ran->exponential(1.0);
		time *= mean;
		return time;
	}
	coalescent& deactivate_node(const int id)
	{
		inactive_node[n_inactive]=active_node[id];
		inactive_node[n_inactive]->recycle();
		++n_inactive;
		active_node[id]=active_node[ARG_k-1];
		active_node[id]->active_id = id;
		active_node[ARG_k-1]=NULL;
		--ARG_k;
		/*NB no memory reallocation occurs*/
		return *this;
	}
	/*coalescent& deactivate_node(ap_node *id)
	{
		int lin1 = 0;
		while(active_node[lin1]!=id) ++lin1;
		return deactivate_node(lin1);
	}*/
	class ap_node* create_node(double *time)
	{
		if(ARG_k==nodes_reserved)reserve_nodes(2*(ARG_k+1));
		active_node[ARG_k]=inactive_node[n_inactive-1];
		active_node[ARG_k]->activate(time);
		active_node[ARG_k]->active_id = ARG_k;
		inactive_node[n_inactive-1]=NULL;
		++ARG_k;
		--n_inactive;
		return active_node[ARG_k-1];
	}
	coalescent& deactivate_trees()
	{
		int new_seg_tree_id=1-seg_tree_id;
		int i,j;
		for(i=0,j=0;i<n_segregating;i++) /*i counts along old tree, j along new tree*/
		{
			int tree_id=(*segregating_tree)[i];
			if(tree[tree_id].get_k()!=1)
			{
				internal_seg_tree[new_seg_tree_id][j]=tree_id;
				j++;
			}
		}
		int new_n_segregating = j;
		for(;j<n_segregating;j++) internal_seg_tree[new_seg_tree_id][j]=-1;
		/*if(new_n_segregating<n_segregating) {
			double before,after;
			for(i=0;i<ARG_k;i++) {
				before = active_node[i]->L;
				calc_node_rlen(active_node[i]);
				after = active_node[i]->L;
				if(before!=after) {
					warning("weird");
				}
			}
			calc_rlen();
		}*/
		n_segregating = new_n_segregating;
		seg_tree_id=new_seg_tree_id;
		segregating_tree=&(internal_seg_tree[seg_tree_id]);
		return *this;
	}
	coalescent& deactivate_trees2()
	{
		int new_seg_tree_id=seg_tree_id;
		int i,j;
		for(i=0,j=0;i<n_segregating;i++) /*i counts along old tree, j along new tree*/
		{
			int tree_id=(*segregating_tree)[i];
			if(tree[tree_id].get_k()!=1)
			{
				internal_seg_tree[new_seg_tree_id][j]=tree_id;
				j++;
			}
		}
		int new_n_segregating = j;
		for(;j<n_segregating;j++) internal_seg_tree[new_seg_tree_id][j]=-1;
		/*if(new_n_segregating<n_segregating) {
			double before,after;
			for(i=0;i<ARG_k;i++) {
				before = active_node[i]->L;
				calc_node_rlen(active_node[i]);
				after = active_node[i]->L;
				if(before!=after) {
					warning("weird");
				}
			}
			calc_rlen();
		}*/
		n_segregating = new_n_segregating;
		seg_tree_id=new_seg_tree_id;
		segregating_tree=&(internal_seg_tree[seg_tree_id]);
		return *this;
	}
	coalescent& deactivate_tree(const int id)
	{
		int new_seg_tree_id=1-seg_tree_id;
		//segregating_tree[id]=segregating_tree[n_segregating-1];
		int j=0;
		int i;
		for(i=0;i<n_segregating;i++) /*i counts along old tree, j along new tree*/
		{
			int tree_id=(*segregating_tree)[i];
			if(tree_id!=id)
			{
				internal_seg_tree[new_seg_tree_id][j]=tree_id;
				j++;
			}
		}
		internal_seg_tree[new_seg_tree_id][j]=id;
		--n_segregating;
		seg_tree_id=new_seg_tree_id;
		segregating_tree=&(internal_seg_tree[seg_tree_id]);
		return *this;
	}
	virtual coalescent& calc_node_rlen(class ap_node* id)
	{
		/*In this function 'seg' or 'segregating' means that the
		marginal tree at that site has not found its mrca*/
		int left,right;
		/*left is the first non-NULL seg site from the left*/
		int i;
		for(i=0;i<n_segregating;i++)
		{
			left=(*segregating_tree)[i];
			if(id->AMP.ptr[left]!=NULL)break;
		}
		/*right is the first non-NULL seg site from the right*/
		for(i=n_segregating-1;i>=0;i--)
		{
			right=(*segregating_tree)[i];
			if(id->AMP.ptr[right]!=NULL)break;
		}
		if(i<0)/*This occurs when there are no non-NULL seg sites*/
		{
			id->rlen=0.0;
			id->L=0.0;
			id->ltr=0;
			id->rtr=0;
		}
		else
		{
			/*L is the total number of sites bounded by non-NULL seg sites*/
			id->L = (double)(right - left + 1);
			if(no_gene_conversion)
				id->rlen = con->r * (double)(id->L - 1);
			else {
				id->rlen = 0.5 * con->r * ((double)(id->L - 1) + 1./con->lambda * (1. - pow(1.-con->lambda,(double)(id->L-1))));
				if(id->flag == ap_node::FIXED_NODE && id->L>0) id->rlen += 0.5 * con->r/con->lambda*pow(1.-con->lambda,(double)(id->L-1));
				/* Before 11.08.06, 
				id->rlen = con->r * (double)(id->L - 1);
				if(id->flag == ap_node::FIXED_NODE && id->L>0) id->rlen += con->r/con->lambda*pow(1.-con->lambda,(double)(id->L-1));*/
				/* Even earlier, id->rlen = con->r*con->lambda*(1.-pow(1.-1./con->lambda,(double)(id->L-1))); */			
			}
//					+ (con->r)/(con->lambda) * exp(-con->lambda * (double)id->L);
//			if((id->L==0.0)&&(con->r!=0.0)) 
//			{
//				error("calc_rlen: non-zero rlen value when L=0");
//			}
			id->ltr = left;
			id->rtr = right+1;
		}

		return *this;
	}
	coalescent& calc_rlen()
	{
		total_rlen=0.0;
		int i;
		for(i=0;i<ARG_k;i++)
			total_rlen+=active_node[i]->rlen;
		rho=2.0*total_rlen/(double)ARG_k;

		return *this;
	}
	coalescent& calc_rlen(class ap_node* id)
	{
		calc_node_rlen(id);
		total_rlen=0.0;
		int i;
		for(i=0;i<ARG_k;i++)
			total_rlen+=active_node[i]->rlen;
		rho=2.0*total_rlen/(double)ARG_k;

		return *this;
	}
	coalescent& calc_rlen(class ap_node* id1, class ap_node* id2)
	{
		calc_node_rlen(id1);
		calc_node_rlen(id2);
		total_rlen=0.0;
		int i;
		for(i=0;i<ARG_k;i++)
			total_rlen+=active_node[i]->rlen;
		rho=2.0*total_rlen/(double)ARG_k;

		return *this;
	}
	coalescent& migrate_calc_rlen()
	{
		total_rlen=0.0;
		int i;
		for(i=0;i<con->ndemes;i++) rho_deme[i] = 0.0;
		for(i=0;i<ARG_k;i++) {
			total_rlen += active_node[i]->rlen;
			rho_deme[active_node[i]->deme] += active_node[i]->rlen;
		}
		rho=2.0*total_rlen/(double)ARG_k;
		for(i=0;i<con->ndemes;i++) rho_deme[i] *= 2.0 / con->N_deme_over_D[i];

		return *this;
	}
	coalescent& migrate_calc_rlen(class ap_node* id)
	{
		calc_node_rlen(id);
		total_rlen=0.0;
		int i;
		for(i=0;i<con->ndemes;i++) rho_deme[i] = 0.0;
		for(i=0;i<ARG_k;i++) {
			total_rlen += active_node[i]->rlen;
			rho_deme[active_node[i]->deme] += active_node[i]->rlen;
		}
		rho=2.0*total_rlen/(double)ARG_k;
		for(i=0;i<con->ndemes;i++) rho_deme[i] *= 2.0 / con->N_deme_over_D[i];

		return *this;
	}
	coalescent& migrate_calc_rlen(class ap_node* id1, class ap_node* id2)
	{
		calc_node_rlen(id1);
		calc_node_rlen(id2);
		total_rlen=0.0;
		int i;
		for(i=0;i<con->ndemes;i++) rho_deme[i] = 0.0;
		for(i=0;i<ARG_k;i++) {
			total_rlen+=active_node[i]->rlen;
			rho_deme[active_node[i]->deme] += active_node[i]->rlen;
		}
		rho=2.0*total_rlen/(double)ARG_k;
		for(i=0;i<con->ndemes;i++) rho_deme[i] *= 2.0 / con->N_deme_over_D[i];

		return *this;
	}
	coalescent& perform_single_crossover(int *ltr, int *rtr)
	{
		/*ltr and rtr are modified so that when they are		*/
		/*returned they dictate the recombination boundaries	*/

		/*Fragment length, 1<=f<=L-1, where L=rtr-ltr			*/
		/*Simulate X=(f-1) Truncated exponential, mean=1/lambda,*/
		/*truncation point=L-1 s.t. Pr(L-1)=0.					*/

		//printf("\nLTR: %3d RTR: %3d ",*ltr,*rtr);
		int L = (*rtr)-(*ltr);
		int X;
		if(no_gene_conversion) X = ran->discrete(0,L-2);
		else
		{
			//double mean = 1.0/(con->lambda);	
			//X = floor(ran->trunc_exponential(mean,L-1));
			X = ran->trunc_geometric(con->lambda,L-1) - 1;
		}

		/*Implement it according to direction*/	
		bool dir = ran->bernoulliTF(0.5);
		if (dir) /*left to right*/
		{
			(*rtr) = (*ltr) + (X+1);
		}
		else /*right to left*/
		{
			(*ltr) = (*rtr) - (X+1);
		}
		
		//printf("X: %3d Dir: %d LTR: %3d RTR: %3d\n",X,dir,*ltr,*rtr);

		return *this;
	}
	coalescent& perform_double_crossover(int *ltr, int *rtr)
	{
		/*First part as for single cross-over					*/
		int L = (*rtr)-(*ltr);
		//double mean = 1.0/(con->lambda);
		//int X = floor(ran->trunc_exponential(mean,L-1));
		int X = ran->trunc_geometric(con->lambda,L-1) - 1;

		/*Implement it according to direction*/	
		int dir = ran->discrete(0,1);
		if (dir) /*left to right*/
		{
			(*ltr) += X+1;		/*except now this is ltr*/
		}
		else /*right to left*/
		{
			(*rtr) -= (X+1);	/*and this is now rtr	*/
		}

		/*And repeat with re-defined boundaries					*/
		L = (*rtr)-(*ltr);
		//X = floor(ran->trunc_exponential(mean,L-1));
		X = ran->trunc_geometric(con->lambda,L-1) - 1;

		if (dir) /*left to right*/
		{
			(*rtr) = (*ltr)+X+1;
		}
		else /*right to left*/
		{
			(*ltr) = (*rtr)-X-1;
		}

		return *this;
	}
	coalescent& mutate_tree(const int tid, const int nid, const int state)
	{
		int my_state=mutate_edge(state,tree[tid].node[nid].edge_time);
		if(tree[tid].node[nid].descendant[0]==NULL)
		{	/*I.e. a base node*/
			if(tree[tid].node[nid].descendant[0]!=NULL)error("mutate_tree() err1: node has one, not two descendants");
			genotype.element[nid][tid]=my_state;
		}
		else
		{
			if(tree[tid].node[nid].descendant[1]==NULL)error("mutate_tree() err2: node has one, not two descendants");

			int desc=tree[tid].node[nid].descendant[0]->id;
			mutate_tree(tid, desc, my_state);

			desc=tree[tid].node[nid].descendant[1]->id;
			mutate_tree(tid, desc, my_state);
		}
		return *this;
	}
	coalescent& mutate_tree(const int tid, const int nid, const int state, Mutation_Matrix *M)
	{
		int my_state=mutate_edge(state,tree[tid].node[nid].edge_time,M);
		if(tree[tid].node[nid].descendant[0]==NULL)
		{	/*I.e. a base node*/
			if(tree[tid].node[nid].descendant[0]!=NULL)error("mutate_tree() err1: node has one, not two descendants");
			genotype.element[nid][tid]=my_state;
		}
		else
		{
			if(tree[tid].node[nid].descendant[1]==NULL)error("mutate_tree() err2: node has one, not two descendants");

			int desc=tree[tid].node[nid].descendant[0]->id;
			mutate_tree(tid, desc, my_state, M);

			desc=tree[tid].node[nid].descendant[1]->id;
			mutate_tree(tid, desc, my_state, M);
		}
		return *this;
	}
	coalescent& mutate_tree_and_record(const int tid, const int nid, const int state, Mutation_Matrix *M, vector<int> &mutLog)
	{
		int my_state=mutate_edge_and_record(state,tree[tid].node[nid].edge_time,M,mutLog);
		if(tree[tid].node[nid].descendant[0]==NULL)
		{	/*I.e. a base node*/
			if(tree[tid].node[nid].descendant[0]!=NULL)error("mutate_tree_and_record() err1: node has one, not two descendants");
			genotype.element[nid][tid]=my_state;
		}
		else
		{
			if(tree[tid].node[nid].descendant[1]==NULL)error("mutate_tree_and_record() err2: node has one, not two descendants");

			int desc=tree[tid].node[nid].descendant[0]->id;
			mutate_tree_and_record(tid, desc, my_state, M, mutLog);

			desc=tree[tid].node[nid].descendant[1]->id;
			mutate_tree_and_record(tid, desc, my_state, M, mutLog);
		}
		return *this;
	}
	int mutate_mrca()
	{
		double rnum1 = ran->U();
		int i;
		for(i=0;i<con->n_states;i++)
		{
			if (rnum1<=con->state_freq[i]) break;
			rnum1 -= con->state_freq[i];
		}
		if (i>=con->n_states) error("mutate_mrca(): initial state chosen incorrectly");
		//It's important state_freq sums to one
		
		return i;
	}
	int mutate_edge(const int state, const double edge_time)
	{
		int gt=state;
		if (gt==-1)
		{
			error("mutate_edge(): genotype does not exist");
		}
		double time_left=edge_time;
		double next_mut=ran->exponential(con->state_M[gt]);
		while (next_mut<time_left)
		{
			double rnum1=ran->U();
			int i;
			for(i=0;i<con->n_states;i++)
			{
				if(rnum1<=con->mut_matrix.element[gt][i])
					break;
				rnum1-=con->mut_matrix.element[gt][i];
			}
			if(i>=con->n_states) 
			{
				printf("\nCurrent state: %d Random uniform[0,1] deviate: %g",gt,rnum1);
				printf("\nError by-passed, mutation anulled\n");
				//error("mutate_edge(): transition incorrectly chosen");
				i=gt;
			}

			//node[st_node]->genotype[site]=i;
			gt=i;
			++nmut;

			//++mutation_spectrum[count_descendants(site,st_node)];

			time_left-=next_mut;
			next_mut=ran->exponential(con->state_M[gt]);
		}
		
		/*Make descendants inherit state*/
		//inherit(site,st_node);
		return gt;
	}
	int mutate_edge(const int state, const double edge_time, Mutation_Matrix *M)
	{
		int gt=state;
		if (gt<0||gt>=M->n_states)
		{
			error("mutate_edge(): genotype does not exist");
		}
		double time_left=edge_time;
		double next_mut=ran->exponential(M->mutation_mean[gt]);
		double rp;
		while (next_mut<time_left)
		{
			double rnum1=ran->U();
			int i;
			for(i=0;i<M->n_states-1;i++)
			{
				rp = M->D[gt][i];
				if(rnum1 <= rp)
					break;
				rnum1 -= rp;
			}

			gt=i;
			++nmut;

			time_left-=next_mut;
			next_mut=ran->exponential(M->mutation_mean[gt]);
		}
		return gt;
	}
	int mutate_edge_and_record(const int state, const double edge_time, Mutation_Matrix *M, vector<int> &mutLog)
	{
		int gt=state;
		if (gt<0||gt>=M->n_states)
		{
			error("mutate_edge(): genotype does not exist");
		}
		double time_left=edge_time;
		double next_mut=ran->exponential(M->mutation_mean[gt]);
		double rp;
		while (next_mut<time_left)
		{
			mutLog.push_back(gt);
			double rnum1=ran->U();
			int i;
			for(i=0;i<M->n_states-1;i++)
			{
				rp = M->D[gt][i];
				if(rnum1 <= rp)
					break;
				rnum1 -= rp;
			}

			gt=i;
			++nmut;
			mutLog.push_back(gt);

			time_left-=next_mut;
			next_mut=ran->exponential(M->mutation_mean[gt]);
		}
		return gt;
	}

public:
	/* Number of segregating sites */
	double S() {
		double result = 0.0;
		if(con->nsamp==0) return 0.0;

		int i,j;
		for(j=0;j<con->seq_len;j++) {
			double hap = genotype[0][j];
			for(i=1;i<con->nsamp;i++)
				if(genotype[i][j]!=hap) {
					++result;
					break;
				}
		}

		return result;
	}

	/* Number of unique haplotypes */
	double H() {
		int result = 1;

		if(con->nsamp==0) return 0.0;
		_uniqueHaps = vector<int>(con->nsamp,-1);
		_uniqueHaps[0] = 0;

		int i,ii,j;
		bool unique;
		for(i=1;i<con->nsamp;i++) {
			unique = true;
			for(ii=0;ii<result;ii++) {
				for(j=0;j<con->seq_len;j++) {
					if(genotype[i][j]!=genotype[_uniqueHaps[ii]][j]) break;
				}
				if(j==con->seq_len) unique = false;
			}
			if(unique==true) {
				_uniqueHaps[result] = i;
				++result;
			}
		}
		return (double)result;
	}

	/* Average number of pairwise differences */
	double pi() {
		double result = 0.0;

		int i,j,k;
		for(i=0;i<con->nsamp;i++)
			for(j=0;j<i;j++)
				for(k=0;k<con->seq_len;k++)
					result += (genotype[i][k]==genotype[j][k]) ? 0.0 : 1.0;
		result *= 2.0/(double)(con->nsamp)/(double)(con->nsamp-1);

		return result;
	}

	/* Variance in number of pairwise differences */
	double Varpi() {
		double E,EE,pi;
		int i,j,k;

		E = EE = 0.0;
		for(i=0;i<con->nsamp;i++)
			for(j=0;j<i;j++) {
				pi = 0.0;
				for(k=0;k<con->seq_len;k++)
					pi += (genotype[i][k]==genotype[j][k]) ? 0.0 : 1.0;
				E += pi;
				EE += pi*pi;
			}
		E *= 2.0/(double)(con->nsamp)/(double)(con->nsamp-1);
		EE *= 2.0/(double)(con->nsamp)/(double)(con->nsamp-1);

		double result = EE - E*E;
		return result;
	}

	double Tajima() {
		double D = 0.0;
		int i,j,k,n,L;
		n = con->nsamp;
		L = con->seq_len;
		double a1,a2,b1,b2,c1,c2,e1,e2,khat,S;
		bool segregating;
		khat = S = 0.0;
		for(k=0;k<L;k++) {
			segregating = false;
			for(i=0;i<n;i++)
				for(j=0;j<i;j++)
					if(genotype[i][k]!=genotype[j][k]) {
						++khat;
						segregating = true;
					}
			if(segregating) ++S;
		}
		if(S==0) return 0.0;
		khat /= (double)(n*(n-1)/2);
		a1 = a2 = 0.0;
		for(i=1;i<=n-1;i++) {
			a1 += 1./(double)i;
			a2 += 1./(double)(i*i);
		}
		b1 = (double)(n+1)/(double)(3*(n-1));
		b2 = (double)(2*(n*n+n+3))/(double)(9*n*(n-1));
		c1 = b1 - 1./a1;
		c2 = b2 - (double)(n+2)/a1/(double)(n) + a2/a1/a1;
		e1 = c1/a1;
		e2 = c2/(a1*a1+a2);
		D = (khat - S/a1)/sqrt(e1*S+e2*S*(S-1.));
		return D;
	}
	/* This function counts the average number of pairwise differences, where the matrix
	   diff defines those differences using 0 = different or 1 = identical.
	   If diff is the identity matrix then this function is equivalent to double pi(). */
	double pi(Matrix<int> &diff) {
		double result = 0.0;

		int i,j,k;
		for(i=0;i<con->nsamp;i++)
			for(j=0;j<i;j++)
				for(k=0;k<con->seq_len;k++)
					result += (diff[(int)genotype[i][k]][(int)genotype[j][k]]==0);
		result *= 2.0/(double)(con->nsamp)/(double)(con->nsamp-1);

		return result;
	}
	/* Hudson and Kaplan's Rm, the minimum # recombinations. See Myers and Griffiths(2003)*/
	double Rm() {
		if(con->nsamp==0) return 0.0;
		if(con->seq_len==0) return 0.0;

		/* Determine which sites are biallelic segregating */
		_sites = vector<int>(con->seq_len,0);
		int i,j,k;
		int S = 0;
		double hap0,hap1;
		bool segregating;
		for(j=0;j<con->seq_len;j++) {
			segregating = false;
			hap0 = genotype[0][j];
			for(i=1;i<con->nsamp;i++) {
				if(!segregating && genotype[i][j]!=hap0) {
					segregating = true;
					hap1 = genotype[i][j];
				}
				else if(segregating && genotype[i][j]!=hap0 && genotype[i][j]!=hap1) {
					segregating = false;	// define segregating only for biallelic sites
					break;
				}
			}
			if(segregating) {
				_sites[S] = j;
				++S;
			}
		}
		if(S<2) return 0.0;

		/* Calculate the compatibility matrix */
		____B = LowerTriangularMatrix<int>(S,0);	// so j>=k always
		// ____B[j][k] = 0 for compatible, 1 for incompatible
		bool comb[3];
		for(j=0;j<S;j++)
			for(k=0;k<j;k++)
			{
				hap0 = genotype[0][_sites[j]];
				hap1 = genotype[0][_sites[k]];
				comb[0] = false;				// hap0  hap1'
				comb[1] = false;				// hap0' hap1
				comb[2] = false;				// hap0' hap1'
				for(i=1;i<con->nsamp;i++) {
					if(genotype[i][_sites[j]]==hap0 && genotype[i][_sites[k]]!=hap1) comb[0] = true;
					if(genotype[i][_sites[j]]!=hap0 && genotype[i][_sites[k]]==hap1) comb[1] = true;
					if(genotype[i][_sites[j]]!=hap0 && genotype[i][_sites[k]]!=hap1) comb[2] = true;
					if(comb[0] && comb[1] && comb[2]) break;
				}
				____B[j][k] = (comb[0] && comb[1] && comb[2]) ? 1 : 0;			
			}

		/* Calculate the dynamic programming partition matrix */
		_M = vector<int>(S,0);
		int maxM = 0;
		_M[S-1] = 0;
		_M[S-2] = ____B[S-1][S-2];
		for(i=S-3;i>=0;i--) {
			_M[i] = ____B[i+1][i] + _M[i+1];
			for(k=i+2;k<S;k++) if(____B[k][i]+_M[k]>_M[i]) _M[i] = ____B[k][i]+_M[k];
		}

		return (double)_M[0];
	}
	void RecCorrelations(double* result) {
		RecCorrelations(result,true);
	}
	void RecCovariances(double* result) {
		RecCorrelations(result,false);
	}
	void RecCorrelations(double* result, bool normalize) {
		result[0] = result[1] = result[2] = 0.0;

		if(con->nsamp==0) return;
		if(con->seq_len==0) return;

		/* Determine which sites are biallelic segregating */
		_sites = vector<int>(con->seq_len,0);
		int i,j,k;
		int S = 0;
		double hap0,hap1;
		bool segregating;
		for(j=0;j<con->seq_len;j++) {
			segregating = false;
			hap0 = genotype[0][j];
			for(i=1;i<con->nsamp;i++) {
				if(!segregating && genotype[i][j]!=hap0) {
					segregating = true;
					hap1 = genotype[i][j];
				}
				else if(segregating && genotype[i][j]!=hap0 && genotype[i][j]!=hap1) {
					segregating = false;	// define segregating only for biallelic sites
					break;
				}
			}
			if(segregating) {
				_sites[S] = j;
				++S;
			}
		}
		if(S<3) return;
		
		/* Calculate frequency statistics */
		_F = vector<double>(S,1.0);							/* _F is the marginal frequency of hap0 at site j */
		for(j=0;j<S;j++) {
			hap0 = genotype[0][_sites[j]];
			for(i=1;i<con->nsamp;i++)
				if(genotype[i][_sites[j]]==hap0) _F[j]++;
			_F[j] /= (double)con->nsamp;
		}

		_four = vector<double>(4,0.0);							/* _G[j][k] is the frequency of AB (_G[j][k][0]),	*/
		_G = LowerTriangularMatrix< vector<double> >(S,_four);	/* Ab (1), aB (2), ab (3) for sites j and k		*/
		for(j=0;j<S;j++)
		  for(k=0;k<j;k++) {
			  hap0 = genotype[0][_sites[j]];
			  hap1 = genotype[0][_sites[k]];
			  for(i=0;i<con->nsamp;i++) {
				  if(genotype[i][_sites[j]]==hap0 && genotype[i][_sites[k]]==hap1) ++_G[j][k][0];
				  else if(genotype[i][_sites[j]]==hap0 && genotype[i][_sites[k]]!=hap1) ++_G[j][k][1];
				  else if(genotype[i][_sites[j]]!=hap0 && genotype[i][_sites[k]]==hap1) ++_G[j][k][2];
				  else if(genotype[i][_sites[j]]!=hap0 && genotype[i][_sites[k]]!=hap1) ++_G[j][k][3];
				  else warning("Unexpected choice");
			  }
			  for(i=0;i<4;i++) _G[j][k][i] /= (double)con->nsamp;
		  }
  
		/* Calculate LD statistics for pairs of sites */
		_A = LowerTriangularMatrix<double>(S,0.0);			//	rsq
		___B = LowerTriangularMatrix<double>(S,0.0);			//	Dprime
		___C = LowerTriangularMatrix<double>(S,0.0);			//	G4
		_D = Matrix<double>(S,S,0.0);

		double temp;
		for(i=0;i<S;i++) {
			for(j=0;j<i;j++) {
				temp = _G[i][j][0] - _F[i]*_F[j];
				_A[i][j] = pow(temp,2.0)/(_F[i]*(1.-_F[i])*_F[j]*(1.-_F[j]));
				___B[i][j] = (temp < 0.0) ? -temp/MIN(_F[i]*_F[j],(1.-_F[i])*(1.-_F[j])) : temp/MIN(_F[i]*(1.-_F[j]),(1.-_F[i])*_F[j]);
				___C[i][j] = (_G[i][j][0]>0.0 && _G[i][j][1]>0.0 && _G[i][j][2]>0.0 && _G[i][j][3]>0.0) ? 1.0 : 0.0;
				_D[i][j] = _D[j][i] = _sites[i] - _sites[j];
			}
		}

		double  E[4] = {0.0,0.0,0.0,0.0};
		double EE[4] = {0.0,0.0,0.0,0.0};
		double ED[3] = {0.0,0.0,0.0};
		int ctr;
//		ofstream out("ld.txt");
		for(i=0,ctr=0;i<S;i++)
			for(j=0;j<i;j++,ctr++) {
//				out << _A[i][j] << "\t" << _D[i][j] << endl;
				E[0] += _A[i][j]; E[1] += ___B[i][j]; E[2] += ___C[i][j]; E[3] += _D[i][j];
				EE[0] += _A[i][j]*_A[i][j]; EE[1] += ___B[i][j]*___B[i][j]; EE[2] += ___C[i][j]*___C[i][j]; EE[3] += _D[i][j]*_D[i][j];
				ED[0] += _A[i][j]*_D[i][j]; ED[1] += ___B[i][j]*_D[i][j]; ED[2] += ___C[i][j]*_D[i][j];
			}
//		out.close();

		if(normalize)	// Calculate correlation
			for(k=0;k<3;k++)
				result[k] = (ED[k]-E[k]*E[3]/(double)(ctr))/sqrt((EE[k]-E[k]*E[k]/(double)(ctr))*(EE[3]-E[3]*E[3]/(double)(ctr)));
		else			// Calculate covariance
			for(k=0;k<3;k++)
				result[k] = (ED[k]-E[k]*E[3]/(double)(ctr))/(double)(ctr);

		return;
	}
};

#endif // ___COALESCENT_PROCESS_H_
