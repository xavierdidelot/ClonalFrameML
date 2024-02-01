/*
 *  main.h
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
#ifndef _MAIN_H_
#define _MAIN_H_
#include <iostream>
#include <string.h>
#include "myutils/newick.h"
#include "coalesce/coalescent_record.h"
#include <sstream>
#include "xmfa.h"
#include <fstream>
#include <algorithm>
#include "myutils/DNA.h"
#include "myutils/mydouble.h"
#include "powell.h"
#include "myutils/argumentwizard.h"
#include <time.h>
#include "myutils/random.h"
#include <limits>
#include <iomanip>
#define ClonalFrameML_version "v1.13"

using std::cout;
using myutils::NewickTree;
using std::stringstream;
using myutils::error;
using myutils::ArgumentWizard;
using myutils::DATA_TYPE;

// Global definition of random number generator
Random ran;

enum Nucleotide {Adenine=0, Guanine, Cytosine, Thymine, N_ambiguous};
enum ImportationState {Unimported=0, Imported};

marginal_tree convert_rooted_NewickTree_to_marginal_tree(NewickTree &newick, vector<string> &tip_labels, vector<string> &all_node_labels);
marginal_tree convert_unrooted_NewickTree_to_marginal_tree(NewickTree &newick, vector<string> &tip_labels, vector<string> &all_node_labels);
vector<int> compute_compatibility(DNA &fa, marginal_tree &tree, vector<bool> &anyN, bool purge_singletons=true);
NewickTree read_Newick(const char* newick_file);
Matrix<Nucleotide> FASTA_to_nucleotide(DNA &fa, vector<double> &empirical_nucleotide_frequencies, vector<bool> usesite);
void find_alignment_patterns(Matrix<Nucleotide> &nuc, vector<bool> &iscompat, vector<string> &pat, vector<int> &pat1, vector<int> &cpat, vector<int> &ipat);
vector< Matrix<double> > compute_HKY85_ptrans(const marginal_tree &ctree, const double kappa, const vector<double> &pi);
Matrix<mydouble> compute_HKY85_ptrans(const double x, const double k, const vector<double> &pi);
Matrix<double> dcompute_HKY85_ptrans(const double x, const double kappa, const vector<double> &pi);
double HKY85_expected_rate(const vector<double> &n, const double kappa, const vector<double> &pi);
mydouble maximum_likelihood_ancestral_sequences(Matrix<Nucleotide> &nuc, marginal_tree &ctree, const double kappa, const vector<double> &pi, vector<int> &pat1, vector<int> &cpat, Matrix<Nucleotide> &node_sequence);
void write_newick(const marginal_tree &ctree, const vector<string> &all_node_names, const char* file_name);
void write_newick(const marginal_tree &ctree, const vector<string> &all_node_names, ofstream &fout);
void write_newick_node(const mt_node *node, const vector<string> &all_node_names, ofstream &fout);
void write_ancestral_fasta(Matrix<Nucleotide> &nuc, vector<string> &all_node_names, const char* file_name);
void write_filtered_fasta(vector< vector<ImportationState> > &imported, DNA * fa,vector<bool> & ignore_site, const char* file_name);
void write_position_cross_reference(vector<bool> &iscompat, vector<int> &ipat, const char* file_name);
void write_position_cross_reference(vector<bool> &iscompat, vector<int> &ipat, ofstream &fout);
mydouble likelihood_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<int> &pat1, const vector<int> &cpat, const double kappa, const vector<double> &pinuc, const double branch_length);
bool string_to_bool(const string s, const string label="");
void write_importation_status_intervals(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, const char* file_name, const int root_node,const char* chr_name);
double Baum_Welch(const marginal_tree &tree, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const vector<bool> &informative, const vector<double> &prior_a, const vector<double> &prior_b, vector<double> &full_param, vector<double> &posterior_a, int &neval, const bool coutput, double &priorL);
double Baum_Welch0(const marginal_tree &tree, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const vector<bool> &informative, const vector<double> &prior_a, const vector<double> &prior_b, const vector<double> &full_param, const vector<double> &posterior_a, const bool coutput);
double gamma_loglikelihood(const double x, const double a, const double b);
Matrix<double> Baum_Welch_simulate_posterior(const marginal_tree &tree, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const vector<bool> &informative, const vector<double> &prior_a, const vector<double> &prior_b, const vector<double> &full_param, int &neval, const bool coutput, const int nsim);
double Baum_Welch_Rho_Per_Branch(const marginal_tree &tree, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const vector<bool> &informative, const vector<double> &prior_a, const vector<double> &prior_b, vector<double> &mean_param, Matrix<double> &full_param, Matrix<double> &posterior_a, int &neval, const bool coutput);
mydouble maximum_likelihood_ClonalFrame_branch_allsites(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pi, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, vector<ImportationState> &is_imported);

class orderNewickNodesByStatusLabelAndAge {
public:
  using first_argument_type = size_t;
  using second_argument_type = size_t;
  using result_type = bool;
	
	const vector<NewickNode*> &root2tip;	// temporary ordering of Newick nodes from root to tips
	const vector<double> &ageroot2tip;		// corresponding age of each node in root2tip
	const vector<size_t> &labelorder;		// The position where each node comes in the label order (for tips; the label is ignored for internal nodes)
	orderNewickNodesByStatusLabelAndAge(const vector<NewickNode*> &root2tip_in, const vector<double> &ageroot2tip_in, const vector<size_t> &labelorder_in) :
	root2tip(root2tip_in), ageroot2tip(ageroot2tip_in), labelorder(labelorder_in) {
	}
	// Test if i is less than j
	bool operator()(size_t i, size_t j) const {
		if(root2tip[i]->dec.size()==0 && root2tip[j]->dec.size()!=0) {
			// If i is a tip and j is not
			return true;
		} else if(root2tip[i]->dec.size()==0 && root2tip[j]->dec.size()==0) {
			// If i and j are both tips
			// Then order by label
			if(labelorder[i]==labelorder[j]) {
				stringstream errTxt;
				errTxt << "orderNewickNodesByStatusLabelAndAge::operator(): ";
				errTxt << "tips cannot have the same label order";
				error(errTxt.str().c_str());
			}
			return labelorder[i] < labelorder[j];
		} else if(root2tip[i]->dec.size()!=0 && root2tip[j]->dec.size()==0) {
			// If i is not a tip but j is
			return false;
		} else {
			// If neither are tips
			// Then order by age
			return ageroot2tip[i] < ageroot2tip[j];
		}
	}
};

class ClonalFrameRescaleBranchFunction : public PowellFunction {
public:
	// References to non-member variables
	const mt_node &node;
	const Matrix<Nucleotide> &node_nuc;
	const vector<int> &pat1;
	const vector<int> &cpat;
	const double kappa;
	const vector<double> &pi;
	// True member variable
	mydouble ML;
	int neval;
	const bool multithread;
	double crude_branch_length;
	double min_branch_length;
public:
	ClonalFrameRescaleBranchFunction(const mt_node &_node, const Matrix<Nucleotide> &_node_nuc, const vector<int> &_pat1, const vector<int> &_cpat, const double _kappa,
									const vector<double> &_pi, const bool _multithread, const double _crude_branch_length, const double _min_branch_length) :
	node(_node), node_nuc(_node_nuc), pat1(_pat1), cpat(_cpat), kappa(_kappa), pi(_pi), neval(0),
	multithread(_multithread), crude_branch_length(_crude_branch_length), min_branch_length(_min_branch_length) {};
	double f(const vector<double>& x) {
		++neval;
		// Process parameters
		if(!(x.size()==1)) error("ClonalFrameRescaleBranchFunction::f(): 1 argument required");
		double branch_length = pow(10.,x[0]);
		if(branch_length<min_branch_length) branch_length = min_branch_length;
		const int dec_id = node.id;
		const int anc_id = node.ancestor->id;
		// Calculate likelihood
		ML = likelihood_branch(dec_id,anc_id,node_nuc,pat1,cpat,kappa,pi,branch_length);
		return -ML.LOG();
	}
};

/*	Maximum likelihood routine based on the Baum-Welch EM algorithm for estimating
	a single set of recombination parameters (R/M, import length, import divergence)
	and an independent branch length per branch. Note that the approach is classical
	and the priors act through pseudocounts - i.e. a form of data augmentation prior */
class ClonalFrameBaumWelch {
public:
	// References to non-member variables
	const marginal_tree &tree;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	vector< vector<ImportationState> > &is_imported;
	// True member variable
	double ML,ML0,priorL;
	double PR;
	int neval;
	const vector<double> prior_a;
	const vector<double> prior_b;
	vector<double> which_compat;
	const int root_node;
	vector<bool> informative;
	vector<double> initial_branch_length;
	vector<double> full_param;
	vector<double> posterior_a;
	bool guess_initial_m;
	bool coutput;
public:
	ClonalFrameBaumWelch(const marginal_tree &_tree, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
							   const vector<double> &_pi, vector< vector<ImportationState> > &_is_imported,
							   const vector<double> &_prior_a, const vector<double> &_prior_b, const int _root_node, const bool _guess_initial_m, const bool _coutput=false) :
	tree(_tree), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa),
	pi(_pi), neval(0), is_imported(_is_imported),
	prior_a(_prior_a), prior_b(_prior_b), root_node(_root_node), initial_branch_length(_root_node), informative(_root_node), guess_initial_m(_guess_initial_m),
	coutput(_coutput) {
		if(prior_a.size()!=4) error("ClonalFrameBaumWelch: prior a must have length 4");
		if(prior_b.size()!=4) error("ClonalFrameBaumWelch: prior b must have length 4");
		int i;
		// Precompute which sites are compatible
		which_compat = vector<double>(0);
		for(i=0;i<iscompat.size();i++) {
			if(iscompat[i]) {
				which_compat.push_back((double)i);
			}
		}
		int j,k;
		for(i=0;i<root_node;i++) {
			// Crudely re-estimate branch length: use this as the mean of the prior on branch length ????
			double pd = 1.0, pd_den = 2.0;
			const int dec_id = tree.node[i].id;
			const int anc_id = tree.node[i].ancestor->id;
			for(j=0,k=0;j<iscompat.size();j++) {
				if(iscompat[j]) {
					Nucleotide dec = node_nuc[dec_id][ipat[k]];
					Nucleotide anc = node_nuc[anc_id][ipat[k]];
					if(dec!=anc) ++pd;
					++pd_den;
					++k;
				}
			}
			initial_branch_length[i] = pd/pd_den;
//			initial_branch_length[i] = tree.node[i].edge_time;
			informative[i] = (pd>=2.0) ? true : false;
		}
	}
	vector<double> maximize_likelihood(const vector<double> &param) {
		if(!(param.size()==3)) error("ClonalFrameBaumWelch::maximize_likelihood(): 3 arguments required");
		// Starting points for the shared parameters
		full_param = vector<double>(0);
		posterior_a = vector<double>(0);
		full_param.push_back(param[0]);		// rho_over_theta
		full_param.push_back(param[1]);		// mean_import_length: may need to invert
		full_param.push_back(param[2]);		// import_divergence
		int i;
		for(i=0;i<initial_branch_length.size();i++) {
			// Standard approach: initially equate the expected number of mutations and substitutions
			double ibl = initial_branch_length[i];
			if(guess_initial_m) {
				// Alternative approach: apportion the expected number of mutations and substitutions proportionally among M and nu according to the prior
				// E(S) = 1/(1+R*delta)*M + R*delta/(1+R*delta)*nu
				//      = M * ( 1/(1+M*R/M*delta) + R/M*delta*nu/(1+M*R/M*delta) ) = M * ( 1 + nu*R/M*delta )/( 1 + M*R/M*delta )
				// So let, and possibly iterate a few times
				// M = S * ( 1 + M*R/M*delta )/( 1 + nu*R/M*delta )
				// log(M) = log(S) + log(1+10^(logM+logR/M+logdelta)) - log(1+10^(lognu+logR/M+logdelta))
				double log10ibl = log10(ibl);
				int j;
				for(j=0;j<3;j++) log10ibl = log10(initial_branch_length[i]) + log(1.0+pow(10.,log10ibl+param[0]+param[1])) - log(1.0+pow(10.,param[2]+param[0]+param[1]));
				// Make sure it hasn't gone wrong for some reason before substituting...
				if(!(log10ibl!=log10ibl)) ibl = pow(10.,log10ibl);
			}
			full_param.push_back(ibl);
		}
		// Iterate
		ML = Baum_Welch(tree,node_nuc,which_compat,ipat,kappa,pi,informative,prior_a,prior_b,full_param,posterior_a,neval,coutput,priorL);
		// Update importation status for all branches **for ALL SITES**, including uninformative ones
		for(i=0;i<initial_branch_length.size();i++) {
			const int dec_id = tree.node[i].id;
			const int anc_id = tree.node[i].ancestor->id;
			const double rho_over_theta = full_param[0];
			const double mean_import_length = full_param[1];
			const double import_divergence = full_param[2];
			const double branch_length = (informative[i]) ? full_param[3+i] : initial_branch_length[i];
			maximum_likelihood_ClonalFrame_branch_allsites(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,import_divergence,is_imported[i]);
		}
		ML0 = Baum_Welch0(tree,node_nuc,which_compat,ipat,kappa,pi,informative,prior_a,prior_b,full_param,posterior_a,coutput);
		return full_param;
	}
	Matrix<double> simulate_posterior(const vector<double> &param, const int nsim) {
		if(!(param.size()==3+informative.size())) error("ClonalFrameBaumWelch::simulate_posterior(): 3 arguments required");
		return Baum_Welch_simulate_posterior(tree,node_nuc,which_compat,ipat,kappa,pi,informative,prior_a,prior_b,param,neval,coutput,nsim);
	}
};

/*	In this version, the Baum-Welch algorithm is used to maximize the likelihood of
	all four parameters (R/M, import length, import divergence, branch length)
	for each branch. As for ClonalFrameBaumWelch, the prior acts through pseudocounts
	i.e. a data augmentation prior, and there is an extra parameter whose prior
	determines the variance in estimates of the recombination parameters per branch.
	This parameter needs to be set fairly stringently to prevent wild estimates in
	the absence of strong information per branch.									*/
class ClonalFrameBaumWelchRhoPerBranch {
public:
	// References to non-member variables
	const marginal_tree &tree;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	vector< vector<ImportationState> > &is_imported;
	// True member variable
	double ML;
	double PR;
	int neval;
	const vector<double> prior_a;
	const vector<double> prior_b;
	vector<double> which_compat;
	const int root_node;
	vector<bool> informative;
	vector<double> initial_branch_length;
	vector<double> mean_param;					//	Mean recombination parameters
	Matrix<double> full_param;					//	Branch-specific recombination parameters and branch length
	Matrix<double> posterior_a;
	bool guess_initial_m;
	bool coutput;
public:
	ClonalFrameBaumWelchRhoPerBranch(const marginal_tree &_tree, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
						 const vector<double> &_pi, vector< vector<ImportationState> > &_is_imported,
						 const vector<double> &_prior_a, const vector<double> &_prior_b, const int _root_node, const bool _guess_initial_m, const bool _coutput=false) :
	tree(_tree), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa),
	pi(_pi), neval(0), is_imported(_is_imported),
	prior_a(_prior_a), prior_b(_prior_b), root_node(_root_node), initial_branch_length(_root_node), informative(_root_node), guess_initial_m(_guess_initial_m),
	coutput(_coutput) {
		if(prior_a.size()!=5) error("ClonalFrameBaumWelchRhoPerBranch: prior a must have length 5");
		if(prior_b.size()!=5) error("ClonalFrameBaumWelchRhoPerBranch: prior b must have length 5");
		int i;
		// Precompute which sites are compatible
		which_compat = vector<double>(0);
		for(i=0;i<iscompat.size();i++) {
			if(iscompat[i]) {
				which_compat.push_back((double)i);
			}
		}
		int j,k;
		for(i=0;i<root_node;i++) {
			// Crudely re-estimate branch length: use this as the mean of the prior on branch length ????
			double pd = 1.0, pd_den = 2.0;
			const int dec_id = tree.node[i].id;
			const int anc_id = tree.node[i].ancestor->id;
			for(j=0,k=0;j<iscompat.size();j++) {
				if(iscompat[j]) {
					Nucleotide dec = node_nuc[dec_id][ipat[k]];
					Nucleotide anc = node_nuc[anc_id][ipat[k]];
					if(dec!=anc) ++pd;
					++pd_den;
					++k;
				}
			}
			initial_branch_length[i] = pd/pd_den;
			informative[i] = (pd>=2.0) ? true : false;
		}
	}
	void maximize_likelihood(const vector<double> &param) {
		if(!(param.size()==4)) error("ClonalFrameBaumWelchRhoPerBranch::maximize_likelihood(): 4 arguments required");
		// Starting points for the shared parameters
		mean_param = vector<double>(0);
		mean_param.push_back(param[0]);		// rho_over_theta
		// NB:- **internally** define second parameter to be INVERSE mean import length
		mean_param.push_back(1.0/param[1]);		// 1/mean_import_length
		mean_param.push_back(param[2]);		// import_divergence
		// Specially for the mean branch length, set it to the crudely estimated value assuming no recombnation
		mean_param.push_back(0.0);		// mean branch length
		int i;
		for(i=0;i<initial_branch_length.size();i++) {
			mean_param[3] += initial_branch_length[i];
		}
		mean_param[3] /= (double)initial_branch_length.size();
		full_param = Matrix<double>(initial_branch_length.size(),4);
		posterior_a = Matrix<double>(initial_branch_length.size(),4,0.0);
		for(i=0;i<initial_branch_length.size();i++) {
			full_param[i][0] = 1.0;		// factor rho_over_theta
			full_param[i][1] = 1.0;		// factor mean_import_length
			full_param[i][2] = 1.0;		// factor import_divergence
			// Standard approach: initially equate the expected number of mutations and substitutions
			double ibl = initial_branch_length[i];
			if(guess_initial_m) {
				// Alternative approach: apportion the expected number of mutations and substitutions proportionally among M and nu according to the prior
				// E(S) = 1/(1+R*delta)*M + R*delta/(1+R*delta)*nu
				//      = M * ( 1/(1+M*R/M*delta) + R/M*delta*nu/(1+M*R/M*delta) ) = M * ( 1 + nu*R/M*delta )/( 1 + M*R/M*delta )
				// So let, and possibly iterate a few times
				// M = S * ( 1 + M*R/M*delta )/( 1 + nu*R/M*delta )
				// log(M) = log(S) + log(1+10^(logM+logR/M+logdelta)) - log(1+10^(lognu+logR/M+logdelta))
				double log10ibl = log10(ibl);
				int j;
				for(j=0;j<3;j++) log10ibl = log10(initial_branch_length[i]) + log(1.0+pow(10.,log10ibl+param[0]+param[1])) - log(1.0+pow(10.,param[2]+param[0]+param[1]));
				// Make sure it hasn't gone wrong for some reason before substituting...
				if(!(log10ibl!=log10ibl)) ibl = pow(10.,log10ibl);
			}
			full_param[i][3] = ibl/mean_param[3];
		}
		// Iterate
		ML = Baum_Welch_Rho_Per_Branch(tree,node_nuc,which_compat,ipat,kappa,pi,informative,prior_a,prior_b,mean_param,full_param,posterior_a,neval,coutput);
		// Update importation status for all branches **for ALL SITES**, including uninformative ones
		for(i=0;i<initial_branch_length.size();i++) {
			const int dec_id = tree.node[i].id;
			const int anc_id = tree.node[i].ancestor->id;
			const double rho_over_theta = mean_param[0]*full_param[i][0];
			const double mean_import_length = 1.0/(mean_param[1]*full_param[i][1]);
			const double import_divergence = mean_param[2]*full_param[i][2];
			const double branch_length = (informative[i]) ? mean_param[3]*full_param[i][3] : initial_branch_length[i];
			maximum_likelihood_ClonalFrame_branch_allsites(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,import_divergence,is_imported[i]);
		}
		return;
	}
	Matrix<double> simulate_posterior(const vector<double> &param, const int nsim) {
		error("Not implemented yet");
//		if(!(param.size()==3+informative.size())) error("ClonalFrameBaumWelchRhoPerBranch::simulate_posterior(): 3 arguments required");
//		return Baum_Welch_simulate_posterior(tree,node_nuc,which_compat,ipat,kappa,pi,informative,prior_a,prior_b,param,neval,coutput,nsim);
		return Matrix<double>(0,0,0);
	}
};

#endif // _MAIN_H_
