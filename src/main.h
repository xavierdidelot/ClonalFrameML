/*  Copyright 2013 Daniel Wilson and Xavier Didelot.
 *
 *  main.h
 *  Part of ClonalFrameML
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
 *
 */
#ifndef _MAIN_H_
#define _MAIN_H_
#include <iostream>
#include <newick.h>
#include <coalesce.h>
#include <coalescent_record.h>
#include <sstream>
#include <myutils.h>
#include <fstream>
#include <algorithm>
#include <DNA.h>
#include <mydouble.h>
#include <mutation.h>
#include <powell.h>
#include <argumentwizard.h>
#include <time.h>
#include <omp.h>
#include <random.h>
#include <limits>
#include <bfgs.h>
#include <version.h>

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
vector< Matrix<double> > compute_HKY85_ptrans(const marginal_tree &ctree, const double kappa, const vector<double> &pi, const double import_divergence, const bool excess_divergence_model);
Matrix<mydouble> compute_HKY85_ptrans(const double x, const double k, const vector<double> &pi);
Matrix<double> dcompute_HKY85_ptrans(const double x, const double kappa, const vector<double> &pi);
double HKY85_expected_rate(const vector<double> &n, const double kappa, const vector<double> &pi);
mydouble maximum_likelihood_ancestral_sequences(Matrix<Nucleotide> &nuc, marginal_tree &ctree, const double kappa, const vector<double> &pi, vector<int> &pat1, vector<int> &cpat, Matrix<Nucleotide> &node_sequence);
void write_newick(const marginal_tree &ctree, const vector<string> &all_node_names, const char* file_name);
void write_newick(const marginal_tree &ctree, const vector<string> &all_node_names, ofstream &fout);
void write_newick_node(const mt_node *node, const vector<string> &all_node_names, ofstream &fout);
void write_ancestral_fasta(Matrix<Nucleotide> &nuc, vector<string> &all_node_names, const char* file_name);
void write_ancestral_fasta(Matrix<Nucleotide> &nuc, vector<string> &all_node_names, ofstream &fout);
void write_position_cross_reference(vector<bool> &iscompat, vector<int> &ipat, const char* file_name);
void write_position_cross_reference(vector<bool> &iscompat, vector<int> &ipat, ofstream &fout);
mydouble maximum_likelihood_ClonalFrame_branch_allsites(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pi, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, vector<ImportationState> &is_imported);
mydouble maximum_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pi, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, vector<ImportationState> &is_imported);
mydouble maximum_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<double> &which_compat, const vector<int> &ipat, const double kappa, const vector<double> &pi, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, vector<ImportationState> &is_imported);
double marginal_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pi, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence);
mydouble mydouble_marginal_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pi, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence);
mydouble mydouble_marginal_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<double> &which_compat, const vector<int> &ipat, const double kappa, const vector<double> &pi, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence);
mydouble likelihood_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<int> &pat1, const vector<int> &cpat, const double kappa, const vector<double> &pinuc, const double branch_length);
bool string_to_bool(const string s, const string label="");
void write_importation_status(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, const char* file_name, const int root_node);
void write_importation_status(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, ofstream &fout, const int root_node);
void write_importation_status_intervals(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, const char* file_name, const int root_node);
void write_importation_status_intervals(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, ofstream &fout, const int root_node);
void maximum_likelihood_parameters_given_path(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const vector<ImportationState> &is_imported, vector<double> &MLE);
double Viterbi_training(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, double &branch_length, double &rho_over_theta, double &mean_import_length, double &import_divergence, vector<ImportationState> &is_imported, int &neval);

class orderNewickNodesByStatusAndAge : public std::binary_function<size_t,size_t,bool> {
public:
	const vector<NewickNode*> &root2tip;	// temporary ordering of Newick nodes from root to tips
	const vector<double> &ageroot2tip;		// corresponding age of each node in root2tip
	orderNewickNodesByStatusAndAge(const vector<NewickNode*> &root2tip_in, const vector<double> &ageroot2tip_in) : root2tip(root2tip_in), ageroot2tip(ageroot2tip_in) {
	}
	// Test if i is less than j
	bool operator()(size_t i, size_t j) const {
		// If i is a tip
		if(root2tip[i]->dec.size()==0) {
			// And j is not
			if(root2tip[j]->dec.size()!=0) return true;
			// If both are tips
			return ageroot2tip[i] < ageroot2tip[j];
		} else {
			// If j is
			if(root2tip[j]->dec.size()==0) return false;
			// Neither are tips
			return ageroot2tip[i] < ageroot2tip[j];
		}
	}
};

class orderNewickNodesByStatusLabelAndAge : public std::binary_function<size_t,size_t,bool> {
public:
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

class ClonalFrameBranchFunction : public PowellFunction {
public:
	// References to non-member variables
	const mt_node &node;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	// True member variable
	vector<ImportationState> is_imported;
	mydouble ML;
	int neval;
	bool excess_divergence_model;
public:
	ClonalFrameBranchFunction(const mt_node &_node, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
							  const vector<double> &_pi, bool _excess_divergence_model) : 
	node(_node), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), pi(_pi), neval(0), excess_divergence_model(_excess_divergence_model) {};
	double f(const vector<double>& x) {
		++neval;
		if(x.size()!=4) error("ClonalFrameBranchFunction::f(): 4 arguments required");
		const int dec_id = node.id;
		const int anc_id = node.ancestor->id;
		const double branch_length = pow(10.,x[0]);
		const double rho_over_theta = pow(10.,x[1]);
		const double mean_import_length = pow(10.,x[2]);
		const double import_divergence = (excess_divergence_model) ? branch_length + pow(10.,x[3]) : pow(10.,x[3]);
		ML = maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,import_divergence,is_imported);
		return -ML.LOG();
	}
};

class ClonalFrameFunction : public PowellFunction {
public:
	// References to non-member variables
	const marginal_tree &ctree;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	// True member variable
	vector< vector<ImportationState> > is_imported;
	vector< mydouble > ML;
	int neval;
	bool excess_divergence_model;
	clock_t last_update;
	const bool multithread;
	const int root_node;
public:
	ClonalFrameFunction(const marginal_tree &_ctree, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
							  const vector<double> &_pi, bool _excess_divergence_model, const bool _multithread, const int _root_node) : 
	ctree(_ctree), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), pi(_pi), is_imported(_root_node), ML(_root_node), neval(0), 
	excess_divergence_model(_excess_divergence_model), last_update(clock()), multithread(_multithread), root_node(_root_node) {};
	double f(const vector<double>& x) {
		++neval;
		if((clock()-last_update)/CLOCKS_PER_SEC>60.0) {
			cout << "Done " << neval << " iterations" << endl;
			last_update = clock();
		}
		// Process parameters
		const int nparam = root_node+3;
		if(x.size()!=nparam) {
			stringstream errTxt;
			errTxt << "ClonalFrameFunction::f(): " << nparam << " arguments required";
			error(errTxt.str().c_str());
		}
		const double rho_over_theta = pow(10.,x[0]);
		const double mean_import_length = pow(10.,x[1]);
		const double import_divergence_base = pow(10.,x[2]);
		// Calculate likelihood
		int i;
		for(i=0;i<root_node;i++) {
			const double branch_length = pow(10.,x[3+i]);
			const double import_divergence = (excess_divergence_model) ? import_divergence_base + branch_length : import_divergence_base;
			const mt_node &node = ctree.node[i];
			const int dec_id = node.id;
			const int anc_id = node.ancestor->id;
			ML[i] = maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,import_divergence,is_imported[i]);
		}
		// Separate: not included in the parallelizable code
		mydouble treeML = 1.0;
		for(i=0;i<root_node;i++) {
			treeML *= ML[i];
		}
		// The penultimate node (i==ctree.size-2) is a special case: do not try to optimize its branch length nor calculate its likelihood if it is a root node (will be 1)
		// Return minus log-likelihood
		return -treeML.LOG();
	}
};

class ClonalFrameJointBranchLengthFunction : public PowellFunction {
public:
	// References to non-member variables
	const mt_node &node;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	vector<ImportationState> &is_imported;
	// True member variable
	double ML;
	int neval;
	bool excess_divergence_model;
	double rho_over_theta, mean_import_length, import_divergence;
	bool use_viterbi;
public:
	ClonalFrameJointBranchLengthFunction(const bool _use_viterbi, const mt_node &_node, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
									const vector<double> &_pi, bool _excess_divergence_model, const double _rho_over_theta, const double _mean_import_length, const double _import_divergence, vector< ImportationState > &_is_imported) : 
	use_viterbi(_use_viterbi), node(_node), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), pi(_pi), neval(0), excess_divergence_model(_excess_divergence_model),
	rho_over_theta(_rho_over_theta), mean_import_length(_mean_import_length), import_divergence(_import_divergence), is_imported(_is_imported) {};
	double f(const vector<double>& x) {
		++neval;
		if(x.size()!=1) error("ClonalFrameJointBranchLengthFunction::f(): 1 argument required");
		const int dec_id = node.id;
		const int anc_id = node.ancestor->id;
		const double branch_length = pow(10.,x[0]);
		// NB:- in this joint model, excess divergence is additive
		const double final_import_divergence = (excess_divergence_model) ? branch_length + import_divergence : import_divergence;
		if(use_viterbi) {
			ML = maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence,is_imported).LOG();
		} else {
			ML = marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence);
		}
		return -ML;
	}
};

class ClonalFrameFunctionJoint : public PowellFunction {
public:
	// References to non-member variables
	const marginal_tree &ctree;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	vector< vector<ImportationState> > &is_imported;
	// True member variable
	vector< double > ML;
	vector< double > branch_length;
	int neval;
	bool excess_divergence_model;
	clock_t last_update;
	vector<double> &substitutions_per_branch;
	double min_branch_length;
	bool use_viterbi;
	bool coutput;
	double brent_tolerance;
	double powell_tolerance;
	const int root_node;
public:
	ClonalFrameFunctionJoint(const bool _use_viterbi, const marginal_tree &_ctree, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
						const vector<double> &_pi, bool _excess_divergence_model, vector< vector<ImportationState> > &_is_imported, vector<double> &_substitutions_per_branch, const double _min_branch_length, const bool _coutput, const double _brent_tolerance, const double _powell_tolerance, const int _root_node) : 
	use_viterbi(_use_viterbi), ctree(_ctree), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), pi(_pi), ML(_root_node), branch_length(_root_node), neval(0), 
	excess_divergence_model(_excess_divergence_model), last_update(clock()), is_imported(_is_imported), substitutions_per_branch(_substitutions_per_branch), min_branch_length(_min_branch_length), coutput(_coutput), brent_tolerance(_brent_tolerance), powell_tolerance(_powell_tolerance), root_node(_root_node) {};
	double f(const vector<double>& x) {
		++neval;
		if(coutput && (clock()-last_update)/CLOCKS_PER_SEC>60.0) {
			cout << "Done " << neval << " iterations" << endl;
			last_update = clock();
		}
		// Process parameters
		if(!(x.size()==3)) error("ClonalFrameFunctionJoint::f(): 3 arguments required");
		const double rho_over_theta = get_rho_over_theta(x);
		const double mean_import_length = get_mean_import_length(x);
		const double import_divergence = get_import_divergence(x);
		// Calculate likelihood
		int i;
		// The penultimate node (i==ctree.size-2) is a special case: do not try to optimize its branch length nor calculate its likelihood if a root node (will be 1)
		for(i=0;i<root_node;i++) {
			const mt_node &node = ctree.node[i];
			ClonalFrameJointBranchLengthFunction cfb(use_viterbi,node,node_nuc,iscompat,ipat,kappa,pi,excess_divergence_model,rho_over_theta,mean_import_length,import_divergence,is_imported[i]);
			// Maximize the branch-specific likelihood with respect to the branch length
			Powell Pow(cfb);
			Pow.coutput = Pow.brent.coutput = coutput;
			Pow.TOL = brent_tolerance;
			vector<double> param(1,log10(get_branch_length(substitutions_per_branch[i],x)));
			param = Pow.minimize(param,powell_tolerance);
			ML[i] = -Pow.function_minimum;
			branch_length[i] = pow(10.,param[0]);
			// Ensure importation status is updated correctly (force use of viterbi algorithm)
			const bool old_use_viterbi = cfb.use_viterbi;
			cfb.use_viterbi = true;
			cfb.f(param);
			cfb.use_viterbi = old_use_viterbi;
		}
		double treeML = 0.0;
		for(i=0;i<root_node;i++) {
			treeML += ML[i];
		}
		// Return minus log-likelihood
		return -treeML;
	}
	double get_rho_over_theta(const vector<double>& x) {
		return pow(10.,x[0]);
	}
	double get_mean_import_length(const vector<double>& x) {
		return pow(10.,x[1]);
	}
	double get_import_divergence(const vector<double>& x) {
		return pow(10.,x[2]);
	}
	double get_branch_length(const double expected_number_of_substitutions, const vector<double>& x) {
		// Start each branch length so that the expected number of substitutions matches the ancestral state reconstruction
		// SUBJECT to a minimum (positive) branch length
		const double rho_over_theta = get_rho_over_theta(x);
		const double mean_import_length = get_mean_import_length(x);
		const double import_divergence = get_import_divergence(x);
		double b;
		if(excess_divergence_model) {
			// Additive model of excess divergence:
			// Pr(import) = (rho*tau*T/2)/(1+rho*tau*T/2) = (branch_length*rho/theta*tau)/(1+branch_length*rho/theta*tau)
			// where tau = mean import length
			// E(# substitutions) = E(# mutations) + pr(import) * import_divergence
			//                    = branch_length * [1 + (rho/theta*tau*import_divergence)/(1+branch_length*rho/theta*tau)]
			// so   branch_length =approx= E(# substitutions)/(1 + rho/theta*tau*import_divergence)
			b = expected_number_of_substitutions/(1.0+rho_over_theta*mean_import_length*import_divergence);
			if(b!=b || b<min_branch_length) b = min_branch_length;
		} else {
			// Standard model of divergence: equal for every branch
			// Pr(import) = (rho*tau*T/2)/(1+rho*tau*T/2) = (branch_length*rho/theta*tau)/(1+branch_length*rho/theta*tau)
			// where tau = mean import length
			// E(# substitutions) = pr(no import) * E(# mutations) + pr(import) * import_divergence
			//                    = branch_length * (1 + rho/theta*tau*import_divergence)/(1+branch_length*rho/theta*tau)
			// so   branch_length = E(# substitutions)/(1 + rho/theta*tau*(import_divergence-E(# substitutions)))
			b = expected_number_of_substitutions/(1.0+rho_over_theta*mean_import_length*(import_divergence-expected_number_of_substitutions));
			if(expected_number_of_substitutions>=import_divergence) b = expected_number_of_substitutions;
			if(b!=b || b<min_branch_length) b = min_branch_length;
		}
		return b;
	}
};

class ClonalFrameApproxBranchLengthFunction : public PowellFunction {
public:
	// References to non-member variables
	const marginal_tree &ctree;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	// True member variable
	vector< vector<ImportationState> > is_imported;
	vector< mydouble > ML;
	int neval;
	bool excess_divergence_model;
	clock_t last_update;
	vector<double> adjusted_pmut;
	const double min_branch_length;
	vector<double> branch_length_hat;
	const int root_node;
public:
	ClonalFrameApproxBranchLengthFunction(const marginal_tree &_ctree, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
						const vector<double> &_pi, bool _excess_divergence_model, const int _root_node, const double _min_branch_length=0.0000001) : 
	ctree(_ctree), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), pi(_pi), is_imported(_root_node), ML(_root_node), neval(0), 
	excess_divergence_model(_excess_divergence_model), last_update(clock()), adjusted_pmut(_root_node), min_branch_length(_min_branch_length), root_node(_root_node), branch_length_hat(_root_node) {
		int i,j;
		for(i=0;i<root_node;i++) {
			const int dec = i;
			const int anc = ctree.node[i].ancestor->id;
			vector<double> n(4,0.0);
			double nmut = 0;
			for(j=0;j<ipat.size();j++) {
				Nucleotide from = node_nuc[anc][ipat[j]];
				Nucleotide to = node_nuc[dec][ipat[j]];
				n[(int)from]++;
				if(from!=to) ++nmut;
			}
			// Branch-specific observed proportion of mutant sites adjusted for base composition of parent branch. Generally, the denominator should be close to the number of sites.
			adjusted_pmut[i] = nmut/HKY85_expected_rate(n,kappa,pi);
			if(excess_divergence_model) {
				cout << "WARNING: excess divergence model not available for ClonalFrameApproxBranchLengthFunction, ignoring.";
			}
		}
	};
	double f(const vector<double>& x) {
		++neval;
		if((clock()-last_update)/CLOCKS_PER_SEC>60.0) {
			cout << "Done " << neval << " iterations" << endl;
			last_update = clock();
		}
		// Process parameters
		const int nparam = 3;
		if(x.size()!=nparam) {
			stringstream errTxt;
			errTxt << "ClonalFrameApproxBranchLengthFunction::f(): " << nparam << " arguments required";
			error(errTxt.str().c_str());
		}
		const double rho_over_theta = pow(10.,x[0]);
		const double mean_import_length = pow(10.,x[1]);
		const double import_divergence_base = pow(10.,x[2]);
		// Calculate likelihood
		int i;
		for(i=0;i<root_node;i++) {
			// Fast approximate estimate of branch length
			branch_length_hat[i] = (1/mean_import_length*adjusted_pmut[i])/(1/mean_import_length+rho_over_theta*(import_divergence_base-adjusted_pmut[i]));
			if(branch_length_hat[i]<min_branch_length) branch_length_hat[i]=min_branch_length;
			const mt_node &node = ctree.node[i];
			const int dec_id = node.id;
			const int anc_id = node.ancestor->id;
			ML[i] = maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length_hat[i],rho_over_theta,mean_import_length,import_divergence_base,is_imported[i]);
		}
		// Separate: not included in the parallelizable code
		mydouble treeML = 1.0;
		for(i=0;i<root_node;i++) {
			treeML *= ML[i];
		}
		// The penultimate node (i==ctree.size-2) is a special case: do not try to optimize its branch length nor calculate its likelihood (will be 1)
		// Return minus log-likelihood
		return -treeML.LOG();
	}
};


class ClonalFrameParameterFunction : public PowellFunction {
public:
	// References to non-member variables
	const marginal_tree &ctree;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	// True member variable
	vector< vector<ImportationState> > is_imported;
	vector< mydouble > ML;
	int neval;
	bool excess_divergence_model;
	clock_t last_update;
	bool multithread;
	const int root_node;
public:
	ClonalFrameParameterFunction(const marginal_tree &_ctree, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
								 const vector<double> &_pi, bool _excess_divergence_model, const bool _multithread, const int _root_node) : 
	ctree(_ctree), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), pi(_pi), is_imported(_root_node), ML(_root_node), neval(0), 
	excess_divergence_model(_excess_divergence_model), last_update(clock()), multithread(_multithread), root_node(_root_node) {};
	double f(const vector<double>& x) {
		++neval;
		if((clock()-last_update)/CLOCKS_PER_SEC>60.0) {
			cout << "Done " << neval << " iterations" << endl;
			last_update = clock();
		}
		// Process parameters
		const int nparam = 3;
		if(x.size()!=nparam) {
			stringstream errTxt;
			errTxt << "ClonalFrameParameterFunction::f(): " << nparam << " arguments required";
			error(errTxt.str().c_str());
		}
		const double rho_over_theta = pow(10.,x[0]);
		const double mean_import_length = pow(10.,x[1]);
		const double import_divergence_base = pow(10.,x[2]);
		// Calculate likelihood
		int i;
		for(i=0;i<root_node;i++) {
			const mt_node &node = ctree.node[i];
			const double branch_length = node.edge_time;
			const double import_divergence = (excess_divergence_model) ? import_divergence_base + branch_length : import_divergence_base;
			const int dec_id = node.id;
			const int anc_id = node.ancestor->id;
			ML[i] = maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,import_divergence,is_imported[i]);
		}
		// Separate: not included in the parallelizable code
		mydouble treeML = 1.0;
		for(i=0;i<root_node;i++) {
			treeML *= ML[i];
		}
		// The penultimate node (i==ctree.size-2) is a special case: do not try to optimize its branch length nor calculate its likelihood if a root node (will be 1)
		// Return minus log-likelihood
		return -treeML.LOG();
	}
};

class ClonalFrameBranchLengthFunction : public PowellFunction {
public:
	// References to non-member variables
	const mt_node &node;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	// True member variable
	vector<ImportationState> is_imported;
	mydouble ML;
	int neval;
	bool excess_divergence_model;
	double rho_over_theta, mean_import_length, import_divergence;
	const bool multithread;
public:
	ClonalFrameBranchLengthFunction(const mt_node &_node, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
							  const vector<double> &_pi, bool _excess_divergence_model, const double _rho_over_theta, const double _mean_import_length, const double _import_divergence, const bool _multithread) : 
	node(_node), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), pi(_pi), neval(0), excess_divergence_model(_excess_divergence_model),
	rho_over_theta(_rho_over_theta), mean_import_length(_mean_import_length), import_divergence(_import_divergence), multithread(_multithread) {};
	double f(const vector<double>& x) {
		++neval;
		if(x.size()!=1) error("ClonalFrameBranchLengthFunction::f(): 1 argument required");
		const int dec_id = node.id;
		const int anc_id = node.ancestor->id;
		const double branch_length = pow(10.,x[0]);
		const double final_import_divergence = (excess_divergence_model) ? branch_length + import_divergence : import_divergence;
		ML = maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence,is_imported);
		return -ML.LOG();
	}
};

class ClonalFrameRhoPerBranchFunction : public PowellFunction {
public:
	// References to non-member variables
	const mt_node &node;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	vector<ImportationState> &is_imported;
	// True member variable
	double ML;
	int neval;
	bool excess_divergence_model;
	const bool multithread;
	double crude_branch_length;
	double min_branch_length;
	vector<double> which_compat;
public:
	ClonalFrameRhoPerBranchFunction(const mt_node &_node, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
									const vector<double> &_pi, bool _excess_divergence_model, const bool _multithread, vector<ImportationState> &_is_imported, const double _crude_branch_length, const double _min_branch_length) : 
	node(_node), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), pi(_pi), neval(0),
	excess_divergence_model(_excess_divergence_model), multithread(_multithread), is_imported(_is_imported), crude_branch_length(_crude_branch_length), min_branch_length(_min_branch_length) {
		if(!excess_divergence_model) error("ClonalFrameRhoPerBranchFunction::f(): excess_divergence_model is mandatory");
		// Precompute which sites are compatible
		which_compat = vector<double>(0);
		int i;
		for(i=0;i<iscompat.size();i++) {
			if(iscompat[i]) {
				which_compat.push_back((double)i);
			}
		}		
	};
	double f(const vector<double>& x) {
		++neval;
		// Process parameters
		if(!(x.size()==3 || x.size()==4)) error("ClonalFrameRhoPerBranchFunction::f(): 3 or 4 arguments required");
		const double rho_over_theta = pow(10.,x[0]);
		const double import_ratio = 1.0/(1.0+pow(10.,-x[1]));
		const double import_divergence = pow(10.,x[2]);
		double branch_length;
		if(x.size()==3) {
			// Constrain so that the expected number of substitutions equals crude_branch_length
			// crude_branch_length = branch_length + import_ratio/(1+import_ratio)*branch_length*(2+import_divergence)
			//                     = branch_length*(1 + import_ratio/(1+import_ratio)*(2+import_divergence))
			// so    branch_length = crude_branch_length/(1 + import_ratio/(1+import_ratio)*(2+import_divergence))
			branch_length = crude_branch_length/(1.0+import_ratio/(1.0+import_ratio)*(2.0+import_divergence));
		} else {
			branch_length = pow(10.,x[3]);
		}
		if(branch_length<min_branch_length) branch_length = min_branch_length;
		const int dec_id = node.id;
		const int anc_id = node.ancestor->id;
		const double mean_import_length = import_ratio/branch_length/rho_over_theta;
		const double final_import_divergence = (excess_divergence_model) ? branch_length*(2.0 + import_divergence) : import_divergence;
		// Calculate likelihood
		//ML = maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence,is_imported);
		//return -ML.LOG();
		//ML = marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence);
		ML = mydouble_marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,which_compat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence).LOG();
		//cout << "l(" << branch_length << "," <<  rho_over_theta << "," <<  mean_import_length << "," <<  final_import_divergence << ") = " <<  ML << endl;
		return -ML;
	}
};

// Based on ClonalFrameRhoPerBranchFunction this class implements the posterior density function for a single branch under the ClonalFrame model
// A parameterization is implemented that helps impose the original ClonalFrame constraint that the parameters are the same over branches
// (except the branch length which is different for every branch).
//			x[0]	=	log10{ rho/theta (need to check this is exactly the same as in ClonalFrame) }
//			x[1]	=	log10{ import_length (aka 1/delta) }
//			x[2]	=	log10{ import_divergence (aka nu) }
//			x[3]	=	log10{ branch_length }
// In this class, a "driving prior" is implemented. The idea is that this helps with the maximization, which is started from the mode of the driving
// prior, and that post hoc the effect of the prior can be removed to obtain a Normal (Laplace) approximation to the likelihood. The prior used is a multi-
// variate normal distribution with specified mean and precision matrix.
class ClonalFrameLaplacePerBranchFunction : public PowellFunction, public BFGSFunction {
public:
	// References to non-member variables
	const mt_node &node;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	vector<ImportationState> &is_imported;
	// True member variable
	double ML;
	double PR;
	int neval;
	const bool multithread;
	const vector<double> prior_mean;
	const vector<double> prior_precision;
	int parameterization;
	vector<double> which_compat;
public:
	ClonalFrameLaplacePerBranchFunction(const mt_node &_node, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
									const vector<double> &_pi, const bool _multithread, vector<ImportationState> &_is_imported, 
									const vector<double> &_prior_mean, const vector<double> &_prior_precision) : 
	node(_node), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), 
	pi(_pi), neval(0), multithread(_multithread), is_imported(_is_imported),
	prior_mean(_prior_mean), prior_precision(_prior_precision), parameterization(1) {
		if(prior_mean.size()!=4) error("ClonalFrameLaplacePerBranchFunction: prior mean must have length 4");
		if(prior_precision.size()!=4) error("ClonalFrameLaplacePerBranchFunction: prior precision must have length 4");
		// Precompute which sites are compatible
		which_compat = vector<double>(0);
		int i;
		for(i=0;i<iscompat.size();i++) {
			if(iscompat[i]) {
				which_compat.push_back((double)i);
			}
		}		
	};
	// Return function to be minimized
	double f(const vector<double>& x) {
		++neval;
		// Process parameters
		if(x.size()!=4) error("ClonalFrameLaplacePerBranchFunction::f(): 4 arguments required");
		double rho_over_theta, mean_import_length, final_import_divergence, branch_length;
		if(parameterization==0) {
			rho_over_theta = pow(10.,x[0]);
			mean_import_length = pow(10.,x[1]);
			final_import_divergence = pow(10.,x[2]);
			branch_length = pow(10.,x[3]);
			// The following constraint may be important to avoid inverting the signal of recombinant and non-recombinant sites, but for consistency with the original CF parameterization is not applied
			// final_import_divergence = (excess_divergence_model) ? branch_length*(2.0 + import_divergence) : import_divergence;
		} else if(parameterization==1) {
			// 4/6/14 new parameterization (with hard constraint that the branch length be positive)
			rho_over_theta = pow(10.,x[0]);
			mean_import_length = pow(10.,x[1]);
			final_import_divergence = pow(10.,x[2]);
			branch_length = pow(10.,x[3])/(1.0+rho_over_theta*mean_import_length*(final_import_divergence-pow(10.,x[3])));
			if(branch_length<=0.0) return numeric_limits<double>::max();
			// 4/6/14 new constraint: branch_length*rho_over_theta*mean_import_length <= 1. This ensures that the majority of sites are unimported for every branch
			if(branch_length*rho_over_theta*mean_import_length>1.0) return numeric_limits<double>::max();
		} else if(parameterization==2) {
			// 6/6/14 new parameterization informed by mixture distribution with mixing proportion pUnimported and divergences branch_length and final_import_length
			// 10^x[0] = rho*mean_import_length. pUnimported = 1/(1+rho*mean_import_length)
			// 10^x[1] = mean_import_length
			// 10^x[2] = final_import_divergence
			// 10^x[3] = pUnimported*branch_length + (1-pUnimported)*final_import_divergence
			const double pUnimported = 1.0/(1.0+pow(10.,x[0]));
			mean_import_length = pow(10.,x[1]);
			final_import_divergence = pow(10.,x[2]);
			branch_length = (pow(10.,x[3])-(1.0-pUnimported)*final_import_divergence)/pUnimported;
			if(branch_length<=0.0) return numeric_limits<double>::max();
			rho_over_theta = 1.0/branch_length/mean_import_length/pUnimported;
			if(rho_over_theta<=0.0) return numeric_limits<double>::max();
		} else if(parameterization==3) {
			// See "possible reparameterization.docx"
			const double pUnimported = 1.0/(1.0+pow(10.,x[1]));
			final_import_divergence = pow(10.,x[0]+x[3])/(1.0-pUnimported)/(1.0+pow(10.,x[3]));
			branch_length = pow(10.,x[1]-x[3])*final_import_divergence;
			mean_import_length = pow(10.,x[2]);
			rho_over_theta = pow(10.,x[1])/mean_import_length/branch_length;
			// = pow(10.,x[3]-x[2])/final_import_divergence;
		} else {
			error("ClonalFrameLaplacePerBranchFunction::f(): parameterization code not recognized");
		}
		// Test for NaNs
		if(rho_over_theta!=rho_over_theta || mean_import_length!=mean_import_length || final_import_divergence!=final_import_divergence || branch_length!=branch_length) return numeric_limits<double>::max();
		// Identify the nodes
		const int dec_id = node.id;
		const int anc_id = node.ancestor->id;
		// Calculate likelihood
		//ML = marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence);
		ML = mydouble_marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,which_compat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence).LOG();
		// Calculate prior (note that this is only calculated up to a normalizing constant)
		if(parameterization==3) {
			vector<double> xx(4);
			xx[0] = log10(rho_over_theta);
			xx[1] = log10(mean_import_length);
			xx[2] = log10(final_import_divergence);
			xx[3] = log10(branch_length);
			PR = log_prior(xx);
		} else {
			PR = log_prior(x);
		}
		// Print results to screen
		// cout << "Node " << dec_id << " params " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " gave lik " << ML << " prior " << PR << " total " << ML+PR << endl;
		// The optimization routine is assumed to be a minimization routine, hence minus the posterior density is returned
		const double ret = -(ML+PR);
		// Test for NaNs
		if(ret!=ret) return numeric_limits<double>::max();
		return ret;
	}
	double log_prior(const vector<double>& x) {
		double ret = 0.0;
		int i;
		for(i=0;i<4;i++) {
			ret -= 0.5*pow(prior_mean[i]-x[i],2.0)*prior_precision[i];
		}
		return ret;
	}
	vector<double> convert_parameterization_1_to_3(vector<double> &initial_values, const double min_branch_length) {
		// Input parameters, log10-scaled, are rho_over_theta, mean_import_length, final_import_divergence, branch_length (expected divergence, B)
		vector<double> param(4);
		// Expected divergence
		param[0] = initial_values[3];
		// Mean import length
		param[2] = initial_values[1];
		// r/m
		param[3] = initial_values[0]+initial_values[1]+initial_values[2];
		double M = pow(10.,initial_values[3])/(1.0+pow(10.,param[0]+initial_values[1])*(pow(10.,initial_values[2])-pow(10.,initial_values[3])));
		if(M<min_branch_length) {
			M = min_branch_length;
		}
		param[1] = initial_values[0]+initial_values[1]+log10(M);
		return param;
	}
	vector<double> convert_parameterization_0_to_3(vector<double> &initial_values, const double min_branch_length) {
		// Input parameters, log10-scaled, are rho_over_theta, mean_import_length, final_import_divergence, branch_length (divergence at unimported sites, M)
		vector<double> param(4);
		const double nu = pow(10., initial_values[2]);
		const double M = pow(10., initial_values[3]);
		const double Rdelta = pow(10, initial_values[0]+initial_values[1]+initial_values[3]);
		// Expected divergence
		param[0] = log10(nu + (M-nu)/(1.0+Rdelta));
		// Ratio of imported to unimported sites
		param[1] = initial_values[0]+initial_values[1]+initial_values[3];
		// Mean import length
		param[2] = initial_values[1];
		// r/m
		param[3] = initial_values[0]+initial_values[1]+initial_values[2];
		return param;
	}
	vector<double> convert_parameterization_3_to_0(vector<double> &x) {
		const double pUnimported = 1.0/(1.0+pow(10.,x[1]));
		const double final_import_divergence = pow(10.,x[0]+x[3])/(1.0-pUnimported)/(1.0+pow(10.,x[3]));
		const double branch_length = pow(10.,x[1]-x[3])*final_import_divergence;
		const double mean_import_length = pow(10.,x[2]);
		const double rho_over_theta = pow(10.,x[1])/mean_import_length/branch_length;
		vector<double> param(4);
		param[0] = log10(rho_over_theta);
		param[1] = log10(mean_import_length);
		param[2] = log10(final_import_divergence);
		param[3] = log10(branch_length);
		return param;
	}
};

// Based on ClonalFrameLaplacePerBranchFunction this class implements the posterior density function for a single branch under the ClonalFrame model
// A parameterization is implemented that helps impose the original ClonalFrame constraint that the parameters are the same over branches
// (except the branch length which is different for every branch).
//			x[0]	=	log10{ rho/theta (need to check this is exactly the same as in ClonalFrame) }
//			x[1]	=	log10{ import_length (aka 1/delta) }
//			x[2]	=	log10{ import_divergence (aka nu) }
//			x[3]	=	log10{ branch_length }
// In this class, a "driving prior" is implemented. The idea is that this helps with the maximization, which is started from the mode of the driving
// prior, and that post hoc the effect of the prior can be removed to obtain a Normal (Laplace) approximation to the likelihood. The prior used is a multi-
// variate normal distribution with specified mean and precision matrix.
class ClonalFrameMCMCJointFunction {
public:
	// References to non-member variables
	const mt_node *node;
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
	const bool multithread;
	const double min_branch_length;
	const vector<double> prior_mean;
	const vector<double> prior_precision;
	int parameterization;
	vector<double> which_compat;
	bool improper_prior;
	int root_node;
public:
	ClonalFrameMCMCJointFunction(const mt_node *_node, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
										const vector<double> &_pi, const bool _multithread, vector< vector<ImportationState> > &_is_imported, const double _min_branch_length,
										const vector<double> &_prior_mean, const vector<double> &_prior_precision, const bool _improper_prior, const int _root_node) : 
	node(_node), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), 
	pi(_pi), neval(0), multithread(_multithread), is_imported(_is_imported), min_branch_length(_min_branch_length),
	prior_mean(_prior_mean), prior_precision(_prior_precision), improper_prior(_improper_prior), 
	root_node(_root_node), parameterization(0) {
		if(prior_mean.size()!=4) error("ClonalFrameMCMCJointFunction: prior mean must have length 4");
		if(prior_precision.size()!=4) error("ClonalFrameMCMCJointFunction: prior precision must have length 4");
		// Precompute which sites are compatible
		which_compat = vector<double>(0);
		int i;
		for(i=0;i<iscompat.size();i++) {
			if(iscompat[i]) {
				which_compat.push_back((double)i);
			}
		}
	};
	// Return function to be minimized
	double log_posterior(const vector<double>& x) {
		++neval;
		// Process parameters
		double rho_over_theta, mean_import_length, final_import_divergence, branch_length;
		if(parameterization==0) {
			if(x.size()!=3+root_node) error("ClonalFrameMCMCJointFunction::f(): wrong number of arguments for parameterization 0");
			rho_over_theta = pow(10.,x[0]);
			mean_import_length = pow(10.,x[1]);
			final_import_divergence = pow(10.,x[2]);
		} else if(parameterization==1) {
			if(x.size()!=3+root_node) error("ClonalFrameMCMCJointFunction::f(): wrong number of arguments for parameterization 1");
			// Automatic branch lengths
			rho_over_theta = pow(10.,x[0]);
			mean_import_length = pow(10.,x[1]);
			final_import_divergence = pow(10.,x[2]);
			// branch_length = pow(10.,x[3])/(1.0+rho_over_theta*mean_import_length*(final_import_divergence-pow(10.,x[3])));
			// if(branch_length<=min_branch_length) branch_length = min_branch_length;
		} else {
			error("ClonalFrameMCMCJointFunction::f(): parameterization code not recognized");
		}
				
		// Calculate likelihood
		int i;
		ML = 0.0;
		for(i=0;i<root_node;i++) {
			// Identify the nodes
			const int dec_id = node[i].id;
			const int anc_id = node[i].ancestor->id;
			// Branch length
			branch_length = (parameterization==0) ? pow(10.,x[3+i]) : MAX(min_branch_length,pow(10.,x[3+i])/(1.0+rho_over_theta*mean_import_length*(final_import_divergence-pow(10.,x[3+i]))));
			// Update the likelihood
			ML += mydouble_marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,which_compat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence).LOG();
		}
		// Calculate prior (note that this is only calculated up to a normalizing constant)
		PR = log_prior(x);
		const double ret = (ML+PR);
		// Test for NaNs
		if(ret!=ret) return numeric_limits<double>::min();
		return ret;
	}
	// Return function to be minimized
	double log_posterior(const vector<double>& x, vector<double> &partial_post) {
		++neval;
		// Process parameters
		double rho_over_theta, mean_import_length, final_import_divergence, branch_length;
		if(parameterization==0) {
			if(x.size()!=3+root_node) error("ClonalFrameMCMCJointFunction::f(): wrong number of arguments for parameterization 0");
			rho_over_theta = pow(10.,x[0]);
			mean_import_length = pow(10.,x[1]);
			final_import_divergence = pow(10.,x[2]);
		} else if(parameterization==1) {
			if(x.size()!=3+root_node) error("ClonalFrameMCMCJointFunction::f(): wrong number of arguments for parameterization 1");
			// Automatic branch lengths
			rho_over_theta = pow(10.,x[0]);
			mean_import_length = pow(10.,x[1]);
			final_import_divergence = pow(10.,x[2]);
			// branch_length = pow(10.,x[3])/(1.0+rho_over_theta*mean_import_length*(final_import_divergence-pow(10.,x[3])));
			// if(branch_length<=min_branch_length) branch_length = min_branch_length;
		} else {
			error("ClonalFrameMCMCJointFunction::f(): parameterization code not recognized");
		}
		
		// Calculate likelihood
		int i;
		ML = 0.0;
		for(i=0;i<root_node;i++) {
			// Identify the nodes
			const int dec_id = node[i].id;
			const int anc_id = node[i].ancestor->id;
			// Branch length
			branch_length = (parameterization==0) ? pow(10.,x[3+i]) : MAX(min_branch_length,pow(10.,x[3+i])/(1.0+rho_over_theta*mean_import_length*(final_import_divergence-pow(10.,x[3+i]))));
			// Update the likelihood
			partial_post[i] = mydouble_marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,which_compat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence).LOG();
			ML += partial_post[i];
		}
		// Calculate prior (note that this is only calculated up to a normalizing constant)
		PR = log_prior(x);
		const double ret = (ML+PR);
		// Test for NaNs
		if(ret!=ret) return numeric_limits<double>::min();
		return ret;
	}
	// Return function to be minimized
	double log_posterior(const vector<double>& x, vector<double> &partial_post, const int branch_update) {
		++neval;
		// Process parameters
		double rho_over_theta, mean_import_length, final_import_divergence, branch_length;
		if(parameterization==0) {
			if(x.size()!=3+root_node) error("ClonalFrameMCMCJointFunction::f(): wrong number of arguments for parameterization 0");
			rho_over_theta = pow(10.,x[0]);
			mean_import_length = pow(10.,x[1]);
			final_import_divergence = pow(10.,x[2]);
		} else if(parameterization==1) {
			if(x.size()!=3+root_node) error("ClonalFrameMCMCJointFunction::f(): wrong number of arguments for parameterization 1");
			// Automatic branch lengths
			rho_over_theta = pow(10.,x[0]);
			mean_import_length = pow(10.,x[1]);
			final_import_divergence = pow(10.,x[2]);
			// branch_length = pow(10.,x[3])/(1.0+rho_over_theta*mean_import_length*(final_import_divergence-pow(10.,x[3])));
			// if(branch_length<=min_branch_length) branch_length = min_branch_length;
		} else {
			error("ClonalFrameMCMCJointFunction::f(): parameterization code not recognized");
		}
		
		// Calculate likelihood
		int i;
		ML = 0.0;
		for(i=0;i<root_node;i++) {
			if(branch_update==i) {
				// Identify the nodes
				const int dec_id = node[i].id;
				const int anc_id = node[i].ancestor->id;
				// Branch length
				branch_length = (parameterization==0) ? pow(10.,x[3+i]) : MAX(min_branch_length,pow(10.,x[3+i])/(1.0+rho_over_theta*mean_import_length*(final_import_divergence-pow(10.,x[3+i]))));
				// Update the likelihood
				partial_post[i] = mydouble_marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,which_compat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,final_import_divergence).LOG();
				ML += partial_post[i];
			} else {
				ML += partial_post[i];
			}
		}
		// Calculate prior (note that this is only calculated up to a normalizing constant)
		PR = log_prior(x);
		const double ret = (ML+PR);
		// Test for NaNs
		if(ret!=ret) return numeric_limits<double>::min();
		return ret;
	}
	double log_prior(const vector<double>& x) {
		if(improper_prior) return 0.0;
		error("ClonalFrameMCMCJointFunction::log_prior(): normal prior not yet implemented");
		return 0.0;
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

class ClonalFrameSingleRho : public PowellFunction {
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
	const bool use_viterbi;
	double ML;
	int neval;
	bool excess_divergence_model;
	const bool multithread;
	vector<double> &substitutions_per_branch;
	double min_branch_length;
	const int root_node;
public:
	ClonalFrameSingleRho(const bool _use_viterbi, const marginal_tree &_tree, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
									const vector<double> &_pi, bool _excess_divergence_model, const bool _multithread, vector< vector<ImportationState> > &_is_imported, vector<double> &_substitutions_per_branch, const double _min_branch_length, const int _root_node) : 
	use_viterbi(_use_viterbi), tree(_tree), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), pi(_pi), neval(0),
	excess_divergence_model(_excess_divergence_model), multithread(_multithread), is_imported(_is_imported), substitutions_per_branch(_substitutions_per_branch), min_branch_length(_min_branch_length), root_node(_root_node) {
		if(excess_divergence_model) error("ClonalFrameSingleRho::f(): excess_divergence_model not implemented");	
	};
	double f(const vector<double>& x) {
		++neval;
		// Process parameters
		if(!(x.size()==3)) error("ClonalFrameSingleRho::f(): 3 arguments required");
		const double rho_over_theta = get_rho_over_theta(x);
		const double mean_import_length = get_mean_import_length(x);
		const double import_divergence = get_import_divergence(x);
		ML = 0.0;
		int i;
		for(i=0;i<root_node;i++) {
			const mt_node &node = tree.node[i];
			const double branch_length = get_branch_length(substitutions_per_branch[i],x);
			const int dec_id = node.id;
			const int anc_id = node.ancestor->id;
			if(use_viterbi) {
				// Calculate likelihood: Viterbi
				ML += maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,import_divergence,is_imported[i]).LOG();
			} else {
				// Calculate likelihood: Forward algorithm
				ML += marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,import_divergence);
			}
		}
		return -ML;
	}
	double get_rho_over_theta(const vector<double>& x) {
		return pow(10.,x[0]);
	}
	double get_mean_import_length(const vector<double>& x) {
		return pow(10.,x[1]);
	}
	double get_import_divergence(const vector<double>& x) {
		return pow(10.,x[2]);
	}
	double get_branch_length(const double expected_number_of_substitutions, const vector<double>& x) {
		if(excess_divergence_model) error("ClonalFrameSingleRho::f(): excess_divergence_model not implemented");
		// Constrain EACH BRANCH so that the expected number of substitutions matches the ancestral state reconstruction
		// SUBJECT to a minimum (positive) branch length
		// Standard model of divergence: equal for every branch
		// Pr(import) = (rho*tau*T/2)/(1+rho*tau*T/2) = (branch_length*rho/theta*tau)/(1+branch_length*rho/theta*tau)
		// where tau = mean import length
		// E(# substitutions) = pr(no import) * E(# mutations) + pr(import) * import_divergence
		//                    = branch_length * (1 + rho/theta*tau*import_divergence)/(1+branch_length*rho/theta*tau)
		// so   branch_length = E(# substitutions)/(1 + rho/theta*tau*(import_divergence-E(# substitutions)))
		const double rho_over_theta = get_rho_over_theta(x);
		const double mean_import_length = get_mean_import_length(x);
		const double import_divergence = get_import_divergence(x);
		double b = expected_number_of_substitutions/(1.0+rho_over_theta*mean_import_length*(import_divergence-expected_number_of_substitutions));
		if(expected_number_of_substitutions>=import_divergence) b = expected_number_of_substitutions;
		if(b!=b || b<min_branch_length) b = min_branch_length;
		return b;
	}
};

class ClonalFrameViterbiTrainingPerBranch {
public:
	// References to non-member variables
	const mt_node &node;
	const Matrix<Nucleotide> &node_nuc;
	const vector<bool> &iscompat;
	const vector<int> &ipat;
	const double kappa;
	const vector<double> &pi;
	vector<ImportationState> &is_imported;
	// True member variable
	double ML;
	double PR;
	int neval;
	const vector<double> prior_mean;
	const vector<double> prior_precision;
	vector<double> which_compat;
public:
	ClonalFrameViterbiTrainingPerBranch(const mt_node &_node, const Matrix<Nucleotide> &_node_nuc, const vector<bool> &_iscompat, const vector<int> &_ipat, const double _kappa,
										const vector<double> &_pi, vector<ImportationState> &_is_imported, 
										const vector<double> &_prior_mean, const vector<double> &_prior_precision) : 
	node(_node), node_nuc(_node_nuc), iscompat(_iscompat), ipat(_ipat), kappa(_kappa), 
	pi(_pi), neval(0), is_imported(_is_imported),
	prior_mean(_prior_mean), prior_precision(_prior_precision) {
		if(prior_mean.size()!=4) error("ClonalFrameLaplacePerBranchFunction: prior mean must have length 4");
		if(prior_precision.size()!=4) error("ClonalFrameLaplacePerBranchFunction: prior precision must have length 4");
		// Precompute which sites are compatible
		which_compat = vector<double>(0);
		int i;
		for(i=0;i<iscompat.size();i++) {
			if(iscompat[i]) {
				which_compat.push_back((double)i);
			}
		}		
	}
	vector<double> maximize_likelihood(const vector<double> &param) {
		const int dec_id = node.id;
		const int anc_id = node.ancestor->id;
		double branch_length = pow(10.,param[3]);
		double rho_over_theta = pow(10.,param[0]);
		double mean_import_length = pow(10.,param[1]);
		double import_divergence = pow(10.,param[2]);
		ML = Viterbi_training(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,import_divergence,is_imported,neval);
		vector<double> ret(4);
		ret[0] = log10(rho_over_theta);
		ret[1] = log10(mean_import_length);
		ret[2] = log10(import_divergence);
		ret[3] = log10(branch_length);
		return(ret);
	}
	double log_likelihood(const vector<double> &param) {
		const int dec_id = node.id;
		const int anc_id = node.ancestor->id;
		double branch_length = pow(10.,param[3]);
		double rho_over_theta = pow(10.,param[0]);
		double mean_import_length = pow(10.,param[1]);
		double import_divergence = pow(10.,param[2]);
		return maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,which_compat,ipat,kappa,pi,branch_length,rho_over_theta,mean_import_length,import_divergence,is_imported).LOG();
	}
};

#endif // _MAIN_H_
