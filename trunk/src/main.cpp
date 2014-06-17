/*  Copyright 2013 Daniel Wilson and Xavier Didelot.
 *
 *  main.cpp
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
#include <main.h>

int main (const int argc, const char* argv[]) {
	clock_t start_time = clock();
	// Output version number
	cout << "ClonalFrameML version 0." << ClonalFrameML_SVNRevision << endl;
	// Process the command line arguments	
	if(argc<5) {
		stringstream errTxt;
		errTxt << "Syntax: ClonalFrameML newick_file fasta_file kappa output_file [OPTIONS]" << endl;
		errTxt << endl;
		errTxt << "The following options are available:" << endl;
		errTxt << "-fasta_file_list               true or false (default)   Take fasta_file to be a white-space separated file list." << endl;
		errTxt << "-correct_branch_lengths        true (default) or false   Correct branch lengths using ClonalFrame model." << endl;
		errTxt << "-excess_divergence_model       true or false (default)   Use the 'excess divergence' model. Mandatory for two sequences." << endl;
		errTxt << "-ignore_incomplete_sites       true or false (default)   Ignore sites with any ambiguous bases." << endl;
		errTxt << "-ignore_user_sites             sites_file                Ignore sites listed in whitespace-separated sites_file." << endl;
		errTxt << "-reconstruct_invariant_sites   true or false (default)   Reconstruct the ancestral states at invariant sites." << endl;
		errTxt << "-use_incompatible_sites        true or false (default)   Use homoplasious and multiallelic sites to correct branch lengths." << endl;
		errTxt << "-brent_tolerance               tolerance (default .001)  Set the tolerance of the Brent routine." << endl;
		errTxt << "-powell_tolerance              tolerance (default .001)  Set the tolerance of the Powell routine." << endl;
		errTxt << "-joint_branch_param            true or false (default)   Jointly optimize branch lengths and recombination parameters." << endl;
		errTxt << "-rho_per_branch                true or false (default)   Estimate recombination parameters separately for each branch." << endl;
		errTxt << "-rho_per_branch_no_lrt         true or false (default)   As above but suppress likelihood ratio test for recombination." << endl;
		errTxt << "-single_rho_viterbi            true or false (default)   Jointly optimize recombination parameters using Viterbi algorithm." << endl;
		errTxt << "-single_rho_forward            true or false (default)   Jointly optimize recombination parameters using forward algorithm." << endl;
		errTxt << "-rescale_no_recombination      true or false (default)   Rescale branch lengths for given sites with no recombination model." << endl;		
		errTxt << "-multithread                   true or false (default)   Enable OpenMP parallel code. Overhead may cancel out gains." << endl;
		errTxt << "-show_progress                 true or false (default)   Output the progress of the maximum likelihood routines." << endl;
		errTxt << "-compress_reconstructed_sites  true (default) or false   Reduce the number of columns in the output FASTA file." << endl;
		errTxt << "-initial_rho_over_theta        value > 0 (default 0.1)   Initial value of rho/theta used in the search." << endl;
		errTxt << "-initial_import_divergence     value > 0 (default 0.1)   Initial value of import divergence used in the search." << endl;
		errTxt << "-initial_mean_import_length    value > 1 (default 500)   Initial value of mean import length used in the search." << endl;
		errTxt << "-min_branch_length             value > 0 (default 1e-7)  Minimum branch length." << endl;
		errTxt << "-mcmc_per_branch               true or false (default)   Estimate by MCMC recombination parameters for each branch." << endl;
		errTxt << "-laplace_approx                true or false (default)   rho_per_branch model with approximation of the joint posterior." << endl;
		errTxt << "-driving_prior_mean            4 values (df \"0 0 0 0\")   Mean of the driving prior used by the Laplace approximation." << endl;
		errTxt << "-driving_prior_precision       4 values (df \"1 1 1 1\")   Precision of the driving prior used by the Laplace approximation." << endl;
		errTxt << "-grid_approx                   value 0/2+  (default 0)   Number of points for grid approximation (0 = off)." << endl;
		errTxt << "-mcmc                          true or false (default)   Estimate by MCMC recombination parameters for all branches." << endl;
		errTxt << "-mcmc_infer_branch_lengths     true or false (default)   Estimate by MCMC branch lengths for all branches." << endl;
		errTxt << "-initial_values                3 values/empty (def \"\")   Initial values used by the Laplace approximation only." << endl;
		error(errTxt.str().c_str());
	}
	// Process required arguments
	const char* newick_file = argv[1];
	const char* fasta_file = argv[2];
	const double kappa = atof(argv[3]);
	const char* out_file = argv[4];
	string tree_out_file = string(out_file) + ".labelled_tree.newick";
	string fasta_out_file = string(out_file) + ".ML_sequence.fasta";
	string xref_out_file = string(out_file) + ".position_cross_reference.txt";
	string import_out_file = string(out_file) + ".importation_status.txt";
	string laplace_out_file = string(out_file) + ".laplace.txt";
	string mcmc0_out_file = string(out_file) + ".mcmc0.txt";
	string mcmc1_out_file = string(out_file) + ".mcmc1.txt";
	string mcmc2_out_file = string(out_file) + ".mcmc2.txt";
	string mcmc3_out_file = string(out_file) + ".mcmc3.txt";
	string prop0_out_file = string(out_file) + ".prop0.txt";
	string prop1_out_file = string(out_file) + ".prop1.txt";
	string prop2_out_file = string(out_file) + ".prop2.txt";
	string prop3_out_file = string(out_file) + ".prop3.txt";
	string loglik_out_file = string(out_file) + ".loglik.txt";
	string loglik_prop_out_file = string(out_file) + ".loglik_prop.txt";
	// Set default options
	ArgumentWizard arg;
	arg.case_sensitive = false;
	string fasta_file_list="false", correct_branch_lengths="true", excess_divergence_model="false", ignore_incomplete_sites="false", ignore_user_sites="", reconstruct_invariant_sites="false";
	string use_incompatible_sites="false", joint_branch_param="false", rho_per_branch="false", rho_per_branch_no_LRT="false", rescale_no_recombination="false";
	string single_rho_viterbi="false", single_rho_forward="false", multithread="false", show_progress="false", compress_reconstructed_sites="true", laplace_approx="false", mcmc_per_branch = "false";
	string string_driving_prior_mean="0 0 0 0", string_driving_prior_precision="1 1 1 1", mcmc_joint = "false", mcmc_infer_branch_lengths = "false", string_initial_values = "";
	double brent_tolerance = 1.0e-3, powell_tolerance = 1.0e-3, initial_rho_over_theta = 0.1, initial_mean_import_length = 500.0, initial_import_divergence = 0.1, global_min_branch_length = 1.0e-7;
	double grid_approx = 0.0;
	// Process options
	arg.add_item("fasta_file_list",				TP_STRING, &fasta_file_list);
	arg.add_item("correct_branch_lengths",		TP_STRING, &correct_branch_lengths);
	arg.add_item("excess_divergence_model",		TP_STRING, &excess_divergence_model);
	arg.add_item("ignore_incomplete_sites",		TP_STRING, &ignore_incomplete_sites);
	arg.add_item("ignore_user_sites",			TP_STRING, &ignore_user_sites);
	arg.add_item("reconstruct_invariant_sites", TP_STRING, &reconstruct_invariant_sites);
	arg.add_item("use_incompatible_sites",		TP_STRING, &use_incompatible_sites);
	arg.add_item("brent_tolerance",				TP_DOUBLE, &brent_tolerance);
	arg.add_item("powell_tolerance",			TP_DOUBLE, &powell_tolerance);
	arg.add_item("joint_branch_param",			TP_STRING, &joint_branch_param);
	arg.add_item("rho_per_branch",				TP_STRING, &rho_per_branch);
	arg.add_item("rho_per_branch_no_lrt",		TP_STRING, &rho_per_branch_no_LRT);
	arg.add_item("rescale_no_recombination",	TP_STRING, &rescale_no_recombination);
	arg.add_item("single_rho_viterbi",			TP_STRING, &single_rho_viterbi);
	arg.add_item("single_rho_forward",			TP_STRING, &single_rho_forward);
	arg.add_item("multithread",			        TP_STRING, &multithread);
	arg.add_item("show_progress",				TP_STRING, &show_progress);
	arg.add_item("compress_reconstructed_sites",TP_STRING, &compress_reconstructed_sites);	
	arg.add_item("initial_rho_over_theta",		TP_DOUBLE, &initial_rho_over_theta);	
	arg.add_item("initial_mean_import_length",	TP_DOUBLE, &initial_mean_import_length);	
	arg.add_item("initial_import_divergence",	TP_DOUBLE, &initial_import_divergence);	
	arg.add_item("min_branch_length",			TP_DOUBLE, &global_min_branch_length);
	arg.add_item("mcmc_per_branch",				TP_STRING, &mcmc_per_branch);
	arg.add_item("laplace_approx",				TP_STRING, &laplace_approx);
	arg.add_item("driving_prior_mean",			TP_STRING, &string_driving_prior_mean);
	arg.add_item("driving_prior_precision",		TP_STRING, &string_driving_prior_precision);
	arg.add_item("grid_approx",					TP_DOUBLE, &grid_approx);
	arg.add_item("mcmc",						TP_STRING, &mcmc_joint);
	arg.add_item("mcmc_infer_branch_lengths",	TP_STRING, &mcmc_infer_branch_lengths);
	arg.add_item("initial_values",				TP_STRING, &string_initial_values);
	arg.read_input(argc-4,argv+4);
	bool FASTA_FILE_LIST				= string_to_bool(fasta_file_list,				"fasta_file_list");
	bool CORRECT_BRANCH_LENGTHS			= string_to_bool(correct_branch_lengths,		"correct_branch_lengths");
	bool EXCESS_DIVERGENCE_MODEL		= string_to_bool(excess_divergence_model,		"excess_divergence_model");
	bool IGNORE_INCOMPLETE_SITES		= string_to_bool(ignore_incomplete_sites,		"ignore_incomplete_sites");
	bool RECONSTRUCT_INVARIANT_SITES	= string_to_bool(reconstruct_invariant_sites,	"reconstruct_invariant_sites");
	bool USE_INCOMPATIBLE_SITES			= string_to_bool(use_incompatible_sites,		"use_incompatible_sites");
	bool OLD_JOINT_BRANCH_PARAM			= false;
	bool JOINT_BRANCH_PARAM				= string_to_bool(joint_branch_param,			"joint_branch_param");
	bool RHO_PER_BRANCH					= string_to_bool(rho_per_branch,				"rho_per_branch");
	bool RHO_PER_BRANCH_NO_LRT			= string_to_bool(rho_per_branch_no_LRT,			"rho_per_branch_no_lrt");
	bool RESCALE_NO_RECOMBINATION		= string_to_bool(rescale_no_recombination,		"rescale_no_recombination");
	bool SINGLE_RHO_VITERBI				= string_to_bool(single_rho_viterbi,			"single_rho_viterbi");
	bool SINGLE_RHO_FORWARD				= string_to_bool(single_rho_forward,			"single_rho_forward");
	bool MULTITHREAD					= string_to_bool(multithread,					"multithread");
	bool SHOW_PROGRESS					= string_to_bool(show_progress,					"show_progress");
	bool COMPRESS_RECONSTRUCTED_SITES	= string_to_bool(compress_reconstructed_sites,	"compress_reconstructed_sites");
	bool MCMC_PER_BRANCH				= string_to_bool(mcmc_per_branch,				"mcmc_per_branch");
	bool LAPLACE_APPROX					= string_to_bool(laplace_approx,				"laplace_approx");
	bool GRID_APPROX = (grid_approx != 0.0);
	bool MCMC_JOINT						= string_to_bool(mcmc_joint,					"mcmc");
	bool MCMC_INFER_BRANCH_LENGTHS		= string_to_bool(mcmc_infer_branch_lengths,		"mcmc_infer_branch_lengths");
	if(brent_tolerance<=0.0 || brent_tolerance>=0.1) {
		stringstream errTxt;
		errTxt << "brent_tolerance value out of range (0,0.1], default 0.001";
		error(errTxt.str().c_str());
	}
	if(powell_tolerance<=0.0 || powell_tolerance>=0.1) {
		stringstream errTxt;
		errTxt << "powell_tolerance value out of range (0,0.1], default 0.001";
		error(errTxt.str().c_str());
	}
	if(((int)JOINT_BRANCH_PARAM + (int)RHO_PER_BRANCH + (int)RHO_PER_BRANCH_NO_LRT + (int)RESCALE_NO_RECOMBINATION) + (int)SINGLE_RHO_VITERBI + (int)SINGLE_RHO_FORWARD + (int)MCMC_PER_BRANCH + (int)LAPLACE_APPROX + (int)GRID_APPROX + (int)MCMC_JOINT>1) {
		stringstream errTxt;
		errTxt << "joint_branch_param, rho_per_branch, rho_per_branch_no_lrt, rescale_no_recombination, single_rho_viterbi, single_rho_forward, mcmc_per_branch, laplace_approx, grid_approx and mcmc are mutually incompatible";
		error(errTxt.str().c_str());
	}
	if((EXCESS_DIVERGENCE_MODEL || JOINT_BRANCH_PARAM || RHO_PER_BRANCH || RESCALE_NO_RECOMBINATION || SINGLE_RHO_VITERBI || SINGLE_RHO_FORWARD || MCMC_PER_BRANCH || LAPLACE_APPROX || GRID_APPROX || MCMC_JOINT) && !CORRECT_BRANCH_LENGTHS) {
		stringstream wrnTxt;
		wrnTxt << "branch correction options will be ignored because correct_branch_lengths=false";
		warning(wrnTxt.str().c_str());
	}
	if(!EXCESS_DIVERGENCE_MODEL && RHO_PER_BRANCH) {
		stringstream wrnTxt;
		wrnTxt << "rho_per_branch implies excess_divergence_model";
		warning(wrnTxt.str().c_str());
		EXCESS_DIVERGENCE_MODEL = true;
	}
	if(!EXCESS_DIVERGENCE_MODEL && MCMC_PER_BRANCH) {
		stringstream wrnTxt;
		wrnTxt << "mcmc_per_branch implies excess_divergence_model";
		warning(wrnTxt.str().c_str());
		EXCESS_DIVERGENCE_MODEL = true;
	}
	if(EXCESS_DIVERGENCE_MODEL && LAPLACE_APPROX) {
		stringstream wrnTxt;
		wrnTxt << "laplace_approx implies no excess_divergence_model";
		warning(wrnTxt.str().c_str());
		EXCESS_DIVERGENCE_MODEL = false;
	}
	if(EXCESS_DIVERGENCE_MODEL && (SINGLE_RHO_FORWARD || SINGLE_RHO_VITERBI)) {
		stringstream errTxt;
		errTxt << "single_rho not implemented under excess_divergence_model";
		error(errTxt.str().c_str());
	}
	if(EXCESS_DIVERGENCE_MODEL && MCMC_JOINT) {
		stringstream wrnTxt;
		wrnTxt << "mcmc currently implies no excess_divergence_model";
		warning(wrnTxt.str().c_str());
		EXCESS_DIVERGENCE_MODEL = false;
	}
	if(MULTITHREAD) {
		cout << "WARNING: multithreaded version not implemented, ignoring." << endl;
	}
	if(initial_mean_import_length<1) {
		stringstream errTxt;
		errTxt << "initial_mean_import_length must be greater than 1";
		error(errTxt.str().c_str());
	}
	if(global_min_branch_length<=0.0) {
		error("Minimum branch length must be positive");
	}
//	if(LAPLACE_APPROX && !(RHO_PER_BRANCH)) {
//		error("Laplace approximation only applied under rho per branch model");
//	}
	// Process the driving prior mean and precision
	vector<double> driving_prior_mean(0), driving_prior_precision(0);
	stringstream sstream_driving_prior_mean;
	sstream_driving_prior_mean << string_driving_prior_mean;
	int i;
	for(i=0;i<1000;i++) {
		if(sstream_driving_prior_mean.eof()) break;
		double driving_prior_mean_elem;
		sstream_driving_prior_mean >> driving_prior_mean_elem;
		if(sstream_driving_prior_mean.fail()) error("Could not interpret value specified by driving_prior_mean");
		driving_prior_mean.push_back(driving_prior_mean_elem);
	}
	if(i==1000) error("Maximum length of vector exceeded by driving_prior_mean");
	stringstream sstream_driving_prior_precision;
	sstream_driving_prior_precision << string_driving_prior_precision;
	for(i=0;i<1000;i++) {
		if(sstream_driving_prior_precision.eof()) break;
		double driving_prior_precision_elem;
		sstream_driving_prior_precision >> driving_prior_precision_elem;
		if(sstream_driving_prior_precision.fail()) error("Could not interpret value specified by driving_prior_precision");
		driving_prior_precision.push_back(driving_prior_precision_elem);
	}
	if(driving_prior_mean.size()!=4) error("driving_prior_mean must have 4 values separated by spaces");
	if(driving_prior_precision.size()!=4) error("driving_prior_precision must have 4 values separated by spaces");
	if(fabs(grid_approx-(double)(int)(grid_approx))>1.0e-12) error("grid_approx must have an integer value");
	if(GRID_APPROX && grid_approx<2.0) error("grid_approx must be 0 (off) or 2 or more");
	if(!MCMC_JOINT & MCMC_INFER_BRANCH_LENGTHS) warning("mcmc_infer_branch_lengths will be ignored because -mcmc is false");
	// Process the initial values for the Laplace approximation
	vector<double> initial_values(0);
	stringstream sstream_initial_values;
	sstream_initial_values << string_initial_values;
	for(i=0;i<1000;i++) {
		if(sstream_initial_values.eof()) break;
		double initial_values_elem;
		sstream_initial_values >> initial_values_elem;
		if(sstream_initial_values.fail()) error("Could not interpret value specified by initial_values");
		initial_values.push_back(initial_values_elem);
	}
	if(i==1000) error("Maximum length of vector exceeded by initial_values");
	if(!(initial_values.size()==0 || initial_values.size()==3)) error("initial values must have 0 or 3 values separated by spaces");
	if(initial_values.size()>0 && !LAPLACE_APPROX) warning("-initial_values only used by -laplace_approx currently");

	
	// Open the FASTA file(s)
	DNA fa;
	if(FASTA_FILE_LIST) {
		ifstream file_list(fasta_file);
		if(!file_list.is_open()) {
			stringstream errTxt;
			errTxt << "could not find file " << fasta_file;
			error(errTxt.str().c_str());
		}
		int n = 0;
		int L = -1;
		while(!file_list.eof()) {
			string filename;
			file_list >> filename;
			// Pre-check: does it exist?
			ifstream file_list1(filename.c_str());
			if(!file_list1.is_open()) {
				stringstream errTxt;
				errTxt << "could not find listed file " << fasta_file;
				error(errTxt.str().c_str());
			}
			// Read the file
			DNA fa1(filename.c_str());
			n += fa1.nseq;
			if(L==-1) L = fa1.lseq;
			if(fa1.lseq!=L) {
				stringstream errTxt;
				errTxt << "listed file " << fasta_file << " had sequence length " << fa1.lseq << " expecting " << L;
				error(errTxt.str().c_str());
			}
			// Add to list
			int ni;
			for(ni=0;ni<fa1.nseq;ni++) {
				fa.label.push_back(fa1.label[ni]);
				fa.sequence.push_back(fa1.sequence[ni]);
				fa.nseq++;
				fa.ntimes.push_back(fa1.ntimes[ni]);
			}
		}
	} else {
		fa.readFASTA_1pass(fasta_file);
	}
	cout << "Read " << fa.nseq << " sequences of length " << fa.lseq << " sites from " << fasta_file << endl;
	if(fa.nseq==2) {
		if(!EXCESS_DIVERGENCE_MODEL) cout << "WARNING: with only two sequences, the excess divergence model is mandatory." << endl;
		EXCESS_DIVERGENCE_MODEL = true;
	}
	// Open the Newick file and convert to internal rooted tree format, outputting the names of the tips and internal nodes
	NewickTree newick = read_Newick(newick_file);
	vector<string> ctree_node_labels;
	const bool is_rooted = (newick.root.dec.size()==2);
	marginal_tree ctree = (is_rooted) ? convert_rooted_NewickTree_to_marginal_tree(newick,fa.label,ctree_node_labels) : convert_unrooted_NewickTree_to_marginal_tree(newick,fa.label,ctree_node_labels);
	const int root_node = (is_rooted) ? ctree.size-1 : ctree.size-2;
	// Open the list of sites to ignore
	vector<bool> ignore_site(fa.lseq,false);
	if(ignore_user_sites!="") {
		ifstream user_sites(ignore_user_sites.c_str());
		int debug_last_elem = -1;
		while(!user_sites.eof()) {
			int elem;
			user_sites >> elem;
			elem--;
			if(!(elem>=0 && elem<fa.lseq)) {
				stringstream errTxt;
				errTxt << "In file " << ignore_user_sites << " value " << elem+1 << " not in range 1-" << fa.lseq;
				error(errTxt.str().c_str());
			}
			ignore_site[elem] = true;
			debug_last_elem = elem;
		}
	}
	
	// Compute compatibility and test every site for any sequences with 'N','-','X' or '?'
	// Key to results: -1: invariant, 0: compatible biallelic (including singletons), 1: incompatible biallelic, 2: more than two alleles
	vector<bool> anyN;
	vector<int> compat = compute_compatibility(fa,ctree,anyN,false);
	if(IGNORE_INCOMPLETE_SITES) {
		for(i=0;i<fa.lseq;i++) {
			if(anyN[i]) ignore_site[i] = true;
		}
	}
	
	// Flag sites for Imputation and Reconstruction of Ancestral States
	int nIRAS = 0;
	vector<bool> isIRAS(fa.lseq,false);
	for(i=0;i<fa.lseq;i++) {
		if(ignore_site[i]==false && (compat[i]!=-1 || RECONSTRUCT_INVARIANT_SITES==true)) {
			isIRAS[i] = true;
			++nIRAS;
		}
	}
	
	// Flag sites for Branch Length Correction
	int nBLC = 0;
	vector<bool> isBLC(fa.lseq,false);
	for(i=0;i<fa.lseq;i++) {
		if(ignore_site[i]==false && (compat[i]<=0 || USE_INCOMPATIBLE_SITES==true)) {
			isBLC[i] = true;
			++nBLC;
		}
	}
	
	// IMPUTATION AND RECONSTRUCTION OF ANCESTRAL STATES
	if(true) {
		// Convert FASTA file to internal representation of nucleotides for Branch Length Correction
		vector<double> empirical_nucleotide_frequencies(4,0.25);
		Matrix<Nucleotide> nuc = FASTA_to_nucleotide(fa,empirical_nucleotide_frequencies,isIRAS);
		// Identify and count unique patterns
		vector<string> pat;				// Pattern as string of AGCTNs
		vector<int> pat1, cpat, ipat;	// First example of each pattern, number of sites with that pattern, the pattern at each (compatible) site (-1 otherwise)
		vector<bool> nuc_ispoly(nuc.ncols(),true);
		find_alignment_patterns(nuc,nuc_ispoly,pat,pat1,cpat,ipat);
		// Storage for the MLE of the nucleotide sequence at every node
		Matrix<Nucleotide> node_nuc;
		// Begin by computing the joint maximum likelihood ancestral sequences
		mydouble ML = maximum_likelihood_ancestral_sequences(nuc,ctree,kappa,empirical_nucleotide_frequencies,pat1,cpat,node_nuc);

		cout << "IMPUTATION AND RECONSTRUCTION OF ANCESTRAL STATES:" << endl;
		cout << "Analysing " << nIRAS << " sites" << endl;
		// Report the estimated equilibrium frequencies
		cout << "Empirical nucleotide frequencies:   A " << round(1000*empirical_nucleotide_frequencies[Adenine])/10 << "%   C " << round(1000*empirical_nucleotide_frequencies[Cytosine])/10;
		cout << "%   G " << round(1000*empirical_nucleotide_frequencies[Guanine])/10 << "%   T " << round(1000*empirical_nucleotide_frequencies[Thymine])/10 << "%" << endl;
		// Report the ML
		cout << "Maximum log-likelihood for imputation and ancestral state reconstruction = " << ML.LOG() << endl;

		if(!COMPRESS_RECONSTRUCTED_SITES) cout << "WARNING: -compress_reconstructed_sites=false not yet implemented, ignoring." << endl;
		// Output the ML reconstructed sequences
		write_ancestral_fasta(node_nuc, ctree_node_labels, fasta_out_file.c_str());
		// For every position in the original FASTA file, output the corresponding position in the output FASTA file, or -1 (not included)
		write_position_cross_reference(isIRAS, ipat, xref_out_file.c_str());
		cout << "Wrote imputed and reconstructed ancestral states to " << fasta_out_file << endl;
		cout << "Wrote position cross-reference file to " << xref_out_file << endl;
	}
	
	// BRANCH LENGTH CORRECTION
	if(CORRECT_BRANCH_LENGTHS) {
		// Convert FASTA file to internal representation of nucleotides for Branch Length Correction
		vector<double> empirical_nucleotide_frequencies(4,0.25);
		Matrix<Nucleotide> nuc = FASTA_to_nucleotide(fa,empirical_nucleotide_frequencies,isBLC);
		// Identify and count unique patterns
		vector<string> pat;				// Pattern as string of AGCTNs
		vector<int> pat1, cpat, ipat;	// First example of each pattern, number of sites with that pattern, the pattern at each (compatible) site (-1 otherwise)
		vector<bool> nuc_ispoly(nuc.ncols(),true);
		find_alignment_patterns(nuc,nuc_ispoly,pat,pat1,cpat,ipat);
		// Storage for the MLE of the nucleotide sequence at every node
		Matrix<Nucleotide> node_nuc;
		// Begin by computing the joint maximum likelihood ancestral sequences
		mydouble ML = maximum_likelihood_ancestral_sequences(nuc,ctree,kappa,empirical_nucleotide_frequencies,pat1,cpat,node_nuc);
		
		cout << "BRANCH LENGTH CORRECTION:" << endl;
		cout << "Analysing " << nBLC << " sites" << endl;
		// Report the estimated equilibrium frequencies
		cout << "Empirical nucleotide frequencies:   A " << round(1000*empirical_nucleotide_frequencies[Adenine])/10 << "%   C " << round(1000*empirical_nucleotide_frequencies[Cytosine])/10;
		cout << "%   G " << round(1000*empirical_nucleotide_frequencies[Guanine])/10 << "%   T " << round(1000*empirical_nucleotide_frequencies[Thymine])/10 << "%" << endl;
		
		if(OLD_JOINT_BRANCH_PARAM) {
			// Prepare to correct branch lengths
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) under the ClonalFrame model
			ClonalFrameFunction cff(ctree,node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,EXCESS_DIVERGENCE_MODEL,MULTITHREAD,root_node);
			vector<double> param(3+root_node);
			param[0] = log10(initial_rho_over_theta); param[1] = log10(initial_mean_import_length); param[2] = log10(initial_import_divergence);
			for(i=0;i<root_node;i++) param[3+i] = log10(ctree.node[i].edge_time); 
			
			// Maximize the likelihood of the parameters
			Powell Pow(cff);
			Pow.coutput = Pow.brent.coutput = SHOW_PROGRESS;
			Pow.TOL = brent_tolerance;
			param = Pow.minimize(param,powell_tolerance);
			// Ensure importation status is updated correctly
			cff.f(param);
			
			cout << "Maximum log-likelihood of " << -Pow.function_minimum << " found in " << cff.neval << " iterations" << endl;
			cout << "rho/theta         = " << pow(10.,param[0]) << endl;
			cout << "Import length     = " << pow(10.,param[1]) << endl;
			cout << "Import divergence = " << pow(10.,param[2]) << endl;
			for(i=0;i<root_node;i++) {
				cout << "Branch length " << ctree_node_labels[i] << " = " << pow(10.,param[3+i]) << endl;
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = pow(10.,param[3+i]);
			}

			// Output the importation status
			write_importation_status(cff.is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);	
			cout << "Wrote inferred importation status to " << import_out_file << endl;
		} else if(JOINT_BRANCH_PARAM) {
			cout << "Beginning parameter estimation" << endl;
			vector< vector<ImportationState> > is_imported(root_node);
			// Calculate the expected number of substitutions per branch according to the ancestral state reconstruction
			vector<double> substitutions_per_branch(root_node,0);
			for(i=0;i<root_node;i++) {
				double pd = 1.0, pd_den = 2.0;
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				int j,k;
				for(j=0,k=0;j<isBLC.size();j++) {
					if(isBLC[j]) {
						Nucleotide dec = node_nuc[dec_id][ipat[k]];
						Nucleotide anc = node_nuc[anc_id][ipat[k]];
						if(dec!=anc) ++pd;
						++pd_den;
						++k;
					}
				}
				substitutions_per_branch[i] = pd/pd_den;
			}
			// Prepare to correct branch lengths
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) under the ClonalFrame model
			// Minimum branch length
			const double min_branch_length = global_min_branch_length;
			const bool USE_VITERBI = false;
			ClonalFrameFunctionJoint cff(USE_VITERBI,ctree,node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,EXCESS_DIVERGENCE_MODEL,is_imported,substitutions_per_branch,min_branch_length,SHOW_PROGRESS,brent_tolerance,powell_tolerance,root_node);
			vector<double> param(3);
			param[0] = log10(initial_rho_over_theta); param[1] = log10(initial_mean_import_length); param[2] = log10(initial_import_divergence);
			// Maximize the likelihood of the parameters
			Powell Pow(cff);
			Pow.coutput = Pow.brent.coutput = SHOW_PROGRESS;
			Pow.TOL = brent_tolerance;
			param = Pow.minimize(param,powell_tolerance);
			// Ensure importation status is updated correctly
			cff.f(param);
			
			cout << "Maximum log-likelihood of " << -Pow.function_minimum << " found in " << cff.neval << " iterations" << endl;
			cout << "rho/theta         = " << pow(10.,param[0]) << endl;
			cout << "Import length     = " << pow(10.,param[1]) << endl;
			cout << "Import divergence = " << pow(10.,param[2]) << endl;
			for(i=0;i<root_node;i++) {
				double blen = cff.branch_length[i];
				if(blen<min_branch_length) blen = min_branch_length;
				cout << "Branch length " << ctree_node_labels[i] << " = " << blen << endl;
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = blen;
			}
			
			// Output the importation status
			write_importation_status(cff.is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			cout << "Wrote inferred importation status to " << import_out_file << endl;
		} else if(RHO_PER_BRANCH_NO_LRT) {
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) AND recombination parameters under the ClonalFrame model
			// SUBJECT to the constraints that the importation state have frequency < 0.5 AND the recombination divergence exceeds the branch length divergence
			// Now estimate branch lengths
			cout << "Beginning branch optimization. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L   maximum log-likelihood per branch" << endl;
			cout << "M   expected number of mutations per branch                  (> 0)" << endl;
			cout << "R   expected number of recombination imports per branch      (> 0)" << endl;
			cout << "F   ratio of imported versus unimported sites                (0-1)" << endl;
			cout << "E   excess substitutions imported by recombination           (> M)" << endl;			
			double ML = 0.0;
			vector< vector<ImportationState> > is_imported(root_node);
			for(i=0;i<root_node;i++) {
				// Crudely re-estimate branch length
				double pd = 1.0, pd_den = 2.0;
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				int j,k;
				for(j=0,k=0;j<isBLC.size();j++) {
					if(isBLC[j]) {
						Nucleotide dec = node_nuc[dec_id][ipat[k]];
						Nucleotide anc = node_nuc[anc_id][ipat[k]];
						if(dec!=anc) ++pd;
						++pd_den;
						++k;
					}
				}
				const double initial_branch_length = pd/pd_den;
				// Initial values for the other parameters
				const double initial_import_ratio = 1.0/initial_mean_import_length;			// constrained to be less than 1
				// Minimum branch length
				const double min_branch_length = global_min_branch_length;
				ClonalFrameRhoPerBranchFunction cff(ctree.node[i],node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,EXCESS_DIVERGENCE_MODEL,MULTITHREAD,is_imported[i],initial_branch_length,min_branch_length);
				// Setup optimization function
				Powell Pow(cff);
				Pow.coutput = Pow.brent.coutput = SHOW_PROGRESS;
				Pow.TOL = brent_tolerance;
				// Initially, estimate parameters with total substitution rate fixed at the crude estimate
				vector<double> param(3);
				param[0] = log10(initial_rho_over_theta); param[1] = log10(initial_import_ratio/(1.0-initial_import_ratio)); param[2] = log10(initial_import_divergence);
				param = Pow.minimize(param,powell_tolerance);
				// Then refine
				const double import_ratio = 1.0/(1.0+pow(10.,-param[1]));
				const double import_divergence = pow(10.,param[2]);
				param.push_back(log10(initial_branch_length/(1.0+import_ratio/(1.0+import_ratio)*(2.0+import_divergence))));
				param = Pow.minimize(param,powell_tolerance);
				// Ensure importation status is updated correctly
				double final_branch_length = pow(10.,param[3]);
				if(final_branch_length<min_branch_length) final_branch_length = min_branch_length;
				const double final_rho_over_theta = pow(10.,param[0]);
				const double final_mean_import_length = (1.0/(1.0+pow(10.,-param[1])))/final_branch_length/final_rho_over_theta;
				const double final_import_divergence = final_branch_length*(2.0+pow(10.,param[2]));
				maximum_likelihood_ClonalFrame_branch_allsites(dec_id, anc_id, node_nuc, isBLC, ipat, kappa, empirical_nucleotide_frequencies, final_branch_length, final_rho_over_theta, final_mean_import_length, final_import_divergence, is_imported[i]);
				// Update branch length in the tree
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = final_branch_length;
				cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L = " << -Pow.function_minimum << " M = " << final_branch_length << " R = " << pow(10.,param[0]+param[3]) << " F = " << 1.0/(1.0+pow(10.,-param[1])) << " E = " << final_branch_length*(1.0+pow(10.,param[2])) << endl;
				ML += -Pow.function_minimum;
			}
			cout << "Log-likelihood after branch optimization is " << ML << endl;
						
			// Output the importation status
			write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			cout << "Wrote inferred importation status to " << import_out_file << endl;
		} else if(RESCALE_NO_RECOMBINATION) {
			// Rescale the branch lengths using given sites without a model of recombination
			cout << "Beginning branch optimization. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L   maximum log-likelihood per branch" << endl;
			cout << "M   corrected branch length/expected number of mutations     (> 0)" << endl;
			double ML = 0.0;
			for(i=0;i<root_node;i++) {
				// Crudely re-estimate branch length
				double pd = 1.0, pd_den = 2.0;
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				int j,k;
				for(j=0,k=0;j<isBLC.size();j++) {
					if(isBLC[j]) {
						Nucleotide dec = node_nuc[dec_id][ipat[k]];
						Nucleotide anc = node_nuc[anc_id][ipat[k]];
						if(dec!=anc) ++pd;
						++pd_den;
						++k;
					}
				}
				const double initial_branch_length = pd/pd_den;
				// Minimum branch length
				const double min_branch_length = global_min_branch_length;
				ClonalFrameRescaleBranchFunction cff(ctree.node[i],node_nuc,pat1,cpat,kappa,empirical_nucleotide_frequencies,MULTITHREAD,initial_branch_length,min_branch_length);
				// Setup optimization function
				Powell Pow(cff);
				Pow.coutput = Pow.brent.coutput = SHOW_PROGRESS;
				Pow.TOL = brent_tolerance;
				// Estimate parameter
				vector<double> param(1,log10(initial_branch_length));
				param = Pow.minimize(param,powell_tolerance);
				double final_branch_length = pow(10.,param[0]);
				if(final_branch_length<min_branch_length) final_branch_length = min_branch_length;
				// Update branch length in the tree
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = final_branch_length;
				cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L = " << -Pow.function_minimum << " M = " << final_branch_length << endl;
				ML += -Pow.function_minimum;
			}
			cout << "Log-likelihood after branch optimization is " << ML << endl;
		} else if(RHO_PER_BRANCH) {
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) AND recombination parameters under the ClonalFrame model
			// SUBJECT to the constraints that the importation state have frequency < 0.5 AND the recombination divergence exceeds the branch length divergence
			// Also impose a likelihood ratio test, so only report simpler model if no signifcant improvement with the recombination model.
			// Now estimate branch lengths
			cout << "Beginning branch optimization. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L0  maximum log-likelihood per branch (no recombination model)" << endl;
			cout << "L   maximum log-likelihood per branch" << endl;
			cout << "*   preferred model based on likelihood ratio test" << endl;
			cout << "M   expected number of mutations per branch                  (> 0)" << endl;
			cout << "R   expected number of recombination imports per branch      (> 0)" << endl;
			cout << "F   ratio of imported versus unimported sites                (0-1)" << endl;
			cout << "E   excess substitutions imported by recombination           (> M)" << endl;			
			double ML = 0.0;
			vector< vector<ImportationState> > is_imported(root_node);
			for(i=0;i<root_node;i++) {
				// Crudely re-estimate branch length
				double pd = 1.0, pd_den = 2.0;
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				int j,k;
				for(j=0,k=0;j<isBLC.size();j++) {
					if(isBLC[j]) {
						Nucleotide dec = node_nuc[dec_id][ipat[k]];
						Nucleotide anc = node_nuc[anc_id][ipat[k]];
						if(dec!=anc) ++pd;
						++pd_den;
						++k;
					}
				}
				const double initial_branch_length = pd/pd_den;
				// Initial values for the other parameters
				const double initial_import_ratio = 1.0/initial_mean_import_length;			// constrained to be less than 1
				// Minimum branch length
				const double min_branch_length = global_min_branch_length;
				// Object for the no-recombination model
				ClonalFrameRescaleBranchFunction cff0(ctree.node[i],node_nuc,pat1,cpat,kappa,empirical_nucleotide_frequencies,MULTITHREAD,initial_branch_length,min_branch_length);
				// Setup optimization function
				Powell Pow0(cff0);
				Pow0.coutput = Pow0.brent.coutput = SHOW_PROGRESS;
				Pow0.TOL = brent_tolerance;
				// Initially, estimate parameter for the no-recombination model
				vector<double> param0(1,log10(initial_branch_length)), param(3);
				param0 = Pow0.minimize(param0,powell_tolerance);
				// Object for the per-branch recombination model: initial branch length set from no-recombination model estimate
				ClonalFrameRhoPerBranchFunction cff(ctree.node[i],node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,EXCESS_DIVERGENCE_MODEL,MULTITHREAD,is_imported[i],pow(10.,param0[0]),min_branch_length);
				// Setup optimization function
				Powell Pow(cff);
				Pow.coutput = Pow.brent.coutput = SHOW_PROGRESS;
				Pow.TOL = brent_tolerance;
				// Now estimate parameters for the recombination model with total substitution rate fixed at this estimate
				param[0] = log10(initial_rho_over_theta); param[1] = log10(initial_import_ratio/(1.0-initial_import_ratio)); param[2] = log10(initial_import_divergence);
				param = Pow.minimize(param,powell_tolerance);
				// Then refine
				const double import_ratio = 1.0/(1.0+pow(10.,-param[1]));
				const double import_divergence = pow(10.,param[2]);
				param.push_back(log10(pow(10.,param0[0])/(1.0+import_ratio/(1.0+import_ratio)*(2.0+import_divergence))));
				param = Pow.minimize(param,powell_tolerance);
				// Test if the recombination model is a significant improvement (2 log-likelihoods per additional parameter)
				const bool LRT_pass = (-Pow.function_minimum>-Pow0.function_minimum+6.0);
				// Ensure importation status is updated correctly
				double final_branch_length = (LRT_pass) ? pow(10.,param[3]) : pow(10.,param0[0]);
				if(final_branch_length<min_branch_length) final_branch_length = min_branch_length;
				const double final_rho_over_theta = (LRT_pass) ? pow(10.,param[0]) : 1e-20;
				const double final_mean_import_length = (LRT_pass) ? (1.0/(1.0+pow(10.,-param[1])))/final_branch_length/final_rho_over_theta : 1e20;
				const double final_import_divergence = (LRT_pass) ? final_branch_length*(2.0+pow(10.,param[2])) : 1e20;
				maximum_likelihood_ClonalFrame_branch_allsites(dec_id, anc_id, node_nuc, isBLC, ipat, kappa, empirical_nucleotide_frequencies, final_branch_length, final_rho_over_theta, final_mean_import_length, final_import_divergence, is_imported[i]);
				// Update branch length in the tree
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = final_branch_length;
				if(LRT_pass) {
					cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L0 = " << -Pow0.function_minimum << "  L = " << -Pow.function_minimum << "* M = " << final_branch_length << " R = " << pow(10.,param[0]+param[3]) << " F = " << 1.0/(1.0+pow(10.,-param[1])) << " E = " << final_branch_length*(1.0+pow(10.,param[2])) << endl;
					ML += -Pow.function_minimum;
				} else {
					cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L0 = " << -Pow0.function_minimum << "* L = " << -Pow.function_minimum << "  M = " << final_branch_length << endl;
					ML += -Pow0.function_minimum;
				}
			}
			cout << "Log-likelihood after branch optimization is " << ML << endl;
			
			// Output the importation status
			write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			cout << "Wrote inferred importation status to " << import_out_file << endl;
		} else if(LAPLACE_APPROX) {
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) AND recombination parameters under the ClonalFrame model
			// using the Laplace approximation to the likelihood and a "driving prior" to assist with the optimization (this ensures a proper posterior)
			cout << "Beginning branch optimization. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L   maximum unnormalized log-posterior per branch" << endl;
			cout << "R   rho/theta per branch                                     (> 0)" << endl;
			cout << "I   mean DNA import length per branch                        (> 0)" << endl;
			cout << "D   divergence of DNA imported by recombination              (> 0)" << endl;			
			cout << "M   expected number of mutations per branch                  (> 0)" << endl;
			double ML = 0.0;
			vector< vector<ImportationState> > is_imported(root_node);
			vector< vector<double> > laplaceMLE(0);
			vector< Matrix<double> > laplaceQ(0);
			laplaceMLE = vector< vector<double> >(root_node,vector<double>(4));
			laplaceQ = vector< Matrix<double> >(root_node,Matrix<double>(4,4));
			int j,k;
			for(i=0;i<root_node;i++) {
				// Crudely re-estimate branch length: use this as the mean of the prior on branch length ????
				double pd = 1.0, pd_den = 2.0;
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				for(j=0,k=0;j<isBLC.size();j++) {
					if(isBLC[j]) {
						Nucleotide dec = node_nuc[dec_id][ipat[k]];
						Nucleotide anc = node_nuc[anc_id][ipat[k]];
						if(dec!=anc) ++pd;
						++pd_den;
						++k;
					}
				}
				const double initial_branch_length = pd/pd_den;
				// Object for the per-branch recombination model: initial branch length set from no-recombination model estimate
				const int PARAMETERIZATION = 3;		// Do not let this take other values without reviewing code
				ClonalFrameLaplacePerBranchFunction cff(ctree.node[i],node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,MULTITHREAD,is_imported[i],driving_prior_mean,driving_prior_precision);
				cff.parameterization = PARAMETERIZATION;
				// Setup optimization function
				Powell Pow(cff);
				Pow.coutput = Pow.brent.coutput = SHOW_PROGRESS;
				Pow.TOL = brent_tolerance;
				// Now estimate parameters for the recombination model starting at the mean of the prior (or initial values), except the branch length
				vector<double> param;
				if(initial_values.size()==0) {
					param = driving_prior_mean;
					param[3] = log10(initial_branch_length);
				} else {
					param = initial_values;
					param.push_back(log10(initial_branch_length));
				}
				if(PARAMETERIZATION==3) {
					param = cff.convert_parameterization_1_to_3(param,global_min_branch_length);
				}
				clock_t pow_start_time = clock();
				int neval = cff.neval;
				param = Pow.minimize(param,powell_tolerance);
				if(PARAMETERIZATION==3) {
					param = cff.convert_parameterization_3_to_0(param);
				}
				cout << "Powell gave param = " << param[0] << " " << param[1] << " " << param[2] << " " << param[3] << " post = " << -Pow.function_minimum << " in " << (double)(clock()-pow_start_time)/CLOCKS_PER_SEC << " s and " << cff.neval-neval << " evaluations" << endl;
				if(false) {
					// Attempt to refine using BFGS
					param = driving_prior_mean;
					param[3] = log10(initial_branch_length);
					BFGS bfgs(cff);
					bfgs.coutput = SHOW_PROGRESS;
					clock_t bfgs_start_time = clock();
					neval = cff.neval;
					param = bfgs.minimize(param,powell_tolerance);
					cout << "BFGS gave param = " << param[0] << " " << param[1] << " " << param[2] << " " << param[3] << " post = " << -bfgs.function_minimum << " in " << (double)(clock()-bfgs_start_time)/CLOCKS_PER_SEC << " s and " << cff.neval-neval << " evaluations" << endl;
					// Get the approximate inverse Hessian
					// laplaceQ[i] = bfgs.hessin;
					cout << "BFGS gave marginal st devs = " << sqrt(laplaceQ[i][0][0]) << " " << sqrt(laplaceQ[i][1][1]) << " " << sqrt(laplaceQ[i][2][2]) << " " << sqrt(laplaceQ[i][3][3]) << endl;
				}
				// Approximate the likelihood by a multivariate Gaussian
				laplaceMLE[i] = param;
				// Interval over which to numerically compute second derivatives
				// Assumes a log-likelihood accuracy calculation of 0.001 and a curvature scale of 1
				const double h = 0.1;
				vector<double> paramQ;
				double calcQ;
				// The maximum log-likelihood
				if(PARAMETERIZATION!=3) {
					const double calcQ0 = -cff.f(param);
					for(j=0;j<4;j++) {
						for(k=0;k<j;k++) {
							paramQ = param; paramQ[j] += h; paramQ[k] += h;
							calcQ = -cff.f(paramQ);
							paramQ = param; paramQ[j] += h; paramQ[k] -= h;
							calcQ -= -cff.f(paramQ);
							paramQ = param; paramQ[j] -= h; paramQ[k] += h;
							calcQ -= -cff.f(paramQ);
							paramQ = param; paramQ[j] -= h; paramQ[k] -= h;
							calcQ += -cff.f(paramQ);
							laplaceQ[i][j][k] = laplaceQ[i][k][j] = -calcQ/4.0/h/h;
						}
						paramQ = param; paramQ[j] += 2.0*h;
						calcQ = -cff.f(paramQ);
						paramQ = param; paramQ[j] -= 2.0*h;
						calcQ += -cff.f(paramQ);
						calcQ -= 2.0*calcQ0;
						laplaceQ[i][j][j] = -calcQ/4.0/h/h;
					}
				} else {
					const double calcQ0 = -cff.f(cff.convert_parameterization_0_to_3(param,global_min_branch_length));
					for(j=0;j<4;j++) {
						for(k=0;k<j;k++) {
							paramQ = param; paramQ[j] += h; paramQ[k] += h;
							calcQ = -cff.f(cff.convert_parameterization_0_to_3(paramQ,global_min_branch_length));
							paramQ = param; paramQ[j] += h; paramQ[k] -= h;
							calcQ -= -cff.f(cff.convert_parameterization_0_to_3(paramQ,global_min_branch_length));
							paramQ = param; paramQ[j] -= h; paramQ[k] += h;
							calcQ -= -cff.f(cff.convert_parameterization_0_to_3(paramQ,global_min_branch_length));
							paramQ = param; paramQ[j] -= h; paramQ[k] -= h;
							calcQ += -cff.f(cff.convert_parameterization_0_to_3(paramQ,global_min_branch_length));
							laplaceQ[i][j][k] = laplaceQ[i][k][j] = -calcQ/4.0/h/h;
						}
						paramQ = param; paramQ[j] += 2.0*h;
						calcQ = -cff.f(cff.convert_parameterization_0_to_3(paramQ,global_min_branch_length));
						paramQ = param; paramQ[j] -= 2.0*h;
						calcQ += -cff.f(cff.convert_parameterization_0_to_3(paramQ,global_min_branch_length));
						calcQ -= 2.0*calcQ0;
						laplaceQ[i][j][j] = -calcQ/4.0/h/h;
					}
				}
				// Ensure importation status is updated at the MAP parameter estimate (for now, this is affected by the driving prior)
				const double final_rho_over_theta = pow(10.,param[0]);
				const double final_mean_import_length = pow(10.,param[1]);
				const double final_import_divergence = pow(10.,param[2]);
				const double final_branch_length = (PARAMETERIZATION==1) ? pow(10.,param[3])/(1+final_rho_over_theta*final_mean_import_length*(final_import_divergence-pow(10.,param[3]))) : pow(10.,param[3]);
				maximum_likelihood_ClonalFrame_branch_allsites(dec_id, anc_id, node_nuc, isBLC, ipat, kappa, empirical_nucleotide_frequencies, final_branch_length, final_rho_over_theta, final_mean_import_length, final_import_divergence, is_imported[i]);
				// Update branch length in the tree
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = final_branch_length;
				// Output results to screen
				cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L = " << -Pow.function_minimum << " R = " << final_rho_over_theta << " I = " << final_mean_import_length << " D = " << final_import_divergence << " M = " << final_branch_length << endl;
				ML += -Pow.function_minimum;
			}
			cout << "Unnormalized log-posterior after branch optimization is " << ML << endl;
			
			// Output the importation status
			write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			cout << "Wrote inferred importation status to " << import_out_file << endl;
			// Output the Laplace approximation
			ofstream lout(laplace_out_file.c_str());
			lout << "Branch";
			for(j=0;j<4;j++) lout << "\t" << "P" << j;
			for(j=0;j<4;j++) {
				for(k=0;k<4;k++) {
					lout << "\t" << "Q" << j << k;
				}
			}
			lout << endl;
			for(i=0;i<root_node;i++) {
				lout << ctree_node_labels[i];
				for(j=0;j<4;j++) lout << "\t" << laplaceMLE[i][j];
				for(j=0;j<4;j++) {
					for(k=0;k<4;k++) {
						lout << "\t" << laplaceQ[i][j][k];
					}
				}
				lout << endl;
			}
			lout.close();
			cout << "Wrote Laplace approximation information to " << laplace_out_file << endl;
		} else if(MCMC_PER_BRANCH) {
			// For a given branch, infer by MCMC importation state (unimported vs imported) AND recombination parameters under the ClonalFrame model
			cout << "Beginning branch optimization. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L   maximum unnormalized log-posterior per branch" << endl;
			cout << "R   rho/theta per branch                                     (> 0)" << endl;
			cout << "I   mean DNA import length per branch                        (> 0)" << endl;
			cout << "D   divergence of DNA imported by recombination              (> 0)" << endl;			
			cout << "M   expected number of mutations per branch                  (> 0)" << endl;
			vector< vector<ImportationState> > is_imported(root_node);
			// Open files for writing
			vector< ofstream* > mout(4);
			vector< ofstream* > pout(4);
			mout[0] = new ofstream(mcmc0_out_file.c_str());
			mout[1] = new ofstream(mcmc1_out_file.c_str());
			mout[2] = new ofstream(mcmc2_out_file.c_str());
			mout[3] = new ofstream(mcmc3_out_file.c_str());
			pout[0] = new ofstream(prop0_out_file.c_str());
			pout[1] = new ofstream(prop1_out_file.c_str());
			pout[2] = new ofstream(prop2_out_file.c_str());
			pout[3] = new ofstream(prop3_out_file.c_str());
			ofstream lout(loglik_out_file.c_str());
			ofstream lpout(loglik_prop_out_file.c_str());
			// For each branch
			int j,k;
			for(i=0;i<root_node;i++) {
				// Crudely re-estimate branch length
				double pd = 1.0, pd_den = 2.0;
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				for(j=0,k=0;j<isBLC.size();j++) {
					if(isBLC[j]) {
						Nucleotide dec = node_nuc[dec_id][ipat[k]];
						Nucleotide anc = node_nuc[anc_id][ipat[k]];
						if(dec!=anc) ++pd;
						++pd_den;
						++k;
					}
				}
				const double initial_branch_length = pd/pd_den;
				// Object for the per-branch recombination model
				ClonalFrameLaplacePerBranchFunction cff(ctree.node[i],node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,MULTITHREAD,is_imported[i],driving_prior_mean,driving_prior_precision);
				// 4/6/14: start MCMC at prior mean except 4th element which is now the branch length (not expected divergence at unimported sites)
				vector<double> param = driving_prior_mean;
				param[3] = log10(initial_branch_length);
				// Output preamble
				for(j=0;j<4;j++) *mout[j] << ctree.node[i].id;
				for(j=0;j<4;j++) *pout[j] << ctree.node[i].id;
				lout << ctree.node[i].id;
				lpout << ctree.node[i].id;
				// Do MCMC
				int iter;
				const int niter = 1000;
				double loglik = -cff.f(param);
				double proposal_sd[4] = {0.5,0.5,0.5,0.02};
				for(iter=0;iter<niter;iter++) {
					for(j=0;j<4;j++) {
						vector<double> new_param = param;
						new_param[j] += ran.normal(0,proposal_sd[j]);
						const double new_loglik = -cff.f(new_param);
						const double alpha = new_loglik-loglik;
						const bool accept = alpha>=0.0 || alpha>=log(ran.U());
						if(accept) {
							// Accept
							param = new_param;
							loglik = new_loglik;
						} else {
							// Reject (do nothing)
						}
						// Record
						if(iter>=0 && (iter%1)==0) {
							*mout[j] << "\t" << param[j];
							*pout[j] << "\t" << new_param[j];
							lout << "\t" << loglik;
							lpout << "\t" << new_loglik;
						}
					}
				}
				// End of line
				for(j=0;j<4;j++) *mout[j] << endl;
				for(j=0;j<4;j++) *pout[j] << endl;
				lout << endl;
				lpout << endl;
				
				cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << endl;
/*				// NOT YET IMPLEMENTED
				// Ensure importation status is updated correctly
				double final_branch_length = (LRT_pass) ? pow(10.,param[3]) : pow(10.,param0[0]);
				if(final_branch_length<min_branch_length) final_branch_length = min_branch_length;
				const double final_rho_over_theta = (LRT_pass) ? pow(10.,param[0]) : 1e-20;
				const double final_mean_import_length = (LRT_pass) ? (1.0/(1.0+pow(10.,-param[1])))/final_branch_length/final_rho_over_theta : 1e20;
				const double final_import_divergence = (LRT_pass) ? final_branch_length*(2.0+pow(10.,param[2])) : 1e20;
				maximum_likelihood_ClonalFrame_branch_allsites(dec_id, anc_id, node_nuc, isBLC, ipat, kappa, empirical_nucleotide_frequencies, final_branch_length, final_rho_over_theta, final_mean_import_length, final_import_divergence, is_imported[i]);
				// Update branch length in the tree
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = final_branch_length;
				if(LRT_pass) {
					cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L0 = " << -Pow0.function_minimum << "  L = " << -Pow.function_minimum << "* M = " << final_branch_length << " R = " << pow(10.,param[0]+param[3]) << " F = " << 1.0/(1.0+pow(10.,-param[1])) << " E = " << final_branch_length*(1.0+pow(10.,param[2])) << endl;
					ML += -Pow.function_minimum;
				} else {
					cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L0 = " << -Pow0.function_minimum << "* L = " << -Pow.function_minimum << "  M = " << final_branch_length << endl;
					ML += -Pow0.function_minimum;
				}
*/
			}
			// Close files
			for(j=0;j<4;j++) {
				mout[j]->close();
				delete mout[j];
				pout[j]->close();
				delete pout[j];
			}
			lout.close();
			lpout.close();
			
			// NOT YET IMPLEMENTED
/*			// Output the importation status
			write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			cout << "Wrote inferred importation status to " << import_out_file << endl;
 */
		} else if(MCMC_JOINT) {
			// For all branches jointly, infer by MCMC recombination parameters under the ClonalFrame model
			cout << "Beginning MCMC. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L   maximum unnormalized log-posterior per branch" << endl;
			cout << "R   rho/theta per branch                                     (> 0)" << endl;
			cout << "I   mean DNA import length per branch                        (> 0)" << endl;
			cout << "D   divergence of DNA imported by recombination              (> 0)" << endl;			
			cout << "M   expected number of mutations per branch                  (> 0)" << endl;
			vector< vector<ImportationState> > is_imported(root_node);
			// Open files for writing
			vector< ofstream* > mout(4);
			vector< ofstream* > pout(4);
			mout[0] = new ofstream(mcmc0_out_file.c_str());
			mout[1] = new ofstream(mcmc1_out_file.c_str());
			mout[2] = new ofstream(mcmc2_out_file.c_str());
			mout[3] = new ofstream(mcmc3_out_file.c_str());
			pout[0] = new ofstream(prop0_out_file.c_str());
			pout[1] = new ofstream(prop1_out_file.c_str());
			pout[2] = new ofstream(prop2_out_file.c_str());
			pout[3] = new ofstream(prop3_out_file.c_str());
			ofstream lout(loglik_out_file.c_str());
			ofstream lpout(loglik_prop_out_file.c_str());
			// Minimum branch length
			const double min_branch_length = global_min_branch_length;
			// For each branch find good starting values for the branch lengths
			int j,k;
			vector<double> initial_branch_length(root_node);
			for(i=0;i<root_node;i++) {
				// Crudely re-estimate branch length (if !MCMC_INFER_BRANCH_LENGTHS these values will be fixed)
				double pd = 1.0, pd_den = 2.0;
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				for(j=0,k=0;j<isBLC.size();j++) {
					if(isBLC[j]) {
						Nucleotide dec = node_nuc[dec_id][ipat[k]];
						Nucleotide anc = node_nuc[anc_id][ipat[k]];
						if(dec!=anc) ++pd;
						++pd_den;
						++k;
					}
				}
				initial_branch_length[i] = pd/pd_den;
			}
			// Initialize the parameters at the driving prior mean or, in the case of the branch lengths, their crude estimates
			vector<double> param(root_node+3), new_param;
			for(i=0;i<3;i++) param[i] = driving_prior_mean[i];
			for(i=0;i<root_node;i++) param[3+i] = log10(initial_branch_length[i]);
			vector<double> partial_post(root_node), new_partial_post(root_node);
			// Initialize the log likelihood class
			ClonalFrameMCMCJointFunction cff(ctree.node,node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,MULTITHREAD,is_imported,min_branch_length,driving_prior_mean,driving_prior_precision,true,root_node);
			cff.parameterization = (MCMC_INFER_BRANCH_LENGTHS) ? 1 : 0;
			// Output preamble to files
			*mout[0] << "rho_over_theta" << endl;
			*mout[1] << "import_length" << endl;
			*mout[2] << "import_divergence" << endl;
			for(i=0;i<root_node;i++) {
				if(i>0) *mout[3] << "\t";
				*mout[3] << ctree.node[i].id;
			}
			*mout[3] << endl;
			*pout[0] << "rho_over_theta" << endl;
			*pout[1] << "import_length" << endl;
			*pout[2] << "import_divergence" << endl;
			for(i=0;i<root_node;i++) {
				if(i>0) *pout[3] << "\t";
				*pout[3] << ctree.node[i].id;
			}
			*pout[3] << endl;
			lout << "posterior" << endl;
			lpout << "posterior" << endl;
			// Do MCMC			
			int iter;
			const int niter = 1000;
			double loglik = cff.log_posterior(param,partial_post), new_loglik;
			vector<double> proposal_sd(param.size(),0.1);
			const int nmix_param = (MCMC_INFER_BRANCH_LENGTHS) ? param.size() : 3;
			for(iter=0;iter<niter;iter++) {
				for(j=0;j<nmix_param;j++) {
					new_param = param;
					new_param[j] += ran.normal(0,proposal_sd[j]);
					if(j<3) {
						new_loglik = cff.log_posterior(new_param,new_partial_post);
					} else {
						new_partial_post = partial_post;
						new_loglik = cff.log_posterior(new_param,new_partial_post,j-3);
					}
					const double alpha = new_loglik-loglik;
					const bool accept = alpha>=0.0 || alpha>=log(ran.U());
					if(accept) {
						// Accept
						param = new_param;
						loglik = new_loglik;
						partial_post = new_partial_post;
					} else {
						// Reject (do nothing)
					}
				}
				// Record
				if(iter>=0 && (iter%1)==0) {
					*mout[0] << param[0] << endl;
					*mout[1] << param[1] << endl;
					*mout[2] << param[2] << endl;
					for(i=0;i<root_node;i++) {
						if(i>0) *mout[3] << "\t";
						*mout[3] << param[3+i];
					}
					*mout[3] << endl;
					lout << loglik << endl;
				}
			}
			
			//cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << endl;
			/*				// NOT YET IMPLEMENTED
			 // Ensure importation status is updated correctly
			 double final_branch_length = (LRT_pass) ? pow(10.,param[3]) : pow(10.,param0[0]);
			 if(final_branch_length<min_branch_length) final_branch_length = min_branch_length;
			 const double final_rho_over_theta = (LRT_pass) ? pow(10.,param[0]) : 1e-20;
			 const double final_mean_import_length = (LRT_pass) ? (1.0/(1.0+pow(10.,-param[1])))/final_branch_length/final_rho_over_theta : 1e20;
			 const double final_import_divergence = (LRT_pass) ? final_branch_length*(2.0+pow(10.,param[2])) : 1e20;
			 maximum_likelihood_ClonalFrame_branch_allsites(dec_id, anc_id, node_nuc, isBLC, ipat, kappa, empirical_nucleotide_frequencies, final_branch_length, final_rho_over_theta, final_mean_import_length, final_import_divergence, is_imported[i]);
			 // Update branch length in the tree
			 // Note this is unsafe in general because the corresponding node times are not adjusted
			 ctree.node[i].edge_time = final_branch_length;
			 if(LRT_pass) {
			 cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L0 = " << -Pow0.function_minimum << "  L = " << -Pow.function_minimum << "* M = " << final_branch_length << " R = " << pow(10.,param[0]+param[3]) << " F = " << 1.0/(1.0+pow(10.,-param[1])) << " E = " << final_branch_length*(1.0+pow(10.,param[2])) << endl;
			 ML += -Pow.function_minimum;
			 } else {
			 cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L0 = " << -Pow0.function_minimum << "* L = " << -Pow.function_minimum << "  M = " << final_branch_length << endl;
			 ML += -Pow0.function_minimum;
			 }
			 */
			// Close files
			for(j=0;j<4;j++) {
				mout[j]->close();
				delete mout[j];
				pout[j]->close();
				delete pout[j];
			}
			lout.close();
			lpout.close();

			// NOT YET IMPLEMENTED
			/*			// Output the importation status
			 write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			 cout << "Wrote inferred importation status to " << import_out_file << endl;
			 */
		} else if(GRID_APPROX) {
			// For a given branch, approximate the mean and covariance matrix of the posterior distribution by evaluating over a grid of points
			// using the driving prior
			const int nparam = 4;
			const int ngrid = (int)grid_approx;
			Matrix<double> grid(nparam,ngrid);
			int i,j,k,l,m;
			for(i=0;i<nparam;i++) {
				// Create a uniform grid within the marginal 95% prior quantiles of each parameter
				const double grid_beg = driving_prior_mean[i]-2.0/sqrt(driving_prior_precision[i]);
				const double grid_end = driving_prior_mean[i]+2.0/sqrt(driving_prior_precision[i]);
				const double grid_inc = (grid_end-grid_beg)/((double)ngrid-1.0);
				for(j=0;j<ngrid;j++) {
					grid[i][j] = grid_beg + ((double)j)*grid_inc;
				}
			}
			cout << "Beginning grid approximation. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L   maximum unnormalized log-posterior per branch" << endl;
			cout << "R   rho/theta per branch                                     (> 0)" << endl;
			cout << "I   mean DNA import length per branch                        (> 0)" << endl;
			cout << "D   divergence of DNA imported by recombination              (> 0)" << endl;			
			cout << "M   expected number of mutations per branch                  (> 0)" << endl;
			vector< vector<ImportationState> > is_imported(root_node);
			vector< vector<double> > Egrid(0);
			vector< Matrix<double> > Vgrid(0);
			Egrid = vector< vector<double> >(root_node,vector<double>(nparam,0.0));
			Vgrid = vector< Matrix<double> >(root_node,Matrix<double>(nparam,nparam,0.0));
			// Cycle over branches
			for(i=0;i<root_node;i++) {
				// Crudely re-estimate branch length: use this as the mean of the prior on branch length ????
				double pd = 1.0, pd_den = 2.0;
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				for(j=0,k=0;j<isBLC.size();j++) {
					if(isBLC[j]) {
						Nucleotide dec = node_nuc[dec_id][ipat[k]];
						Nucleotide anc = node_nuc[anc_id][ipat[k]];
						if(dec!=anc) ++pd;
						++pd_den;
						++k;
					}
				}
				const double initial_branch_length = pd/pd_den;
				// Object for the per-branch recombination model: initial branch length set from no-recombination model estimate
				ClonalFrameLaplacePerBranchFunction cff(ctree.node[i],node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,MULTITHREAD,is_imported[i],driving_prior_mean,driving_prior_precision);
				// Force parameterization 1
				cff.parameterization = 1;
				// Evaluate the grid
				vector< mydouble > evalpost(ngrid*ngrid*ngrid*ngrid);
				mydouble denom = 0.0;
				vector<double> param(4);
				ofstream gout;
				if(i==0) gout.open("/Users/wilson/temp/grid.txt");
				char tab = '\t';
				int ctr = 0;
				for(j=0;j<ngrid;j++) {
					param[0] = grid[0][j];
					for(k=0;k<ngrid;k++) {
						param[1] = grid[1][k];
						for(l=0;l<ngrid;l++) {
							param[2] = grid[2][l];
							for(m=0;m<ngrid;m++,ctr++) {
								param[3] = grid[3][m];
								// Evaluate the unnormalized posterior at the grid point
								evalpost[ctr].setlog(-cff.f(param));
								// Update the running total
								denom += evalpost[ctr];
								if(i==0) gout << param[0] << tab << param[1] << tab << param[2] << tab << param[3] << tab << evalpost[ctr].LOG() << endl;
							}
						}
					}
				}
				if(i==0) gout.close();
				// Approximately normalize the posterior density evaluations
				for(ctr=0;ctr<evalpost.size();ctr++) evalpost[ctr] /= denom;
				// Calculate the mean and variance-covariance matrix for the branch				
				ctr = 0;		
				for(j=0;j<ngrid;j++) {
					param[0] = grid[0][j];
					for(k=0;k<ngrid;k++) {
						param[1] = grid[1][k];
						for(l=0;l<ngrid;l++) {
							param[2] = grid[2][l];
							for(m=0;m<ngrid;m++,ctr++) {
								param[3] = grid[3][m];
								// Cycle through the elements of the mean and covariance matrix
								int ii,jj;
								double unlog_evalpost = evalpost[ctr].todouble();
								for(ii=0;ii<nparam;ii++) {
									for(jj=ii+1;jj<nparam;jj++) {
										double tmp = param[ii]*param[jj]*unlog_evalpost;
										Vgrid[i][ii][jj] += tmp;
										Vgrid[i][jj][ii] += tmp;
									}
									Vgrid[i][ii][ii] += param[ii]*param[ii]*unlog_evalpost;
									Egrid[i][ii] += param[ii]*unlog_evalpost;
								}
							}
						}
					}
				}
				// Ensure importation status is updated at the MAP parameter estimate (for now, this is affected by the driving prior)
				//cff.f(Egrid[i]);
				// Update branch length in the tree
				// Note this is unsafe in general because the corresponding node times are not adjusted
				//ctree.node[i].edge_time = final_branch_length;
				// Output results to screen
				cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " R = " << Egrid[i][0] << " I = " << Egrid[i][1] << " D = " << Egrid[i][2] << " M = " << Egrid[i][3] << endl;
			}
			
			// Output the importation status
			//write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			//cout << "Wrote inferred importation status to " << import_out_file << endl;
			// Output the grid approximation
			ofstream lout(laplace_out_file.c_str());
			lout << "Branch";
			for(j=0;j<4;j++) lout << "\t" << "E" << j;
			for(j=0;j<4;j++) {
				for(k=0;k<4;k++) {
					lout << "\t" << "V" << j << k;
				}
			}
			lout << endl;
			for(i=0;i<root_node;i++) {
				lout << ctree_node_labels[i];
				for(j=0;j<4;j++) lout << "\t" << Egrid[i][j];
				for(j=0;j<4;j++) {
					for(k=0;k<4;k++) {
						lout << "\t" << Vgrid[i][j][k];
					}
				}
				lout << endl;
			}
			lout.close();
			cout << "Wrote grid approximation information to " << laplace_out_file << endl;
		} else if(SINGLE_RHO_VITERBI || SINGLE_RHO_FORWARD) {
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) AND recombination parameters under the ClonalFrame model
			// SUBJECT to the constraints that the importation state have frequency < 0.5 AND the recombination divergence exceeds the branch length divergence
			// For computational efficiency, branch lengths are set equal to the number of substitutions implied by the ancestral state reconstruction
			cout << "Beginning parameter estimation" << endl;
			vector< vector<ImportationState> > is_imported(root_node);
			// Constrain the expected number of substitutions per branch according to the ancestral state reconstruction
			vector<double> substitutions_per_branch(root_node,0);
			for(i=0;i<root_node;i++) {
				double pd = 1.0, pd_den = 2.0;
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				int j,k;
				for(j=0,k=0;j<isBLC.size();j++) {
					if(isBLC[j]) {
						Nucleotide dec = node_nuc[dec_id][ipat[k]];
						Nucleotide anc = node_nuc[anc_id][ipat[k]];
						if(dec!=anc) ++pd;
						++pd_den;
						++k;
					}
				}
				substitutions_per_branch[i] = pd/pd_den;
			}
			// Jointly estimate the recombination parameters across branches, subject to the fixed expected number of substitutions per branch
			// Minimum branch length
			const double min_branch_length = global_min_branch_length;
			ClonalFrameSingleRho cff(SINGLE_RHO_VITERBI,ctree,node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,EXCESS_DIVERGENCE_MODEL,MULTITHREAD,is_imported,substitutions_per_branch,min_branch_length,root_node);
			// Setup optimization function
			Powell Pow(cff);
			Pow.coutput = Pow.brent.coutput = SHOW_PROGRESS;
			Pow.TOL = brent_tolerance;
			vector<double> param(3);
			param[0] = log10(initial_rho_over_theta);
			param[1] = log10(initial_mean_import_length);
			param[2] = log10(initial_import_divergence);
			param = Pow.minimize(param,powell_tolerance);
			// Output results to screen
			const double ML = -Pow.function_minimum;
			const double final_rho_over_theta = cff.get_rho_over_theta(param);
			const double final_mean_import_length = cff.get_mean_import_length(param);
			const double final_import_divergence = cff.get_import_divergence(param);
			cout << "Maximum log-likelihood = " << ML << " with parameter estimates: " << endl;
			cout << "rho/theta              = " << final_rho_over_theta << endl;
			cout << "mean import length     = " << final_mean_import_length << endl;
			cout << "mean import divergence = " << final_import_divergence << endl;			
			// Ensure importation status is updated correctly
			for(i=0;i<root_node;i++) {
				const int dec_id = ctree.node[i].id;
				const int anc_id = ctree.node[i].ancestor->id;
				double final_branch_length = cff.get_branch_length(substitutions_per_branch[i],param);
				if(final_branch_length<min_branch_length) final_branch_length = min_branch_length;
				// This part of the algorithm always uses the Viterbi algorithm
				maximum_likelihood_ClonalFrame_branch_allsites(dec_id, anc_id, node_nuc, isBLC, ipat, kappa, empirical_nucleotide_frequencies, final_branch_length, final_rho_over_theta, final_mean_import_length, final_import_divergence, is_imported[i]);
				// Update branch length in the tree - this will set it equal to the number of substitutions implied by the ancestral state reconstruction
				// Note this operation is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = substitutions_per_branch[i];
			}
			// Output the importation status
			write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			cout << "Wrote inferred importation status to " << import_out_file << endl;
		} else {
			// Estimate parameters with branch lengths fixed
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) under the ClonalFrame model
			ClonalFrameParameterFunction cff(ctree,node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,EXCESS_DIVERGENCE_MODEL,MULTITHREAD,root_node);
			vector<double> param(3);
			param[0] = log10(initial_rho_over_theta); param[1] = log10(initial_mean_import_length); param[2] = log10(initial_import_divergence);
			
			// Maximize the likelihood of the parameters
			Powell Pow(cff);
			Pow.coutput = Pow.brent.coutput = SHOW_PROGRESS;
			Pow.TOL = brent_tolerance;
			param = Pow.minimize(param,powell_tolerance);
			// Ensure importation status is updated correctly
			cff.f(param);
			double rho_over_theta = pow(10.,param[0]);
			double mean_import_length = pow(10.,param[1]);
			double import_divergence = pow(10.,param[2]);
			
			cout << "Maximum log-likelihood of " << -Pow.function_minimum << " found in " << cff.neval << " iterations" << endl;
			cout << "rho/theta         = " << rho_over_theta << endl;
			cout << "Import length     = " << mean_import_length << endl;
			cout << "Import divergence = " << import_divergence << endl;
			
			// Now estimate branch lengths
			cout << "Beginning branch optimization:" << endl;
			for(i=0;i<root_node;i++) {
				ClonalFrameBranchLengthFunction cfblf(ctree.node[i],node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,EXCESS_DIVERGENCE_MODEL,rho_over_theta,mean_import_length,import_divergence,MULTITHREAD);
				vector<double> paramBL(1,log10(ctree.node[i].edge_time));
				Powell PowBL(cfblf);
				PowBL.coutput = PowBL.brent.coutput = SHOW_PROGRESS;
				PowBL.TOL = brent_tolerance;
				paramBL = PowBL.minimize(paramBL,powell_tolerance);
				// Update branch length in the tree
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = pow(10.,paramBL[0]);
				cout << "Branch length " << ctree_node_labels[i] << " = " << ctree.node[i].edge_time << endl;
			}
			cout << "Log-likelihood after branch optimization is " << -cff.f(param) << endl;
			
			// Output the importation status
			write_importation_status(cff.is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			cout << "Wrote inferred importation status to " << import_out_file << endl;
		} 
	}
	
	// Output the tree with internal nodes automatically labelled, for cross-referencing with the reconstructed sequences
	write_newick(ctree,ctree_node_labels,tree_out_file.c_str());
	cout << "Wrote processed tree to " << tree_out_file << endl;
	
	cout << "All done in " << (double)(clock()-start_time)/CLOCKS_PER_SEC/60.0 << " minutes." << endl;
	return 0;
}

marginal_tree convert_rooted_NewickTree_to_marginal_tree(NewickTree &newick, vector<string> &tip_labels, vector<string> &all_node_labels) {
	size_t i;
	vector<string> order = tip_labels;
	const int n = tip_labels.size();
	for(i=0;i<n-1;i++) {
		order.push_back("");
	}

	marginal_tree tree;
	// Identify the tips in the NewickTree
	vector<NewickNode*> &allnodes = newick.allnodes;
	size_t nnode = allnodes.size();
	vector<NewickNode*> tips(0);
	vector<NewickNode*> coals(0);
	NewickNode* root = 0;
	for(i=0;i<nnode;i++) {
		// Test for a tip
		if(allnodes[i]->dec.size()==0) {
			tips.push_back(allnodes[i]);
		} else {
			coals.push_back(allnodes[i]);
		}
		// Test for multifurcations
		if(allnodes[i]->dec.size()>2) {
			stringstream errTxt;
			errTxt << "convert_NewickTree_to_marginal_tree(): ";
			errTxt << "Number of descendant nodes (" << allnodes[i]->dec.size();
			errTxt << ") incompatible with a strictly bifurcating rooted tree";
			error(errTxt.str().c_str());
		}
		// Test for the root
		if(allnodes[i]->anc==0) {
			if(root==0) {
				root = allnodes[i];
			} else {
				stringstream errTxt;
				errTxt << "convert_NewickTree_to_marginal_tree(): ";
				errTxt << "Found multiple roots in Newick tree";
				error(errTxt.str().c_str());
			}
		}
	}
	size_t ntips = tips.size();
	// Make sure the number of tips equals that specified by the tip labels
	if(ntips!=tip_labels.size()) {
		stringstream errTxt;
		errTxt << "Number of nodes in Newick tree inconsistent with that expected";
		error(errTxt.str().c_str());
	}
	// Check that the Newick tree is strictly bifurcating, and assume it is correctly rooted
	if(nnode!=2*ntips-1) {
		stringstream errTxt;
		errTxt << "convert_NewickTree_to_marginal_tree(): ";
		errTxt << "Number of nodes (" << nnode << ") and number of tips (" << ntips;
		errTxt << ") incompatible with a strictly bifurcating rooted tree";
		error(errTxt.str().c_str());
	}
	// Calculate node times. Ensure all branches have non-zero length
	const double minbranchlength = 1e-12;
	vector<NewickNode*> root2tip(1,root);	// temporary ordering of nodes from root to tips
	vector<double> ageroot2tip(1,0.0);		// corresponding age of each node in root2tip
	double youngest_node = 0.0;
	size_t iroot2tip;
	for(iroot2tip=0;iroot2tip<nnode;iroot2tip++) {
		if(iroot2tip>=root2tip.size()) {
			stringstream errTxt;
			errTxt << "convert_NewickTree_to_marginal_tree(): ";
			errTxt << "iroot2tip exceeded size of root2tip";
			error(errTxt.str().c_str());
		}
		// Add descendants of current node to list and calculate node times
		// counting with age increasing backwards in time, but the root node at time 0
		int idec;
		for(idec=0;idec<root2tip[iroot2tip]->dec.size();idec++) {
			root2tip.push_back(root2tip[iroot2tip]->dec[idec]);
			// Ensure the descendant is always younger than its ancestor
			double branchlength = root2tip[root2tip.size()-1]->len;
			if(branchlength<minbranchlength) branchlength = minbranchlength;
			ageroot2tip.push_back(ageroot2tip[iroot2tip]-branchlength);
			if(ageroot2tip[ageroot2tip.size()-1]<youngest_node) {
				youngest_node = ageroot2tip[ageroot2tip.size()-1];
			}
		}
	}
	// Give the nodes ids for use in the marginal tree, using the convention
	// that the tips comprise the first ntips nodes, and the remainder are in height order.
	// The order of the tips must be the same for all marginal trees for cross-referencing.
	// To achieve this, set the ordering on the first call, when order==vector<string>(0)
	// and then impose this ordering thereafter. Note that the time-ordering of the tips is
	// unimportant.
	
	vector<size_t> ixroot2tip(0);
	for(iroot2tip=0;iroot2tip<nnode;iroot2tip++) ixroot2tip.push_back(iroot2tip);
	if(order.size()==0) {
		stringstream errTxt;
		errTxt << "convert_NewickTree_to_marginal_tree(): ";
		errTxt << "Node labels in Newick file incompatible with expected tip labels";
		error(errTxt.str().c_str());
	} else {
		// Calculate the position (first instance) in "order" of the label of each node in root2tip
		vector<size_t> labelorder;
		for(iroot2tip=0;iroot2tip<nnode;iroot2tip++) {
			vector<string>::iterator _find = std::find(order.begin(),order.end(),root2tip[iroot2tip]->str);
			if(_find==order.end()) {
				stringstream errTxt;
				errTxt << "convert_NewickTree_to_marginal_tree(): ";
				errTxt << "Newick tree tip label " << root2tip[iroot2tip]->str << " was not expected";
				error(errTxt.str().c_str());
			}
			labelorder.push_back(_find-order.begin());
		}
		// Re-order root2tip and ageroot2tip by (1) label (tips only) (2) age (coalescences only)
		std::stable_sort(ixroot2tip.begin(),ixroot2tip.end(),orderNewickNodesByStatusLabelAndAge(root2tip,ageroot2tip,labelorder));
	}
	// Assign each node in root2tip an index by calculating the rank of each element in root2tip in ixroot2tip
	map<const NewickNode*,int> nodeIndex;
	for(iroot2tip=0;iroot2tip<nnode;iroot2tip++) nodeIndex[root2tip[ixroot2tip[iroot2tip]]] = (int)iroot2tip;
	
	// Construct the internal representation of the tree
	tree.initialize(0,(int)ntips);
	// Add the tips and then the coalescence events via root2tip in the order stored in ixroot2tip
	bool internal_nodes_begun = false;
	// Record the tip names
	all_node_labels = vector<string>(0);
	for(iroot2tip=0;iroot2tip<nnode;iroot2tip++) {
		size_t ix = ixroot2tip[iroot2tip];
		const NewickNode *node = root2tip[ix];
		// Sanity check
		if(nodeIndex[node]!=iroot2tip) {
			stringstream errTxt;
			errTxt << "convert_NewickTree_to_marginal_tree(): ";
			errTxt << "Inconsistency in internal node numbering";
			error(errTxt.str().c_str());
		}
		if(node->dec.size()==0) {
			// If tip
			if(internal_nodes_begun) {
				stringstream errTxt;
				errTxt << "convert_NewickTree_to_marginal_tree(): ";
				errTxt << "internal nodes added to marginal tree before all tips";
				error(errTxt.str().c_str());
			}
			double age = ageroot2tip[ix]-youngest_node;
			if(fabs(age)<1e-6) age = 0.0;
			tree.add_base_node(&age,nodeIndex[node]);
		} else if(node->dec.size()==2) {
			// If internal node
			internal_nodes_begun = true;
			double age = ageroot2tip[ix]-youngest_node;
			if(fabs(age)<1e-6) age = 0.0;
			tree.coalesce(age,nodeIndex[node->dec[0]],nodeIndex[node->dec[1]]);
		} else {
			stringstream errTxt;
			errTxt << "convert_NewickTree_to_marginal_tree(): ";
			errTxt << "only tips or bifurcating nodes expected. " << node->dec.size() << " descendants not allowed.";
			error(errTxt.str().c_str());
		}
		if(node->str!="") {
			all_node_labels.push_back(node->str);
		} else {
			stringstream autolab;
			autolab << "NODE_" << iroot2tip+1;
			all_node_labels.push_back(autolab.str());
		}
	}
	return tree;
}

marginal_tree convert_unrooted_NewickTree_to_marginal_tree(NewickTree &newick, vector<string> &tip_labels, vector<string> &all_node_labels) {
	size_t i;
	vector<string> order = tip_labels;
	const int n = tip_labels.size();
	for(i=0;i<n-1;i++) {
		order.push_back("");
	}
	
	marginal_tree tree;
	// Identify the tips in the NewickTree
	vector<NewickNode*> &allnodes = newick.allnodes;
	size_t nnode = allnodes.size();
	vector<NewickNode*> tips(0);
	vector<NewickNode*> coals(0);
	NewickNode* root = 0;
	for(i=0;i<nnode;i++) {
		// Test for a tip
		if(allnodes[i]->dec.size()==0) {
			tips.push_back(allnodes[i]);
		} else {
			coals.push_back(allnodes[i]);
		}
		// Test for multifurcations
		if(allnodes[i]->anc!=0 && allnodes[i]->dec.size()==3) {
			stringstream errTxt;
			errTxt << "convert_NewickTree_to_marginal_tree(): ";
			errTxt << "Only the root is allowed 3 descendant nodes";
			error(errTxt.str().c_str());
		}
		if(allnodes[i]->dec.size()>3) {
			stringstream errTxt;
			errTxt << "convert_NewickTree_to_marginal_tree(): ";
			errTxt << "Number of descendant nodes (" << allnodes[i]->dec.size();
			errTxt << ") incompatible with a bifurcating unrooted tree";
			error(errTxt.str().c_str());
		}
		// Test for the root
		if(allnodes[i]->anc==0) {
			if(root==0) {
				root = allnodes[i];
				if(root->dec.size()!=3) {
					stringstream errTxt;
					errTxt << "convert_NewickTree_to_marginal_tree(): ";
					errTxt << "Deepest node in unrooted Newick tree expected to have 3 descendants";
					error(errTxt.str().c_str());
				}
			} else {
				stringstream errTxt;
				errTxt << "convert_NewickTree_to_marginal_tree(): ";
				errTxt << "Found multiple roots in Newick tree";
				error(errTxt.str().c_str());
			}
		}
	}
	size_t ntips = tips.size();
	// Make sure the number of tips equals that specified by the tip labels
	if(ntips!=tip_labels.size()) {
		stringstream errTxt;
		errTxt << "Number of nodes in Newick tree inconsistent with that expected";
		error(errTxt.str().c_str());
	}
	// Check that the Newick tree is consistent with an unrooted strictly bifurcating tree
	if(nnode!=2*ntips-2) {
		stringstream errTxt;
		errTxt << "convert_NewickTree_to_marginal_tree(): ";
		errTxt << "Number of nodes (" << nnode << ") and number of tips (" << ntips;
		errTxt << ") incompatible with an unrooted bifurcating tree";
		error(errTxt.str().c_str());
	}
	// Calculate node times. Ensure all branches have non-zero length
	const double minbranchlength = 1e-12;
	vector<NewickNode*> root2tip(1,root);	// temporary ordering of nodes from root to tips
	vector<double> ageroot2tip(1,0.0);		// corresponding age of each node in root2tip
	double youngest_node = 0.0;
	size_t iroot2tip;
	for(iroot2tip=0;iroot2tip<nnode;iroot2tip++) {
		if(iroot2tip>=root2tip.size()) {
			stringstream errTxt;
			errTxt << "convert_NewickTree_to_marginal_tree(): ";
			errTxt << "iroot2tip exceeded size of root2tip";
			error(errTxt.str().c_str());
		}
		// Add descendants of current node to list and calculate node times
		// counting with age increasing backwards in time, but the root node at time 0
		int idec;
		for(idec=0;idec<root2tip[iroot2tip]->dec.size();idec++) {
			root2tip.push_back(root2tip[iroot2tip]->dec[idec]);
			// Ensure the descendant is always younger than its ancestor
			double branchlength = root2tip[root2tip.size()-1]->len;
			if(branchlength<minbranchlength) branchlength = minbranchlength;
			ageroot2tip.push_back(ageroot2tip[iroot2tip]-branchlength);
			if(ageroot2tip[ageroot2tip.size()-1]<youngest_node) {
				youngest_node = ageroot2tip[ageroot2tip.size()-1];
			}
		}
	}
	// Give the nodes ids for use in the marginal tree, using the convention
	// that the tips comprise the first ntips nodes, and the remainder are in height order.
	// The order of the tips must be the same for all marginal trees for cross-referencing.
	// To achieve this, set the ordering on the first call, when order==vector<string>(0)
	// and then impose this ordering thereafter. Note that the time-ordering of the tips is
	// unimportant.
	
	vector<size_t> ixroot2tip(0);
	for(iroot2tip=0;iroot2tip<nnode;iroot2tip++) ixroot2tip.push_back(iroot2tip);
	if(order.size()==0) {
		stringstream errTxt;
		errTxt << "convert_NewickTree_to_marginal_tree(): ";
		errTxt << "Node labels in Newick file incompatible with expected tip labels";
		error(errTxt.str().c_str());
	} else {
		// Calculate the position (first instance) in "order" of the label of each node in root2tip
		vector<size_t> labelorder;
		for(iroot2tip=0;iroot2tip<nnode;iroot2tip++) {
			vector<string>::iterator _find = std::find(order.begin(),order.end(),root2tip[iroot2tip]->str);
			if(_find==order.end()) {
				stringstream errTxt;
				errTxt << "convert_NewickTree_to_marginal_tree(): ";
				errTxt << "Newick tree tip label " << root2tip[iroot2tip]->str << " was not expected";
				error(errTxt.str().c_str());
			}
			labelorder.push_back(_find-order.begin());
		}
		// Re-order root2tip and ageroot2tip by (1) label (tips only) (2) age (coalescences only)
		std::stable_sort(ixroot2tip.begin(),ixroot2tip.end(),orderNewickNodesByStatusLabelAndAge(root2tip,ageroot2tip,labelorder));
	}
	// Assign each node in root2tip an index by calculating the rank of each element in root2tip in ixroot2tip
	map<const NewickNode*,int> nodeIndex;
	for(iroot2tip=0;iroot2tip<nnode;iroot2tip++) nodeIndex[root2tip[ixroot2tip[iroot2tip]]] = (int)iroot2tip;
	
	// Construct the internal representation of the tree
	tree.initialize(0,(int)ntips);
	// Add the tips and then the coalescence events via root2tip in the order stored in ixroot2tip
	bool internal_nodes_begun = false;
	// Record the tip names
	all_node_labels = vector<string>(0);
	for(iroot2tip=0;iroot2tip<nnode-1;iroot2tip++) {
		size_t ix = ixroot2tip[iroot2tip];
		const NewickNode *node = root2tip[ix];
		// Sanity check
		if(nodeIndex[node]!=iroot2tip) {
			stringstream errTxt;
			errTxt << "convert_NewickTree_to_marginal_tree(): ";
			errTxt << "Inconsistency in internal node numbering";
			error(errTxt.str().c_str());
		}
		if(node->dec.size()==0) {
			// If tip
			if(internal_nodes_begun) {
				stringstream errTxt;
				errTxt << "convert_NewickTree_to_marginal_tree(): ";
				errTxt << "internal nodes added to marginal tree before all tips";
				error(errTxt.str().c_str());
			}
			double age = ageroot2tip[ix]-youngest_node;
			if(fabs(age)<1e-6) age = 0.0;
			tree.add_base_node(&age,nodeIndex[node]);
		} else if(node->dec.size()==2) {
			// If internal node
			internal_nodes_begun = true;
			double age = ageroot2tip[ix]-youngest_node;
			if(fabs(age)<1e-6) age = 0.0;
			tree.coalesce(age,nodeIndex[node->dec[0]],nodeIndex[node->dec[1]]);
		} else {
			stringstream errTxt;
			errTxt << "convert_NewickTree_to_marginal_tree(): ";
			errTxt << "only tips or bifurcating nodes expected. " << node->dec.size() << " descendants not allowed.";
			error(errTxt.str().c_str());
		}
		if(node->str!="") {
			all_node_labels.push_back(node->str);
		} else {
			stringstream autolab;
			autolab << "NODE_" << iroot2tip+1;
			all_node_labels.push_back(autolab.str());
		}
	}
	// Deal with the root separately
	iroot2tip = nnode-1;
	size_t ix = ixroot2tip[iroot2tip];
	const NewickNode *node = root2tip[ix];
	// Sanity check
	if(nodeIndex[node]!=iroot2tip) {
		stringstream errTxt;
		errTxt << "convert_NewickTree_to_marginal_tree(): ";
		errTxt << "Inconsistency in internal node numbering";
		error(errTxt.str().c_str());
	}
	if(node->dec.size()!=3) {
		stringstream errTxt;
		errTxt << "convert_NewickTree_to_marginal_tree(): ";
		errTxt << "expected trifurcating root node";
		error(errTxt.str().c_str());
	}
	if(node->anc!=NULL) {
		stringstream errTxt;
		errTxt << "convert_NewickTree_to_marginal_tree(): ";
		errTxt << "expected orphan root node";
		error(errTxt.str().c_str());
	}
	double age = ageroot2tip[ix]-youngest_node;
	if(fabs(age)<1e-6) age = 0.0;
	// Coalesce the first two descendants
	tree.coalesce(age,nodeIndex[node->dec[0]],nodeIndex[node->dec[1]]);
	if(node->str!="") {
		all_node_labels.push_back(node->str);
	} else {
		stringstream autolab;
		autolab << "NODE_" << iroot2tip+1;
		all_node_labels.push_back(autolab.str());
	}
	// Coalesce the resulting node with the third descendant to make the absolute root (this branch has exactly zero length)
	int penultimate_nodeid = nnode-1;
	tree.coalesce(age,penultimate_nodeid,nodeIndex[node->dec[2]]);
	stringstream autolab;
	autolab << "NODE_" << iroot2tip+2;
	all_node_labels.push_back(autolab.str());
	
	return tree;
}


vector<int> compute_compatibility(DNA &fa, marginal_tree &ctree, vector<bool> &anyN, bool purge_singletons) {
	// Sample size
	const int n = fa.nseq;
	// Sequence length
	const int L = fa.lseq;
	
	// Results of initial incompatibility test: -1 (invariant or singleton, compatible), 0 (2 alleles, not tested), 2 (>2 alleles, incompatible)
	vector<int> iscompat(L,0);
	anyN = vector<bool>(L,false);
	// Convert FASTA file to binary: if more than two alleles mark as incompatible: -1 (uninitialized), 0 (reference allele), 1 (first non-reference allele), 2 (second non-reference allele)
	// Let -2 be a no-call (N)
	Matrix<int> bip(n,L,-1);
	int i,pos;
	for(pos=0;pos<L;pos++) {		
		char* allele0 = NULL;
		int nallele0 = 0;
		char* allele1 = NULL;
		int nallele1 = 0;
		for(i=0;i<n;i++) {
			if(fa[i][pos]!='N' && fa[i][pos]!='-' && fa[i][pos]!='X' && fa[i][pos]!='?') {
				// If not an N
				if(allele0==NULL) {
					allele0 = new char(fa[i][pos]);
					bip[i][pos] = 0;
					++nallele0;
				} else {
					if(fa[i][pos]==*allele0) {
						bip[i][pos] = 0;
						++nallele0;
					} else {
						if(allele1==NULL) {
							allele1 = new char(fa[i][pos]);
							bip[i][pos] = 1;
							++nallele1;
						} else {
							if(fa[i][pos]==*allele1) {
								bip[i][pos] = 1;
								++nallele1;
							} else {
								bip[i][pos] = 2;
								iscompat[pos] = 2;
								break;
							}
						}
					}
				}
			} else {
				// If an N
				anyN[pos] = true;
				bip[i][pos] = -2;
			}
		}
		if(iscompat[pos]==0) {
			if(allele0==NULL || allele1==NULL) {
				// Invariant site (or all Ns): must be compatible
				iscompat[pos] = -1;
			} else if(nallele0==1 || nallele1==1) {
				// Singleton: must be compatible
				iscompat[pos] = (purge_singletons) ? -1 : -2;
			}
		}
		if(allele0!=NULL) delete allele0;
		if(allele1!=NULL) delete allele1;
	}
	// Create a bip file for the internal branches of the tree, of which there will be n-2
	Matrix<int> treebip(n,n-2,-1);
	
	// Add "mutations" encoding the branches of the clonal frame
	// The first index is for the sequence (including internal sequences) and the second is for the branch encoded (equivalent to the site)
	Matrix<int> cstate(2*n-1,2*n-1,-1);
	int j,k;
	// Assign 0 to the root node for every site
	for(k=0;k<2*n-1;k++) cstate[2*n-2][k] = 0;
	// Work from root to tips inheriting the state or, if the focal branch, introducing the mutated state
	for(j=2*n-3;j>=0;j--) {
		for(k=0;k<2*n-1;k++) {
			if(j==k) {
				cstate[j][k] = 1;
			} else {
				const mt_node *node = &(ctree.node[j]);
				const mt_node *parent = node->ancestor;
				const int parentState = cstate[parent->id][k];
				cstate[j][k] = parentState;
			}
		}
	}
	
	// Determine compatibility with the clonal frame
	// Test whether the observed partitions in the FASTA file are incompatible with any branches in the Newick tree
	// by tracking whether each of the four possible "haplotypes" has been observed. 
	// pos is the position in the FASTA file, j is the individual in the FASTA file and k is the branch in the Newick tree
	for(pos=0;pos<L;pos++) {
		// First index is for the branches in the tree. Second and third indices are for the four possible
		// "haplotypes" (00, 01, 10 and 11).
		vector< Matrix<bool> > hap(2*n-1, Matrix<bool>(2,2,false));
		if(iscompat[pos]==0) {
			for(j=0;j<n;j++) {
				const int jallele = bip[j][pos];
				if(jallele!=-2) {
					for(k=0;k<2*n-1;k++) {
						const int kallele = cstate[j][k];
						if(kallele!=-2) {
							hap[k][jallele][kallele] = true;
						}
					}
				}
			}
			bool allcompat = true;
			for(k=0;k<2*n-1;k++) {
				if(hap[k][0][0] && hap[k][0][1] && hap[k][1][0] && hap[k][1][1]) {
					allcompat = false;
					break;
				}
			}
			if(!allcompat) {
				iscompat[pos] = 1;
			}
		} else if(iscompat[pos]==-2) {
			iscompat[pos] = 0;
		}
	}
	
	return iscompat;
}

NewickTree read_Newick(const char* newick_file) {
	std::ifstream fnewick(newick_file);
	if(!fnewick.is_open()) {
		stringstream errTxt;
		errTxt << "Could not open " << newick_file;
		error(errTxt.str().c_str());
	}
	string snewick;
	std::getline(fnewick,snewick);	
	fnewick.close();
	return NewickTree(snewick);
}

Matrix<Nucleotide> FASTA_to_nucleotide(DNA &fa, vector<double> &empirical_nucleotide_frequencies, vector<bool> usesite) {
	int i,j,k;
	int nsites = 0;
	for(j=0;j<usesite.size();j++) {
		if(usesite[j]) ++nsites;
	}
	Matrix<Nucleotide> nuc(fa.nseq,nsites,N_ambiguous);
	empirical_nucleotide_frequencies = vector<double>(4,0.0);
	double total_empirical_count = 0.0;
	for(j=0,k=0;j<fa.lseq;j++) {
		if(usesite[j]) {
			for(i=0;i<fa.nseq;i++) {
				switch(toupper(fa[i][j])) {
					case 'A':
						nuc[i][k] = Adenine;
						++empirical_nucleotide_frequencies[Adenine];
						++total_empirical_count;
						break;
					case 'G':
						nuc[i][k] = Guanine;
						++empirical_nucleotide_frequencies[Guanine];
						++total_empirical_count;
						break;					
					case 'C':
						nuc[i][k] = Cytosine;
						++empirical_nucleotide_frequencies[Cytosine];
						++total_empirical_count;
						break;					
					case 'T': case 'U':
						nuc[i][k] = Thymine;
						++empirical_nucleotide_frequencies[Thymine];
						++total_empirical_count;
						break;
					case 'N': case 'X': case '-': case '?':
						nuc[i][k] = N_ambiguous;
						break;
					default:
						stringstream errTxt;
						errTxt << "FASTA_to_nucleotide(): unsupported base " << fa[i][j] << " in sequence " << i;
						errTxt << " ("<< fa.label[i] << ") position " << j;
						error(errTxt.str().c_str());
				}
			}
			++k;
		} else {
			for(i=0;i<fa.nseq;i++) {
				switch(toupper(fa[i][j])) {
					case 'A':
						++empirical_nucleotide_frequencies[Adenine];
						++total_empirical_count;
						break;
					case 'G':
						++empirical_nucleotide_frequencies[Guanine];
						++total_empirical_count;
						break;					
					case 'C':
						++empirical_nucleotide_frequencies[Cytosine];
						++total_empirical_count;
						break;					
					case 'T': case 'U':
						++empirical_nucleotide_frequencies[Thymine];
						++total_empirical_count;
						break;
					case 'N': case 'X': case '-': case '?':
						break;
					default:
						stringstream errTxt;
						errTxt << "FASTA_to_nucleotide(): unsupported base " << fa[i][j] << " in sequence " << i;
						errTxt << " ("<< fa.label[i] << ") position " << j;
						error(errTxt.str().c_str());
				}
			}
		}
	}
	for(i=0;i<4;i++) empirical_nucleotide_frequencies[i] /= total_empirical_count;
	return nuc;
}

void find_alignment_patterns(Matrix<Nucleotide> &nuc, vector<bool> &iscompat, vector<string> &pat, vector<int> &pat1, vector<int> &cpat, vector<int> &ipat) {
	pat = vector<string>(0);
	pat1 = vector<int>(0);
	cpat = vector<int>(0);
	ipat = vector<int>(nuc.ncols());
	static const char AGCTN[5] = {'A','G','C','T','N'};
	int i,j,pos;
	for(pos=0;pos<nuc.ncols();pos++) {
		if(iscompat[pos]) {
			string pospat = "";
			for(i=0;i<nuc.nrows();i++) {
				pospat += AGCTN[nuc[i][pos]];
			}
			for(j=0;j<pat.size();j++) {
				if(pospat==pat[j]) break;
			}
			if(j==pat.size()) {
				pat.push_back(pospat);
				pat1.push_back(pos);
				cpat.push_back(1);
				ipat[pos] = j;
			} else {
				++cpat[j];
				ipat[pos] = j;
			}
		} else {
			// If not a compatible site
			ipat[pos] = -1;
		}
	}
}

/* Use the following in Maple to generate this code (k is 1/transition:transversion ratio, i.e. k=1/kappa)

M := Matrix([
[-g-k*(c+t),g,k*c,k*t],
[a,-a-k*(c+t),k*c,k*t],
[k*a,k*g,-k*(a+g)-t,t],
[k*a,k*g,c,-k*(a+g)-c]]);
R:= simplify(-a*M[1,1]-g*M[2,2]-c*M[3,3]-t*M[4,4]);
M:=simplify(M/R);
 
P:=simplify(LinearAlgebra:-MatrixExponential(M,x));
CodeGeneration:-C(subs([a=pi[1],g=pi[2],c=pi[3],t=pi[4]],P),optimize,resultname="ptrans[i]");

 */ 
vector< Matrix<double> > compute_HKY85_ptrans(const marginal_tree &ctree, const double kappa, const vector<double> &pi) {
	const double k = 1.0/kappa;
	const int nnodes = ctree.size;
	Matrix<double> ptrans_element(4,4,0.0);
	vector< Matrix<double> > ptrans(nnodes,ptrans_element);
	int i;
	for(i=0;i<nnodes-1;i++) {
		const double x = ctree.node[i].edge_time;
		double t1 = pi[2] + pi[3];
		double t2 = t1 * pi[0];
		double t3 = t1 * pi[1];
		double t4 = pi[0] * pi[1] + pi[2] * pi[3] + (t2 + t3) * k;
		t4 = 0.1e1 / t4;
		double t5 = -0.1e1 / 0.2e1;
		double t6 = exp(t5 * (t1 * k + pi[0] + pi[1]) * x * t4);
		double t7 = pi[2] + pi[3] + pi[0] + pi[1];
		double t8 = exp(t5 * k * t7 * x * t4);
		double t9 = pow(pi[1], 0.2e1);
		double t10 = pow(pi[0], 0.2e1);
		double t11 = pi[0] + pi[1];
		double t12 = t7 * t6 - t1 * t8 - pi[0] - pi[1];
		double t13 = t8 - 0.1e1;
		double t14 = 0.1e1 / t11;
		double t15 = 0.1e1 / t7;
		double t16 = t13 * pi[2] * t15;
		double t17 = t13 * pi[3] * t15;
		t4 = exp(t5 * (t11 * k + pi[2] + pi[3]) * x * t4);
		t5 = pow(pi[3], 0.2e1);
		double t18 = pow(pi[2], 0.2e1);
		t11 = t11 * t8;
		t7 = -t11 + t7 * t4 - pi[3] - pi[2];
		t1 = 0.1e1 / t1;
		double t19 = t13 * pi[0] * t15;
		t13 = t13 * pi[1] * t15;
		ptrans[i][0][0] = (t6 * t9 + ((pi[0] + pi[3] + pi[2]) * t6 + pi[0]) * pi[1] + t2 * t8 + t10) * t14 * t15;
		ptrans[i][0][1] = -pi[1] * t12 * t14 * t15;
		ptrans[i][0][2] = -t16;
		ptrans[i][0][3] = -t17;
		ptrans[i][1][0] = -pi[0] * t12 * t14 * t15;
		ptrans[i][1][1] = (t6 * t10 + ((pi[2] + pi[1] + pi[3]) * t6 + pi[1]) * pi[0] + t3 * t8 + t9) * t14 * t15;
		ptrans[i][1][2] = -t16;
		ptrans[i][1][3] = -t17;
		ptrans[i][2][0] = -t19;
		ptrans[i][2][1] = -t13;
		ptrans[i][2][2] = (t4 * t5 + ((pi[0] + pi[2] + pi[1]) * t4 + pi[2]) * pi[3] + t11 * pi[2] + t18) * t1 * t15;
		ptrans[i][2][3] = -t7 * pi[3] * t1 * t15;
		ptrans[i][3][0] = -t19;
		ptrans[i][3][1] = -t13;
		ptrans[i][3][2] = -t7 * pi[2] * t1 * t15;
		ptrans[i][3][3] = (t4 * t18 + ((pi[0] + pi[1] + pi[3]) * t4 + pi[3]) * pi[2] + t11 * pi[3] + t5) * t1 * t15;		
	}
	// For the root node, whose edge length is undefined/infinite
	ptrans[nnodes-1][0][0] = pi[0];
	ptrans[nnodes-1][1][0] = pi[0];
	ptrans[nnodes-1][2][0] = pi[0];
	ptrans[nnodes-1][3][0] = pi[0];
	ptrans[nnodes-1][0][1] = pi[1];
	ptrans[nnodes-1][1][1] = pi[1];
	ptrans[nnodes-1][2][1] = pi[1];
	ptrans[nnodes-1][3][1] = pi[1];
	ptrans[nnodes-1][0][2] = pi[2];
	ptrans[nnodes-1][1][2] = pi[2];
	ptrans[nnodes-1][2][2] = pi[2];
	ptrans[nnodes-1][3][2] = pi[2];
	ptrans[nnodes-1][0][3] = pi[3];
	ptrans[nnodes-1][1][3] = pi[3];
	ptrans[nnodes-1][2][3] = pi[3];
	ptrans[nnodes-1][3][3] = pi[3];
	// Prevent any underflow problems
	int j,l;
	for(i=0;i<nnodes;i++) {
		for(j=0;j<4;j++) {
			for(l=0;l<4;l++) {
				if(ptrans[i][j][l]>1.0) {
					ptrans[i][j][l] = 1.0;
				} else if(ptrans[i][j][l]<1.0e-100) {
					ptrans[i][j][l] = 1.0e-100;
				}
			}
		}
	}
	return ptrans;
}

vector< Matrix<double> > compute_HKY85_ptrans(const marginal_tree &ctree, const double kappa, const vector<double> &pi, const double import_divergence, const bool excess_divergence_model) {
	const double k = 1.0/kappa;
	const int nnodes = ctree.size;
	Matrix<double> ptrans_element(4,4,0.0);
	vector< Matrix<double> > ptrans(nnodes,ptrans_element);
	int i;
	for(i=0;i<nnodes-1;i++) {
		double x = import_divergence;
		if(excess_divergence_model) x += ctree.node[i].edge_time;
		double t1 = pi[2] + pi[3];
		double t2 = t1 * pi[0];
		double t3 = t1 * pi[1];
		double t4 = pi[0] * pi[1] + pi[2] * pi[3] + (t2 + t3) * k;
		t4 = 0.1e1 / t4;
		double t5 = -0.1e1 / 0.2e1;
		double t6 = exp(t5 * (t1 * k + pi[0] + pi[1]) * x * t4);
		double t7 = pi[2] + pi[3] + pi[0] + pi[1];
		double t8 = exp(t5 * k * t7 * x * t4);
		double t9 = pow(pi[1], 0.2e1);
		double t10 = pow(pi[0], 0.2e1);
		double t11 = pi[0] + pi[1];
		double t12 = t7 * t6 - t1 * t8 - pi[0] - pi[1];
		double t13 = t8 - 0.1e1;
		double t14 = 0.1e1 / t11;
		double t15 = 0.1e1 / t7;
		double t16 = t13 * pi[2] * t15;
		double t17 = t13 * pi[3] * t15;
		t4 = exp(t5 * (t11 * k + pi[2] + pi[3]) * x * t4);
		t5 = pow(pi[3], 0.2e1);
		double t18 = pow(pi[2], 0.2e1);
		t11 = t11 * t8;
		t7 = -t11 + t7 * t4 - pi[3] - pi[2];
		t1 = 0.1e1 / t1;
		double t19 = t13 * pi[0] * t15;
		t13 = t13 * pi[1] * t15;
		ptrans[i][0][0] = (t6 * t9 + ((pi[0] + pi[3] + pi[2]) * t6 + pi[0]) * pi[1] + t2 * t8 + t10) * t14 * t15;
		ptrans[i][0][1] = -pi[1] * t12 * t14 * t15;
		ptrans[i][0][2] = -t16;
		ptrans[i][0][3] = -t17;
		ptrans[i][1][0] = -pi[0] * t12 * t14 * t15;
		ptrans[i][1][1] = (t6 * t10 + ((pi[2] + pi[1] + pi[3]) * t6 + pi[1]) * pi[0] + t3 * t8 + t9) * t14 * t15;
		ptrans[i][1][2] = -t16;
		ptrans[i][1][3] = -t17;
		ptrans[i][2][0] = -t19;
		ptrans[i][2][1] = -t13;
		ptrans[i][2][2] = (t4 * t5 + ((pi[0] + pi[2] + pi[1]) * t4 + pi[2]) * pi[3] + t11 * pi[2] + t18) * t1 * t15;
		ptrans[i][2][3] = -t7 * pi[3] * t1 * t15;
		ptrans[i][3][0] = -t19;
		ptrans[i][3][1] = -t13;
		ptrans[i][3][2] = -t7 * pi[2] * t1 * t15;
		ptrans[i][3][3] = (t4 * t18 + ((pi[0] + pi[1] + pi[3]) * t4 + pi[3]) * pi[2] + t11 * pi[3] + t5) * t1 * t15;		
	}
	// For the root node, whose edge length is undefined/infinite
	ptrans[nnodes-1][0][0] = pi[0];
	ptrans[nnodes-1][1][0] = pi[0];
	ptrans[nnodes-1][2][0] = pi[0];
	ptrans[nnodes-1][3][0] = pi[0];
	ptrans[nnodes-1][0][1] = pi[1];
	ptrans[nnodes-1][1][1] = pi[1];
	ptrans[nnodes-1][2][1] = pi[1];
	ptrans[nnodes-1][3][1] = pi[1];
	ptrans[nnodes-1][0][2] = pi[2];
	ptrans[nnodes-1][1][2] = pi[2];
	ptrans[nnodes-1][2][2] = pi[2];
	ptrans[nnodes-1][3][2] = pi[2];
	ptrans[nnodes-1][0][3] = pi[3];
	ptrans[nnodes-1][1][3] = pi[3];
	ptrans[nnodes-1][2][3] = pi[3];
	ptrans[nnodes-1][3][3] = pi[3];
	// Prevent any underflow problems
	int j,l;
	for(i=0;i<nnodes;i++) {
		for(j=0;j<4;j++) {
			for(l=0;l<4;l++) {
				if(ptrans[i][j][l]>1.0) {
					ptrans[i][j][l] = 1.0;
				} else if(ptrans[i][j][l]<1.0e-100) {
					ptrans[i][j][l] = 1.0e-100;
				}
			}
		}
	}
	return ptrans;
}

Matrix<mydouble> compute_HKY85_ptrans(const double x, const double kappa, const vector<double> &pi) {
	const double k = 1.0/kappa;
	Matrix<mydouble> ptrans(4,4,0.0);
	double t1 = pi[2] + pi[3];
	double t2 = t1 * pi[0];
	double t3 = t1 * pi[1];
	double t4 = pi[0] * pi[1] + pi[2] * pi[3] + (t2 + t3) * k;
	t4 = 0.1e1 / t4;
	double t5 = -0.1e1 / 0.2e1;
	double t6 = exp(t5 * (t1 * k + pi[0] + pi[1]) * x * t4);
	double t7 = pi[2] + pi[3] + pi[0] + pi[1];
	double t8 = exp(t5 * k * t7 * x * t4);
	double t9 = pow(pi[1], 0.2e1);
	double t10 = pow(pi[0], 0.2e1);
	double t11 = pi[0] + pi[1];
	double t12 = t7 * t6 - t1 * t8 - pi[0] - pi[1];
	double t13 = t8 - 0.1e1;
	double t14 = 0.1e1 / t11;
	double t15 = 0.1e1 / t7;
	double t16 = t13 * pi[2] * t15;
	double t17 = t13 * pi[3] * t15;
	t4 = exp(t5 * (t11 * k + pi[2] + pi[3]) * x * t4);
	t5 = pow(pi[3], 0.2e1);
	double t18 = pow(pi[2], 0.2e1);
	t11 = t11 * t8;
	t7 = -t11 + t7 * t4 - pi[3] - pi[2];
	t1 = 0.1e1 / t1;
	double t19 = t13 * pi[0] * t15;
	t13 = t13 * pi[1] * t15;
	double temp;
	temp = (t6 * t9 + ((pi[0] + pi[3] + pi[2]) * t6 + pi[0]) * pi[1] + t2 * t8 + t10) * t14 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[0][0] = temp;
	temp = -pi[1] * t12 * t14 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[0][1] = temp;
	temp = -t16;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[0][2] = temp;
	temp = -t17;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[0][3] = temp;
	temp = -pi[0] * t12 * t14 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[1][0] = temp;
	temp = (t6 * t10 + ((pi[2] + pi[1] + pi[3]) * t6 + pi[1]) * pi[0] + t3 * t8 + t9) * t14 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[1][1] = temp;
	temp = -t16;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[1][2] = temp;
	temp = -t17;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[1][3] = temp;
	temp = -t19;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[2][0] = temp;
	temp = -t13;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[2][1] = temp;
	temp = (t4 * t5 + ((pi[0] + pi[2] + pi[1]) * t4 + pi[2]) * pi[3] + t11 * pi[2] + t18) * t1 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[2][2] = temp;
	temp = -t7 * pi[3] * t1 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[2][3] = temp;
	temp = -t19;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[3][0] = temp;
	temp = -t13;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[3][1] = temp;
	temp = -t7 * pi[2] * t1 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[3][2] = temp;
	temp = (t4 * t18 + ((pi[0] + pi[1] + pi[3]) * t4 + pi[3]) * pi[2] + t11 * pi[3] + t5) * t1 * t15;		
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[3][3] = temp;
	return ptrans;
}

Matrix<double> dcompute_HKY85_ptrans(const double x, const double kappa, const vector<double> &pi) {
	const double k = 1.0/kappa;
	Matrix<double> ptrans(4,4,0.0);
	double t1 = pi[2] + pi[3];
	double t2 = t1 * pi[0];
	double t3 = t1 * pi[1];
	double t4 = pi[0] * pi[1] + pi[2] * pi[3] + (t2 + t3) * k;
	t4 = 0.1e1 / t4;
	double t5 = -0.1e1 / 0.2e1;
	double t6 = exp(t5 * (t1 * k + pi[0] + pi[1]) * x * t4);
	double t7 = pi[2] + pi[3] + pi[0] + pi[1];
	double t8 = exp(t5 * k * t7 * x * t4);
	double t9 = pow(pi[1], 0.2e1);
	double t10 = pow(pi[0], 0.2e1);
	double t11 = pi[0] + pi[1];
	double t12 = t7 * t6 - t1 * t8 - pi[0] - pi[1];
	double t13 = t8 - 0.1e1;
	double t14 = 0.1e1 / t11;
	double t15 = 0.1e1 / t7;
	double t16 = t13 * pi[2] * t15;
	double t17 = t13 * pi[3] * t15;
	t4 = exp(t5 * (t11 * k + pi[2] + pi[3]) * x * t4);
	t5 = pow(pi[3], 0.2e1);
	double t18 = pow(pi[2], 0.2e1);
	t11 = t11 * t8;
	t7 = -t11 + t7 * t4 - pi[3] - pi[2];
	t1 = 0.1e1 / t1;
	double t19 = t13 * pi[0] * t15;
	t13 = t13 * pi[1] * t15;
	double temp;
	temp = (t6 * t9 + ((pi[0] + pi[3] + pi[2]) * t6 + pi[0]) * pi[1] + t2 * t8 + t10) * t14 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[0][0] = temp;
	temp = -pi[1] * t12 * t14 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[0][1] = temp;
	temp = -t16;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[0][2] = temp;
	temp = -t17;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[0][3] = temp;
	temp = -pi[0] * t12 * t14 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[1][0] = temp;
	temp = (t6 * t10 + ((pi[2] + pi[1] + pi[3]) * t6 + pi[1]) * pi[0] + t3 * t8 + t9) * t14 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[1][1] = temp;
	temp = -t16;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[1][2] = temp;
	temp = -t17;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[1][3] = temp;
	temp = -t19;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[2][0] = temp;
	temp = -t13;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[2][1] = temp;
	temp = (t4 * t5 + ((pi[0] + pi[2] + pi[1]) * t4 + pi[2]) * pi[3] + t11 * pi[2] + t18) * t1 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[2][2] = temp;
	temp = -t7 * pi[3] * t1 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[2][3] = temp;
	temp = -t19;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[3][0] = temp;
	temp = -t13;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[3][1] = temp;
	temp = -t7 * pi[2] * t1 * t15;
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[3][2] = temp;
	temp = (t4 * t18 + ((pi[0] + pi[1] + pi[3]) * t4 + pi[3]) * pi[2] + t11 * pi[3] + t5) * t1 * t15;		
	if(temp>1.0) temp = 1.0; if(temp<1e-100) temp=1e-100;
	ptrans[3][3] = temp;
	return ptrans;
}

/* Use the following in Maple to generate this code: (k is 1/transition:transversion ratio, i.e. k=1/kappa)
 
 M := Matrix([
 [-g-k*(c+t),g,k*c,k*t],
 [a,-a-k*(c+t),k*c,k*t],
 [k*a,k*g,-k*(a+g)-t,t],
 [k*a,k*g,c,-k*(a+g)-c]]);
 R:= simplify(-a*M[1,1]-g*M[2,2]-c*M[3,3]-t*M[4,4]);
 M:=simplify(M/R);
 
 CodeGeneration:-C(subs([a=pi[1],g=pi[2],c=pi[3],t=pi[4]],-n[1]*M[1,1]-n[2]*M[2,2]-n[3]*M[3,3]-n[4]*M[4,4]),optimize,resultname="ptrans");
 
 */ 
double HKY85_expected_rate(const vector<double> &n, const double kappa, const vector<double> &pi) {
	const double k = 1.0/kappa;
	double t2 = pi[1];
	double t3 = pi[2];
	double t4 = k * t3;
	double t5 = pi[3];
	double t6 = k * t5;
	double t9 = pi[0];
	double t11 = t9 * k;
	double t14 = t2 * k;
	double t19 = 0.1e1 / (t9 * t2 + t11 * t3 + t11 * t5 + t14 * t3 + t14 * t5 + t3 * t5);
	return n[0] * (t2 + t4 + t6) * t19 / 0.2e1 + n[1] * (t9 + t4 + t6) * t19 / 0.2e1 + n[2] * (t11 + t14 + t5) * t19 / 0.2e1 + n[3] * (t11 + t14 + t3) * t19 / 0.2e1;
}


/*	For a full description of this algorithm:
	A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences
	Tal Pupko, Itsik Peer, Ron Shamir, and Dan Graur. Mol. Biol. Evol. 17(6):890–896. 2000
 */
mydouble maximum_likelihood_ancestral_sequences(Matrix<Nucleotide> &nuc, marginal_tree &ctree, const double kappa, const vector<double> &pi, vector<int> &pat1, vector<int> &cpat, Matrix<Nucleotide> &node_sequence) {
	mydouble ML(1.0);
	// Every node in the tree has a likelihood attached of the best subtree likelihood, and the sequence eventually identified as the global maximum likelihood estimate
	const int nseq = nuc.nrows();
	const int nnodes = 2*nseq-1;
	const int npat = pat1.size();
	node_sequence = Matrix<Nucleotide>(nnodes,npat,N_ambiguous);
	// subtree_ML[i][j][k] is, for node i, pattern j, the subtree maximum likelihood given the parent node has state k = {A,G,C,T}
	Matrix<mydouble> subtree_ML_element(npat,4,0.0);
	vector< Matrix<mydouble> > subtree_ML(nnodes,subtree_ML_element);
	// path_ML[i][j][k] is, for node i, pattern j, the state of node i that maximizes the subtree likelihood given the parent node has state k = {A,G,C,T}
	Matrix<Nucleotide> path_ML_element(npat,4,N_ambiguous);
	vector< Matrix<Nucleotide> > path_ML(nnodes,path_ML_element);
	// For each node (except the root node), define an HKY85 transition probability matrix
	vector< Matrix<double> > ptrans = compute_HKY85_ptrans(ctree,kappa,pi);
	// Nodes are ordered in the tree first in tip order (0..n-1) then in ascending time order towards the root node (2*n-2)
	// First, do the tips
	int i,j,k,l;
	for(i=0;i<nseq;i++) {
		for(j=0;j<npat;j++) {
			const Nucleotide obs = nuc[i][pat1[j]];
			for(k=0;k<4;k++) {
				// If the parent node's state is k, what is the maximum likelihood of the subtree?
				// And what is the state of the node that achieves that maximum value?
				if(obs==Adenine || obs==Guanine || obs==Cytosine || obs==Thymine) {
					subtree_ML[i][j][k] = ptrans[i][k][obs];
					path_ML[i][j][k] = obs;
				} else if(obs==N_ambiguous) {
					// If multiple equally good paths are possible, the path is chosen in the following order of decreasing preference: A, G, C, T
					subtree_ML[i][j][k] = ptrans[i][k][0];
					path_ML[i][j][k] = (Nucleotide)0;
					for(l=1;l<4;l++) {
						double subtree_ML_l = ptrans[i][k][l];
						if(subtree_ML_l>subtree_ML[i][j][k]) {
							subtree_ML[i][j][k] = subtree_ML_l;
							path_ML[i][j][k] = (Nucleotide)l;
						}
					}
				} else {
					stringstream errTxt;
					errTxt << "maximum_likelihood_ancestral_sequences(): unexpected base " << obs << " (out of range 0-5) in sequence " << i << " pattern " << j;
					error(errTxt.str().c_str());
				}
			}
		}
	}
	// Now the internal nodes, all of which are bifurcating
	for(;i<nnodes;i++) {
		for(j=0;j<npat;j++) {
			const mt_node* d0 = ctree.node[i].descendant[0];
			const mt_node* d1 = ctree.node[i].descendant[1];
			// Check the descendant nodes exist
			if(d0==NULL || d1==NULL) {
				stringstream errTxt;
				errTxt << "maximum_likelihood_ancestral_sequences(): null pointer during Viterbi-like algorithm";
				error(errTxt.str().c_str());
			}
			const int i0 = d0->id;
			const int i1 = d1->id;
			if(i0<0 || i0>=nnodes || i1<0 || i1>=nnodes) {
				stringstream errTxt;
				errTxt << "maximum_likelihood_ancestral_sequences(): node index during Viterbi-like algorithm";
				error(errTxt.str().c_str());
			}
			// Check subtree ML has been computed
			for(l=0;l<4;l++) {
				if(subtree_ML[i0][j][l].iszero() || subtree_ML[i1][j][l].iszero()) {
					stringstream errTxt;
					errTxt << "maximum_likelihood_ancestral_sequences(): uninitialized subtree ML during Viterbi-like algorithm";
					error(errTxt.str().c_str());
				}
			}
			for(k=0;k<4;k++) {
				// If the parent node's state is k, what is the maximum likelihood of the subtree?
				// And what is the state of the node that achieves that maximum value?
				// If multiple equally good paths are possible, the path is chosen in the following order of decreasing preference: A, G, C, T
				subtree_ML[i][j][k] = ptrans[i][k][0]*subtree_ML[i0][j][0]*subtree_ML[i1][j][0];
				path_ML[i][j][k] = (Nucleotide)0;
				for(l=1;l<4;l++) {
					const mydouble subtree_ML_l = ptrans[i][k][l]*subtree_ML[i0][j][l]*subtree_ML[i1][j][l];
					if(subtree_ML_l > subtree_ML[i][j][k]) {
						subtree_ML[i][j][k] = subtree_ML_l;
						path_ML[i][j][k] = (Nucleotide)l;
					}
				}
			}
		}
	}
	// Now work back from root to tips choosing the ML path
	// Start at the root (this is redundant as the root's ancestor has no bearing so subtree_ML[nnodes-1][j][l] and path_ML[nnodes-1][j][l] are the same for different l's)
	for(j=0;j<npat;j++) {
		int best_state = 0;
		mydouble ML_temp = subtree_ML[nnodes-1][j][0];
		for(l=1;l<4;l++) {
			if(subtree_ML[nnodes-1][j][l]>ML) {
				best_state = l;
				ML_temp = subtree_ML[nnodes-1][j][l];
			}
		}
		node_sequence[nnodes-1][j] = path_ML[nnodes-1][j][best_state];
		ML *= pow(ML_temp,cpat[j]);
	}
	for(i=nnodes-2;i>=0;i--) {
		const mt_node* anc = ctree.node[i].ancestor;
		// Check the descendant nodes exist
		if(anc==NULL) {
			stringstream errTxt;
			errTxt << "maximum_likelihood_ancestral_sequences(): null pointer during Viterbi-like algorithm second pass";
			error(errTxt.str().c_str());
		}
		const int ianc = anc->id;
		if(ianc<0 || ianc>=nnodes) {
			stringstream errTxt;
			errTxt << "maximum_likelihood_ancestral_sequences(): node index during Viterbi-like algorithm second pass";
			error(errTxt.str().c_str());
		}
		for(j=0;j<npat;j++) {
			// Check ancestor path has been computed
			const int parent_state = node_sequence[ianc][j];
			if(parent_state==N_ambiguous) {
				stringstream errTxt;
				errTxt << "maximum_likelihood_ancestral_sequences(): uninitialized ML ancestor during Viterbi-like algorithm second pass";
				error(errTxt.str().c_str());
			}
			node_sequence[i][j] = path_ML[i][j][parent_state];
		}
	}
	
	// At the end, convert patterns to sites working from the end of the sequence forwards
	return ML;
}

void write_newick(const marginal_tree &ctree, const vector<string> &all_node_names, const char* file_name) {
	ofstream fout(file_name);
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_newick(): could not open file " << file_name << " for writing";
		error(errTxt.str().c_str());
	}
	write_newick(ctree,all_node_names,fout);
	fout.close();
}

void write_newick(const marginal_tree &ctree, const vector<string> &all_node_names, ofstream &fout) {
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_newick(): could not open file stream for writing";
		error(errTxt.str().c_str());
	}
	const int nnodes = ctree.size;
	if(all_node_names.size()!=nnodes) {
		stringstream errTxt;
		errTxt << "write_newick(): length of node names vector does not equal number of nodes";
		error(errTxt.str().c_str());
	}
	const mt_node* root = (const mt_node*)(&ctree.node[nnodes-1]);
	if(root==NULL) {
		stringstream errTxt;
		errTxt << "write_newick(): null pointer to root";
		error(errTxt.str().c_str());
	}
	const int id = root->id;
	const mt_node* d0 = root->descendant[0];
	const mt_node* d1 = root->descendant[1];
	// Check the descendant nodes exist
	if(d0==NULL || d1==NULL) {
		stringstream errTxt;
		errTxt << "write_newick(): null pointer to root descendant";
		error(errTxt.str().c_str());
	}
	// Write to Newick
	fout << "(";
	write_newick_node(d0,all_node_names,fout);
	fout << ",";
	write_newick_node(d1,all_node_names,fout);
	fout << ")" << all_node_names[id] << ";" << endl;
}

void write_newick_node(const mt_node *node, const vector<string> &all_node_names, ofstream &fout) {
	const int id = node->id;
	const mt_node* d0 = node->descendant[0];
	const mt_node* d1 = node->descendant[1];
	// Check the descendant nodes exist
	if(d0==NULL && d1==NULL) {
		// Node is a tip
		fout << all_node_names[id] << ":" << node->edge_time;
	} else if(d0!=NULL && d1!=NULL) {
		// Node is internal
		fout << "(";
		write_newick_node(d0,all_node_names,fout);
		fout << ",";
		write_newick_node(d1,all_node_names,fout);
		fout << ")" << all_node_names[id] << ":" << node->edge_time;
	} else {
		stringstream errTxt;
		errTxt << "write_newick_node(): node has unexpectedly just one descendant";
		error(errTxt.str().c_str());
	}
}

void write_ancestral_fasta(Matrix<Nucleotide> &nuc, vector<string> &all_node_names, const char* file_name) {
	ofstream fout(file_name);
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_ancestral_fasta(): could not open file " << file_name << " for writing";
		error(errTxt.str().c_str());
	}
	write_ancestral_fasta(nuc,all_node_names,fout);
	fout.close();
}

void write_ancestral_fasta(Matrix<Nucleotide> &nuc, vector<string> &all_node_names, ofstream &fout) {
	static const char AGCTN[5] = {'A','G','C','T','N'};
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_ancestral_fasta(): could not open file stream for writing";
		error(errTxt.str().c_str());
	}
	if(nuc.nrows()!=all_node_names.size()) {
		stringstream errTxt;
		errTxt << "write_ancestral_fasta(): number of sequences (" << nuc.nrows() << ") does not equal number of node labels (" << all_node_names.size() << ")";
		error(errTxt.str().c_str());
	}
	int i,pos;
	for(i=0;i<nuc.nrows();i++) {
		fout << ">" << all_node_names[i] << endl;
		for(pos=0;pos<nuc.ncols();pos++) {
			fout << AGCTN[nuc[i][pos]];
		}
		fout << endl;
	}	
}

void write_position_cross_reference(vector<bool> &iscompat, vector<int> &ipat, const char* file_name) {
	ofstream fout(file_name);
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_position_cross_reference(): could not open file " << file_name << " for writing";
		error(errTxt.str().c_str());
	}
	write_position_cross_reference(iscompat,ipat,fout);
	fout.close();
}

void write_position_cross_reference(vector<bool> &iscompat, vector<int> &ipat, ofstream &fout) {
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_position_cross_reference(): could not open file stream for writing";
		error(errTxt.str().c_str());
	}
	int i,j,pat;
	for(i=0,j=0;i<iscompat.size();i++) {
		pat = -1;
		if(iscompat[i]) {
			if(j>=ipat.size()) {
				stringstream errTxt;
				errTxt << "write_position_cross_reference(): internal inconsistency in number of compatible sizes (" << j+1 << " or more) and number of patterns (" << ipat.size() << ")";
				error(errTxt.str().c_str());
			}
			pat = ipat[j];
			++j;
		}
		if(i>0) fout << ',';
		fout << pat+1;
	}
	fout << endl;
}

mydouble maximum_likelihood_ClonalFrame_branch_allsites(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, vector<ImportationState> &is_imported) {
	mydouble ML(0.0);
	// Store the positions of **all** sites
	is_imported = vector<ImportationState>(iscompat.size(),Unimported);
	// subseq_ML[i][j] is, for position i, the subsequence maximum likelihood given the next position has state j = {Unimported,Imported}
	static Matrix<mydouble> subseq_ML(iscompat.size(),2);
	// path_ML[i][j] is, for position i, the state of position i that maximizes the subsequence likelihood given the next position has state j = {Unimported,Imported}
	static Matrix<ImportationState> path_ML(iscompat.size(),2);
	// Define an HKY85 emission probability matrix for Unimported sites
	static Matrix<mydouble> pemisUnimported;
	pemisUnimported = compute_HKY85_ptrans(branch_length,kappa,pinuc);
	// Define an HKY85 emission probability matrix for Imported sites
	static Matrix<mydouble> pemisImported;
	pemisImported = compute_HKY85_ptrans(import_divergence,kappa,pinuc);
	// Recombination parameters
	const double recrate = rho_over_theta*branch_length;
	const double endrecrate = 1.0/mean_import_length;
	const double totrecrate = recrate+endrecrate;	
	// Equilibrium frequency of unimported and imported sites respectively
	const double pi[2] = {endrecrate/totrecrate,recrate/totrecrate};
	// Define a transition probability matrix
	static Matrix<mydouble> ptrans(2,2,0.0);
	// These probabilities do not change until (i==0)
	ptrans[0][0] = (mydouble)(exp(-totrecrate)+pi[0]*(1-exp(-totrecrate)));
	ptrans[0][1] = (mydouble)(pi[1]*(1-exp(-totrecrate)));
	ptrans[1][1] = (mydouble)(exp(-totrecrate)+pi[1]*(1-exp(-totrecrate)));
	ptrans[1][0] = (mydouble)(pi[0]*(1-exp(-totrecrate)));
	// Beginning at the last variable site, calculate the subsequence maximum likelihood
	int i,j;
	for(i=iscompat.size()-1,j=ipat.size();i>=0;i--) {
		if(i==0) {
			ptrans[0][0] = pi[0];
			ptrans[1][0] = pi[0];
			ptrans[0][1] = pi[1];
			ptrans[1][1] = pi[1];
		}
		// If the previous position's state (leftwards) is j, what is the maximum likelihood of the subsequence from the current position to the last (rightwards)?
		// And what is the state k of the position that achieves that maximum value?
		mydouble UU,UI,IU,II;
		if(iscompat[i]) {
			j--;
			if(j<0) {
				stringstream errTxt;
				errTxt << "maximum_likelihood_ClonalFrame_branch_allsites(): internal inconsistency in tracking informative sites";
				error(errTxt.str().c_str());
			}
			Nucleotide dec = node_nuc[dec_id][ipat[j]];
			Nucleotide anc = node_nuc[anc_id][ipat[j]];
			if(i<iscompat.size()-1) {
				UU = ptrans[0][0]*pemisUnimported[anc][dec]*subseq_ML[i+1][0];
				UI = ptrans[0][1]*  pemisImported[anc][dec]*subseq_ML[i+1][1];
				IU = ptrans[1][0]*pemisUnimported[anc][dec]*subseq_ML[i+1][0];
				II = ptrans[1][1]*  pemisImported[anc][dec]*subseq_ML[i+1][1];
			} else {
				UU = ptrans[0][0]*pemisUnimported[anc][dec];
				UI = ptrans[0][1]*  pemisImported[anc][dec];
				IU = ptrans[1][0]*pemisUnimported[anc][dec];
				II = ptrans[1][1]*  pemisImported[anc][dec];
			}					
		} else {
			if(i<iscompat.size()-1) {
				UU = ptrans[0][0]*subseq_ML[i+1][0];
				UI = ptrans[0][1]*subseq_ML[i+1][1];
				IU = ptrans[1][0]*subseq_ML[i+1][0];
				II = ptrans[1][1]*subseq_ML[i+1][1];
			} else {
				UU = ptrans[0][0];
				UI = ptrans[0][1];
				IU = ptrans[1][0];
				II = ptrans[1][1];
			}					
		}
		subseq_ML[i][0] = (UU>=UI) ? UU : UI;
		path_ML[i][0]   = (UU>=UI) ? Unimported : Imported;
		subseq_ML[i][1] = (IU>=II) ? IU : II;
		path_ML[i][1]   = (IU>=II) ? Unimported : Imported;
	}
	// Beginning at the first variable site, identify the most likely path
	// Sanity check
	if(path_ML[0][0]!=path_ML[0][1]) {
		stringstream errTxt;
		errTxt << "maximum_likelihood_ClonalFrame_branch_allsites(): internal inconsistency when choosing the first importation state in the best path";
		error(errTxt.str().c_str());
	}
	is_imported[0] = path_ML[0][0];
	ML = subseq_ML[0][0];
	for(i=1;i<iscompat.size();i++) {
		const int prev = is_imported[i-1];
		is_imported[i] = path_ML[i][prev];
	}
	
	return ML;
}

mydouble maximum_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, vector<ImportationState> &is_imported) {
	mydouble ML(0.0);
	// Store the positions of compatible sites
	vector<double> position(0);
	int i;
	for(i=0;i<iscompat.size();i++) {
		if(iscompat[i]) {
			position.push_back((double)i);
		}
	}
	return maximum_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,position,ipat,kappa,pinuc,branch_length,rho_over_theta,mean_import_length,import_divergence,is_imported);
}

mydouble maximum_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, vector<ImportationState> &is_imported) {
	mydouble ML(0.0);
	const int npos = position.size();
	is_imported = vector<ImportationState>(npos,Unimported);
	// subseq_ML[i][j] is, for position i, the subsequence maximum likelihood given the next position has state j = {Unimported,Imported}
	static Matrix<mydouble> subseq_ML(npos,2);
	// path_ML[i][j] is, for position i, the state of position i that maximizes the subsequence likelihood given the next position has state j = {Unimported,Imported}
	static Matrix<ImportationState> path_ML(npos,2);
	// Define an HKY85 emission probability matrix for Unimported sites
	static Matrix<mydouble> pemisUnimported;
	pemisUnimported = compute_HKY85_ptrans(branch_length,kappa,pinuc);
	// Define an HKY85 emission probability matrix for Imported sites
	static Matrix<mydouble> pemisImported;
	pemisImported = compute_HKY85_ptrans(import_divergence,kappa,pinuc);
	// Recombination parameters
	const double recrate = rho_over_theta*branch_length;
	const double endrecrate = 1.0/mean_import_length;
	const double totrecrate = recrate+endrecrate;	
	// Equilibrium frequency of unimported and imported sites respectively
	const double pi[2] = {endrecrate/totrecrate,recrate/totrecrate};
	// Define a transition probability matrix
	static Matrix<mydouble> ptrans(2,2,0.0);
	// Beginning at the last variable site, calculate the subsequence maximum likelihood
	int i;
	for(i=npos-1;i>=0;i--) {
		if(i>0) {
			const double expterm = exp(-totrecrate*(position[i]-position[i-1]));
			ptrans[0][0] = (mydouble)(expterm+pi[0]*(1-expterm));
			ptrans[0][1] = (mydouble)(pi[1]*(1-expterm));
			ptrans[1][1] = (mydouble)(expterm+pi[1]*(1-expterm));
			ptrans[1][0] = (mydouble)(pi[0]*(1-expterm));
		} else {
			ptrans[0][0] = pi[0];
			ptrans[1][0] = pi[0];
			ptrans[0][1] = pi[1];
			ptrans[1][1] = pi[1];
		}
		Nucleotide dec = node_nuc[dec_id][ipat[i]];
		Nucleotide anc = node_nuc[anc_id][ipat[i]];
		// If the previous position's state (leftwards) is j, what is the maximum likelihood of the subsequence from the current position to the last (rightwards)?
		// And what is the state k of the position that achieves that maximum value?
		mydouble UU,UI,IU,II;
		if(i<npos-1) {
			UU = ptrans[0][0]*pemisUnimported[anc][dec]*subseq_ML[i+1][0];
			UI = ptrans[0][1]*  pemisImported[anc][dec]*subseq_ML[i+1][1];
			IU = ptrans[1][0]*pemisUnimported[anc][dec]*subseq_ML[i+1][0];
			II = ptrans[1][1]*  pemisImported[anc][dec]*subseq_ML[i+1][1];
		} else {
			UU = ptrans[0][0]*pemisUnimported[anc][dec];
			UI = ptrans[0][1]*  pemisImported[anc][dec];
			IU = ptrans[1][0]*pemisUnimported[anc][dec];
			II = ptrans[1][1]*  pemisImported[anc][dec];
		}
		subseq_ML[i][0] = (UU>=UI) ? UU : UI;
		path_ML[i][0]   = (UU>=UI) ? Unimported : Imported;
		subseq_ML[i][1] = (IU>=II) ? IU : II;
		path_ML[i][1]   = (IU>=II) ? Unimported : Imported;
	}
	// Beginning at the first variable site, identify the most likely path
	// Sanity check
	if(path_ML[0][0]!=path_ML[0][1]) {
		stringstream errTxt;
		errTxt << "maximum_likelihood_ClonalFrame_branch(): internal inconsistency when choosing the first importation state in the best path";
		error(errTxt.str().c_str());
	}
	is_imported[0] = path_ML[0][0];
	ML = subseq_ML[0][0];
	for(i=1;i<npos;i++) {
		const int prev = is_imported[i-1];
		is_imported[i] = path_ML[i][prev];
	}
	
	return ML;
}

double marginal_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence) {
	// Force the mydouble version of the function to be used (defined below)
	mydouble ret = mydouble_marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,iscompat,ipat,kappa,pinuc,branch_length,rho_over_theta,mean_import_length,import_divergence);
	return ret.LOG();
	// Store the positions of compatible sites
	vector<double> position(0);
	int i,npos=0;
	for(i=0;i<iscompat.size();i++) {
		if(iscompat[i]) {
			position.push_back((double)i);
			++npos;
		}
	}
	// Define an HKY85 emission probability matrix for Unimported sites
	static Matrix<double> pemisUnimported;
	pemisUnimported = dcompute_HKY85_ptrans(branch_length,kappa,pinuc);
	// Define an HKY85 emission probability matrix for Imported sites
	static Matrix<double> pemisImported;
	pemisImported = dcompute_HKY85_ptrans(import_divergence,kappa,pinuc);
	// Recombination parameters
	const double recrate = rho_over_theta*branch_length;
	const double endrecrate = 1.0/mean_import_length;
	const double totrecrate = recrate+endrecrate;	
	// Storage
	double aprev[2];
	double a[2];
	// Use these factorizations to avoid under/overflow
	const double ftr = 1.0;
	double sumlogsuma = 0.0;
	// Equilibrium frequency of unimported and imported sites respectively
	const double pi[2] = {endrecrate/totrecrate,recrate/totrecrate};
	// Beginning at the first variable site, calculate the subsequence marginal likelihood
	for(i=0;i<npos;i++) {
		Nucleotide dec = node_nuc[dec_id][ipat[i]];
		Nucleotide anc = node_nuc[anc_id][ipat[i]];
		if(i==0) {
			a[0] = pi[0]*pemisUnimported[anc][dec];
			a[1] = pi[1]*  pemisImported[anc][dec];
		} else {
			aprev[0] = a[0];
			aprev[1] = a[1];
			const double sumaprev = aprev[0]+aprev[1];
			sumlogsuma += log(sumaprev);
			const double prnotrans = exp(-totrecrate*(position[i]-position[i-1]));
			const double prtrans = 1.0-prnotrans;
			a[0] = (aprev[0]/sumaprev*prnotrans+pi[0]*prtrans)*pemisUnimported[anc][dec];
			a[1] = (aprev[1]/sumaprev*prnotrans+pi[1]*prtrans)*  pemisImported[anc][dec];
		}
	}
	// Output
	const double loglik = log(a[0]+a[1])+sumlogsuma-(double)npos*log(ftr);
	return loglik;
}

mydouble mydouble_marginal_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence) {
	// Store the positions of compatible sites
	vector<double> position(0);
	int i;
	for(i=0;i<iscompat.size();i++) {
		if(iscompat[i]) {
			position.push_back((double)i);
		}
	}
	return mydouble_marginal_likelihood_ClonalFrame_branch(dec_id,anc_id,node_nuc,position,ipat,kappa,pinuc,branch_length,rho_over_theta,mean_import_length,import_divergence);
}
	
mydouble mydouble_marginal_likelihood_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence) {
	const int npos = position.size();
	// Define an HKY85 emission probability matrix for Unimported sites
	static Matrix<mydouble> pemisUnimported;
	pemisUnimported = compute_HKY85_ptrans(branch_length,kappa,pinuc);
	// Define an HKY85 emission probability matrix for Imported sites
	static Matrix<mydouble> pemisImported;
	pemisImported = compute_HKY85_ptrans(import_divergence,kappa,pinuc);
	// Recombination parameters
	const double recrate = rho_over_theta*branch_length;
	const double endrecrate = 1.0/mean_import_length;
	const double totrecrate = recrate+endrecrate;	
	// Storage
	mydouble aprev[2];
	mydouble a[2];
	// Equilibrium frequency of unimported and imported sites respectively
	const mydouble pi[2] = {endrecrate/totrecrate,recrate/totrecrate};
	// Beginning at the first variable site, calculate the subsequence marginal likelihood
	int i;
	for(i=0;i<npos;i++) {
		Nucleotide dec = node_nuc[dec_id][ipat[i]];
		Nucleotide anc = node_nuc[anc_id][ipat[i]];
		if(i==0) {
			a[0] = pi[0]*pemisUnimported[anc][dec];
			a[1] = pi[1]*  pemisImported[anc][dec];
		} else {
			aprev[0] = a[0];
			aprev[1] = a[1];
			const mydouble sumaprev = aprev[0]+aprev[1];
			mydouble prnotrans;
			prnotrans.setlog(-totrecrate*(position[i]-position[i-1]));
			const mydouble prtrans = mydouble(1.0)-prnotrans;
			a[0] = (aprev[0]*prnotrans+sumaprev*pi[0]*prtrans)*pemisUnimported[anc][dec];
			a[1] = (aprev[1]*prnotrans+sumaprev*pi[1]*prtrans)*  pemisImported[anc][dec];
		}
	}
	// Output
	const mydouble ML = (a[0]+a[1]);
	return ML;
}

mydouble likelihood_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<int> &pat1, const vector<int> &cpat, const double kappa, const vector<double> &pinuc, const double branch_length) {
	mydouble ML(1.0);
	const int npat = pat1.size();
	// Define an HKY85 emission probability matrix for Unimported sites
	static Matrix<mydouble> pemis;
	pemis = compute_HKY85_ptrans(branch_length,kappa,pinuc);
	// Cycle through the patterns calculating the likelihood
	int i;
	for(i=0;i<npat;i++) {
		Nucleotide dec = node_nuc[dec_id][i];
		Nucleotide anc = node_nuc[anc_id][i];
		ML *= pow(pemis[anc][dec],(double)cpat[i]);
	}
	return ML;
}

bool string_to_bool(const string s, const string label) {
	int i;
	string S = s;
	for(i=0;i<s.length();i++) S[i] = tolower(s[i]);
	if(S=="true") return true;
	if(S=="false") return false;
	stringstream errTxt;
	errTxt << "string_to_bool(): string " << label << " has value " << S << " when true or false expected";
	error(errTxt.str().c_str());
	return false;
}

void write_importation_status(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, const char* file_name, const int root_node) {
	ofstream fout(file_name);
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_importation_status(): could not open file " << file_name << " for writing";
		error(errTxt.str().c_str());
	}
	write_importation_status(imported,all_node_names,isBLC,compat,fout,root_node);
	fout.close();
}

void write_importation_status(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, ofstream &fout, const int root_node) {
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_importation_status(): could not open file stream for writing";
		error(errTxt.str().c_str());
	}
	if(imported.size()!=root_node) {
		stringstream errTxt;
		errTxt << "write_importation_status(): number of lineages (" << imported.size() << ") does not equal the number of non-root node labels (" << root_node << ")";
		error(errTxt.str().c_str());
	}
	if(all_node_names.size()<root_node) {
		stringstream errTxt;
		errTxt << "write_importation_status(): number of non-root lineages (" << root_node << ") exceeds the number of node labels (" << all_node_names.size() << ")";
		error(errTxt.str().c_str());
	}
	int i,pos;
	for(i=0;i<root_node;i++) {
		fout << ">" << all_node_names[i] << endl;
		int k = 0;
		for(pos=0;pos<isBLC.size();pos++) {
			if(isBLC[pos]) {
				// If used in branch length correction, 0 (unimported), 1 (imported), 2 (homoplasy/multiallelic unimported), 3 (homoplasy/multiallelic imported)
				int out = 2*(compat[pos]>0) + (int)imported[i][k];
				fout << out;
				++k;
			} else if(compat[pos]<=0) {
				// If compatible but not used in branch length correction, 4
				fout << 4;
			} else {
				// If homoplasy/multiallelic and not used in branch length correction, 5
				fout << 5;
			}
		}
		fout << endl;
	}	
}

void write_importation_status_intervals(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, const char* file_name, const int root_node) {
	ofstream fout(file_name);
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_importation_status_intervals(): could not open file " << file_name << " for writing";
		error(errTxt.str().c_str());
	}
	write_importation_status_intervals(imported,all_node_names,isBLC,compat,fout,root_node);
	fout.close();
}

void write_importation_status_intervals(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, ofstream &fout, const int root_node) {
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_importation_status_intervals(): could not open file stream for writing";
		error(errTxt.str().c_str());
	}
	if(imported.size()!=root_node) {
		stringstream errTxt;
		errTxt << "write_importation_status_intervals(): number of lineages (" << imported.size() << ") does not equal the number of non-root node labels (" << root_node << ")";
		error(errTxt.str().c_str());
	}
	if(all_node_names.size()<root_node) {
		stringstream errTxt;
		errTxt << "write_importation_status_intervals(): number of non-root lineages (" << root_node << ") exceeds the number of node labels (" << all_node_names.size() << ")";
		error(errTxt.str().c_str());
	}
	const char tab = '\t';
	fout << "Node" << tab << "Beg" << tab << "End" << endl;
	int i,pos;
	for(i=0;i<root_node;i++) {
		// Identify intervals
		bool in_interval = (imported[i][0]==Imported);
		int interval_beg = 0;
		for(pos=1;pos<imported[i].size();pos++) {
			if(in_interval) {
				if(imported[i][pos]==Unimported) {
					fout << all_node_names[i] << tab << interval_beg+1 << tab << pos << endl;
					in_interval = false;
				}
			} else {
				if(imported[i][pos]==Imported) {
					interval_beg = pos;
					in_interval = true;
				}
			}
		}
		if(in_interval) {
			fout << all_node_names[i] << tab << interval_beg+1 << tab << pos << endl;
		}
	}
}



