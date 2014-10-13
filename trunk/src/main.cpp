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
	// Output version number (There were 178 SVN revisions up to version 1.0)
	cout << "ClonalFrameML version 1." << ClonalFrameML_SVNRevision-178 << endl;
	// Process the command line arguments	
	if(argc<5) {
		stringstream errTxt;
		errTxt << "Syntax: ClonalFrameML newick_file fasta_file output_file [OPTIONS]" << endl;
		errTxt << endl;
		errTxt << "Options specifying the analysis type:" << endl;
		errTxt << "-em                            true (default) or false   Estimate parameters by a Baum-Welch expectation maximization algorithm." << endl;
		errTxt << "-embranch                      true or false (default)   Estimate parameters for each branch using the EM algorithm." << endl;
		errTxt << "-rescale_no_recombination      true or false (default)   Rescale branch lengths for given sites with no recombination model." << endl;		
		errTxt << "-imputation_only               true or false (default)   Perform only ancestral state reconstruction and imputation." << endl;
		errTxt << "Options affecting all analyses:" << endl;
		errTxt << "-kappa                         value > 0 (default 2.0)   Relative rate of transitions vs transversions in substitution model" << endl;
		errTxt << "-fasta_file_list               true or false (default)   Take fasta_file to be a white-space separated file list." << endl;
		errTxt << "-ignore_user_sites             sites_file                Ignore sites listed in whitespace-separated sites_file." << endl;
		errTxt << "-ignore_incomplete_sites       true or false (default)   Ignore sites with any ambiguous bases." << endl;
		errTxt << "-use_incompatible_sites        true (default) or false   Use homoplasious and multiallelic sites to correct branch lengths." << endl;
		errTxt << "-show_progress                 true or false (default)   Output the progress of the maximum likelihood routines." << endl;
		errTxt << "-min_branch_length             value > 0 (default 1e-7)  Minimum branch length." << endl;
		errTxt << "-reconstruct_invariant_sites   true or false (default)   Reconstruct the ancestral states at invariant sites." << endl;
//		errTxt << "-compress_reconstructed_sites  true (default) or false   Reduce the number of columns in the output FASTA file." << endl;	// Alternative not currently implemented, so not optional
		errTxt << "Options affecting -em and -embranch:" << endl;
		errTxt << "-prior_mean                    df \"0.1 0.001 0.1 0.0001\" Prior mean for R/theta, 1/delta, nu and M." << endl;
		errTxt << "-prior_sd                      df \"0.1 0.001 0.1 0.0001\" Prior standard deviation for R/theta, 1/delta, nu and M." << endl;
		errTxt << "-initial_values                default \"0.1 0.001 0.05\"  Initial values for R/theta, 1/delta and nu." << endl;
		errTxt << "-guess_initial_m               true (default) or false   Initialize M and nu jointly in the EM algorithms." << endl;
		errTxt << "-emsim                         value >= 0  (default 0)   Number of simulations to estimate uncertainty in the EM results." << endl;
		errTxt << "-embranch_dispersion           value > 0 (default .01)   Dispersion in parameters among branches in the -embranch model." << endl;
		errTxt << "Options affecting -rescale_no_recombination:" << endl;
		errTxt << "-brent_tolerance               tolerance (default .001)  Set the tolerance of the Brent routine for -rescale_no_recombination." << endl;
		errTxt << "-powell_tolerance              tolerance (default .001)  Set the tolerance of the Powell routine for -rescale_no_recombination." << endl;
		error(errTxt.str().c_str());
	}
	// Process required arguments
	const char* newick_file = argv[1];
	const char* fasta_file = argv[2];
	const char* out_file = argv[3];
	string tree_out_file = string(out_file) + ".labelled_tree.newick";
	string fasta_out_file = string(out_file) + ".ML_sequence.fasta";
	string xref_out_file = string(out_file) + ".position_cross_reference.txt";
	string import_out_file = string(out_file) + ".importation_status.txt";
	string em_out_file = string(out_file) + ".em.txt";
	string emsim_out_file = string(out_file) + ".emsim.txt";
	// Set default options
	ArgumentWizard arg;
	arg.case_sensitive = false;
	string fasta_file_list="false", imputation_only="false", ignore_incomplete_sites="false", ignore_user_sites="", reconstruct_invariant_sites="false";
	string use_incompatible_sites="true", rescale_no_recombination="false";
	string show_progress="false", compress_reconstructed_sites="true";
	string string_prior_mean="0.1 0.001 0.1 0.0001", string_prior_sd="0.1 0.001 0.1 0.0001", string_initial_values = "0.1 0.001 0.05";
	string guess_initial_m="true", em="true", embranch="false";
	double brent_tolerance = 1.0e-3, powell_tolerance = 1.0e-3, global_min_branch_length = 1.0e-7;
	double embranch_dispersion = 0.01, kappa = 2.0;
	int emsim = 0;
	// Process options
	arg.add_item("fasta_file_list",				TP_STRING, &fasta_file_list);
	arg.add_item("imputation_only",				TP_STRING, &imputation_only);
	arg.add_item("ignore_incomplete_sites",		TP_STRING, &ignore_incomplete_sites);
	arg.add_item("ignore_user_sites",			TP_STRING, &ignore_user_sites);
	arg.add_item("reconstruct_invariant_sites", TP_STRING, &reconstruct_invariant_sites);
	arg.add_item("use_incompatible_sites",		TP_STRING, &use_incompatible_sites);
	arg.add_item("brent_tolerance",				TP_DOUBLE, &brent_tolerance);
	arg.add_item("powell_tolerance",			TP_DOUBLE, &powell_tolerance);
	arg.add_item("rescale_no_recombination",	TP_STRING, &rescale_no_recombination);
	arg.add_item("show_progress",				TP_STRING, &show_progress);
	arg.add_item("compress_reconstructed_sites",TP_STRING, &compress_reconstructed_sites);	
	arg.add_item("min_branch_length",			TP_DOUBLE, &global_min_branch_length);
	arg.add_item("prior_mean",					TP_STRING, &string_prior_mean);
	arg.add_item("prior_sd",					TP_STRING, &string_prior_sd);
	arg.add_item("initial_values",				TP_STRING, &string_initial_values);
	arg.add_item("guess_initial_m",				TP_STRING, &guess_initial_m);
	arg.add_item("em",							TP_STRING, &em);
	arg.add_item("emsim",						TP_INT,	   &emsim);
	arg.add_item("embranch",					TP_STRING, &embranch);
	arg.add_item("embranch_dispersion",			TP_DOUBLE, &embranch_dispersion);
	arg.add_item("kappa",						TP_DOUBLE, &kappa);
	arg.read_input(argc-3,argv+3);
	bool FASTA_FILE_LIST				= string_to_bool(fasta_file_list,				"fasta_file_list");
	bool CORRECT_BRANCH_LENGTHS			= !string_to_bool(imputation_only,				"imputation_only");
	bool IGNORE_INCOMPLETE_SITES		= string_to_bool(ignore_incomplete_sites,		"ignore_incomplete_sites");
	bool RECONSTRUCT_INVARIANT_SITES	= string_to_bool(reconstruct_invariant_sites,	"reconstruct_invariant_sites");
	bool USE_INCOMPATIBLE_SITES			= string_to_bool(use_incompatible_sites,		"use_incompatible_sites");
	bool RESCALE_NO_RECOMBINATION		= string_to_bool(rescale_no_recombination,		"rescale_no_recombination");
	bool SHOW_PROGRESS					= string_to_bool(show_progress,					"show_progress");
	bool COMPRESS_RECONSTRUCTED_SITES	= string_to_bool(compress_reconstructed_sites,	"compress_reconstructed_sites");
	bool GUESS_INITIAL_M				= string_to_bool(guess_initial_m,				"guess_initial_m");
	bool EM								= string_to_bool(em,							"em");
	bool EMBRANCH						= string_to_bool(embranch,						"embranch");
	bool MULTITHREAD = false;
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
	if(((int)RESCALE_NO_RECOMBINATION) + (int)EM +(int)EMBRANCH>1) {
		stringstream errTxt;
		errTxt << "rescale_no_recombination, em and embranch are mutually incompatible";
		error(errTxt.str().c_str());
	}
	if((RESCALE_NO_RECOMBINATION || EM || EMBRANCH) && !CORRECT_BRANCH_LENGTHS) {
		stringstream wrnTxt;
		wrnTxt << "advanced options will be ignored because imputation_only=true";
		warning(wrnTxt.str().c_str());
	}
	if(CORRECT_BRANCH_LENGTHS && !(RESCALE_NO_RECOMBINATION || EM || EMBRANCH)) {
		error("One of -em, -embranch or -rescale_no_recombination must be specified when imputation_only=false");
	}
	if(MULTITHREAD) {
		cout << "WARNING: multithreaded version not implemented, ignoring." << endl;
	}
	if(global_min_branch_length<=0.0) {
		error("Minimum branch length must be positive");
	}
	// Process the prior mean and standard deviation
	vector<double> prior_mean(0), prior_sd(0);
	stringstream sstream_prior_mean;
	sstream_prior_mean << string_prior_mean;
	int i;
	for(i=0;i<1000;i++) {
		if(sstream_prior_mean.eof()) break;
		double prior_mean_elem;
		sstream_prior_mean >> prior_mean_elem;
		if(sstream_prior_mean.fail()) error("Could not interpret value specified by prior_mean");
		prior_mean.push_back(prior_mean_elem);
	}
	if(i==1000) error("Maximum length of vector exceeded by prior_mean");
	stringstream sstream_prior_sd;
	sstream_prior_sd << string_prior_sd;
	for(i=0;i<1000;i++) {
		if(sstream_prior_sd.eof()) break;
		double prior_sd_elem;
		sstream_prior_sd >> prior_sd_elem;
		if(sstream_prior_sd.fail()) error("Could not interpret value specified by prior_sd");
		prior_sd.push_back(prior_sd_elem);
	}
	if(prior_mean.size()!=4) error("prior_mean must have 4 values separated by spaces");
	if(prior_sd.size()!=4) error("prior_sd must have 4 values separated by spaces");
	// Process the initial values
	vector<double> initial_values(0);
	if(string_initial_values!="") {
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
		if(!(initial_values.size()==3)) error("initial values must have 3 values separated by spaces");
	}
	if(emsim<0) error("-emsim cannot be negative");
	if(emsim>0 && !(EM || EMBRANCH)) error("-emsim only applicable with -em or -embranch");
	if(embranch_dispersion<=0.0) error("-embranch_dispersion must be positive");
	if(kappa<=0.0) error("-kappa must be positive");
	
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
		// Sanity check: are all branch lengths non-negative
		for(i=0;i<ctree_node_labels.size();i++) {
			if(ctree.node[i].edge_time<0.0) {
				stringstream errTxt;
				errTxt << "Negative branch length of " << ctree.node[i].edge_time << " found for branch " << ctree_node_labels[i];
				error(errTxt.str().c_str());
			}
		}
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
		
		cout << "BRANCH LENGTH CORRECTION/RECOMBINATION ANALYSIS:" << endl;
		cout << "Analysing " << nBLC << " sites" << endl;
		// Report the estimated equilibrium frequencies
		cout << "Empirical nucleotide frequencies:   A " << round(1000*empirical_nucleotide_frequencies[Adenine])/10 << "%   C " << round(1000*empirical_nucleotide_frequencies[Cytosine])/10;
		cout << "%   G " << round(1000*empirical_nucleotide_frequencies[Guanine])/10 << "%   T " << round(1000*empirical_nucleotide_frequencies[Thymine])/10 << "%" << endl;
		
		if(RESCALE_NO_RECOMBINATION) {
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
				// Update branch length in the output tree
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = final_branch_length;
				cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L = " << -Pow.function_minimum << " M = " << final_branch_length << endl;
				ML += -Pow.function_minimum;
			}
			cout << "Log-likelihood after branch optimization is " << ML << endl;
		} else if(EM) {
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) AND recombination parameters under the ClonalFrame model
			// using Baum-Welch EM algorithm
			cout << "Beginning branch optimization. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L   maximum unnormalized log-posterior per branch" << endl;
			cout << "R   R/theta per branch                                       (> 0)" << endl;
			cout << "I   mean DNA import length per branch                        (> 0)" << endl;
			cout << "D   divergence of DNA imported by recombination              (> 0)" << endl;			
			cout << "M   expected number of mutations per branch                  (> 0)" << endl;
			double ML = 0.0;
			vector< vector<ImportationState> > is_imported(root_node);
			// Calculate the a and b parameters of the priors
			vector<double> prior_a(4), prior_b(4);
			for(i=0;i<4;i++) {
				// Mean = a/b and variance = a/b/b so sd = sqrt(a)/b
				// So b = mean/sd/sd and a = b*mean
				if(prior_mean[i]<=0.0) error("EM: prior_mean must be positive");
				if(prior_sd[i]<=0.0) error("EM: prior_sd must be positive");
				prior_b[i] = prior_mean[i]/prior_sd[i]/prior_sd[i];
				prior_a[i] = prior_b[i]*prior_mean[i];
			}
			// Initial values for R_over_theta, mean_import_length and import_divergence from prior
			vector<double> param(3);
			param[0] = initial_values[0];
			param[1] = 1.0/initial_values[1];
			param[2] = initial_values[2];
			// Do inference
			clock_t pow_start_time = clock();
			ClonalFrameBaumWelch cff(ctree,node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,is_imported,prior_a,prior_b,root_node,GUESS_INITIAL_M,SHOW_PROGRESS);
			param = cff.maximize_likelihood(param);
			ML = cff.ML;
			cout << " L = " << ML << " R = " << param[0] << " I = " << param[1] << " D = " << param[2] << " in " << (double)(clock()-pow_start_time)/CLOCKS_PER_SEC << " s and " << cff.neval << " evaluations" << endl;
			cout << " Posterior alphas: R = " << cff.posterior_a[0] << " I = " << cff.posterior_a[1] << " D = " << cff.posterior_a[2] << endl;
			for(i=0;i<root_node;i++) {
				if(cff.informative[i]) {
					cout << "Branch " << ctree_node_labels[i] << " B = " << cff.initial_branch_length[i] << " M = " << param[3+i] << endl;
					// Update branch length in the output tree
					// Note this is unsafe in general because the corresponding node times are not adjusted
					ctree.node[i].edge_time = param[3+i];
				} else {
					ctree.node[i].edge_time = global_min_branch_length;
				}
			}
			// Output the point estimates, 95% credible intervals, posterior_a and posterior_b parameters
			ofstream vout(em_out_file.c_str());
			char tab = '\t';
			vout << "Parameter" << tab << "Posterior Mean" << tab << "Posterior Variance" << tab << "a_post" << tab << "b_post" << endl;
			vout << "R/theta"	<< tab << param[0] << tab << param[0]*param[0]/cff.posterior_a[0] << tab << cff.posterior_a[0] << tab << cff.posterior_a[0]/param[0] << endl;
			vout << "1/delta"	<< tab << 1./param[1] << tab << 1./param[1]/param[1]/cff.posterior_a[1] << tab << cff.posterior_a[1] << tab << cff.posterior_a[1]*param[1] << endl;
			vout << "nu"		<< tab << param[2] << tab << param[2]*param[2]/cff.posterior_a[2] << tab << cff.posterior_a[2] << tab << cff.posterior_a[2]/param[2] << endl;
			for(i=0;i<root_node;i++) {
				if(cff.informative[i]) {
					vout << ctree_node_labels[i] << tab << param[3+i] << tab << param[3+i]*param[3+i]/cff.posterior_a[3+i] << tab << cff.posterior_a[3+i] << tab << cff.posterior_a[3+i]/param[3+i] << endl;
				}
			}
			vout.close();
			// Output the importation status
			write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			cout << "Wrote inferred importation status to " << import_out_file << endl;
			
			// If required, simulate under the point estimates and output posterior samples of the parameters
			if(emsim>0) {
				Matrix<double> sim = cff.simulate_posterior(param,emsim);
				if(sim.nrows()!=3 || sim.ncols()!=emsim) error("ClonalFrameBaumWelch::simulate_posterior() produced unexpected results");
				ofstream eout(emsim_out_file.c_str());
				eout << "R/theta" << tab << "delta" << tab << "nu" << endl;
				for(i=0;i<emsim;i++) {
					eout << sim[0][i] << tab << sim[1][i] << tab << sim[2][i] << endl;
				}
				eout.close();				
			}
			
		} else if(EMBRANCH) {
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) AND recombination parameters
			// for each branch under the ClonalFrame model using Baum-Welch EM algorithm
			cout << "Beginning branch optimization. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L   maximum unnormalized log-posterior per branch" << endl;
			cout << "R   R/theta per branch                                       (> 0)" << endl;
			cout << "I   mean DNA import length per branch                        (> 0)" << endl;
			cout << "D   divergence of DNA imported by recombination              (> 0)" << endl;			
			cout << "M   expected number of mutations per branch                  (> 0)" << endl;
			double ML = 0.0;
			vector< vector<ImportationState> > is_imported(root_node);
			// Calculate the a and b parameters of the priors
			vector<double> prior_a(5), prior_b(5);
			for(i=0;i<4;i++) {
				// Mean = a/b and variance = a/b/b so sd = sqrt(a)/b
				// So b = mean/sd/sd and a = b*mean
				if(prior_mean[i]<=0.0) error("EMBRANCH: prior_mean must be positive");
				if(prior_sd[i]<=0.0) error("EMBRANCH: prior_sd must be positive");
				prior_b[i] = prior_mean[i]/prior_sd[i]/prior_sd[i];
				prior_a[i] = prior_b[i]*prior_mean[i];
			}
			// Set the prior on the fifth parameter
			prior_a[4] = prior_b[4] = 1.0/embranch_dispersion;
			// Initial values for rho_over_theta, mean_import_length and import_divergence from prior
			// Note that the fourth value (mean branch length) is ignored and computed from the tree
			vector<double> param(4);
			param[0] = initial_values[0];
			param[1] = 1.0/initial_values[1];
			param[2] = initial_values[2];
			param.push_back(1.0e-5);
			// Do inference
			clock_t pow_start_time = clock();
			ClonalFrameBaumWelchRhoPerBranch cff(ctree,node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,is_imported,prior_a,prior_b,root_node,GUESS_INITIAL_M,SHOW_PROGRESS);
			cff.maximize_likelihood(param);
			ML = cff.ML;
			cout << "Mean parameters:" << endl;
			cout << " L = " << ML << " R = " << cff.mean_param[0] << " I = " << 1.0/cff.mean_param[1] << " D = " << cff.mean_param[2] << " M = " << cff.mean_param[3] << " in " << (double)(clock()-pow_start_time)/CLOCKS_PER_SEC << " s and " << cff.neval << " evaluations" << endl;
			cout << "Parameters per branch:" << endl;
			for(i=0;i<root_node;i++) {
				if(cff.informative[i]) {
					cout << "Branch " << ctree_node_labels[i] << " B = " << cff.initial_branch_length[i] << " R = " << cff.mean_param[0]*cff.full_param[i][0] << " I = " << 1.0/(cff.mean_param[1]*cff.full_param[i][1]) << " D = " << cff.mean_param[2]*cff.full_param[i][2] << " M = " << cff.mean_param[3]*cff.full_param[i][3] << endl;
					// Update branch length in the output tree
					// Note this is unsafe in general because the corresponding node times are not adjusted
					ctree.node[i].edge_time = cff.mean_param[3]*cff.full_param[i][3];
				} else {
					ctree.node[i].edge_time = global_min_branch_length;
				}
			}
			// Output the point estimates, 95% credible intervals, posterior_a and posterior_b parameters
			ofstream vout(em_out_file.c_str());
			char tab = '\t';
			vout << "Parameter" << tab << "Branch" << tab << "Posterior Mean" << tab << "Posterior Variance" << tab << "a_post" << tab << "b_post" << endl;
			int p;
			string pname[4] = {"R/theta","1/delta","nu","M"};
			for(p=0;p<4;p++) {
				vout << pname[p] << tab << "Mean" << tab <<cff.mean_param[p] << tab << "NA" << tab << "NA" << tab << "NA" << endl;
			}
			for(i=0;i<root_node;i++) {
				if(cff.informative[i]) {
					for(p=0;p<4;p++) {
						vout << pname[p] << tab << ctree_node_labels[i] << tab << cff.full_param[i][p] << tab << cff.full_param[i][p]*cff.full_param[i][p]/cff.posterior_a[i][p] << tab << cff.posterior_a[i][p] << tab << cff.posterior_a[i][p]/cff.full_param[i][p] << endl;
					}
				}
			}
			vout.close();
			// Output the importation status
			write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str(),root_node);
			cout << "Wrote inferred importation status to " << import_out_file << endl;
			
			// If required, simulate under the point estimates and output posterior samples of the parameters
			if(emsim>0) {
				warning("-emsim not yet implemented for -embranch");
//				Matrix<double> sim = cff.simulate_posterior(param,emsim);
//				if(sim.nrows()!=3 || sim.ncols()!=emsim) error("ClonalFrameBaumWelch::simulate_posterior() produced unexpected results");
//				ofstream eout(emsim_out_file.c_str());
//				eout << "R/theta" << tab << "delta" << tab << "nu" << endl;
//				for(i=0;i<emsim;i++) {
//					eout << sim[0][i] << tab << sim[1][i] << tab << sim[2][i] << endl;
//				}
//				eout.close();				
			}
		} else {
			error("No option specified");
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

mydouble likelihood_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<int> &pat1, const vector<int> &cpat, const double kappa, const vector<double> &pinuc, const double branch_length) {
	mydouble ML(1.0);
	const int npat = pat1.size();
	// Define an HKY85 emission probability matrix for Unimported sites
	Matrix<mydouble> pemis;
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

mydouble maximum_likelihood_ClonalFrame_branch_allsites(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<bool> &iscompat, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, vector<ImportationState> &is_imported) {
	mydouble ML(0.0);
	// Store the positions of **all** sites
	is_imported = vector<ImportationState>(iscompat.size(),Unimported);
	// subseq_ML[i][j] is, for position i, the subsequence maximum likelihood given the next position has state j = {Unimported,Imported}
	Matrix<mydouble> subseq_ML(iscompat.size(),2);
	// path_ML[i][j] is, for position i, the state of position i that maximizes the subsequence likelihood given the next position has state j = {Unimported,Imported}
	Matrix<ImportationState> path_ML(iscompat.size(),2);
	// Define an HKY85 emission probability matrix for Unimported sites
	Matrix<mydouble> pemisUnimported;
	pemisUnimported = compute_HKY85_ptrans(branch_length,kappa,pinuc);
	// Define an HKY85 emission probability matrix for Imported sites
	Matrix<mydouble> pemisImported;
	pemisImported = compute_HKY85_ptrans(import_divergence,kappa,pinuc);
	// Recombination parameters
	const double recrate = rho_over_theta*branch_length;
	const double endrecrate = 1.0/mean_import_length;
	const double totrecrate = recrate+endrecrate;	
	// Equilibrium frequency of unimported and imported sites respectively
	const double pi[2] = {endrecrate/totrecrate,recrate/totrecrate};
	// Define a transition probability matrix
	Matrix<mydouble> ptrans(2,2,0.0);
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

// The following function calculates, for a particular branch of the tree, the expected number of transitions from state i to state j and emissions from state i to observation j
// This requires storage for the forward algorithm calculations and a second pass using the backward algorithm to calculate the marginal expectations
// The marginal likelihood for the branch is returned
mydouble mydouble_forward_backward_expectations_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, Matrix<double> &numEmis, vector<double> &denEmis, Matrix<double> &numTrans, vector<double> &denTrans) {
	const int npos = position.size();
	// Define an HKY85 emission probability matrix for Unimported sites
	Matrix<mydouble> pemisUnimported;
	pemisUnimported = compute_HKY85_ptrans(branch_length,kappa,pinuc);
	// Define an HKY85 emission probability matrix for Imported sites
	Matrix<mydouble> pemisImported;
	pemisImported = compute_HKY85_ptrans(import_divergence,kappa,pinuc);
	// Define storage space for the intermediate forward calculations
	Matrix<mydouble> A;
	A = Matrix<mydouble>(npos,2);
	// Resize if necessary and zero the output objects
	numEmis = Matrix<double>(2,2,0.0);
	denEmis = vector<double>(2,0.0);
	numTrans = Matrix<double>(2,2,0.0);
	denTrans = vector<double>(2,0.0);
	//	cout << "numTrans = " << numTrans[0][0].todouble() << " " << numTrans[0][1].todouble() << " " << numTrans[1][0].todouble() << " " << numTrans[0][0].todouble() << endl;
	// Recombination parameters
	const double recrate = rho_over_theta*branch_length;
	const double endrecrate = 1.0/mean_import_length;
	const double totrecrate = recrate+endrecrate;	
	// Transient storage
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
		// Store for the second pass
		A[i][0] = a[0];
		A[i][1] = a[1];
	}
	// Record the marginal likelihood for output later
	const mydouble ML = (a[0]+a[1]);
	// Second pass: backward algorithm
	mydouble bnext[2];
	mydouble b[2];
	// Beginning at the last variable site, calculate the marginal likelihood of the 3prime sites
	for(i=npos-1;i>=0;i--) {
		if(i==(npos-1)) {
			b[0] = mydouble(1.0);
			b[1] = mydouble(1.0);
			
			// Update the expected number of emissions
			mydouble pU = A[i][0]*b[0];
			mydouble pI = A[i][1]*b[1];
			// NB:- pU+pI should always equal ML but just in case it introduces small errors
			const mydouble MLi = pU + pI;
			pU /= MLi;
			pI /= MLi;
			const double ppost[2]  = {pU.todouble(),1.0-pU.todouble()};
			// Increment the numerator and denominator of the expected number of emissions from state j to observation k
			int j;
			// NB:- *** obs refers to the PRESENT site !!! ***
			const int obs = (int)(node_nuc[dec_id][ipat[i]]!=node_nuc[anc_id][ipat[i]]);		// 0 = same, 1 = different
			for(j=0;j<2;j++) {
				// Total number of emissions from j to k equals indicator of actual observation k (0 or 1) weighted by probability the site was in state j
				numEmis[j][obs] += ppost[j];
				// Total number of possible emissions from j to k equals the number of sites, each weighted by probability the site was in state j
				denEmis[j]      += ppost[j];		// NB:- the denominator is the same for both observation states
			}			
		} else {
			bnext[0] = b[0];
			bnext[1] = b[1];
			// Note that these retrieve the ancestral and descendant nucleotides at the 3prime adjacent site
			Nucleotide dec = node_nuc[dec_id][ipat[i+1]];
			Nucleotide anc = node_nuc[anc_id][ipat[i+1]];
			const mydouble pemisU = pemisUnimported[anc][dec];
			const mydouble pemisI = pemisImported[anc][dec];
			mydouble prnotrans;
			prnotrans.setlog(-totrecrate*(position[i+1]-position[i]));
			const mydouble prtrans = mydouble(1.0)-prnotrans;
			const mydouble sumbnext = prtrans*(pi[0]*pemisU*bnext[0] + pi[1]*pemisI*bnext[1]);
			b[0] = prnotrans*pemisU*bnext[0]+sumbnext;
			b[1] = prnotrans*pemisI*bnext[1]+sumbnext;
			
			// Update the expected number of transitions and emissions
			// Calculate the marginal probabilities that the hidden state is Unimported or Imported
			//			if(fabs((A[i][0]*b[0]+A[i][1]*b[1]).LOG()-ML.LOG())>1e-6) {
			//				cout << ML.LOG() << "\t" << (A[i][0]*b[0]+A[i][1]*b[1]).LOG() << endl;
			//			}
			mydouble pU = A[i][0]*b[0];
			mydouble pI = A[i][1]*b[1];
			// NB:- pU+pI should always equal ML but just in case it introduces small errors
			const mydouble MLi = pU + pI;
			pU /= MLi;
			pI /= MLi;
			const double ppost[2]  = {pU.todouble(),1.0-pU.todouble()};
			// Increment the numerator and denominator of the expected number of emissions from state j to observation k
			int j;
			// NB:- *** obs refers to the PRESENT site !!! ***
			const int obs = (int)(node_nuc[dec_id][ipat[i]]!=node_nuc[anc_id][ipat[i]]);		// 0 = same, 1 = different
			for(j=0;j<2;j++) {
				// Total number of emissions from j to k equals indicator of actual observation k (0 or 1) weighted by probability the site was in state j
				numEmis[j][obs] += ppost[j];
				// Total number of possible emissions from j to k equals the number of sites, each weighted by probability the site was in state j
				denEmis[j]      += ppost[j];		// NB:- the denominator is the same for both observation states
			}
			// Increment the numerator and denominator of the expected number of transitions from state j to state k
			// Impose maximum adjacent site distance of 1kb (needed for small-p Poisson approximation to heterogeneous bernoulli)
			const mydouble pemis[2]  = {pemisU,pemisI};
			const double dist = position[i+1]-position[i];
			if(dist<=1000.) {
				int k;
				for(j=0;j<2;j++) {
					for(k=0;k<2;k++) {
						const int istrans = (int)(j!=k);
						// Probability of transition from j to k given the data equals the joint likelihood of the data and transition from j to k, divided by marginal likelihood of the data
						if(istrans) {
							numTrans[j][k] += (A[i][j]*prtrans*pi[k]*pemis[k]*bnext[k]/MLi).todouble();		// Note the use of bnext, not b
							//							if(j==0 && k==1) cout << "pos = " << i << " numTrans[0][1] = " << numTrans[j][k].todouble() << endl; //(A[i][j]*ptrans[istrans]*pemis[k]*bnext[k]/ML).LOG() << endl;
						} else {
							numTrans[j][k] += (A[i][j]*(prnotrans+prtrans*pi[k])*pemis[k]*bnext[k]/MLi).todouble();		// Note the use of bnext, not b
						}
					}
					// Expected distance between sites equals actual distance weighted by the probability the 5prime site was in state j
					denTrans[j] += dist*ppost[j];											// NB:- the denominator is the same for both destination states
				}
			}
		}
	}
	// Return the marginal likelihood
	//	cout << "numTrans = " << numTrans[0][0].todouble() << " " << numTrans[0][1].todouble() << " " << numTrans[1][0].todouble() << " " << numTrans[0][0].todouble() << endl;
	return ML;
}

double Baum_Welch(const marginal_tree &tree, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const vector<bool> &informative, const vector<double> &prior_a, const vector<double> &prior_b, vector<double> &full_param, vector<double> &posterior_a, int &neval, const bool coutput) {
	int i;
	if(coutput) cout << setprecision(9);
	// Initial parameters
	double rho_over_theta = full_param[0];
	double mean_import_length = full_param[1];
	double import_divergence = full_param[2];
	posterior_a = vector<double>(3+informative.size());
	// Storage for the expected number of transitions and emissions in the HMM
	Matrix<double> numEmiss(2,2), numTrans(2,2);
	vector<double> denEmiss(2),   denTrans(2);
	// Counters
	double mutI=0.0;			// Running total divergence at imported sites
	double numU=0.0, numI=0.0;	// Running total number of transitions *to* unimported, imported regions
	double nsiI=0.0;			// Running total number of imported sites
	double lenU=0.0, lenI=0.0;	// Running total length of unimported, imported regions
	// Calculate the marginal likelihood and expected number of transitions and emissions by the forward-backward algorithm
	// Include the effect of the prior
	double ML = gamma_loglikelihood(full_param[0], prior_a[0], prior_b[0]) + gamma_loglikelihood(1.0/full_param[1], prior_a[1], prior_b[1]) + gamma_loglikelihood(full_param[2], prior_a[2], prior_b[2]);
	for(i=0;i<informative.size();i++) {
		if(informative[i]) {
			ML += gamma_loglikelihood(full_param[3+i], prior_a[3], prior_b[3]);
			const int dec_id = tree.node[i].id;
			const int anc_id = tree.node[i].ancestor->id;
			const double branch_length = full_param[3+i];
			ML += mydouble_forward_backward_expectations_ClonalFrame_branch(dec_id,anc_id,node_nuc,position,ipat,kappa,pinuc,branch_length,rho_over_theta,mean_import_length,import_divergence,numEmiss,denEmiss,numTrans,denTrans).LOG();
			// Update estimate of the branch length
			const double mutU_br = numEmiss[0][1];
			const double nsiU_br = denEmiss[0];
			full_param[3+i] = (prior_a[3]+mutU_br)/(prior_b[3]+nsiU_br);
			posterior_a[3+i] = (prior_a[3]+mutU_br);
			// Increment counters for the other expectations
			mutI += numEmiss[1][1];
			nsiI += denEmiss[1];
			const double numI_br = numTrans[0][1];
			const double lenU_br = denTrans[0];
			numI += numI_br;
			lenU += full_param[3+i]*lenU_br;
			numU += numTrans[1][0];
			lenI += denTrans[1];
			if(coutput) {
				cout << "nmut = " << mutU_br << " nU = " << nsiU_br << " nsub = " << numEmiss[1][1] << " nI = " << denEmiss[1] << endl;
				cout << "nU>I = " << numI_br << " dU = " << lenU_br << " nI>U = " << numTrans[1][0] << " dI = " << denTrans[1] << endl;
				cout << "numTrans = " << numTrans[0][0] << " " << numTrans[0][1] << " " << numTrans[1][0] << " " << numTrans[0][0] << endl;
			}
		}
	}
	++neval;
	// Update estimates of the recombination parameters
	full_param[0] = (prior_a[0]+numI)/(prior_b[0]+lenU);
	full_param[1] = (prior_b[1]+lenI)/(prior_a[1]+numU);
	full_param[2] = (prior_a[2]+mutI)/(prior_b[2]+nsiI);
	posterior_a[0] = (prior_a[0]+numI);
	posterior_a[1] = (prior_a[1]+numU);
	posterior_a[2] = (prior_a[2]+mutI);
	if(coutput) {
		cout << "params =";
		for(int j=0;j<full_param.size();j++) cout << " " << full_param[j];
		cout << " ML = " << ML << endl;
	}
	// Iterate until the maximum likelihood improves by less than some threshold
	const int maxit = 200;
	const double threshold = 1.0e-2;
	vector<double> MLE;
	int it;
	for(it=0;it<maxit;it++) {
		// Identify the model parameters
		rho_over_theta = full_param[0];
		mean_import_length = full_param[1];
		import_divergence = full_param[2];
		// Update the likelihood
		mutI=0.0; numU=0.0; numI=0.0; nsiI=0.0; lenU=0.0; lenI=0.0;
		double new_ML = gamma_loglikelihood(full_param[0], prior_a[0], prior_b[0]) + gamma_loglikelihood(1.0/full_param[1], prior_a[1], prior_b[1]) + gamma_loglikelihood(full_param[2], prior_a[2], prior_b[2]);
		for(i=0;i<informative.size();i++) {
			if(informative[i]) {
				new_ML += gamma_loglikelihood(full_param[3+i], prior_a[3], prior_b[3]);
				const int dec_id = tree.node[i].id;
				const int anc_id = tree.node[i].ancestor->id;
				const double branch_length = full_param[3+i];
				new_ML += mydouble_forward_backward_expectations_ClonalFrame_branch(dec_id,anc_id,node_nuc,position,ipat,kappa,pinuc,branch_length,rho_over_theta,mean_import_length,import_divergence,numEmiss,denEmiss,numTrans,denTrans).LOG();
				// Update estimate of the branch length
				const double mutU_br = numEmiss[0][1];
				const double nsiU_br = denEmiss[0];
				full_param[3+i] = (prior_a[3]+mutU_br)/(prior_b[3]+nsiU_br);
				posterior_a[3+i] = (prior_a[3]+mutU_br);
				// Increment counters for the other expectations
				mutI += numEmiss[1][1];
				nsiI += denEmiss[1];
				const double numI_br = numTrans[0][1];
				const double lenU_br = denTrans[0];
				numI += numI_br;
				lenU += full_param[3+i]*lenU_br;
				numU += numTrans[1][0];
				lenI += denTrans[1];
				if(coutput) {
					cout << "nmut = " << mutU_br << " nU = " << nsiU_br << " nsub = " << numEmiss[1][1] << " nI = " << denEmiss[1] << endl;
					cout << "nU>I = " << numI_br << " dU = " << lenU_br << " nI>U = " << numTrans[1][0] << " dI = " << denTrans[1] << endl;
					cout << "numTrans = " << numTrans[0][0] << " " << numTrans[0][1] << " " << numTrans[1][0] << " " << numTrans[0][0] << endl;
				}
			}
		}
		++neval;
		// Update estimates of the recombination parameters
		full_param[0] = (prior_a[0]+numI)/(prior_b[0]+lenU);
		full_param[1] = (prior_b[1]+lenI)/(prior_a[1]+numU);
		full_param[2] = (prior_a[2]+mutI)/(prior_b[2]+nsiI);
		posterior_a[0] = (prior_a[0]+numI);
		posterior_a[1] = (prior_a[1]+numU);
		posterior_a[2] = (prior_a[2]+mutI);
		if(coutput) {
			cout << "params =";
			for(int j=0;j<full_param.size();j++) cout << " " << full_param[j];
			cout << " ML = " << new_ML << endl;
		}
		// Test for no further improvement
		if(new_ML-ML< -threshold) {
			//cout << "Old likelihood = " << ML << " delta = " << new_ML-ML << endl;
			//warning("Likelihood got worse in Baum_Welch");
		} else if(fabs(new_ML-ML)<threshold) {
			break;
		}
		// Otherwise continue
		ML = new_ML;
	}
	if(it==maxit) warning("Baum_Welch(): maximum number of iterations reached");
	// Once more for debugging purposes
	// mydouble_forward_backward_expectations_ClonalFrame_branch(dec_id,anc_id,node_nuc,position,ipat,kappa,pinuc,branch_length,rho_over_theta,mean_import_length,import_divergence,numEmiss,denEmiss,numTrans,denTrans);
	return ML;
}

double gamma_loglikelihood(const double x, const double a, const double b) {
	return a*log(b)-lgamma(a)+(a-1)*log(x)-b*x;
}

void forward_backward_simulate_expectations_ClonalFrame_branch(const int dec_id, const int anc_id, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const double branch_length, const double rho_over_theta, const double mean_import_length, const double import_divergence, const int nsim, vector<double> &mutU, vector<double> &nsiU, vector<double> &mutI, vector<double> &nsiI, vector<double> &numUI, vector<double> &lenU, vector<double> &numIU, vector<double> &lenI) {
	const int npos = position.size();
	// Define an HKY85 emission probability matrix for Unimported sites
	Matrix<mydouble> pemisUnimported;
	pemisUnimported = compute_HKY85_ptrans(branch_length,kappa,pinuc);
	// Define an HKY85 emission probability matrix for Imported sites
	Matrix<mydouble> pemisImported;
	pemisImported = compute_HKY85_ptrans(import_divergence,kappa,pinuc);
	// Define storage space for the intermediate forward calculations and counters
	Matrix<mydouble> A;
	A = Matrix<mydouble>(npos,2);
	Matrix<double> numEmis, numTrans;
	vector<double> denEmis, denTrans;
	// Define storage space for the observation at every site
	vector<int> emittedState;
	emittedState = vector<int>(npos);
	// Recombination parameters
	const double recrate = rho_over_theta*branch_length;
	const double endrecrate = 1.0/mean_import_length;
	const double totrecrate = recrate+endrecrate;	
	// Transient storage
	mydouble aprev[2];
	mydouble a[2];
	// Equilibrium frequency of unimported and imported sites respectively
	const mydouble pi[2] = {endrecrate/totrecrate,recrate/totrecrate};
	// Beginning at the first variable site, do the forward algorithm
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
		// Store for the second pass
		A[i][0] = a[0];
		A[i][1] = a[1];
		// Store the emitted state for simulation later
		emittedState[i] = (int)(anc!=dec);
	}
	// Second pass: backward algorithm
	// Define storage for the backward simulation probabilities
	Matrix<double> P;
	P = Matrix<double>(npos,2); // P[i][j] is the probability of going from position (i+1) state j to position i state 1
	mydouble bnext[2];
	mydouble b[2];
	// Beginning at the last variable site, do the backward algorithm and calculate backward simulation probabilities
	for(i=npos-1;i>=0;i--) {
		if(i==(npos-1)) {
			// Backward algorithm
			b[0] = mydouble(1.0);
			b[1] = mydouble(1.0);
			// Calculate the backwards simulation probability
			// A[npos-1][j]*b[j] is the joint probability of the data and state j at the final position
			const mydouble num = A[npos-1][1]*b[1];
			const mydouble den = A[npos-1][0]*b[0] + num;
			P[npos-1][0] = P[npos-1][1] = (num/den).todouble();
		} else {
			// Backward algorithm
			bnext[0] = b[0];
			bnext[1] = b[1];
			// Note that these retrieve the ancestral and descendant nucleotides at the 3prime adjacent site
			Nucleotide dec = node_nuc[dec_id][ipat[i+1]];
			Nucleotide anc = node_nuc[anc_id][ipat[i+1]];
			const mydouble pemisU = pemisUnimported[anc][dec];
			const mydouble pemisI = pemisImported[anc][dec];
			mydouble prnotrans;
			prnotrans.setlog(-totrecrate*(position[i+1]-position[i]));
			const mydouble prtrans = mydouble(1.0)-prnotrans;
			const mydouble sumbnext = prtrans*(pi[0]*pemisU*bnext[0] + pi[1]*pemisI*bnext[1]);
			b[0] = prnotrans*pemisU*bnext[0]+sumbnext;
			b[1] = prnotrans*pemisI*bnext[1]+sumbnext;
			// Calculate the backwards simulation probability
			const mydouble pemis[2]  = {pemisU,pemisI};
			// numjk is proportional to the probability of going from state j at position (i+1) to state k at position i
			mydouble num00 = A[i][0]*(prnotrans+prtrans*pi[0])*pemis[0]*bnext[0];
			mydouble num01 = A[i][1]*prtrans*pi[0]*pemis[0]*bnext[0];
			mydouble num10 = A[i][0]*prtrans*pi[1]*pemis[1]*bnext[1];
			mydouble num11 = A[i][1]*(prnotrans+prtrans*pi[1])*pemis[1]*bnext[1];
			P[i][0] = (num01/(num00+num01)).todouble();
			P[i][1] = (num11/(num10+num11)).todouble();
		}
	}
	// Simulate the number of transitions and emissions
	int sim;
	for(sim=0;sim<nsim;sim++) {
		// Resize if necessary and zero the counters
		numEmis = Matrix<double>(2,2,0.0);
		denEmis = vector<double>(2,0.0);
		numTrans = Matrix<double>(2,2,0.0);
		denTrans = vector<double>(2,0.0);
		// Cycle from 3prime to 5prime
		int last;		// Last hidden state
		for(i=npos-1;i>=0;i--) {
			if(i==(npos-1)) {
				// Start by simulating the 3prime-most position
				last = ran.bernoulli(P[i][0]);
				// Update relevant counters
				++numEmis[last][emittedState[i]];
				++denEmis[last];
			} else {
				// Simulate the 5prime-next position
				const int next = ran.bernoulli(P[i][last]);
				// Update all the counters
				++numEmis[next][emittedState[i]];
				++denEmis[next];
				const double dist = position[i+1]-position[i];
				if(dist<=1000.0) {
					++numTrans[next][last];
					denTrans[next] += dist;
				}
				last = next;
			}
		}
		mutU[sim] = numEmis[0][1];
		nsiU[sim] = denEmis[0];
		mutI[sim] = numEmis[1][1];
		nsiI[sim] = denEmis[1];
		numUI[sim] = numTrans[0][1];
		lenU[sim] = denTrans[0];
		numIU[sim] = numTrans[1][0];
		lenI[sim] = denTrans[1];
	}
}

Matrix<double> Baum_Welch_simulate_posterior(const marginal_tree &tree, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const vector<bool> &informative, const vector<double> &prior_a, const vector<double> &prior_b, const vector<double> &full_param, int &neval, const bool coutput, const int nsim) {
	// Storage for output: for each parameter, simulated values
	Matrix<double> post(3,nsim,0.0);
	// Storage for the simulated counts of transitions and emissions
	vector<double> /*mutU(nsim,0.0), nsiU(nsim,0.0),*/ mutI(nsim,0.0), nsiI(nsim,0.0);
	vector<double> numUI(nsim,0.0), lenU(nsim,0.0), numIU(nsim,0.0), lenI(nsim,0.0);
	vector<double> mutU_br(nsim,0.0), nsiU_br(nsim,0.0), mutI_br(nsim,0.0), nsiI_br(nsim,0.0);
	vector<double> numUI_br(nsim,0.0), lenU_br(nsim,0.0), numIU_br(nsim,0.0), lenI_br(nsim,0.0);
	// Estimated parameters
	double rho_over_theta = full_param[0];
	double mean_import_length = full_param[1];
	double import_divergence = full_param[2];
	// Do all the simulations for each branch individually, and combine
	int i;
	for(i=0;i<informative.size();i++) {
		if(informative[i]) {
			const int dec_id = tree.node[i].id;
			const int anc_id = tree.node[i].ancestor->id;
			const double branch_length = full_param[3+i];
			forward_backward_simulate_expectations_ClonalFrame_branch(dec_id,anc_id,node_nuc,position,ipat,kappa,pinuc,branch_length,rho_over_theta,mean_import_length,import_divergence,nsim,mutU_br,nsiU_br,mutI_br,nsiI_br,numUI_br,lenU_br,numIU_br,lenI_br);
			// Update the running totals for each simulation
			int sim;
			for(sim=0;sim<nsim;sim++) {
				// For each branch, necessary to simulate the branch length
				const double a = prior_a[3]+mutU_br[sim];
				const double b = prior_b[3]+nsiU_br[sim];
				const double sim_branch_length = ran.gamma(1.0/b,a);
				numUI[sim] += numUI_br[sim];
				lenU[sim] += sim_branch_length*lenU_br[sim];
				numIU[sim] += numIU_br[sim];
				lenI[sim] += lenI_br[sim];
				mutI[sim] += mutI_br[sim];
				nsiI[sim] += nsiI_br[sim];
			}
		}
	}
	// Simulate the recombination parameters
	int sim;
	for(sim=0;sim<nsim;sim++) {
		// rho over theta
		const double a0 = prior_a[0]+numUI[sim];
		const double b0 = prior_b[0]+lenU[sim];
		post[0][sim] = ran.gamma(1.0/b0,a0);
		// Mean import length
		const double a1 = prior_a[1]+numIU[sim];
		const double b1 = prior_b[1]+lenI[sim];
		post[1][sim] = 1.0/ran.gamma(1.0/b1,a1);
		// Mean import divergence
		const double a2 = prior_a[2]+mutI[sim];
		const double b2 = prior_b[2]+nsiI[sim];
		post[2][sim] = ran.gamma(1.0/b2,a2);
	}
	return post;
}


double Baum_Welch_Rho_Per_Branch(const marginal_tree &tree, const Matrix<Nucleotide> &node_nuc, const vector<double> &position, const vector<int> &ipat, const double kappa, const vector<double> &pinuc, const vector<bool> &informative, const vector<double> &prior_a, const vector<double> &prior_b, vector<double> &mean_param, Matrix<double> &full_param, Matrix<double> &posterior_a, int &neval, const bool coutput) {
	int i;
	if(coutput) cout << setprecision(9);
	// Resize as necessary
	posterior_a = Matrix<double>(informative.size(),4);
	// Storage for the expected number of transitions and emissions in the HMM per branch
	Matrix<double> numEmiss(2,2), numTrans(2,2);
	vector<double> denEmiss(2),   denTrans(2);
	// Counters per branch
	vector<double> mutU_br(informative.size(),0.0), mutI_br(informative.size(),0.0);
	vector<double> nsiU_br(informative.size(),0.0), nsiI_br(informative.size(),0.0);
	vector<double> numI_br(informative.size(),0.0), numU_br(informative.size(),0.0);
	vector<double> lenU_br(informative.size(),0.0), lenI_br(informative.size(),0.0);
	// Calculate the marginal likelihood and expected number of transitions and emissions by the forward-backward algorithm
	// Include the effect of the prior (this is dubious - should instead compute loglikelihood of the pseudocounts)
	double ML = gamma_loglikelihood(mean_param[0], prior_a[0], prior_b[0]) + gamma_loglikelihood(mean_param[1], prior_a[1], prior_b[1]) + gamma_loglikelihood(mean_param[2], prior_a[2], prior_b[2]) + gamma_loglikelihood(mean_param[3], prior_a[3], prior_b[3]);
	for(i=0;i<informative.size();i++) {
		if(informative[i]) {
			// Include the effect of the prior (this is dubious - should instead compute loglikelihood of the pseudocounts)
			ML += gamma_loglikelihood(full_param[i][0], prior_a[4], prior_b[4]) + gamma_loglikelihood(full_param[i][1], prior_a[4], prior_b[4])
			+ gamma_loglikelihood(full_param[i][2], prior_a[4], prior_b[4]) + gamma_loglikelihood(full_param[i][3], prior_a[4], prior_b[4]);
			const int dec_id = tree.node[i].id;
			const int anc_id = tree.node[i].ancestor->id;
			// Initial parameters
			const double rho_over_theta = mean_param[0]*full_param[i][0];
			const double mean_import_length = 1.0/(mean_param[1]*full_param[i][1]);	// NB internal definition
			const double import_divergence = mean_param[2]*full_param[i][2];
			const double branch_length = mean_param[3]*full_param[i][3];
			ML += mydouble_forward_backward_expectations_ClonalFrame_branch(dec_id,anc_id,node_nuc,position,ipat,kappa,pinuc,branch_length,rho_over_theta,mean_import_length,import_divergence,numEmiss,denEmiss,numTrans,denTrans).LOG();
			// Store counters per branch
			mutU_br[i] = numEmiss[0][1];
			nsiU_br[i] = denEmiss[0];
			mutI_br[i] = numEmiss[1][1];
			nsiI_br[i] = denEmiss[1];
			numI_br[i] = numTrans[0][1];
			lenU_br[i] = denTrans[0];
			numU_br[i] = numTrans[1][0];
			lenI_br[i] = denTrans[1];
			//			if(coutput) {
			//				cout << "nmut = " << mutU_br << " nU = " << nsiU_br << " nsub = " << numEmiss[1][1] << " nI = " << denEmiss[1] << endl;
			//				cout << "nU>I = " << numI_br << " dU = " << lenU_br << " nI>U = " << numTrans[1][0] << " dI = " << denTrans[1] << endl;
			//				cout << "numTrans = " << numTrans[0][0] << " " << numTrans[0][1] << " " << numTrans[1][0] << " " << numTrans[0][0] << endl;
			//			}
		}
	}
	++neval;
	// Update estimates of all the parameters: start with the branch lengths
	double mean_param_num, mean_param_den;
	// First, iterate to update the mean branch length parameter (max 3 times)
	int j;
	for(j=0;j<3;j++) {
		mean_param_num = prior_a[3];
		mean_param_den = prior_b[3];
		for(i=0;i<informative.size();i++) {
			if(informative[i]) {
				mean_param_num += mutU_br[i];
				mean_param_den += (prior_a[4]+mutU_br[i])*nsiU_br[i]/(prior_b[4]+mean_param[3]*nsiU_br[i]);
			}
		}
		mean_param[3] = mean_param_num/mean_param_den;
	}
	// Second, update the individual branch length factors
	for(i=0;i<informative.size();i++) {
		if(informative[i]) {
			full_param[i][3] = (prior_a[4]+mutU_br[i])/(prior_b[4]+mean_param[3]*nsiU_br[i]);
			posterior_a[i][3] = (prior_a[4]+mutU_br[i]);
		}
	}
	// Next update the recombination parameters
	int p;
	for(p=0;p<3;p++) {
		// First, iterate to update the mean branch length parameter (max 3 times)
		for(j=0;j<3;j++) {
			mean_param_num = prior_a[p];
			mean_param_den = prior_b[p];
			for(i=0;i<informative.size();i++) {
				if(informative[i]) {
					double num, den;
					if(p==0) {
						num = numI_br[i];
						// Adjust the counts for R/M according to the estimated branch lengths
						den = lenU_br[i]*mean_param[3]*full_param[i][3];
					} else if(p==1) {
						num = numU_br[i];
						den = lenI_br[i];
					} else {
						num = mutI_br[i];
						den = nsiI_br[i];
					}
					mean_param_num += num;
					mean_param_den += (prior_a[4]+num)*den/(prior_b[4]+mean_param[p]*den);
				}
			}
			mean_param[p] = mean_param_num/mean_param_den;
		}
		// Second, update the individual per branch factors
		for(i=0;i<informative.size();i++) {
			if(informative[i]) {
				double num, den;
				if(p==0) {
					num = numI_br[i];
					// Adjust the counts for R/M according to the estimated branch lengths
					den = lenU_br[i]*mean_param[3]*full_param[i][3];
				} else if(p==1) {
					num = numU_br[i];
					den = lenI_br[i];
				} else {
					num = mutI_br[i];
					den = nsiI_br[i];
				}
				full_param[i][p] = (prior_a[4]+num)/(prior_b[4]+mean_param[p]*den);
				posterior_a[i][p] = (prior_a[4]+num);
			}
		}
	}
	if(coutput) {
		cout << "mean params =";
		for(int j=0;j<mean_param.size();j++) cout << " " << mean_param[j];
		cout << " ML = " << ML << endl;
	}
	
	// Iterate until the maximum likelihood improves by less than some threshold
	const int maxit = 200;
	const double threshold = 1.0e-2;
	int it;
	for(it=0;it<maxit;it++) {
		// Update the likelihood
		double new_ML = gamma_loglikelihood(mean_param[0], prior_a[0], prior_b[0]) + gamma_loglikelihood(mean_param[1], prior_a[1], prior_b[1]) + gamma_loglikelihood(mean_param[2], prior_a[2], prior_b[2]) + gamma_loglikelihood(mean_param[3], prior_a[3], prior_b[3]);
		for(i=0;i<informative.size();i++) {
			if(informative[i]) {
				// Include the effect of the prior (this is dubious - should instead compute loglikelihood of the pseudocounts)
				new_ML += gamma_loglikelihood(full_param[i][0], prior_a[4], prior_b[4]) + gamma_loglikelihood(full_param[i][1], prior_a[4], prior_b[4])
				+ gamma_loglikelihood(full_param[i][2], prior_a[4], prior_b[4]) + gamma_loglikelihood(full_param[i][3], prior_a[4], prior_b[4]);
				const int dec_id = tree.node[i].id;
				const int anc_id = tree.node[i].ancestor->id;
				// Initial parameters
				const double rho_over_theta = mean_param[0]*full_param[i][0];
				const double mean_import_length = 1.0/(mean_param[1]*full_param[i][1]);	// NB internal definition
				const double import_divergence = mean_param[2]*full_param[i][2];
				const double branch_length = mean_param[3]*full_param[i][3];
				new_ML += mydouble_forward_backward_expectations_ClonalFrame_branch(dec_id,anc_id,node_nuc,position,ipat,kappa,pinuc,branch_length,rho_over_theta,mean_import_length,import_divergence,numEmiss,denEmiss,numTrans,denTrans).LOG();
				// Store counters per branch
				mutU_br[i] = numEmiss[0][1];
				nsiU_br[i] = denEmiss[0];
				mutI_br[i] = numEmiss[1][1];
				nsiI_br[i] = denEmiss[1];
				numI_br[i] = numTrans[0][1];
				lenU_br[i] = denTrans[0];
				numU_br[i] = numTrans[1][0];
				lenI_br[i] = denTrans[1];
			}
		}
		++neval;
		// Update estimates of all the parameters: start with the branch lengths
		// First, iterate to update the mean branch length parameter (max 3 times)
		int j;
		for(j=0;j<3;j++) {
			mean_param_num = prior_a[3];
			mean_param_den = prior_b[3];
			for(i=0;i<informative.size();i++) {
				if(informative[i]) {
					mean_param_num += mutU_br[i];
					mean_param_den += (prior_a[4]+mutU_br[i])*nsiU_br[i]/(prior_b[4]+mean_param[3]*nsiU_br[i]);
				}
			}
			mean_param[3] = mean_param_num/mean_param_den;
		}
		// Second, update the individual branch length factors
		for(i=0;i<informative.size();i++) {
			if(informative[i]) {
				full_param[i][3] = (prior_a[4]+mutU_br[i])/(prior_b[4]+mean_param[3]*nsiU_br[i]);
				posterior_a[i][3] = (prior_a[4]+mutU_br[i]);
			}
		}
		// Next update the recombination parameters
		for(p=0;p<3;p++) {
			// First, iterate to update the mean branch length parameter (max 3 times)
			for(j=0;j<3;j++) {
				mean_param_num = prior_a[p];
				mean_param_den = prior_b[p];
				for(i=0;i<informative.size();i++) {
					if(informative[i]) {
						double num, den;
						if(p==0) {
							num = numI_br[i];
							// Adjust the counts for R/M according to the estimated branch lengths
							den = lenU_br[i]*mean_param[3]*full_param[i][3];
						} else if(p==1) {
							num = numU_br[i];
							den = lenI_br[i];
						} else {
							num = mutI_br[i];
							den = nsiI_br[i];
						}
						mean_param_num += num;
						mean_param_den += (prior_a[4]+num)*den/(prior_b[4]+mean_param[p]*den);
					}
				}
				mean_param[p] = mean_param_num/mean_param_den;
			}
			// Second, update the individual branch length factors
			for(i=0;i<informative.size();i++) {
				if(informative[i]) {
					double num, den;
					if(p==0) {
						num = numI_br[i];
						// Adjust the counts for R/M according to the estimated branch lengths
						den = lenU_br[i]*mean_param[3]*full_param[i][3];
					} else if(p==1) {
						num = numU_br[i];
						den = lenI_br[i];
					} else {
						num = mutI_br[i];
						den = nsiI_br[i];
					}
					full_param[i][p] = (prior_a[4]+num)/(prior_b[4]+mean_param[p]*den);
					posterior_a[i][p] = (prior_a[4]+num);
				}
			}
		}
		if(coutput) {
			cout << "mean params =";
			for(int j=0;j<mean_param.size();j++) cout << " " << mean_param[j];
			cout << " ML = " << new_ML << endl;
		}
		// Test for no further improvement
		if(new_ML-ML< -threshold) {
			//cout << "Old likelihood = " << ML << " delta = " << new_ML-ML << endl;
			//warning("Likelihood got worse in Baum_Welch");
		} else if(fabs(new_ML-ML)<threshold) {
			break;
		}
		// Otherwise continue
		ML = new_ML;
	}
	if(it==maxit) warning("Baum_Welch_Rho_Per_Branch(): maximum number of iterations reached");
	return ML;
}
