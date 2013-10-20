#include "main.h"

int main (const int argc, const char* argv[]) {
	clock_t start_time = clock();
	// Process the command line arguments	
	if(argc<5) {
		stringstream errTxt;
		errTxt << "Syntax: ClonalFrameML newick_file fasta_file kappa output_file [OPTIONS]" << endl;
		errTxt << endl;
		errTxt << "The following options are available:" << endl;
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
		errTxt << "-rescale_no_recombination      true or false (default)   Rescale branch lengths for given sites with no recombination model." << endl;		
		errTxt << "-multithread                   true or false (default)   Enable OpenMP parallel code. Overhead may cancel out gains." << endl;
		errTxt << "-show_progress                 true or false (default)   Output the progress of the maximum likelihood routines." << endl;
		errTxt << "-compress_reconstructed_sites  true (default) or false   Reduce the number of columns in the output FASTA file." << endl;
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
	// Set default options
	ArgumentWizard arg;
	arg.case_sensitive = false;
	string correct_branch_lengths="true", excess_divergence_model="false", ignore_incomplete_sites="false", ignore_user_sites="", reconstruct_invariant_sites="false";
	string use_incompatible_sites="false", joint_branch_param="false", rho_per_branch="false", rescale_no_recombination="false", multithread="false", show_progress="false", compress_reconstructed_sites="true";
	double brent_tolerance = 1.0e-3, powell_tolerance = 1.0e-3;
	// Process options
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
	arg.add_item("rescale_no_recombination",	TP_STRING, &rescale_no_recombination);
	arg.add_item("multithread",			        TP_STRING, &multithread);
	arg.add_item("show_progress",				TP_STRING, &show_progress);
	arg.add_item("compress_reconstructed_sites",TP_STRING, &compress_reconstructed_sites);	
	arg.read_input(argc-4,argv+4);
	bool CORRECT_BRANCH_LENGTHS			= string_to_bool(correct_branch_lengths,		"correct_branch_lengths");
	bool EXCESS_DIVERGENCE_MODEL		= string_to_bool(excess_divergence_model,		"excess_divergence_model");
	bool IGNORE_INCOMPLETE_SITES		= string_to_bool(ignore_incomplete_sites,		"ignore_incomplete_sites");
	bool RECONSTRUCT_INVARIANT_SITES	= string_to_bool(reconstruct_invariant_sites,	"reconstruct_invariant_sites");
	bool USE_INCOMPATIBLE_SITES			= string_to_bool(use_incompatible_sites,		"use_incompatible_sites");
	bool JOINT_BRANCH_PARAM				= string_to_bool(joint_branch_param,			"joint_branch_param");
	bool RHO_PER_BRANCH					= string_to_bool(rho_per_branch,				"rho_per_branch");
	bool RESCALE_NO_RECOMBINATION		= string_to_bool(rescale_no_recombination,		"rescale_no_recombination");
	bool MULTITHREAD					= string_to_bool(multithread,					"multithread");
	bool SHOW_PROGRESS					= string_to_bool(show_progress,					"show_progress");
	bool COMPRESS_RECONSTRUCTED_SITES	= string_to_bool(compress_reconstructed_sites,	"compress_reconstructed_sites");
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
	if(((int)JOINT_BRANCH_PARAM + (int)RHO_PER_BRANCH + (int)RESCALE_NO_RECOMBINATION)>1) {
		stringstream errTxt;
		errTxt << "joint_branch_param, rho_per_branch and rescale_no_recombination are mutually incompatible";
		error(errTxt.str().c_str());
	}
	if((EXCESS_DIVERGENCE_MODEL || JOINT_BRANCH_PARAM || RHO_PER_BRANCH || RESCALE_NO_RECOMBINATION) && !CORRECT_BRANCH_LENGTHS) {
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
	if(MULTITHREAD) {
		cout << "WARNING: multithreaded version not implemented, ignoring." << endl;
	}
	
	// Open the FASTA file
	DNA fa(fasta_file);
	cout << "Read " << fa.nseq << " sequences of length " << fa.lseq << " sites from " << fasta_file << endl;
	if(fa.nseq==2) {
		if(!EXCESS_DIVERGENCE_MODEL) cout << "WARNING: with only two sequences, the excess divergence model is mandatory." << endl;
		EXCESS_DIVERGENCE_MODEL = true;
	}
	// Open the Newick file and convert to internal rooted tree format, outputting the names of the tips and internal nodes
	NewickTree newick = read_Newick(newick_file);
	vector<string> ctree_node_labels;
	marginal_tree ctree = convert_unrooted_NewickTree_to_marginal_tree(newick,fa.label,ctree_node_labels);
	// Open the list of sites to ignore
	vector<bool> ignore_site(fa.lseq,false);
	if(ignore_user_sites!="") {
		ifstream user_sites(ignore_user_sites.c_str());
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
		}
	}
	
	// Compute compatibility and test every site for any sequences with 'N','-','X' or '?'
	// Key to results: -1: invariant, 0: compatible biallelic (including singletons), 1: incompatible biallelic, 2: more than two alleles
	vector<bool> anyN;
	vector<int> compat = compute_compatibility(fa,ctree,anyN,false);
	int i;
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
		
		if(JOINT_BRANCH_PARAM) {
			// Prepare to correct branch lengths
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) under the ClonalFrame model
			const double initial_rho_over_theta = 0.1;
			const double initial_mean_import_length = 500.0;
			const double initial_import_divergence = 0.01;
			ClonalFrameFunction cff(ctree,node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,EXCESS_DIVERGENCE_MODEL,MULTITHREAD);
			vector<double> param(3+ctree.size-2);
			param[0] = log10(initial_rho_over_theta); param[1] = log10(initial_mean_import_length); param[2] = log10(initial_import_divergence);
			for(i=0;i<ctree.size-2;i++) param[3+i] = log10(ctree.node[i].edge_time); 
			
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
			for(i=0;i<ctree.size-2;i++) {
				cout << "Branch length " << ctree_node_labels[i] << " = " << pow(10.,param[3+i]) << endl;
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = pow(10.,param[3+i]);
			}

			// Output the importation status
			write_importation_status(cff.is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str());	
			cout << "Wrote inferred importation status to " << import_out_file << endl;
		} else if(RHO_PER_BRANCH) {
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
			vector< vector<ImportationState> > is_imported(ctree.size-2);
			for(i=0;i<ctree.size-2;i++) {
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
				const double initial_rho_over_theta = 0.1;
				const double initial_import_ratio = 0.1;			// constrained to be less than 1
				const double initial_import_divergence = 0.1;		// multiplicative excess (excess model is mandatory), so subst_rate = mut_rate * (2+import_divergence)
				// Minimum branch length
				const double min_branch_length = 1.0e-7;
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
				const double final_branch_length = pow(10.,param[3]);
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
			write_importation_status_intervals(is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str());	
			cout << "Wrote inferred importation status to " << import_out_file << endl;
		} else if(RESCALE_NO_RECOMBINATION) {
			// Rescale the branch lengths using given sites without a model of recombination
			cout << "Beginning branch optimization. Key to parameters (and constraints):" << endl;
			cout << "B   uncorrected branch length" << endl;
			cout << "L   maximum log-likelihood per branch" << endl;
			cout << "M   corrected branch length/expected number of mutations     (> 0)" << endl;
			double ML = 0.0;
			for(i=0;i<ctree.size-2;i++) {
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
				const double min_branch_length = 1.0e-7;
				ClonalFrameRescaleBranchFunction cff(ctree.node[i],node_nuc,pat1,cpat,kappa,empirical_nucleotide_frequencies,MULTITHREAD,initial_branch_length,min_branch_length);
				// Setup optimization function
				Powell Pow(cff);
				Pow.coutput = Pow.brent.coutput = SHOW_PROGRESS;
				Pow.TOL = brent_tolerance;
				// Estimate parameter
				vector<double> param(1,log10(initial_branch_length));
				param = Pow.minimize(param,powell_tolerance);
				const double final_branch_length = pow(10.,param[0]);
				// Update branch length in the tree
				// Note this is unsafe in general because the corresponding node times are not adjusted
				ctree.node[i].edge_time = final_branch_length;
				cout << "Branch " << ctree_node_labels[i] << " B = " << initial_branch_length << " L = " << -Pow.function_minimum << " M = " << final_branch_length << endl;
				ML += -Pow.function_minimum;
			}
			cout << "Log-likelihood after branch optimization is " << ML << endl;
		} else {
			// Estimate parameters with branch lengths fixed
			// For a given branch, compute the maximum likelihood importation state (unimported vs imported) under the ClonalFrame model
			const double initial_rho_over_theta = 0.1;
			const double initial_mean_import_length = 500.0;
			const double initial_import_divergence = 0.01;
			ClonalFrameParameterFunction cff(ctree,node_nuc,isBLC,ipat,kappa,empirical_nucleotide_frequencies,EXCESS_DIVERGENCE_MODEL,MULTITHREAD);
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
			for(i=0;i<ctree.size-2;i++) {
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
			write_importation_status(cff.is_imported,ctree_node_labels,isBLC,compat,import_out_file.c_str());	
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
	int i,npos=0;
	for(i=0;i<iscompat.size();i++) {
		if(iscompat[i]) {
			position.push_back((double)i);
			++npos;
		}
	}
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

void write_importation_status(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, const char* file_name) {
	ofstream fout(file_name);
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_importation_status(): could not open file " << file_name << " for writing";
		error(errTxt.str().c_str());
	}
	write_importation_status(imported,all_node_names,isBLC,compat,fout);
	fout.close();
}

void write_importation_status(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, ofstream &fout) {
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_importation_status(): could not open file stream for writing";
		error(errTxt.str().c_str());
	}
	if(imported.size()!=all_node_names.size()-2) {
		stringstream errTxt;
		errTxt << "write_importation_status(): number of lineages (" << imported.size() << ") does not equal the number of node labels (" << all_node_names.size() << ") minus two";
		error(errTxt.str().c_str());
	}
	int i,pos;
	for(i=0;i<imported.size();i++) {
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

void write_importation_status_intervals(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, const char* file_name) {
	ofstream fout(file_name);
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_importation_status_intervals(): could not open file " << file_name << " for writing";
		error(errTxt.str().c_str());
	}
	write_importation_status_intervals(imported,all_node_names,isBLC,compat,fout);
	fout.close();
}

void write_importation_status_intervals(vector< vector<ImportationState> > &imported, vector<string> &all_node_names, vector<bool> &isBLC, vector<int> &compat, ofstream &fout) {
	if(!fout) {
		stringstream errTxt;
		errTxt << "write_importation_status_intervals(): could not open file stream for writing";
		error(errTxt.str().c_str());
	}
	if(imported.size()!=all_node_names.size()-2) {
		stringstream errTxt;
		errTxt << "write_importation_status_intervals(): number of lineages (" << imported.size() << ") does not equal the number of node labels (" << all_node_names.size() << ") minus two";
		error(errTxt.str().c_str());
	}
	const char tab = '\t';
	fout << "Node" << tab << "Beg" << tab << "End" << endl;
	int i,pos;
	for(i=0;i<imported.size();i++) {
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



