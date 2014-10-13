ClonalFrameML
Xavier Didelot and Daniel Wilson. 2013

This program reads in a Newick tree and FASTA file and, for all variable sites, reconstructs
the joint maximum likelihood sequences at all nodes (including, for the purposes of imputation, tips)
using the HKY85 nucleotide substitution model and an algorithm described in:

    A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences
    Tal Pupko, Itsik Peer, Ron Shamir, and Dan Graur. Mol. Biol. Evol. 17(6):890–896. 2000
	
Branch lengths of the tree are corrected for heterospecific horizontal gene transfer using a new maximum-
likelihood algorithm implementing the ClonalFrame model that was described in:

    Inference of Bacterial Microevolution Using Multilocus Sequence Data
	Xavier Didelot, and Daniel Falush. Genetics 175(3):1251-1266. 2007

Syntax: ClonalFrameML newick_file fasta_file output_file [OPTIONS]

newick_file  The tree specified in Newick format. It must be an unrooted bifurcating tree. All
             tips should be uniquely labelled and the internal nodes must not be labelled. Note that the
             branch lengths must be scaled in units of expected number of substitutions per site.
             Failure to provide appropriately scaled branch lengths will adversely affect results.
fasta_file   The nucleotide sequences specified in FASTA format, with labels exactly matching those in
			 the newick_file. The letter codes A, C, G and T are interpreted directly, U is converted
			 to T, and N, -, ? and X are treated equivalently as ambiguity codes. No other codes are
			 allowed.
output_file  The prefix for the output files, described below.
[OPTIONS]    Run ClonalFrameML with no arguments to see the options available.

The program reports the empirical nucleotide frequencies and the joint log-likelihood of the reconstructed
sequences for variable sites. Files are output with the following suffixes:

.labelled_tree.newick          The corrected Newick tree is ouput with internal nodes labelled so that they
                               correspond with the reconstructed ancestral sequence file.
.ML_sequence.fasta             The reconstructed sequences (ancestral and, for the purposes of imputation,
							   observed) in FASTA format with letter codes A, C, G and T only. The labels
							   match exactly those in the output Newick tree.
.position_cross_reference.txt  A vector of comma-separated values equal in length to the input FASTA file
							   relating the positions of (variable) sites in the input FASTA file to the
							   positions of their reconstructed sequences in the output FASTA file, starting
							   with position 1. Sites in the input file not reconstructed are assigned a 0.
.importation_status.txt        A FASTA file representing the inferred importation status of every site
                               coded as 0 (unimported) 1 (imported) 2 (unimported homoplasy/multiallelic)
							   3 (imported homoplasy/multiallelic) 4 (untested compatible) 5 (untested
							   homoplasy).
