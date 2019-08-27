# ClonalFrameML

# Introduction #

This is the homepage of ClonalFrameML, a software package that performs efficient inference of recombination in bacterial genomes. ClonalFrameML was created by [Xavier Didelot](http://xavierdidelot.github.io) and [Daniel Wilson](http://www.danielwilson.me.uk/). ClonalFrameML can be applied to any type of aligned sequence data, but is especially aimed at analysis of whole genome sequences. It is able to compare hundreds of whole genomes in a matter of hours on a standard Desktop computer. There are three main outputs from a run of ClonalFrameML: a phylogeny with branch lengths corrected to account for recombination, an estimation of the key parameters of the recombination process, and a genomic map of where recombination took place for each branch of the phylogeny.

ClonalFrameML is a maximum likelihood implementation of the Bayesian software [ClonalFrame](http://xavierdidelot.github.io/clonalframe.html) which was previously described by [Didelot and Falush (2007)](http://www.genetics.org/cgi/content/abstract/175/3/1251). The recombination model underpinning ClonalFrameML is exactly the same as for ClonalFrame, but this new implementation is a lot faster, is able to deal with much larger genomic dataset, and does not suffer from MCMC convergence issues. A scientific paper describing ClonalFrameML in detail has been published, see [Didelot X, Wilson DJ (2015) ClonalFrameML: Efficient Inference of Recombination in Whole Bacterial Genomes. PLoS Comput Biol 11(2): e1004041. doi:10.1371/journal.pcbi.1004041](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004041).

# Download and Installation #

You can obtain the most up to date version of ClonalFrameML by downloading and compiling the C++ source code via GIT using the command:
```
git clone https://github.com/xavierdidelot/ClonalFrameML
```

Please note that the code for ClonalFrameML is distributed under the terms of the GNU GPL v3 license, for more details see https://www.gnu.org/copyleft/gpl.html

You can compile the code using the following command:

```
cd ClonalFrameML/src
./make.sh
```

Compilation requires a C++ compiler, such as [GCC](https://gcc.gnu.org/), to be installed. Running the bundled R scripts requires [R](http://cran.r-project.org/) to be installed with the ape and phangorn packages.

# User guide #

The user guide for ClonalFrameML is available [here](https://github.com/xavierdidelot/clonalframeml/wiki).

# Getting help #

If you need assistance using ClonalFrameML, you can get in touch by emailing either [Xavier Didelot](http://xavierdidelot.github.io/contact.html) or [Daniel Wilson](http://www.danielwilson.me.uk/contact.html).
