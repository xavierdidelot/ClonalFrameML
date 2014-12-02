# ClonalFrameML results
# Planned usage: Rscript cfml_results.R ...
help = paste(
"cfml_results.R summarizes the results of a ClonalFrameML analysis",
"Daniel Wilson (2014)",
"",
"Usage: Rscript cfml_results.R prefix [coresites_list]",
sep="\n")

# Preliminaries
library(ape)
library(phangorn)
source("fasta.R")
# Alternative method of plotting using lines. Assume m>0 is interesting
alt.image = function(m,col=heat.colors(1+max(m,na.rm=TRUE)),xpos=NULL,ypos=NULL,length=1,background.fun=NULL,...) {
	if(is.null(xpos)) xpos = 1:nrow(m)
	if(is.null(ypos)) ypos = 1:ncol(m)
	x = matrix(rep(xpos,ncol(m)),nrow=nrow(m))
	y = matrix(rep(ypos,each=nrow(m)),ncol=ncol(m))
	plot(range(xpos),range(ypos)+c(-length,length)/2,type="n",...)
	rect(min(xpos),min(ypos)-length/2,max(xpos),max(ypos)+length/2,col=col[1],border="NA")
	if(!is.null(background.fun)) background.fun()
	gd = m>0
	COL = matrix(col[1+m],nrow=nrow(m))
	arrows(x[gd],y[gd]-length/2,x[gd],y[gd]+length/2,col=COL[gd],len=0)
}

# Read options from command line
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=1 & length(args)!=2) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")
}
prefix = args[1]
coresites_list = ifelse(length(args)==2,args[2],NA)

if(!is.na(coresites_list)) {
	coresites = scan(coresites_list)
} else {
	coresites = NA
}

# Automatically set
treefile = paste(prefix,".labelled_tree.newick",sep="")
xreffile = paste(prefix,".position_cross_reference.txt",sep="")
ML_seqfile = paste(prefix,".ML_sequence.fasta",sep="")
istatefile = paste(prefix,".importation_status.txt",sep="")
if(!file.exists(istatefile)) istatefile = NA

# Load the phyML tree estimated from all core variant and invariant sites
#tree0 = read.tree(treefile); tree = midpoint(tree0); tree$node.label = c(tree$node.label,setdiff(tree0$node.label,tree$node.label))
tree = read.tree(treefile)

# Load a list cross-referencing patterns in the original data to the output FASTA file
xref = scan(xreffile,sep=",")
genome_length = length(xref)
if(is.na(coresites_list)) {
	coresites = 1:genome_length
} else if(any(coresites>genome_length)) stop("Core site ",which(coresites>genome_length)[1]," exceeds genome length ",genome_length)
if(any(coresites<1)) stop("Core sites must be positive")

# Load the imputed and reconstructed ancestral sequences
ML_seq=scan(ML_seqfile,what=character(0))
tp = substr(ML_seq[seq(1,length(ML_seq),by=2)],2,1000)
ML_seq = ML_seq[seq(2,length(ML_seq),by=2)]; names(ML_seq) = tp
# M is a matrix containing the FASTA file base calls
M = matrix("",length(ML_seq),nchar(ML_seq[1]))
for(i in 1:length(ML_seq)) {
	v = unlist(strsplit(ML_seq[i],""))
	M[i,] = v
	gc()
}
rownames(M) = names(ML_seq)

# Precompute various mappings
# Combine the tip and node labels
treelabels = c(tree$tip.label,tree$node.label)
# For each row of M, identify the node index
M_node_index = match(rownames(M),treelabels)
# And the reverse operation
rev_M_node_index = match(treelabels,rownames(M))
# For each row of M, identify the node index of its ancestor
# To do this, identify the node index in tree$edge[,2] and read tree$edge[,1]
M_anc_node_index = tree$edge[match(M_node_index,tree$edge[,2]),1]
# Find, by name, the ancestor
M_anc_node = treelabels[M_anc_node_index]
# Find its position in M
M_anc_node_M_index = match(M_anc_node,rownames(M))
# Not-root
nonroot = !is.na(M_anc_node_index)
# Map edge order on to M order, and vice versa
edge2M = match(tree$edge[,2],M_node_index)
M2edge = match(M_node_index,tree$edge[,2])

# Precompute the positions of mutations on branches of the tree
# For each pattern, record the mutated nodes
# wh.mut is a matrix, in the same order as M, recording whether the base represents a mutation
wh.mut = apply(M,2,function(m) 1*(m!=m[M_anc_node_M_index])); wh.mut[nrow(wh.mut),] = 0
# Weight of each pattern
wpat = as.vector(table(factor(xref,levels=1:ncol(M))))
# For each node, what proportion of mutations are shared with each other node?
#tp = sapply(1:nrow(wh.mut),function(i) apply(t(t(wh.mut[,wh.mut[i,]==1,drop=FALSE])*wpat[wh.mut[i,]==1]),1,sum)/sum(wpat[wh.mut[i,]==1]))

# A homoplasy is a mutation that occurs on multiple branches. Count the number of homoplasic mutations per branch
# Exclude reference sequences from the count
gd = !is.na(as.numeric(rownames(wh.mut))) | substr(rownames(wh.mut),1,4)=="NODE"
n.mut = apply(wh.mut[gd,],2,sum)
is.homoplasy = n.mut>1
is.core = !is.na(match(1:genome_length,coresites))

# A homoplasy is a mutation that occurs on multiple branches. Count the number of homoplasic mutations per branch
# Exclude reference sequences from the count
#gd = !is.na(as.numeric(rownames(wh.mut))) | substr(rownames(wh.mut),1,4)=="NODE"
#plot.mut = t(wh.mut[,xref[xref>0]])*(1+is.homoplasy[xref[xref>0]])
spectrum.mut = t(wh.mut[,xref[xref>0]])*(n.mut[xref[xref>0]])

# Identify contiguous non-core regions
noncore.beg = 1+which(is.core[2:length(is.core)]==0 & (is.core[2:length(is.core)]!=is.core[1:(length(is.core)-1)])); if(!is.core[1]) noncore.beg = c(1,noncore.beg)
noncore.end = which(is.core[2:length(is.core)]==1 & (is.core[2:length(is.core)]!=is.core[1:(length(is.core)-1)])); if(!is.core[length(is.core)]) noncore.end = c(noncore.end,length(is.core))
noncore.len = noncore.end-noncore.beg+1
noncore.plot = noncore.len>=1000

# Plot "raw" mutations/homoplasies
#f = function() rect(noncore.beg[noncore.plot],0,noncore.end[noncore.plot],ncol(wh.mut),col="white",border=NA)
#noncore.plot = noncore.len>=1000
#alt.image(plot.mut,col=c("skyblue","yellow","yellow"),xlab="Position",ylab="Branch",axes=FALSE,xaxs="i",yaxs="i",xpos=which(xref>0),background.fun=f)
#axis(1); axis(2,1:nrow(wh.mut),rownames(wh.mut),las=2,cex.axis=.4); box()
# Plot the recombination intervals
#ypos = match(itv2$Node,rownames(wh.mut))
#arrows(itv2$Beg,ypos,itv2$End,ypos,len=0,lwd=2,col="blue",lend=2)

if(!is.na(istatefile)) itv2 = read.table(istatefile,h=T,as.is=T,sep="\t")

if(FALSE){
# Histogram of recombination tract lengths
tlen = itv2$End-itv2$Beg
# Identify ones that straddle the original
wh = which(itv2$End==genome_length)
for(i in wh) {
	if(any(itv2$Beg[itv2$Node==itv2$Node[i]]==1)) {
		wh2 = which(itv2$Beg[itv2$Node==itv2$Node[i]]==1)
		tlen[i] = tlen[i]+tlen[itv2$Node==itv2$Node[i]][wh2]
		tlen[itv2$Node==itv2$Node[i]][wh2] = NA
	}
}

hist(tlen,100,col="orange3",prob=T)
hist(log10(tlen),100,col="orange3",prob=T)
plot.ecdf(log10(tlen),col="orange3")
}

# Make all branch lengths equal
#tree.bkp = tree
#tree$edge.length = rep(1,length(tree$edge.length))

tree$comid = ifelse(is.na(as.numeric(tree$tip.label)),tree$tip.label,paste("C0000",as.numeric(tree$tip.label),sep=""))
wh.mlst = ifelse(is.na(as.numeric(rownames(wh.mut))),rownames(wh.mut),paste("C0000",as.numeric(rownames(wh.mut)),sep=""))
#wh.mlst_or_ref = ifelse(is.na(as.numeric(rownames(wh.mut))),rownames(wh.mut),mlst[paste(">",ifelse(is.na(as.numeric(rownames(wh.mut))),rownames(wh.mut),paste("C0000",as.numeric(rownames(wh.mut)),sep="")),"_n1",sep="")]); wh.mlst_or_ref[(1+ceiling(nrow(wh.mut)/2)):nrow(wh.mut)] = ""

pdf(file="/dev/null",width=14,height=7)
par(mfrow=c(1,2))
plot(tree,type="phylogram")
dev.off()
# Based on the phylogram tree plot, find the vertical positions and horizontal end-points of every branch
vpos = get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[M_node_index]
lpos = get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[M_anc_node_index]
rpos = get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[M_node_index]

# Manipulate the vertical positions
new_plot.phylo = get("last_plot.phylo", envir = .PlotPhyloEnv); new_plot.phylo$yy = rank(vpos)[rev_M_node_index]
assign("last_plot.phylo",new_plot.phylo,envir=.PlotPhyloEnv)
vpos = get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[M_node_index]

pdf(file=paste0(prefix,".cfml.pdf"),width=14,height=7)
par(mfrow=c(1,2))
xrg = range(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)
plot(xrg+c(0,diff(xrg)/20),range(get("last_plot.phylo", envir = .PlotPhyloEnv)$yy)+c(-0.5,0.5),type="n",axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
# Plot the horizontal branches
arrows(lpos,vpos,rpos,vpos,col=1,len=0)
# Plot the vertical branches
sapply(sort(union(M_anc_node_index,c())),function(i) {
	vpos = get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[M_node_index[!is.na(M_anc_node_index) & M_anc_node_index==i]]
	hpos = get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[M_node_index[!is.na(M_node_index) & M_node_index==i]]
	if(length(vpos)==2)	arrows(hpos,vpos[1],hpos,vpos[2],len=0,col=1)
})
#text(max(get("last_plot.phylo", envir = .PlotPhyloEnv)$xx)+diff(xrg)/20/2,get("last_plot.phylo", envir = .PlotPhyloEnv)$yy,wh.mlst_or_ref[rev_M_node_index],cex=.4)
# Draw lines from the nodes
arrows(rpos,vpos,rep(xrg[2],length(vpos)),vpos,lty=2,len=0,col="grey")

# Plot "raw" mutations/homoplasies
od = order(vpos)
if(length(noncore.beg)>0) background.noncore = function() rect(noncore.beg[noncore.plot],0,noncore.end[noncore.plot],ncol(wh.mut),col="grey",border=NA) else background.noncore = function() {}
noncore.plot = noncore.len>=10000
alt.image(spectrum.mut[,od],col=c("skyblue","white","yellow",colorRampPalette(c("orange","red"))(pmax(0,max(spectrum.mut)-2))),xlab="Position",ylab="Branch",axes=FALSE,xaxs="i",yaxs="i",xpos=which(xref>0),background.fun=background.noncore)
axis(1); axis(2,1:nrow(wh.mut),ifelse((1:nrow(wh.mut))<=ceiling(nrow(wh.mut)/2),rownames(wh.mut),"")[od],las=2,cex.axis=.4); box()
# Plot the recombination intervals
if(!is.na(istatefile)) {
	ypos = match(itv2$Node,rownames(wh.mut)[od])
	arrows(itv2$Beg,ypos,itv2$End,ypos,len=0,lwd=2,col="blue",lend=2)
}
dev.off()
