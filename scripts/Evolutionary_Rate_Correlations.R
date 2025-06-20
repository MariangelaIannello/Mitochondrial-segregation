#!/usr/bin/env Rscript

'
Evolutionary_rate_correlations.R

Usage: 
	Evolutionary_Rate_Correlations.R (--terminalbranches | --root2tip) --norm=<norm_tree> --folder=<og_folder> --output=<out_name> --reftree=<tree>

Options: 
	-h --help                       Show this screen

	--terminalbranches              Calculate the correlations between the terminal branches of the trees
	--root2tip                      Calculate the correlations between the root-to-tip distances of each OTU

Arguments:

	--norm=<normtree>               Tree for correction of phylogenetic singnal before fitting a linear model between the tested branch length distributions (Phylogenetic Independet Contrasts)

	--tree_folder=<og_folder>       Folder containing single files for each tree for which to test the correlations against the reference tree

	--reftree=<tree>                Tree to test correlations against

	-o --output=<out_name>          Name for tabular output ("OG","rho","pvalue","BH_fdr","BY_fdr","Hommel_cpvalue","Holm_cpvalue")
' -> doc

library(docopt)
arguments <- docopt(doc)


## Directories
suppressWarnings(dir.create("ERC_analyses"))

## Libraries
library(ape) # read.tree / drop.tip / pic
suppressMessages(library(adephylo)) # distRoot (option "root2tip")


#### ERC ####

## Create empty table to be filled with ERC values
columns <- c("OG","rho","pvalue","BH_fdr","BY_fdr","Hommel_cpvalue","Holm_cpvalue")
erc <- data.frame(matrix(nrow=0,ncol=length(columns)))
colnames(erc) = columns

## Read reference and normalization trees

ref <- read.tree(arguments$reftree)

norm <- read.tree(arguments$norm)

## Perform Evolutionary Rate Correlations

for (i in list.files(path=arguments$folder, full.names=FALSE, recursive=FALSE)) {
	
	# Input the trees and prune them to have the same tips
	og <- read.tree(paste(arguments$folder,i,sep="/"))
	og.pr <- drop.tip(og,setdiff(og$tip.label,ref$tip.label))

	ref.pr <- drop.tip(ref,setdiff(ref$tip.label,og$tip.label))		

	norm.pr <- drop.tip(norm,setdiff(norm$tip.label,og.pr$tip.label))
	
	if (arguments$terminalbranches) {
	# This method calculates correlations between terminal branches of the trees

		# Extract and order terminal branch lengths for og and reference tree
		og.tblen <- as.matrix(setNames(og.pr$edge.length[sapply(1:length(og.pr$tip.label),function(x,y) which (y==x),y=og.pr$edge[,2])],og.pr$tip.label))
		og.tblen.order <- as.matrix(og.tblen[match(norm.pr$tip.label,row.names(og.tblen)),])
		#og.tblen.order <- og.tblen[order(row.names(og.tblen)),,drop=F]
		ref.tblen <- as.matrix(setNames(ref.pr$edge.length[sapply(1:length(ref.pr$tip.label),function(x,y) which (y==x),y=ref.pr$edge[,2])],ref.pr$tip.label))
		ref.tblen.order <- as.matrix(ref.tblen[match(norm.pr$tip.label,row.names(ref.tblen)),])
		#ref.tblen.order <- ref.tblen[order(row.names(ref.tblen)),,drop=F]
		
		# Fit a linear model and output p-value and r coefficient (adjusting the distribution for phylogenesis, based on the normalization tree)
		y <- pic(og.tblen.order,norm.pr)
		x <- pic(ref.tblen.order,norm.pr)
		model <- summary(lm( y ~ x - 1)) # "-1" is necessary since these are pic-transformed values
		erc[nrow(erc)+1,1] = gsub("\\.aln.*", "",i)
		erc[nrow(erc),2] = sqrt(model$r.squared)
		erc[nrow(erc),3] = model$coefficients[,4]
		
	} else { if (arguments$root2tip) {
	# The methods specified calculates correlations between the root-to-tip distances of each tip		

		# Calculate and order root-to-tip distances for og and reference trees
		og.dist <- as.matrix(distRoot(og.pr))
		og.order <- as.matrix(og.dist[match(norm.pr$tip.label,row.names(og.dist)),])
		#og.order <- og.dist[order(row.names(og.dist)),,drop=F]		
		ref.dist <- as.matrix(distRoot(ref.pr))
		ref.order <- as.matrix(ref.dist[match(norm.pr$tip.label,row.names(ref.dist)),])
		#ref.order <- ref.dist[order(row.names(ref.dist)),,drop=F]

		# Fit a linear model and output p-value and r coefficient (adjusting the distribution for phylogenesis, based on the normalization tree)
		y <- pic(og.order,norm.pr)
		x <- pic(ref.order,norm.pr)
		model <- summary(lm( y ~ x - 1)) # "-1" is necessary since these are pic-transformed values
		erc[nrow(erc)+1,1] = gsub("\\.aln.*", "",i)
		erc[nrow(erc),2] = sqrt(model$r.squared)
		erc[nrow(erc),3] = model$coefficients[,4]

		} 
	} 

}

## Adjust pvalues with different methods and output final table
erc$BH_fdr <- p.adjust(erc$pvalue,method="BH")
erc$BY_fdr <- p.adjust(erc$pvalue,method="BY")
erc$Hommel_cpvalue <- p.adjust(erc$pvalue,method="hommel")
erc$Holm_cpvalue <- p.adjust(erc$pvalue,method="holm")

write.table(erc,file=paste("ERC_analyses/",arguments$output,sep=""),quote=F,row.names=F,col.names=T,sep="\t")

