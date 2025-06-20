#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(ape)
#library(bnpsd) #tree_reoder modified to avoid checking for blen>1
source("../tree_order_modify.R") #from bnpsd - sourced directly due to problems in library installation

inalt <- paste(args[1],"_backbone_tagged.nwk",sep="")
outalt <- paste(args[1],"_alternative.nwk",sep="")
innull <- paste(args[1],"_blen.nwk",sep="")
outnull <- paste(args[1],"_null.nwk",sep="")

alt <- read.tree(inalt)
null <- read.tree(innull)

null.r <- root(null,sub("#1","",head(alt$tip.label,1)),resolve.root=TRUE)
alt.r <- root(alt,head(alt$tip.label,1),resolve.root=TRUE)
null.r <- tree_reorder_mod(null.r,sub("#1","",alt.r$tip.label))
alt.r.order <- ladderize(alt.r)
alt.r.order$edge.length <- ladderize(null.r)$edge.length


write.tree(unroot(alt.r.order),file=outalt)
write.tree(unroot(null.r),file=outnull)

